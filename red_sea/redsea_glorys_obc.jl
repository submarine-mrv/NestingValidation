# redsea_glorys_obc.jl
# GLORYS-driven open boundary for the Red Sea + Gulf of Aden.
# Adapted from salish_sea/salish_glorys_obc.jl; see docs/obc-implementation-recipe.md for the
# reusable pattern and salish_sea/README.md for why each piece is here.
#
# The point of this domain: ONE open boundary, in deep open ocean the parent resolves well, with
# Bab el-Mandeb (~36 km, ~4 GLORYS cells) INSIDE the domain so the model resolves the strait its
# parent cannot. That is the nesting question the Salish Sea configuration could never test.

using CairoMakie
using Oceananigans
using NumericalEarth
using CopernicusMarine   # loads NumericalEarthCopernicusMarineExt -> the GLORYS download backend
using Oceananigans: architecture
using Oceananigans.Advection: cell_advection_timescale
using Oceananigans.Units
using CFTime
using Dates
using Printf

const DW = NumericalEarth.DataWrangling
const SWP = Oceananigans.BuoyancyFormulations.SeawaterPolynomials

arch = CPU()

# --- Run window -------------------------------------------------------------
run_start   = DateTime(1993, 1, 1)
# run_days sizes the PARENT DATA window; REDSEA_STOP_DAYS sizes the simulation. They were decoupled,
# so a 2-day arm still loaded 10 days of GLORYS and held all of it resident (trap 3 forbids a smaller
# time_indices_in_memory here -- this script repairs the parent in memory with fill_down_fts!, and a
# rolling window silently reloads unrepaired slices). That footprint is what the OS killed twice.
# Shortening it is safe for comparison between arms: the daily samples covering t are the same files
# either way, so boundary values for t <= run_days are identical to a longer window's.
run_days    = parse(Int, get(ENV, "REDSEA_RUN_DAYS", "7"))
margin_days = 2
boundary_dates = run_start : Day(1) : (run_start + Day(run_days + margin_days))

# --- Grid --------------------------------------------------------------------
# Edges verified against ETOPO: south (10.0N), north (30.0N) and west (32.0E) are ALL LAND, so the
# EAST edge at 50.0E is the only open boundary. Suez and Aqaba terminate naturally inside.
# The 50E section is a single segment, 334 cells over 11.52-14.85N, 2848 m deep, 719 km^2 -- about
# 37 GLORYS cells across, against 1-2 for the Salish Sea channels.
λ₁, λ₂ = (32.0, 50.0)
φ₁, φ₂ = (10.0, 30.0)

# Resolution. Measured: 0.05 deg = 5.23 km, 360x400 = 144k columns, ~28 h wall for 7 days on CPU.
# 0.025 deg = 2.62 km, 576k columns, ~20x the Salish Sea cell count -> ~8-9 DAYS for 7 days on CPU.
# 2.5 km is the right target but needs GPU/HPC; develop at 0.05 and scale up.
#   REDSEA_DELTA=0.025 julia --project=<repo> red_sea/redsea_glorys_obc.jl
Δdeg = parse(Float64, get(ENV, "REDSEA_DELTA", "0.05"))
Nx = round(Int, (λ₂ - λ₁) / Δdeg)
Ny = round(Int, (φ₂ - φ₁) / Δdeg)

# Vertical. Measured max depth 3176 m, and GLORYS has valid data to 3226 m in the boundary box --
# only ~50 m of margin, so do NOT round `depth` up: model cells below the parent's valid data are
# what destabilised the Salish Sea configuration.
# `depth` stays 3200 even though the RED SEA reaches only 2478 m (median 425 m), because the 3176 m
# maximum is in the Gulf of Aden beside the open boundary; clipping it would truncate the deep water
# the boundary exists to supply. Only 0.7% of columns exceed 2500 m, so those levels are cheap to
# leave coarse -- GLORYS salinity is 34.776-34.802 over 2871-3226 m, i.e. homogeneous.
# Nz/scale protect the TOP 150 m: the Hanish sill is 136 m, so the entire Bab el-Mandeb exchange
# lives above it. 32/480 puts 17 levels there with a 0.94 m surface cell.
Nz = 32
depth = 3200
z_scale = 480          # prints Δz near the surface below so it is checkable
z_reference = ExponentialDiscretization(Nz, -depth, 0; scale = z_scale)

# z* (free-surface-following). Wrapping the reference discretization in MutableVerticalDiscretization
# is the only change needed: ocean_simulation does not pass `vertical_coordinate`, so the model takes
# default_vertical_coordinate(grid) = ZStarCoordinate() iff grid.z is mutable.
#
# ON by default here because of the 0.94 m surface cell, NOT because of tides. A static grid in
# Oceananigans is a LINEAR free surface: thicknesses are fixed and η never changes them
# (reconcile_vertical_coordinate! is a no-op, ∂t_σ = 0). It does not fail when η exceeds the top
# cell -- it just never lets η alter the geometry. The cost is accuracy: surface tracer concentration
# and the depth over which surface fluxes are spread each carry an O(max|η| / Δz_top) error. Here
# max|η| is 0.4-1.0 m against a 0.94 m top cell, so that ratio is O(1).
# Evaporation is this basin's driver, which sharpens it further: a surface freshwater flux over a
# fixed-volume thin cell cannot concentrate salt in a thinning layer, which is the physics.
# REDSEA_ZSTAR=0 falls back to static z -- if you do that, raise z_scale to ~600 (2.8 m top cell),
# or the surface-layer budget is not meaningful.
use_zstar = get(ENV, "REDSEA_ZSTAR", "1") == "1"
z = use_zstar ? MutableVerticalDiscretization(z_reference) : z_reference

@info @sprintf("Grid: %d x %d x %d  (Δ = %.4f deg ≈ %.2f km at 20N)  |  vertical: %s",
               Nx, Ny, Nz, Δdeg, Δdeg * 111.32 * cosd(20.0),
               use_zstar ? "ZStar (mutable)" : "Z (static)")

underlying_grid = LatitudeLongitudeGrid(arch;
                             size = (Nx, Ny, Nz),
                             latitude  = (φ₁, φ₂),
                             longitude = (λ₁, λ₂),
                             z,
                             halo = (7, 7, 7),
                             topology = (Bounded, Bounded, Bounded))

# Fewer smoothing passes than the Salish config's 10. Bab el-Mandeb's main channel is ~36 km, i.e.
# ~7 cells at 0.05 deg; 10 passes closed a comparable Salish channel entirely. Check the strait stays
# open (the diagnostic below reports its cross-section) before trusting a run.
interpolation_passes = parse(Int, get(ENV, "REDSEA_PASSES", "5"))
bottom_height = regrid_bathymetry(underlying_grid;
                                  minimum_depth = 10,
                                  interpolation_passes = interpolation_passes)
grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bottom_height))

# GLORYS subset: fixed, not derived from the domain edges, so resolution changes do not invalidate
# the download cache (the cached filenames encode the region).
glorys_region = DW.BoundingBox(longitude = (31.5, 50.5), latitude = (9.5, 30.5))
@assert λ₁ >= 31.5 && λ₂ <= 50.5 && φ₁ >= 9.5 && φ₂ <= 30.5 "domain outside glorys_region"
glorys = DW.GLORYSDaily()

# --- Boundary conditions -----------------------------------------------------
using Oceananigans.BoundaryConditions: NormalRadiation,
                                       GravityWaveRadiationBoundaryCondition,
                                       SurfaceWaveRadiationBoundaryCondition
using Oceananigans.Operators: Δzᶜᶜᶜ, Δzᶠᶜᶜ, Axᶠᶜᶜ, Ayᶜᶠᶜ, Azᶜᶜᶜ, Vᶜᶜᶜ
using Oceananigans.ImmersedBoundaries: immersed_inactive_node, immersed_peripheral_node
using Oceananigans.Units: Time

# Scheme selection. Defaults reproduce the validated 7-day runs exactly: NormalRadiation on
# everything, no reservoir. The alternatives are the Mid-Atlantic Bight pattern that canary/ settled
# on -- oblique on VELOCITIES, TracerReservoir on T/S -- see docs/obc-implementation-recipe.md.
const HAS_OBLIQUE   = isdefined(Oceananigans.BoundaryConditions, :ObliqueRadiation)
const HAS_RESERVOIR = isdefined(Oceananigans.BoundaryConditions, :TracerReservoir)

const RAD_NORMAL = NormalRadiation(inflow_timescale = 6hours, outflow_timescale = Inf)

# REDSEA_RADIATION=oblique selects Raymond & Kuo (1984) two-dimensional radiation (Oceananigans
# >= 0.113). Only ever applied to velocities here: oblique on TRACERS at a boundary that terminates
# in land drove salinity to 45.7 psu in canary/, and the east boundary is cut by the Yemeni and
# Somali coasts at both ends.
radiation_kind = get(ENV, "REDSEA_RADIATION", "normal")
radiation_kind in ("normal", "oblique") || error("REDSEA_RADIATION must be \"normal\" or \"oblique\"")
radiation_kind == "oblique" && !HAS_OBLIQUE &&
    error("REDSEA_RADIATION=oblique needs Oceananigans >= 0.113; this environment has no ObliqueRadiation")
radiation = radiation_kind == "oblique" ?
    getfield(Oceananigans.BoundaryConditions, :ObliqueRadiation)(inflow_timescale = 6hours, outflow_timescale = Inf) :
    RAD_NORMAL

# T/S carry their own scheme. The reservoir (MOM6-style, metres) gives the boundary a memory: a
# memoryless Value condition hands re-entering water the exterior value regardless of what just
# left. Bab el-Mandeb is a two-layer exchange -- surface inflow, deep outflow -- so water leaving
# and returning is the NORMAL state here, not an edge case.
const res_in  = parse(Float64, get(ENV, "REDSEA_RESERVOIR_IN",  "0"))
const res_out = parse(Float64, get(ENV, "REDSEA_RESERVOIR_OUT", "0"))
ts_kind = get(ENV, "REDSEA_TS_SCHEME", "same")
ts_kind in ("same", "normal", "reservoir") || error("REDSEA_TS_SCHEME must be same/normal/reservoir")
ts_kind == "reservoir" && !HAS_RESERVOIR &&
    error("REDSEA_TS_SCHEME=reservoir needs Oceananigans >= 0.113.2")
ts_scheme = ts_kind == "normal"    ? RAD_NORMAL :
            ts_kind == "reservoir" ? getfield(Oceananigans.BoundaryConditions, :TracerReservoir)(
                                         inflow_length_scale = res_in, outflow_length_scale = res_out) :
            radiation
# REDSEA_TANGENTIAL adds a boundary-condition for the ALONG-boundary velocity, which this domain has
# never had: v was left to a default no-flux wall on an open edge. OFF by default, like every other
# knob here, so the defaults still reproduce the validated runs exactly -- the +32-39% retained
# enstrophy that justified it was measured in canary/, not in this domain, and a baseline arm is only
# worth running if it changes one thing at a time.
use_tangential = get(ENV, "REDSEA_TANGENTIAL", "0") == "1"
@info "Open-boundary scheme: velocities $(radiation_kind), T/S $(ts_kind), tangential $(use_tangential ? "on" : "off")" *
      (ts_kind == "reservoir" ? @sprintf(" (reservoir in %.0f m, out %.0f m)", res_in, res_out) : "")

@inline wetcell(i, j, k, grid, ℓx, ℓy, ℓz) =
    !immersed_peripheral_node(i, j, k, grid, ℓx, ℓy, ℓz) & !immersed_inactive_node(i, j, k, grid, ℓx, ℓy, ℓz)

# --- Nudging band ------------------------------------------------------------
# Single band, along the east edge only. 0.5 deg (~10 cells at 0.05 deg) sits far from Bab el-Mandeb
# (43.3E), so the strait is left free -- the Salish lesson was that a wide sponge pins the very
# passage the domain exists to resolve.
sponge_width        = 0.5
restoring_timescale = 12hours
restoring_rate      = 1 / restoring_timescale

@inline gaussian_edge(d, w) = exp(-(d / w)^2)
@inline boundary_mask(x, y, z, t) = gaussian_edge(λ₂ - x, sponge_width)

restoring_meta(var) = DW.Metadata(var; dataset = glorys, dates = boundary_dates, region = glorys_region)
bc_inpainting = DW.NearestNeighborInpainting(100)

# Repair parent levels with no valid data (they come back zero-filled, and interpolation then blends
# a valid value with zero to give anything in between -- an absolute floor cannot catch that).
#
# In THIS domain it is insurance, not load-bearing: glorys_interior_check.jl re-ran the whole IC with
# fill_down! disabled and got a bit-identical model state. It does repair ~630k native cells, but all
# of them lie below -3226 m or in columns the model never samples, which is why zeros survive in the
# printed range below and are still harmless. Keep it -- the next domain may need it, and the Salish
# Sea did. But the real interior risk here is HORIZONTAL: inpainting reaches up to 100 cells sideways
# and imports the wrong water mass. A vertical repair cannot fix that. See red_sea/README.md.
function fill_down!(a::AbstractArray{<:Any,3}, floor_value)
    nx, ny, nz = size(a)
    @inbounds for j in 1:ny, i in 1:nx
        last = NaN
        for k in nz:-1:1                    # k = nz is the surface
            v = a[i, j, k]
            if isfinite(v) && v > floor_value
                last = v
            elseif isfinite(last)
                a[i, j, k] = last
            end
        end
    end
    return a
end
fill_down_fts!(fts, floor) = (for n in 1:length(fts.times); fill_down!(interior(fts[n]), floor); end; fts)

# All time indices resident: the default of 2 reloads slices from disk as the clock advances, which
# would silently undo the repair above.
n_times = length(boundary_dates)
# time_indexing: DatasetRestoring defaults to Cyclical(), and these series are not only the interior
# sponge -- their field_time_series is what the open-boundary exterior values are sampled from. Under
# NumericalEarth 0.8 GLORYSDaily is stamped at the CENTRE of its daily-mean window (12:00), so the
# first sample lands 12 h AFTER the run starts, and Cyclical fills that gap by wrapping to the FAR END
# of the record: the last day of the series would drive the first half-day of the run, at both the
# boundary and the sponge. Clamp holds the first sample instead (persistence), matching what the
# hand-rolled kernels in canary/canary_glorys_obc.jl do. A lead-in day is not available: GLORYS
# itself begins 1993-01-01, which is run_start.
restoring_kw = (; rate = restoring_rate, mask = boundary_mask, inpainting = bc_inpainting,
                  time_indices_in_memory = n_times,
                  time_indexing = Oceananigans.OutputReaders.Clamp())

T_restoring = DW.DatasetRestoring(restoring_meta(:temperature), arch; restoring_kw...)
S_restoring = DW.DatasetRestoring(restoring_meta(:salinity),    arch; restoring_kw...)
u_restoring = DW.DatasetRestoring(restoring_meta(:u_velocity),  arch; restoring_kw...)
v_restoring = DW.DatasetRestoring(restoring_meta(:v_velocity),  arch; restoring_kw...)

# Red Sea S runs 35-41 (hypersaline in the north), T 21-31 -- floors well below anything physical.
fill_down_fts!(T_restoring.field_time_series, 1.0)
fill_down_fts!(S_restoring.field_time_series, 5.0)
# Zeros here are expected (see above) -- what matters is that none reach a wet model cell, which the
# "Initial condition (wet)" line below reports.
let Ti = interior(T_restoring.field_time_series[1]), Si = interior(S_restoring.field_time_series[1])
    @info @sprintf("GLORYS native after fill_down: T ∈ [%.3f, %.3f], S ∈ [%.3f, %.3f] | residual zeros %.2f%% (below -3226 m, unsampled)",
                   minimum(Ti), maximum(Ti), minimum(Si), maximum(Si), 100count(==(0), Si) / length(Si))
end

T_fts, T_native = T_restoring.field_time_series, T_restoring.native_grid
S_fts, S_native = S_restoring.field_time_series, S_restoring.native_grid
u_fts, u_native = u_restoring.field_time_series, u_restoring.native_grid
v_fts, v_native = v_restoring.field_time_series, v_restoring.native_grid

# --- GLORYS section residual (diagnostic only) --------------------------------
# There is NO net-transport correction. The raw GLORYS exterior carries a nonzero net transport
# through the east section, but correcting it turned out not to matter: two 2-day runs whose
# prescribed exterior nets differed by 3.0e5 m³/s were indistinguishable (strait exchange within
# 0.6%, max|η| and max|u| unchanged, volume drift within 1.5%), and the no-correction case lies
# BETWEEN those two, so it is already bracketed.
#
# The reason is that Flather is a restoring feedback, not an integrator: U = Uᵉˣᵗ + √(gH)(ηᵇ − ηᵉˣᵗ),
# so η settles ~4 mm below the exterior and sources the residual indefinitely without accumulating.
# The same holds for evaporation, whose demand (~1.3e4 m³/s) is an order of magnitude SMALLER than
# the residual and needs ~0.2 mm. Removing the correction also removes a global Ref mutated by a
# callback and read inside boundary-condition closures -- machinery that carried a sign error through
# an entire validation run with no detectable symptom. See docs/obc-implementation-recipe.md.
#
# Q_raw is still worth measuring: it is a PARENT/CHILD MISMATCH, not a physical signal, and it flips
# sign with resolution (+3.15e5 m³/s at 0.2°, -1.49e5 at 0.05°) because parent and child disagree
# about bathymetry and cell geometry. `progress` reports it.

# East boundary: positive u is OUT of the domain, so inflow is -u.
@inline u_east_exterior(j, k, grid, clock, fields) =
    DW.interpolate_to_grid(u_fts, grid.Nx, j, k, u_native, grid, Time(clock.time))

# Tangential (along-boundary) velocity. Absent until now: the boundary-normal component was
# prescribed and v left to a default no-flux wall, which is wrong for an open edge and measurably
# so -- in an idealised truncation test adding it was worth +32-39% of retained enstrophy.
# Switched by REDSEA_TANGENTIAL; see the scheme block above.
@inline v_east_tangential(j, k, grid, clock, fields) =
    DW.interpolate_to_grid(v_fts, grid.Nx, j, k, v_native, grid, Time(clock.time))

@inline T_east_exterior(j, k, grid, clock, fields) =
    DW.interpolate_to_grid(T_fts, grid.Nx, j, k, T_native, grid, Time(clock.time))
@inline S_east_exterior(j, k, grid, clock, fields) =
    DW.interpolate_to_grid(S_fts, grid.Nx, j, k, S_native, grid, Time(clock.time))

# GLORYS sea surface height as the Flather exterior η. Held fully in memory: unlike the T/S/u/v
# series this one is NOT reachable from model.forcing, so update_model_field_time_series! cannot
# advance its window -- a closure is opaque to that traversal.
zos_meta = DW.Metadata(:free_surface; dataset = glorys, dates = boundary_dates, region = glorys_region)
zos_fts  = FieldTimeSeries(zos_meta, arch;
                           time_indices_in_memory = n_times,
                           # Clamp, NOT Cyclical. NumericalEarth 0.8 stamps GLORYSDaily at the CENTRE of
                           # its daily-mean window (12:00), so the first sample lands 12 h AFTER the run
                           # starts. Cyclical() fills that gap by wrapping to the FAR END of the record --
                           # the last day of the series forcing the first half-day of the run. Clamp holds
                           # the first sample instead, i.e. persistence, which is what the hand-rolled
                           # kernels in canary/canary_glorys_obc.jl already do. A lead-in day is not an
                           # option: GLORYS itself begins 1993-01-01, the usual run_start.
                           time_indexing = Oceananigans.OutputReaders.Clamp(),
                           inpainting = bc_inpainting)
zos_native = zos_fts.grid
@inline zos_at(i, j, grid, t) = DW.interpolate_to_grid(zos_fts, i, j, 1, zos_native, grid, t)

# Flather at an EAST boundary:  U = Uᵉˣᵗ + √(gH)(ηᵇ − ηᵉˣᵗ)   (sign is opposite the west case).
# Returns (exterior transport, exterior η).
@inline function east_U_obc(j, k, grid, clock, fields)
    i = grid.Nx
    t = Time(clock.time)
    U = zero(grid)
    @inbounds for kk in 1:grid.Nz
        wet = wetcell(i, j, kk, grid, Face(), Center(), Center())
        uk  = DW.interpolate_to_grid(u_fts, i, j, kk, u_native, grid, t)
        U  += ifelse(wet, uk * Δzᶠᶜᶜ(i, j, kk, grid), zero(grid))
    end
    return (U, zos_at(i, j, grid, t))
end

u_bcs = FieldBoundaryConditions(east = NormalFlowBoundaryCondition(u_east_exterior; scheme = radiation, discrete_form = true))
v_bcs = use_tangential ?
    FieldBoundaryConditions(east = ValueBoundaryCondition(v_east_tangential; scheme = radiation, discrete_form = true)) :
    nothing
U_bcs = FieldBoundaryConditions(east = GravityWaveRadiationBoundaryCondition(east_U_obc; discrete_form = true))
# The Chapman condition radiates interior η outward and never reads an exterior value -- see the
# recipe. Present for the barotropic mode to be well posed, not to impose anything.
η_bcs = FieldBoundaryConditions(grid, (Center(), Center(), Face()); east = SurfaceWaveRadiationBoundaryCondition())
T_bcs = FieldBoundaryConditions(east = ValueBoundaryCondition(T_east_exterior; scheme = ts_scheme, discrete_form = true))
S_bcs = FieldBoundaryConditions(east = ValueBoundaryCondition(S_east_exterior; scheme = ts_scheme, discrete_form = true))

# --- Ocean simulation --------------------------------------------------------
weno_order = parse(Int, get(ENV, "REDSEA_WENO_ORDER", "7"))
tracer_advection = WENO(order = weno_order,
                        time_discretization = AdaptiveVerticallyImplicitDiscretization(cfl = 0.5))

ocean = ocean_simulation(grid;
                         tracers = (:T, :S),
                         tracer_advection,
                         forcing = (T = T_restoring, S = S_restoring,
                                    u = u_restoring, v = v_restoring),
                         boundary_conditions = merge((; u = u_bcs, U = U_bcs, η = η_bcs, T = T_bcs, S = S_bcs),
                                                     isnothing(v_bcs) ? (;) : (; v = v_bcs)))

# --- Atmospheric forcing -----------------------------------------------------
# API note: this differs from older ClimaOcean/NumericalEarth examples. In 0.6.1 radiation is a
# SEPARATE top-level component (SW/LW left the atmosphere), there is no `JRA55NetCDFBackend` --
# `time_indices_in_memory` replaces it -- and `Simulation` is Oceananigans', not NumericalEarth's.
# Rivers/icebergs are simply not built by the constructor, so there is no kwarg to switch them off.
#
# Dataset: RepeatYearJRA55 by default because it is already cached (~12 GB) and its year runs
# 1990-01-01..12-31, so mapping run_start onto the same month/day keeps the SEASON aligned with the
# 1993 GLORYS ocean. It is a repeating climatological year, so the synoptic weather does NOT
# correspond to 1993. REDSEA_JRA_DATASET=multiyear uses the real 1993 fields instead -- correct, but
# it is another ~12 GB download and the disk had only ~12 GB free, so check before switching.
atmos_source = get(ENV, "REDSEA_ATMOS", "JRA")     # "JRA" | "none"
use_atmos = atmos_source == "JRA"
jra_multiyear = get(ENV, "REDSEA_JRA_DATASET", "repeat") == "multiyear"

atmosphere = nothing
radiation  = nothing
if use_atmos
    jra_dataset = jra_multiyear ? NumericalEarth.MultiYearJRA55() : NumericalEarth.RepeatYearJRA55()
    # Repeat-year data lives in 1990; shift the run window onto it so t = 0 is the same day of year.
    # Qualify: Oceananigans.Units also exports `day`/`year`, which shadows the Dates accessors.
    jra_year   = jra_multiyear ? Dates.year(run_start) : 1990
    jra_start  = DateTime(jra_year, Dates.month(run_start), Dates.day(run_start))
    jra_stop   = jra_start + Day(run_days + margin_days)
    @info @sprintf("JRA55 %s: %s .. %s", jra_multiyear ? "MultiYear" : "RepeatYear", jra_start, jra_stop)

    atmos_kw = (; dataset = jra_dataset, start_date = jra_start, end_date = jra_stop,
                  time_indices_in_memory = 24)     # 3-hourly, so 24 slices = 3 days resident
    atmosphere = NumericalEarth.JRA55PrescribedAtmosphere(arch; atmos_kw...)
    radiation  = NumericalEarth.JRA55PrescribedRadiation(arch; atmos_kw...)
end

# --- Initial conditions ------------------------------------------------------
# Sampled from the SAME repaired FieldTimeSeries the nudging and the boundary use, so IC and forcing
# agree by construction. The alternative (`set!(field, ::Metadatum)`) re-reads the cached file and
# never sees the repair -- that mismatch was the root cause of the Salish Sea instability.
# Written to scratch fields then copied: assigning to model tracers directly triggers a halo fill
# through the clock-dependent OBCs, which a plain set! cannot satisfy.
function ic_from_fts!(field, fts, native_grid, t)
    a = interior(field)
    @inbounds for k in 1:Nz, j in 1:Ny, i in 1:Nx
        a[i, j, k] = DW.interpolate_to_grid(fts, i, j, k, native_grid, grid, t)
    end
    return field
end

T_ic = CenterField(grid); ic_from_fts!(T_ic, T_fts, T_native, Time(0.0))
S_ic = CenterField(grid); ic_from_fts!(S_ic, S_fts, S_native, Time(0.0))
interior(ocean.model.tracers.T) .= interior(T_ic)
interior(ocean.model.tracers.S) .= interior(S_ic)

free_surface_field(fs) = hasproperty(fs, :η) ? fs.η : fs.displacement
η_field = free_surface_field(ocean.model.free_surface)

# η from GLORYS zos rather than rest: starting at 0 against a boundary saying ~0.5 m manufactures a
# basin-filling transient that swamps any bounded test.
let ηi = interior(η_field), t₀ = Time(0.0)
    for i in 1:Nx, j in 1:Ny
        wet = !immersed_inactive_node(i, j, Nz, grid, Center(), Center(), Center())
        ηi[i, j, 1] = ifelse(wet, zos_at(i, j, grid, t₀), 0.0)
    end
end
@info @sprintf("η from GLORYS zos: [%.4f, %.4f] m", minimum(interior(η_field)), maximum(interior(η_field)))

# --- Coupled model -----------------------------------------------------------
# ocean_simulation merges user boundary conditions SIDE BY SIDE with its defaults (every side left
# unspecified inherits the default), so the east-only T/S conditions above keep the default TOP flux
# boundary condition and the atmospheric fluxes do reach the ocean. Verified in
# Oceans/ocean_simulation.jl:merge_boundary_conditions -- worth re-checking on any version bump,
# because a replace-rather-than-merge would silently disable surface forcing.
if use_atmos
    coupled_model = NumericalEarth.OceanOnlyModel(ocean; atmosphere, radiation)
    simulation = Simulation(coupled_model; Δt = 10, stop_iteration = 100)
else
    simulation = ocean
end

# `coupled_model.ocean` is the ocean SIMULATION, not the model; unforced runs hand `progress` the
# ocean model directly. Everything below goes through this so both paths share one code path.
ocean_model_of(m) = hasproperty(m, :ocean) ? (hasproperty(m.ocean, :model) ? m.ocean.model : m.ocean) : m

# --- Diagnostics -------------------------------------------------------------
z_centers = znodes(grid, Center(), Center(), Center())
λ_centers = λnodes(grid, Center(), Center(), Center())
φ_centers = φnodes(grid, Center(), Center(), Center())

@info @sprintf("Δz near surface: %s m (top of %d levels, domain %d m)",
               join((@sprintf("%.1f", d) for d in reverse(diff(z_centers))[1:min(6, Nz-1)]), ", "), Nz, depth)

wet3d = [!immersed_inactive_node(i, j, k, grid, Center(), Center(), Center()) for i in 1:Nx, j in 1:Ny, k in 1:Nz]
wet_idx = findall(wet3d)
@info "Wet cells: $(length(wet_idx)) of $(Nx*Ny*Nz) ($(round(100*length(wet_idx)/(Nx*Ny*Nz), digits=1))%)"

wet_surface_area = sum(Azᶜᶜᶜ(i, j, Nz, grid) for i in 1:Nx, j in 1:Ny
                       if !immersed_inactive_node(i, j, Nz, grid, Center(), Center(), Center()); init = 0.0)

# Bab el-Mandeb cross-section, so a run reports whether the strait survived the bathymetry smoothing
# and what it is carrying. j index nearest the narrowest channel (12.69N).
j_strait = argmin(abs.(φ_centers .- 12.69))
strait_cells = [i for i in 1:Nx if λ_centers[i] < 44.5 &&
                wetcell(i, j_strait, Nz, grid, Center(), Center(), Center())]
if isempty(strait_cells)
    @warn "Bab el-Mandeb section at 12.69N is CLOSED at this resolution/smoothing -- reduce REDSEA_PASSES"
else
    @info @sprintf("Bab el-Mandeb (12.69N): %d wet cells, lon %.2f..%.2f, deepest %.0f m",
                   length(strait_cells), λ_centers[strait_cells[1]], λ_centers[strait_cells[end]],
                   maximum(-Array(interior(bottom_height))[i, j_strait, 1] for i in strait_cells))
end

function wet_extrema(field)
    a = interior(field)[wet_idx]
    nmin, nmax = argmin(a), argmax(a)
    return a[nmin], wet_idx[nmin], a[nmax], wet_idx[nmax]
end

function whereis(I)
    i, j, k = Tuple(I)
    zz = z_centers[clamp(k, 1, Nz)]
    tag = i >= Nx - 1 ? "EAST-bnd" :
          i >= Nx - 10 ? "near-east" :
          (abs(φ_centers[clamp(j,1,Ny)] - 12.69) < 1.0 && λ_centers[clamp(i,1,Nx)] < 44.5) ? "BAB-EL-MANDEB" : "interior"
    return @sprintf("(%3d,%3d,%2d) z=%7.1fm %s", i, j, k, zz, tag)
end

# East section: positive u is OUT of the domain, so inflow is -u.
function east_transport(model)
    u = model.velocities.u
    Φ = 0.0
    for k in 1:Nz, j in 1:Ny
        wetcell(Nx, j, k, grid, Face(), Center(), Center()) || continue
        @inbounds Φ -= u[Nx, j, k] * Axᶠᶜᶜ(Nx, j, k, grid)
    end
    return Φ
end

# Exchange across Bab el-Mandeb (zonal section at j_strait, so the normal component is v).
# THIS is the number to compare with observations -- roughly 0.3-0.4 Sv each way with a seasonally
# varying 2-3 layer structure. The boundary exchange below measures boundary activity instead: at
# 50E nothing has been transformed yet, so its ΔS is near zero by construction.
function strait_exchange(model)
    v, S = model.velocities.v, model.tracers.S
    Qin = Qout = SQin = SQout = 0.0
    for k in 1:Nz, i in strait_cells
        wetcell(i, j_strait, k, grid, Center(), Face(), Center()) || continue
        @inbounds q = v[i, j_strait, k] * Ayᶜᶠᶜ(i, j_strait, k, grid)   # + = northward, into the Red Sea
        @inbounds s = S[i, j_strait, k]
        q > 0 ? (Qin += q; SQin += q*s) : (Qout += q; SQout += q*s)
    end
    return Qin, Qout, (Qin > 0 ? SQin/Qin : 0.0), (Qout < 0 ? SQout/Qout : 0.0)
end

function east_exchange(model)
    u, S = model.velocities.u, model.tracers.S
    Qin = Qout = SQin = SQout = 0.0
    for k in 1:Nz, j in 1:Ny
        wetcell(Nx, j, k, grid, Face(), Center(), Center()) || continue
        @inbounds q = -u[Nx, j, k] * Axᶠᶜᶜ(Nx, j, k, grid)   # + into the domain
        @inbounds s = S[Nx, j, k]
        q > 0 ? (Qin += q; SQin += q*s) : (Qout += q; SQout += q*s)
    end
    return Qin, Qout, (Qin > 0 ? SQin/Qin : 0.0), (Qout < 0 ? SQout/Qout : 0.0)
end

# ∫ Jʷ dA over wet surface cells. Jʷ lives on the coupled interfaces, so it is zero without an
# atmosphere. The basin-mean rate this implies is printed by `progress` in m/yr, which is the number
# to sanity-check: the Red Sea proper evaporates ~2 m/yr, and a domain half Gulf of Aden should still
# come out clearly negative. A POSITIVE value means the domain is gaining freshwater and the
# overturning will run backwards -- see red_sea/README.md.
function surface_freshwater_flux(m)
    use_atmos || return 0.0
    Jη = m.interfaces.net_fluxes.ocean.η
    Q = 0.0
    @inbounds for j in 1:Ny, i in 1:Nx
        immersed_inactive_node(i, j, Nz, grid, Center(), Center(), Center()) && continue
        Q += Jη[i, j, 1] * Azᶜᶜᶜ(i, j, Nz, grid)
    end
    return Q
end

# Net inflow the RAW GLORYS exterior would carry through the east section, and the section area.
# Diagnostic only -- nothing is corrected with it.
function glorys_section_residual(t)
    Q = A = 0.0
    for j in 1:Ny, k in 1:Nz
        wetcell(Nx, j, k, grid, Face(), Center(), Center()) || continue
        a = Axᶠᶜᶜ(Nx, j, k, grid)
        Q -= DW.interpolate_to_grid(u_fts, Nx, j, k, u_native, grid, t) * a
        A += a
    end
    return Q, A
end

let (Q, A) = glorys_section_residual(Time(0.0))
    @info @sprintf("GLORYS section residual at t=0: Q_raw = %+.4e m³/s over A = %.4e m² (Δu-equivalent %+.3e m/s); NOT corrected",
                   Q, A, A > 0 ? Q / A : 0.0)
end

# Excursion is referenced to the FORCING DATA range, not the IC: the boundary can legitimately import
# water outside the initial range, and referencing the IC reports physical adjustment as spurious.
# Valid-only: fill_down! propagates a value downward, so a column with NO valid value anywhere
# (deep inland, beyond inpainting reach across this much land) keeps its zeros. Those are harmless
# in themselves but would make the reference -- and therefore the excursion -- meaningless.
function native_range(fts, floor_value)
    a = filter(x -> x > floor_value, vec(interior(fts[1])))
    isempty(a) ? (NaN, NaN) : (minimum(a), maximum(a))
end
T_ref = native_range(T_fts, 1.0); S_ref = native_range(S_fts, 5.0)
T_ic_range = let (lo,_,hi,_) = wet_extrema(ocean.model.tracers.T); (lo,hi) end
S_ic_range = let (lo,_,hi,_) = wet_extrema(ocean.model.tracers.S); (lo,hi) end
@info @sprintf("Reference (GLORYS native): T ∈ [%.3f, %.3f] | S ∈ [%.3f, %.3f]", T_ref[1], T_ref[2], S_ref[1], S_ref[2])
@info @sprintf("Initial condition (wet):    T ∈ [%.3f, %.3f] | S ∈ [%.3f, %.3f]",
               T_ic_range[1], T_ic_range[2], S_ic_range[1], S_ic_range[2])

# --- Conservation and stability monitors -------------------------------------
# An INTEGRATED metric, because min/max S cannot see a slow drift. Volume is the other half of the
# same question: under ZStar the thicknesses move with η, so wet volume should track the net boundary
# transport and nothing else. This is what cleared z* on these package versions -- over 7 days mean S
# went 37.5030 -> 37.5028 with a non-monotonic 1.4e-3 psu excursion, and salt tracked volume
# throughout, i.e. concentration conserved under the moving grid. See red_sea/README.md.
function salt_and_volume(model)
    S = model.tracers.S
    V = 0.0; SV = 0.0
    @inbounds for k in 1:Nz, j in 1:Ny, i in 1:Nx
        wetcell(i, j, k, grid, Center(), Center(), Center()) || continue
        v = Vᶜᶜᶜ(i, j, k, grid)
        V += v; SV += v * S[i, j, k]
    end
    return SV, V
end

# Static stability. glorys_interior_check.jl found 10558 unstable pairs in the IC, worst 0.92 kg/m³,
# shallow and coastal, from horizontal inpainting rather than from fill_down!. Over the 7-day run they
# RISE first (to 13394 during the adjustment transient, which stirs faster than the closure mixes)
# and then decay to 7325 at worst 0.44 -- so an early sample reads as growth. They do not reach zero:
# ~7300 persist as a stirring-against-mixing balance. Reference pressure is the cell midpoint from the
# rest-state grid; σ departs from 1 by O(1e-3), far below what this comparison resolves.
const eos = SWP.TEOS10.TEOS10EquationOfState()

function static_instability(model)
    T, S = model.tracers.T, model.tracers.S
    n = 0; worst = 0.0
    @inbounds for j in 1:Ny, i in 1:Nx, k in 1:(Nz-1)
        (wetcell(i, j, k,   grid, Center(), Center(), Center()) &&
         wetcell(i, j, k+1, grid, Center(), Center(), Center())) || continue
        zr = 0.5 * (z_centers[k] + z_centers[k+1])
        Δρ = SWP.ρ(T[i, j, k+1], S[i, j, k+1], zr, eos) - SWP.ρ(T[i, j, k], S[i, j, k], zr, eos)
        if Δρ > 1e-4
            n += 1
            worst = max(worst, Δρ)
        end
    end
    return n, worst
end

# The drift baseline is taken on the FIRST PROGRESS CALL, not here. Under z* the grid still carries
# σ = 1 at this point while η already holds ~0.4 m of GLORYS zos; σ is reconciled inside the first
# time_step! (update_state! does not do it). Baselining at rest therefore books that one-time
# reconciliation as a +3.7e-4 volume jump, which is the same size as the drift being looked for.
const salt₀ = Ref(0.0)
const vol₀  = Ref(0.0)

let (SV, V) = salt_and_volume(ocean.model)
    n, w = static_instability(ocean.model)
    @info @sprintf("At rest: salt %.6e psu·m³ | wet volume %.6e m³ | mean S %.4f | unstable pairs %d (worst Δρ %.3f kg/m³)",
                   SV, V, SV/V, n, w)
end

const wall_t0   = Ref(Base.time())
const wall_last = Ref(Base.time())
const sim_last  = Ref(0.0)

# Divisor for the pacing line below. Must be the SAME value the callback is registered with, or
# the reported s/iter is wrong by their ratio -- it read 2.17 s/iter at interval 50 when the truth
# was 0.43, because the divisor was hardcoded to 10.
const PROGRESS_ITERS = parse(Int, get(ENV, "REDSEA_PROGRESS_ITERS", "10"))

function progress(sim)
    model = ocean_model_of(sim.model)
    u, v, w = model.velocities
    T, S = model.tracers.T, model.tracers.S
    η = free_surface_field(model.free_surface)

    Tmin, _, Tmax, _       = wet_extrema(T)
    Smin, Imin, Smax, Imax = wet_extrema(S)
    _, Iw = findmax(abs, interior(w))
    _, Iu = findmax(abs, interior(u))

    Φ = east_transport(model)
    @info @sprintf("iter %d | t = %s | Δt = %s | max|u,v,w| = (%.2e, %.2e, %.2e) | max|η| = %.3e | T ∈ [%.3f, %.3f] | S ∈ [%.3f, %.3f] (excursion %+.3f / %+.3f)",
                   iteration(sim), prettytime(time(sim)), prettytime(sim.Δt),
                   maximum(abs, u), maximum(abs, v), maximum(abs, w), maximum(abs, η),
                   Tmin, Tmax, Smin, Smax, Smin - S_ref[1], Smax - S_ref[2])
    @info @sprintf("        where: Smax %s | Smin %s | max|w| %s | max|u| %s",
                   whereis(Imax), whereis(Imin), whereis(Iw), whereis(Iu))
    Qraw, _ = glorys_section_residual(Time(model.clock.time))
    @info @sprintf("        transport: east %+.4e m³/s (+ = inflow)  ->  %+.3e m/hour mean sea level | GLORYS residual %+.4e (uncorrected)",
                   Φ, Φ / wet_surface_area * 3600, Qraw)

    Qin, Qout, Sin, Sout = east_exchange(model)
    @info @sprintf("        east boundary: Qin %+.3e  Qout %+.3e m³/s | ΔS %.3f  (boundary activity, not transformation)",
                   Qin, Qout, Sin - Sout)

    if !isempty(strait_cells)
        qi, qo, si, so = strait_exchange(model)
        # One limb can be empty early in spin-up; ΔS and the Knudsen residual are meaningless then.
        if qi > 0 && qo < 0
            salt_in = qi * si
            @info @sprintf("        BAB EL-MANDEB: Qin %+.3e  Qout %+.3e m³/s (lit. ~0.3-0.4 Sv) | Sin %.3f  Sout %.3f  ΔS %.3f | Knudsen %+.1f%%",
                           qi, qo, si, so, si - so, 100 * (salt_in + qo*so) / abs(salt_in))
        else
            @info @sprintf("        BAB EL-MANDEB: Qin %+.3e  Qout %+.3e m³/s | one-way flow, ΔS n/a", qi, qo)
        end
    end

    SV, V = salt_and_volume(model)
    ninv, worstΔρ = static_instability(model)
    # iteration > 0 is load-bearing: σ is reconciled inside the first time_step!, and under the
    # coupled model progress fires at iteration 0 BEFORE that happens. Baselining there books the
    # reconciliation as a +3.7e-4 volume jump -- the artifact this whole scheme exists to avoid.
    if vol₀[] == 0.0 && iteration(sim) > 0
        salt₀[], vol₀[] = SV, V
        @info @sprintf("        drift baseline set at t = %s (σ reconciled): salt %.6e | volume %.6e",
                       prettytime(time(sim)), SV, V)
    end
    if use_atmos
        Qfw = surface_freshwater_flux(sim.model)
        @info @sprintf("        surface freshwater: %+.4e m³/s (%+.3f m/yr basin mean)",
                       Qfw, Qfw / wet_surface_area * 3.15576e7)
    end
    have_base = vol₀[] > 0
    @info @sprintf("        conservation: salt %s | volume %s (relative drift vs baseline) | mean S %.4f | unstable pairs %d (worst Δρ %.3f)",
                   have_base ? @sprintf("%+.3e", (SV - salt₀[]) / salt₀[]) : "  (baseline pending)",
                   have_base ? @sprintf("%+.3e", (V  - vol₀[])  / vol₀[])  : "  (pending)",
                   SV / V, ninv, worstΔρ)

    now  = Base.time(); wall = now - wall_t0[]
    Δwall, Δsim = now - wall_last[], time(sim) - sim_last[]
    wall_last[], sim_last[] = now, time(sim)
    rate = Δwall > 0 ? Δsim / Δwall : 0.0
    eta  = (isfinite(sim.stop_time) && rate > 0) ? (sim.stop_time - time(sim)) / rate : NaN
    @info @sprintf("        pacing: wall %s | %.0f× realtime (recent) | %.2f s/iter | ETA %s",
                   prettytime(wall), rate, Δwall / PROGRESS_ITERS, isnan(eta) ? "n/a" : prettytime(eta))
    flush(stderr)
    return nothing
end
# progress walks every wet cell on the host -- wet_extrema on T and S, two findmax, the east and
# Bab el-Mandeb transports, and a conservation check that compares density pairs over all ~497k wet
# cells. At IterationInterval(10) that dominates wall time: the cost per iteration stays ~32 s whether
# Δt is 10 s or 15 min, which is the signature of a fixed per-call cost rather than of the dynamics.
# Thinning it does not touch the solution -- the diagnostics are read-only and the output writers keep
# their own schedule. Mirrors CANARY_PROGRESS_ITERS.
add_callback!(simulation, progress, IterationInterval(PROGRESS_ITERS))

# --- Spin-up then the run ----------------------------------------------------
simulation.Δt = 10
simulation.stop_iteration = 100
run!(simulation)

ocean_timescale(model) = cell_advection_timescale(ocean_model_of(model))
wizard = TimeStepWizard(; cfl = 0.15, max_Δt = 15minutes, max_change = 1.1, cell_advection_timescale = ocean_timescale)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))
simulation.stop_iteration = Inf
stop_days = parse(Float64, get(ENV, "REDSEA_STOP_DAYS", string(run_days)))
simulation.stop_time = stop_days * days

output_fields   = merge(ocean.model.velocities, ocean.model.tracers)
output_schedule = TimeInterval(6hours)
# Atmosphere is tagged so forced runs never overwrite the unforced 7-day validation baseline.
# REDSEA_TAG appends an arbitrary label so an A/B run cannot overwrite a baseline whose flags are
# otherwise identical -- an unforced re-run tags `_zstar` exactly like the 7-day validation does.
tag = (weno_order == 7 ? "" : "_weno$(weno_order)") * (use_zstar ? "_zstar" : "") *
      (use_atmos ? (jra_multiyear ? "_jra$(Dates.year(run_start))" : "_jraryf") : "") *
      (Δdeg == 0.05 ? "" : "_d$(replace(string(Δdeg), "." => ""))") *
      get(ENV, "REDSEA_TAG", "")
out(name) = joinpath(@__DIR__, "redsea_glorys_$(name)$(tag)")

# Oceananigans renamed this writer keyword in 0.113 (overwrite_existing -> overwrite_files), for
# Checkpointer as well as JLD2Writer. Picking it from the method's own kwarg list keeps ONE code path
# across both versions; splat it into every writer below rather than naming the keyword directly.
const OVERWRITE = (:overwrite_files in Base.kwarg_decl(which(JLD2Writer, Tuple{Any, Any}))) ?
                  (; overwrite_files = true) : (; overwrite_existing = true)

simulation.output_writers[:surface] = JLD2Writer(ocean.model, output_fields;
    schedule = output_schedule, OVERWRITE..., indices = (:, :, Nz), filename = out("surface"))
simulation.output_writers[:east_section] = JLD2Writer(ocean.model, output_fields;
    schedule = output_schedule, OVERWRITE..., indices = (Nx, :, :), filename = out("east_section"))
simulation.output_writers[:strait_section] = JLD2Writer(ocean.model, output_fields;
    schedule = output_schedule, OVERWRITE..., indices = (:, j_strait, :), filename = out("strait_section"))
simulation.output_writers[:free_surface] = JLD2Writer(ocean.model, (; η = free_surface_field(ocean.model.free_surface));
    schedule = output_schedule, OVERWRITE..., filename = out("free_surface"))

# A 28 h run should survive a crash. cleanup keeps only the most recent checkpoint (~300 MB).
simulation.output_writers[:checkpointer] = Checkpointer(ocean.model;
    schedule = TimeInterval(12hours), dir = @__DIR__, prefix = "redsea_checkpoint$(tag)",
    OVERWRITE..., cleanup = true)

run!(simulation)
