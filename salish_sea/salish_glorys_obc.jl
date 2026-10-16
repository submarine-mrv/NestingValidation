# salish_glorys_obc.jl
# GLORYS-driven open boundaries for the Salish Sea. NumericalEarth v0.6.1 / Oceananigans v0.110.17
# (the ss/open-boundary-conditions branch, now merged). Validated by a stable 7-day run at
# res_factor=1; see docs/obc-implementation-recipe.md.

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

arch = CPU()

# --- Run window -------------------------------------------------------------
# GLORYS boundary/IC data only needs to cover the simulated window (+ a small margin
run_start   = DateTime(1993, 1, 1)
run_days    = 7
margin_days = 2
boundary_dates = run_start : Day(1) : (run_start + Day(run_days + margin_days))

# --- Grid --------------------------------------------------------------------
# Domain retargeted on Juan de Fuca.
# φ₁ 47.0 -> 47.28: rows 1-30 of the old domain were 100% land in the model's own bathymetry
#   (first water at 47.305, the southern tip of Puget Sound) -- 11% of the domain doing nothing.
# φ₂ 49.75 -> 49.05: the old north edge was closed by bathymetry smoothing anyway (0 wet cells at
#   interpolation_passes=10). 49.05 cuts the Strait of Georgia as a single clean channel, 81 wet
#   cells and 139 m deep at 10 passes, ~60 km clear of the west boundary so the two open boundaries
#   no longer meet at a corner. Critically it stays above the 320 m depth below which GLORYS has NO
#   valid data anywhere in this region -- a 49.30 cut (367 m) would prescribe extrapolated values
#   across the deepest 13% of the section, in the deep-renewal inflow layer.
# Cost: the Fraser (mouth ~49.12) is now outside the domain and its freshwater arrives through the
#   northern boundary as GLORYS represents it. Failure signal to watch: JdF surface layer too salty
#   / estuarine exchange too weak -> revisit 49.30 with the OBC tapered below 320 m.
λ₁, λ₂  = (235.5, 238.0)  # domain in longitude (0-360° convention)
φ₁, φ₂  = (47.28, 49.05)  # domain in latitude

Nz = 36
depth = 400
z_scale = 150
z_reference = ExponentialDiscretization(Nz, -depth, 0; scale = z_scale)

# z* (free-surface-following) vertical coordinate. Wrapping the reference discretization in a
# MutableVerticalDiscretization is the ONLY change needed: ocean_simulation does not pass
# `vertical_coordinate`, so the model takes `default_vertical_coordinate(grid)`, which returns
# ZStarCoordinate() iff grid.z is mutable (and ZCoordinate() otherwise). Asking for ZStar on a
# static grid throws.
# A static grid here is a LINEAR free surface: thicknesses are fixed and η never changes them, so it
# does not fail when η is large -- the cost is an O(max|η| / Δz_top) error in surface tracer
# concentration and in the depth over which surface fluxes are spread. At Δz_top = 2.3 m with η of
# 0.2-0.5 m that ratio is 0.09-0.22, which is why static z is fine here. Tides (3-5 m) or a thinner
# surface cell would change that. See docs/salish_sea_notes.md. Off by default until validated
# against the stable static-z baseline.
#   SALISH_ZSTAR=1 julia --project=. salish_glorys_obc.jl
use_zstar = get(ENV, "SALISH_ZSTAR", "0") == "1"
z = use_zstar ? MutableVerticalDiscretization(z_reference) : z_reference
@info "Vertical coordinate: $(use_zstar ? "ZStar (mutable)" : "Z (static)")"

# Resolution knob. res_factor = 1 is the production grid: 0.01° cells = Δx ≈ 0.74 km (at 48.4°N)
# by Δy ≈ 1.11 km -- i.e. ~8× FINER than the 1/12° GLORYS (≈6.2 × 9.3 km here) that drives the
# boundaries and the nudging. That 8:1 refinement ratio is itself a suspect (see notes at the
# diagnostics block below), so keep this knob handy: res_factor = 3 gives ~2.2 × 3.3 km cells,
# a 3:1 ratio and ~9× fewer columns, which turns the 1-hour boundary test into a short run.
# Changing it regenerates the bathymetry cache once (a few minutes), then it's cached.
# Override without editing:  SALISH_RES_FACTOR=3 julia --project=. salish_glorys_obc.jl
res_factor = parse(Int, get(ENV, "SALISH_RES_FACTOR", "1"))
Δdeg = 0.01 * res_factor            # cell size in degrees; derive Nx, Ny so domain edits stay consistent
Nx = round(Int, (λ₂ - λ₁) / Δdeg)
Ny = round(Int, (φ₂ - φ₁) / Δdeg)
@info "Grid resolution: res_factor = $res_factor -> $Nx × $Ny × $Nz  (lon $(λ₁) to $(λ₂), lat $(φ₁) to $(φ₂))"

underlying_grid = LatitudeLongitudeGrid(arch;
                             size = (Nx, Ny, Nz),
                             latitude  = (φ₁, φ₂),
                             longitude = (λ₁, λ₂),
                             z,
                             halo = (7, 7, 7),
                             topology = (Bounded, Bounded, Bounded))

bottom_height = regrid_bathymetry(underlying_grid;
                                  minimum_depth = 10,
                                  interpolation_passes = 10)

# Bottom representation. GridFittedBottom makes every cell fully wet or fully dry, so a sloping
# basin becomes a staircase -- the remaining suspect for the deep-basin salinity growth, since the
# step slope Δz/Δx is ~0.038 at res_factor=1 against ~0.013 at res_factor=3, matching where the
# problem appears. PartialCellBottom instead shortens the bottommost cell to the true depth.
# ϵ = minimum fractional cell height guards against vanishingly thin (CFL-hostile) bottom cells.
#   SALISH_PARTIAL_CELLS=1 julia --project=. salish_glorys_obc.jl
use_partial_cells = get(ENV, "SALISH_PARTIAL_CELLS", "0") == "1"
bottom = use_partial_cells ? PartialCellBottom(bottom_height; minimum_fractional_cell_height = 0.2) :
                             GridFittedBottom(bottom_height)
@info "Bottom: $(use_partial_cells ? "PartialCellBottom(ϵ=0.2)" : "GridFittedBottom")"
grid = ImmersedBoundaryGrid(underlying_grid, bottom)

# GLORYS region: the model domain plus a small pad so exterior/boundary interpolation
# has data on all sides (GLORYS is 1/12° ≈ 0.083°).
# Deliberately NOT derived from the domain edges. The cached GLORYS filenames encode the region, so
# tying it to λ₁/φ₁ means every domain tweak invalidates the cache and forces a re-download (which
# needs Copernicus credentials). It only has to CONTAIN the domain plus a pad; this box covers the
# original 47.0–49.75 extent, so it still serves any domain inside it from the existing cache.
glorys_region = DW.BoundingBox(longitude = (235.25, 238.25), latitude = (46.75, 50.00))
@assert λ₁ >= 235.25 && λ₂ <= 238.25 && φ₁ >= 46.75 && φ₂ <= 50.00 "domain outside glorys_region"
glorys = DW.GLORYSDaily()

# --- Boundary conditions -----------------------------------------------------
# Names changed when ss/open-boundary-conditions merged (Oceananigans 0.110.8 branch -> 0.110.17):
#   Radiation                                     -> NormalRadiation
#   ImplicitGravityWaveRadiation(BoundaryCondition) -> SurfaceWaveRadiation(BoundaryCondition)
# GravityWaveRadiation (Flather) and NormalFlowBoundaryCondition are unchanged, the Flather value is
# still a 2-tuple (U, η), and SurfaceWaveRadiation's halo fill is byte-identical to the old Chapman.
using Oceananigans.BoundaryConditions: NormalRadiation,
                                       GravityWaveRadiationBoundaryCondition,
                                       SurfaceWaveRadiationBoundaryCondition
using Oceananigans.Operators: Δzᶠᶜᶜ, Δzᶜᶠᶜ, Axᶠᶜᶜ, Ayᶜᶠᶜ, Azᶜᶜᶜ
using Oceananigans.ImmersedBoundaries: immersed_inactive_node, immersed_peripheral_node

# Scheme selection, mirroring red_sea/. Defaults reproduce the validated run exactly:
# NormalRadiation everywhere, no reservoir, no tangential BC.
#
# This domain is the EXTREME case for the refinement hypothesis in
# docs/obc-implementation-recipe.md. Oblique radiation blew up at Canary (4.2x refinement vs the
# 1/12 deg parent) but ran clean at Red Sea (1.7x), both on boundaries terminated by land --
# evidence the failure tracks refinement rather than land. Salish runs 8.3x, and BOTH its open
# boundaries are narrow channels walled on both sides, so if refinement is the driver this is
# where it should fail hardest.
const HAS_OBLIQUE   = isdefined(Oceananigans.BoundaryConditions, :ObliqueRadiation)
const HAS_RESERVOIR = isdefined(Oceananigans.BoundaryConditions, :TracerReservoir)

const RAD_NORMAL = NormalRadiation(inflow_timescale = 6hours, outflow_timescale = Inf)

radiation_kind = get(ENV, "SALISH_RADIATION", "normal")
radiation_kind in ("normal", "oblique") || error("SALISH_RADIATION must be \"normal\" or \"oblique\"")
radiation_kind == "oblique" && !HAS_OBLIQUE &&
    error("SALISH_RADIATION=oblique needs Oceananigans >= 0.113; this environment has no ObliqueRadiation")
radiation = radiation_kind == "oblique" ?
    getfield(Oceananigans.BoundaryConditions, :ObliqueRadiation)(inflow_timescale = 6hours, outflow_timescale = Inf) :
    RAD_NORMAL

# T/S carry their own scheme. Oblique is never offered here: oblique on TRACERS at a
# land-terminated boundary is what drove salinity to 45.7 psu at Canary, and both boundaries
# here are channels walled on both sides.
const res_in  = parse(Float64, get(ENV, "SALISH_RESERVOIR_IN",  "0"))
const res_out = parse(Float64, get(ENV, "SALISH_RESERVOIR_OUT", "0"))
ts_kind = get(ENV, "SALISH_TS_SCHEME", "same")
ts_kind in ("same", "normal", "reservoir") || error("SALISH_TS_SCHEME must be same/normal/reservoir")
ts_kind == "reservoir" && !HAS_RESERVOIR &&
    error("SALISH_TS_SCHEME=reservoir needs Oceananigans >= 0.113.2")
ts_scheme = ts_kind == "normal"    ? RAD_NORMAL :
            ts_kind == "reservoir" ? getfield(Oceananigans.BoundaryConditions, :TracerReservoir)(
                                         inflow_length_scale = res_in, outflow_length_scale = res_out) :
            radiation

# SALISH_TANGENTIAL adds the ALONG-boundary velocity component, which this domain has never had:
# v was left to a default no-flux wall at the west boundary and u at the north one. In a channel
# the walls impose strong along-boundary structure, so this should matter more here than anywhere.
use_tangential = get(ENV, "SALISH_TANGENTIAL", "0") == "1"
@info "Open-boundary scheme: velocities $(radiation_kind), T/S $(ts_kind), tangential $(use_tangential ? "on" : "off")" *
      (ts_kind == "reservoir" ? @sprintf(" (reservoir in %.0f m, out %.0f m)", res_in, res_out) : "")

@inline wetcell(i, j, k, grid, ℓx, ℓy, ℓz) = !immersed_peripheral_node(i, j, k, grid, ℓx, ℓy, ℓz) & !immersed_inactive_node(i, j, k, grid, ℓx, ℓy, ℓz)

# All open boundary conditions (u/v normal flow, U/V Flather transport, η Chapman, T/S value)
# are now defined AFTER the restoring block below (STEP B): they reuse the GLORYS
# FieldTimeSeries built there for the nudging, so the exterior comes from GLORYS.

# --- STEP A: GLORYS nudging band (replaces the constant Relaxation sponge) ---
# DatasetRestoring adds a forcing  F = rate * mask * (GLORYS - model)  per field,
# with GLORYS linearly interpolated in space+time onto the model grid. We restrict it
# to a Gaussian band along the two open boundaries (west = min λ, north = max φ).
# WAS 1.0° per the Barents recipe (PR #5351), but that domain is far larger than this one: measured,
# a 1° band nudged the deep Strait of Georgia at 0.80 and the domain centre at 0.21 -- near-global
# restoring rather than a boundary sponge.
# Now 0.15° (15 cells) because the north boundary sits at 49.05, much closer to the priority region:
# at 0.3° it would restore Boundary Pass at mask 0.445 (a 27 h timescale), pinning the tidal-mixing
# passage that feeds Juan de Fuca to GLORYS. At 0.15° that drops to 0.039 (306 h), with Haro Strait
# and JdF at zero. NOTE this width is shared by the west and north bands in boundary_mask.
sponge_width      = 0.15       # degrees
restoring_timescale = 12hours  # TODO tune: strength of the boundary nudge
restoring_rate    = 1 / restoring_timescale

# Gaussian decay from each open boundary; band is the max over the two boundaries.
@inline gaussian_edge(d, w) = exp(-(d / w)^2)
@inline function boundary_mask(x, y, z, t)
    dw = x - λ₁     # distance inward from the west boundary (≥ 0 in-domain)
    dn = φ₂ - y     # distance inward from the north boundary (≥ 0 in-domain)
    return max(gaussian_edge(dw, sponge_width), gaussian_edge(dn, sponge_width))
end

# One DatasetRestoring per field. Pass `arch` (NOT `grid`): the FieldTimeSeries then stays
# on GLORYS's tiny native grid (fast native inpaint) and interpolates onto the model grid
# ON-THE-FLY during timestepping. Passing `grid` instead pre-interpolates every loaded slice
# onto the fine 250×275×36 immersed grid at setup -- very slow here, and looks like a hang.
# inpainting fills GLORYS land so the on-the-fly sample near the coast never hits NaN.
restoring_meta(var) = DW.Metadata(var; dataset = glorys, dates = boundary_dates, region = glorys_region)

# Inpainting fills GLORYS land cells (nearest-neighbor) so the regrid onto the model
# grid never samples NaN. The DEFAULT for T/S is NearestNeighborInpainting(Inf), which
# loops until zero NaNs remain across ALL levels -- but this coastal bbox is ~92% land and
# its deepest GLORYS levels have NO wet cells at all, so Inf never converges (infinite hang).
# A finite maxiter terminates; the loop still exits early on any level that fully fills, so
# results are identical to Inf on every fillable level (max fill distance on the 38×41 native
# grid is ≤ 79). Levels with NO wet cells at all cannot be filled and are ZEROED -- and those zeros
# ARE sampled at wet model points (see fill_down! below), so they must be repaired.
bc_inpainting = DW.NearestNeighborInpainting(100)

# Replace invalid (zero-filled) values with the deepest valid value above them, per column.
# GLORYS has no wet cells below ~300 m in the Strait of Georgia, so those levels come back as 0.0;
# the model bathymetry reaches -353 m, so interpolating onto the deepest wet model cell blended a
# valid level with a zeroed one and produced S = 30.16 sitting UNDER S = 33.92 -- a 3.8 psu density
# inversion, continuously re-imposed by the nudging. That was the instability seed.
# Only T and S: zero is a legitimate velocity, and a zero-filled deep u/v is harmless.
function fill_down!(a::AbstractArray{<:Any,3}, floor_value)
    nx, ny, nz = size(a)
    @inbounds for j in 1:ny, i in 1:nx
        last = NaN
        for k in nz:-1:1              # k = nz is the surface, k = 1 the bottom
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

# Floors are well below anything physical here (GLORYS S bottoms out ~19 in the Fraser plume,
# T ~4.4) and well above the 0.0 that marks an unfilled level.
function fill_down_fts!(fts, floor_value)
    for n in 1:length(fts.times)
        fill_down!(interior(fts[n]), floor_value)
    end
    return fts
end

# time_indices_in_memory = all: the default of 2 reloads slices from the cached inpainted files as
# the clock advances, which would silently undo the repair below. The whole series is ~25 MB.
n_times = length(boundary_dates)
# time_indexing: DatasetRestoring defaults to Cyclical(), and these series are not only the interior
# sponge -- their field_time_series is what the open-boundary exterior values are sampled from. Under
# NumericalEarth 0.8 GLORYSDaily is stamped at the CENTRE of its daily-mean window (12:00), so the
# first sample lands 12 h AFTER the run starts, and Cyclical fills that gap by wrapping to the FAR END
# of the record: the last day of the series would drive the first half-day of the run, at both the
# boundary and the sponge. Clamp holds the first sample instead (persistence), matching what the
# hand-rolled kernels in benchmarking/canary_template.jl do. A lead-in day is not available: GLORYS
# itself begins 1993-01-01, which is run_start.
restoring_kw = (; rate = restoring_rate, mask = boundary_mask, inpainting = bc_inpainting,
                  time_indices_in_memory = n_times,
                  time_indexing = Oceananigans.OutputReaders.Clamp())

T_restoring = DW.DatasetRestoring(restoring_meta(:temperature), arch; restoring_kw...)
S_restoring = DW.DatasetRestoring(restoring_meta(:salinity),    arch; restoring_kw...)
u_restoring = DW.DatasetRestoring(restoring_meta(:u_velocity),  arch; restoring_kw...)
v_restoring = DW.DatasetRestoring(restoring_meta(:v_velocity),  arch; restoring_kw...)

fill_down_fts!(T_restoring.field_time_series, 1.0)
fill_down_fts!(S_restoring.field_time_series, 5.0)
@info @sprintf("Filled down GLORYS deep levels: T ∈ [%.3f, %.3f], S ∈ [%.3f, %.3f] (native, t₀)",
               minimum(interior(T_restoring.field_time_series[1])), maximum(interior(T_restoring.field_time_series[1])),
               minimum(interior(S_restoring.field_time_series[1])), maximum(interior(S_restoring.field_time_series[1])))

# --- STEP B: GLORYS-driven T/S open-boundary exteriors -----------------------
# Replace Stage-A's constant exteriors (T=8/S=31 at every depth) with GLORYS values sampled at
# the boundary faces -- depth- AND time-varying. This removes the constant-vs-GLORYS conflict
# that drove the subsurface salinity overshoot (worst at depth, where GLORYS S≈34 but the
# constant BC said 31). We REUSE the native-grid FieldTimeSeries built above for the nudging
# (T_restoring/S_restoring): those live in model.forcing so Oceananigans advances them in time
# automatically, and the closures below read that same (always-current) GLORYS state.
# getbc calls the closure as (j,k,…) on west and (i,k,…) on north; DW.interpolate_to_grid maps
# native GLORYS data onto the model boundary node (i=1 west, j=Ny north) at Time(clock.time) --
# the exact machinery DatasetRestoring uses internally.
using Oceananigans.Units: Time

T_glorys_fts, T_glorys_native = T_restoring.field_time_series, T_restoring.native_grid
S_glorys_fts, S_glorys_native = S_restoring.field_time_series, S_restoring.native_grid

@inline T_west_exterior(j, k, grid, clock, fields)  = DW.interpolate_to_grid(T_glorys_fts, 1,       j, k, T_glorys_native, grid, Time(clock.time))
@inline T_north_exterior(i, k, grid, clock, fields) = DW.interpolate_to_grid(T_glorys_fts, i, grid.Ny, k, T_glorys_native, grid, Time(clock.time))
@inline S_west_exterior(j, k, grid, clock, fields)  = DW.interpolate_to_grid(S_glorys_fts, 1,       j, k, S_glorys_native, grid, Time(clock.time))
@inline S_north_exterior(i, k, grid, clock, fields) = DW.interpolate_to_grid(S_glorys_fts, i, grid.Ny, k, S_glorys_native, grid, Time(clock.time))

T_bcs = FieldBoundaryConditions(west  = ValueBoundaryCondition(T_west_exterior;  scheme = ts_scheme, discrete_form = true),
                                north = ValueBoundaryCondition(T_north_exterior; scheme = ts_scheme, discrete_form = true))
S_bcs = FieldBoundaryConditions(west  = ValueBoundaryCondition(S_west_exterior;  scheme = ts_scheme, discrete_form = true),
                                north = ValueBoundaryCondition(S_north_exterior; scheme = ts_scheme, discrete_form = true))

# --- STEP B: GLORYS-driven velocity + barotropic-transport open boundaries ---
# Stage A forced a steady 0.1 m/s inflow at BOTH open boundaries with no matching outflow -- a
# net mass source that piled water up (subsurface salinity overshoot + a max|w| jump after
# ~30 min, unfixed by the tighter wizard or the GLORYS T/S BCs). GLORYS velocities carry inflow
# AND outflow, balancing the transport. Reuse the u/v restoring FTS (auto-advanced via forcing).
u_glorys_fts, u_glorys_native = u_restoring.field_time_series, u_restoring.native_grid
v_glorys_fts, v_glorys_native = v_restoring.field_time_series, v_restoring.native_grid

# --- Net-transport correction ------------------------------------------------
# Measured: the raw depth-integrated GLORYS transport at the west boundary is a ~1.5e5 m³/s NET
# INFLOW with nothing balancing it (the north section carries ~100x less), which is comparable to
# the ENTIRE real one-way exchange through Juan de Fuca. GLORYS spans the strait in only 2-3 cells
# at 1/12°, so it cannot represent the estuarine exchange profile; depth-integrating that over the
# model's deeper, better-resolved column yields a smeared, mostly-inflow signal. Left alone it
# fills the basin at ~5 cm/hour.
#
# Fix (standard ROMS/NEMO regional practice): subtract a SPATIALLY UNIFORM normal-velocity
# correction Δu from the boundary exteriors so the total net flux equals a prescribed target,
#
#     Δu = (Q_raw - Q_target) / A          Q = net inflow (m³/s), A = total wet open-boundary area
#
# applied as u_west -= Δu and v_north += Δu (positive v at the north face is OUT of the domain), so
#     Q_corrected = Q_raw - Δu·A_west - Δu·A_north = Q_raw - Δu·A = Q_target.   ✓
#
# Applied to BOTH the 3-D normal-flow exteriors and the depth-integrated Flather transport, so the
# baroclinic and barotropic boundary data stay mutually consistent.
#
# Q_target = 0 is correct FOR THIS MODEL because it currently has no river forcing at all. In the
# real system the net outflow through the open boundaries is set by freshwater input: the Salish Sea
# receives O(5e3) m³/s in January (Fraser is at winter low, ~1e3 m³/s, but the coastal rain-fed
# rivers are high; the annual-mean total is ~1e4 m³/s). If river forcing is added later -- GloFAS is
# available in this NumericalEarth version, see src/DataWrangling/GloFAS -- then Q_target must become
# NEGATIVE (net OUTFLOW) and equal to that discharge, since rivers add the volume the boundaries must
# then export. Note the scale: even the annual-mean freshwater input is ~6% of the ~1.5e5 m³/s error
# being corrected here, so zero is the right target now and the sign subtlety only matters later.
net_transport_target = 0.0

# Δu is a GLOBAL quantity (it needs integrals over the whole boundary) but the BC closures are
# called per-(j,k), so it cannot be computed inside them. It is recomputed in a callback below and
# read from this Ref. `const` matters: a non-const global would make every closure call type-unstable.
const transport_correction = Ref(0.0)

# 3-D (baroclinic) boundary-normal velocity: Radiation with a GLORYS exterior.
@inline u_west_exterior(j, k, grid, clock, fields)  = DW.interpolate_to_grid(u_glorys_fts, 1,       j, k, u_glorys_native, grid, Time(clock.time)) - transport_correction[]
@inline v_north_exterior(i, k, grid, clock, fields) = DW.interpolate_to_grid(v_glorys_fts, i, grid.Ny, k, v_glorys_native, grid, Time(clock.time)) + transport_correction[]

# Tangential (along-boundary) components: v along the WEST boundary, u along the NORTH one.
# No transport correction on these -- that correction fixes the net volume flux THROUGH the
# boundary, which is carried entirely by the normal component.
@inline v_west_tangential(j, k, grid, clock, fields)  = DW.interpolate_to_grid(v_glorys_fts, 1,       j, k, v_glorys_native, grid, Time(clock.time))
@inline u_north_tangential(i, k, grid, clock, fields) = DW.interpolate_to_grid(u_glorys_fts, i, grid.Ny, k, u_glorys_native, grid, Time(clock.time))

nf(f) = NormalFlowBoundaryCondition(f; scheme = radiation, discrete_form = true)
vb(f) = ValueBoundaryCondition(f; scheme = radiation, discrete_form = true)

u_bcs = use_tangential ?
    FieldBoundaryConditions(west = nf(u_west_exterior), north = vb(u_north_tangential)) :
    FieldBoundaryConditions(west = nf(u_west_exterior))
v_bcs = use_tangential ?
    FieldBoundaryConditions(north = nf(v_north_exterior), west = vb(v_west_tangential)) :
    FieldBoundaryConditions(north = nf(v_north_exterior))

# GLORYS free surface (zos) as the Flather exterior η. This is the ONLY place an exterior η
# enters the model -- see the note on the Chapman BC below, which despite appearances never
# reads one. The Flather halo fill (gravity_wave_schemes.jl in the pinned Oceananigans) is:
#
#     west:   U = Uᵉˣᵗ - √(gH) (ηᵇ - ηᵉˣᵗ)        north:  V = Vᵉˣᵗ + √(gH) (ηᵇ - ηᵉˣᵗ)
#
# with ηᵇ the face-average of the adjacent interior η. Leaving ηᵉˣᵗ = 0 while the model's own η
# sits at ~0.2 m therefore injected a spurious ±√(gH)·0.2 ≈ ±8.8 m²/s per unit width of
# transport. GLORYS zos in this region is 0.12-0.23 m, i.e. very close to the model's η, so
# using it removes that term rather than merely rescaling it.
#
# The series is 2-D ((Center, Center, Nothing) for :free_surface) and tiny, so we hold ALL time
# indices in memory. That matters: unlike the T/S/u/v series this one is NOT reachable from
# model.forcing, so Oceananigans' update_model_field_time_series! cannot advance its window --
# a closure is opaque to that traversal. With the whole series resident, time interpolation
# needs no advancing at all.
zos_meta = DW.Metadata(:free_surface; dataset = glorys, dates = boundary_dates, region = glorys_region)
zos_fts  = FieldTimeSeries(zos_meta, arch;
                           time_indices_in_memory = length(boundary_dates),
                           # Clamp, NOT Cyclical. NumericalEarth 0.8 stamps GLORYSDaily at the CENTRE of
                           # its daily-mean window (12:00), so the first sample lands 12 h AFTER the run
                           # starts. Cyclical() fills that gap by wrapping to the FAR END of the record --
                           # the last day of the series forcing the first half-day of the run. Clamp holds
                           # the first sample instead, i.e. persistence, which is what the hand-rolled
                           # kernels in benchmarking/canary_template.jl already do. A lead-in day is not an
                           # option: GLORYS itself begins 1993-01-01, the usual run_start.
                           time_indexing = Oceananigans.OutputReaders.Clamp(),
                           inpainting = bc_inpainting)
zos_native = zos_fts.grid

# Sampled at the boundary-adjacent interior node (west i=1, north j=Ny), matching how the T/S and
# u/v exteriors above sample. GLORYS is ~8x coarser than the model here, so the offset between
# that node and the boundary face is far inside the interpolation error.
@inline zos_at(i, j, grid, t) = DW.interpolate_to_grid(zos_fts, i, j, 1, zos_native, grid, t)

# Barotropic transport (Flather): depth-integrate GLORYS u/v over the wet column, paired with the
# GLORYS zos exterior η above. Returns (exterior transport, exterior η).
@inline function west_U_obc(j, k, grid, clock, fields)
    i = 1
    t = Time(clock.time)
    U = zero(grid)
    @inbounds for kk in 1:grid.Nz
        wet = wetcell(i, j, kk, grid, Face(), Center(), Center())
        uk  = DW.interpolate_to_grid(u_glorys_fts, i, j, kk, u_glorys_native, grid, t) - transport_correction[]
        U  += ifelse(wet, uk * Δzᶠᶜᶜ(i, j, kk, grid), zero(grid))
    end
    return (U, zos_at(i, j, grid, t))
end

@inline function north_V_obc(i, k, grid, clock, fields)
    j = grid.Ny + 1
    t = Time(clock.time)
    V = zero(grid)
    @inbounds for kk in 1:grid.Nz
        wet = wetcell(i, j, kk, grid, Center(), Face(), Center())
        vk  = DW.interpolate_to_grid(v_glorys_fts, i, j, kk, v_glorys_native, grid, t) + transport_correction[]
        V  += ifelse(wet, vk * Δzᶜᶠᶜ(i, j, kk, grid), zero(grid))
    end
    return (V, zos_at(i, grid.Ny, grid, t))
end

U_bcs = FieldBoundaryConditions(grid, (Face(), Center(), nothing);
                                west  = GravityWaveRadiationBoundaryCondition(west_U_obc; discrete_form = true))
V_bcs = FieldBoundaryConditions(grid, (Center(), Face(), nothing);
                                north = GravityWaveRadiationBoundaryCondition(north_V_obc; discrete_form = true))

# Free surface (Chapman): deliberately left as-is, and it does NOT take a GLORYS exterior.
# Checked the halo fill in the pinned Oceananigans (gravity_wave_schemes.jl, _fill_west_halo! for
# IGWVBC): it computes  η[0,j,k] = (ηᵇ + C η₁) / (1 + C),  with C = √(gH) Δt / Δx -- pure implicit
# radiation of the interior η outward. It never calls getbc, so the `0` inside
# SurfaceWaveRadiationBoundaryCondition() is not an exterior value that is being applied;
# it is simply ignored. Passing GLORYS zos here would change nothing. The Flather 2-tuple above is
# the only route by which an exterior η reaches the solution.
η_bcs = FieldBoundaryConditions(grid, (Center(), Center(), Face());
                                west  = SurfaceWaveRadiationBoundaryCondition(),
                                north = SurfaceWaveRadiationBoundaryCondition())

# --- Ocean simulation -------------------------------------------------------
# ocean_simulation merges these forcings with its own barotropic/radiative forcing per field.
#
# Tracer advection order is a knob because it is the leading suspect for the residual deep-basin
# salinity growth: WENO is not strictly monotone, and with the nudging and zero-fill seeds now
# eliminated, a rising maximum can only come from dispersive overshoot or from topography error.
# Order 5 is more dissipative than the ocean_simulation default of 7. Everything else about the
# scheme matches that default, so the comparison isolates the order.
#   SALISH_WENO_ORDER=5 julia --project=. salish_glorys_obc.jl
weno_order = parse(Int, get(ENV, "SALISH_WENO_ORDER", "7"))
tracer_advection = WENO(order = weno_order,
                        time_discretization = AdaptiveVerticallyImplicitDiscretization(cfl = 0.5))
@info "Tracer advection: WENO(order = $weno_order)"

ocean = ocean_simulation(grid;
                        tracers = (:T, :S),
                        tracer_advection,
                        forcing = (T = T_restoring, S = S_restoring,
                                   u = u_restoring, v = v_restoring),
                        boundary_conditions = (; u = u_bcs, U = U_bcs, v = v_bcs, V = V_bcs, η = η_bcs, T = T_bcs, S = S_bcs),
)

# --- STEP A: GLORYS initial conditions (stratified) --------------------------
# Sampled from the SAME fill_down!-repaired FieldTimeSeries the nudging and the boundary exteriors
# use, so the IC and the forcing agree by construction.
#
# The previous `set!(field, ::Metadatum)` route re-read the cached *_inpainted.jld2 from disk and so
# never saw the repair, and an absolute floor could not clean it afterwards: interpolation blends a
# valid 33.9 with a zeroed deep level to give ANY value in between (observed 30.16, and 15.24 once
# PartialCellBottom activated deeper cells -- below the 19.098 minimum that exists anywhere in
# GLORYS, i.e. 14 psu of fresh water sitting under salty water). Sampling the repaired series has no
# such failure mode: the result is bounded by the repaired data.
#
# Written into scratch fields, then copied: assigning to the model tracers directly triggers a halo
# fill through the GLORYS Value+Radiation OBCs, whose discrete closures need a clock that a plain
# set! cannot supply. Model-tracer halos fill at run! start, when the clock exists.
function ic_from_fts!(field, fts, native_grid, t)
    a = interior(field)
    @inbounds for k in 1:Nz, j in 1:Ny, i in 1:Nx
        a[i, j, k] = DW.interpolate_to_grid(fts, i, j, k, native_grid, grid, t)
    end
    return field
end

T_ic = CenterField(grid); ic_from_fts!(T_ic, T_glorys_fts, T_glorys_native, Time(0.0))
S_ic = CenterField(grid); ic_from_fts!(S_ic, S_glorys_fts, S_glorys_native, Time(0.0))
interior(ocean.model.tracers.T) .= interior(T_ic)
interior(ocean.model.tracers.S) .= interior(S_ic)

# The free-surface field is `displacement` on SplitExplicitFreeSurface (what ocean_simulation
# selects here); older/other free surfaces call it `η`. Accept either.
free_surface_field(fs) = hasproperty(fs, :η) ? fs.η : fs.displacement

# --- Free-surface initial condition from GLORYS zos --------------------------
# Previously η started at REST (0) while GLORYS zos at the open boundary is +0.211 m. That
# manufactured a basin-filling spin-up: the Flather term -√(gH)(ηᵇ - ηᵉˣᵗ) correctly drove inflow
# to close a ~0.2 m deficit, and at the observed ~4 cm/hour that adjustment takes ~5 HOURS --
# far longer than the 1-hour test, so every bounded run so far has been measuring the transient
# rather than the boundary balance. Starting from GLORYS zos removes it.
#
# Sampled with the same zos_at / interpolate_to_grid path the Flather BC uses, so the initial
# interior and the boundary exterior are consistent by construction rather than by coincidence.
# Land columns are zeroed. Assigning to `interior` deliberately avoids a halo fill: the η BC's
# halo fill needs a clock (anchored_fill/stage_Δt) that a plain set! cannot supply -- the same
# trap that broke the tracer ICs once the clock-dependent OBCs went in.
η_field = free_surface_field(ocean.model.free_surface)
let ηi = interior(η_field), t₀ = Time(0.0)
    for i in 1:Nx, j in 1:Ny
        wet = !immersed_inactive_node(i, j, Nz, grid, Center(), Center(), Center())
        ηi[i, j, 1] = ifelse(wet, zos_at(i, j, grid, t₀), 0.0)
    end
end
@info @sprintf("Initialized η from GLORYS zos: η ∈ [%.4f, %.4f] m", minimum(interior(η_field)), maximum(interior(η_field)))

# Velocities START FROM REST (the conventional choice: T/S from the parent, momentum spun up by the
# model). Seeding u/v from GLORYS was tested at res_factor=3 and deliberately turned back off --
# it did NOT improve the volume budget, it slightly worsened it:
#
#                            u,v from REST      u,v from GLORYS
#   net at t = 0             0                  +5.50e4 m³/s   (raw GLORYS convergence in the IC)
#   net, mean last 10        +8.5e2 m³/s        +3.7e3 m³/s
#   drift                    -4.0e-4 m/hr       -9.3e-4 m/hr
#   max|η|                   0.343 m            0.503 m        (looks like a seiche off the t=0 shock)
#
# The only argument for seeding was spin-up speed: from rest the model reaches ~0.15 m/s in an hour
# against a real Juan de Fuca exchange of 0.3-0.4 m/s, so the estuarine circulation takes several
# hours to develop. That is a transient to wait out, not a reason to inject an inconsistent field.
#
# Set to `true` to seed from GLORYS. NOTE if you do: it is seeded with RAW GLORYS, without the Δu
# net-transport correction -- Δu is a boundary construct (a uniform normal-velocity adjustment to the
# BC exteriors) and smearing it across the interior would not be meaningful. Cost is Nx*Ny*Nz*2
# interpolations: trivial at res_factor=3, ~5M at res_factor=1.
initialize_velocities_from_glorys = false
if initialize_velocities_from_glorys
    let ui = interior(ocean.model.velocities.u), vi = interior(ocean.model.velocities.v), t₀ = Time(0.0)
        for i in 1:Nx, j in 1:Ny, k in 1:Nz
            uwet = wetcell(i, j, k, grid, Face(), Center(), Center())
            vwet = wetcell(i, j, k, grid, Center(), Face(), Center())
            ui[i, j, k] = ifelse(uwet, DW.interpolate_to_grid(u_glorys_fts, i, j, k, u_glorys_native, grid, t₀), 0.0)
            vi[i, j, k] = ifelse(vwet, DW.interpolate_to_grid(v_glorys_fts, i, j, k, v_glorys_native, grid, t₀), 0.0)
        end
    end
    @info "Initialized u, v from GLORYS"
end

# --- Progress monitoring + blow-up localization ------------------------------
# The last three runs all went unstable at t ≈ 29-35 min, and the previous progress line could
# not say WHERE: it reported whole-array extrema (land included) and the JLD2 output was
# surface-only, while the actual excursion is SUBSURFACE (3-D log hit S = 36.28 / max|w| = 5.5e-2
# while the saved surface stayed flat at S = 31.82 / 9e-3). So this block reports, per callback:
#
#   (1) WHERE the extrema are -- argmax (i,j,k) + depth + whether that cell sits on an open
#       boundary. Distinguishes a boundary-formulation problem (west i=1 / north j=Ny) from a
#       dispersive advection overshoot in the interior, which is the fork the next fix hangs on.
#   (2) The NET VOLUME FLUX through each open section, and the domain-mean sea-level rise it
#       implies. Tests directly whether the GLORYS boundary velocities -- interpolated onto a
#       bathymetry ~8× better resolved than GLORYS's own -- leave a mass imbalance.
#   (3) max|η| and the west-boundary η vs GLORYS zos mismatch, which drives the Flather term.
#
# Extrema are WET-ONLY: land cells hold zeros that corrupted the old minima (the "S ∈ [0, ...]"
# and "0.000" artifacts). Velocity fields are already masked to zero on land by Oceananigans, so
# only the tracers need the mask.
z_centers = znodes(grid, Center(), Center(), Center())

# Vertical spacing, top-down -- confirms the z-grid is spending its levels in the water column.
@info @sprintf("Δz near surface: %s m (top of %d levels, domain %d m)",
               join((@sprintf("%.1f", d) for d in reverse(diff(z_centers))[1:min(6, Nz-1)]), ", "),
               Nz, depth)

wet3d = [!immersed_inactive_node(i, j, k, grid, Center(), Center(), Center())
         for i in 1:Nx, j in 1:Ny, k in 1:Nz]
wet_idx = findall(wet3d)
@info "Diagnostics: $(length(wet_idx)) wet cells of $(Nx*Ny*Nz) ($(round(100*length(wet_idx)/(Nx*Ny*Nz), digits=1))%)"

# Total wet surface area, for converting a net volume imbalance into a sea-level tendency.
wet_surface_area = sum(Azᶜᶜᶜ(i, j, Nz, grid) for i in 1:Nx, j in 1:Ny
                       if !immersed_inactive_node(i, j, Nz, grid, Center(), Center(), Center()); init = 0.0)

# Wet extremum + its location. Returns (min, argmin, max, argmax).
function wet_extrema(field)
    a = interior(field)[wet_idx]
    nmin, nmax = argmin(a), argmax(a)
    return a[nmin], wet_idx[nmin], a[nmax], wet_idx[nmax]
end

# Label a cell by its position relative to the two open boundaries (west = i 1, north = j Ny).
function whereis(I)
    i, j, k = Tuple(I)
    z = z_centers[clamp(k, 1, Nz)]   # w is on Faces (k can be Nz+1); u/v extend one index in i/j
    tag = i == 1        ? "WEST-bnd"   :
          j == Ny       ? "NORTH-bnd"  :
          i <= 3        ? "near-west"  :
          j >= Ny - 2   ? "near-north" : "interior"
    return @sprintf("(%3d,%3d,%2d) z=%6.1fm %s", i, j, k, z, tag)
end

# Net volume flux through each open section (m³/s). Sign convention: POSITIVE = into the domain,
# so west inflow is +u and north inflow is -v (the north face normal points out of the domain).
function section_transports(model)
    u, v = model.velocities.u, model.velocities.v
    Φw = 0.0
    for k in 1:Nz, j in 1:Ny
        wetcell(1, j, k, grid, Face(), Center(), Center()) || continue
        @inbounds Φw += u[1, j, k] * Axᶠᶜᶜ(1, j, k, grid)
    end
    Φn = 0.0
    for k in 1:Nz, i in 1:Nx
        wetcell(i, Ny + 1, k, grid, Center(), Face(), Center()) || continue
        @inbounds Φn -= v[i, Ny + 1, k] * Ayᶜᶠᶜ(i, Ny + 1, k, grid)
    end
    return Φw, Φn
end

# EXCHANGE flow through the west (Juan de Fuca) section: inflow and outflow accumulated separately,
# rather than the net. This is the number the literature reports and the one that governs JdF water
# properties -- roughly 1e5 m³/s each way, i.e. ~15x the net. The net can be right while the exchange
# is wrong by a factor of several, and the net-only diagnostic would not show it.
# Also returns the flux-weighted salinity of each limb, which with Qin/Qout is the standard
# Knudsen-style check: Qin·Sin ≈ Qout·Sout, and Sin − Sout measures how much mixing the model does.
function west_exchange(model)
    u, S = model.velocities.u, model.tracers.S
    Qin = Qout = 0.0
    SQin = SQout = 0.0
    for k in 1:Nz, j in 1:Ny
        wetcell(1, j, k, grid, Face(), Center(), Center()) || continue
        @inbounds q = u[1, j, k] * Axᶠᶜᶜ(1, j, k, grid)
        @inbounds s = S[1, j, k]
        if q > 0
            Qin  += q; SQin  += q * s
        else
            Qout += q; SQout += q * s
        end
    end
    Sin  = Qin  > 0 ? SQin / Qin   : 0.0
    Sout = Qout < 0 ? SQout / Qout : 0.0
    return Qin, Qout, Sin, Sout
end

# Recompute the uniform normal-velocity correction Δu = (Q_raw - Q_target) / A from the RAW GLORYS
# exteriors (no correction applied inside -- that is what makes it idempotent rather than a feedback
# loop). Depends only on the GLORYS data and the clock, never on model state.
function update_transport_correction!(sim)
    t = Time(sim.model.clock.time)
    Q = 0.0   # net inflow implied by the raw GLORYS boundary data, m³/s
    A = 0.0   # total wet open-boundary cross-sectional area, m²
    for j in 1:Ny, k in 1:Nz
        wetcell(1, j, k, grid, Face(), Center(), Center()) || continue
        a  = Axᶠᶜᶜ(1, j, k, grid)
        Q += DW.interpolate_to_grid(u_glorys_fts, 1, j, k, u_glorys_native, grid, t) * a
        A += a
    end
    for i in 1:Nx, k in 1:Nz
        wetcell(i, Ny + 1, k, grid, Center(), Face(), Center()) || continue
        a  = Ayᶜᶠᶜ(i, Ny + 1, k, grid)
        Q -= DW.interpolate_to_grid(v_glorys_fts, i, Ny + 1, k, v_glorys_native, grid, t) * a
        A += a
    end
    transport_correction[] = A > 0 ? (Q - net_transport_target) / A : 0.0
    return nothing
end

# Every 10 iterations is ample: GLORYS is daily, so Δu drifts on a ~day timescale while 10 iterations
# is at most a few minutes. Updating every step would re-run ~19k interpolations per iteration at
# res_factor = 1, comparable to a boundary halo fill, for no accuracy gain.
update_transport_correction!(ocean)
add_callback!(ocean, update_transport_correction!, IterationInterval(10))
@info @sprintf("Initial net-transport correction: Δu = %+.4e m/s (target net = %.1f m³/s)",
               transport_correction[], net_transport_target)

# West-boundary mean interior η vs mean GLORYS zos there. This is the mismatch the Flather η_ext
# fix targets, so it should now be small; a persistent gap means the two have different reference
# levels (GLORYS zos is height above the geoid, the model's η is relative to z = 0) and that an
# offset correction is needed rather than the raw field.
function west_boundary_eta(model)
    ηi = interior(free_surface_field(model.free_surface))
    t  = Time(model.clock.time)
    n, ηs, zs = 0, 0.0, 0.0
    for j in 1:Ny
        wetcell(1, j, Nz, grid, Center(), Center(), Center()) || continue
        n += 1
        ηs += @inbounds ηi[1, j, 1]
        zs += zos_at(1, j, grid, t)
    end
    return n == 0 ? (0.0, 0.0) : (ηs / n, zs / n)
end

# Reference range for the excursion diagnostic. Bound against the FORCING DATA, not the initial
# condition: the boundary and the nudging can legitimately bring in water saltier than anything the
# regridded IC contained (native max 33.952 vs IC max 33.920 here), so referencing the IC reports
# ~0.03 of physical adjustment as if it were spurious. Only an excursion beyond the native range
# is unphysical. The native FTS is post-fill_down, so its extrema are the true available bounds.
native_range(fts) = (minimum(interior(fts[1])), maximum(interior(fts[1])))
T_ref = native_range(T_restoring.field_time_series)
S_ref = native_range(S_restoring.field_time_series)
T_ic_range = let (lo, _, hi, _) = wet_extrema(ocean.model.tracers.T); (lo, hi) end
S_ic_range = let (lo, _, hi, _) = wet_extrema(ocean.model.tracers.S); (lo, hi) end
@info @sprintf("Reference (GLORYS native): T ∈ [%.3f, %.3f] | S ∈ [%.3f, %.3f]", T_ref[1], T_ref[2], S_ref[1], S_ref[2])
@info @sprintf("Initial condition (model grid, wet): T ∈ [%.3f, %.3f] | S ∈ [%.3f, %.3f]",
               T_ic_range[1], T_ic_range[2], S_ic_range[1], S_ic_range[2])

# Wall-clock pacing, so a long run reports its own ETA instead of needing to be extrapolated.
# The ETA uses the rate over the LAST callback interval, not the cumulative rate: Δt ramps from 2 s
# to ~280 s over the first ~520 iterations, so a cumulative rate is dominated by the slow early
# steps and would give a wildly pessimistic ETA that keeps falling for hours.
const wall_t0    = Ref(Base.time())
const wall_last  = Ref(Base.time())
const sim_last   = Ref(0.0)

function progress(sim)
    model = sim.model
    u, v, w = model.velocities
    T, S = model.tracers.T, model.tracers.S
    η = free_surface_field(model.free_surface)

    Tmin, _, Tmax, _       = wet_extrema(T)
    Smin, Imin, Smax, Imax = wet_extrema(S)
    _, Iw = findmax(abs, interior(w))
    _, Iu = findmax(abs, interior(u))

    Φw, Φn = section_transports(model)
    net    = Φw + Φn                          # m³/s into the domain
    dηdt   = net / wet_surface_area * 3600    # implied domain-mean sea-level rise, m/hour

    @info @sprintf("iter %d | t = %s | Δt = %s | max|u,v,w| = (%.2e, %.2e, %.2e) | max|η| = %.3e | T ∈ [%.3f, %.3f] | S ∈ [%.3f, %.3f] (excursion %+.3f / %+.3f)",
                   iteration(sim), prettytime(time(sim)), prettytime(sim.Δt),
                   maximum(abs, u), maximum(abs, v), maximum(abs, w), maximum(abs, η),
                   Tmin, Tmax, Smin, Smax, Smin - S_ref[1], Smax - S_ref[2])
    @info @sprintf("        where: Smax %s | Smin %s | max|w| %s | max|u| %s",
                   whereis(Imax), whereis(Imin), whereis(Iw), whereis(Iu))
    η̄, z̄ = west_boundary_eta(model)
    @info @sprintf("        transport: west %+.4e  north %+.4e  net %+.4e m³/s  ->  %+.3e m/hour mean sea level | west-bnd η̄ = %+.4f m, GLORYS zos = %+.4f m, mismatch %+.4f m",
                   Φw, Φn, net, dηdt, η̄, z̄, η̄ - z̄)

    Qin, Qout, Sin, Sout = west_exchange(model)
    salt_in = Qin * Sin
    knudsen = salt_in + Qout * Sout
    @info @sprintf("        JdF exchange: Qin %+.3e  Qout %+.3e m³/s  (lit. TEF ~1.4e5 at the mouth) | Sin %.3f  Sout %.3f  ΔS %.3f | Knudsen residual %+.1f%% of Qin·Sin",
                   Qin, Qout, Sin, Sout, Sin - Sout, 100 * knudsen / max(abs(salt_in), eps()))

    now  = Base.time()
    wall = now - wall_t0[]
    Δwall, Δsim = now - wall_last[], time(sim) - sim_last[]
    wall_last[], sim_last[] = now, time(sim)
    rate = Δwall > 0 ? Δsim / Δwall : 0.0        # simulated seconds per wall second, recent window
    eta  = (isfinite(sim.stop_time) && rate > 0) ? (sim.stop_time - time(sim)) / rate : NaN
    @info @sprintf("        pacing: wall %s | %.0f× realtime (recent) | %.2f s/iter | ETA %s",
                   prettytime(wall), rate, Δwall / 10, isnan(eta) ? "n/a" : prettytime(eta))
    flush(stderr)
    return nothing
end
add_callback!(ocean, progress, IterationInterval(10))

# --- Spin-up: short fixed-Δt run to check the GLORYS ICs + nudging are stable ----
ocean.Δt = 2
ocean.stop_iteration = 100
run!(ocean)

# --- Post-spin-up snapshot: confirm the GLORYS regrid is spatially sensible ---
# The spun-up state is still ≈ the GLORYS ICs, so this checks the regrid landed in the
# right place (Fraser plume in Georgia Strait, oceanic inflow in Juan de Fuca) -- something
# a min/max can't catch. Land (immersed) columns are masked to NaN so the coastline shows.
using Oceananigans.Fields: location
# 1-D lon/lat in degrees. NB: xnodes/ynodes on a LatLonGrid return metres; use λnodes/φnodes.
snap_λ(f) = λnodes(f.grid, map(L -> L(), location(f))...)
snap_φ(f) = φnodes(f.grid, map(L -> L(), location(f))...)
# surface (top level) wet mask: NaN over land, 1 over ocean
wet2d = [immersed_inactive_node(i, j, Nz, grid, Center(), Center(), Center()) ? NaN : 1.0
         for i in 1:Nx, j in 1:Ny]
surftop(f) = Array(interior(f))[1:Nx, 1:Ny, Nz]   # crop Face (u/v) fields to Nx×Ny centres

try
    Tm  = surftop(ocean.model.tracers.T) .* wet2d
    Sm  = surftop(ocean.model.tracers.S) .* wet2d
    spd = sqrt.(surftop(ocean.model.velocities.u).^2 .+ surftop(ocean.model.velocities.v).^2) .* wet2d
    λ, φ = snap_λ(ocean.model.tracers.T), snap_φ(ocean.model.tracers.T)

    fig = Figure(size = (1500, 520))
    for (n, (data, label, cmap)) in enumerate(((Tm,  "surface T (°C)",     :thermal),
                                               (Sm,  "surface S (psu)",    :haline),
                                               (spd, "surface speed (m/s)", :speed)))
        ax = Axis(fig[1, n]; title = label, xlabel = "λ", ylabel = "φ", aspect = DataAspect())
        hm = heatmap!(ax, λ, φ, data; nan_color = :gray70, colormap = cmap)
        Colorbar(fig[2, n], hm; vertical = false, flipaxis = false)
    end
    out = "salish_glorys_ic_snapshot.png"
    save(out, fig)
    @info "Saved model snapshot -> $out"
catch err
    @error "Snapshot plotting failed (the run itself was fine)" exception = (err, catch_backtrace())
end

# --- Adaptive-timestep run: wizard + JLD2 surface output ---------------------
# Continues from the spun-up state. The wizard ramps Δt up from 2s under a CFL limit.
#
# CURRENTLY: bounded 1-hour VALIDATION test (output every 5 min, run! active) -- confirms the
# nudging + adaptive stepping stay stable past spin-up and that JLD2 output writes.
# For the real 7-day run, swap the two TEST lines to the PRODUCTION values noted beside them
# and launch in your own terminal (hours of wall-clock; keep a REPL open to avoid recompiling).
ocean_timescale(model) = cell_advection_timescale(model)
# Conservative wizard, kept from the instability hunt. CFL caps Δt near 280 s here (Δx ≈ 740 m,
# max|u| ~0.4 m/s), just under the 5-minute ceiling, after a ~520-iteration ramp from 2 s that
# itself covers ~7.8 h of simulation. A 7-day run is therefore ~2580 iterations.
wizard = TimeStepWizard(; cfl = 0.15, max_Δt = 5minutes, max_change = 1.1, cell_advection_timescale = ocean_timescale)
ocean.callbacks[:wizard] = Callback(wizard, IterationInterval(10))
ocean.stop_iteration = Inf
# Override for a short test:  SALISH_STOP_DAYS=0.04 julia --project=. salish_glorys_obc.jl
stop_days = parse(Float64, get(ENV, "SALISH_STOP_DAYS", string(run_days)))
ocean.stop_time = stop_days * days
@info @sprintf("Production run: stop_time = %.3g days (~%.0f iterations, est. %.1f h wall at ~8.7 s/iter)",
               stop_days, 520 + max(0, stop_days*86400 - 28200)/280, (520 + max(0, stop_days*86400 - 28200)/280)*8.7/3600)

output_fields   = merge(ocean.model.velocities, ocean.model.tracers)
# 3-hourly: 56 frames over 7 days, ~250 MB across the five writers. Finer buys nothing -- the
# forcing is daily GLORYS and there are no tides, so there is no sub-3-hourly signal to resolve.
output_schedule = TimeInterval(3hours)
# Tag non-default runs so a comparison run does not overwrite the baseline outputs.
run_tag = (weno_order == 7 ? "" : "_weno$(weno_order)") * (use_partial_cells ? "_pcb" : "") *
          (use_zstar ? "_zstar" : "")

# Oceananigans renamed this writer keyword in 0.113 (overwrite_existing -> overwrite_files), for
# Checkpointer as well as JLD2Writer. Picking it from the method's own kwarg list keeps ONE code path
# across both versions; splat it into every writer below rather than naming the keyword directly.
const OVERWRITE = (:overwrite_files in Base.kwarg_decl(which(JLD2Writer, Tuple{Any, Any}))) ?
                  (; overwrite_files = true) : (; overwrite_existing = true)

# Surface slice, as before. NOTE: this alone hid the blow-up on the last two runs -- the excursion
# is subsurface, and the surface stayed flat for the whole hour. Keep it, but do not judge from it.
ocean.output_writers[:surface] = JLD2Writer(ocean.model, output_fields;
                                            schedule = output_schedule,
                                            OVERWRITE...,
                                            indices = (:, :, Nz),
                                            filename = "salishsea_glorys_openboundaries_test1$(run_tag)")

# Open-boundary sections: the FULL depth structure at each open boundary. These are the cheap,
# high-value ones -- 275×36 and 250×36 cells, ~1/250th of a full 3-D dump -- and they show
# directly whether the deep inflow at the west (Juan de Fuca, where GLORYS S ≈ 34) and the north
# (Georgia Strait) is where the salinity overshoot and the w spike are born.
ocean.output_writers[:west_section] = JLD2Writer(ocean.model, output_fields;
                                                 schedule = output_schedule,
                                                 OVERWRITE...,
                                                 indices = (1, :, :),
                                                 filename = "salishsea_glorys_west_section$(run_tag)")

ocean.output_writers[:north_section] = JLD2Writer(ocean.model, output_fields;
                                                  schedule = output_schedule,
                                                  OVERWRITE...,
                                                  indices = (:, Ny, :),
                                                  filename = "salishsea_glorys_north_section$(run_tag)")

# Mid-depth horizontal slice (~150 m): catches an excursion that is subsurface but NOT on a
# boundary, i.e. the dispersive-advection branch of the diagnosis.
k_mid = argmin(abs.(z_centers .+ 150))
@info @sprintf("Subsurface output slice at k = %d, z = %.1f m", k_mid, z_centers[k_mid])
ocean.output_writers[:subsurface] = JLD2Writer(ocean.model, output_fields;
                                               schedule = output_schedule,
                                               OVERWRITE...,
                                               indices = (:, :, k_mid),
                                               filename = "salishsea_glorys_subsurface$(run_tag)")

# Free surface, in its own writer: η is a single-level field, so it cannot share the `indices`
# of the slice writers above. η_ext is still hardcoded to zero in the Flather and Chapman BCs,
# so watching η grow (and where) is the direct test of that mismatch.
ocean.output_writers[:free_surface] = JLD2Writer(ocean.model, (; η = free_surface_field(ocean.model.free_surface));
                                                 schedule = output_schedule,
                                                 OVERWRITE...,
                                                 filename = "salishsea_glorys_free_surface$(run_tag)")

run!(ocean)   # bounded 1-hour test; for production keep this and use the PRODUCTION values above
