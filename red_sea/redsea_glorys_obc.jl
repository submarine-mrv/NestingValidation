# redsea_glorys_obc.jl
# GLORYS-driven open boundary for the Red Sea + Gulf of Aden.
# Adapted from salish_sea/salish_glorys_obc.jl; see docs/obc-implementation-recipe.md for the
# reusable pattern and docs/salish_sea_notes.md for why each piece is here.
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

arch = CPU()

# --- Run window -------------------------------------------------------------
run_start   = DateTime(1993, 1, 1)
run_days    = 7
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
using Oceananigans.Operators: Δzᶜᶜᶜ, Δzᶠᶜᶜ, Axᶠᶜᶜ, Ayᶜᶠᶜ, Azᶜᶜᶜ
using Oceananigans.ImmersedBoundaries: immersed_inactive_node, immersed_peripheral_node
using Oceananigans.Units: Time

radiation = NormalRadiation(inflow_timescale = 6hours, outflow_timescale = Inf)

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
# GLORYS is valid to 3226 m here against a 3176 m domain, so this should be a no-op; the printed
# range below is the check that it is.
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
restoring_kw = (; rate = restoring_rate, mask = boundary_mask, inpainting = bc_inpainting,
                  time_indices_in_memory = n_times)

T_restoring = DW.DatasetRestoring(restoring_meta(:temperature), arch; restoring_kw...)
S_restoring = DW.DatasetRestoring(restoring_meta(:salinity),    arch; restoring_kw...)
u_restoring = DW.DatasetRestoring(restoring_meta(:u_velocity),  arch; restoring_kw...)
v_restoring = DW.DatasetRestoring(restoring_meta(:v_velocity),  arch; restoring_kw...)

# Red Sea S runs 35-41 (hypersaline in the north), T 21-31 -- floors well below anything physical.
fill_down_fts!(T_restoring.field_time_series, 1.0)
fill_down_fts!(S_restoring.field_time_series, 5.0)
@info @sprintf("GLORYS native after fill_down: T ∈ [%.3f, %.3f], S ∈ [%.3f, %.3f]",
               minimum(interior(T_restoring.field_time_series[1])), maximum(interior(T_restoring.field_time_series[1])),
               minimum(interior(S_restoring.field_time_series[1])), maximum(interior(S_restoring.field_time_series[1])))

T_fts, T_native = T_restoring.field_time_series, T_restoring.native_grid
S_fts, S_native = S_restoring.field_time_series, S_restoring.native_grid
u_fts, u_native = u_restoring.field_time_series, u_restoring.native_grid
v_fts, v_native = v_restoring.field_time_series, v_restoring.native_grid

# --- Net-transport correction ------------------------------------------------
# One boundary, so Δu = (Q_raw - Q_target)/A applies to the east section alone.
# Q_target = 0 while there is no surface buoyancy forcing. NOTE the sign convention differs from a
# river-forced domain: the Red Sea is evaporation-driven, so once evaporation is added the target
# becomes a net INFLOW (positive) equal to E, where a river-forced domain needs a net outflow.
net_transport_target = 0.0
const transport_correction = Ref(0.0)

# East boundary: positive u is OUT of the domain, so inflow is -u.
@inline u_east_exterior(j, k, grid, clock, fields) =
    DW.interpolate_to_grid(u_fts, grid.Nx, j, k, u_native, grid, Time(clock.time)) + transport_correction[]

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
                           time_indexing = Oceananigans.OutputReaders.Cyclical(),
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
        uk  = DW.interpolate_to_grid(u_fts, i, j, kk, u_native, grid, t) + transport_correction[]
        U  += ifelse(wet, uk * Δzᶠᶜᶜ(i, j, kk, grid), zero(grid))
    end
    return (U, zos_at(i, j, grid, t))
end

u_bcs = FieldBoundaryConditions(east = NormalFlowBoundaryCondition(u_east_exterior; scheme = radiation, discrete_form = true))
U_bcs = FieldBoundaryConditions(east = GravityWaveRadiationBoundaryCondition(east_U_obc; discrete_form = true))
# The Chapman condition radiates interior η outward and never reads an exterior value -- see the
# recipe. Present for the barotropic mode to be well posed, not to impose anything.
η_bcs = FieldBoundaryConditions(grid, (Center(), Center(), Face()); east = SurfaceWaveRadiationBoundaryCondition())
T_bcs = FieldBoundaryConditions(east = ValueBoundaryCondition(T_east_exterior; scheme = radiation, discrete_form = true))
S_bcs = FieldBoundaryConditions(east = ValueBoundaryCondition(S_east_exterior; scheme = radiation, discrete_form = true))

# --- Ocean simulation --------------------------------------------------------
weno_order = parse(Int, get(ENV, "REDSEA_WENO_ORDER", "7"))
tracer_advection = WENO(order = weno_order,
                        time_discretization = AdaptiveVerticallyImplicitDiscretization(cfl = 0.5))

ocean = ocean_simulation(grid;
                         tracers = (:T, :S),
                         tracer_advection,
                         forcing = (T = T_restoring, S = S_restoring,
                                    u = u_restoring, v = v_restoring),
                         boundary_conditions = (; u = u_bcs, U = U_bcs, η = η_bcs, T = T_bcs, S = S_bcs))

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

function update_transport_correction!(sim)
    t = Time(sim.model.clock.time)
    Q = A = 0.0
    for j in 1:Ny, k in 1:Nz
        wetcell(Nx, j, k, grid, Face(), Center(), Center()) || continue
        a = Axᶠᶜᶜ(Nx, j, k, grid)
        Q -= DW.interpolate_to_grid(u_fts, Nx, j, k, u_native, grid, t) * a
        A += a
    end
    transport_correction[] = A > 0 ? -(Q - net_transport_target) / A : 0.0
    return nothing
end
update_transport_correction!(ocean)
add_callback!(ocean, update_transport_correction!, IterationInterval(10))
@info @sprintf("Initial net-transport correction: Δu = %+.4e m/s (target net = %.1f m³/s)",
               transport_correction[], net_transport_target)

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

const wall_t0   = Ref(Base.time())
const wall_last = Ref(Base.time())
const sim_last  = Ref(0.0)

function progress(sim)
    model = sim.model
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
    @info @sprintf("        transport: east %+.4e m³/s (+ = inflow)  ->  %+.3e m/hour mean sea level",
                   Φ, Φ / wet_surface_area * 3600)

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

    now  = Base.time(); wall = now - wall_t0[]
    Δwall, Δsim = now - wall_last[], time(sim) - sim_last[]
    wall_last[], sim_last[] = now, time(sim)
    rate = Δwall > 0 ? Δsim / Δwall : 0.0
    eta  = (isfinite(sim.stop_time) && rate > 0) ? (sim.stop_time - time(sim)) / rate : NaN
    @info @sprintf("        pacing: wall %s | %.0f× realtime (recent) | %.2f s/iter | ETA %s",
                   prettytime(wall), rate, Δwall / 10, isnan(eta) ? "n/a" : prettytime(eta))
    flush(stderr)
    return nothing
end
add_callback!(ocean, progress, IterationInterval(10))

# --- Spin-up then the run ----------------------------------------------------
ocean.Δt = 10
ocean.stop_iteration = 100
run!(ocean)

ocean_timescale(model) = cell_advection_timescale(model)
wizard = TimeStepWizard(; cfl = 0.15, max_Δt = 15minutes, max_change = 1.1, cell_advection_timescale = ocean_timescale)
ocean.callbacks[:wizard] = Callback(wizard, IterationInterval(10))
ocean.stop_iteration = Inf
stop_days = parse(Float64, get(ENV, "REDSEA_STOP_DAYS", string(run_days)))
ocean.stop_time = stop_days * days

output_fields   = merge(ocean.model.velocities, ocean.model.tracers)
output_schedule = TimeInterval(6hours)
tag = (weno_order == 7 ? "" : "_weno$(weno_order)") * (use_zstar ? "_zstar" : "") *
      (Δdeg == 0.05 ? "" : "_d$(replace(string(Δdeg), "." => ""))")
out(name) = joinpath(@__DIR__, "redsea_glorys_$(name)$(tag)")

ocean.output_writers[:surface] = JLD2Writer(ocean.model, output_fields;
    schedule = output_schedule, overwrite_existing = true, indices = (:, :, Nz), filename = out("surface"))
ocean.output_writers[:east_section] = JLD2Writer(ocean.model, output_fields;
    schedule = output_schedule, overwrite_existing = true, indices = (Nx, :, :), filename = out("east_section"))
ocean.output_writers[:strait_section] = JLD2Writer(ocean.model, output_fields;
    schedule = output_schedule, overwrite_existing = true, indices = (:, j_strait, :), filename = out("strait_section"))
ocean.output_writers[:free_surface] = JLD2Writer(ocean.model, (; η = free_surface_field(ocean.model.free_surface));
    schedule = output_schedule, overwrite_existing = true, filename = out("free_surface"))

run!(ocean)
