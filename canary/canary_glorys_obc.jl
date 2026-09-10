# canary_glorys_obc.jl
# GLORYS-driven open boundaries for the Canary upwelling region (NW Africa).
# Adapted from red_sea/redsea_glorys_obc.jl (monitors, atmosphere, diagnostics) and
# salish_sea/salish_glorys_obc.jl (the multi-boundary pattern). See docs/obc-implementation-recipe.md.
#
# What this domain tests that the others could not:
#   * THREE open boundaries, all in deep open ocean the parent resolves well -- the Red Sea had one,
#     the Salish Sea had two in badly-resolved channels.
#   * A WIND-DRIVEN interior. Upwelling is the signal, so this is the first domain where the
#     atmosphere is the point rather than a correction. It was rejected as a candidate earlier
#     precisely because the configuration had no wind forcing; JRA55 changed that.
#   * Islands in the flow: the Canary archipelago sits inside the domain and sheds wakes the parent
#     cannot resolve, which is the nesting question in a different guise from Bab el-Mandeb.

using CairoMakie
using Oceananigans
using Oceananigans.Architectures: on_architecture
using NumericalEarth
using CopernicusMarine
using Oceananigans.Advection: cell_advection_timescale
using Oceananigans.Units
using CFTime
using Dates
using Printf
using Statistics

const DW  = NumericalEarth.DataWrangling
const SWP = Oceananigans.BuoyancyFormulations.SeawaterPolynomials

# CANARY_ARCH=GPU for the long tracer runs. Everything below is written so the same script runs on
# either: the boundary conditions and the tracer forcing are device kernels already, and every
# DIAGNOSTIC operates on CPU copies (see grid_c / arr below) because GPU disallows the scalar indexing
# they all use.
arch = get(ENV, "CANARY_ARCH", "CPU") == "GPU" ? GPU() : CPU()
@info "Architecture: $(typeof(arch).name.name)"

# --- Run window ---------------------------------------------------------------
# January, matching salish_sea/ and red_sea/ so the three domains share a window. Note the cost:
# Canary upwelling is quasi-permanent but peaks in boreal SUMMER at these latitudes, so a January
# run exercises the machinery on a weaker signal than the domain can produce. CANARY_START=1993-07-01
# for the strong-upwelling case.
run_start   = DateTime(get(ENV, "CANARY_START", "1993-01-01"))
run_days    = parse(Int, get(ENV, "CANARY_RUN_DAYS", "7"))
margin_days = 2
boundary_dates = run_start : Day(1) : (run_start + Day(run_days + margin_days))

# --- Grid ---------------------------------------------------------------------
# Measured from ETOPO, not assumed: the EAST edge at -10.0 is entirely land (the African coast), so
# WEST, SOUTH and NORTH are the open boundaries. West is 100% wet at 4578-4885 m; north 99% wet;
# south 58% wet, its land being the African coast at the SE corner. Combined open area ~7.6e3 km²
# against 719 km² for the Red Sea.
λ₁, λ₂ = (-21.5, -10.0)
φ₁, φ₂ = ( 25.0,  29.5)

# 0.02 deg = 2.06 km at 27.25N, 575 x 225 x 32 = 4.1M cells, i.e. about the Red Sea's cost per step.
# 0.01 deg (1.03 km) is the target for resolving island wakes and the upwelling filaments, but at
# 16.6M cells and half the timestep it is a GPU/HPC run, not a CPU one.
Δdeg = parse(Float64, get(ENV, "CANARY_DELTA", "0.02"))
Nx = round(Int, (λ₂ - λ₁) / Δdeg)
Ny = round(Int, (φ₂ - φ₁) / Δdeg)

# Vertical. GLORYS's deepest valid level over this region is -4836.6 m while raw ETOPO reaches
# 4885 m, so unlike the Red Sea the margin runs the WRONG WAY and `depth` is a hard cap rather than
# a convenience: 4800 keeps every model cell inside the parent's valid data, at the cost of clipping
# 1.2% of columns by at most 85 m of abyssal plain.
# scale=700 spends the 32 levels where upwelling lives: 1.2 m surface cell, 11 levels in the top
# 50 m, 16 in the top 150 m. The price is ~930 m bottom cells, which the quiescent abyss can carry.
Nz = 32
depth   = parse(Float64, get(ENV, "CANARY_DEPTH",  "4800"))
z_scale = parse(Float64, get(ENV, "CANARY_ZSCALE", "700"))
z_reference = ExponentialDiscretization(Nz, -depth, 0; scale = z_scale)

# z* on by default, as in red_sea/: with a 1.2 m surface cell the O(max|η|/Δz_top) error of a static
# (linear free surface) grid is not small. CANARY_ZSTAR=0 falls back to static z.
use_zstar = get(ENV, "CANARY_ZSTAR", "1") == "1"
z = use_zstar ? MutableVerticalDiscretization(z_reference) : z_reference

@info @sprintf("Grid: %d x %d x %d (Δ = %.3f deg ≈ %.2f km at 27.25N) | depth %.0f m | vertical: %s",
               Nx, Ny, Nz, Δdeg, Δdeg * 111.32 * cosd(27.25), depth,
               use_zstar ? "ZStar (mutable)" : "Z (static)")

underlying_grid = LatitudeLongitudeGrid(arch;
                             size = (Nx, Ny, Nz),
                             latitude = (φ₁, φ₂), longitude = (λ₁, λ₂),
                             z, halo = (7, 7, 7),
                             topology = (Bounded, Bounded, Bounded))

# The islands are the reason to keep smoothing light: at 0.02 deg El Hierro is only a few cells
# across, and over-smoothing would dissolve the wake-generating topography this domain exists to
# resolve. The startup check below reports each island's land-cell count.
interpolation_passes = parse(Int, get(ENV, "CANARY_PASSES", "5"))
bottom_height = regrid_bathymetry(underlying_grid;
                                  minimum_depth = 10,
                                  interpolation_passes = interpolation_passes)
grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bottom_height))

# CPU mirror. The diagnostics all walk the grid cell by cell, and the metric operators (Vᶜᶜᶜ, Axᶠᶜᶜ,
# Azᶜᶜᶜ …) read the grid's own coordinate arrays, so BOTH the fields and the grid must be on CPU for
# them. On a CPU run this is the same grid and `arr` is a plain copy.
const grid_c = on_architecture(CPU(), grid)
@inline arr(f) = Array(interior(f))

# On CPU the mirror SHARES the z* state arrays (σ, η, ∂t_σ), so its metrics track the live grid for
# free. On GPU it is a genuine copy, and without this every diagnostic would silently use FROZEN
# metrics -- including salt_and_volume, which exists to measure z* drift and would then read ~0 no
# matter what z* did. Called once per progress block; a no-op on CPU.
const _ZSTAR_FIELDS = (:ηⁿ, :σᶜᶜⁿ, :σᶠᶜⁿ, :σᶜᶠⁿ, :σᶠᶠⁿ, :σᶜᶜ⁻, :∂t_σ)
function sync_grid_c!()
    (arch isa CPU || !(grid.underlying_grid.z isa MutableVerticalDiscretization)) && return nothing
    zg, zc = grid.underlying_grid.z, grid_c.underlying_grid.z
    for f in _ZSTAR_FIELDS
        hasfield(typeof(zg), f) || continue
        copyto!(parent(getfield(zc, f)), Array(parent(getfield(zg, f))))
    end
    return nothing
end

# Fixed subset so resolution changes do not invalidate the download cache.
glorys_region = DW.BoundingBox(longitude = (-22.0, -9.5), latitude = (24.5, 30.0))
@assert λ₁ >= -22.0 && λ₂ <= -9.5 && φ₁ >= 24.5 && φ₂ <= 30.0 "domain outside glorys_region"
glorys = DW.GLORYSDaily()

# --- Passive tracer release (optional) ----------------------------------------
# For pathway and residence-time work: a continuous release at a point on Gran Canaria, of the kind
# an ocean alkalinity enhancement outfall would represent.
#
# TWO tracers, following Deleersnijder et al.'s age formulation:
#   c  -- concentration.      ∂c/∂t  + ∇·(uc)  = S      with S the release
#   ac -- age concentration.  ∂ac/∂t + ∇·(uac) = c      (grows at rate c wherever tracer is present)
# so the MEAN AGE of the water carrying the tracer is ac/c. A single "ideal age" tracer cannot do
# this for a point source -- it measures time since surface contact, not time since release.
#
# S is normalised to ∫S dV = 1 unit/s over the release volume, which makes the budget checkable by
# inspection: total inventory should equal elapsed seconds minus whatever has left through the open
# boundaries. Default site is Las Palmas (-15.38, 28.16): wet, 535 m, 47 of 49 neighbouring cells wet,
# and the island's main port.
use_tracer  = get(ENV, "CANARY_TRACER", "0") == "1"
λ_release   = parse(Float64, get(ENV, "CANARY_RELEASE_LON",   "-15.38"))
φ_release   = parse(Float64, get(ENV, "CANARY_RELEASE_LAT",    "28.16"))
R_release   = parse(Float64, get(ENV, "CANARY_RELEASE_R",       "0.05"))   # degrees, ~5 km
D_release   = parse(Float64, get(ENV, "CANARY_RELEASE_DEPTH",  "20.0"))    # metres below surface

# --- Boundary conditions ------------------------------------------------------
using Oceananigans.BoundaryConditions: NormalRadiation,
                                       GravityWaveRadiationBoundaryCondition,
                                       SurfaceWaveRadiationBoundaryCondition
using Oceananigans.Operators: Δzᶜᶜᶜ, Δzᶠᶜᶜ, Δzᶜᶠᶜ, Axᶠᶜᶜ, Ayᶜᶠᶜ, Azᶜᶜᶜ, Vᶜᶜᶜ
using Oceananigans.ImmersedBoundaries: immersed_inactive_node, immersed_peripheral_node
using Oceananigans.Units: Time

radiation = NormalRadiation(inflow_timescale = 6hours, outflow_timescale = Inf)

@inline wetcell(i, j, k, grid, ℓx, ℓy, ℓz) =
    !immersed_peripheral_node(i, j, k, grid, ℓx, ℓy, ℓz) & !immersed_inactive_node(i, j, k, grid, ℓx, ℓy, ℓz)

# --- Nudging band -------------------------------------------------------------
# Bands on all three open edges, combined with max() so the corners do not double-count. Kept narrow
# (0.25 deg ≈ 12 cells at 0.02) so the shelf and the island wakes run free: the Salish lesson was
# that a wide sponge pins the very features the domain exists to resolve, and here the upwelling
# front sits against the eastern coast, far from every sponge.
sponge_width        = parse(Float64, get(ENV, "CANARY_SPONGE", "0.25"))
restoring_timescale = 12hours
restoring_rate      = 1 / restoring_timescale

@inline gaussian_edge(d, w) = exp(-(d / w)^2)
@inline function boundary_mask(x, y, z, t)
    dw = x - λ₁        # west
    ds = y - φ₁        # south
    dn = φ₂ - y        # north
    return max(gaussian_edge(dw, sponge_width),
               gaussian_edge(ds, sponge_width),
               gaussian_edge(dn, sponge_width))
end

restoring_meta(var) = DW.Metadata(var; dataset = glorys, dates = boundary_dates, region = glorys_region)
bc_inpainting = DW.NearestNeighborInpainting(100)

# Vertical constant extrapolation below the parent's deepest valid level. Insurance in the Red Sea;
# here it is expected to do real work, because `depth` sits only ~37 m above the parent's deepest
# valid centre and the abyssal columns are sparsely covered (2.9% of columns have data at -4836 m).
function fill_down!(a::AbstractArray{<:Any,3}, floor_value)
    nx, ny, nz = size(a)
    @inbounds for j in 1:ny, i in 1:nx
        last = NaN
        for k in nz:-1:1
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
# via a host copy: fill_down! walks the array cell by cell, which a device array will not allow.
function fill_down_fts!(fts, floor)
    for n in 1:length(fts.times)
        a = Array(interior(fts[n])); fill_down!(a, floor); copyto!(interior(fts[n]), a)
    end
    return fts
end

n_times = length(boundary_dates)
restoring_kw = (; rate = restoring_rate, mask = boundary_mask, inpainting = bc_inpainting,
                  time_indices_in_memory = n_times)

T_restoring = DW.DatasetRestoring(restoring_meta(:temperature), arch; restoring_kw...)
S_restoring = DW.DatasetRestoring(restoring_meta(:salinity),    arch; restoring_kw...)
u_restoring = DW.DatasetRestoring(restoring_meta(:u_velocity),  arch; restoring_kw...)
v_restoring = DW.DatasetRestoring(restoring_meta(:v_velocity),  arch; restoring_kw...)

# Canary Current water is S 35.5-36.8, T 15-24 at the surface; floors sit well below anything real.
fill_down_fts!(T_restoring.field_time_series, 1.0)
fill_down_fts!(S_restoring.field_time_series, 5.0)
let Ti = interior(T_restoring.field_time_series[1]), Si = interior(S_restoring.field_time_series[1])
    @info @sprintf("GLORYS native after fill_down: T ∈ [%.3f, %.3f], S ∈ [%.3f, %.3f] | residual zeros %.2f%%",
                   minimum(Ti), maximum(Ti), minimum(Si), maximum(Si), 100count(==(0), Si) / length(Si))
end

T_fts, T_native = T_restoring.field_time_series, T_restoring.native_grid
S_fts, S_native = S_restoring.field_time_series, S_restoring.native_grid
u_fts, u_native = u_restoring.field_time_series, u_restoring.native_grid
v_fts, v_native = v_restoring.field_time_series, v_restoring.native_grid

# --- GLORYS section residual (diagnostic only) --------------------------------
# No net-transport correction: measured in red_sea/ to change nothing a Flather boundary does not
# already absorb (see docs/obc-implementation-recipe.md). Q_raw is reported per boundary instead,
# because with three sections it also says which edge the parent disagrees with most.

# --- Exteriors ----------------------------------------------------------------
# Inflow is +u at west, +v at south, -v at north. Oceananigans applies the side-dependent Flather
# sign itself (west/south: U = Uᵉˣᵗ − √(gH)(ηᵇ − ηᵉˣᵗ); north/east: +), so these closures return the
# exterior value only. Do NOT pre-flip signs here -- that double inversion is exactly the bug the
# Red Sea configuration carried through its whole validation run.
@inline u_west_exterior(j, k, grid, clock, fields) =
    DW.interpolate_to_grid(u_fts, 1, j, k, u_native, grid, Time(clock.time))
@inline v_south_exterior(i, k, grid, clock, fields) =
    DW.interpolate_to_grid(v_fts, i, 1, k, v_native, grid, Time(clock.time))
@inline v_north_exterior(i, k, grid, clock, fields) =
    DW.interpolate_to_grid(v_fts, i, grid.Ny + 1, k, v_native, grid, Time(clock.time))

@inline T_west_exterior(j, k, grid, clock, fields)  = DW.interpolate_to_grid(T_fts, 1, j, k, T_native, grid, Time(clock.time))
@inline S_west_exterior(j, k, grid, clock, fields)  = DW.interpolate_to_grid(S_fts, 1, j, k, S_native, grid, Time(clock.time))
@inline T_south_exterior(i, k, grid, clock, fields) = DW.interpolate_to_grid(T_fts, i, 1, k, T_native, grid, Time(clock.time))
@inline S_south_exterior(i, k, grid, clock, fields) = DW.interpolate_to_grid(S_fts, i, 1, k, S_native, grid, Time(clock.time))
@inline T_north_exterior(i, k, grid, clock, fields) = DW.interpolate_to_grid(T_fts, i, grid.Ny, k, T_native, grid, Time(clock.time))
@inline S_north_exterior(i, k, grid, clock, fields) = DW.interpolate_to_grid(S_fts, i, grid.Ny, k, S_native, grid, Time(clock.time))

zos_meta = DW.Metadata(:free_surface; dataset = glorys, dates = boundary_dates, region = glorys_region)
zos_fts  = FieldTimeSeries(zos_meta, arch;
                           time_indices_in_memory = n_times,
                           time_indexing = Oceananigans.OutputReaders.Cyclical(),
                           inpainting = bc_inpainting)
zos_native = zos_fts.grid
@inline zos_at(i, j, grid, t) = DW.interpolate_to_grid(zos_fts, i, j, 1, zos_native, grid, t)

# Depth-integrated exterior transport for the Flather condition, one per open edge. Velocity is
# sampled at the boundary FACE index (west i=1, south j=1, north j=Ny+1) and η at the adjacent
# interior node, matching salish_sea/.
@inline function west_U_obc(j, k, grid, clock, fields)
    i = 1; t = Time(clock.time); U = zero(grid)
    @inbounds for kk in 1:grid.Nz
        wet = wetcell(i, j, kk, grid, Face(), Center(), Center())
        U += ifelse(wet, DW.interpolate_to_grid(u_fts, i, j, kk, u_native, grid, t) * Δzᶠᶜᶜ(i, j, kk, grid), zero(grid))
    end
    return (U, zos_at(i, j, grid, t))
end

@inline function south_V_obc(i, k, grid, clock, fields)
    j = 1; t = Time(clock.time); V = zero(grid)
    @inbounds for kk in 1:grid.Nz
        wet = wetcell(i, j, kk, grid, Center(), Face(), Center())
        V += ifelse(wet, DW.interpolate_to_grid(v_fts, i, j, kk, v_native, grid, t) * Δzᶜᶠᶜ(i, j, kk, grid), zero(grid))
    end
    return (V, zos_at(i, j, grid, t))
end

@inline function north_V_obc(i, k, grid, clock, fields)
    j = grid.Ny + 1; t = Time(clock.time); V = zero(grid)
    @inbounds for kk in 1:grid.Nz
        wet = wetcell(i, j, kk, grid, Center(), Face(), Center())
        V += ifelse(wet, DW.interpolate_to_grid(v_fts, i, j, kk, v_native, grid, t) * Δzᶜᶠᶜ(i, j, kk, grid), zero(grid))
    end
    return (V, zos_at(i, grid.Ny, grid, t))
end

u_bcs = FieldBoundaryConditions(west  = NormalFlowBoundaryCondition(u_west_exterior;  scheme = radiation, discrete_form = true))
v_bcs = FieldBoundaryConditions(south = NormalFlowBoundaryCondition(v_south_exterior; scheme = radiation, discrete_form = true),
                                north = NormalFlowBoundaryCondition(v_north_exterior; scheme = radiation, discrete_form = true))
U_bcs = FieldBoundaryConditions(grid, (Face(), Center(), nothing);
                                west  = GravityWaveRadiationBoundaryCondition(west_U_obc;  discrete_form = true))
V_bcs = FieldBoundaryConditions(grid, (Center(), Face(), nothing);
                                south = GravityWaveRadiationBoundaryCondition(south_V_obc; discrete_form = true),
                                north = GravityWaveRadiationBoundaryCondition(north_V_obc; discrete_form = true))
# Chapman: radiates interior η outward and never reads an exterior value. Present so the barotropic
# mode is well posed, not to impose anything -- see the recipe.
η_bcs = FieldBoundaryConditions(grid, (Center(), Center(), Face());
                                west  = SurfaceWaveRadiationBoundaryCondition(),
                                south = SurfaceWaveRadiationBoundaryCondition(),
                                north = SurfaceWaveRadiationBoundaryCondition())
T_bcs = FieldBoundaryConditions(west  = ValueBoundaryCondition(T_west_exterior;  scheme = radiation, discrete_form = true),
                                south = ValueBoundaryCondition(T_south_exterior; scheme = radiation, discrete_form = true),
                                north = ValueBoundaryCondition(T_north_exterior; scheme = radiation, discrete_form = true))
S_bcs = FieldBoundaryConditions(west  = ValueBoundaryCondition(S_west_exterior;  scheme = radiation, discrete_form = true),
                                south = ValueBoundaryCondition(S_south_exterior; scheme = radiation, discrete_form = true),
                                north = ValueBoundaryCondition(S_north_exterior; scheme = radiation, discrete_form = true))

# --- Tracer release forcing and boundaries ------------------------------------
# Release volume is measured on the grid so the source can be normalised; a continuous-form forcing
# (a pure function of x, y, z, t) is used rather than a masked array so the same code runs on GPU.
# In a function, not a top-level loop: accumulating into a global from inside `for` hits Julia's
# soft-scope rule and fails with UndefVarError.
function compute_release_volume()
    λc = λnodes(grid, Center(), Center(), Center())
    φc = φnodes(grid, Center(), Center(), Center())
    zc = znodes(grid, Center(), Center(), Center())
    V = 0.0
    for k in 1:Nz, j in 1:Ny, i in 1:Nx
        wetcell(i, j, k, grid_c, Center(), Center(), Center()) || continue
        dx = (λc[i] - λ_release) * cosd(φ_release)
        dy = φc[j] - φ_release
        (dx*dx + dy*dy < R_release^2 && zc[k] > -D_release) || continue
        V += Vᶜᶜᶜ(i, j, k, grid_c)
    end
    return V
end

release_volume = use_tracer ? compute_release_volume() : 0.0
if use_tracer
    release_volume <= 0 ?
        @warn("Release site ($(λ_release), $(φ_release)) contains NO wet cells -- tracer will never enter the domain") :
        @info @sprintf("Tracer release at %.3fE %.3fN, r=%.3f deg, top %.0f m: volume %.4e m³",
                       λ_release, φ_release, R_release, D_release, release_volume)
end

const RELEASE_RATE = release_volume > 0 ? 1 / release_volume : 0.0   # ∫S dV = 1 unit/s

@inline function c_source(x, y, z, t)
    dx = (x - λ_release) * cosd(φ_release)
    dy = y - φ_release
    inside = (dx*dx + dy*dy < R_release^2) & (z > -D_release)
    return ifelse(inside, RELEASE_RATE, zero(x))
end
@inline age_source(x, y, z, t, c) = c        # ∂ac/∂t = c, so ac/c is the mean age

# Zero exterior at every open boundary: there is no released tracer outside the domain, so this lets
# it leave freely. A default no-flux condition would trap it and make residence time meaningless.
@inline zero_west(j, k, grid, clock, fields)  = zero(grid)
@inline zero_ns(i, k, grid, clock, fields)    = zero(grid)
tracer_bcs() = FieldBoundaryConditions(
    west  = ValueBoundaryCondition(zero_west; scheme = radiation, discrete_form = true),
    south = ValueBoundaryCondition(zero_ns;   scheme = radiation, discrete_form = true),
    north = ValueBoundaryCondition(zero_ns;   scheme = radiation, discrete_form = true))

tracer_names   = use_tracer ? (:T, :S, :c, :ac) : (:T, :S)
tracer_forcing = use_tracer ? (; c = Forcing(c_source), ac = Forcing(age_source, field_dependencies = :c)) : (;)
tracer_bc_pair = use_tracer ? (; c = tracer_bcs(), ac = tracer_bcs()) : (;)

# --- Ocean simulation ---------------------------------------------------------
weno_order = parse(Int, get(ENV, "CANARY_WENO_ORDER", "7"))
tracer_advection = WENO(order = weno_order,
                        time_discretization = AdaptiveVerticallyImplicitDiscretization(cfl = 0.5))

ocean = ocean_simulation(grid;
                         tracers = tracer_names,
                         tracer_advection,
                         forcing = merge((T = T_restoring, S = S_restoring,
                                          u = u_restoring, v = v_restoring), tracer_forcing),
                         boundary_conditions = merge((; u = u_bcs, U = U_bcs, v = v_bcs, V = V_bcs,
                                                        η = η_bcs, T = T_bcs, S = S_bcs), tracer_bc_pair))

# --- Atmospheric forcing ------------------------------------------------------
# The point of this domain, not a correction to it: upwelling is wind-driven, so CANARY_ATMOS=none
# gives a boundary-only control run rather than a usable configuration.
# API notes and the ocean_simulation BC-merge trap are documented in red_sea/README.md.
atmos_source = get(ENV, "CANARY_ATMOS", "JRA")
use_atmos = atmos_source == "JRA"
jra_multiyear = get(ENV, "CANARY_JRA_DATASET", "repeat") == "multiyear"

atmosphere = nothing
radiation_component = nothing
if use_atmos
    jra_dataset = jra_multiyear ? NumericalEarth.MultiYearJRA55() : NumericalEarth.RepeatYearJRA55()
    jra_year  = jra_multiyear ? Dates.year(run_start) : 1990
    jra_start = DateTime(jra_year, Dates.month(run_start), Dates.day(run_start))
    jra_stop  = jra_start + Day(run_days + margin_days)
    @info @sprintf("JRA55 %s: %s .. %s", jra_multiyear ? "MultiYear" : "RepeatYear", jra_start, jra_stop)
    atmos_kw = (; dataset = jra_dataset, start_date = jra_start, end_date = jra_stop,
                  time_indices_in_memory = 24)
    atmosphere = NumericalEarth.JRA55PrescribedAtmosphere(arch; atmos_kw...)
    radiation_component = NumericalEarth.JRA55PrescribedRadiation(arch; atmos_kw...)
end

# --- Initial conditions -------------------------------------------------------
# Sampled from the SAME repaired series the nudging and boundaries use, so IC and forcing agree by
# construction; the alternative re-reads the cache and never sees the repair.
#
# On GPU the restoring series live on the device and this samples them point by point, which the
# device disallows -- so CPU copies are loaded for the purpose, with the same fill_down! repair. They
# are only needed at setup and cost a few tens of MB. On CPU the existing series are reused.
function cpu_parent(var, floor_value)
    f = FieldTimeSeries(restoring_meta(var), CPU(); time_indices_in_memory = n_times, inpainting = bc_inpainting)
    fill_down_fts!(f, floor_value)
    return f
end
T_fts_c = arch isa CPU ? T_fts : cpu_parent(:temperature, 1.0)
S_fts_c = arch isa CPU ? S_fts : cpu_parent(:salinity,    5.0)

function ic_array(fts, native_grid, t)
    a = zeros(eltype(interior(fts[1])), Nx, Ny, Nz)
    @inbounds for k in 1:Nz, j in 1:Ny, i in 1:Nx
        a[i, j, k] = DW.interpolate_to_grid(fts, i, j, k, native_grid, grid_c, t)
    end
    return a
end

set!(ocean.model.tracers.T, ic_array(T_fts_c, T_fts_c.grid, Time(0.0)))
set!(ocean.model.tracers.S, ic_array(S_fts_c, S_fts_c.grid, Time(0.0)))

free_surface_field(fs) = hasproperty(fs, :η) ? fs.η : fs.displacement
η_field = free_surface_field(ocean.model.free_surface)
# zos_fts itself must stay on the device -- the Flather closures sample it there. Building η samples
# it point by point on the host instead, so that needs a host copy.
zos_fts_c = arch isa CPU ? zos_fts :
    FieldTimeSeries(zos_meta, CPU(); time_indices_in_memory = n_times,
                    time_indexing = Oceananigans.OutputReaders.Cyclical(), inpainting = bc_inpainting)
let t₀ = Time(0.0)
    ηa = zeros(Nx, Ny, 1)
    for j in 1:Ny, i in 1:Nx
        wet = !immersed_inactive_node(i, j, Nz, grid_c, Center(), Center(), Center())
        ηa[i, j, 1] = ifelse(wet, DW.interpolate_to_grid(zos_fts_c, i, j, 1, zos_fts_c.grid, grid_c, t₀), 0.0)
    end
    set!(η_field, ηa)
end
@info @sprintf("η from GLORYS zos: [%.4f, %.4f] m", minimum(arr(η_field)), maximum(arr(η_field)))

if use_atmos
    coupled_model = NumericalEarth.OceanOnlyModel(ocean; atmosphere, radiation = radiation_component)
    simulation = Simulation(coupled_model; Δt = 10, stop_iteration = 100)
else
    simulation = ocean
end
ocean_model_of(m) = hasproperty(m, :ocean) ? (hasproperty(m.ocean, :model) ? m.ocean.model : m.ocean) : m

# --- Diagnostics --------------------------------------------------------------
z_centers = znodes(grid, Center(), Center(), Center())
λ_centers = λnodes(grid, Center(), Center(), Center())
φ_centers = φnodes(grid, Center(), Center(), Center())

@info @sprintf("Δz near surface: %s m (top of %d levels, domain %.0f m)",
               join((@sprintf("%.1f", d) for d in reverse(diff(z_centers))[1:min(6, Nz-1)]), ", "), Nz, depth)

wet3d = [!immersed_inactive_node(i, j, k, grid_c, Center(), Center(), Center()) for i in 1:Nx, j in 1:Ny, k in 1:Nz]
wet_idx = findall(wet3d)
@info "Wet cells: $(length(wet_idx)) of $(Nx*Ny*Nz) ($(round(100*length(wet_idx)/(Nx*Ny*Nz), digits=1))%)"
wet_surface_area = sum(Azᶜᶜᶜ(i, j, Nz, grid_c) for i in 1:Nx, j in 1:Ny
                       if !immersed_inactive_node(i, j, Nz, grid_c, Center(), Center(), Center()); init = 0.0)

# Islands: the topography this domain exists to resolve. Over-smoothing dissolves them silently, so
# report the land-cell count in a 0.4 deg box around each before trusting a run.
for (nm, il, ip) in (("La Palma", -17.87, 28.68), ("El Hierro", -18.01, 27.74), ("Tenerife", -16.60, 28.29),
                     ("Gran Canaria", -15.60, 27.96), ("Fuerteventura", -14.02, 28.36), ("Lanzarote", -13.62, 29.03))
    box = [(i, j) for i in 1:Nx, j in 1:Ny if abs(λ_centers[i]-il) < 0.2 && abs(φ_centers[j]-ip) < 0.2]
    isempty(box) && continue
    nl = count(!wet3d[i, j, Nz] for (i, j) in box)
    nl == 0 ? (@warn "$nm is NOT RESOLVED at this resolution/smoothing -- reduce CANARY_PASSES") :
              @info @sprintf("  island %-14s %3d land of %3d cells", nm, nl, length(box))
end

function wet_extrema(a3)          # a3 is a CPU array, not a Field
    a = a3[wet_idx]
    nmin, nmax = argmin(a), argmax(a)
    return a[nmin], wet_idx[nmin], a[nmax], wet_idx[nmax]
end

# Tag by the feature that matters here: which open edge, the upwelling shelf, or the island zone.
function whereis(I)
    i, j, k = Tuple(I)
    λi, φj = λ_centers[clamp(i,1,Nx)], φ_centers[clamp(j,1,Ny)]
    tag = i <= 2                     ? "WEST-bnd"  :
          j <= 2                     ? "SOUTH-bnd" :
          j >= Ny-1                  ? "NORTH-bnd" :
          (λi > -14.0)               ? "SHELF/upwelling" :
          (-18.3 < λi < -13.2 && 27.4 < φj < 29.3) ? "ISLANDS" : "interior"
    return @sprintf("(%3d,%3d,%2d) z=%7.1fm %s", i, j, k, z_centers[clamp(k,1,Nz)], tag)
end

# Volume transport through each open edge, positive INTO the domain.
function boundary_transports(u, v)
    W = S_ = N = 0.0
    @inbounds for k in 1:Nz, j in 1:Ny
        wetcell(1, j, k, grid_c, Face(), Center(), Center()) && (W += u[1, j, k] * Axᶠᶜᶜ(1, j, k, grid_c))
    end
    @inbounds for k in 1:Nz, i in 1:Nx
        wetcell(i, 1, k, grid_c, Center(), Face(), Center())  && (S_ += v[i, 1, k]  * Ayᶜᶠᶜ(i, 1, k, grid_c))
        # Ny+1, not Ny: for Ny cells a y-Face field has Ny+1 points and the north DOMAIN EDGE is the
        # last one. Sampling Ny measures the flux through an interior face instead, which is not the
        # domain's volume budget. (north_V_obc already uses Ny+1.)
        wetcell(i, Ny+1, k, grid_c, Center(), Face(), Center()) && (N -= v[i, Ny+1, k] * Ayᶜᶠᶜ(i, Ny+1, k, grid_c))
    end
    return W, S_, N
end

# What the raw GLORYS exterior would carry through each edge: a parent/child mismatch measure, not a
# physical signal. Reported, never corrected.
#
# Computed live, not precomputed. An earlier version evaluated it once at setup and interpolated in
# time; the parent interpolation is exactly linear so that part was sound, but it also froze the
# SECTION AREAS at their rest state, and under z* those breathe with the free surface. That shifted
# the diagnostic by a few tenths of a percent against the validated baseline. Host copies of the
# parent series keep this off the device without changing what is measured.
u_fts_c = arch isa CPU ? u_fts : FieldTimeSeries(restoring_meta(:u_velocity), CPU(); time_indices_in_memory = n_times, inpainting = bc_inpainting)
v_fts_c = arch isa CPU ? v_fts : FieldTimeSeries(restoring_meta(:v_velocity), CPU(); time_indices_in_memory = n_times, inpainting = bc_inpainting)

function glorys_residuals(t)
    W = S_ = N = 0.0
    @inbounds for k in 1:Nz, j in 1:Ny
        wetcell(1, j, k, grid_c, Face(), Center(), Center()) &&
            (W += DW.interpolate_to_grid(u_fts_c, 1, j, k, u_fts_c.grid, grid_c, t) * Axᶠᶜᶜ(1, j, k, grid_c))
    end
    @inbounds for k in 1:Nz, i in 1:Nx
        wetcell(i, 1, k, grid_c, Center(), Face(), Center()) &&
            (S_ += DW.interpolate_to_grid(v_fts_c, i, 1, k, v_fts_c.grid, grid_c, t) * Ayᶜᶠᶜ(i, 1, k, grid_c))
        wetcell(i, Ny+1, k, grid_c, Center(), Face(), Center()) &&
            (N -= DW.interpolate_to_grid(v_fts_c, i, Ny+1, k, v_fts_c.grid, grid_c, t) * Ayᶜᶠᶜ(i, Ny+1, k, grid_c))
    end
    return W, S_, N
end

# THE signal for this domain. Upwelling is offshore Ekman transport drawing cold water up against the
# African coast, so the diagnostic is the cross-shore SST contrast: coastal minus offshore at MATCHED
# latitude, plus the temperature at the ~100 m source depth over the shelf.
#
# The coast runs diagonally -- λ = -14.87 at 25N to -10.23 at 29.4N -- so a fixed-longitude strip
# cannot follow it. The first version used λ > -13.5 and sampled 0 wet cells at 25-26N and 7 at
# 26-27N, missing Cape Bojador (~26N) entirely, i.e. one of the two upwelling centres in the domain.
# `upwelling_index` therefore walks in from the coast ROW BY ROW. `upwelling_index_fixedλ` keeps the
# original definition so July stays comparable with the January run, which used it.
const k_src = argmin(abs.(z_centers .+ 100))     # ~100 m, the source depth for coastal upwelling
# Bands are defined in DEGREES and converted to cells, not fixed cell counts: a hard-coded 150-cell
# offshore offset silently exceeded the grid_c width at 0.1° and returned NaN for every row.
const N_COAST     = max(1, round(Int, 0.50 / Δdeg))   # ~50 km inshore band
const N_OFFSHORE  = round(Int, 3.00 / Δdeg)           # reference centred ~300 km offshore
const N_OFF_HALF  = max(1, round(Int, 0.25 / Δdeg))   # ~25 km half-width on that reference
# Easternmost wet surface cell per row = the coast. NaN where the row has no water.
const i_coast = [findlast(i -> wet3d[i, j, Nz], 1:Nx) for j in 1:Ny]

function upwelling_index(T)
    cT = oT = sT = 0.0; nc = no = ns = 0
    @inbounds for j in 1:Ny
        ic = i_coast[j]
        ic === nothing && continue
        for i in max(1, ic-N_COAST+1):ic
            wet3d[i, j, Nz] || continue
            cT += T[i, j, Nz]; nc += 1
            if wet3d[i, j, k_src]; sT += T[i, j, k_src]; ns += 1; end
        end
        lo, hi = ic - N_OFFSHORE - N_OFF_HALF, ic - N_OFFSHORE + N_OFF_HALF
        lo < 1 && continue
        for i in lo:hi
            wet3d[i, j, Nz] || continue
            oT += T[i, j, Nz]; no += 1
        end
    end
    Tc = nc > 0 ? cT/nc : NaN
    To = no > 0 ? oT/no : NaN
    return Tc, To, Tc - To, ns > 0 ? sT/ns : NaN
end

# Original fixed-longitude definition, kept only so the January run stays comparable.
function upwelling_index_fixedλ(T)
    cT = oT = 0.0; nc = no = 0
    @inbounds for j in 1:Ny, i in 1:Nx
        wet3d[i, j, Nz] || continue
        λi = λ_centers[i]
        if λi > -13.5;              cT += T[i, j, Nz]; nc += 1
        elseif -20.5 < λi < -19.0;  oT += T[i, j, Nz]; no += 1 end
    end
    return (nc > 0 ? cT/nc : NaN) - (no > 0 ? oT/no : NaN)
end

# --- Tracer budget ------------------------------------------------------------
# The source is normalised to 1 unit/s, so cumulative input is exactly the elapsed time in seconds:
#   inventory + cumulative export == t
# which makes the budget checkable by inspection rather than by trusting the advection scheme.
# Mean age is ac/c over the whole tracer-bearing volume (Deleersnijder); it is only meaningful once
# there is tracer to speak of, so it is reported as NaN until the inventory is non-trivial.
function tracer_budget(c, ac)
    C = 0.0; AC = 0.0; cmax = -Inf; imax = CartesianIndex(1, 1, 1)
    @inbounds for k in 1:Nz, j in 1:Ny, i in 1:Nx
        wetcell(i, j, k, grid_c, Center(), Center(), Center()) || continue
        v = Vᶜᶜᶜ(i, j, k, grid_c)
        C += v * c[i, j, k]; AC += v * ac[i, j, k]
        if c[i, j, k] > cmax; cmax = c[i, j, k]; imax = CartesianIndex(i, j, k); end
    end
    return C, AC, cmax, imax
end

# Advective tracer flux through each open edge; + = into the domain, so a negative sum is export.
function tracer_export(c, u, v)
    W = 0.0; So = 0.0; N = 0.0
    @inbounds for k in 1:Nz, j in 1:Ny
        wetcell(1, j, k, grid_c, Face(), Center(), Center()) &&
            (W += u[1, j, k] * Axᶠᶜᶜ(1, j, k, grid_c) * c[1, j, k])
    end
    @inbounds for k in 1:Nz, i in 1:Nx
        wetcell(i, 1, k, grid_c, Center(), Face(), Center()) &&
            (So += v[i, 1, k] * Ayᶜᶠᶜ(i, 1, k, grid_c) * c[i, 1, k])
        wetcell(i, Ny+1, k, grid_c, Center(), Face(), Center()) &&
            (N -= v[i, Ny+1, k] * Ayᶜᶠᶜ(i, Ny+1, k, grid_c) * c[i, Ny, k])
    end
    return W, So, N
end

# --- Conservation and stability monitors --------------------------------------
function salt_and_volume(S)
    V = 0.0; SV = 0.0
    @inbounds for k in 1:Nz, j in 1:Ny, i in 1:Nx
        wetcell(i, j, k, grid_c, Center(), Center(), Center()) || continue
        v = Vᶜᶜᶜ(i, j, k, grid_c)
        V += v; SV += v * S[i, j, k]
    end
    return SV, V
end

const eos = SWP.TEOS10.TEOS10EquationOfState()
# Transport split by POTENTIAL density, for comparison with the Canary Current literature.
# The inverse box models (Hernández-Guerra 2005; Pérez-Hernández 2013; Casanova-Masjoan 2020) report
# transports over NEUTRAL density (γⁿ) layers -- thermocline γⁿ < 27.38, intermediate 27.38–27.922,
# deep > 27.922. γⁿ needs the Jackett & McDougall lookup, which this stack does not have, so σθ =
# ρ(T,S,0) − 1000 stands in at the same thresholds. Good at the 27.38 thermocline boundary; it drifts
# with depth, so treat the intermediate/deep split as indicative only.
#
# Fixed-DEPTH bins would be the obvious alternative and are wrong for this purpose: isopycnals heave
# seasonally, so a fixed split aliases that heave into apparent transport change -- corrupting exactly
# the January-vs-July comparison this diagnostic exists to make.
#
# Literature check for the thermocline layer, our sign convention (+ = into the domain): the Canary
# Current flows SOUTHWARD, so expect north > 0 (inflow) and south < 0 (outflow), a few Sv each, weaker
# in January than July. Hernández-Guerra 2005 give 4.7 ± 0.8 Sv southward north of the islands.
# Thresholds match the DEPTH the literature's layers occupy, not their γⁿ numbers -- σθ ≠ γⁿ, so
# copying the numbers copies the label and not the layer. Measured in the January run: σθ=27.38 sits
# at 1195 m against the literature's ~700 m thermocline base, while σθ=27.10 sits at ~500 m. The
# 27.38 bin therefore straddled TWO OPPOSING CURRENTS -- the Canary Current above ~500 m (north
# inflow) and the poleward undercurrent at 700-1500 m (north outflow) -- and as the undercurrent spun
# up it flipped the sum from +8.5 to -0.5 Sv, which looked like the model losing the Canary Current
# when it was the diagnostic mixing two signals. Offline re-binning of the saved sections found it.
#
# The deep bound also cannot use the literature's number: σθ maxes at 27.798 in this domain (coldest,
# saltiest corner T=1.97 S=34.944), so γⁿ=27.922 is unreachable and that bin would always read zero.
const σ_THERMOCLINE  = 27.10     # ≈ 500 m here; the Canary Current layer
const σ_INTERMEDIATE = 27.70     # ≈ 1800 m; NADW
@inline σθ_of(T, S) = SWP.ρ(T, S, 0.0, eos) - 1000.0
@inline density_layer(σ) = σ < σ_THERMOCLINE ? 1 : (σ < σ_INTERMEDIATE ? 2 : 3)

function layered_boundary_transports(T, S, u, v)
    Q = zeros(3, 3)                     # [boundary W/S/N, layer thermocline/intermediate/deep]
    @inbounds for k in 1:Nz, j in 1:Ny
        wetcell(1, j, k, grid_c, Face(), Center(), Center()) || continue
        Q[1, density_layer(σθ_of(T[1, j, k], S[1, j, k]))] += u[1, j, k] * Axᶠᶜᶜ(1, j, k, grid_c)
    end
    @inbounds for k in 1:Nz, i in 1:Nx
        if wetcell(i, 1, k, grid_c, Center(), Face(), Center())
            Q[2, density_layer(σθ_of(T[i, 1, k], S[i, 1, k]))] += v[i, 1, k] * Ayᶜᶠᶜ(i, 1, k, grid_c)
        end
        if wetcell(i, Ny+1, k, grid_c, Center(), Face(), Center())
            Q[3, density_layer(σθ_of(T[i, Ny, k], S[i, Ny, k]))] -= v[i, Ny+1, k] * Ayᶜᶠᶜ(i, Ny+1, k, grid_c)
        end
    end
    return Q
end

function static_instability(T, S)
    n = 0; worst = 0.0
    @inbounds for j in 1:Ny, i in 1:Nx, k in 1:(Nz-1)
        (wetcell(i, j, k,   grid_c, Center(), Center(), Center()) &&
         wetcell(i, j, k+1, grid_c, Center(), Center(), Center())) || continue
        zr = 0.5 * (z_centers[k] + z_centers[k+1])
        Δρ = SWP.ρ(T[i, j, k+1], S[i, j, k+1], zr, eos) - SWP.ρ(T[i, j, k], S[i, j, k], zr, eos)
        if Δρ > 1e-4; n += 1; worst = max(worst, Δρ); end
    end
    return n, worst
end

function surface_freshwater_flux(m)
    use_atmos || return 0.0
    Jη = arr(m.interfaces.net_fluxes.ocean.η)
    Q = 0.0
    @inbounds for j in 1:Ny, i in 1:Nx
        immersed_inactive_node(i, j, Nz, grid_c, Center(), Center(), Center()) && continue
        Q += Jη[i, j, 1] * Azᶜᶜᶜ(i, j, Nz, grid_c)
    end
    return Q
end

function native_range(fts, floor_value)
    a = filter(x -> x > floor_value, vec(Array(interior(fts[1]))))
    isempty(a) ? (NaN, NaN) : (minimum(a), maximum(a))
end
T_ref = native_range(T_fts, 1.0); S_ref = native_range(S_fts, 5.0)
@info @sprintf("Reference (GLORYS native): T ∈ [%.3f, %.3f] | S ∈ [%.3f, %.3f]", T_ref[1], T_ref[2], S_ref[1], S_ref[2])
let (lo,_,hi,_) = wet_extrema(arr(ocean.model.tracers.T)), (slo,_,shi,_) = wet_extrema(arr(ocean.model.tracers.S))
    @info @sprintf("Initial condition (wet):    T ∈ [%.3f, %.3f] | S ∈ [%.3f, %.3f]", lo, hi, slo, shi)
end
let (W, S_, N) = glorys_residuals(Time(0.0))
    @info @sprintf("GLORYS residuals at t=0 (+ = inflow): west %+.3e  south %+.3e  north %+.3e  net %+.3e m³/s; NOT corrected",
                   W, S_, N, W + S_ + N)
end

# Calibrate the proxy: report where the σθ surfaces sit in depth, so the offset from the literature's
# γⁿ surfaces is visible rather than assumed. The papers put γⁿ=26.85 at roughly 80-320 m; if our
# σθ=27.38 surface sits much deeper than their ~700 m thermocline base, our "thermocline" layer is
# including intermediate water and will read high against theirs.
let T = arr(ocean.model.tracers.T), S = arr(ocean.model.tracers.S)
    for σ in (26.85, σ_THERMOCLINE, σ_INTERMEDIATE)
        d = Float64[]
        for j in 1:Ny, i in 1:Nx
            k = findlast(k -> wet3d[i, j, k] && σθ_of(T[i, j, k], S[i, j, k]) >= σ, 1:Nz)
            k === nothing && continue
            push!(d, -z_centers[k])
        end
        isempty(d) ? @info(@sprintf("σθ = %.2f surface: absent from this domain", σ)) :
            @info @sprintf("σθ = %.2f surface sits at %.0f m (median over wet columns; range %.0f–%.0f)",
                           σ, median(d), minimum(d), maximum(d))
    end
end

const salt₀ = Ref(0.0); const vol₀ = Ref(0.0)
let Ta0 = arr(ocean.model.tracers.T), Sa0 = arr(ocean.model.tracers.S)
    SV, V = salt_and_volume(Sa0)
    n, w = static_instability(Ta0, Sa0)
    Tc, To, ΔT, Ts = upwelling_index(Ta0)
    @info @sprintf("At rest: salt %.6e | volume %.6e | mean S %.4f | unstable pairs %d (worst Δρ %.3f)",
                   SV, V, SV/V, n, w)
    @info @sprintf("At rest: SST coastal %.3f, offshore %.3f, ΔT %+.3f °C | T at %.0f m over shelf %.3f °C",
                   Tc, To, ΔT, z_centers[k_src], Ts)
end

const wall_t0 = Ref(Base.time()); const wall_last = Ref(Base.time()); const sim_last = Ref(0.0)

function progress(sim)
    model = ocean_model_of(sim.model)
    sync_grid_c!()
    # ONE device->host copy per field per block, then everything below is plain array work. This is
    # what makes the script GPU-safe: the diagnostics all index cell by cell, which a device array
    # will not allow. ~17 MB per field at 0.02°, every 10 iterations -- negligible against the step.
    ua, va, wa = arr(model.velocities.u), arr(model.velocities.v), arr(model.velocities.w)
    Ta, Sa = arr(model.tracers.T), arr(model.tracers.S)
    ηa = arr(free_surface_field(model.free_surface))
    Tmin, _, Tmax, _ = wet_extrema(Ta)
    Smin, Imin, Smax, Imax = wet_extrema(Sa)
    _, Iw = findmax(abs, wa); _, Iu = findmax(abs, ua)

    @info @sprintf("iter %d | t = %s | Δt = %s | max|u,v,w| = (%.2e, %.2e, %.2e) | max|η| = %.3e | T ∈ [%.3f, %.3f] | S ∈ [%.3f, %.3f] (excursion %+.3f / %+.3f)",
                   iteration(sim), prettytime(time(sim)), prettytime(sim.Δt),
                   maximum(abs, ua), maximum(abs, va), maximum(abs, wa), maximum(abs, ηa),
                   Tmin, Tmax, Smin, Smax, Smin - S_ref[1], Smax - S_ref[2])
    @info @sprintf("        where: Smax %s | Smin %s | max|w| %s | max|u| %s",
                   whereis(Imax), whereis(Imin), whereis(Iw), whereis(Iu))

    W, So, N = boundary_transports(ua, va)
    Wr, Sr, Nr = glorys_residuals(Time(model.clock.time))
    @info @sprintf("        transport (+ = inflow): west %+.3e  south %+.3e  north %+.3e | net %+.3e m³/s -> %+.3e m/hour",
                   W, So, N, W+So+N, (W+So+N) / wet_surface_area * 3600)
    @info @sprintf("        GLORYS residual (uncorrected): west %+.3e  south %+.3e  north %+.3e | net %+.3e m³/s",
                   Wr, Sr, Nr, Wr+Sr+Nr)

    Q = layered_boundary_transports(Ta, Sa, ua, va)
    for (li, lab) in ((1, "Canary Cur σθ<27.10"), (2, "undercurrent  <27.70"), (3, "deep      σθ≥27.70"))
        @info @sprintf("        %s: west %+7.3f  south %+7.3f  north %+7.3f | net %+7.3f Sv",
                       lab, Q[1,li]/1e6, Q[2,li]/1e6, Q[3,li]/1e6, sum(Q[:,li])/1e6)
    end

    Tc, To, ΔT, Ts = upwelling_index(Ta)
    @info @sprintf("        UPWELLING: SST coastal %.3f  offshore %.3f  ΔT %+.3f °C | T at %.0f m over shelf %.3f °C | ΔT(fixed-λ, Jan-comparable) %+.3f",
                   Tc, To, ΔT, z_centers[k_src], Ts, upwelling_index_fixedλ(Ta))

    if use_tracer
        ca, aca = arr(model.tracers.c), arr(model.tracers.ac)
        C, AC, cmax, imax = tracer_budget(ca, aca)
        Wt, St, Nt = tracer_export(ca, ua, va)
        elapsed = max(time(sim), eps())
        @info @sprintf("        TRACER: inventory %.4e of %.4e released (%.2f%%) | mean age %s | max c %.3e %s",
                       C, elapsed, 100C/elapsed,
                       C > 1e-6*elapsed ? @sprintf("%.3f d", AC/C/86400) : "n/a (too little tracer)",
                       cmax, whereis(imax))
        @info @sprintf("        TRACER export (+ = in): west %+.3e  south %+.3e  north %+.3e | net %+.3e /s",
                       Wt, St, Nt, Wt+St+Nt)
    end

    SV, V = salt_and_volume(Sa)
    ninv, worstΔρ = static_instability(Ta, Sa)
    if vol₀[] == 0.0 && iteration(sim) > 0
        salt₀[], vol₀[] = SV, V
        @info @sprintf("        drift baseline set at t = %s (σ reconciled)", prettytime(time(sim)))
    end
    if use_atmos
        Qfw = surface_freshwater_flux(sim.model)
        @info @sprintf("        surface freshwater: %+.4e m³/s (%+.3f m/yr basin mean)",
                       Qfw, Qfw / wet_surface_area * 3.15576e7)
    end
    have = vol₀[] > 0
    @info @sprintf("        conservation: salt %s | volume %s | mean S %.4f | unstable pairs %d (worst Δρ %.3f)",
                   have ? @sprintf("%+.3e", (SV-salt₀[])/salt₀[]) : "(pending)",
                   have ? @sprintf("%+.3e", (V-vol₀[])/vol₀[])   : "(pending)",
                   SV/V, ninv, worstΔρ)

    now = Base.time(); wall = now - wall_t0[]
    Δwall, Δsim = now - wall_last[], time(sim) - sim_last[]
    wall_last[], sim_last[] = now, time(sim)
    rate = Δwall > 0 ? Δsim/Δwall : 0.0
    eta = (isfinite(sim.stop_time) && rate > 0) ? (sim.stop_time - time(sim))/rate : NaN
    @info @sprintf("        pacing: wall %s | %.0f× realtime | %.2f s/iter | ETA %s",
                   prettytime(wall), rate, Δwall/10, isnan(eta) ? "n/a" : prettytime(eta))
    flush(stderr)
    return nothing
end
add_callback!(simulation, progress, IterationInterval(10))

# --- Spin-up then the run -----------------------------------------------------
simulation.Δt = 10
simulation.stop_iteration = 100
run!(simulation)

ocean_timescale(model) = cell_advection_timescale(ocean_model_of(model))
wizard = TimeStepWizard(; cfl = 0.15, max_Δt = 15minutes, max_change = 1.1, cell_advection_timescale = ocean_timescale)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))
simulation.stop_iteration = Inf
simulation.stop_time = parse(Float64, get(ENV, "CANARY_STOP_DAYS", string(run_days))) * days

output_fields   = merge(ocean.model.velocities, ocean.model.tracers)
output_schedule = TimeInterval(6hours)
tag = (weno_order == 7 ? "" : "_weno$(weno_order)") * (use_zstar ? "_zstar" : "") *
      (use_atmos ? (jra_multiyear ? "_jra$(Dates.year(run_start))" : "_jraryf") : "") *
      (Δdeg == 0.02 ? "" : "_d$(replace(string(Δdeg), "." => ""))") * get(ENV, "CANARY_TAG", "")
out(name) = joinpath(@__DIR__, "canary_glorys_$(name)$(tag)")

simulation.output_writers[:surface] = JLD2Writer(ocean.model, output_fields;
    schedule = output_schedule, overwrite_existing = true, indices = (:, :, Nz), filename = out("surface"))
simulation.output_writers[:west_section] = JLD2Writer(ocean.model, output_fields;
    schedule = output_schedule, overwrite_existing = true, indices = (1, :, :), filename = out("west_section"))
# Zonal section at 26.5N, cutting the shelf, the upwelling front and the open ocean.
const j_shelf = argmin(abs.(φ_centers .- 26.5))
simulation.output_writers[:shelf_section] = JLD2Writer(ocean.model, output_fields;
    schedule = output_schedule, overwrite_existing = true, indices = (:, j_shelf, :), filename = out("shelf_section"))
# All three open boundaries are saved, so the layered transports can be RE-BINNED offline: the layer
# definition (σθ as a γⁿ proxy, and where to put the thresholds) is the least settled part of this
# diagnostic, and re-deriving it from file is free where re-running is not. T, S and the grid all
# come back from these files, so σθ and area-weighted transport are fully recoverable.
#
# The north needs TWO writers. Face-in-y indexing puts the domain edge at Ny+1 for v while T and S
# live at Ny, and one slice cannot hold both -- `(:, Ny, :)` would silently give the last INTERIOR
# face, which is the off-by-one this script had until it was caught. West and south are exempt:
# index 1 is both the boundary face and the first interior cell.
simulation.output_writers[:south_section] = JLD2Writer(ocean.model, output_fields;
    schedule = output_schedule, overwrite_existing = true, indices = (:, 1, :), filename = out("south_section"))
simulation.output_writers[:north_section] = JLD2Writer(ocean.model, output_fields;
    schedule = output_schedule, overwrite_existing = true, indices = (:, Ny, :), filename = out("north_section"))
simulation.output_writers[:north_face_v] = JLD2Writer(ocean.model, (; v = ocean.model.velocities.v);
    schedule = output_schedule, overwrite_existing = true, indices = (:, Ny+1, :), filename = out("north_face_v"))

# Daily 3-D snapshot (~83 MB/frame, ~0.6 GB for 7 days). Cheap against the cost of re-running to
# answer a question this run could already have answered: interior sections, water-mass
# transformation, volume-integrated budgets.
simulation.output_writers[:snapshot3d] = JLD2Writer(ocean.model, output_fields;
    schedule = TimeInterval(24hours), overwrite_existing = true, filename = out("snapshot3d"))  # 24hours, not 1day: Oceananigans.Units and Dates both export `day`

simulation.output_writers[:free_surface] = JLD2Writer(ocean.model, (; η = free_surface_field(ocean.model.free_surface));
    schedule = output_schedule, overwrite_existing = true, filename = out("free_surface"))
simulation.output_writers[:checkpointer] = Checkpointer(ocean.model;
    schedule = TimeInterval(12hours), dir = @__DIR__, prefix = "canary_checkpoint$(tag)",
    overwrite_existing = true, cleanup = true)

run!(simulation)
