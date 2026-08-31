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
using NumericalEarth
using CopernicusMarine
using Oceananigans.Advection: cell_advection_timescale
using Oceananigans.Units
using CFTime
using Dates
using Printf

const DW  = NumericalEarth.DataWrangling
const SWP = Oceananigans.BuoyancyFormulations.SeawaterPolynomials

arch = CPU()

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

# Fixed subset so resolution changes do not invalidate the download cache.
glorys_region = DW.BoundingBox(longitude = (-22.0, -9.5), latitude = (24.5, 30.0))
@assert λ₁ >= -22.0 && λ₂ <= -9.5 && φ₁ >= 24.5 && φ₂ <= 30.0 "domain outside glorys_region"
glorys = DW.GLORYSDaily()

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
fill_down_fts!(fts, floor) = (for n in 1:length(fts.times); fill_down!(interior(fts[n]), floor); end; fts)

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

# --- Ocean simulation ---------------------------------------------------------
weno_order = parse(Int, get(ENV, "CANARY_WENO_ORDER", "7"))
tracer_advection = WENO(order = weno_order,
                        time_discretization = AdaptiveVerticallyImplicitDiscretization(cfl = 0.5))

ocean = ocean_simulation(grid;
                         tracers = (:T, :S),
                         tracer_advection,
                         forcing = (T = T_restoring, S = S_restoring,
                                    u = u_restoring, v = v_restoring),
                         boundary_conditions = (; u = u_bcs, U = U_bcs, v = v_bcs, V = V_bcs,
                                                  η = η_bcs, T = T_bcs, S = S_bcs))

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
let ηi = interior(η_field), t₀ = Time(0.0)
    for i in 1:Nx, j in 1:Ny
        wet = !immersed_inactive_node(i, j, Nz, grid, Center(), Center(), Center())
        ηi[i, j, 1] = ifelse(wet, zos_at(i, j, grid, t₀), 0.0)
    end
end
@info @sprintf("η from GLORYS zos: [%.4f, %.4f] m", minimum(interior(η_field)), maximum(interior(η_field)))

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

wet3d = [!immersed_inactive_node(i, j, k, grid, Center(), Center(), Center()) for i in 1:Nx, j in 1:Ny, k in 1:Nz]
wet_idx = findall(wet3d)
@info "Wet cells: $(length(wet_idx)) of $(Nx*Ny*Nz) ($(round(100*length(wet_idx)/(Nx*Ny*Nz), digits=1))%)"
wet_surface_area = sum(Azᶜᶜᶜ(i, j, Nz, grid) for i in 1:Nx, j in 1:Ny
                       if !immersed_inactive_node(i, j, Nz, grid, Center(), Center(), Center()); init = 0.0)

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

function wet_extrema(field)
    a = interior(field)[wet_idx]
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
function boundary_transports(model)
    u, v = model.velocities.u, model.velocities.v
    W = S_ = N = 0.0
    @inbounds for k in 1:Nz, j in 1:Ny
        wetcell(1, j, k, grid, Face(), Center(), Center()) && (W += u[1, j, k] * Axᶠᶜᶜ(1, j, k, grid))
    end
    @inbounds for k in 1:Nz, i in 1:Nx
        wetcell(i, 1, k, grid, Center(), Face(), Center())  && (S_ += v[i, 1, k]  * Ayᶜᶠᶜ(i, 1, k, grid))
        wetcell(i, Ny, k, grid, Center(), Face(), Center()) && (N  -= v[i, Ny, k] * Ayᶜᶠᶜ(i, Ny, k, grid))
    end
    return W, S_, N
end

# What the raw GLORYS exterior would carry through each edge: a parent/child mismatch measure, not a
# physical signal. Reported, never corrected.
function glorys_residuals(t)
    W = S_ = N = 0.0
    @inbounds for k in 1:Nz, j in 1:Ny
        wetcell(1, j, k, grid, Face(), Center(), Center()) &&
            (W += DW.interpolate_to_grid(u_fts, 1, j, k, u_native, grid, t) * Axᶠᶜᶜ(1, j, k, grid))
    end
    @inbounds for k in 1:Nz, i in 1:Nx
        wetcell(i, 1, k, grid, Center(), Face(), Center()) &&
            (S_ += DW.interpolate_to_grid(v_fts, i, 1, k, v_native, grid, t) * Ayᶜᶠᶜ(i, 1, k, grid))
        wetcell(i, Ny, k, grid, Center(), Face(), Center()) &&
            (N -= DW.interpolate_to_grid(v_fts, i, Ny, k, v_native, grid, t) * Ayᶜᶠᶜ(i, Ny, k, grid))
    end
    return W, S_, N
end

# THE signal for this domain. Upwelling is offshore Ekman transport drawing cold water up against the
# African coast, so the diagnostic is the cross-shore SST contrast in the eastern shelf strip against
# the open ocean at the same latitude, plus the mean temperature at the upwelling source depth.
const k_src = argmin(abs.(z_centers .+ 100))     # ~100 m, the source depth for coastal upwelling
function upwelling_index(model)
    T = model.tracers.T
    coastalT = offshoreT = 0.0; nc = no = 0
    srcT = 0.0; ns = 0
    @inbounds for j in 1:Ny, i in 1:Nx
        wet3d[i, j, Nz] || continue
        λi = λ_centers[i]
        if λi > -13.5                       # shelf strip against the coast
            coastalT += T[i, j, Nz]; nc += 1
        elseif -20.5 < λi < -19.0           # open-ocean reference, clear of the islands
            offshoreT += T[i, j, Nz]; no += 1
        end
        if λi > -13.5 && wet3d[i, j, k_src]
            srcT += T[i, j, k_src]; ns += 1
        end
    end
    Tc = nc > 0 ? coastalT/nc : NaN
    To = no > 0 ? offshoreT/no : NaN
    return Tc, To, Tc - To, ns > 0 ? srcT/ns : NaN
end

# --- Conservation and stability monitors --------------------------------------
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

const eos = SWP.TEOS10.TEOS10EquationOfState()
function static_instability(model)
    T, S = model.tracers.T, model.tracers.S
    n = 0; worst = 0.0
    @inbounds for j in 1:Ny, i in 1:Nx, k in 1:(Nz-1)
        (wetcell(i, j, k,   grid, Center(), Center(), Center()) &&
         wetcell(i, j, k+1, grid, Center(), Center(), Center())) || continue
        zr = 0.5 * (z_centers[k] + z_centers[k+1])
        Δρ = SWP.ρ(T[i, j, k+1], S[i, j, k+1], zr, eos) - SWP.ρ(T[i, j, k], S[i, j, k], zr, eos)
        if Δρ > 1e-4; n += 1; worst = max(worst, Δρ); end
    end
    return n, worst
end

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

function native_range(fts, floor_value)
    a = filter(x -> x > floor_value, vec(interior(fts[1])))
    isempty(a) ? (NaN, NaN) : (minimum(a), maximum(a))
end
T_ref = native_range(T_fts, 1.0); S_ref = native_range(S_fts, 5.0)
@info @sprintf("Reference (GLORYS native): T ∈ [%.3f, %.3f] | S ∈ [%.3f, %.3f]", T_ref[1], T_ref[2], S_ref[1], S_ref[2])
let (lo,_,hi,_) = wet_extrema(ocean.model.tracers.T), (slo,_,shi,_) = wet_extrema(ocean.model.tracers.S)
    @info @sprintf("Initial condition (wet):    T ∈ [%.3f, %.3f] | S ∈ [%.3f, %.3f]", lo, hi, slo, shi)
end
let (W, S_, N) = glorys_residuals(Time(0.0))
    @info @sprintf("GLORYS residuals at t=0 (+ = inflow): west %+.3e  south %+.3e  north %+.3e  net %+.3e m³/s; NOT corrected",
                   W, S_, N, W + S_ + N)
end

const salt₀ = Ref(0.0); const vol₀ = Ref(0.0)
let (SV, V) = salt_and_volume(ocean.model)
    n, w = static_instability(ocean.model)
    Tc, To, ΔT, Ts = upwelling_index(ocean.model)
    @info @sprintf("At rest: salt %.6e | volume %.6e | mean S %.4f | unstable pairs %d (worst Δρ %.3f)",
                   SV, V, SV/V, n, w)
    @info @sprintf("At rest: SST coastal %.3f, offshore %.3f, ΔT %+.3f °C | T at %.0f m over shelf %.3f °C",
                   Tc, To, ΔT, z_centers[k_src], Ts)
end

const wall_t0 = Ref(Base.time()); const wall_last = Ref(Base.time()); const sim_last = Ref(0.0)

function progress(sim)
    model = ocean_model_of(sim.model)
    u, v, w = model.velocities
    T, S = model.tracers.T, model.tracers.S
    η = free_surface_field(model.free_surface)
    Tmin, _, Tmax, _ = wet_extrema(T)
    Smin, Imin, Smax, Imax = wet_extrema(S)
    _, Iw = findmax(abs, interior(w)); _, Iu = findmax(abs, interior(u))

    @info @sprintf("iter %d | t = %s | Δt = %s | max|u,v,w| = (%.2e, %.2e, %.2e) | max|η| = %.3e | T ∈ [%.3f, %.3f] | S ∈ [%.3f, %.3f] (excursion %+.3f / %+.3f)",
                   iteration(sim), prettytime(time(sim)), prettytime(sim.Δt),
                   maximum(abs, u), maximum(abs, v), maximum(abs, w), maximum(abs, η),
                   Tmin, Tmax, Smin, Smax, Smin - S_ref[1], Smax - S_ref[2])
    @info @sprintf("        where: Smax %s | Smin %s | max|w| %s | max|u| %s",
                   whereis(Imax), whereis(Imin), whereis(Iw), whereis(Iu))

    W, So, N = boundary_transports(model)
    Wr, Sr, Nr = glorys_residuals(Time(model.clock.time))
    @info @sprintf("        transport (+ = inflow): west %+.3e  south %+.3e  north %+.3e | net %+.3e m³/s -> %+.3e m/hour",
                   W, So, N, W+So+N, (W+So+N) / wet_surface_area * 3600)
    @info @sprintf("        GLORYS residual (uncorrected): west %+.3e  south %+.3e  north %+.3e | net %+.3e m³/s",
                   Wr, Sr, Nr, Wr+Sr+Nr)

    Tc, To, ΔT, Ts = upwelling_index(model)
    @info @sprintf("        UPWELLING: SST coastal %.3f  offshore %.3f  ΔT %+.3f °C | T at %.0f m over shelf %.3f °C",
                   Tc, To, ΔT, z_centers[k_src], Ts)

    SV, V = salt_and_volume(model)
    ninv, worstΔρ = static_instability(model)
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
simulation.output_writers[:free_surface] = JLD2Writer(ocean.model, (; η = free_surface_field(ocean.model.free_surface));
    schedule = output_schedule, overwrite_existing = true, filename = out("free_surface"))
simulation.output_writers[:checkpointer] = Checkpointer(ocean.model;
    schedule = TimeInterval(12hours), dir = @__DIR__, prefix = "canary_checkpoint$(tag)",
    overwrite_existing = true, cleanup = true)

run!(simulation)
