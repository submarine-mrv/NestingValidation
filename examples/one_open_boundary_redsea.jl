# # One open boundary: the Red Sea driven by GLORYS
#
# The Red Sea and Gulf of Aden. South, north and west are coast; the east edge at 50°E is the only
# open boundary, and Bab el-Mandeb sits inside the domain so the model resolves it.
#
# Each open side carries five conditions: the boundary-normal and boundary-tangential velocities,
# the barotropic transport, the free surface, and the tracers. `examples/README.md` explains what
# each one carries and the traps in the parent data; `docs/obc-implementation-recipe.md` has the
# experiments behind the choices.
#
# `three_open_boundaries_canary.jl` adds multiple boundaries, per-side scheme selection, and the
# face precomputation that makes the condition kernels GPU-compilable.

using NumericalEarth
using Oceananigans
using Oceananigans.Units
using Oceananigans.Units: Time                        # not exported by Units
using CopernicusMarine                                   # activates the GLORYS download backend
using Oceananigans.BoundaryConditions: ObliqueRadiation, TracerReservoir,
                                       SurfaceWaveRadiationBoundaryCondition
using Oceananigans.Operators: Δzᶠᶜᶜ
using Oceananigans.ImmersedBoundaries: immersed_inactive_node, immersed_peripheral_node
using Dates, Printf

DW = NumericalEarth.DataWrangling

# ## Configuration

architecture = CPU()
start_date = DateTime(1993, 1, 1)
run_duration = 7days

boundary_dates = start_date : Day(1) : (start_date + Day(9))

λ_west, λ_east = 32.0, 50.0
φ_south, φ_north = 10.0, 30.0
Δλ = Δφ = 0.05
Nx = round(Int, (λ_east - λ_west) / Δλ)
Ny = round(Int, (φ_north - φ_south) / Δφ)
Nz = 40

depth = 3200   # GLORYS carries no valid data below 3226 m here
z = MutableVerticalDiscretization(ExponentialDiscretization(Nz, -depth, 0; scale = 1800))

grid = LatitudeLongitudeGrid(architecture; size = (Nx, Ny, Nz),
                             longitude = (λ_west, λ_east), latitude = (φ_south, φ_north),
                             z, halo = (7, 7, 7),
                             topology = (Bounded, Bounded, Bounded))

bottom_height = regrid_bathymetry(grid; minimum_depth = 10, interpolation_passes = 5)
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height))

@inline wet(i, j, k, grid, ℓx, ℓy, ℓz) =
    !immersed_peripheral_node(i, j, k, grid, ℓx, ℓy, ℓz) &
    !immersed_inactive_node(i, j, k, grid, ℓx, ℓy, ℓz)

# ## Parent data

glorys = DW.GLORYSDaily()
glorys_region = DW.BoundingBox(longitude = (31.5, 50.5), latitude = (9.5, 30.5))
parent_metadata(name) = DW.Metadata(name; dataset = glorys, dates = boundary_dates, region = glorys_region)

sponge_width = 0.5
@inline sponge_mask(λ, φ, z, t) = exp(-((λ_east - λ) / sponge_width)^2)

function fill_down!(field, floor)
    nx, ny, nz = size(field)
    @inbounds for j in 1:ny, i in 1:nx
        deepest_valid = NaN
        for k in nz:-1:1
            value = field[i, j, k]
            value > floor ? (deepest_valid = value) :
                            (isfinite(deepest_valid) && (field[i, j, k] = deepest_valid))
        end
    end
    return field
end

fill_down_series!(series, floor) =
    (for n in 1:length(series.times); fill_down!(interior(series[n]), floor); end; series)

restoring_kw = (; rate = 1/12hours,
                  mask = sponge_mask,
                  inpainting = DW.NearestNeighborInpainting(50),
                  time_indices_in_memory = length(boundary_dates),
                  time_indexing = Oceananigans.OutputReaders.Clamp())

T_restoring = DW.DatasetRestoring(parent_metadata(:temperature), architecture; restoring_kw...)
S_restoring = DW.DatasetRestoring(parent_metadata(:salinity),    architecture; restoring_kw...)
u_restoring = DW.DatasetRestoring(parent_metadata(:u_velocity),  architecture; restoring_kw...)
v_restoring = DW.DatasetRestoring(parent_metadata(:v_velocity),  architecture; restoring_kw...)

fill_down_series!(T_restoring.field_time_series, 1.0)
fill_down_series!(S_restoring.field_time_series, 5.0)

T_series, T_grid = T_restoring.field_time_series, T_restoring.native_grid
S_series, S_grid = S_restoring.field_time_series, S_restoring.native_grid
u_series, u_grid = u_restoring.field_time_series, u_restoring.native_grid
v_series, v_grid = v_restoring.field_time_series, v_restoring.native_grid

zos_series = FieldTimeSeries(parent_metadata(:free_surface), architecture;
                             time_indices_in_memory = length(boundary_dates),
                             time_indexing = Oceananigans.OutputReaders.Clamp(),
                             inpainting = DW.NearestNeighborInpainting(50))
zos_grid = zos_series.grid
@inline zos(i, j, grid, t) = DW.interpolate_to_grid(zos_series, i, j, 1, zos_grid, grid, t)

# ## Boundary conditions

velocity_scheme = ObliqueRadiation(inflow_timescale = 6hours, outflow_timescale = Inf)
tracer_scheme = TracerReservoir(inflow_length_scale = 5000, outflow_length_scale = 0)

@inline u_east(j, k, grid, clock, fields) =
    DW.interpolate_to_grid(u_series, grid.Nx, j, k, u_grid, grid, Time(clock.time))
@inline v_east(j, k, grid, clock, fields) =
    DW.interpolate_to_grid(v_series, grid.Nx, j, k, v_grid, grid, Time(clock.time))
@inline T_east(j, k, grid, clock, fields) =
    DW.interpolate_to_grid(T_series, grid.Nx, j, k, T_grid, grid, Time(clock.time))
@inline S_east(j, k, grid, clock, fields) =
    DW.interpolate_to_grid(S_series, grid.Nx, j, k, S_grid, grid, Time(clock.time))

## Flather pairs the depth-integrated exterior transport with the parent's free surface.
@inline function U_east(j, k, grid, clock, fields)
    i = grid.Nx
    t = Time(clock.time)
    U = zero(grid)
    @inbounds for k in 1:grid.Nz
        uk = DW.interpolate_to_grid(u_series, i, j, k, u_grid, grid, t)
        U += ifelse(wet(i, j, k, grid, Face(), Center(), Center()), uk * Δzᶠᶜᶜ(i, j, k, grid), zero(grid))
    end
    return (U, zos(i, j, grid, t))
end

u_bcs = FieldBoundaryConditions(east = NormalFlowBoundaryCondition(u_east; scheme = velocity_scheme, discrete_form = true))
v_bcs = FieldBoundaryConditions(east = ValueBoundaryCondition(v_east; scheme = velocity_scheme, discrete_form = true))
U_bcs = FieldBoundaryConditions(east = GravityWaveRadiationBoundaryCondition(U_east; discrete_form = true))
T_bcs = FieldBoundaryConditions(east = ValueBoundaryCondition(T_east; scheme = tracer_scheme, discrete_form = true))
S_bcs = FieldBoundaryConditions(east = ValueBoundaryCondition(S_east; scheme = tracer_scheme, discrete_form = true))
η_bcs = FieldBoundaryConditions(grid, (Center(), Center(), Face());
                                east = SurfaceWaveRadiationBoundaryCondition())

ocean = ocean_simulation(grid;
                         tracers = (:T, :S),
                         tracer_advection = WENO(order = 7,
                             time_discretization = AdaptiveVerticallyImplicitDiscretization(cfl = 0.5)),
                         forcing = (T = T_restoring, S = S_restoring,
                                    u = u_restoring, v = v_restoring),
                         boundary_conditions = (; u = u_bcs, v = v_bcs, U = U_bcs,
                                                  η = η_bcs, T = T_bcs, S = S_bcs))

# ## Initial conditions
#
# Drawn from the same repaired series the boundaries use, through a scratch field.

function interpolate_parent!(field, series, native_grid, t)
    values = interior(field)
    @inbounds for k in 1:Nz, j in 1:Ny, i in 1:Nx
        values[i, j, k] = DW.interpolate_to_grid(series, i, j, k, native_grid, grid, t)
    end
    return field
end

T₀ = interpolate_parent!(CenterField(grid), T_series, T_grid, Time(0))
S₀ = interpolate_parent!(CenterField(grid), S_series, S_grid, Time(0))
interior(ocean.model.tracers.T) .= interior(T₀)
interior(ocean.model.tracers.S) .= interior(S₀)

free_surface = ocean.model.free_surface
η = hasproperty(free_surface, :η) ? free_surface.η : free_surface.displacement
let values = interior(η)
    for i in 1:Nx, j in 1:Ny
        surface_wet = !immersed_inactive_node(i, j, Nz, grid, Center(), Center(), Center())
        values[i, j, 1] = ifelse(surface_wet, zos(i, j, grid, Time(0)), 0.0)
    end
end

# ## Run

function progress(sim)
    u, v, w = sim.model.velocities
    T, S = sim.model.tracers
    @info @sprintf("iter %d | t = %s | Δt = %s | max|u, v, w| = (%.2e, %.2e, %.2e) | T ∈ [%.2f, %.2f] | S ∈ [%.2f, %.2f]",
                   iteration(sim), prettytime(time(sim)), prettytime(sim.Δt),
                   maximum(abs, u), maximum(abs, v), maximum(abs, w),
                   minimum(T), maximum(T), minimum(S), maximum(S))
    return nothing
end

add_callback!(ocean, progress, IterationInterval(50))

ocean.output_writers[:east_section] =
    JLD2Writer(ocean.model, merge(ocean.model.velocities, ocean.model.tracers);
               schedule = TimeInterval(6hours), indices = (Nx, :, :),
               filename = "redsea_east_section", overwrite_files = true)

ocean.Δt = 10
ocean.stop_iteration = 100
run!(ocean)

ocean.callbacks[:wizard] = Callback(TimeStepWizard(; cfl = 0.15, max_Δt = 15minutes, max_change = 1.1),
                                    IterationInterval(10))
ocean.stop_iteration = Inf
ocean.stop_time = run_duration
run!(ocean)
