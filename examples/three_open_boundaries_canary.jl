# # Three open boundaries: Canary upwelling, with GPU-compilable boundary kernels
#
# The Canary upwelling region off NW Africa. West, south and north are open; the east edge is the
# African coast, and the archipelago sits inside the domain. Read `one_open_boundary_redsea.jl`
# first — it is the same pattern with one boundary and simpler plumbing. This example adds:
#
# - several boundaries at once, and the crossed tangential wiring that implies;
# - per-side scheme selection, because the southern boundary terminates in coast;
# - the parent precomputed onto the boundary faces, which is what lets the condition kernels
#   compile for a GPU: sampling a `FieldTimeSeries` inside a kernel works on the CPU but sends the
#   CUDA compiler's inliner past 31 GB. Every array the kernels touch reaches them through
#   `parameters`, so nothing is captured from global scope.
#
# `examples/README.md` covers the five conditions an open side needs and the traps in the parent
# data; `docs/obc-implementation-recipe.md` has the experiments behind the choices.

using NumericalEarth
using Oceananigans
using Oceananigans.Units
using Oceananigans.Units: Time                        # not exported by Units
using CopernicusMarine
using Oceananigans.Architectures: on_architecture
using Oceananigans.BoundaryConditions: ObliqueRadiation, TracerReservoir,
                                       SurfaceWaveRadiationBoundaryCondition
using Oceananigans.Operators: Δzᶠᶜᶜ, Δzᶜᶠᶜ
using Oceananigans.ImmersedBoundaries: immersed_inactive_node, immersed_peripheral_node
using Dates, Printf

DW = NumericalEarth.DataWrangling

# ## Configuration

architecture = CPU()
start_date = DateTime(1993, 1, 1)
run_duration = 7days
boundary_dates = start_date : Day(1) : (start_date + Day(9))

λ_west, λ_east = -21.5, -10.0
φ_south, φ_north = 25.0, 29.5
Δλ = Δφ = 0.02
Nx = round(Int, (λ_east - λ_west) / Δλ)
Ny = round(Int, (φ_north - φ_south) / Δφ)
Nz = 32

depth = 4800.0   # GLORYS carries no valid data below 4836.6 m here
z = MutableVerticalDiscretization(ExponentialDiscretization(Nz, -depth, 0; scale = 700))

grid = LatitudeLongitudeGrid(architecture; size = (Nx, Ny, Nz),
                             longitude = (λ_west, λ_east), latitude = (φ_south, φ_north),
                             z, halo = (7, 7, 7),
                             topology = (Bounded, Bounded, Bounded))

## Light smoothing keeps the island wakes this domain exists to resolve.
bottom_height = regrid_bathymetry(grid; minimum_depth = 10, interpolation_passes = 5)
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bottom_height))

## The precompute below walks the grid cell by cell, so it needs a host copy.
host_grid = on_architecture(CPU(), grid)

@inline wet(i, j, k, grid, ℓx, ℓy, ℓz) =
    !immersed_peripheral_node(i, j, k, grid, ℓx, ℓy, ℓz) &
    !immersed_inactive_node(i, j, k, grid, ℓx, ℓy, ℓz)

# ## Parent data

glorys = DW.GLORYSDaily()
glorys_region = DW.BoundingBox(longitude = (-22.0, -9.5), latitude = (24.5, 30.0))
parent_metadata(name) = DW.Metadata(name; dataset = glorys, dates = boundary_dates, region = glorys_region)

sponge_width = 0.25
@inline edge(distance) = exp(-(distance / sponge_width)^2)
## max(), not a sum, so the corners do not double-count.
@inline sponge_mask(λ, φ, z, t) = max(edge(λ - λ_west), edge(φ - φ_south), edge(φ_north - φ))

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

inpainting = DW.NearestNeighborInpainting(50)
restoring_kw = (; rate = 1/12hours, mask = sponge_mask, inpainting,
                  time_indices_in_memory = length(boundary_dates),
                  time_indexing = Oceananigans.OutputReaders.Clamp())

T_restoring = DW.DatasetRestoring(parent_metadata(:temperature), architecture; restoring_kw...)
S_restoring = DW.DatasetRestoring(parent_metadata(:salinity),    architecture; restoring_kw...)
u_restoring = DW.DatasetRestoring(parent_metadata(:u_velocity),  architecture; restoring_kw...)
v_restoring = DW.DatasetRestoring(parent_metadata(:v_velocity),  architecture; restoring_kw...)

fill_down_series!(T_restoring.field_time_series, 1.0)
fill_down_series!(S_restoring.field_time_series, 5.0)

# ## Precomputing the parent onto the faces
#
# Each face array is `(along-boundary index, depth, parent time)`; the kernels read it by index and
# interpolate linearly in time.

function host_series(device_series, name, floor)
    architecture isa CPU && return device_series
    series = FieldTimeSeries(parent_metadata(name), CPU();
                             time_indices_in_memory = length(boundary_dates), inpainting)
    floor === nothing || fill_down_series!(series, floor)
    return series
end

T_series = host_series(T_restoring.field_time_series, :temperature, 1.0)
S_series = host_series(S_restoring.field_time_series, :salinity,    5.0)
u_series = host_series(u_restoring.field_time_series, :u_velocity,  nothing)
v_series = host_series(v_restoring.field_time_series, :v_velocity,  nothing)

zos_series = FieldTimeSeries(parent_metadata(:free_surface), CPU();
                             time_indices_in_memory = length(boundary_dates),
                             time_indexing = Oceananigans.OutputReaders.Clamp(), inpainting)

parent_times = Float64.(collect(u_series.times))
parent_Δt = parent_times[2] - parent_times[1]
@assert all(diff(parent_times) .≈ parent_Δt) "the kernels below assume uniformly spaced parent times"

## Carry the parent's origin: GLORYSDaily is stamped at the centre of its daily-mean window.
parent_t₀ = parent_times[1]
parent_Nt = length(parent_times)

function face_values(series, side, j_index)
    FT = eltype(interior(series[1]))
    n_along = side === :west ? Ny : Nx
    values = zeros(FT, n_along, Nz, parent_Nt)
    for n in 1:parent_Nt, k in 1:Nz, a in 1:n_along
        i, j = side === :west ? (1, a) : side === :south ? (a, 1) : (a, j_index)
        values[a, k, n] = DW.interpolate_to_grid(series, i, j, k, series.grid, host_grid, Time(parent_times[n]))
    end
    return on_architecture(architecture, values)
end

function face_free_surface(side)
    n_along = side === :west ? Ny : Nx
    values = zeros(Float64, n_along, parent_Nt)
    for n in 1:parent_Nt, a in 1:n_along
        i, j = side === :west ? (1, a) : side === :south ? (a, 1) : (a, Ny)
        values[a, n] = DW.interpolate_to_grid(zos_series, i, j, 1, zos_series.grid, host_grid, Time(parent_times[n]))
    end
    return on_architecture(architecture, values)
end

## Normal velocities at the boundary face index; tracers, tangential velocities and η at the
## adjacent interior node, half a cell inside.
west = (; u = face_values(u_series, :west, 1),  tangential = face_values(v_series, :west, 1),
          T = face_values(T_series, :west, 1),  S = face_values(S_series, :west, 1),
          η = face_free_surface(:west), Δt = parent_Δt, Nt = parent_Nt, t₀ = parent_t₀)
south = (; v = face_values(v_series, :south, 1), tangential = face_values(u_series, :south, 1),
           T = face_values(T_series, :south, 1), S = face_values(S_series, :south, 1),
           η = face_free_surface(:south), Δt = parent_Δt, Nt = parent_Nt, t₀ = parent_t₀)
north = (; v = face_values(v_series, :north, Ny + 1), tangential = face_values(u_series, :north, Ny),
           T = face_values(T_series, :north, Ny),     S = face_values(S_series, :north, Ny),
           η = face_free_surface(:north), Δt = parent_Δt, Nt = parent_Nt, t₀ = parent_t₀)

@inline function parent_index(t, p)
    τ = (t - p.t₀) / p.Δt
    n = clamp(unsafe_trunc(Int, τ) + 1, 1, p.Nt - 1)
    return n, clamp(τ - (n - 1), zero(τ), one(τ))
end

@inline function at_face(values, a, k, t, p)
    n, w = parent_index(t, p)
    @inbounds return (1 - w) * values[a, k, n] + w * values[a, k, n + 1]
end

@inline function at_face(values, a, t, p)
    n, w = parent_index(t, p)
    @inbounds return (1 - w) * values[a, n] + w * values[a, n + 1]
end

# ## Boundary conditions

@inline u_west(j, k, grid, clock, fields, p) = at_face(p.u, j, k, clock.time, p)
@inline v_south(i, k, grid, clock, fields, p) = at_face(p.v, i, k, clock.time, p)
@inline v_north(i, k, grid, clock, fields, p) = at_face(p.v, i, k, clock.time, p)
@inline tangential(a, k, grid, clock, fields, p) = at_face(p.tangential, a, k, clock.time, p)
@inline T_exterior(a, k, grid, clock, fields, p) = at_face(p.T, a, k, clock.time, p)
@inline S_exterior(a, k, grid, clock, fields, p) = at_face(p.S, a, k, clock.time, p)

## Flather pairs the depth-integrated exterior transport with the parent's free surface.
@inline function U_west(j, k, grid, clock, fields, p)
    U = zero(grid)
    @inbounds for k in 1:grid.Nz
        uk = at_face(p.u, j, k, clock.time, p)
        U += ifelse(wet(1, j, k, grid, Face(), Center(), Center()), uk * Δzᶠᶜᶜ(1, j, k, grid), zero(grid))
    end
    return (U, at_face(p.η, j, clock.time, p))
end

@inline function meridional_transport(i, grid, clock, p, j)
    V = zero(grid)
    @inbounds for k in 1:grid.Nz
        vk = at_face(p.v, i, k, clock.time, p)
        V += ifelse(wet(i, j, k, grid, Center(), Face(), Center()), vk * Δzᶜᶠᶜ(i, j, k, grid), zero(grid))
    end
    return (V, at_face(p.η, i, clock.time, p))
end

@inline V_south(i, k, grid, clock, fields, p) = meridional_transport(i, grid, clock, p, 1)
@inline V_north(i, k, grid, clock, fields, p) = meridional_transport(i, grid, clock, p, grid.Ny)

## Oblique diagnoses an along-boundary phase speed, so it needs a side that does not end in coast.
## West is fully wet and north 99%; south is 58% wet and terminates in the African coast.
velocity_scheme(side) = side === :south ?
    NormalRadiation(inflow_timescale = 6hours, outflow_timescale = Inf) :
    ObliqueRadiation(inflow_timescale = 6hours, outflow_timescale = Inf)

tracer_scheme = TracerReservoir(inflow_length_scale = 5000, outflow_length_scale = 0)

normal_flow(f, p, side) = NormalFlowBoundaryCondition(f; scheme = velocity_scheme(side), discrete_form = true, parameters = p)
along_flow(f, p, side) = ValueBoundaryCondition(f; scheme = velocity_scheme(side), discrete_form = true, parameters = p)
tracer(f, p) = ValueBoundaryCondition(f; scheme = tracer_scheme, discrete_form = true, parameters = p)
transport(f, p) = GravityWaveRadiationBoundaryCondition(f; discrete_form = true, parameters = p)

## With perpendicular boundaries the normal component on one side is the tangential component on
## another, so u and v each carry conditions on all three sides.
u_bcs = FieldBoundaryConditions(west = normal_flow(u_west, west, :west),
                                south = along_flow(tangential, south, :south),
                                north = along_flow(tangential, north, :north))
v_bcs = FieldBoundaryConditions(south = normal_flow(v_south, south, :south),
                                north = normal_flow(v_north, north, :north),
                                west = along_flow(tangential, west, :west))
U_bcs = FieldBoundaryConditions(grid, (Face(), Center(), nothing); west = transport(U_west, west))
V_bcs = FieldBoundaryConditions(grid, (Center(), Face(), nothing); south = transport(V_south, south),
                                                                   north = transport(V_north, north))
η_bcs = FieldBoundaryConditions(grid, (Center(), Center(), Face());
                                west = SurfaceWaveRadiationBoundaryCondition(),
                                south = SurfaceWaveRadiationBoundaryCondition(),
                                north = SurfaceWaveRadiationBoundaryCondition())
T_bcs = FieldBoundaryConditions(west = tracer(T_exterior, west), south = tracer(T_exterior, south),
                                north = tracer(T_exterior, north))
S_bcs = FieldBoundaryConditions(west = tracer(S_exterior, west), south = tracer(S_exterior, south),
                                north = tracer(S_exterior, north))

ocean = ocean_simulation(grid;
                         tracers = (:T, :S),
                         tracer_advection = WENO(order = 7,
                             time_discretization = AdaptiveVerticallyImplicitDiscretization(cfl = 0.5)),
                         forcing = (T = T_restoring, S = S_restoring,
                                    u = u_restoring, v = v_restoring),
                         boundary_conditions = (; u = u_bcs, v = v_bcs, U = U_bcs, V = V_bcs,
                                                  η = η_bcs, T = T_bcs, S = S_bcs))

# ## Initial conditions
#
# Drawn from the same repaired series the boundaries use, through a scratch field.

function interpolate_parent!(field, series, t)
    values = interior(field)
    @inbounds for k in 1:Nz, j in 1:Ny, i in 1:Nx
        values[i, j, k] = DW.interpolate_to_grid(series, i, j, k, series.grid, host_grid, t)
    end
    return field
end

T₀ = interpolate_parent!(CenterField(host_grid), T_series, Time(0))
S₀ = interpolate_parent!(CenterField(host_grid), S_series, Time(0))
interior(ocean.model.tracers.T) .= on_architecture(architecture, interior(T₀))
interior(ocean.model.tracers.S) .= on_architecture(architecture, interior(S₀))

free_surface = ocean.model.free_surface
η = hasproperty(free_surface, :η) ? free_surface.η : free_surface.displacement
let values = zeros(Nx, Ny)
    for i in 1:Nx, j in 1:Ny
        surface_wet = !immersed_inactive_node(i, j, Nz, host_grid, Center(), Center(), Center())
        values[i, j] = ifelse(surface_wet,
                              DW.interpolate_to_grid(zos_series, i, j, 1, zos_series.grid, host_grid, Time(0)), 0.0)
    end
    interior(η)[:, :, 1] .= on_architecture(architecture, values)
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

ocean.output_writers[:surface] =
    JLD2Writer(ocean.model, merge(ocean.model.velocities, ocean.model.tracers);
               schedule = TimeInterval(6hours), indices = (:, :, Nz),
               filename = "canary_surface", overwrite_files = true)

ocean.Δt = 10
ocean.stop_iteration = 100
run!(ocean)

ocean.callbacks[:wizard] = Callback(TimeStepWizard(; cfl = 0.15, max_Δt = 15minutes, max_change = 1.1),
                                    IterationInterval(10))
ocean.stop_iteration = Inf
ocean.stop_time = run_duration
run!(ocean)
