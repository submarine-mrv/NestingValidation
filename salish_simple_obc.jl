using CairoMakie
using Oceananigans
using Oceananigans: architecture
using Oceananigans.Advection: cell_advection_timescale
using Oceananigans.Units
using CFTime
using Dates
using Printf

arch = CPU()

start_date = DateTime(1993, 1, 1)
end_date   = DateTime(1993, 12, 31)
dates = start_date: Month(1): end_date

#define grid
λ₁, λ₂  = (235.5,238) # domain in longitude
φ₁, φ₂  = (47, 49.75) # domain in latitude

Nz = 36
depth = 1000
scale_factor = 0.85
z = ExponentialDiscretization(Nz, -depth, 0; scale = scale_factor*depth)

Nx = 250 #lower resolution for testing
Ny = 275
#Nx = 500 #start with 1/80 degree or ~1.25 km resolution, work toward 200 m (1/200 deg is 500 m, 1/400 is 350 m) ?
#Ny = 550


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

grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bottom_height))

#define boundary conditions
#u: west  open boundary for Juan de Fuca
#v: north open boundary for Georgia Strait
#T/S: west + north GLORYS value boundaries

#first define input dataset fields
physics_bc_dataset = ECCO4Monthly()
temperature = Metadata(:temperature; dates, dataset = physics_bc_dataset)
salinity = Metadata(:salinity; dates, dataset = physics_bc_dataset)
u_velocity = Metadata(:u_velocity; dates, dataset = physics_bc_dataset)
v_velocity = Metadata(:v_velocity; dates, dataset = physics_bc_dataset)


# --- 3-D GLORYS 'NormalFlow' Radiation boundary conditions --
u_west_bc = NormalFlowBoundaryCondition(forcing; scheme = Radiation(), parameters = )
u_north_bc = NormalFlowBoundaryCondition(forcing; scheme = Radiation(), parameters = )
u_bcs     = FieldBoundaryConditions(west = u_west_bc, north = u_north_bc)

# ---- Barotropic Flather boundary conditions ---

using Oceananigans.ImmersedBoundaries: immersed_inactive_node, immersed_peripheral_node

@inline wetcell(i, j, k, grid, ℓx, ℓy, ℓz) = !immersed_peripheral_node(i, j, k, grid, ℓx, ℓy, ℓz) & !immersed_inactive_node(i, j, k, grid, ℓx, ℓy, ℓz)

@inline function west_U_obc(j, k, grid, clock, fields)
    i = 0
    U = zero(grid)
    @inbounds for k in 1:grid.Nz
        wet = wetcell(i, j, k, grid, Face(), Center(), Center())
        U += ifelse(wet, 0.05 * Oceananigans.Operators.Δzᶠᶜᶜ(i, j, k, grid), zero(U))
    end
    return (U, 0.0)
end



U_west_bc = FlatherBoundaryCondition(U_west_glorys;
                                     discrete_form = true,
                                     parameters = glorys_params)

V_north_bc = FlatherBoundaryCondition(V_north_glorys;
                                      discrete_form = true,
                                      parameters = glorys_params)

U_bcs = FieldBoundaryConditions(grid, (Face(), Center(), nothing); west = FlatherBoundaryCondition(west_U_obc, discrete_form = true))
V_bcs = FieldBoundaryConditions(grid, (Face(), Center(), nothing); north = FlatherBoundaryCondition(north_U_obc, discrete_form = true))


#T/S boundary conditions


T = FieldBoundaryConditions(
        west  = ValueBoundaryCondition(T_west_glorys;  scheme = Radiation(), parameters = glorys_params),
        north = ValueBoundaryCondition(T_north_glorys; scheme = Radiation(), parameters = glorys_params),
    )
S = FieldBoundaryConditions(
        west  = ValueBoundaryCondition(S_west_glorys;  scheme = Radiation(), parameters = glorys_params),
        north = ValueBoundaryCondition(S_north_glorys; scheme = Radiation(), parameters = glorys_params),
    )

#define ocean, atmosphere, and coupled model
ocean = ocean_simulation(grid;
                        forcing = (T = FT, S = FS, u = Fu, v = Fv),
                        tracers = (:T, :S, :e),
                        boundary_conditions = (; u = u_bcs, U = U_bcs, v = v_bcs, V = V_bcs, T = T_bcs, S = S_bcs),
)


#atmosphere
atmosphere = NumericalEarth.JRA55PrescribedAtmosphere(arch)

#define the coupled model
coupled_model = NumericalEarth.OceanOnlyModel(ocean; atmosphere)

#define the simulation
simulation = Simulation(coupled_model)

#Initial conditions
date = DateTime(1993, 1, 1)
ecco_temperature = Metadatum(:temperature; date, physics_bc_dataset)
ecco_salinity = Metadatum(:salinity; date, physics_bc_dataset)
set!(ocean.model, T=ecco_temperature,
                  S=ecco_salinity,
                  )

function progress(sim)
    u, v, w = ocean.model.velocities
    T, S, e = ocean.model.tracers

    @info @sprintf("Time: %s, Iteration %d, Δt %s, max(vel): (%.2e, %.2e, %.2e), min(T, S): %+.4f, %+.4f, max(T, S): %+.4f, %+.4f,  \n",
                   prettytime(simulation.model.clock.time),
                   simulation.model.clock.iteration,
                   prettytime(simulation.Δt),
                   maximum(abs, u), maximum(abs, v), maximum(abs, w),
                   minimum(T), minimum(S),
                   maximum(T), maximum(S))
end

simulation.callbacks[:progress] = Callback(progress, IterationInterval(1))


#run spin up
simulation.Δt = 2
simulation.stop_iteration = 100
run!(simulation)

#define full run

ocean_timescale(model) = cell_advection_timescale(coupled_model.ocean.model)

wizard = TimeStepWizard(; cfl = 0.2, max_Δt = 10minutes, max_change = 1.1, cell_advection_timescale = ocean_timescale)

simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))
simulation.stop_iteration = Inf
simulation.stop_time = 7days

ocean.output_writers[:surface_daily] = NetCDFWriter(model, merge(model.velocities, model.tracers);
                                                         schedule = TimeInterval(1days),
                                                         overwrite_existing = true,
                                                         indices = (:, :, Nz), #Nz is the surface level
                                                         filename = "salishsea_physics_openboundaries_test1")
