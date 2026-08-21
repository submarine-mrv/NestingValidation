using CairoMakie
using Oceananigans
using NumericalEarth
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
#u: west  open boundary for Juan de Fuca   (Radiation / Orlanski on 3-D velocity)
#v: north open boundary for Georgia Strait (Radiation / Orlanski on 3-D velocity)
#U/V: barotropic transport -> GravityWaveRadiation (Flather 1976) for split-explicit solver
#η: free surface           -> ImplicitGravityWaveRadiation (Chapman 1985), pairs with U/V
#T/S: west + north Value boundaries with Radiation relaxation
#
# NOTE: this branch renamed the Flather/Chapman conditions. On ss/open-boundary-conditions:
#   Flather (barotropic transport)  == GravityWaveRadiationBoundaryCondition
#   Chapman (free surface)          == ImplicitGravityWaveRadiationBoundaryCondition
# The old `FlatherBoundaryCondition` / `ChapmanBoundaryCondition` names no longer exist.
#
# STAGE 1 (this file): constant scalar inflow, no external data.
# STAGE 2 (later): replace the constants below with GLORYS fields via DataWrangling.

using Oceananigans.BoundaryConditions: Radiation,
                                       GravityWaveRadiationBoundaryCondition,
                                       ImplicitGravityWaveRadiationBoundaryCondition
using Oceananigans.Operators: Δzᶠᶜᶜ, Δzᶜᶠᶜ
using Oceananigans.ImmersedBoundaries: immersed_inactive_node, immersed_peripheral_node

# --- placeholder constants (tune these) ---
inflow_speed = 0.1     # m/s, prescribed normal inflow magnitude
T_inflow     = 8.0     # °C, placeholder inflow temperature
S_inflow     = 31.0    # psu, placeholder inflow salinity

# Sign convention (normal velocity is +x eastward / +y northward):
#   west  boundary inflow (into domain) => +u
#   north boundary inflow (into domain) => -v (flowing southward)
u_west_inflow  = +inflow_speed
v_north_inflow = -inflow_speed

# Orlanski radiation scheme: relax toward the prescribed exterior value on inflow,
# radiate freely on outflow. inflow_timescale small => strong nudging to the exterior.
radiation = Radiation(inflow_timescale = 6hours, outflow_timescale = Inf)

# helper: true only for genuinely wet cells (not immersed / not peripheral)
@inline wetcell(i, j, k, grid, ℓx, ℓy, ℓz) = !immersed_peripheral_node(i, j, k, grid, ℓx, ℓy, ℓz) & !immersed_inactive_node(i, j, k, grid, ℓx, ℓy, ℓz)

# --- 3-D velocity: Radiation / Orlanski on the boundary-normal component ---
u_bcs = FieldBoundaryConditions(west  = NormalFlowBoundaryCondition(u_west_inflow;  scheme = radiation))
v_bcs = FieldBoundaryConditions(north = NormalFlowBoundaryCondition(v_north_inflow; scheme = radiation))

# --- Barotropic transport: GravityWaveRadiation (Flather) ---
# The barotropic transport U/V is a 2-D field (Face,Center,nothing)/(Center,Face,nothing).
# Its halo fill calls getbc(bc, tangential_idx, k, grid, clock, fields), so we integrate the
# constant inflow velocity over the wet water column here to get the exterior transport.
# The condition must return the 2-tuple (transport, exterior free-surface η).

@inline function west_U_obc(j, k, grid, clock, fields)
    i = 1                      # west boundary face (x-Face index 1)
    U = zero(grid)
    @inbounds for kk in 1:grid.Nz
        wet = wetcell(i, j, kk, grid, Face(), Center(), Center())
        U += ifelse(wet, u_west_inflow * Δzᶠᶜᶜ(i, j, kk, grid), zero(grid))
    end
    return (U, zero(grid))     # (exterior transport, exterior η = 0)
end

@inline function north_V_obc(i, k, grid, clock, fields)
    j = grid.Ny + 1            # north boundary face (y-Face index Ny+1)
    V = zero(grid)
    @inbounds for kk in 1:grid.Nz
        wet = wetcell(i, j, kk, grid, Center(), Face(), Center())
        V += ifelse(wet, v_north_inflow * Δzᶜᶠᶜ(i, j, kk, grid), zero(grid))
    end
    return (V, zero(grid))     # (exterior transport, exterior η = 0)
end

U_bcs = FieldBoundaryConditions(grid, (Face(), Center(), nothing);
                                west  = GravityWaveRadiationBoundaryCondition(west_U_obc; discrete_form = true))
V_bcs = FieldBoundaryConditions(grid, (Center(), Face(), nothing);
                                north = GravityWaveRadiationBoundaryCondition(north_V_obc; discrete_form = true))

# --- Free surface: ImplicitGravityWaveRadiation (Chapman), pairs with the U/V transport bcs ---
η_bcs = FieldBoundaryConditions(grid, (Center(), Center(), Face());
                                west  = ImplicitGravityWaveRadiationBoundaryCondition(),
                                north = ImplicitGravityWaveRadiationBoundaryCondition())

# --- T/S: constant exterior value with Radiation relaxation ---
T_bcs = FieldBoundaryConditions(west  = ValueBoundaryCondition(T_inflow; scheme = radiation),
                                north = ValueBoundaryCondition(T_inflow; scheme = radiation))
S_bcs = FieldBoundaryConditions(west  = ValueBoundaryCondition(S_inflow; scheme = radiation),
                                north = ValueBoundaryCondition(S_inflow; scheme = radiation))

# --- Sponge (flow-relaxation) layer at the open boundaries ---
# Complements the OBCs: inside a band next to each open boundary, relax the normal
# velocity toward the prescribed inflow to absorb reflections. This replaces the
# original Fu/Fv sponge forcings. Widths are in degrees (grid is lon/lat).
sponge_width = 0.2                 # ~20 km band inward from each open boundary
sponge_rate  = 1 / 1hours          # max relaxation rate at the boundary (1/s), tunable

# Masks: 1 at the boundary, linearly ->0 over sponge_width inward.
west_mask  = PiecewiseLinearMask{:x}(center = λ₁, width = sponge_width)  # west open boundary (min longitude)
north_mask = PiecewiseLinearMask{:y}(center = φ₂, width = sponge_width)  # north open boundary (max latitude)

u_sponge = Relaxation(rate = sponge_rate, mask = west_mask,  target = u_west_inflow)
v_sponge = Relaxation(rate = sponge_rate, mask = north_mask, target = v_north_inflow)

#define ocean simulation
# ocean_simulation appends these u/v forcings to its own barotropic-potential + radiative
# forcing automatically (it merges user forcing per-field).
ocean = ocean_simulation(grid;
                        tracers = (:T, :S),   # CATKE closure adds :e automatically; listing it explicitly errors
                        forcing = (u = u_sponge, v = v_sponge),
                        boundary_conditions = (; u = u_bcs, U = U_bcs, v = v_bcs, V = V_bcs, η = η_bcs, T = T_bcs, S = S_bcs),
)


# STAGE 1 (validated): run the bare `ocean` simulation -- it is already a Simulation
# returned by ocean_simulation -- driven purely by the open boundaries, no surface forcing.
#
# STAGE 2 (adds JRA55 surface fluxes; untested here, downloads data):
#   atmosphere    = JRA55PrescribedAtmosphere(arch)
#   coupled_model = OceanOnlyModel(ocean; atmosphere)
#   simulation    = Simulation(coupled_model)   # then run! this instead of `ocean`

# --- Initial conditions ---
# STAGE 1: uniform T/S matching the prescribed inflow (no external data).
set!(ocean.model, T = T_inflow, S = S_inflow)
#
# STAGE 2: initialize from ECCO/GLORYS instead, e.g.
#   physics_bc_dataset = ECCO4Monthly()
#   date = DateTime(1993, 1, 1)
#   set!(ocean.model, T = Metadatum(:temperature; date, dataset = physics_bc_dataset),
#                     S = Metadatum(:salinity;    date, dataset = physics_bc_dataset))

# --- Progress monitoring ---
function progress(sim)
    u, v, w = sim.model.velocities
    T, S = sim.model.tracers.T, sim.model.tracers.S
    @info @sprintf("iter %d | t = %s | Δt = %s | max|u,v,w| = (%.2e, %.2e, %.2e) | T ∈ [%.3f, %.3f] | S ∈ [%.3f, %.3f]",
                   iteration(sim), prettytime(time(sim)), prettytime(sim.Δt),
                   maximum(abs, u), maximum(abs, v), maximum(abs, w),
                   minimum(T), maximum(T), minimum(S), maximum(S))
    flush(stderr)   # force output through Julia's block buffer when stderr is redirected
    return nothing
end
add_callback!(ocean, progress, IterationInterval(10))

# --- Spin-up: short run at a small fixed timestep (this is what has been validated) ---
ocean.Δt = 2
ocean.stop_iteration = 100
run!(ocean)

# --- Full run: adaptive timestep out to 7 days with daily surface output ---
# NOTE: only the spin-up above has been tested. The wizard + output writer + long run
# below are configured but not yet validated; uncomment the final run! to launch it.
ocean_timescale(model) = cell_advection_timescale(model)
wizard = TimeStepWizard(; cfl = 0.2, max_Δt = 20minutes, max_change = 1.3, cell_advection_timescale = ocean_timescale)
ocean.callbacks[:wizard] = Callback(wizard, IterationInterval(10))
ocean.stop_iteration = Inf
ocean.stop_time = 7days

ocean.output_writers[:surface_daily] = JLD2Writer(ocean.model, merge(ocean.model.velocities, ocean.model.tracers);
                                                   schedule = TimeInterval(1days),
                                                   overwrite_existing = true,
                                                   indices = (:, :, Nz), # surface level
                                                   filename = "salishsea_physics_openboundaries_test1")
# run!(ocean)   # <-- uncomment to launch the full 7-day run
