# truncation_truth.jl -- the "truth" run for the open-boundary truncation experiment.
#
# A doubly-periodic 2-D turbulence field. Periodic means truth has NO open boundaries anywhere, so
# every point in it is uncontaminated -- including the ring where the child's boundaries will sit.
# That is the whole reason for this design: differences in the child can only come from its open
# boundaries, because nothing else differs.
#
# The child covers truth cells CHILD_I in both directions, on the SAME spacing, so boundary values
# transfer by index rather than interpolation and introduce no error of their own.
#
# Usage: julia --project=obc_v113 obc_v113/truncation_truth.jl

using Oceananigans, Printf, Statistics

const N      = 256                 # truth cells per side
const L      = 1.0
const CHILD_I = 65:192             # truth cells the child covers (128 per side)
const OUT    = joinpath(@__DIR__, "truth.jld2")

grid = RectilinearGrid(CPU(); size = (N, N), x = (0, L), y = (0, L),
                       topology = (Periodic, Periodic, Flat))

# A mean flow past an obstacle placed INSIDE the child's footprint, plus a turbulent background.
#
# The obstacle is the point. A child that is merely initialised from truth and driven by truth has no
# structure of its own, so the best strategy is to trust the boundary and a hard Dirichlet wins --
# measured, with the first version of this script: dirichlet 8.8e-3 against normal 2.4e-2, and
# degrading the parent 8x in space and 10x in time barely moved it (9.4e-3 vs 2.5e-2). Radiation
# earns its keep only when the child carries features the parent CANNOT represent and those features
# have to leave. The wake supplies exactly that: truth resolves it, the coarsened parent does not.
struct Cylinder; x0::Float64; y0::Float64; r::Float64; end
@inline (c::Cylinder)(x, y) = ifelse((x-c.x0)^2 + (y-c.y0)^2 < c.r^2, 1.0, 0.0)
const CYL = Cylinder(0.40, 0.50, 0.035)          # inside the child box, upstream of its east face

model = NonhydrostaticModel(grid; advection = WENO(order = 5), closure = ScalarDiffusivity(ν = 1e-5),
                            forcing = (u = Relaxation(rate = 50.0, mask = CYL, target = 0.0),
                                       v = Relaxation(rate = 50.0, mask = CYL, target = 0.0)))

using Random; Random.seed!(42)
const UBAR = 1.0        # mean flow that carries the wake out through the child's east boundary
uᵢ(x, y) = sum(cos(2π*k*x/L + 0.7k) * sin(2π*m*y/L + 1.3m) / (k^2+m^2) for k in 1:6, m in 1:6)
vᵢ(x, y) = sum(sin(2π*k*x/L + 1.1k) * cos(2π*m*y/L + 0.3m) / (k^2+m^2) for k in 1:6, m in 1:6)
set!(model, u = uᵢ, v = vᵢ)
# normalise to a sensible speed
U0 = 1.0
umax = maximum(abs, interior(model.velocities.u))
set!(model, u = (x,y) -> UBAR + 0.25*uᵢ(x,y)*U0/umax, v = (x,y) -> 0.25*vᵢ(x,y)*U0/umax)

Δt = 0.15 * (L/N) / (UBAR + U0)
simulation = Simulation(model; Δt, stop_time = 4.0)

# Spin-up first so the field is genuinely turbulent before the child ever starts.
const SPINUP = 1.0
@info @sprintf("truth: %d^2, Δ=%.5f, Δt=%.5f, spin-up to t=%.1f then record to t=%.1f",
               N, L/N, Δt, SPINUP, simulation.stop_time)

u, v = model.velocities
simulation.output_writers[:fields] = JLD2Writer(model, (; u, v);
    filename = OUT, schedule = TimeInterval(0.02), overwrite_files = true)

prog(sim) = @info @sprintf("  t=%.3f  max|u|=%.4f", time(sim), maximum(abs, interior(sim.model.velocities.u)))
simulation.callbacks[:p] = Callback(prog, TimeInterval(0.5))

run!(simulation)
@info "wrote $OUT"
