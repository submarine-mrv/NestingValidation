# oblique_tracer_outflow.jl -- does ObliqueRadiation beat NormalRadiation for a tracer leaving an
# open boundary at an angle?
#
# This is the canary question in miniature. The release plume leaves the southern boundary on a
# slant, and the tracer boundary condition there is ValueBoundaryCondition(...; scheme = radiation).
# NormalRadiation only diagnoses the boundary-NORMAL phase speed, so an obliquely-advected feature
# is radiated with the wrong speed; ObliqueRadiation (Raymond & Kuo 1984) adds the tangential term.
#
# The flow is PRESCRIBED and uniform, so truth is analytic: a Gaussian blob translating at
# (U cos θ, U sin θ), wrapping in the periodic y. Any departure is a boundary artefact, so there is
# no reference run to argue about. Two numbers per configuration:
#   rms   -- RMS error against the analytic solution while the blob straddles the open boundary
#   resid -- tracer mass left in the domain after the blob should have fully exited (ideal 0)
#
# NOTE: θ = 0 is NOT a control. The schemes coincide when the tangential GRADIENT OF THE FIELD
# vanishes, not when the tangential velocity does -- a blob has dc/dy != 0 at the boundary even in
# purely normal flow. The exact-agreement control is Oceananigans' own unit test
# ("ObliqueRadiation reduces to NormalRadiation at normal incidence"), which passes.
#
# A hard Dirichlet arm sets the scale: if radiation-vs-Dirichlet dwarfs Oblique-vs-Normal, then the
# choice between the two radiation schemes is second order for this configuration.
#
# Usage: julia --project=obc_v113 obc_v113/oblique_tracer_outflow.jl

using Oceananigans
using Oceananigans.BoundaryConditions: NormalRadiation, ObliqueRadiation
using Printf, Statistics

const N  = (96, 64)
const U  = 1.0
const σ  = 0.08                      # blob width
const x₀, y₀ = 0.35, 0.5             # start well inside, upstream of the east boundary

grid = RectilinearGrid(CPU(); size = (N[1], N[2], 1), x = (0, 1), y = (0, 1), z = (0, 1),
                       topology = (Bounded, Periodic, Bounded))

xc = xnodes(grid, Center()); yc = ynodes(grid, Center())

# Analytic solution: rigid translation, periodic in y.
function analytic(t, θ)
    u, v = U*cos(θ), U*sin(θ)
    [exp(-(((x - (x₀ + u*t))^2 + (mod(y - (y₀ + v*t) + 0.5, 1.0) - 0.5)^2) / (2σ^2)))
     for x in xc, y in yc]
end

# Timing follows the angle: the blob must cross the EAST boundary, and the normal speed is
# U cos θ, so a fixed step count would leave the steep angles never reaching it at all.
function run_case(scheme, θ; Δt = 2e-3)
    u, v = U*cos(θ), U*sin(θ)
    t_cross = (1.0 - x₀) / u                 # centre reaches the east boundary
    sample  = round(Int, t_cross / Δt)       # straddling it
    nsteps  = round(Int, (t_cross + 4σ/u) / Δt)  # well past: blob fully outside
    bc = isnothing(scheme) ? ValueBoundaryCondition(0) : ValueBoundaryCondition(0; scheme)
    model = HydrostaticFreeSurfaceModel(grid;
                velocities = PrescribedVelocityFields(u = u, v = v),
                tracer_advection = Centered(), buoyancy = nothing, tracers = :c,
                boundary_conditions = (; c = FieldBoundaryConditions(east = bc, west = bc)))
    set!(model, c = (x, y, z) -> exp(-(((x-x₀)^2 + (y-y₀)^2) / (2σ^2))))

    rms = NaN
    for n in 1:nsteps
        time_step!(model, Δt)
        if n == sample                                   # blob straddling the east boundary
            c = Array(interior(model.tracers.c, :, :, 1))
            rms = sqrt(mean((c .- analytic(n*Δt, θ)).^2))
        end
        _ = n
    end
    c = Array(interior(model.tracers.c, :, :, 1))
    return rms, sum(abs, c) / (N[1]*N[2])                # residual after full exit
end

nrm = NormalRadiation(inflow_timescale = 0, outflow_timescale = Inf)
obl = ObliqueRadiation(inflow_timescale = 0, outflow_timescale = Inf)

@printf("%6s | %-11s %-11s %-11s | %-11s %-11s %-11s | %s\n", "theta",
        "rms Dirich", "rms Normal", "rms Obliq", "res Dirich", "res Normal", "res Obliq",
        "Obl vs Nrm (resid)")
println(repeat("-", 112))
for θdeg in (0, 15, 30, 45, 60)
    θ = deg2rad(θdeg)
    (rd, dd) = run_case(nothing, θ)
    (rn, dn) = run_case(nrm, θ)
    (ro, do_) = run_case(obl, θ)
    @printf("%5d° | %-11.3e %-11.3e %-11.3e | %-11.3e %-11.3e %-11.3e | %+.0f%%\n",
            θdeg, rd, rn, ro, dd, dn, do_, 100*(do_-dn)/dn)
end
