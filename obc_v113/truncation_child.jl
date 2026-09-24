# truncation_child.jl -- the child half of the truncation experiment.
#
# Runs a sub-box of the truth domain with OPEN boundaries fed from truth, and measures how far the
# child drifts from truth in its own interior. Truth is doubly periodic, so it has no boundary error
# anywhere: every difference the child shows is caused by its open boundaries and nothing else.
#
# Arms compared (see obc_v113/README.md):
#   dirichlet  -- hard Value/NormalFlow, no radiation scheme: the "do nothing clever" control
#   normal     -- NormalRadiation on the boundary-NORMAL components only (what canary/ does today)
#   normal+tan -- NormalRadiation on normal AND TANGENTIAL components
#   oblique    -- ObliqueRadiation on normal components only
#   oblique+tan-- ObliqueRadiation on normal and tangential
#
# The +tan arms matter because Oceananigans PR #5962 attributes its headline gain to radiating the
# tangential velocity, and canary/ currently sets no tangential open boundary condition at all.
#
# THE PARENT MUST BE DEGRADED, or the experiment is vacuous. With exact truth on the child's own grid
# at output cadence, a hard Dirichlet boundary is optimal by construction and every radiation scheme
# loses to it -- measured, before this degradation existed: dirichlet 8.8e-3 against normal 2.4e-2.
# Radiation exists to let interior features leave when the exterior is WRONG. So the parent here is
# coarsened in space by PARENT_R and in time by PARENT_K, mimicking canary/'s real mismatch: GLORYS
# at 1/12 deg (~4.5x the child's 0.02 deg) sampled daily.
#
# Usage: julia --project=obc_v113 obc_v113/truncation_child.jl

using Oceananigans, Printf, Statistics
using Oceananigans.BoundaryConditions: NormalRadiation, ObliqueRadiation, NormalFlowBoundaryCondition

const N, L    = 256, 1.0
const OFF     = 64                    # child cell/face index + OFF = truth index
const Nc      = 128                   # child cells per side
const Δ       = L / N
const SPINUP  = 1.0
const PARENT_R = parse(Int, get(ENV, "PARENT_R", "8"))   # parent coarser than child by this factor
const PARENT_K = parse(Int, get(ENV, "PARENT_K", "10"))  # parent sampled every K-th truth frame
const TRUTH   = joinpath(@__DIR__, "truth.jld2")

ut = FieldTimeSeries(TRUTH, "u"); vt = FieldTimeSeries(TRUTH, "v")
times = Float64.(collect(ut.times))
n0 = findfirst(t -> t >= SPINUP - 1e-9, times)
const TT   = times[n0:end]                     # truth times the child will use
const NT   = length(TT)
const T0   = TT[1]
const DT   = TT[2] - TT[1]
@assert all(diff(TT) .≈ DT) "truth output must be uniformly spaced"

U(n) = Array(interior(ut[n0 + n - 1]))[:, :, 1]
V(n) = Array(interior(vt[n0 + n - 1]))[:, :, 1]

# --- boundary arrays, by exact index mapping (grids are aligned, so no interpolation error) -------
# normal components sit ON the child's boundary faces; tangential components sit half a cell off it,
# so those are the average of the two truth cells straddling the boundary.
UW = zeros(Nc, NT); UE = zeros(Nc, NT); VS = zeros(Nc, NT); VN = zeros(Nc, NT)   # normal
VW = zeros(Nc+1, NT); VE = zeros(Nc+1, NT); US = zeros(Nc+1, NT); UN = zeros(Nc+1, NT) # tangential
for n in 1:NT
    u = U(n); v = V(n)
    for j in 1:Nc
        UW[j, n] = u[OFF+1,    OFF+j]          # child west face  -> truth u face 65
        UE[j, n] = u[OFF+Nc+1, OFF+j]          # child east face  -> truth u face 193
        VS[j, n] = v[OFF+j,    OFF+1]
        VN[j, n] = v[OFF+j,    OFF+Nc+1]
    end
    for j in 1:Nc+1
        VW[j, n] = 0.5 * (v[OFF,      OFF+j] + v[OFF+1,    OFF+j])   # v at x = child west face
        VE[j, n] = 0.5 * (v[OFF+Nc,   OFF+j] + v[OFF+Nc+1, OFF+j])
        US[j, n] = 0.5 * (u[OFF+j,    OFF]   + u[OFF+j,    OFF+1])   # u at y = child south face
        UN[j, n] = 0.5 * (u[OFF+j,    OFF+Nc]+ u[OFF+j,    OFF+Nc+1])
    end
end

# Block-average along the boundary then interpolate back: a coarse parent seen on the child's grid.
function coarsen_space(A, R)
    R <= 1 && return copy(A)
    na, nt = size(A); nb = cld(na, R)
    B = [mean(@view A[((b-1)*R+1):min(b*R, na), n]) for b in 1:nb, n in 1:nt]
    xc = [(b - 0.5) * R + 0.5 for b in 1:nb]
    out = similar(A)
    for n in 1:nt, a in 1:na
        if a <= xc[1];        out[a, n] = B[1, n]
        elseif a >= xc[end];  out[a, n] = B[end, n]
        else
            b = searchsortedlast(xc, a)
            w = (a - xc[b]) / (xc[b+1] - xc[b])
            out[a, n] = (1 - w) * B[b, n] + w * B[b+1, n]
        end
    end
    return out
end

# Degraded parent: coarse in space, sparse in time. Carries its own time metadata so the child still
# interpolates linearly between the samples it is actually given.
struct Parent{A}; data::A; t0::Float64; dt::Float64; nt::Int; end
function Parent(A, R, K)
    C = coarsen_space(A, R)
    idx = 1:K:size(C, 2)
    return Parent(C[:, idx], T0, K * DT, length(idx))
end
@inline function samp(p::Parent, a, t)
    τ = (t - p.t0) / p.dt
    n = clamp(unsafe_trunc(Int, τ) + 1, 1, p.nt - 1)
    w = clamp(τ - (n - 1), 0.0, 1.0)
    @inbounds return (1 - w) * p.data[a, n] + w * p.data[a, n+1]
end

# discrete_form closures: (i or j, k, grid, clock, fields) -> exterior value
const pUW = Parent(UW, PARENT_R, PARENT_K); const pUE = Parent(UE, PARENT_R, PARENT_K)
const pVS = Parent(VS, PARENT_R, PARENT_K); const pVN = Parent(VN, PARENT_R, PARENT_K)
const pVW = Parent(VW, PARENT_R, PARENT_K); const pVE = Parent(VE, PARENT_R, PARENT_K)
const pUS = Parent(US, PARENT_R, PARENT_K); const pUN = Parent(UN, PARENT_R, PARENT_K)

@inline ew(j, k, grid, clock, f) = samp(pUW, j, clock.time)
@inline ee(j, k, grid, clock, f) = samp(pUE, j, clock.time)
@inline es(i, k, grid, clock, f) = samp(pVS, i, clock.time)
@inline en(i, k, grid, clock, f) = samp(pVN, i, clock.time)
@inline tvw(j, k, grid, clock, f) = samp(pVW, j, clock.time)
@inline tve(j, k, grid, clock, f) = samp(pVE, j, clock.time)
@inline tus(i, k, grid, clock, f) = samp(pUS, i, clock.time)
@inline tun(i, k, grid, clock, f) = samp(pUN, i, clock.time)

child = RectilinearGrid(CPU(); size = (Nc, Nc),
                        x = (OFF*Δ, (OFF+Nc)*Δ), y = (OFF*Δ, (OFF+Nc)*Δ),
                        topology = (Bounded, Bounded, Flat))

function build_bcs(scheme, tangential::Bool)
    nf(f) = isnothing(scheme) ? NormalFlowBoundaryCondition(f; discrete_form = true) :
                                NormalFlowBoundaryCondition(f; scheme, discrete_form = true)
    vb(f) = isnothing(scheme) ? ValueBoundaryCondition(f; discrete_form = true) :
                                ValueBoundaryCondition(f; scheme, discrete_form = true)
    u_kw = Dict{Symbol,Any}(:west => nf(ew), :east => nf(ee))
    v_kw = Dict{Symbol,Any}(:south => nf(es), :north => nf(en))
    if tangential
        u_kw[:south] = vb(tus); u_kw[:north] = vb(tun)
        v_kw[:west]  = vb(tvw); v_kw[:east]  = vb(tve)
    end
    return (; u = FieldBoundaryConditions(; u_kw...), v = FieldBoundaryConditions(; v_kw...))
end

# truth restricted to the child box, for the initial condition and the error metric
truth_u(n) = U(n)[OFF+1:OFF+Nc+1, OFF+1:OFF+Nc]
truth_v(n) = V(n)[OFF+1:OFF+Nc,   OFF+1:OFF+Nc+1]

# Enstrophy <ζ²> over the child interior. Pointwise RMS is useless once the wake goes chaotic -- the
# child decorrelates from truth whatever the boundary does, and that swamps the boundary signal
# (measured: every arm ~0.25, i.e. a quarter of the mean flow). The enstrophy RATIO against truth is
# a STATISTICAL measure and survives decorrelation: >1 means the boundary is manufacturing vorticity,
# <1 means it is damping it. This is the diagnostic Oceananigans PR #5962 reports.
function enstrophy(u, v, m)
    s = 0.0; n = 0
    @inbounds for j in m, i in m
        ζ = (v[i, j] - v[i-1, j]) / Δ - (u[i, j] - u[i, j-1]) / Δ
        s += ζ^2; n += 1
    end
    return s / n
end

function run_arm(name, scheme, tangential)
    model = NonhydrostaticModel(child; advection = WENO(order = 5),
                                closure = ScalarDiffusivity(ν = 1e-5),
                                boundary_conditions = build_bcs(scheme, tangential))
    u0, v0 = truth_u(1), truth_v(1)
    set!(model, u = (x,y) -> 0.0, v = (x,y) -> 0.0)
    interior(model.velocities.u) .= reshape(u0, size(u0)..., 1)
    interior(model.velocities.v) .= reshape(v0, size(v0)..., 1)
    model.clock.time = T0

    Δt = 0.2 * Δ / 1.0
    nsteps = round(Int, (TT[end] - T0) / Δt)
    errs = Float64[]; ratios = Float64[]
    # interior only: drop 8 cells next to each boundary, so we measure contamination that has
    # propagated INWARD rather than the boundary cells themselves
    m = 9:Nc-8
    for s in 1:nsteps
        time_step!(model, Δt)
        if s % 50 == 0
            t = model.clock.time
            n = clamp(round(Int, (t - T0)/DT) + 1, 1, NT)
            uc = Array(interior(model.velocities.u))[:, :, 1]
            vc = Array(interior(model.velocities.v))[:, :, 1]
            push!(errs, sqrt(mean((uc[m, m] .- truth_u(n)[m, m]).^2)))
            push!(ratios, enstrophy(uc, vc, m) / enstrophy(truth_u(n), truth_v(n), m))
        end
    end
    @printf("%-12s  enstrophy ratio: final %.4f  mean %.4f  |dev| %.4f   (RMS(u) %.4e)\n",
            name, ratios[end], mean(ratios), abs(mean(ratios) - 1), mean(errs))
    return mean(ratios)
end

@printf("truth frames %d, t=%.2f..%.2f, Δt_out=%.3f\n", NT, T0, TT[end], DT)
@printf("PARENT degraded: %dx coarser in space (%.4f -> %.4f), %dx sparser in time (%.3f -> %.3f)\n\n",
        PARENT_R, Δ, PARENT_R*Δ, PARENT_K, DT, PARENT_K*DT)
run_arm("dirichlet",  nothing, false)
run_arm("normal",     NormalRadiation(inflow_timescale=0, outflow_timescale=Inf), false)
run_arm("normal+tan", NormalRadiation(inflow_timescale=0, outflow_timescale=Inf), true)
run_arm("oblique",    ObliqueRadiation(inflow_timescale=0, outflow_timescale=Inf), false)
run_arm("oblique+tan",ObliqueRadiation(inflow_timescale=0, outflow_timescale=Inf), true)
