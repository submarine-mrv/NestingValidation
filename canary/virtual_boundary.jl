# virtual_boundary.jl -- how far south must the domain reach to capture a majority of the CDR?
#
# The run only holds data inside its own boundaries, so a larger domain cannot be measured directly.
# Instead measure a family of NESTED sub-domains: move a virtual southern boundary north, and for each
# one record the tracer inventory that remains north of it. That gives tau(L) over a measured range of
# L, which can be fitted and extrapolated south.
#
# Bias to keep in mind: tracer that crosses a virtual boundary southward can return, and here it is
# counted as retained again. A real boundary there would have flushed it (the exterior carries no
# tracer), so these inventories are an UPPER bound on retention, hence tau(L) is an upper bound too.
#
# Writes one CSV row per 3-D snapshot frame. Usage:
#   julia --project=<repo> canary/virtual_boundary.jl <snapshot3d.jld2> <out.csv>

using Oceananigans, Printf
using Oceananigans.Operators: Vᶜᶜᶜ
using Oceananigans.ImmersedBoundaries: immersed_inactive_node

file = length(ARGS) > 0 ? ARGS[1] : error("need snapshot3d path")
out  = length(ARGS) > 1 ? ARGS[2] : "virtual_boundary.csv"

cs = FieldTimeSeries(file, "c")
grid = cs.grid
Nx, Ny, Nz = size(grid)
times = cs.times
φ = φnodes(grid, Center(), Center(), Center())

# Candidate virtual southern boundaries, kept at least 0.75 deg north of the domain's OWN southern
# edge so none sits inside the real sponge (0.25 deg Gaussian, felt to ~0.75 deg). The six latitudes
# the phi_S = 25.0 run measured are always offered, so the two domains can be compared at identical
# boundaries; a domain that reaches further south picks up the extras below as well.
φ_edge   = φ[1] - (φ[2] - φ[1]) / 2
φ_min_ok = φ_edge + 0.75 - 1e-9
φ_bounds = sort(filter(>=(φ_min_ok), [23.0, 23.5, 24.0, 24.5, 25.0, 25.25, 25.5,
                                      25.75, 26.0, 26.5, 27.0, 27.25, 27.5]))
@info "domain southern edge $(round(φ_edge, digits=3)) -> virtual boundaries $(φ_bounds)"
j_bounds = [findfirst(>=(φb), φ) for φb in φ_bounds]

# Static cell volumes. z* makes these time-dependent at the ~0.1% level (free-surface displacement
# over column depth), which is far below the signal here.
vol = zeros(Nx, Ny, Nz)
for k in 1:Nz, j in 1:Ny, i in 1:Nx
    vol[i, j, k] = immersed_inactive_node(i, j, k, grid, Center(), Center(), Center()) ?
                   0.0 : Vᶜᶜᶜ(i, j, k, grid)
end

open(out, "w") do io
    println(io, "t_days,total," * join(["north_$(φb)" for φb in φ_bounds], ","))
    for n in 1:length(times)
        c = Array(interior(cs[n]))
        cv = c .* vol
        total = sum(cv)
        norths = [sum(@view cv[:, j:end, :]) for j in j_bounds]
        @printf(io, "%.4f,%.8e,%s\n", times[n]/86400, total,
                join((@sprintf("%.8e", x) for x in norths), ","))
        # Label the probe from φ_bounds, not a literal: the list is now domain-dependent, so a
        # hard-coded "27.0" silently mislabelled norths[4] as 24.5 on the phi_S = 22 run.
        @printf("frame %2d  t=%6.2f d  total=%.4e  north(%.2f)=%.4e\n",
                n, times[n]/86400, total, φ_bounds[end], norths[end])
        flush(stdout)
    end
end
@info "wrote $out"
