# boundary_excursion_west.jl -- Salish west (Juan de Fuca) boundary. Same method as
# red_sea/boundary_excursion.jl; the section is the WEST face, so the normal component is u at i=1.
#
# Sets the TracerReservoir length scale from this domain's own data rather than from canary/'s.
# A memoryless Value boundary destroys tracer that leaves and returns; the reservoir's length scale
# is how much of that return it remembers, so the scale that matters is the distance a parcel
# actually excursions across the boundary.
#
# Per wet face cell, integrate the boundary-normal velocity in time to get displacement, and take
# its peak-to-peak range. Report the VOLUME-WEIGHTED median: a plain median is dominated by the many
# near-stagnant deep cells and by noise sign-flips, which is the error that had to be corrected when
# this was measured for canary/.
#
#   julia --project=v113_staging salish_sea/boundary_excursion_west.jl <west_section.jld2>

using Oceananigans, Printf, Statistics
using Oceananigans.Operators: Δzᶠᶜᶜ, Δyᶠᶜᶜ
using Oceananigans.ImmersedBoundaries: immersed_inactive_node

file = length(ARGS) > 0 ? ARGS[1] : error("need an east_section .jld2")
us = FieldTimeSeries(file, "u")
grid = us.grid
Nx, Ny, Nz = size(grid)
t = us.times
@printf("%d frames, %.3f -> %.3f days (Δ %.2f h)\n", length(t), t[1]/86400, t[end]/86400,
        (t[2]-t[1])/3600)

# The section is the west face, so the stored i index is the only one.
u = [Array(interior(us[n]))[1, :, :] for n in 1:length(t)]

exc = Float64[]; wgt = Float64[]
for k in 1:Nz, j in 1:Ny
    immersed_inactive_node(1, j, k, grid, Center(), Center(), Center()) && continue
    # Cumulative displacement, trapezoidal in time.
    x = 0.0; xmin = 0.0; xmax = 0.0
    for n in 2:length(t)
        x += 0.5 * (u[n][j, k] + u[n-1][j, k]) * (t[n] - t[n-1])
        xmin = min(xmin, x); xmax = max(xmax, x)
    end
    push!(exc, xmax - xmin)
    push!(wgt, Δyᶠᶜᶜ(1, j, k, grid) * Δzᶠᶜᶜ(1, j, k, grid))
end

function wmedian(v, w)
    p = sortperm(v); v, w = v[p], w[p]
    c = cumsum(w) ./ sum(w)
    return v[findfirst(>=(0.5), c)]
end

@printf("wet face cells: %d\n", length(exc))
@printf("excursion  plain median %.0f m | VOLUME-WEIGHTED median %.0f m | mean %.0f m | p90 %.0f m\n",
        median(exc), wmedian(exc, wgt), mean(exc), quantile(exc, 0.9))
# TracerReservoir relaxes as exp(-E/L) over an excursion E, so L sets how much inflow memory
# survives one excursion. No factor relates E to L -- report the trade-off and let the
# configuration choose. The "0.3x" once printed here was back-derived to describe Canary's
# 5000 m, not a result. Cross-domain spread: ../parent_excursion_survey.jl.
E = wmedian(exc, wgt)
for L in (5e3, 10e3, 20e3, 40e3, 60e3)
    @printf("   L = %5.1f km -> retains %4.1f%% of inflow memory over E = %.1f km\n",
            L/1e3, 100exp(-E/L), E/1e3)
end
