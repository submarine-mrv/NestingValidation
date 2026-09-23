# animate_tracer.jl -- the passive-tracer plume from the Gran Canaria point release, 90 days.
#
# This is the run's reason for existing: where water from a surface point source on Gran Canaria
# goes, and how long it stays. Four views of the same plume:
#   log10(c)  -- pathway over four decades, because the plume spans that much between core and edge
#   mean age  -- ac/c, the Deleersnijder age: how old the water at each point is
#   c on T    -- plume outline over surface temperature, to read it against the upwelling front
#   surface speed with the plume contour, to see which current is carrying it
#
# Usage: julia --project=<repo> canary/animate_tracer.jl [tag] [dir]

using Oceananigans, CairoMakie, Printf, Statistics
using Oceananigans.ImmersedBoundaries: immersed_inactive_node

tag = length(ARGS) > 0 ? ARGS[1] : "_zstar_jraryf_tracer3mo"
dir = length(ARGS) > 1 ? ARGS[2] : joinpath(@__DIR__, "out3mo")
P   = joinpath(dir, "canary_glorys_")

S = P*"surface"*tag*".jld2"
cs = FieldTimeSeries(S, "c"); acs = FieldTimeSeries(S, "ac")
Ts = FieldTimeSeries(S, "T"); us = FieldTimeSeries(S, "u"); vs = FieldTimeSeries(S, "v")
grid = Ts.grid; Nx, Ny, Nz = size(grid); times = Ts.times; nt = length(times)
λ = λnodes(grid, Center(), Center(), Center()); φ = φnodes(grid, Center(), Center(), Center())
wet = [immersed_inactive_node(i, j, Nz, grid, Center(), Center(), Center()) ? NaN : 1.0 for i in 1:Nx, j in 1:Ny]
sfc(fts, n) = Array(interior(fts[n]))[1:Nx, 1:Ny, 1] .* wet
spd(n) = sqrt.(sfc(us, n).^2 .+ sfc(vs, n).^2)

# log10(c) with a floor: c spans ~4 decades and a linear scale shows only the release cell.
const CFLOOR = 1e-10
logc(n) = (a = sfc(cs, n); @. log10(max(a, CFLOOR)))
# Age is only meaningful in the plume CORE. Unbounded WENO undershoots to c < 0 in the dilute edge
# (6.8% of surface cells at day 90), and ac/c there is numerical noise: at c > 1e-6*cmax the age
# reaches 29000 days in a 90-day run. Masking at 1% of the instantaneous max keeps the map physical
# (max 117 d at day 90, 0.07% of cells above elapsed time). The mask stays: bound-preserving WENO
# would be the principled fix but is non-conservative under z* -- see docs/canary-domain-sizing.md.
const AGE_MASK_FRAC = 0.01
function age_days(n)
    c = sfc(cs, n); ac = sfc(acs, n)
    cmax = maximum(filter(isfinite, c))
    thr = AGE_MASK_FRAC * cmax
    @. ifelse(c > thr, ac / c / 86400, NaN)
end

allv(f, sub) = filter(isfinite, vcat((vec(f(n)) for n in sub)...))
sub = unique(round.(Int, range(max(1, nt÷10), nt, length = min(nt, 12))))
cq = quantile(allv(logc, sub), (0.60, 0.9999))
aq = (0.0, quantile(allv(age_days, sub), 0.995))
vq = quantile(allv(spd, sub), 0.99)
Tq = quantile(allv(n -> sfc(Ts, n), sub), (0.005, 0.995))
@info @sprintf("log10 c %.2f..%.2f | age 0..%.1f d | speed <= %.3f | %d frames", cq..., aq[2], vq, nt)

n = Observable(1)
title = @lift(@sprintf("Gran Canaria tracer release, 0.02° (2 km), JRA55-forced  —  day %.1f", times[$n]/86400))
fig = Figure(size = (1700, 1020))
Label(fig[0, 1:4], title; fontsize = 21, font = :bold)
panels = ((@lift(logc($n)),     "log₁₀ tracer concentration — the pathway",       :turbo,   cq,      λ, φ),
          (@lift(age_days($n)), "mean age (days) — plume core only (c > 1% of max)", :viridis, aq,      λ, φ),
          (@lift(sfc(Ts, $n)),  "surface T (°C) with plume outline",              :thermal, Tq,      λ, φ),
          (@lift(spd($n)),      "surface speed (m/s) with plume outline",         :speed,   (0, vq), λ, φ))
for (m, (data, lab, cmap, rng, xs, ys)) in enumerate(panels)
    r, c = fld1(m, 2), mod1(m, 2)
    ax = Axis(fig[r, 2c-1]; title = lab, titlesize = 13, xlabel = "λ", ylabel = "φ", aspect = DataAspect())
    # highclip/lowclip: values off the ends keep the end colour rather than going white.
    hm = heatmap!(ax, xs, ys, data; colormap = cmap, colorrange = rng, nan_color = :gray70,
                  highclip = Makie.to_colormap(cmap)[end], lowclip = Makie.to_colormap(cmap)[1])
    m >= 3 && contour!(ax, λ, φ, @lift(logc($n)); levels = [cq[1] + 0.5(cq[2]-cq[1])],
                       color = :white, linewidth = 1.2)
    Colorbar(fig[r, 2c], hm; width = 12)
end
colgap!(fig.layout, 6)
out = joinpath(@__DIR__, "canary_tracer_90day.mp4")
record(fig, out, 1:nt; framerate = 6) do i; n[] = i; end
n[] = nt; save(joinpath(@__DIR__, "canary_tracer_day90.png"), fig)
@info "Saved $out"
