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
# Free surface lives in its own file (same output schedule, so the frames line up). Guarded because
# a run can be interrupted between writers, leaving one file a frame or two short of the other.
FS = P*"free_surface"*tag*".jld2"
ηs = isfile(FS) ? FieldTimeSeries(FS, "η") : nothing
cs = FieldTimeSeries(S, "c"); acs = FieldTimeSeries(S, "ac")
Ts = FieldTimeSeries(S, "T"); us = FieldTimeSeries(S, "u"); vs = FieldTimeSeries(S, "v")
grid = Ts.grid; Nx, Ny, Nz = size(grid); times = Ts.times
nt = isnothing(ηs) ? length(times) : min(length(times), length(ηs.times))
λ = λnodes(grid, Center(), Center(), Center()); φ = φnodes(grid, Center(), Center(), Center())
wet = [immersed_inactive_node(i, j, Nz, grid, Center(), Center(), Center()) ? NaN : 1.0 for i in 1:Nx, j in 1:Ny]
sfc(fts, n) = Array(interior(fts[n]))[1:Nx, 1:Ny, 1] .* wet
spd(n) = sqrt.(sfc(us, n).^2 .+ sfc(vs, n).^2)
ssh(n) = isnothing(ηs) ? fill(NaN, Nx, Ny) : Array(interior(ηs[n]))[1:Nx, 1:Ny, 1] .* wet

# log10(c) with a floor: c spans ~4 decades, so a linear scale cannot show the pathway.
const CFLOOR = 1e-10
logc(n) = (a = sfc(cs, n); @. log10(max(a, CFLOOR)))
# Absolute concentration alongside it, because the log panel deliberately hides magnitude: it makes a
# 1e-9 filament and a 1e-5 core look comparably bright. The colour range is capped at a high quantile
# rather than the max, which is the release cell itself and would flatten everything else to zero.
linc(n) = sfc(cs, n)

allv(f, sub) = filter(isfinite, vcat((vec(f(n)) for n in sub)...))
sub = unique(round.(Int, range(max(1, nt÷10), nt, length = min(nt, 12))))
cq = quantile(allv(logc, sub), (0.60, 0.9999))
lq = (0.0, quantile(allv(linc, sub), 0.99995))
vq = quantile(allv(spd, sub), 0.99)
Tq = quantile(allv(n -> sfc(Ts, n), sub), (0.005, 0.995))
# SSH is signed, so centre the diverging range on zero rather than on the data.
ηm = isnothing(ηs) ? 1.0 : quantile(abs.(allv(ssh, sub)), 0.995); ηq = (-ηm, ηm)
@info @sprintf("log10 c %.2f..%.2f | linear c <= %.3e | speed <= %.3f | |SSH| <= %.3f | %d frames",
               cq..., lq[2], vq, ηm, nt)

n = Observable(1)
title = @lift(@sprintf("Gran Canaria tracer release, 0.02° (2 km), JRA55-forced  —  day %.1f", times[$n]/86400))
fig = Figure(size = (1700, 1020))
Label(fig[0, 1:4], title; fontsize = 21, font = :bold)
# last field: whether to overlay the plume outline. log10(c) is no longer a panel -- a linear scale
# shows where the MASS is, which is the honest view -- but it still defines the outline contour.
panels = ((@lift(ssh($n)),     "sea surface height (m) with plume outline",  :balance, ηq,      true),
          (@lift(sfc(Ts, $n)), "surface T (°C) with plume outline",          :thermal, Tq,      true),
          (@lift(spd($n)),     "surface speed (m/s) with plume outline",     :speed,   (0, vq), true),
          (@lift(linc($n)),    "tracer concentration (absolute)",            :turbo,   lq,      false))
for (m, (data, lab, cmap, rng, outline)) in enumerate(panels)
    r, c = fld1(m, 2), mod1(m, 2)
    ax = Axis(fig[r, 2c-1]; title = lab, titlesize = 13, xlabel = "λ", ylabel = "φ", aspect = DataAspect())
    # highclip/lowclip: values off the ends keep the end colour rather than going white.
    hm = heatmap!(ax, λ, φ, data; colormap = cmap, colorrange = rng, nan_color = :gray70,
                  highclip = Makie.to_colormap(cmap)[end], lowclip = Makie.to_colormap(cmap)[1])
    outline && contour!(ax, λ, φ, @lift(logc($n)); levels = [cq[1] + 0.5(cq[2]-cq[1])],
                        color = :white, linewidth = 1.2)
    Colorbar(fig[r, 2c], hm; width = 12)
end
colgap!(fig.layout, 6)
# Name the outputs from the non-default part of the tag, so a domain-sizing run does not overwrite
# the phi_S = 25 animation. The 25N tag reduces to "", keeping its original filenames.
suffix = replace(tag, "_zstar_jraryf" => "", "_tracer3mo" => "")
out = joinpath(@__DIR__, "canary_tracer_90day$(suffix).mp4")
record(fig, out, 1:nt; framerate = 6) do i; n[] = i; end
n[] = nt; save(joinpath(@__DIR__, "canary_tracer_day90$(suffix).png"), fig)
@info "Saved $out"
