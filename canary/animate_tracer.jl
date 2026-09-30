# animate_tracer.jl -- surface view of the Gran Canaria release and the flow carrying it.
#
# Four panels, all on a linear scale:
#   sea surface height  -- diverging map centred on zero, not on the data
#   surface temperature -- the upwelling front
#   surface speed       -- which current is doing the transport
#   tracer c (absolute) -- where the MASS is, rather than where the plume can be detected
#
# The tracer panel is deliberately linear. log10(c) renders a 1e-9 filament and a 1e-5 core at
# comparable brightness, which makes the plume look as though it fills the domain when nearly all of
# its mass sits near the source. `ac` is still written by the run, so the Deleersnijder age can be
# recovered; it is just not plotted here.
#
# Usage: julia --project=<repo> canary/animate_tracer.jl [tag] [dir] [max_days]
#
# max_days truncates the render at a given simulated day, so a run still in flight can be animated
# to a round number rather than to wherever it happens to have reached.

using Oceananigans, CairoMakie, Printf, Statistics
using Oceananigans.ImmersedBoundaries: immersed_inactive_node

tag = length(ARGS) > 0 ? ARGS[1] : "_zstar_jraryf_tracer3mo"
dir = length(ARGS) > 1 ? ARGS[2] : joinpath(@__DIR__, "out3mo")
max_days = length(ARGS) > 2 ? parse(Float64, ARGS[3]) : Inf
# Fraction of the peak concentration at which the tracer colour scale saturates. The plume spans
# ~4 decades and nearly all of its mass sits in a few cells at the source, so a range set by the
# peak renders everything else as the bottom pixel. Saturating the core is the point.
c_sat = length(ARGS) > 3 ? parse(Float64, ARGS[4]) : 0.05
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
# Truncate at max_days if asked. Frames are counted, not interpolated, so the last frame is the last
# one at or before the requested day.
if isfinite(max_days)
    keep = count(t -> t / 86400 <= max_days + 1e-9, @view times[1:nt])
    keep == 0 && error("no frames at or before day $(max_days); first frame is day $(times[1]/86400)")
    nt = keep
end
@info @sprintf("rendering %d frames, day %.2f to %.2f", nt, times[1]/86400, times[nt]/86400)
λ = λnodes(grid, Center(), Center(), Center()); φ = φnodes(grid, Center(), Center(), Center())
wet = [immersed_inactive_node(i, j, Nz, grid, Center(), Center(), Center()) ? NaN : 1.0 for i in 1:Nx, j in 1:Ny]
sfc(fts, n) = Array(interior(fts[n]))[1:Nx, 1:Ny, 1] .* wet
spd(n) = sqrt.(sfc(us, n).^2 .+ sfc(vs, n).^2)
# SSH ANOMALY: the domain-mean offset is barotropic adjustment and carries no information, while
# the spatial gradients are what drive the geostrophic flow. Centring a diverging map on zero when
# the field is almost entirely negative wastes half the colour range.
function ssh(n)
    isnothing(ηs) && return fill(NaN, Nx, Ny)
    a = Array(interior(ηs[n]))[1:Nx, 1:Ny, 1] .* wet
    a .- mean(filter(isfinite, a))
end

# log10(c) with a floor: c spans ~4 decades, so a linear scale cannot show the pathway.
# Absolute concentration alongside it, because the log panel deliberately hides magnitude: it makes a
# 1e-9 filament and a 1e-5 core look comparably bright. The colour range is capped at a high quantile
# rather than the max, which is the release cell itself and would flatten everything else to zero.
linc(n) = sfc(cs, n)

allv(f, sub) = filter(isfinite, vcat((vec(f(n)) for n in sub)...))
sub = unique(round.(Int, range(max(1, nt÷10), nt, length = min(nt, 12))))
cpk = maximum(allv(linc, sub))
lq  = (0.0, c_sat * cpk)
vq = quantile(allv(spd, sub), 0.99)
Tq = quantile(allv(n -> sfc(Ts, n), sub), (0.005, 0.995))
ηm = isnothing(ηs) ? 1.0 : quantile(abs.(allv(ssh, sub)), 0.99); ηq = (-ηm, ηm)
@info @sprintf("c peak %.3e, saturating at %.0f%% = %.3e | speed <= %.3f | |SSH'| <= %.3f | %d frames",
               cpk, 100c_sat, lq[2], vq, ηm, nt)

n = Observable(1)
title = @lift(@sprintf("Gran Canaria tracer release, 0.02° (2 km), JRA55-forced  —  day %.1f", times[$n]/86400))
fig = Figure(size = (1700, 1020))
Label(fig[0, 1:4], title; fontsize = 21, font = :bold)
panels = ((@lift(ssh($n)),     "sea surface height anomaly (m)",  :balance, ηq),
          (@lift(sfc(Ts, $n)), "surface T (°C)",                  :thermal, Tq),
          (@lift(spd($n)),     "surface speed (m/s)",             :speed,   (0, vq)),
          (@lift(linc($n)),    "tracer concentration (absolute)", :turbo,   lq))
for (m, (data, lab, cmap, rng)) in enumerate(panels)
    r, c = fld1(m, 2), mod1(m, 2)
    ax = Axis(fig[r, 2c-1]; title = lab, titlesize = 13, xlabel = "λ", ylabel = "φ", aspect = DataAspect())
    # highclip/lowclip: values off the ends keep the end colour rather than going white.
    hm = heatmap!(ax, λ, φ, data; colormap = cmap, colorrange = rng, nan_color = :gray70,
                  highclip = Makie.to_colormap(cmap)[end], lowclip = Makie.to_colormap(cmap)[1])
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
