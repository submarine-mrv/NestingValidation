# animate_w.jl -- surface view with VERTICAL velocity at depth in place of the tracer.
#
# w at the surface is ~1e-4 m/s (the kinematic condition makes it nearly zero) against ~2e-2 m/s in
# the interior, so a surface-level w panel shows nothing. This reads w from the 3-D snapshots at a
# chosen depth -- default ~100 m, which is the depth the run already treats as the upwelling source.
#
# Snapshots are written every CANARY_SNAPSHOT_HOURS (72 h here), so this animation is far coarser in
# time than animate_tracer.jl. SSH, SST and speed come from the same snapshots for consistency
# rather than from the 6-hourly surface file.
#
# Usage: julia --project=<env> canary/animate_w.jl [tag] [dir] [max_days] [depth_m]

using Oceananigans, CairoMakie, Printf, Statistics
using Oceananigans.ImmersedBoundaries: immersed_inactive_node

tag      = length(ARGS) > 0 ? ARGS[1] : "_zstar_jraryf_oblique_res5000_tsreservoir_ref90"
dir      = length(ARGS) > 1 ? ARGS[2] : joinpath(@__DIR__, "ref90")
max_days = length(ARGS) > 2 ? parse(Float64, ARGS[3]) : Inf
depth_m  = length(ARGS) > 3 ? parse(Float64, ARGS[4]) : 100.0

P  = joinpath(dir, "canary_glorys_")
S3 = P*"snapshot3d"*tag*".jld2"
FS = P*"free_surface"*tag*".jld2"

ws = FieldTimeSeries(S3, "w"); Ts = FieldTimeSeries(S3, "T")
us = FieldTimeSeries(S3, "u"); vs = FieldTimeSeries(S3, "v")
ηs = isfile(FS) ? FieldTimeSeries(FS, "η") : nothing
grid = Ts.grid; Nx, Ny, Nz = size(grid); times = Ts.times
nt = length(times)
if isfinite(max_days)
    nt = count(t -> t/86400 <= max_days + 1e-9, @view times[1:nt])
end

λ = λnodes(grid, Center(), Center(), Center()); φ = φnodes(grid, Center(), Center(), Center())
zc = znodes(grid, Center(), Center(), Center())
kw = argmin(abs.(zc .+ depth_m))                 # nearest cell centre to the requested depth
wet2 = [immersed_inactive_node(i, j, Nz, grid, Center(), Center(), Center()) ? NaN : 1.0 for i in 1:Nx, j in 1:Ny]
wetw = [immersed_inactive_node(i, j, kw, grid, Center(), Center(), Center()) ? NaN : 1.0 for i in 1:Nx, j in 1:Ny]
@info @sprintf("w at z = %.1f m (level %d of %d); %d frames, day %.2f to %.2f",
               zc[kw], kw, Nz, nt, times[1]/86400, times[nt]/86400)

sfc(fts, n) = Array(interior(fts[n]))[1:Nx, 1:Ny, Nz] .* wet2
wlev(n)     = Array(interior(ws[n]))[1:Nx, 1:Ny, kw] .* wetw
spd(n)      = sqrt.(sfc(us, n).^2 .+ sfc(vs, n).^2)
# η is on the 6-hourly schedule, so pick its frame by time rather than by index. Plot the ANOMALY:
# the domain-mean offset is barotropic adjustment and carries no information, and centring a
# diverging map on zero when the field is almost entirely negative wastes half the colour range.
function ssh(n)
    isnothing(ηs) && return fill(NaN, Nx, Ny)
    j = argmin(abs.(Float64.(collect(ηs.times)) .- times[n]))
    a = Array(interior(ηs[j]))[1:Nx, 1:Ny, 1] .* wet2
    a .- mean(filter(isfinite, a))
end

allv(f, sub) = filter(isfinite, vcat((vec(f(n)) for n in sub)...))
sub = unique(round.(Int, range(max(1, nt÷10), nt, length = min(nt, 10))))
wm = quantile(abs.(allv(wlev, sub)), 0.995); wq = (-wm, wm)
vq = quantile(allv(spd, sub), 0.99)
Tq = quantile(allv(n -> sfc(Ts, n), sub), (0.005, 0.995))
ηm = isnothing(ηs) ? 1.0 : quantile(abs.(allv(ssh, sub)), 0.99); ηq = (-ηm, ηm)
@info @sprintf("|w| <= %.2e m/s | speed <= %.3f | |SSH| <= %.3f", wm, vq, ηm)

n = Observable(1)
title = @lift(@sprintf("Canary 0.02° — day %.1f", times[$n]/86400))
fig = Figure(size = (1700, 1020))
Label(fig[0, 1:4], title; fontsize = 21, font = :bold)
panels = ((@lift(ssh($n)),     "sea surface height anomaly (m)",                  :balance, ηq),
          (@lift(sfc(Ts, $n)), "surface T (°C)",                                  :thermal, Tq),
          (@lift(spd($n)),     "surface speed (m/s)",                             :speed,   (0, vq)),
          (@lift(wlev($n)),    @sprintf("vertical velocity at %.0f m (m/s)", -zc[kw]), :balance, wq))
for (m, (data, lab, cmap, rng)) in enumerate(panels)
    r, c = fld1(m, 2), mod1(m, 2)
    ax = Axis(fig[r, 2c-1]; title = lab, titlesize = 13, xlabel = "λ", ylabel = "φ", aspect = DataAspect())
    hm = heatmap!(ax, λ, φ, data; colormap = cmap, colorrange = rng, nan_color = :gray70,
                  highclip = Makie.to_colormap(cmap)[end], lowclip = Makie.to_colormap(cmap)[1])
    Colorbar(fig[r, 2c], hm; width = 12)
end
colgap!(fig.layout, 6)
out = joinpath(@__DIR__, "canary_w_$(round(Int, depth_m))m$(replace(tag, "_zstar_jraryf" => "")).mp4")
record(fig, out, 1:nt; framerate = 3) do i; n[] = i; end
n[] = nt; save(replace(out, ".mp4" => ".png"), fig)
@info "Saved $out"
