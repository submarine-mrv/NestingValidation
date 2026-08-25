# animate_surface.jl -- surface T and S over the 7-day validation run.
# Frames come from the `surface` JLD2 writer: 6-hourly, t = 0.25 .. 7.00 d (TimeInterval writes at the
# first interval, so there is no t = 0 frame and the sub-6 h adjustment transient is not resolved).
# Usage: julia --project=<repo> red_sea/animate_surface.jl

using Oceananigans, CairoMakie, Printf, Statistics
using Oceananigans.ImmersedBoundaries: immersed_inactive_node

infile  = joinpath(@__DIR__, "redsea_glorys_surface_zstar.jld2")
outfile = joinpath(@__DIR__, "redsea_surface_7day.mp4")

Ts = FieldTimeSeries(infile, "T")
Ss = FieldTimeSeries(infile, "S")
grid = Ts.grid
Nx, Ny, Nz = size(grid)
times = Ts.times
nt = length(times)

λ = λnodes(grid, Center(), Center(), Center())
φ = φnodes(grid, Center(), Center(), Center())
wet = [immersed_inactive_node(i, j, Nz, grid, Center(), Center(), Center()) ? NaN : 1.0 for i in 1:Nx, j in 1:Ny]
frame(fts, n) = Array(interior(fts[n]))[1:Nx, 1:Ny, 1] .* wet

# Fixed colour ranges across the animation, or apparent change is just the scale moving. Percentiles
# rather than extrema so one coastal cell cannot wash out the basin.
allvals(fts) = filter(isfinite, vcat((vec(frame(fts, n)) for n in 1:nt)...))
Tv, Sv = allvals(Ts), allvals(Ss)
Trange = (quantile(Tv, 0.01), quantile(Tv, 0.99))
Srange = (quantile(Sv, 0.01), quantile(Sv, 0.99))
@info @sprintf("T range %.2f..%.2f | S range %.2f..%.2f | %d frames", Trange..., Srange..., nt)

n = Observable(1)
Tdata = @lift(frame(Ts, $n))
Sdata = @lift(frame(Ss, $n))
title = @lift(@sprintf("Red Sea + Gulf of Aden, 0.05° z*  —  day %.2f", times[$n] / 86400))

fig = Figure(size = (1500, 760))
Label(fig[0, 1:4], title; fontsize = 22, font = :bold)
axT = Axis(fig[1, 1]; title = "surface T (°C)", xlabel = "λ", ylabel = "φ", aspect = DataAspect())
hmT = heatmap!(axT, λ, φ, Tdata; colormap = :thermal, colorrange = Trange, nan_color = :gray70)
Colorbar(fig[1, 2], hmT; width = 14)
axS = Axis(fig[1, 3]; title = "surface S (psu)", xlabel = "λ", ylabel = "φ", aspect = DataAspect())
hmS = heatmap!(axS, λ, φ, Sdata; colormap = :haline, colorrange = Srange, nan_color = :gray70)
Colorbar(fig[1, 4], hmS; width = 14)

record(fig, outfile, 1:nt; framerate = 4) do i
    n[] = i
end
@info "Saved $outfile ($(nt) frames, $(round(nt/4, digits=1)) s at 4 fps)"
