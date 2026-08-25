# animate_strait_section.jl -- Bab el-Mandeb exchange over the 7-day validation run.
#
# Zonal-vertical section at 12.69°N (the narrowest main channel). v is the component NORMAL to this
# section, so +v is northward, into the Red Sea. Two framings, because the section carries two
# different stories: the strait exchange itself (left, top 140 m, where the sill at 136 m confines it)
# and the Red Sea Outflow Water descending into the Gulf of Aden (right, full depth).
#
# 6-hourly frames, t = 0.25 .. 7.00 d -- no t = 0 frame, and the sub-6 h adjustment transient is not
# resolved. Usage: julia --project=<repo> red_sea/animate_strait_section.jl

using Oceananigans, CairoMakie, Printf, Statistics
using Oceananigans.ImmersedBoundaries: immersed_inactive_node

infile  = joinpath(@__DIR__, "redsea_glorys_strait_section_zstar.jld2")
outfile = joinpath(@__DIR__, "redsea_strait_section_7day.mp4")

Ss = FieldTimeSeries(infile, "S")
vs = FieldTimeSeries(infile, "v")
grid = Ss.grid
Nx, Ny, Nz = size(grid)
times = Ss.times
nt = length(times)

λ = λnodes(grid, Center(), Center(), Center())
z = znodes(grid, Center(), Center(), Center())
φ = φnodes(grid, Center(), Center(), Center())
j = argmin(abs.(φ .- 12.69))
@info @sprintf("section at %.3f°N (j = %d), %d frames", φ[j], j, nt)

# Mask from the grid, not from zero values: v is legitimately zero in water.
wetC = [immersed_inactive_node(i, j, k, grid, Center(), Center(), Center()) ? NaN : 1.0 for i in 1:Nx, k in 1:Nz]
wetF = [immersed_inactive_node(i, j, k, grid, Center(), Face(),   Center()) ? NaN : 1.0 for i in 1:Nx, k in 1:Nz]
slabS(n) = Array(interior(Ss[n]))[1:Nx, 1, 1:Nz] .* wetC
slabv(n) = Array(interior(vs[n]))[1:Nx, 1, 1:Nz] .* wetF

# S uses the FULL range, not percentiles: the dense outflow on the sill (S ≈ 37.8) is a small
# fraction of the section's cells, so a 99th-percentile cap (36.26) saturates precisely the signal
# this animation exists to show. v keeps a percentile cap -- velocity has a long tail and clipping a
# few extreme cells preserves the two-layer structure.
allv(f) = filter(isfinite, vcat((vec(f(n)) for n in 1:nt)...))
Srange = extrema(allv(slabS))
vq = quantile(abs.(allv(slabv)), 0.995)
vrange = (-vq, vq)
@info @sprintf("S %.2f..%.2f | |v| ≤ %.3f m/s", Srange..., vq)

n = Observable(1)
dS = @lift(slabS($n)); dv = @lift(slabv($n))
title = @lift(@sprintf("Bab el-Mandeb section, 12.69°N  —  day %.2f   (+v = northward, into the Red Sea)",
                       times[$n] / 86400))

fig = Figure(size = (1600, 860))
Label(fig[0, 1:4], title; fontsize = 21, font = :bold)

# left column: the strait and its sill; right column: the outflow plume at full depth
for (row, (data, lab, cmap, rng)) in enumerate(((dS, "S (psu)", :haline, Srange),
                                                (dv, "v (m/s)", :balance, vrange)))
    axz = Axis(fig[row, 1]; title = "strait — $lab", xlabel = "λ", ylabel = "z (m)")
    hm = heatmap!(axz, λ, z, data; colormap = cmap, colorrange = rng, nan_color = :gray70)
    xlims!(axz, 43.12, 43.85); ylims!(axz, -140, 0)   # main channel only; 44.5E is a shallow side pocket
    hlines!(axz, [-136]; color = :black, linestyle = :dash, linewidth = 1)   # Hanish sill
    Colorbar(fig[row, 2], hm; width = 13)

    axf = Axis(fig[row, 3]; title = "full section — $lab", xlabel = "λ", ylabel = "z (m)")
    hm2 = heatmap!(axf, λ, z, data; colormap = cmap, colorrange = rng, nan_color = :gray70)
    xlims!(axf, 43.1, 50.0); ylims!(axf, -1600, 0)
    Colorbar(fig[row, 4], hm2; width = 13)
end

n[] = nt
save(joinpath(@__DIR__, "redsea_strait_section_day7.png"), fig)

record(fig, outfile, 1:nt; framerate = 4) do i
    n[] = i
end
@info "Saved $outfile ($nt frames, $(round(nt/4, digits=1)) s at 4 fps)"
