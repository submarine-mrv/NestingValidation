# animate_canary.jl -- surface fields and the cross-shore upwelling section over the 7-day run.
#
# Surface T is the signal this domain exists to produce: cold water drawn up against the African
# coast, and the filaments and island wakes that a 2 km child resolves and a 9 km parent cannot.
# Speed shows the wakes directly.
#
# 6-hourly frames, t = 0.25 .. 7.00 d (TimeInterval writes at the first interval, so there is no t = 0
# frame). Usage: julia --project=<repo> canary/animate_canary.jl [tag]

using Oceananigans, CairoMakie, Printf, Statistics
using Oceananigans.ImmersedBoundaries: immersed_inactive_node

tag = length(ARGS) > 0 ? ARGS[1] : "_zstar_jraryf"
sfx = tag == "_zstar_jraryf" ? "" : tag
P = joinpath(@__DIR__, "canary_glorys_")

Ts = FieldTimeSeries(P*"surface"*tag*".jld2", "T")
Ss = FieldTimeSeries(P*"surface"*tag*".jld2", "S")
us = FieldTimeSeries(P*"surface"*tag*".jld2", "u")
vs = FieldTimeSeries(P*"surface"*tag*".jld2", "v")
grid = Ts.grid; Nx, Ny, Nz = size(grid); times = Ts.times; nt = length(times)
λ = λnodes(grid, Center(), Center(), Center()); φ = φnodes(grid, Center(), Center(), Center())
wet = [immersed_inactive_node(i, j, Nz, grid, Center(), Center(), Center()) ? NaN : 1.0 for i in 1:Nx, j in 1:Ny]
sfc(fts, n) = Array(interior(fts[n]))[1:Nx, 1:Ny, 1] .* wet
spd(n) = sqrt.(sfc(us, n).^2 .+ sfc(vs, n).^2)

# Fixed ranges across the animation, or apparent change is just the scale moving. Percentiles for T
# so the coastal upwelling core is not clipped; speed capped at the 99th so eddies stay legible.
allv(f) = filter(isfinite, vcat((vec(f(n)) for n in 1:nt)...))
Tq = quantile(allv(n -> sfc(Ts, n)), (0.005, 0.995))
Sq = quantile(allv(n -> sfc(Ss, n)), (0.01, 0.99))
vq = quantile(allv(spd), 0.99)
@info @sprintf("T %.2f..%.2f | S %.3f..%.3f | speed ≤ %.3f m/s | %d frames", Tq..., Sq..., vq, nt)

n = Observable(1)
title = @lift(@sprintf("Canary upwelling, 0.02° (2 km), JRA55-forced  —  day %.2f", times[$n]/86400))
fig = Figure(size = (1500, 1000))
Label(fig[0, 1:2], title; fontsize = 21, font = :bold)
for (row, (data, lab, cmap, rng)) in enumerate((
        (@lift(sfc(Ts, $n)), "surface T (°C) — cold water against the coast is the upwelling", :thermal, Tq),
        (@lift(spd($n)),     "surface speed (m/s) — island wakes and filaments",               :speed,   (0, vq)),
        (@lift(sfc(Ss, $n)), "surface S (psu)",                                                :haline,  Sq)))
    ax = Axis(fig[row, 1]; title = lab, titlesize = 14, xlabel = "λ", ylabel = "φ", aspect = DataAspect())
    # No highclip: Makie then clamps out-of-range values to the colormap's end colours. An explicit
    # highclip (white here, before) reads as missing data rather than "off the top of the scale",
    # which matters because these ranges are percentile-capped and ~1% of cells sit above them.
    hm = heatmap!(ax, λ, φ, data; colormap = cmap, colorrange = rng, nan_color = :gray70)
    Colorbar(fig[row, 2], hm; width = 12)
end
colgap!(fig.layout, 6)
out = joinpath(@__DIR__, "canary_surface_7day$(sfx).mp4")
record(fig, out, 1:nt; framerate = 4) do i; n[] = i; end
n[] = nt; save(joinpath(@__DIR__, "canary_surface_day7$(sfx).png"), fig)
@info "Saved $out"

# --- cross-shore section at 26.5N: the upwelling cell itself ------------------
Tsec = FieldTimeSeries(P*"shelf_section"*tag*".jld2", "T")
vsec = FieldTimeSeries(P*"shelf_section"*tag*".jld2", "v")
j = argmin(abs.(φ .- 26.5))
zc = znodes(grid, Center(), Center(), Center())
wetc = [immersed_inactive_node(i, j, k, grid, Center(), Center(), Center()) ? NaN : 1.0 for i in 1:Nx, k in 1:Nz]
slabT(m) = Array(interior(Tsec[m]))[1:Nx, 1, 1:Nz] .* wetc
slabv(m) = Array(interior(vsec[m]))[1:Nx, 1, 1:Nz] .* wetc
Tsq = quantile(filter(isfinite, vcat((vec(slabT(m)) for m in 1:nt)...)), (0.01, 0.999))
vsq = quantile(abs.(filter(isfinite, vcat((vec(slabv(m)) for m in 1:nt)...))), 0.995)

m = Observable(1)
t2 = @lift(@sprintf("Cross-shore section at %.2f°N  —  day %.2f   (+v = northward)", φ[j], times[$m]/86400))
fig2 = Figure(size = (1500, 780))
Label(fig2[0, 1:2], t2; fontsize = 20, font = :bold)
for (row, (data, lab, cmap, rng)) in enumerate((
        (@lift(slabT($m)), "T (°C) — isotherms tilting up toward the coast", :thermal, Tsq),
        (@lift(slabv($m)), "v (m/s) — Canary Current above, poleward undercurrent below", :balance, (-vsq, vsq))))
    ax = Axis(fig2[row, 1]; title = lab, titlesize = 14, xlabel = "λ", ylabel = "z (m)")
    hm = heatmap!(ax, λ, zc, data; colormap = cmap, colorrange = rng, nan_color = :gray70)
    ylims!(ax, -1500, 0)
    Colorbar(fig2[row, 2], hm; width = 12)
end
out2 = joinpath(@__DIR__, "canary_section_7day$(sfx).mp4")
record(fig2, out2, 1:nt; framerate = 4) do i; m[] = i; end
m[] = nt; save(joinpath(@__DIR__, "canary_section_day7$(sfx).png"), fig2)
@info "Saved $out2"
