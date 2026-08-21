# inspect_glorys_run.jl
#
# Diagnose the instability onset seen in the bounded 1-hour GLORYS run: read the JLD2
# surface output and (a) print per-frame extrema over time, (b) map the final frame's
# surface S, w, and speed with land masked -- to see WHERE the spurious S overshoot and
# the max|w| jump live. If they sit on the west/north OPEN boundaries, it's the Step-A
# constant-vs-GLORYS tracer-BC conflict (Step B fixes it). If interior/general, it's the
# wizard/Δt (tighten max_Δt, max_change, cfl).

using Oceananigans
using CairoMakie
using Printf
using Oceananigans.Fields: location
using Oceananigans.ImmersedBoundaries: immersed_inactive_node

file = "salishsea_glorys_openboundaries_test1.jld2"

Ts = FieldTimeSeries(file, "T")
Ss = FieldTimeSeries(file, "S")
us = FieldTimeSeries(file, "u")
vs = FieldTimeSeries(file, "v")
ws = FieldTimeSeries(file, "w")

times = Ts.times
grid  = Ts.grid
Nx, Ny, Nz = size(grid)
@info "Loaded $(length(times)) frames from $file" t_min=first(times)/60 t_max=last(times)/60

# 1-D lon/lat in degrees (λnodes/φnodes, not xnodes/ynodes which give metres)
snap_λ(f) = λnodes(f.grid, map(L -> L(), location(f))...)
snap_φ(f) = φnodes(f.grid, map(L -> L(), location(f))...)

# surface wet mask (NaN on land)
wet2d = [immersed_inactive_node(i, j, Nz, grid, Center(), Center(), Center()) ? NaN : 1.0
         for i in 1:Nx, j in 1:Ny]

# saved fields are single-surface-level -> interior k index is 1; crop Face fields to Nx×Ny
sfc(fts, n) = Array(interior(fts[n]))[1:Nx, 1:Ny, 1]

# --- (a) extrema over time ---------------------------------------------------
println("\n  t(min)    Tmin    Tmax    Smin    Smax     max|w|    max|u|    max|v|")
for n in 1:length(times)
    T, S = interior(Ts[n]), interior(Ss[n])
    u, v, w = interior(us[n]), interior(vs[n]), interior(ws[n])
    @printf("  %6.1f  %6.3f  %6.3f  %6.3f  %6.3f   %.2e  %.2e  %.2e\n",
            times[n]/60, minimum(T), maximum(T), minimum(S), maximum(S),
            maximum(abs, w), maximum(abs, u), maximum(abs, v))
end

# --- (b) final-frame maps ----------------------------------------------------
n   = length(times)
Sm  = sfc(Ss, n) .* wet2d
wm  = sfc(ws, n) .* wet2d
spd = sqrt.(sfc(us, n).^2 .+ sfc(vs, n).^2) .* wet2d
λ, φ = snap_λ(Ss), snap_φ(Ss)

fig = Figure(size = (1500, 520))
for (i, (d, lab, cmap)) in enumerate(((Sm,  "final surface S (psu)",     :haline),
                                      (wm,  "final surface w (m/s)",     :balance),
                                      (spd, "final surface speed (m/s)", :speed)))
    ax = Axis(fig[1, i]; title = lab, xlabel = "λ", ylabel = "φ", aspect = DataAspect())
    hm = heatmap!(ax, λ, φ, d; nan_color = :gray70, colormap = cmap)
    Colorbar(fig[2, i], hm; vertical = false, flipaxis = false)
end
out = "glorys_run_final.png"
save(out, fig)
@info "Saved final-frame maps -> $out (look for the S overshoot / w spike location)"
