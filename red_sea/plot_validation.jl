# plot_validation.jl -- summarise the 7-day validation run: the four monitored time series from the
# run log, plus the final surface and Bab el-Mandeb sections from the JLD2 output.
# Usage: julia --project=<repo> red_sea/plot_validation.jl [logfile] [tag]
#   tag selects the JLD2 output set and names the figure, e.g. "_zstar_jraryf" for the forced run.

using Oceananigans, CairoMakie, Printf
using Oceananigans.Fields: location
using Oceananigans.ImmersedBoundaries: immersed_inactive_node

logfile = length(ARGS) > 0 ? ARGS[1] : joinpath(@__DIR__, "redsea_validation_7day.log")
tag     = length(ARGS) > 1 ? ARGS[2] : "_zstar"

# --- time series from the log ------------------------------------------------
secs = Dict("second"=>1.0,"seconds"=>1.0,"minute"=>60,"minutes"=>60,"hour"=>3600,"hours"=>3600,"day"=>86400,"days"=>86400)
t = Float64[]; η = Float64[]; umax = Float64[]
ts = Float64[]; meanS = Float64[]; salt = Float64[]; vol = Float64[]; ninv = Float64[]; drho = Float64[]
tq = Float64[]; qin = Float64[]; qout = Float64[]; dS = Float64[]
for ln in eachline(logfile)
    m = match(r"iter \d+ \| t = ([\d.]+) (\w+).*max\|u,v,w\| = \(([\d.e+-]+), [\d.e+-]+, [\d.e+-]+\) \| max\|η\| = ([\d.e+-]+)", ln)
    if m !== nothing
        push!(t, parse(Float64, m[1]) * secs[m[2]] / 86400)
        push!(umax, parse(Float64, m[3])); push!(η, parse(Float64, m[4])); continue
    end
    m = match(r"conservation: salt ([+\-\d.e]+) \| volume ([+\-\d.e]+).*mean S ([\d.]+) \| unstable pairs (\d+) \(worst Δρ ([\d.]+)\)", ln)
    if m !== nothing && !isempty(t)
        push!(ts, t[end]); push!(salt, parse(Float64, m[1])); push!(vol, parse(Float64, m[2]))
        push!(meanS, parse(Float64, m[3])); push!(ninv, parse(Float64, m[4])); push!(drho, parse(Float64, m[5])); continue
    end
    m = match(r"BAB EL-MANDEB: Qin ([+\-\d.e]+)\s+Qout ([+\-\d.e]+).*?ΔS ([+\-\d.]+)", ln)
    if m !== nothing && !isempty(t)
        push!(tq, t[end]); push!(qin, parse(Float64, m[1])); push!(qout, parse(Float64, m[2])); push!(dS, parse(Float64, m[3]))
    end
end
@info "parsed $(length(t)) progress blocks, $(length(tq)) strait samples, t = 0 .. $(round(t[end], digits=3)) d"

fig = Figure(size = (1600, 1050))

ax1 = Axis(fig[1,1]; title = "salt and volume: relative drift", xlabel = "t (days)", ylabel = "drift")
lines!(ax1, ts, salt; label = "salt", linewidth = 2)
lines!(ax1, ts, vol;  label = "volume", linewidth = 2, linestyle = :dash)
hlines!(ax1, [0.0]; color = :gray, linestyle = :dot); axislegend(ax1; position = :lb)

ax2 = Axis(fig[1,2]; title = "mean salinity (psu) — the z* drift test", xlabel = "t (days)", ylabel = "mean S")
lines!(ax2, ts, meanS; linewidth = 2, color = :firebrick)

ax3 = Axis(fig[1,3]; title = "max|η| (m), log scale — adjustment transient then settle",
           xlabel = "t (days)", ylabel = "max|η|", yscale = log10)
lines!(ax3, t, η; linewidth = 2, color = :teal)

ax4 = Axis(fig[2,1]; title = "static instability: count and worst Δρ", xlabel = "t (days)", ylabel = "unstable pairs")
lines!(ax4, ts, ninv; linewidth = 2, color = :darkorange)
ax4b = Axis(fig[2,1]; ylabel = "worst Δρ (kg/m³)", yaxisposition = :right, yticklabelcolor = :purple)
hidespines!(ax4b); hidexdecorations!(ax4b); linkxaxes!(ax4, ax4b)
lines!(ax4b, ts, drho; linewidth = 2, color = :purple, linestyle = :dash)

ax5 = Axis(fig[2,2]; title = "Bab el-Mandeb exchange (m³/s)", xlabel = "t (days)", ylabel = "transport")
lines!(ax5, tq, qin;  linewidth = 2, label = "Qin (into Red Sea)")
lines!(ax5, tq, qout; linewidth = 2, label = "Qout")
hlines!(ax5, [0.0]; color = :gray, linestyle = :dot); axislegend(ax5; position = :rc)

ax6 = Axis(fig[2,3]; title = "Bab el-Mandeb ΔS = Sin − Sout\n(negative = inverse estuary)", xlabel = "t (days)", ylabel = "ΔS (psu)")
lines!(ax6, tq, dS; linewidth = 2, color = :seagreen)
hlines!(ax6, [0.0]; color = :gray, linestyle = :dot); ylims!(ax6, -2, 2)

# --- final fields ------------------------------------------------------------
sfile = joinpath(@__DIR__, "redsea_glorys_surface$(tag).jld2")
if isfile(sfile)
    Ss = FieldTimeSeries(sfile, "S")
    g = Ss.grid; Nx, Ny, Nz = size(g)
    n = length(Ss.times)
    wet = [immersed_inactive_node(i, j, Nz, g, Center(), Center(), Center()) ? NaN : 1.0 for i in 1:Nx, j in 1:Ny]
    λ = λnodes(g, Center(), Center(), Center()); φ = φnodes(g, Center(), Center(), Center())
    for (col, idx, lab) in ((1, 1, "day 0"), (2, n, @sprintf("day %.1f", Ss.times[n]/86400)))
        d = Array(interior(Ss[idx]))[1:Nx, 1:Ny, 1] .* wet
        ax = Axis(fig[3, col]; title = "surface S, $lab", xlabel = "λ", ylabel = "φ", aspect = DataAspect())
        hm = heatmap!(ax, λ, φ, d; nan_color = :gray70, colormap = :haline, colorrange = (35, 41))
        col == 2 && Colorbar(fig[3, 3, Left()], hm; width = 12, label = "S (psu)")
    end
end
outpng = joinpath(@__DIR__, "redsea_validation_7day$(tag == "_zstar" ? "" : tag).png")
save(outpng, fig)
@info "Saved $outpng"
