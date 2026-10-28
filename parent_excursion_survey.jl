# How much does the boundary excursion -- and so the TracerReservoir length scale L -- actually
# vary between domains? Estimated from the PARENT alone, which is what a configure-time utility
# would have available, rather than from child output.
#
# Result, GLORYS daily 1993-01-01..10, unweighted median over wet cells:
#
#   Red Sea east 53.9 km | Canary west 28.8, south 39.0, north 39.9 | Salish west 14.9, north 12.4
#
# Medians span 4.3x, so one fixed L gives an 8x spread in retention exp(-E/L) -- 6.7% at Red Sea
# against 54% at Salish north for the 20 km nested_ocean_model default. CAVEAT: this parent-only
# estimate is biased in BOTH directions against the child-derived value (Red Sea 53.9 vs 16.2 km
# volume-weighted, 3.3x high; Salish west 14.9 vs 17.6, 15% low). Where the parent resolves the
# boundary well, daily sampling aliases sub-daily reversals and the integral accumulates one-way
# displacement; where it resolves it badly (Salish: 9 wet cells of 21) it smears the flow toward
# zero. So the estimate degrades exactly where the nesting is already marginal. Note also that the
# median here is unweighted while the child figures are volume-weighted, and this is one January
# window. The solution's SENSITIVITY to L remains unmeasured -- see red_sea/README.md.
using NumericalEarth, Oceananigans, Printf, Statistics, Dates
using CopernicusMarine     # registers the GLORYS download backend; without it, set! throws
const DW = NumericalEarth.DataWrangling

glorys = DW.GLORYSDaily()

# (name, region, boundary-normal component, boundary coordinate, axis the boundary runs along)
cases = [
 ("Red Sea   east",  DW.BoundingBox(longitude=(31.5,50.5),   latitude=(9.5,30.5)),  :u_velocity, 50.0,  :lon),
 ("Canary    west",  DW.BoundingBox(longitude=(-22.0,-9.5),  latitude=(24.5,30.0)), :u_velocity, -21.5, :lon),
 ("Canary    south", DW.BoundingBox(longitude=(-22.0,-9.5),  latitude=(24.5,30.0)), :v_velocity, 25.0,  :lat),
 ("Canary    north", DW.BoundingBox(longitude=(-22.0,-9.5),  latitude=(24.5,30.0)), :v_velocity, 29.5,  :lat),
 ("Salish    west",  DW.BoundingBox(longitude=(235.2,238.2), latitude=(46.8,50.0)), :u_velocity, 235.5, :lon),
 ("Salish    north", DW.BoundingBox(longitude=(235.2,238.2), latitude=(46.8,50.0)), :v_velocity, 49.05, :lat),
]

dates = DateTime(1993,1,1) : Day(1) : DateTime(1993,1,10)
L_default = 20e3           # the nested_ocean_model default, for the retention column

function excursion(name, region, var, coord, axis)
    meta = DW.Metadata(var; dataset = glorys, dates, region)
    fts  = FieldTimeSeries(meta, CPU(); time_indices_in_memory = length(dates),
                           inpainting = DW.NearestNeighborInpainting(10))
    g = fts.grid
    t = Float64.(collect(fts.times))
    Nx, Ny, Nz = size(g)
    data = [Array(interior(fts[n])) for n in 1:length(t)]

    # The boundary line, and the range running along it.
    if axis === :lon
        b = argmin(abs.(Array(λnodes(g, Center())) .- coord)); along = 1:Ny
    else
        b = argmin(abs.(Array(φnodes(g, Center())) .- coord)); along = 1:Nx
    end

    exc = Float64[]
    for k in 1:Nz, a in along
        series = axis === :lon ? [data[n][b, a, k] for n in 1:length(t)] :
                                 [data[n][a, b, k] for n in 1:length(t)]
        all(isfinite, series) || continue
        maximum(abs, series) < 1e-6 && continue          # land / stagnant
        # Trapezoidal displacement, tracking its peak-to-peak range.
        x = 0.0; lo = 0.0; hi = 0.0
        for n in 2:length(t)
            x += 0.5 * (series[n] + series[n-1]) * (t[n] - t[n-1])
            lo = min(lo, x); hi = max(hi, x)
        end
        push!(exc, hi - lo)
    end
    isempty(exc) && return @printf("%-16s  no wet cells\n", name)
    E = median(exc)
    @printf("%-16s  n=%5d   median %6.1f km   mean %6.1f km   p90 %6.1f km   retains %4.1f%% at L=%.0f km\n",
            name, length(exc), E/1e3, mean(exc)/1e3, quantile(exc, 0.9)/1e3,
            100exp(-E/L_default), L_default/1e3)
end

println("Parent-derived boundary excursion, GLORYS daily, 1993-01-01..10")
println("(peak-to-peak displacement from time-integrated boundary-normal velocity)\n")
for c in cases
    try
        excursion(c...)
    catch e
        @printf("%-16s  FAILED: %s\n", c[1], sprint(showerror, e)[1:min(70,end)])
    end
end
