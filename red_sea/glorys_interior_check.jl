# glorys_interior_check.jl
#
# Does GLORYS actually have valid data everywhere the Red Sea model has water?
#
# The boundary box was checked before the design was committed (README, "GLORYS validity at the
# boundary"): valid to 3226 m against a 3176 m domain. But a box minimum says nothing about any
# individual column, and the Salish Sea damage was done in the deep INTERIOR, not at the boundary.
# The failure mode is per-column: wherever the parent's own bathymetry is shallower than the child's,
# the child has wet cells the parent never resolved.
#
# Three stages, in increasing order of what actually matters:
#   1. Coverage    -- per model column, is there valid parent data down to the model floor?
#   2. Provenance  -- where there is not, what does the pipeline (inpaint -> fill_down!) put there,
#                     and is it the right water mass? Inpainting reaches up to 100 cells sideways, so
#                     it always returns a finite number; the question is whether it borrowed that
#                     number from across a sill.
#   3. Stability   -- does the resulting IC contain a density inversion? This is the one that counts:
#                     odd values did not destabilise the Salish Sea, an inversion did.
#
# Read-only. Builds the same grid and bathymetry as redsea_glorys_obc.jl and runs offline off the
# GLORYS cache (~65 s at 0.05 deg). Honours REDSEA_DELTA / REDSEA_PASSES so it can be re-run at the
# run resolution; REDSEA_CHECK_CONTROL=1 adds the fill_down!-disabled control in stage 3.

using Oceananigans
using NumericalEarth
using CopernicusMarine
using CairoMakie
using Dates
using Printf
using Statistics
using Oceananigans.Fields: location
using Oceananigans.ImmersedBoundaries: immersed_inactive_node
using Oceananigans.Units: Time
const SWP = Oceananigans.BuoyancyFormulations.SeawaterPolynomials

const DW = NumericalEarth.DataWrangling

arch = CPU()
check_date = DateTime(1993, 1, 1)      # run_start in the model script

# --- Model grid, identical to redsea_glorys_obc.jl ----------------------------
λ₁, λ₂ = (32.0, 50.0)
φ₁, φ₂ = (10.0, 30.0)
Δdeg = parse(Float64, get(ENV, "REDSEA_DELTA", "0.05"))
Nx = round(Int, (λ₂ - λ₁) / Δdeg)
Ny = round(Int, (φ₂ - φ₁) / Δdeg)
Nz = 32
depth = 3200
z_scale = 480
z_reference = ExponentialDiscretization(Nz, -depth, 0; scale = z_scale)

interpolation_passes = parse(Int, get(ENV, "REDSEA_PASSES", "5"))

@info @sprintf("Model grid: %d x %d x %d (Δ = %.4f deg), %d smoothing passes", Nx, Ny, Nz, Δdeg, interpolation_passes)

underlying_grid = LatitudeLongitudeGrid(arch;
                             size = (Nx, Ny, Nz),
                             latitude = (φ₁, φ₂), longitude = (λ₁, λ₂),
                             z = z_reference, halo = (7, 7, 7),
                             topology = (Bounded, Bounded, Bounded))
bottom_height = regrid_bathymetry(underlying_grid; minimum_depth = 10,
                                  interpolation_passes = interpolation_passes)
grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bottom_height))

hb = Array(interior(bottom_height))[:, :, 1]        # negative, 0 on land
λc = λnodes(grid, Center(), Center(), Center())
φc = φnodes(grid, Center(), Center(), Center())
zc = znodes(grid, Center(), Center(), Center())

wet2d  = [!immersed_inactive_node(i, j, Nz, grid, Center(), Center(), Center()) for i in 1:Nx, j in 1:Ny]
wet3d  = [!immersed_inactive_node(i, j, k, grid, Center(), Center(), Center()) for i in 1:Nx, j in 1:Ny, k in 1:Nz]
wetcols = findall(wet2d)
@info @sprintf("Wet columns: %d of %d (%.1f%%) | wet cells: %d | deepest column %.0f m",
               length(wetcols), Nx*Ny, 100length(wetcols)/(Nx*Ny), count(wet3d), -minimum(hb))

# --- GLORYS native, UN-inpainted: the true validity mask ----------------------
glorys_region = DW.BoundingBox(longitude = (31.5, 50.5), latitude = (9.5, 30.5))
glorys = DW.GLORYSDaily()
raw_meta(var) = DW.Metadatum(var; dataset = glorys, region = glorys_region, date = check_date)

S_raw = Field(raw_meta(:salinity); inpainting = nothing)
T_raw = Field(raw_meta(:temperature); inpainting = nothing)
Sr = Array(interior(S_raw))
Tr = Array(interior(T_raw))
gλ = λnodes(S_raw.grid, map(L -> L(), location(S_raw))...)
gφ = φnodes(S_raw.grid, map(L -> L(), location(S_raw))...)
gz = znodes(S_raw.grid, map(L -> L(), location(S_raw))...)   # index 1 = deepest, end = surface
gzf = znodes(S_raw.grid, Center(), Center(), Face())         # cell faces; bottom face of level k is gzf[k]
ngx, ngy, ngz = size(Sr)
@info @sprintf("GLORYS native: %d x %d x %d | λ %.3f..%.3f  φ %.3f..%.3f | z %.1f..%.1f m",
               ngx, ngy, ngz, first(gλ), last(gλ), first(gφ), last(gφ), first(gz), last(gz))

# Deepest valid level per native column. NaN where the whole column is invalid (land in GLORYS).
# Two floors per parent column, because they answer different questions:
#   gbot  = deepest valid CENTRE. Below it interpolate_to_grid has nothing to interpolate between,
#           so this is the depth below which the child's values are extrapolated by fill_down!.
#   gbotf = bottom FACE of that same cell. This is how deep the parent's water actually goes, and is
#           the fair comparison for "does the parent resolve this column at all". Level spacing is
#           ~60 m at 300 m and ~350 m at 3 km, so the two differ by a lot in the deep.
gbot  = fill(NaN, ngx, ngy)
gbotf = fill(NaN, ngx, ngy)
gktop = zeros(Int, ngx, ngy)     # index of deepest valid level (0 = none)
for j in 1:ngy, i in 1:ngx
    k = findfirst(k -> isfinite(Sr[i, j, k]), 1:ngz)
    if k !== nothing
        gktop[i, j] = k
        gbot[i, j]  = gz[k]
        gbotf[i, j] = gzf[k]
    end
end
nwet_native = count(>(0), gktop)
@info @sprintf("GLORYS wet columns in region: %d of %d (%.1f%%) | deepest valid level over region: %.1f m",
               nwet_native, ngx*ngy, 100nwet_native/(ngx*ngy),
               minimum(gbot[isfinite.(gbot)]))

# Nearest native column for a model column (native grid is regular at 1/12 deg).
gΔλ = (last(gλ) - first(gλ)) / (ngx - 1)
gΔφ = (last(gφ) - first(gφ)) / (ngy - 1)
nearest_gi(λ) = clamp(round(Int, (λ - first(gλ)) / gΔλ) + 1, 1, ngx)
nearest_gj(φ) = clamp(round(Int, (φ - first(gφ)) / gΔφ) + 1, 1, ngy)

# --- Stage 1: coverage --------------------------------------------------------
# margin > 0  : parent has valid data below the model floor            -> fine
# margin <= 0 : model has wet cells below the parent's deepest datum   -> extrapolated or borrowed
# "dry parent": the nearest parent column is land entirely             -> the worst case
margin   = fill(NaN, Nx, Ny)     # against the deepest valid centre  (extrapolation onset)
marginf  = fill(NaN, Nx, Ny)     # against that cell's bottom face    (true parent coverage)
dryparent = falses(Nx, Ny)
for I in wetcols
    i, j = Tuple(I)
    gi, gj = nearest_gi(λc[i]), nearest_gj(φc[j])
    if gktop[gi, gj] == 0
        dryparent[i, j] = true
    else
        margin[i, j]  = hb[i, j] - gbot[gi, gj]     # both negative; positive => parent deeper
        marginf[i, j] = hb[i, j] - gbotf[gi, gj]
    end
end

# Cells (not columns) below the nearest parent column's deepest valid datum.
below = falses(Nx, Ny, Nz)
for k in 1:Nz, I in wetcols
    i, j = Tuple(I)
    wet3d[i, j, k] || continue
    gi, gj = nearest_gi(λc[i]), nearest_gj(φc[j])
    below[i, j, k] = gktop[gi, gj] == 0 || zc[k] < gbot[gi, gj]
end

# Region split: the Gulf of Aden lies entirely east of Bab el-Mandeb (~43.3E), so longitude alone
# separates the two basins.
redsea(i, j) = λc[i] < 43.4
region_of(i, j) = redsea(i, j) ? "Red Sea" : "Gulf of Aden"

fin  = [margin[I]  for I in wetcols if isfinite(margin[I])]
finf = [marginf[I] for I in wetcols if isfinite(marginf[I])]
nviol  = count(<=(0), fin)
nviolf = count(<=(0), finf)
ndry  = count(dryparent)
println("\n", "="^92)
println("STAGE 1 -- COVERAGE: model floor vs deepest valid GLORYS level (nearest parent column)")
println("="^92)
@printf("  wet model columns              : %d\n", length(wetcols))
@printf("  parent column entirely dry     : %d  (%.2f%%)\n", ndry, 100ndry/length(wetcols))
@printf("  below deepest valid CENTRE     : %d  (%.2f%%)   <- extrapolated by fill_down!\n",
        nviol, 100nviol/length(wetcols))
@printf("     quantiles (m)               : min %.0f | 1%% %.0f | 5%% %.0f | 50%% %.0f | max %.0f\n",
        minimum(fin), quantile(fin, 0.01), quantile(fin, 0.05), median(fin), maximum(fin))
@printf("  below parent's BOTTOM FACE     : %d  (%.2f%%)   <- parent has no water here at all\n",
        nviolf, 100nviolf/length(wetcols))
@printf("     quantiles (m)               : min %.0f | 1%% %.0f | 5%% %.0f | 50%% %.0f | max %.0f\n",
        minimum(finf), quantile(finf, 0.01), quantile(finf, 0.05), median(finf), maximum(finf))
@printf("  wet model cells below parent   : %d of %d (%.2f%%)\n",
        count(below), count(wet3d), 100count(below)/count(wet3d))

for (lab, sel) in (("Red Sea (λ<43.4)", (i,j) -> redsea(i,j)), ("Gulf of Aden (λ>43.4)", (i,j) -> !redsea(i,j)))
    sub  = [margin[I]  for I in wetcols if (sel(Tuple(I)...) && isfinite(margin[I]))]
    subf = [marginf[I] for I in wetcols if (sel(Tuple(I)...) && isfinite(marginf[I]))]
    subdry = count(I -> sel(Tuple(I)...) && dryparent[I], wetcols)
    nb = count(I -> sel(Tuple(I)[1], Tuple(I)[2]), findall(below))
    nw = count(I -> sel(Tuple(I)[1], Tuple(I)[2]), findall(wet3d))
    isempty(sub) && continue
    @printf("  %-22s: n=%6d | dry parent %4d | below centre %5d | below face %5d | min face margin %+7.0f m | cells below %d/%d (%.2f%%)\n",
            lab, length(sub), subdry, count(<=(0), sub), count(<=(0), subf), minimum(subf), nb, nw, 100nb/nw)
end

# Worst offenders, deduplicated onto a coarse grid so the list spans the basin instead of one feature.
println("\n  Worst 15 columns (most negative margin), 0.25-deg dedup:")
cand = sort([(margin[I], Tuple(I)) for I in wetcols if isfinite(margin[I]) && margin[I] <= 0], by = first)
seen = Set{Tuple{Int,Int}}(); shown = Ref(0)
for (m, (i, j)) in cand
    key = (round(Int, λc[i]/0.25), round(Int, φc[j]/0.25))
    key in seen && continue
    push!(seen, key); shown[] += 1
    gi, gj = nearest_gi(λc[i]), nearest_gj(φc[j])
    @printf("    %7.3fE %6.3fN  %-13s model floor %6.0f m | parent centre %6.0f m / face %6.0f m | margin %+6.0f m\n",
            λc[i], φc[j], region_of(i, j), -hb[i, j], -gbot[gi, gj], -gbotf[gi, gj], m)
    shown[] >= 15 && break
end
shown[] == 0 && println("    (none -- every wet model column has valid parent data below its floor)")

if ndry > 0
    println("\n  Wet model columns whose nearest parent column is entirely land, 0.25-deg dedup:")
    seen2 = Set{Tuple{Int,Int}}(); shown2 = Ref(0)
    for I in wetcols
        dryparent[I] || continue
        i, j = Tuple(I)
        key = (round(Int, λc[i]/0.25), round(Int, φc[j]/0.25))
        key in seen2 && continue
        push!(seen2, key); shown2[] += 1
        @printf("    %7.3fE %6.3fN  %-13s model floor %6.0f m\n", λc[i], φc[j], region_of(i, j), -hb[i, j])
        shown2[] >= 15 && break
    end
end

# --- Stage 2: provenance ------------------------------------------------------
# What the model actually initialises with: the SAME pipeline as redsea_glorys_obc.jl --
# NearestNeighborInpainting(100), then fill_down!, then interpolate_to_grid.
bc_inpainting = DW.NearestNeighborInpainting(100)

function fill_down!(a::AbstractArray{<:Any,3}, floor_value)
    nx, ny, nz = size(a)
    @inbounds for j in 1:ny, i in 1:nx
        last = NaN
        for k in nz:-1:1
            v = a[i, j, k]
            if isfinite(v) && v > floor_value
                last = v
            elseif isfinite(last)
                a[i, j, k] = last
            end
        end
    end
    return a
end

# Two dates, not one: a single-element FieldTimeSeries indexes out of bounds under Cyclical time
# indexing. Only index 1 is used.
fts_dates = [check_date, check_date + Day(1)]
fts_meta(var) = DW.Metadata(var; dataset = glorys, dates = fts_dates, region = glorys_region)
S_fts = FieldTimeSeries(fts_meta(:salinity), arch; time_indices_in_memory = 2, inpainting = bc_inpainting)
T_fts = FieldTimeSeries(fts_meta(:temperature), arch; time_indices_in_memory = 2, inpainting = bc_inpainting)
S_native, T_native = S_fts.grid, T_fts.grid

Sp = interior(S_fts[1]); Tp = interior(T_fts[1])
nzero_S = count(==(0), Sp); nzero_T = count(==(0), Tp)
@info @sprintf("After inpainting, before fill_down!: S zeros %d (%.3f%%), T zeros %d (%.3f%%)",
               nzero_S, 100nzero_S/length(Sp), nzero_T, 100nzero_T/length(Tp))
fill_down!(Sp, 5.0); fill_down!(Tp, 1.0)
fill_down!(interior(S_fts[2]), 5.0); fill_down!(interior(T_fts[2]), 1.0)
nzero_S2 = count(==(0), Sp); nzero_T2 = count(==(0), Tp)
@info @sprintf("After fill_down!:                    S zeros %d (%.3f%%), T zeros %d (%.3f%%)",
               nzero_S2, 100nzero_S2/length(Sp), nzero_T2, 100nzero_T2/length(Tp))

# Interpolate onto the model grid exactly as ic_from_fts! does.
S_ic = zeros(Nx, Ny, Nz); T_ic = zeros(Nx, Ny, Nz)
t₀ = Time(0.0)
for k in 1:Nz, j in 1:Ny, i in 1:Nx
    S_ic[i, j, k] = DW.interpolate_to_grid(S_fts, i, j, k, S_native, grid, t₀)
    T_ic[i, j, k] = DW.interpolate_to_grid(T_fts, i, j, k, T_native, grid, t₀)
end

println("\n", "="^92)
println("STAGE 2 -- PROVENANCE: what the IC pipeline puts in the cells stage 1 flagged")
println("="^92)

widx = findall(wet3d)
Sw = [S_ic[I] for I in widx]; Tw = [T_ic[I] for I in widx]
@printf("  IC over wet cells : T ∈ [%.3f, %.3f] | S ∈ [%.3f, %.3f]\n",
        minimum(Tw), maximum(Tw), minimum(Sw), maximum(Sw))

nonphys = [I for I in widx if S_ic[I] < 30 || T_ic[I] < 0 || T_ic[I] > 40]
@printf("  non-physical wet cells (S<30, T<0 or T>40): %d\n", length(nonphys))
if !isempty(nonphys)
    println("    first 10:")
    for I in nonphys[1:min(10, end)]
        i, j, k = Tuple(I)
        @printf("      %7.3fE %6.3fN  z=%7.1f m  %-13s T=%.3f S=%.3f\n",
                λc[i], φc[j], zc[k], region_of(i, j), T_ic[I], S_ic[I])
    end
end

# Measured vs extrapolated. `below` marks cells under the nearest parent column's deepest valid
# CENTRE -- exactly the cells whose value came from fill_down! rather than from parent data. If the
# two populations have the same distribution, extrapolation is doing no harm; if they separate, the
# interior IC is carrying a water mass the parent never had there.
println()
for (lab, sel) in (("whole domain", (i,j) -> true),
                   ("Red Sea",      (i,j) -> redsea(i,j)),
                   ("Gulf of Aden", (i,j) -> !redsea(i,j)))
    meas = [I for I in widx if sel(Tuple(I)[1], Tuple(I)[2]) && !below[I]]
    extr = [I for I in widx if sel(Tuple(I)[1], Tuple(I)[2]) &&  below[I]]
    tot = length(meas) + length(extr)
    tot == 0 && continue
    @printf("  %-13s  extrapolated %6d / %6d cells (%5.2f%%)\n", lab, length(extr), tot, 100length(extr)/tot)
    for (tag, idx) in (("measured", meas), ("extrapolated", extr))
        isempty(idx) && continue
        sv = [S_ic[I] for I in idx]; tv = [T_ic[I] for I in idx]
        @printf("      %-13s S %.3f..%.3f (med %.3f) | T %.3f..%.3f (med %.3f)\n",
                tag, minimum(sv), maximum(sv), median(sv), minimum(tv), maximum(tv), median(tv))
    end
end

# Water-mass test. Red Sea Deep Water below the Hanish sill (136 m) is ~40.5-40.6 psu at ~21.5 degC
# and is nearly uniform -- the Red Sea is warm and salty all the way down. Gulf of Aden water at the
# same depth is ~35.3-36.4 and cold. If inpainting has borrowed across the sill, deep Red Sea cells
# will carry Gulf of Aden properties -- an error in the INTERIOR, where nothing nudges it back.
println()
deep_rs = [I for I in widx if (λc[Tuple(I)[1]] < 43.4 && φc[Tuple(I)[2]] > 15.0 && zc[Tuple(I)[3]] < -300)]
deep_ga = [I for I in widx if (λc[Tuple(I)[1]] > 45.0 && zc[Tuple(I)[3]] < -300)]
for (lab, idx) in (("Red Sea deep (λ<43.4, φ>15N, z<-300m)", deep_rs),
                   ("Gulf of Aden deep (λ>45E, z<-300m)",    deep_ga))
    isempty(idx) && continue
    sv = [S_ic[I] for I in idx]; tv = [T_ic[I] for I in idx]
    @printf("  %-40s n=%6d | S %.3f..%.3f (med %.3f) | T %.3f..%.3f (med %.3f)\n",
            lab, length(idx), minimum(sv), maximum(sv), median(sv), minimum(tv), maximum(tv), median(tv))
end

# Deep Red Sea cells carrying Gulf of Aden salinity: the specific failure this stage is looking for.
# T is the sharper discriminator: Red Sea deep water is ~21.5 degC, Gulf of Aden deep water ~12.
suspect = [I for I in deep_rs if S_ic[I] < 39.0 || T_ic[I] < 18.0]
@printf("  deep Red Sea cells with S < 39.0 or T < 18.0 (not Red Sea deep water): %d of %d (%.2f%%)\n",
        length(suspect), length(deep_rs), 100length(suspect)/length(deep_rs))
if !isempty(suspect)
    ss = sort(suspect, by = I -> S_ic[I])
    println("    freshest 10 (extrapolated? = value came from fill_down!, not parent data):")
    for I in ss[1:min(10, end)]
        i, j, k = Tuple(I)
        @printf("      %7.3fE %6.3fN  z=%7.1f m  T=%.3f S=%.3f  | model floor %.0f m, parent centre %.0f m | %s\n",
                λc[i], φc[j], zc[k], T_ic[I], S_ic[I], -hb[i, j],
                -gbot[nearest_gi(λc[i]), nearest_gj(φc[j])],
                below[I] ? "EXTRAPOLATED" : "from parent data")
    end
end

# --- Stage 3: static stability ------------------------------------------------
# The Salish Sea configuration did not blow up because its values were merely odd -- it blew up
# because the IC contained a DENSITY INVERSION (S = 30.16 under S = 33.92, 14 psu of fresh water
# beneath salty water, convecting). That is the defect to test for directly. Everything above is
# only a way of predicting where to look.
#
# Note TEOS-10 wants conservative temperature and absolute salinity while GLORYS supplies potential
# temperature and practical salinity. The model makes the same substitution, so this measures the
# stability the model will actually see.
eos = SWP.TEOS10.TEOS10EquationOfState()
zf = znodes(grid, Center(), Center(), Face())

# Compare adjacent wet cells at the shared interface depth: an inversion is denser water sitting on
# lighter water, i.e. ρ(cell above) > ρ(cell below) referenced to the same pressure.
function inversions(T_ic, S_ic)
    out = Tuple{Float64,Int,Int,Int}[]     # (Δρ, i, j, k) with k the LOWER cell
    for j in 1:Ny, i in 1:Nx
        for k in 1:(Nz-1)
            (wet3d[i, j, k] && wet3d[i, j, k+1]) || continue
            zr = zf[k+1]
            ρl = SWP.ρ(T_ic[i, j, k],   S_ic[i, j, k],   zr, eos)
            ρu = SWP.ρ(T_ic[i, j, k+1], S_ic[i, j, k+1], zr, eos)
            Δρ = ρu - ρl
            Δρ > 1e-4 && push!(out, (Δρ, i, j, k))
        end
    end
    return sort(out, by = first, rev = true)
end

println("\n", "="^92)
println("STAGE 3 -- STATIC STABILITY: does the IC contain the defect that destabilised the Salish Sea?")
println("="^92)

inv = inversions(T_ic, S_ic)
npairs = count(I -> (Tuple(I)[3] < Nz && wet3d[Tuple(I)[1], Tuple(I)[2], Tuple(I)[3]+1]), widx)
@printf("  unstable adjacent pairs: %d of %d wet pairs (%.4f%%)\n", length(inv), npairs, 100length(inv)/npairs)
if !isempty(inv)
    for (lo, hi) in ((1e-4, 1e-3), (1e-3, 1e-2), (1e-2, 1e-1), (1e-1, 1.0), (1.0, Inf))
        n = count(x -> lo <= x[1] < hi, inv)
        n > 0 && @printf("     Δρ ∈ [%.0e, %.0e) : %6d\n", lo, hi, n)
    end
    isext(x) = below[x[2], x[3], x[4]] || below[x[2], x[3], x[4]+1]
    nx_ext = count(isext, inv)
    @printf("  of those, %d (%.1f%%) involve a cell below the parent's deepest valid centre; %d sit entirely in parent data\n",
            nx_ext, 100nx_ext/length(inv), length(inv) - nx_ext)

    # The tail is what matters. Small inversions are ordinary reanalysis noise and the closure mixes
    # them away in the first steps; the Salish Sea blow-up was ~10 kg/m³ (14 psu of fresh under salty).
    tail = filter(x -> x[1] >= 0.1, inv)
    @printf("\n  TAIL, Δρ ≥ 0.1 kg/m³: %d pairs (%.4f%% of wet pairs)\n", length(tail), 100length(tail)/npairs)
    if !isempty(tail)
        @printf("     below parent's deepest valid centre : %d (%.1f%%)\n",
                count(isext, tail), 100count(isext, tail)/length(tail))
        @printf("     Red Sea / Gulf of Aden              : %d / %d\n",
                count(x -> redsea(x[2], x[3]), tail), count(x -> !redsea(x[2], x[3]), tail))
        @printf("     above the Hanish sill (z > -136 m)  : %d   below it: %d\n",
                count(x -> zc[x[4]] > -136, tail), count(x -> zc[x[4]] <= -136, tail))
        @printf("     depth of lower cell (m)             : %.0f .. %.0f (median %.0f)\n",
                -maximum(x -> -zc[x[4]], tail), -minimum(x -> -zc[x[4]], tail),
                -median([-zc[x[4]] for x in tail]))
    end

    println("\n  worst 10:")
    for (Δρ, i, j, k) in inv[1:min(10, end)]
        gi, gj = nearest_gi(λc[i]), nearest_gj(φc[j])
        @printf("    %7.3fE %6.3fN  z=%7.1f/%7.1f m  %-13s Δρ=%+.4f | upper T=%.2f S=%.3f / lower T=%.2f S=%.3f | parent valid to %.0f m %s\n",
                λc[i], φc[j], zc[k+1], zc[k], region_of(i, j), Δρ,
                T_ic[i, j, k+1], S_ic[i, j, k+1], T_ic[i, j, k], S_ic[i, j, k],
                gktop[gi, gj] == 0 ? NaN : -gbot[gi, gj],
                isext((Δρ, i, j, k)) ? "[below parent]" : "")
    end
else
    println("  none -- the initial condition is statically stable everywhere.")
end

# Control: the same count with fill_down! switched off, which is the state that destabilised the
# Salish Sea. Opt-in (REDSEA_CHECK_CONTROL=1) because it re-interpolates the whole domain.
if get(ENV, "REDSEA_CHECK_CONTROL", "0") == "1"
    S_fts0 = FieldTimeSeries(fts_meta(:salinity), arch; time_indices_in_memory = 2, inpainting = bc_inpainting)
    T_fts0 = FieldTimeSeries(fts_meta(:temperature), arch; time_indices_in_memory = 2, inpainting = bc_inpainting)
    S_ic0 = zeros(Nx, Ny, Nz); T_ic0 = zeros(Nx, Ny, Nz)
    for k in 1:Nz, j in 1:Ny, i in 1:Nx
        S_ic0[i, j, k] = DW.interpolate_to_grid(S_fts0, i, j, k, S_fts0.grid, grid, t₀)
        T_ic0[i, j, k] = DW.interpolate_to_grid(T_fts0, i, j, k, T_fts0.grid, grid, t₀)
    end
    Sw0 = [S_ic0[I] for I in widx]
    inv0 = inversions(T_ic0, S_ic0)
    @printf("\n  CONTROL, fill_down! disabled: wet-cell S ∈ [%.3f, %.3f] | unstable pairs %d (%.4f%%) | worst Δρ %.3f kg/m³\n",
            minimum(Sw0), maximum(Sw0), length(inv0), 100length(inv0)/npairs,
            isempty(inv0) ? 0.0 : inv0[1][1])
end

# --- Figure -------------------------------------------------------------------
nan_on_land(a) = [wet2d[i, j] ? Float64(a[i, j]) : NaN for i in 1:Nx, j in 1:Ny]

# Worst inversion per column, so the stage-3 tail is locatable on a map rather than a list.
worst_inv = zeros(Nx, Ny)
for (Δρ, i, j, k) in inv
    worst_inv[i, j] = max(worst_inv[i, j], Δρ)
end

# S at the sill depth, not at the deepest level: the Hanish sill is 136 m, only 0.7% of columns reach
# 2500 m, and the tail of stage 3 sits in the top 200 m.
k_sill = argmin(abs.(zc .+ 150))
S_sill = [wet3d[i, j, k_sill] ? S_ic[i, j, k_sill] : NaN for i in 1:Nx, j in 1:Ny]

panels = ((nan_on_land(marginf), "model floor − parent bottom face (m)\nblue = parent has no water this deep", :balance, (-400, 400)),
          (nan_on_land(worst_inv), "worst IC density inversion per column (kg/m³)\n0 = stable", :inferno, (0, 0.3)),
          (nan_on_land([S_ic[i, j, Nz] for i in 1:Nx, j in 1:Ny]), "IC surface S (psu)", :haline, nothing),
          (S_sill, @sprintf("IC S at z ≈ %.0f m (psu) — sill depth", zc[k_sill]), :haline, nothing))

fig = Figure(size = (1500, 950))
for (n, (d, lab, cmap, lims)) in enumerate(panels)
    r, c = fld1(n, 2), mod1(n, 2)
    ax = Axis(fig[r, 2c-1]; title = lab, titlesize = 14, xlabel = "λ", ylabel = "φ", aspect = DataAspect())
    hm = lims === nothing ? heatmap!(ax, λc, φc, d; nan_color = :gray70, colormap = cmap) :
                            heatmap!(ax, λc, φc, d; nan_color = :gray70, colormap = cmap,
                                     colorrange = lims, highclip = :magenta)
    Colorbar(fig[r, 2c], hm; width = 12)
end
colgap!(fig.layout, 6)
outpng = joinpath(@__DIR__, "redsea_glorys_interior_check.png")
save(outpng, fig)
@info "Saved $outpng"
