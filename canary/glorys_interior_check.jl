# glorys_interior_check.jl -- Canary upwelling domain.
#
# Adapted from red_sea/glorys_interior_check.jl. The question is the same (does GLORYS actually cover
# everywhere this model has water?) but the answer matters more here, for two reasons:
#
#   1. The margin runs the WRONG WAY. GLORYS's deepest valid level over this region is -4836.6 m
#      against a raw ETOPO maximum of 4885 m, so the model floor would sit BELOW the parent unless
#      `depth` is capped. The Red Sea had +50 m of headroom; this has -48 m.
#   2. There are THREE open boundaries (west, south, north; the east edge is the African coast), all
#      in deep water, and the deepest point sits on the west/north corner -- i.e. the coverage
#      deficit bites exactly where boundary values are prescribed.
#
# Read-only. Honours CANARY_DELTA / CANARY_PASSES / CANARY_DEPTH.

using Oceananigans, NumericalEarth, CopernicusMarine, CairoMakie
using Dates, Printf, Statistics
using Oceananigans.Fields: location
using Oceananigans.ImmersedBoundaries: immersed_inactive_node
using Oceananigans.Units: Time
const SWP = Oceananigans.BuoyancyFormulations.SeawaterPolynomials
const DW  = NumericalEarth.DataWrangling

arch = CPU()
check_date = DateTime(get(ENV, "CANARY_START", "1993-01-01"))

λ₁, λ₂ = (-21.5, -10.0)
φ₁, φ₂ = ( 25.0,  29.5)
Δdeg = parse(Float64, get(ENV, "CANARY_DELTA", "0.02"))
Nx = round(Int, (λ₂-λ₁)/Δdeg); Ny = round(Int, (φ₂-φ₁)/Δdeg)
Nz = 32
depth   = parse(Float64, get(ENV, "CANARY_DEPTH", "4800"))
z_scale = parse(Float64, get(ENV, "CANARY_ZSCALE", "700"))
passes  = parse(Int, get(ENV, "CANARY_PASSES", "5"))
@info @sprintf("Canary grid %d x %d x %d (Δ=%.3f deg ≈ %.2f km), depth %.0f m, scale %.0f, %d passes",
               Nx, Ny, Nz, Δdeg, Δdeg*111.32*cosd(27.25), depth, z_scale, passes)

z_ref = ExponentialDiscretization(Nz, -depth, 0; scale = z_scale)
ug = LatitudeLongitudeGrid(arch; size=(Nx,Ny,Nz), latitude=(φ₁,φ₂), longitude=(λ₁,λ₂),
                           z = z_ref, halo=(7,7,7), topology=(Bounded,Bounded,Bounded))
bottom_height = regrid_bathymetry(ug; minimum_depth=10, interpolation_passes=passes)
grid = ImmersedBoundaryGrid(ug, GridFittedBottom(bottom_height))

hb = Array(interior(bottom_height))[:,:,1]
λc = λnodes(grid,Center(),Center(),Center()); φc = φnodes(grid,Center(),Center(),Center())
zc = znodes(grid,Center(),Center(),Center())
wet2d = [!immersed_inactive_node(i,j,Nz,grid,Center(),Center(),Center()) for i in 1:Nx, j in 1:Ny]
wet3d = [!immersed_inactive_node(i,j,k,grid,Center(),Center(),Center()) for i in 1:Nx, j in 1:Ny, k in 1:Nz]
wetcols = findall(wet2d)
@info @sprintf("wet columns %d/%d (%.1f%%) | wet cells %d | deepest column %.0f m | Δz top %.2f m",
               length(wetcols), Nx*Ny, 100length(wetcols)/(Nx*Ny), count(wet3d), -minimum(hb),
               reverse(diff(znodes(grid,Center(),Center(),Face())))[1])

# --- parent validity mask -----------------------------------------------------
glorys_region = DW.BoundingBox(longitude=(-22.0,-9.5), latitude=(24.5,30.0))
glorys = DW.GLORYSDaily()
# CANARY_GLORYS_DIR points at a gcsfuse mount on GCP; unset falls back to the local download cache.
const GLORYS_KW = let d = get(ENV, "CANARY_GLORYS_DIR", ""); isempty(d) ? (;) : (; dir = d) end
S_raw = Field(DW.Metadatum(:salinity; dataset=glorys, region=glorys_region, date=check_date, GLORYS_KW...); inpainting=nothing)
Sr = Array(interior(S_raw)); ngx,ngy,ngz = size(Sr)
gz  = znodes(S_raw.grid, map(L->L(), location(S_raw))...)
gzf = znodes(S_raw.grid, Center(), Center(), Face())
gλ  = λnodes(S_raw.grid, map(L->L(), location(S_raw))...)
gφ  = φnodes(S_raw.grid, map(L->L(), location(S_raw))...)

gbot = fill(NaN,ngx,ngy); gbotf = fill(NaN,ngx,ngy); gktop = zeros(Int,ngx,ngy)
for j in 1:ngy, i in 1:ngx
    k = findfirst(k->isfinite(Sr[i,j,k]), 1:ngz)
    k === nothing && continue
    gktop[i,j]=k; gbot[i,j]=gz[k]; gbotf[i,j]=gzf[k]
end
@info @sprintf("GLORYS %d x %d x %d | wet columns %d | deepest valid level over region %.1f m",
               ngx,ngy,ngz, count(>(0),gktop), minimum(gbot[isfinite.(gbot)]))

gΔλ=(last(gλ)-first(gλ))/(ngx-1); gΔφ=(last(gφ)-first(gφ))/(ngy-1)
gi(λ)=clamp(round(Int,(λ-first(gλ))/gΔλ)+1,1,ngx); gj(φ)=clamp(round(Int,(φ-first(gφ))/gΔφ)+1,1,ngy)

margin=fill(NaN,Nx,Ny); marginf=fill(NaN,Nx,Ny); dry=falses(Nx,Ny)
for I in wetcols
    i,j = Tuple(I); a,b = gi(λc[i]), gj(φc[j])
    if gktop[a,b]==0; dry[i,j]=true
    else margin[i,j]=hb[i,j]-gbot[a,b]; marginf[i,j]=hb[i,j]-gbotf[a,b] end
end
below = falses(Nx,Ny,Nz)
for k in 1:Nz, I in wetcols
    i,j=Tuple(I); wet3d[i,j,k] || continue
    a,b=gi(λc[i]),gj(φc[j]); below[i,j,k] = gktop[a,b]==0 || zc[k] < gbot[a,b]
end

println("\n","="^92); println("STAGE 1 -- COVERAGE"); println("="^92)
fin=[margin[I] for I in wetcols if isfinite(margin[I])]; finf=[marginf[I] for I in wetcols if isfinite(marginf[I])]
@printf("  wet columns %d | parent column dry %d | below deepest valid CENTRE %d (%.2f%%) | below BOTTOM FACE %d (%.2f%%)\n",
        length(wetcols), count(dry), count(<=(0),fin), 100count(<=(0),fin)/length(fin),
        count(<=(0),finf), 100count(<=(0),finf)/length(finf))
@printf("  centre-margin quantiles (m): min %.0f | 1%% %.0f | 5%% %.0f | 50%% %.0f | max %.0f\n",
        minimum(fin), quantile(fin,0.01), quantile(fin,0.05), median(fin), maximum(fin))
@printf("  wet cells below parent: %d of %d (%.2f%%)\n", count(below), count(wet3d), 100count(below)/count(wet3d))

# Per-boundary, because that is where prescribed values live.
println("\n  OPEN BOUNDARIES (the deficit matters most here):")
for (nm, idx) in (("west  λ=$(λ₁)", [(1,j)  for j in 1:Ny]),
                  ("south φ=$(φ₁)", [(i,1)  for i in 1:Nx]),
                  ("north φ=$(φ₂)", [(i,Ny) for i in 1:Nx]))
    w=[(i,j) for (i,j) in idx if wet2d[i,j]]; isempty(w) && (@printf("    %-14s ALL LAND\n",nm); continue)
    m=[margin[i,j] for (i,j) in w if isfinite(margin[i,j])]
    nb=count(below[i,j,k] for (i,j) in w, k in 1:Nz if wet3d[i,j,k])
    nw=count(wet3d[i,j,k] for (i,j) in w, k in 1:Nz)
    @printf("    %-14s %4d wet cols | below centre %4d (%5.1f%%) | min margin %+7.0f m | cells below %d/%d (%.2f%%)\n",
            nm, length(w), count(<=(0),m), 100count(<=(0),m)/length(m), minimum(m), nb, nw, 100nb/nw)
end

# --- Stage 2/3: what the pipeline puts there, and is it stable ----------------
bc_inpainting = DW.NearestNeighborInpainting(100)
function fill_down!(a,floor_value)
    nx,ny,nz=size(a)
    @inbounds for j in 1:ny, i in 1:nx
        last=NaN
        for k in nz:-1:1
            v=a[i,j,k]
            if isfinite(v) && v>floor_value; last=v elseif isfinite(last); a[i,j,k]=last end
        end
    end
    a
end
dates=[check_date, check_date+Day(1)]
mk(var)=DW.Metadata(var; dataset=glorys, dates, region=glorys_region, GLORYS_KW...)
S_fts=FieldTimeSeries(mk(:salinity),arch;time_indices_in_memory=2,inpainting=bc_inpainting)
T_fts=FieldTimeSeries(mk(:temperature),arch;time_indices_in_memory=2,inpainting=bc_inpainting)
for n in 1:2; fill_down!(interior(S_fts[n]),5.0); fill_down!(interior(T_fts[n]),1.0); end

S_ic=zeros(Nx,Ny,Nz); T_ic=zeros(Nx,Ny,Nz); t₀=Time(0.0)
for k in 1:Nz, j in 1:Ny, i in 1:Nx
    S_ic[i,j,k]=DW.interpolate_to_grid(S_fts,i,j,k,S_fts.grid,grid,t₀)
    T_ic[i,j,k]=DW.interpolate_to_grid(T_fts,i,j,k,T_fts.grid,grid,t₀)
end
widx=findall(wet3d)
println("\n","="^92); println("STAGE 2 -- PROVENANCE"); println("="^92)
@printf("  IC over wet cells: T ∈ [%.3f, %.3f] | S ∈ [%.3f, %.3f]\n",
        minimum(T_ic[I] for I in widx), maximum(T_ic[I] for I in widx),
        minimum(S_ic[I] for I in widx), maximum(S_ic[I] for I in widx))
bad=[I for I in widx if S_ic[I]<30 || T_ic[I]<0 || T_ic[I]>35]
@printf("  non-physical wet cells (S<30, T<0 or T>35): %d\n", length(bad))
println("  water masses by depth band (open ocean, λ < -17):")
for (lo,hi,lab) in ((0,-100,"surface / Ekman"),(-100,-300,"upwelling source"),(-300,-1000,"NACW/AAIW"),(-1000,-3000,"NADW"),(-3000,-5000,"abyssal"))
    idx=[I for I in widx if λc[Tuple(I)[1]] < -17 && hi <= zc[Tuple(I)[3]] < lo]
    isempty(idx) && continue
    s=[S_ic[I] for I in idx]; t=[T_ic[I] for I in idx]
    @printf("    %-18s n=%7d  S %.3f..%.3f (med %.3f)  T %.2f..%.2f (med %.2f)\n",
            lab, length(idx), minimum(s), maximum(s), median(s), minimum(t), maximum(t), median(t))
end

eos = SWP.TEOS10.TEOS10EquationOfState()
inv = Tuple{Float64,Int,Int,Int}[]
for j in 1:Ny, i in 1:Nx, k in 1:(Nz-1)
    (wet3d[i,j,k] && wet3d[i,j,k+1]) || continue
    zr=0.5*(zc[k]+zc[k+1])
    Δρ=SWP.ρ(T_ic[i,j,k+1],S_ic[i,j,k+1],zr,eos)-SWP.ρ(T_ic[i,j,k],S_ic[i,j,k],zr,eos)
    Δρ>1e-4 && push!(inv,(Δρ,i,j,k))
end
sort!(inv,by=first,rev=true)
npairs=count(I->(Tuple(I)[3]<Nz && wet3d[Tuple(I)[1],Tuple(I)[2],Tuple(I)[3]+1]), widx)
println("\n","="^92); println("STAGE 3 -- STATIC STABILITY"); println("="^92)
@printf("  unstable adjacent pairs: %d of %d (%.4f%%)\n", length(inv), npairs, 100length(inv)/npairs)
if !isempty(inv)
    tail=filter(x->x[1]>=0.1,inv)
    @printf("  tail Δρ ≥ 0.1 kg/m³: %d pairs (%.4f%%) | worst %.3f kg/m³\n", length(tail), 100length(tail)/npairs, inv[1][1])
    isext(x)=below[x[2],x[3],x[4]]||below[x[2],x[3],x[4]+1]
    @printf("  of the tail, %d (%.0f%%) involve a cell below the parent's deepest valid centre\n",
            count(isext,tail), isempty(tail) ? 0 : 100count(isext,tail)/length(tail))
    println("  worst 6:")
    for (Δρ,i,j,k) in inv[1:min(6,end)]
        @printf("    %8.3fE %7.3fN z=%8.1f/%8.1f m  Δρ=%+.4f | upper T=%.2f S=%.3f / lower T=%.2f S=%.3f %s\n",
                λc[i],φc[j],zc[k+1],zc[k],Δρ,T_ic[i,j,k+1],S_ic[i,j,k+1],T_ic[i,j,k],S_ic[i,j,k],
                (below[i,j,k]||below[i,j,k+1]) ? "[below parent]" : "")
    end
end

# --- figure -------------------------------------------------------------------
nan_land(a)=[wet2d[i,j] ? Float64(a[i,j]) : NaN for i in 1:Nx, j in 1:Ny]
worst=zeros(Nx,Ny); for (Δρ,i,j,k) in inv; worst[i,j]=max(worst[i,j],Δρ); end
k100=argmin(abs.(zc .+ 100))
fig=Figure(size=(1500,900))
panels=((nan_land(marginf),"model floor − parent bottom face (m)\nblue = parent has no water this deep",:balance,(-400,400)),
        (nan_land(worst),"worst IC density inversion (kg/m³)",:inferno,(0,0.3)),
        (nan_land([S_ic[i,j,Nz] for i in 1:Nx,j in 1:Ny]),"IC surface S (psu)",:haline,nothing),
        ([wet3d[i,j,k100] ? T_ic[i,j,k100] : NaN for i in 1:Nx,j in 1:Ny], @sprintf("IC T at z ≈ %.0f m (°C) — upwelling source",zc[k100]),:thermal,nothing))
for (n,(d,lab,cmap,lims)) in enumerate(panels)
    r,c=fld1(n,2),mod1(n,2)
    ax=Axis(fig[r,2c-1];title=lab,titlesize=13,xlabel="λ",ylabel="φ",aspect=DataAspect())
    hm = lims===nothing ? heatmap!(ax,λc,φc,d;nan_color=:gray70,colormap=cmap) :
                          heatmap!(ax,λc,φc,d;nan_color=:gray70,colormap=cmap,colorrange=lims,highclip=:magenta)
    Colorbar(fig[r,2c],hm;width=12)
end
colgap!(fig.layout,6)
out=joinpath(@__DIR__,"canary_glorys_interior_check.png"); save(out,fig)
@info "Saved $out"
