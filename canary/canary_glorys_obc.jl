# canary_glorys_obc.jl
#
# GLORYS-driven open boundaries for the Canary upwelling region (NW Africa). Three open boundaries
# (west, south, north; the east edge is the African coast), a wind-driven interior, and the Canary
# archipelago shedding wakes inside the domain.
#
# Values below are the 90-day reference configuration. See canary/README.md for the validation runs
# and docs/obc-implementation-recipe.md for the boundary-condition reasoning. A minimal, teaching
# version of the same boundary wiring is in examples/three_open_boundaries_canary.jl.
#
# GLORYS-driven open boundaries for the Canary upwelling region (NW Africa).
# Adapted from red_sea/redsea_glorys_obc.jl (monitors, atmosphere, diagnostics) and
# salish_sea/salish_glorys_obc.jl (the multi-boundary pattern). See docs/obc-implementation-recipe.md.
#
# What this domain tests that the others could not:
#   * THREE open boundaries, all in deep open ocean the parent resolves well -- the Red Sea had one,
#     the Salish Sea had two in badly-resolved channels.
#   * A WIND-DRIVEN interior. Upwelling is the signal, so this is the first domain where the
#     atmosphere is the point rather than a correction. It was rejected as a candidate earlier
#     precisely because the configuration had no wind forcing; JRA55 changed that.
#   * Islands in the flow: the Canary archipelago sits inside the domain and sheds wakes the parent
#     cannot resolve, which is the nesting question in a different guise from Bab el-Mandeb.

using CairoMakie
using Oceananigans
using Oceananigans.Architectures: on_architecture, architecture
using Oceananigans.TimeSteppers: update_state!
const JLD2 = Oceananigans.OutputReaders.JLD2
using NumericalEarth
using CopernicusMarine
using Oceananigans.Advection: cell_advection_timescale
using Oceananigans.Units
using CFTime
using Dates
using Printf
using Statistics

const DW  = NumericalEarth.DataWrangling
const SWP = Oceananigans.BuoyancyFormulations.SeawaterPolynomials

# --- Parameters ---------------------------------------------------------------
# The ones that cannot be substituted straight into the expression that uses them: either they are
# read before the line that needs them (USE_GPU gates `using CUDA`) or they appear in several places.
const USE_GPU = true
const ALLOW_DRY_RELEASE = false

# Open-boundary relaxation timescales, shared by whichever scheme the run selects. Inf on outflow
# lets the interior leave untouched; the finite inflow scale is what pulls returning water back
# towards the parent.
const OBC_INFLOW  = 6hours
const OBC_OUTFLOW = Inf

# :auto reads the parent time origin off the data (see BND_T0 below).
const PARENT_T0 = :auto

# Spin-up: a short fixed-Δt crawl before the wizard takes over, so the initial adjustment shock does
# not set the timestep for the whole run.
const SPINUP_DT    = 10
const SPINUP_ITERS = Int(100)

# Divisor for the pacing line's s/iter. Must match the interval `progress` is registered with: it was
# hardcoded to 10 while progress_iters was already a parameter, so any other value misreported the
# rate by their ratio.
const PROGRESS_ITERS = Int(10)

# Where outputs, checkpoints and the pickup search go. "." is the container's working directory under
# the ModelBenchmarking GCP backend, which mounts the job's own bucket there.
const output_dir = "."
isdir(output_dir) || mkpath(output_dir)

# use_gpu: ["true"] for the long tracer runs. Everything below is written so the same script runs on
# either: the boundary conditions and the tracer forcing are device kernels already, and every
# DIAGNOSTIC operates on CPU copies (see grid_c / arr below) because GPU disallows the scalar indexing
# they all use.
# GPU() is defined by Oceananigans' CUDA extension, which only loads once CUDA itself is loaded.
if USE_GPU
    using CUDA
    # Workaround: ZeroField() defaults to an Int placeholder, and NumericalEarth's
    # _assemble_net_ocean_fluxes! kernel hits CUDA_ERROR_ILLEGAL_ADDRESS when it receives those
    # (verified: identical call with ZeroField{Float64} runs). Adding 0.0 instead of 0 changes nothing.
    @eval Oceananigans.Fields ZeroField() = ZeroField{Float64, 3}()
end
arch = USE_GPU ? GPU() : CPU()
@info "Architecture: $(typeof(arch).name.name)"

# --- Run window ---------------------------------------------------------------
# January, matching salish_sea/ and red_sea/ so the three domains share a window. Note the cost:
# Canary upwelling is quasi-permanent but peaks in boreal SUMMER at these latitudes, so a January
# run exercises the machinery on a weaker signal than the domain can produce. start_date: ["1993-07-01"]
# for the strong-upwelling case.
run_start   = DateTime("1993-01-01")
run_days    = Int(90)
# Parent data is fetched for run_days + margin_days, so the boundary series outlasts the run
# even after the time-step wizard overshoots the last output.
margin_days = Int(2)
boundary_dates = run_start : Day(1) : (run_start + Day(run_days + margin_days))

# --- Grid ---------------------------------------------------------------------
# Measured from ETOPO, not assumed: the EAST edge at -10.0 is entirely land (the African coast), so
# WEST, SOUTH and NORTH are the open boundaries. West is 100% wet at 4578-4885 m; north 99% wet;
# south 58% wet, its land being the African coast at the SE corner. Combined open area ~7.6e3 km²
# against 719 km² for the Red Sea.
# const: read inside GPU kernels (sponge mask), where non-const globals cannot be compiled.
const λ₁, λ₂ = (Float64(-21.5), Float64(-10.0))
# lat_south extends the domain south for the domain-sizing experiments. At 22.0 the release sits
# 685 km from the south edge but only 600 km from the WEST edge, so west becomes the nearest open
# boundary and the dominant tracer leak -- see canary/domain-sizing.md.
const φ₁, φ₂ = (Float64(25.0), Float64(29.5))

# 0.02 deg = 2.06 km at 27.25N, 575 x 225 x 32 = 4.1M cells, i.e. about the Red Sea's cost per step.
# 0.01 deg (1.03 km) is the target for resolving island wakes and the upwelling filaments, but at
# 16.6M cells and half the timestep it is a GPU/HPC run, not a CPU one.
Δdeg = Float64(0.02)
Nx = round(Int, (λ₂ - λ₁) / Δdeg)
Ny = round(Int, (φ₂ - φ₁) / Δdeg)

# Vertical. GLORYS's deepest valid level over this region is -4836.6 m while raw ETOPO reaches
# 4885 m, so unlike the Red Sea the margin runs the WRONG WAY and `depth` is a hard cap rather than
# a convenience: 4800 keeps every model cell inside the parent's valid data, at the cost of clipping
# 1.2% of columns by at most 85 m of abyssal plain.
# scale=700 spends the 32 levels where upwelling lives: 1.2 m surface cell, 11 levels in the top
# 50 m, 16 in the top 150 m. The price is ~930 m bottom cells, which the quiescent abyss can carry.
Nz = Int(32)
depth   = Float64(4800)
z_scale = Float64(700)
z_reference = ExponentialDiscretization(Nz, -depth, 0; scale = z_scale)

# z* on in every run so far, as in red_sea/: with a 1.2 m surface cell the O(max|η|/Δz_top) error of a static
# (linear free surface) grid is not small. use_zstar: ["false"] falls back to static z.
use_zstar = true
z = use_zstar ? MutableVerticalDiscretization(z_reference) : z_reference

@info @sprintf("Grid: %d x %d x %d (Δ = %.3f deg ≈ %.2f km at 27.25N) | depth %.0f m | vertical: %s",
               Nx, Ny, Nz, Δdeg, Δdeg * 111.32 * cosd(27.25), depth,
               use_zstar ? "ZStar (mutable)" : "Z (static)")

underlying_grid = LatitudeLongitudeGrid(arch;
                             size = (Nx, Ny, Nz),
                             latitude = (φ₁, φ₂), longitude = (λ₁, λ₂),
                             z, halo = (7, 7, 7),
                             topology = (Bounded, Bounded, Bounded))

# The islands are the reason to keep smoothing light: at 0.02 deg El Hierro is only a few cells
# across, and over-smoothing would dissolve the wake-generating topography this domain exists to
# resolve. The startup check below reports each island's land-cell count.
interpolation_passes = Int(5)
bottom_height = regrid_bathymetry(underlying_grid;
                                  minimum_depth = 10,
                                  interpolation_passes = interpolation_passes)
grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(bottom_height))

# CPU mirror. The diagnostics all walk the grid cell by cell, and the metric operators (Vᶜᶜᶜ, Axᶠᶜᶜ,
# Azᶜᶜᶜ …) read the grid's own coordinate arrays, so BOTH the fields and the grid must be on CPU for
# them. On a CPU run this is the same grid and `arr` is a plain copy.
const grid_c = on_architecture(CPU(), grid)
@inline arr(f) = Array(interior(f))

# On CPU the mirror SHARES the z* state arrays (σ, η, ∂t_σ), so its metrics track the live grid for
# free. On GPU it is a genuine copy, and without this every diagnostic would silently use FROZEN
# metrics -- including salt_and_volume, which exists to measure z* drift and would then read ~0 no
# matter what z* did. Called once per progress block; a no-op on CPU.
const _ZSTAR_FIELDS = (:ηⁿ, :σᶜᶜⁿ, :σᶠᶜⁿ, :σᶜᶠⁿ, :σᶠᶠⁿ, :σᶜᶜ⁻, :∂t_σ)
function sync_grid_c!()
    (arch isa CPU || !(grid.underlying_grid.z isa MutableVerticalDiscretization)) && return nothing
    zg, zc = grid.underlying_grid.z, grid_c.underlying_grid.z
    for f in _ZSTAR_FIELDS
        hasfield(typeof(zg), f) || continue
        copyto!(parent(getfield(zc, f)), Array(parent(getfield(zg, f))))
    end
    return nothing
end

# Fixed subset so resolution changes do not invalidate the download cache.
# The parent box follows the domain with a 0.5 deg margin. Its extent is part of every cached
# filename, so changing φ₁ means a fresh download -- nothing already staged is reusable.
const φ_glorys_S = φ₁ - 0.5
glorys_region = DW.BoundingBox(longitude = (-22.0, -9.5), latitude = (φ_glorys_S, 30.0))
@assert λ₁ >= -22.0 && λ₂ <= -9.5 && φ₁ >= φ_glorys_S && φ₂ <= 30.0 "domain outside glorys_region"
glorys = DW.GLORYSDaily()

# Where the parent data lives. Unset (the default) falls back to NumericalEarth's own download cache,
# which is what a local run wants. On GCP these point at gcsfuse mounts -- GLORYS is in the obs bucket
# and JRA55 in the forcing bucket, so they are separate roots rather than one.
#   glorys_dir: ["/mnt/obs/GLORYS"]  jra55_dir: ["/mnt/forcing/JRA55"]
glorys_dir = "/mnt/obs/glorys"
jra55_dir  = ""
const GLORYS_KW = isempty(glorys_dir) ? (;) : (; dir = glorys_dir)
const JRA55_KW  = isempty(jra55_dir)  ? (;) : (; dir = jra55_dir)
isempty(glorys_dir) || @info "GLORYS dir: $glorys_dir"
isempty(jra55_dir)  || @info "JRA55 dir:  $jra55_dir"

# --- Passive tracer release (optional) ----------------------------------------
# For pathway and residence-time work: a continuous release at a point on Gran Canaria, of the kind
# an ocean alkalinity enhancement outfall would represent.
#
# TWO tracers, following Deleersnijder et al.'s age formulation:
#   c  -- concentration.      ∂c/∂t  + ∇·(uc)  = S      with S the release
#   ac -- age concentration.  ∂ac/∂t + ∇·(uac) = c      (grows at rate c wherever tracer is present)
# so the MEAN AGE of the water carrying the tracer is ac/c. A single "ideal age" tracer cannot do
# this for a point source -- it measures time since surface contact, not time since release.
#
# S is normalised to ∫S dV = 1 unit/s over the release volume, which makes the budget checkable by
# inspection: total inventory should equal elapsed seconds minus whatever has left through the open
# boundaries. Default site is Las Palmas (-15.38, 28.16): wet, 535 m, 47 of 49 neighbouring cells wet,
# and the island's main port.
use_tracer  = true
# Release site, south coast of Gran Canaria. The two completed 90-day runs (phi_S = 25 and 22) used
# the earlier point at (-15.38, 28.16) off the NORTH tip, so they are NOT directly comparable to runs
# made after this change -- it moves the source ~39 km south, i.e. closer to the southern boundary
# and further from the northern one.
# Coastal release off the south-east of Gran Canaria, 2.5 km from the requested (27.812758,
# -15.422082). That exact point regrids to LAND at 0.02 deg (bottom height +11.1 m), which does not
# fail loudly: the Gaussian would sit on a dry cell and the release would leak out through whatever
# wet cells its tail reached, displacing the effective source. This point is 60 m deep with 12
# vertical levels beneath it and 81% of the 3σ patch wet -- deliberately shallow, because the release
# is meant to be coastal. At 0.01 deg the coastline resolves better and it could move back inshore.
const λ_release   = Float64(-15.397082)
const φ_release   = Float64(27.807758)
# Sharp 2-D Gaussian at the surface: σ = 0.02° is one grid cell at the 0.02° target (1.98 km
# zonal, 2.22 km meridional), so the source spans a few cells and is cut off at 3σ. The depth
# covers the top cell only (1.21 m at z_scale: [700]), i.e. a surface release.
const σ_release   = Float64(0.02)     # degrees
const R_release   = 3σ_release                                                    # hard cutoff
const D_release   = Float64(1.5)       # metres below surface

# --- Boundary conditions ------------------------------------------------------
using Oceananigans.BoundaryConditions: NormalRadiation, PerturbationAdvection,
                                       GravityWaveRadiationBoundaryCondition,
                                       SurfaceWaveRadiationBoundaryCondition
using Oceananigans.Operators: Δzᶜᶜᶜ, Δzᶠᶜᶜ, Δzᶜᶠᶜ, Axᶠᶜᶜ, Ayᶜᶠᶜ, Azᶜᶜᶜ, Vᶜᶜᶜ
using Oceananigans.ImmersedBoundaries: immersed_inactive_node, immersed_peripheral_node
using Oceananigans.Units: Time

# radiation: ["oblique"] selects Raymond & Kuo (1984) two-dimensional radiation, which diagnoses a
# tangential phase speed as well as a normal one. It exists only from Oceananigans 0.113, so the
# lookup is conditional and one code path still runs on 0.110.17.
# Default is normal: in an idealised truncation experiment the oblique scheme was within noise of
# NormalRadiation, while adding the tangential boundary condition (now unconditional, below) was
# worth +32-39% of retained enstrophy. See obc_v113/README.md.
const HAS_OBLIQUE = isdefined(Oceananigans.BoundaryConditions, :ObliqueRadiation)
radiation_kind = "oblique"
radiation_kind in ("normal", "oblique") || error("radiation must be \"normal\" or \"oblique\"")
if radiation_kind == "oblique" && !HAS_OBLIQUE
    error("radiation: oblique needs Oceananigans >= 0.113; this environment has no ObliqueRadiation")
end
const RAD_NORMAL = NormalRadiation(inflow_timescale = OBC_INFLOW, outflow_timescale = OBC_OUTFLOW)
const RAD_OBLIQUE = HAS_OBLIQUE ?
    getfield(Oceananigans.BoundaryConditions, :ObliqueRadiation)(inflow_timescale = OBC_INFLOW, outflow_timescale = OBC_OUTFLOW) :
    nothing

# oblique_sides restricts the oblique scheme to named boundaries (comma list of west/south/
# north); unset means all of them. This exists because oblique diverged at 0.02 deg with salinity
# reaching 45.7 psu, and the blow-up sat on cell (332, 1) -- the LAST WET CELL of the south boundary,
# with land in the next cell along. Oblique's tangential derivative runs ALONG the boundary, so there
# it is evaluated across a land/water interface. West is 100% wet and north 99%; south is 58% wet and
# is the only boundary that terminates in coast. See docs/obc-implementation-recipe.md.
oblique_sides = Symbol.(split("", ',', keepempty = false))
use_oblique(side) = radiation_kind == "oblique" && (isempty(oblique_sides) || side in oblique_sides)
rad(side) = use_oblique(side) ? RAD_OBLIQUE : RAD_NORMAL
radiation = rad(:west)          # tracer BCs and anything not side-aware use this
@info "Open-boundary radiation: $(radiation_kind)" *
      (radiation_kind == "oblique" && !isempty(oblique_sides) ? " on $(join(oblique_sides, ",")) only" : "")

# normal_scheme: ["perturbation"] uses PerturbationAdvection for the boundary-NORMAL velocity
# instead of the radiation scheme, keeping radiation for the tangential velocity and the tracers.
# That is the split the Mid-Atlantic Bight configuration in NumericalEarthSandbox uses; ours has
# always used one scheme throughout. Available on both 0.110.17 and 0.113. See the recipe doc.
normal_kind = "radiation"
normal_kind in ("radiation", "perturbation") || error("normal_scheme must be \"radiation\" or \"perturbation\"")
normal_scheme = normal_kind == "perturbation" ?
    PerturbationAdvection(inflow_timescale = OBC_INFLOW, outflow_timescale = OBC_OUTFLOW) : radiation
normal_kind == "perturbation" && @info "Boundary-normal velocity scheme: PerturbationAdvection"

@inline wetcell(i, j, k, grid, ℓx, ℓy, ℓz) =
    !immersed_peripheral_node(i, j, k, grid, ℓx, ℓy, ℓz) & !immersed_inactive_node(i, j, k, grid, ℓx, ℓy, ℓz)

# --- Nudging band -------------------------------------------------------------
# Bands on all three open edges, combined with max() so the corners do not double-count. Kept narrow
# (0.25 deg ≈ 12 cells at 0.02) so the shelf and the island wakes run free: the Salish lesson was
# that a wide sponge pins the very features the domain exists to resolve, and here the upwelling
# front sits against the eastern coast, far from every sponge.
const sponge_width  = Float64(0.25)
restoring_timescale = 12hours
restoring_rate      = 1 / restoring_timescale

@inline gaussian_edge(d, w) = exp(-(d / w)^2)
@inline function boundary_mask(x, y, z, t)
    dw = x - λ₁        # west
    ds = y - φ₁        # south
    dn = φ₂ - y        # north
    return max(gaussian_edge(dw, sponge_width),
               gaussian_edge(ds, sponge_width),
               gaussian_edge(dn, sponge_width))
end

restoring_meta(var) = DW.Metadata(var; dataset = glorys, dates = boundary_dates, region = glorys_region, GLORYS_KW...)
bc_inpainting = DW.NearestNeighborInpainting(100)

# Vertical constant extrapolation below the parent's deepest valid level. Insurance in the Red Sea;
# here it is expected to do real work, because `depth` sits only ~37 m above the parent's deepest
# valid centre and the abyssal columns are sparsely covered (2.9% of columns have data at -4836 m).
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
# via a host copy: fill_down! walks the array cell by cell, which a device array will not allow.
function fill_down_fts!(fts, floor)
    for n in 1:length(fts.times)
        a = Array(interior(fts[n])); fill_down!(a, floor); copyto!(interior(fts[n]), a)
    end
    return fts
end

n_times = length(boundary_dates)
restoring_kw = (; rate = restoring_rate, mask = boundary_mask, inpainting = bc_inpainting,
                  time_indices_in_memory = n_times)

T_restoring = DW.DatasetRestoring(restoring_meta(:temperature), arch; restoring_kw...)
S_restoring = DW.DatasetRestoring(restoring_meta(:salinity),    arch; restoring_kw...)
u_restoring = DW.DatasetRestoring(restoring_meta(:u_velocity),  arch; restoring_kw...)
v_restoring = DW.DatasetRestoring(restoring_meta(:v_velocity),  arch; restoring_kw...)

# Canary Current water is S 35.5-36.8, T 15-24 at the surface; floors sit well below anything real.
fill_down_fts!(T_restoring.field_time_series, 1.0)
fill_down_fts!(S_restoring.field_time_series, 5.0)
let Ti = interior(T_restoring.field_time_series[1]), Si = interior(S_restoring.field_time_series[1])
    @info @sprintf("GLORYS native after fill_down: T ∈ [%.3f, %.3f], S ∈ [%.3f, %.3f] | residual zeros %.2f%%",
                   minimum(Ti), maximum(Ti), minimum(Si), maximum(Si), 100count(==(0), Si) / length(Si))
end

const T_fts, T_native = T_restoring.field_time_series, T_restoring.native_grid
const S_fts, S_native = S_restoring.field_time_series, S_restoring.native_grid
const u_fts, u_native = u_restoring.field_time_series, u_restoring.native_grid
const v_fts, v_native = v_restoring.field_time_series, v_restoring.native_grid

# --- GLORYS section residual (diagnostic only) --------------------------------
# No net-transport correction: measured in red_sea/ to change nothing a Flather boundary does not
# already absorb (see docs/obc-implementation-recipe.md). Q_raw is reported per boundary instead,
# because with three sections it also says which edge the parent disagrees with most.

# --- Exteriors ----------------------------------------------------------------
# Inflow is +u at west, +v at south, -v at north. Oceananigans applies the side-dependent Flather
# sign itself (west/south: U = Uᵉˣᵗ − √(gH)(ηᵇ − ηᵉˣᵗ); north/east: +), so these closures return the
# exterior value only. Do NOT pre-flip signs -- that double inversion is the bug red_sea/ carried.
#
# The parent is interpolated onto the three boundary faces ONCE, here, into plain (face × depth × time)
# arrays that the boundary kernels read by index and interpolate linearly in time. Sampling the
# NumericalEarth series inside the kernels instead works on CPU but not on GPU: their deeply
# parametric types sent the CUDA compiler's inliner past 31 GB before the first halo fill.

# Host copies of the parent series. On CPU these are the series themselves; on GPU the device series
# cannot be sampled point by point, so CPU copies are loaded with the same repair.
function host_series(dev, var, floor)
    arch isa CPU && return dev
    f = FieldTimeSeries(restoring_meta(var), CPU(); time_indices_in_memory = n_times, inpainting = bc_inpainting)
    floor === nothing || fill_down_fts!(f, floor)
    return f
end
const T_fts_c = host_series(T_fts, :temperature, 1.0)
const S_fts_c = host_series(S_fts, :salinity,    5.0)
const u_fts_c = host_series(u_fts, :u_velocity,  nothing)
const v_fts_c = host_series(v_fts, :v_velocity,  nothing)

zos_meta = DW.Metadata(:free_surface; dataset = glorys, dates = boundary_dates, region = glorys_region, GLORYS_KW...)
const zos_fts_c = FieldTimeSeries(zos_meta, CPU();            # only ever sampled on the host now
                                  time_indices_in_memory = n_times,
                                  time_indexing = Oceananigans.OutputReaders.Cyclical(),
                                  inpainting = bc_inpainting)

const BND_TIMES = Float64.(collect(u_fts_c.times))
const BND_Δt    = BND_TIMES[2] - BND_TIMES[1]
# The kernels need UNIFORM parent spacing, but not a zero origin. NumericalEarth 0.8 stamps
# GLORYSDaily at the CENTRE of its daily-mean window (12:00) rather than at 00:00, so the series
# starts at 43200 s there and at 0 s on 0.6.1. Carrying the origin keeps one code path correct on
# both -- see obc_v113/README.md.
@assert all(diff(BND_TIMES) .≈ BND_Δt) "boundary kernels assume uniformly spaced parent times"
# parent_t0: "auto" reads the origin off the data, which is always what a science run wants. A
# number overrides it, and exists only for migration checking: forcing 0 on NumericalEarth 0.8
# reproduces 0.6.1's indexing (same daily files, old assumed timestamps), isolating the
# time-convention change from every other difference between the versions.
const BND_T0    = PARENT_T0 === :auto ? BND_TIMES[1] : Float64(PARENT_T0)
const BND_Nt    = length(BND_TIMES)

function face_array(fts, which, loc_i, loc_j)
    FT = eltype(interior(fts[1]))
    na = which === :west ? Ny : Nx
    A = zeros(FT, na, Nz, BND_Nt)
    for n in 1:BND_Nt, k in 1:Nz, a in 1:na
        i, j = which === :west ? (1, a) : which === :south ? (a, 1) : (a, loc_j)
        A[a, k, n] = DW.interpolate_to_grid(fts, i, j, k, fts.grid, grid_c, Time(BND_TIMES[n]))
    end
    return on_architecture(arch, A)
end

function face_η(which)
    na = which === :west ? Ny : Nx
    A = zeros(Float64, na, BND_Nt)
    for n in 1:BND_Nt, a in 1:na
        i, j = which === :west ? (1, a) : which === :south ? (a, 1) : (a, Ny)
        A[a, n] = DW.interpolate_to_grid(zos_fts_c, i, j, 1, zos_fts_c.grid, grid_c, Time(BND_TIMES[n]))
    end
    return on_architecture(arch, A)
end

# Velocity at the boundary FACE index (west i=1, south j=1, north j=Ny+1); tracers and η at the
# adjacent interior node, matching salish_sea/.
const pW = (; u = face_array(u_fts_c, :west,  1, 1),    T = face_array(T_fts_c, :west,  1, 1),
              S = face_array(S_fts_c, :west,  1, 1),    η = face_η(:west),
              vt = face_array(v_fts_c, :west,  1, 1),
              Δt = BND_Δt, Nt = BND_Nt, T0 = BND_T0)
const pS = (; v = face_array(v_fts_c, :south, 1, 1),    T = face_array(T_fts_c, :south, 1, 1),
              S = face_array(S_fts_c, :south, 1, 1),    η = face_η(:south),
              ut = face_array(u_fts_c, :south, 1, 1),
              Δt = BND_Δt, Nt = BND_Nt, T0 = BND_T0)
const pN = (; v = face_array(v_fts_c, :north, 1, Ny+1), T = face_array(T_fts_c, :north, 1, Ny),
              S = face_array(S_fts_c, :north, 1, Ny),   η = face_η(:north),
              ut = face_array(u_fts_c, :north, 1, Ny),
              Δt = BND_Δt, Nt = BND_Nt, T0 = BND_T0)
@info @sprintf("Boundary exteriors precomputed: %d parent times at Δt = %.0f s", BND_Nt, BND_Δt)

@inline function bnd_time(t, p)
    # Relative to the FIRST parent time, not to zero: see BND_T0. Before it (the first half-day when
    # the parent is centred at noon) τ < 0, both clamps floor to the first sample, i.e. persistence.
    τ = (t - p.T0) / p.Δt
    n = clamp(unsafe_trunc(Int, τ) + 1, 1, p.Nt - 1)
    return n, clamp(τ - (n - 1), zero(τ), one(τ))
end
@inline function face3(A, a, k, t, p)
    n, w = bnd_time(t, p)
    @inbounds return (1 - w) * A[a, k, n] + w * A[a, k, n+1]
end
@inline function face2(A, a, t, p)
    n, w = bnd_time(t, p)
    @inbounds return (1 - w) * A[a, n] + w * A[a, n+1]
end

@inline u_west_exterior(j, k, grid, clock, fields, p)  = face3(p.u, j, k, clock.time, p)
@inline v_south_exterior(i, k, grid, clock, fields, p) = face3(p.v, i, k, clock.time, p)
@inline v_north_exterior(i, k, grid, clock, fields, p) = face3(p.v, i, k, clock.time, p)
# Tangential velocity exteriors: v at the west face, u at the south and north faces. Sampled at the
# adjacent interior node like T and S, since the tangential component is half a cell off the boundary.
@inline v_west_tangential(j, k, grid, clock, fields, p)  = face3(p.vt, j, k, clock.time, p)
@inline u_south_tangential(i, k, grid, clock, fields, p) = face3(p.ut, i, k, clock.time, p)
@inline u_north_tangential(i, k, grid, clock, fields, p) = face3(p.ut, i, k, clock.time, p)
@inline tracer_T_exterior(a, k, grid, clock, fields, p) = face3(p.T, a, k, clock.time, p)
@inline tracer_S_exterior(a, k, grid, clock, fields, p) = face3(p.S, a, k, clock.time, p)

# Depth-integrated exterior transport for the Flather condition, paired with the parent η.
@inline function west_U_obc(j, k, grid, clock, fields, p)
    U = zero(grid)
    @inbounds for kk in 1:grid.Nz
        wet = wetcell(1, j, kk, grid, Face(), Center(), Center())
        U += ifelse(wet, face3(p.u, j, kk, clock.time, p) * Δzᶠᶜᶜ(1, j, kk, grid), zero(grid))
    end
    return (U, face2(p.η, j, clock.time, p))
end
@inline function south_V_obc(i, k, grid, clock, fields, p)
    V = zero(grid)
    @inbounds for kk in 1:grid.Nz
        wet = wetcell(i, 1, kk, grid, Center(), Face(), Center())
        V += ifelse(wet, face3(p.v, i, kk, clock.time, p) * Δzᶜᶠᶜ(i, 1, kk, grid), zero(grid))
    end
    return (V, face2(p.η, i, clock.time, p))
end
@inline function north_V_obc(i, k, grid, clock, fields, p)
    j = grid.Ny + 1; V = zero(grid)
    @inbounds for kk in 1:grid.Nz
        wet = wetcell(i, j, kk, grid, Center(), Face(), Center())
        V += ifelse(wet, face3(p.v, i, kk, clock.time, p) * Δzᶜᶠᶜ(i, j, kk, grid), zero(grid))
    end
    return (V, face2(p.η, i, clock.time, p))
end

const res_in  = Float64(5000)
const res_out = Float64(0)
const use_reservoir = res_in > 0 || res_out > 0
if use_reservoir && !isdefined(Oceananigans.BoundaryConditions, :TracerReservoir)
    error("reservoir_in_m/_out_m need Oceananigans >= 0.113.2; this environment has no TracerReservoir")
end

# T and S carry their own scheme, separate from the velocities. The Mid-Atlantic Bight configuration
# runs oblique on VELOCITIES and TracerReservoir on tracers; feeding oblique to T/S is what diverged
# here, with salinity reaching 47 psu at the south boundary's coastal end.
#   same      - whatever the velocities use (previous behaviour)
#   normal    - NormalRadiation
#   reservoir - TracerReservoir at the reservoir_in_m/_out_m scales, as MAB does
ts_kind = "reservoir"
ts_kind in ("same", "normal", "reservoir") || error("ts_scheme must be same/normal/reservoir")
if ts_kind == "reservoir" && !isdefined(Oceananigans.BoundaryConditions, :TracerReservoir)
    error("ts_scheme: reservoir needs Oceananigans >= 0.113.2")
end
ts_scheme(side) = ts_kind == "normal"    ? RAD_NORMAL :
                  ts_kind == "reservoir" ? getfield(Oceananigans.BoundaryConditions, :TracerReservoir)(
                                               inflow_length_scale = res_in, outflow_length_scale = res_out) :
                  rad(side)
ts_kind == "same" || @info "T/S boundary scheme: $(ts_kind)"

nf(f, p, side) = NormalFlowBoundaryCondition(f; scheme = (normal_kind == "perturbation" ? normal_scheme : rad(side)),
                                               discrete_form = true, parameters = p)
vb(f, p, side) = ValueBoundaryCondition(f; scheme = rad(side), discrete_form = true, parameters = p)
vbts(f, p, side) = ValueBoundaryCondition(f; scheme = ts_scheme(side), discrete_form = true, parameters = p)
gw(f, p) = GravityWaveRadiationBoundaryCondition(f; discrete_form = true, parameters = p)

# The TANGENTIAL velocity carries an open boundary condition as well as the normal one: v at west,
# u at south and north. Unconditional, on three independent lines of evidence -- Oceananigans PR
# #5962 attributes its headline gain to radiating the tangential component, the truncation experiment
# in obc_v113/README.md measured it at +32-39% of retained enstrophy, and the Mid-Atlantic Bight
# configuration sets it by default (docs/obc-implementation-recipe.md).
# Everything run before this had NO condition here, so the tangential component fell back to the
# default wall behaviour; those runs are not comparable to these.
u_bcs = FieldBoundaryConditions(west  = nf(u_west_exterior,    pW, :west),
                                south = vb(u_south_tangential, pS, :south),
                                north = vb(u_north_tangential, pN, :north))
v_bcs = FieldBoundaryConditions(south = nf(v_south_exterior,   pS, :south),
                                north = nf(v_north_exterior,   pN, :north),
                                west  = vb(v_west_tangential,  pW, :west))
U_bcs = FieldBoundaryConditions(grid, (Face(), Center(), nothing);  west  = gw(west_U_obc, pW))
V_bcs = FieldBoundaryConditions(grid, (Center(), Face(), nothing);  south = gw(south_V_obc, pS),
                                                                    north = gw(north_V_obc, pN))
# Chapman: radiates interior η outward and never reads an exterior value. Present so the barotropic
# mode is well posed, not to impose anything -- see the recipe.
η_bcs = FieldBoundaryConditions(grid, (Center(), Center(), Face());
                                west  = SurfaceWaveRadiationBoundaryCondition(),
                                south = SurfaceWaveRadiationBoundaryCondition(),
                                north = SurfaceWaveRadiationBoundaryCondition())
T_bcs = FieldBoundaryConditions(west  = vbts(tracer_T_exterior, pW, :west), south = vbts(tracer_T_exterior, pS, :south),
                                north = vbts(tracer_T_exterior, pN, :north))
S_bcs = FieldBoundaryConditions(west  = vbts(tracer_S_exterior, pW, :west), south = vbts(tracer_S_exterior, pS, :south),
                                north = vbts(tracer_S_exterior, pN, :north))

# --- Tracer release forcing and boundaries ------------------------------------
# Release volume is measured on the grid so the source can be normalised; a continuous-form forcing
# (a pure function of x, y, z, t) is used rather than a masked array so the same code runs on GPU.
# In a function, not a top-level loop: accumulating into a global from inside `for` hits Julia's
# soft-scope rule and fails with UndefVarError.
@inline function release_weight(x, y, z)
    dx = (x - λ_release) * cosd(φ_release)
    dy = y - φ_release
    r2 = dx*dx + dy*dy
    inside = (r2 < R_release^2) & (z > -D_release)
    return ifelse(inside, exp(-r2 / (2 * σ_release^2)), zero(x))
end

# ∫ w dV over the release patch, so the source can be normalised to exactly 1 unit/s.
function compute_release_volume()
    λc = λnodes(grid_c, Center(), Center(), Center())
    φc = φnodes(grid_c, Center(), Center(), Center())
    zc = znodes(grid_c, Center(), Center(), Center())
    V = 0.0; ncell = 0
    for k in 1:Nz, j in 1:Ny, i in 1:Nx
        wetcell(i, j, k, grid_c, Center(), Center(), Center()) || continue
        w = release_weight(λc[i], φc[j], zc[k])
        w > 0 || continue
        V += w * Vᶜᶜᶜ(i, j, k, grid_c); ncell += 1
    end
    @info @sprintf("Release patch: %d wet cells within 3σ (σ = %.3f° ≈ %.1f km), top %.1f m",
                   ncell, σ_release, σ_release * 111.32 * cosd(φ_release), D_release)
    return V
end

# A release point on LAND does not fail loudly on its own, which is the reason this check exists: the
# Gaussian's tail still reaches wet cells, so tracer enters the domain from displaced cells with a
# truncated weighting and nothing in the output says so. Measured: (27.812758, -15.422082) regrids to
# land at 0.02 deg (bottom height +11.1 m) and still produced a 17-cell release that looked healthy.
function check_release_site()
    λc = λnodes(grid_c, Center(), Center(), Center())
    φc = φnodes(grid_c, Center(), Center(), Center())
    zf = znodes(grid_c, Center(), Center(), Face())
    i0 = argmin(abs.(λc .- λ_release)); j0 = argmin(abs.(φc .- φ_release))
    wet(i, j, k) = wetcell(i, j, k, grid_c, Center(), Center(), Center())
    centre_wet = wet(i0, j0, Nz)                       # k = Nz is the surface

    nwet = ntot = 0                                    # surface footprint of the 3σ patch
    for j in 1:Ny, i in 1:Nx
        dx = (λc[i] - λ_release) * cosd(φ_release); dy = φc[j] - φ_release
        dx^2 + dy^2 < R_release^2 || continue
        ntot += 1; wet(i, j, Nz) && (nwet += 1)
    end
    frac = ntot == 0 ? 0.0 : nwet / ntot

    kbot = findfirst(k -> wet(i0, j0, k), 1:Nz)        # lowest wet level in the release column
    dep  = (centre_wet && kbot !== nothing) ? -zf[kbot] : 0.0
    nlev = (centre_wet && kbot !== nothing) ? Nz - kbot + 1 : 0
    @info @sprintf("Release site %.4fE %.4fN: centre cell %s | depth %.0f m (%d levels) | patch %d/%d wet (%.0f%%)",
                   λ_release, φ_release, centre_wet ? "WET" : "DRY", dep, nlev, nwet, ntot, 100frac)

    if !centre_wet && !ALLOW_DRY_RELEASE
        error("Release site ($(λ_release), $(φ_release)) is LAND at this resolution. The run would still " *
              "emit tracer, from displaced cells with a truncated Gaussian, and look healthy. Move the " *
              "site, or set allow_dry_release: [true] in the config to proceed anyway.")
    end
    frac < 0.5 && @warn @sprintf("Only %.0f%% of the release patch is wet: the Gaussian is truncated by the coast, so the effective source sits off-centre", 100frac)
    return nothing
end

use_tracer && check_release_site()
release_volume = use_tracer ? compute_release_volume() : 0.0
if use_tracer
    release_volume <= 0 ?
        @warn("Release site ($(λ_release), $(φ_release)) contains NO wet cells -- tracer will never enter the domain") :
        @info @sprintf("Tracer release at %.3fE %.3fN: ∫w dV = %.4e m³", λ_release, φ_release, release_volume)
end

const RELEASE_RATE = release_volume > 0 ? 1 / release_volume : 0.0   # ∫S dV = 1 unit/s

@inline c_source(x, y, z, t) = RELEASE_RATE * release_weight(x, y, z)

# Zero exterior at every open boundary: there is no released tracer outside the domain, so this lets
# it leave freely. A default no-flux condition would trap it and make residence time meaningless.
@inline zero_west(j, k, grid, clock, fields)  = zero(grid)
@inline zero_ns(i, k, grid, clock, fields)    = zero(grid)
# reservoir_in_m / _out_m give the release tracers a MOM6-style boundary reservoir (metres),
# available from Oceananigans 0.113.2. A memoryless Value boundary hands re-entering water the
# exterior value regardless of what just left -- and the exterior here is ZERO -- so tracer that
# leaves and returns is destroyed. That is the same mechanism this domain's own data implicated: the
# 22N run retained MORE north of every shared latitude than the 25N run's virtual boundaries
# predicted, because the real boundary at 25.0 flushed tracer that would otherwise have come back
# (canary/domain-sizing.md). Upstream measures ~7-8% lost per oscillation at zero net
# displacement. Both scales 0, the default, reproduces the memoryless condition exactly.
# The release tracers carry their own reservoir scales, set independently of the T/S ones. Splitting
# them is what lets a reservoir on/off experiment vary ONE thing: T and S stay on the reservoir in
# both arms (oblique needs that to be stable) while only c and ac change.
# The fallback is NormalRadiation, never `radiation` -- oblique on a tracer at a land-terminated
# boundary is what drove salinity to 47 psu, and c/ac sit on the same boundaries as T/S.
const rel_in  = Float64(5000)
const rel_out = Float64(0)
const use_release_reservoir = rel_in > 0 || rel_out > 0
release_scheme = use_release_reservoir ?
    getfield(Oceananigans.BoundaryConditions, :TracerReservoir)(inflow_length_scale  = rel_in,
                                                                outflow_length_scale = rel_out) : RAD_NORMAL
use_release_reservoir ?
    @info(@sprintf("Release-tracer reservoir: inflow L = %.0f m, outflow L = %.0f m", rel_in, rel_out)) :
    @info("Release tracers: NormalRadiation (no reservoir)")

tracer_bcs() = FieldBoundaryConditions(
    west  = ValueBoundaryCondition(zero_west; scheme = release_scheme, discrete_form = true),
    south = ValueBoundaryCondition(zero_ns;   scheme = release_scheme, discrete_form = true),
    north = ValueBoundaryCondition(zero_ns;   scheme = release_scheme, discrete_form = true))

tracer_names   = use_tracer ? (:T, :S, :c, :ac) : (:T, :S)
# ac gets its source (∂ac/∂t = c) from a callback, not a Forcing: a forcing with field_dependencies
# makes the tracer kernel read model_fields, which fails to compile on GPU (InvalidIRError in
# free_surface_fields). Applied after each step, so the age source lags advection by one Δt.
tracer_forcing = use_tracer ? (; c = Forcing(c_source)) : (;)
tracer_bc_pair = use_tracer ? (; c = tracer_bcs(), ac = tracer_bcs()) : (;)

# --- Ocean simulation ---------------------------------------------------------
weno_order = Int(7)
weno_kw = (; time_discretization = AdaptiveVerticallyImplicitDiscretization(cfl = 0.5))
weno_scheme = WENO(order = weno_order; weno_kw...)
# tracer_bounds gives the release tracers bound-preserving WENO, to fix the c < 0 undershoot
# in the dilute plume edge (6.8% of surface cells at day 90) that makes the pointwise age ac/c
# meaningless there. Kept OFF: in Oceananigans 0.110.17 that path is NON-CONSERVATIVE under z*,
# creating tracer at +25%/day against +0.07% for plain WENO -- see canary/domain-sizing.md.
# T and S must stay unbounded: 0 is not a physical floor for either.
# NOTE: a NamedTuple here defaults any UNNAMED tracer to CenteredSecondOrder, so every tracer in
# `tracer_names` must appear below.
# Only the LOWER bound is physical; the upper is Inf because ac is a concentration-time product that
# any finite cap would clip. bounds must be a homogeneous NTuple{2}: (0, Inf) is Tuple{Int,Float64}
# and is rejected.
positive_weno = WENO(order = weno_order; bounds = (0.0, Inf), weno_kw...)
use_bounds = false
release_advection = use_bounds ? positive_weno : weno_scheme
tracer_advection = use_tracer ?
    (; T = weno_scheme, S = weno_scheme, c = release_advection, ac = release_advection) : weno_scheme

ocean = ocean_simulation(grid;
                         tracers = tracer_names,
                         tracer_advection,
                         forcing = merge((T = T_restoring, S = S_restoring,
                                          u = u_restoring, v = v_restoring), tracer_forcing),
                         boundary_conditions = merge((; u = u_bcs, U = U_bcs, v = v_bcs, V = V_bcs,
                                                        η = η_bcs, T = T_bcs, S = S_bcs), tracer_bc_pair))

# --- Atmospheric forcing ------------------------------------------------------
# The point of this domain, not a correction to it: upwelling is wind-driven, so atmosphere: none
# gives a boundary-only control run rather than a usable configuration.
# API notes and the ocean_simulation BC-merge trap are documented in red_sea/README.md.
atmos_source = "JRA"
atmos_source in ("JRA", "none") || error("atmosphere must be \"JRA\" or \"none\"")
use_atmos = atmos_source == "JRA"
jra_dataset_kind = "repeat"
jra_dataset_kind in ("repeat", "multiyear") || error("jra_dataset must be \"repeat\" or \"multiyear\"")
jra_multiyear = jra_dataset_kind == "multiyear"

atmosphere = nothing
radiation_component = nothing
if use_atmos
    jra_dataset = jra_multiyear ? NumericalEarth.MultiYearJRA55() : NumericalEarth.RepeatYearJRA55()
    jra_year  = jra_multiyear ? Dates.year(run_start) : 1990
    jra_start = DateTime(jra_year, Dates.month(run_start), Dates.day(run_start))
    jra_stop  = jra_start + Day(run_days + margin_days)
    @info @sprintf("JRA55 %s: %s .. %s", jra_multiyear ? "MultiYear" : "RepeatYear", jra_start, jra_stop)
    atmos_kw = (; dataset = jra_dataset, start_date = jra_start, end_date = jra_stop,
                  time_indices_in_memory = 24, JRA55_KW...)
    atmosphere = NumericalEarth.JRA55PrescribedAtmosphere(arch; atmos_kw...)
    radiation_component = NumericalEarth.JRA55PrescribedRadiation(arch; atmos_kw...)
end

# --- Initial conditions -------------------------------------------------------
# Sampled from the SAME repaired series the nudging and boundaries use, so IC and forcing agree by
# construction; the alternative re-reads the cache and never sees the repair.
#
# On GPU the restoring series live on the device and this samples them point by point, which the
# device disallows -- so CPU copies are loaded for the purpose, with the same fill_down! repair. They
# are only needed at setup and cost a few tens of MB. On CPU the existing series are reused.
function ic_array(fts, native_grid, t)
    a = zeros(eltype(interior(fts[1])), Nx, Ny, Nz)
    @inbounds for k in 1:Nz, j in 1:Ny, i in 1:Nx
        a[i, j, k] = DW.interpolate_to_grid(fts, i, j, k, native_grid, grid_c, t)
    end
    return a
end

set!(ocean.model.tracers.T, ic_array(T_fts_c, T_fts_c.grid, Time(0.0)))
set!(ocean.model.tracers.S, ic_array(S_fts_c, S_fts_c.grid, Time(0.0)))

free_surface_field(fs) = hasproperty(fs, :η) ? fs.η : fs.displacement
η_field = free_surface_field(ocean.model.free_surface)
let t₀ = Time(0.0)
    ηa = zeros(Nx, Ny, 1)
    for j in 1:Ny, i in 1:Nx
        wet = !immersed_inactive_node(i, j, Nz, grid_c, Center(), Center(), Center())
        ηa[i, j, 1] = ifelse(wet, DW.interpolate_to_grid(zos_fts_c, i, j, 1, zos_fts_c.grid, grid_c, t₀), 0.0)
    end
    set!(η_field, ηa)
end
@info @sprintf("η from GLORYS zos: [%.4f, %.4f] m", minimum(arr(η_field)), maximum(arr(η_field)))

if use_atmos
    coupled_model = NumericalEarth.OceanOnlyModel(ocean; atmosphere, radiation = radiation_component)
    simulation = Simulation(coupled_model; Δt = SPINUP_DT, stop_iteration = SPINUP_ITERS)
else
    simulation = ocean
end
ocean_model_of(m) = hasproperty(m, :ocean) ? (hasproperty(m.ocean, :model) ? m.ocean.model : m.ocean) : m

if use_tracer
    function accumulate_age!(sim)
        m = ocean_model_of(sim.model)
        interior(m.tracers.ac) .+= sim.Δt .* interior(m.tracers.c)
        return nothing
    end
    add_callback!(simulation, accumulate_age!, IterationInterval(1))
end

# --- Diagnostics --------------------------------------------------------------
z_centers = znodes(grid_c, Center(), Center(), Center())
λ_centers = λnodes(grid_c, Center(), Center(), Center())
φ_centers = φnodes(grid_c, Center(), Center(), Center())

@info @sprintf("Δz near surface: %s m (top of %d levels, domain %.0f m)",
               join((@sprintf("%.1f", d) for d in reverse(diff(z_centers))[1:min(6, Nz-1)]), ", "), Nz, depth)

wet3d = [!immersed_inactive_node(i, j, k, grid_c, Center(), Center(), Center()) for i in 1:Nx, j in 1:Ny, k in 1:Nz]
wet_idx = findall(wet3d)
@info "Wet cells: $(length(wet_idx)) of $(Nx*Ny*Nz) ($(round(100*length(wet_idx)/(Nx*Ny*Nz), digits=1))%)"
wet_surface_area = sum(Azᶜᶜᶜ(i, j, Nz, grid_c) for i in 1:Nx, j in 1:Ny
                       if !immersed_inactive_node(i, j, Nz, grid_c, Center(), Center(), Center()); init = 0.0)

# Islands: the topography this domain exists to resolve. Over-smoothing dissolves them silently, so
# report the land-cell count in a 0.4 deg box around each before trusting a run.
for (nm, il, ip) in (("La Palma", -17.87, 28.68), ("El Hierro", -18.01, 27.74), ("Tenerife", -16.60, 28.29),
                     ("Gran Canaria", -15.60, 27.96), ("Fuerteventura", -14.02, 28.36), ("Lanzarote", -13.62, 29.03))
    box = [(i, j) for i in 1:Nx, j in 1:Ny if abs(λ_centers[i]-il) < 0.2 && abs(φ_centers[j]-ip) < 0.2]
    isempty(box) && continue
    nl = count(!wet3d[i, j, Nz] for (i, j) in box)
    nl == 0 ? (@warn "$nm is NOT RESOLVED at this resolution/smoothing -- reduce bathymetry_passes") :
              @info @sprintf("  island %-14s %3d land of %3d cells", nm, nl, length(box))
end

function wet_extrema(a3)          # a3 is a CPU array, not a Field
    a = a3[wet_idx]
    nmin, nmax = argmin(a), argmax(a)
    return a[nmin], wet_idx[nmin], a[nmax], wet_idx[nmax]
end

# Tag by the feature that matters here: which open edge, the upwelling shelf, or the island zone.
function whereis(I)
    i, j, k = Tuple(I)
    λi, φj = λ_centers[clamp(i,1,Nx)], φ_centers[clamp(j,1,Ny)]
    tag = i <= 2                     ? "WEST-bnd"  :
          j <= 2                     ? "SOUTH-bnd" :
          j >= Ny-1                  ? "NORTH-bnd" :
          (λi > -14.0)               ? "SHELF/upwelling" :
          (-18.3 < λi < -13.2 && 27.4 < φj < 29.3) ? "ISLANDS" : "interior"
    return @sprintf("(%3d,%3d,%2d) z=%7.1fm %s", i, j, k, z_centers[clamp(k,1,Nz)], tag)
end

# Volume transport through each open edge, positive INTO the domain.
function boundary_transports(u, v)
    W = S_ = N = 0.0
    @inbounds for k in 1:Nz, j in 1:Ny
        wetcell(1, j, k, grid_c, Face(), Center(), Center()) && (W += u[1, j, k] * Axᶠᶜᶜ(1, j, k, grid_c))
    end
    @inbounds for k in 1:Nz, i in 1:Nx
        wetcell(i, 1, k, grid_c, Center(), Face(), Center())  && (S_ += v[i, 1, k]  * Ayᶜᶠᶜ(i, 1, k, grid_c))
        # Ny+1, not Ny: for Ny cells a y-Face field has Ny+1 points and the north DOMAIN EDGE is the
        # last one. Sampling Ny measures the flux through an interior face instead, which is not the
        # domain's volume budget. (north_V_obc already uses Ny+1.)
        wetcell(i, Ny+1, k, grid_c, Center(), Face(), Center()) && (N -= v[i, Ny+1, k] * Ayᶜᶠᶜ(i, Ny+1, k, grid_c))
    end
    return W, S_, N
end

# What the raw GLORYS exterior would carry through each edge: a parent/child mismatch measure, not a
# physical signal. Reported, never corrected.
#
# Computed live, not precomputed. An earlier version evaluated it once at setup and interpolated in
# time; the parent interpolation is exactly linear so that part was sound, but it also froze the
# SECTION AREAS at their rest state, and under z* those breathe with the free surface. That shifted
# the diagnostic by a few tenths of a percent against the validated baseline. Host copies of the
# parent series keep this off the device without changing what is measured.

function glorys_residuals(t)
    W = S_ = N = 0.0
    @inbounds for k in 1:Nz, j in 1:Ny
        wetcell(1, j, k, grid_c, Face(), Center(), Center()) &&
            (W += DW.interpolate_to_grid(u_fts_c, 1, j, k, u_fts_c.grid, grid_c, t) * Axᶠᶜᶜ(1, j, k, grid_c))
    end
    @inbounds for k in 1:Nz, i in 1:Nx
        wetcell(i, 1, k, grid_c, Center(), Face(), Center()) &&
            (S_ += DW.interpolate_to_grid(v_fts_c, i, 1, k, v_fts_c.grid, grid_c, t) * Ayᶜᶠᶜ(i, 1, k, grid_c))
        wetcell(i, Ny+1, k, grid_c, Center(), Face(), Center()) &&
            (N -= DW.interpolate_to_grid(v_fts_c, i, Ny+1, k, v_fts_c.grid, grid_c, t) * Ayᶜᶠᶜ(i, Ny+1, k, grid_c))
    end
    return W, S_, N
end

# THE signal for this domain. Upwelling is offshore Ekman transport drawing cold water up against the
# African coast, so the diagnostic is the cross-shore SST contrast: coastal minus offshore at MATCHED
# latitude, plus the temperature at the ~100 m source depth over the shelf.
#
# The coast runs diagonally -- λ = -14.87 at 25N to -10.23 at 29.4N -- so a fixed-longitude strip
# cannot follow it. The first version used λ > -13.5 and sampled 0 wet cells at 25-26N and 7 at
# 26-27N, missing Cape Bojador (~26N) entirely, i.e. one of the two upwelling centres in the domain.
# `upwelling_index` therefore walks in from the coast ROW BY ROW. `upwelling_index_fixedλ` keeps the
# original definition so July stays comparable with the January run, which used it.
const k_src = argmin(abs.(z_centers .+ 100))     # ~100 m, the source depth for coastal upwelling
# Bands are defined in DEGREES and converted to cells, not fixed cell counts: a hard-coded 150-cell
# offshore offset silently exceeded the grid_c width at 0.1° and returned NaN for every row.
# Both upwelling indices are restricted to φ >= 25.0 so the number stays comparable when the domain
# is extended south: otherwise a larger domain silently averages in the Cape Bojador upwelling too.
const J_UPW_MIN   = something(findfirst(>=(25.0), φ_centers), 1)
const N_COAST     = max(1, round(Int, 0.50 / Δdeg))   # ~50 km inshore band
const N_OFFSHORE  = round(Int, 3.00 / Δdeg)           # reference centred ~300 km offshore
const N_OFF_HALF  = max(1, round(Int, 0.25 / Δdeg))   # ~25 km half-width on that reference
# Easternmost wet surface cell per row = the coast. NaN where the row has no water.
const i_coast = [findlast(i -> wet3d[i, j, Nz], 1:Nx) for j in 1:Ny]

function upwelling_index(T)
    cT = oT = sT = 0.0; nc = no = ns = 0
    @inbounds for j in J_UPW_MIN:Ny
        ic = i_coast[j]
        ic === nothing && continue
        for i in max(1, ic-N_COAST+1):ic
            wet3d[i, j, Nz] || continue
            cT += T[i, j, Nz]; nc += 1
            if wet3d[i, j, k_src]; sT += T[i, j, k_src]; ns += 1; end
        end
        lo, hi = ic - N_OFFSHORE - N_OFF_HALF, ic - N_OFFSHORE + N_OFF_HALF
        lo < 1 && continue
        for i in lo:hi
            wet3d[i, j, Nz] || continue
            oT += T[i, j, Nz]; no += 1
        end
    end
    Tc = nc > 0 ? cT/nc : NaN
    To = no > 0 ? oT/no : NaN
    return Tc, To, Tc - To, ns > 0 ? sT/ns : NaN
end

# Original fixed-longitude definition, kept only so the January run stays comparable.
function upwelling_index_fixedλ(T)
    cT = oT = 0.0; nc = no = 0
    @inbounds for j in J_UPW_MIN:Ny, i in 1:Nx
        wet3d[i, j, Nz] || continue
        λi = λ_centers[i]
        if λi > -13.5;              cT += T[i, j, Nz]; nc += 1
        elseif -20.5 < λi < -19.0;  oT += T[i, j, Nz]; no += 1 end
    end
    return (nc > 0 ? cT/nc : NaN) - (no > 0 ? oT/no : NaN)
end

# --- Tracer budget ------------------------------------------------------------
# The source is normalised to 1 unit/s, so cumulative input is exactly the elapsed time in seconds:
#   inventory + cumulative export == t
# which makes the budget checkable by inspection rather than by trusting the advection scheme.
# Mean age is ac/c over the whole tracer-bearing volume (Deleersnijder); it is only meaningful once
# there is tracer to speak of, so it is reported as NaN until the inventory is non-trivial.
function tracer_budget(c, ac)
    C = 0.0; AC = 0.0; cmax = -Inf; imax = CartesianIndex(1, 1, 1)
    @inbounds for k in 1:Nz, j in 1:Ny, i in 1:Nx
        wetcell(i, j, k, grid_c, Center(), Center(), Center()) || continue
        v = Vᶜᶜᶜ(i, j, k, grid_c)
        C += v * c[i, j, k]; AC += v * ac[i, j, k]
        if c[i, j, k] > cmax; cmax = c[i, j, k]; imax = CartesianIndex(i, j, k); end
    end
    return C, AC, cmax, imax
end

# Advective tracer flux through each open edge; + = into the domain, so a negative sum is export.
function tracer_export(c, u, v)
    W = 0.0; So = 0.0; N = 0.0
    @inbounds for k in 1:Nz, j in 1:Ny
        wetcell(1, j, k, grid_c, Face(), Center(), Center()) &&
            (W += u[1, j, k] * Axᶠᶜᶜ(1, j, k, grid_c) * c[1, j, k])
    end
    @inbounds for k in 1:Nz, i in 1:Nx
        wetcell(i, 1, k, grid_c, Center(), Face(), Center()) &&
            (So += v[i, 1, k] * Ayᶜᶠᶜ(i, 1, k, grid_c) * c[i, 1, k])
        wetcell(i, Ny+1, k, grid_c, Center(), Face(), Center()) &&
            (N -= v[i, Ny+1, k] * Ayᶜᶠᶜ(i, Ny+1, k, grid_c) * c[i, Ny, k])
    end
    return W, So, N
end

# --- Conservation and stability monitors --------------------------------------
function salt_and_volume(S)
    V = 0.0; SV = 0.0
    @inbounds for k in 1:Nz, j in 1:Ny, i in 1:Nx
        wetcell(i, j, k, grid_c, Center(), Center(), Center()) || continue
        v = Vᶜᶜᶜ(i, j, k, grid_c)
        V += v; SV += v * S[i, j, k]
    end
    return SV, V
end

const eos = SWP.TEOS10.TEOS10EquationOfState()
# Transport split by POTENTIAL density, for comparison with the Canary Current literature.
# The inverse box models (Hernández-Guerra 2005; Pérez-Hernández 2013; Casanova-Masjoan 2020) report
# transports over NEUTRAL density (γⁿ) layers -- thermocline γⁿ < 27.38, intermediate 27.38–27.922,
# deep > 27.922. γⁿ needs the Jackett & McDougall lookup, which this stack does not have, so σθ =
# ρ(T,S,0) − 1000 stands in at the same thresholds. Good at the 27.38 thermocline boundary; it drifts
# with depth, so treat the intermediate/deep split as indicative only.
#
# Fixed-DEPTH bins would be the obvious alternative and are wrong for this purpose: isopycnals heave
# seasonally, so a fixed split aliases that heave into apparent transport change -- corrupting exactly
# the January-vs-July comparison this diagnostic exists to make.
#
# Literature check for the thermocline layer, our sign convention (+ = into the domain): the Canary
# Current flows SOUTHWARD, so expect north > 0 (inflow) and south < 0 (outflow), a few Sv each, weaker
# in January than July. Hernández-Guerra 2005 give 4.7 ± 0.8 Sv southward north of the islands.
# Thresholds match the DEPTH the literature's layers occupy, not their γⁿ numbers -- σθ ≠ γⁿ, so
# copying the numbers copies the label and not the layer. Measured in the January run: σθ=27.38 sits
# at 1195 m against the literature's ~700 m thermocline base, while σθ=27.10 sits at ~500 m. The
# 27.38 bin therefore straddled TWO OPPOSING CURRENTS -- the Canary Current above ~500 m (north
# inflow) and the poleward undercurrent at 700-1500 m (north outflow) -- and as the undercurrent spun
# up it flipped the sum from +8.5 to -0.5 Sv, which looked like the model losing the Canary Current
# when it was the diagnostic mixing two signals. Offline re-binning of the saved sections found it.
#
# The deep bound also cannot use the literature's number: σθ maxes at 27.798 in this domain (coldest,
# saltiest corner T=1.97 S=34.944), so γⁿ=27.922 is unreachable and that bin would always read zero.
const σ_THERMOCLINE  = 27.10     # ≈ 500 m here; the Canary Current layer
const σ_INTERMEDIATE = 27.70     # ≈ 1800 m; NADW
@inline σθ_of(T, S) = SWP.ρ(T, S, 0.0, eos) - 1000.0
@inline density_layer(σ) = σ < σ_THERMOCLINE ? 1 : (σ < σ_INTERMEDIATE ? 2 : 3)

function layered_boundary_transports(T, S, u, v)
    Q = zeros(3, 3)                     # [boundary W/S/N, layer thermocline/intermediate/deep]
    @inbounds for k in 1:Nz, j in 1:Ny
        wetcell(1, j, k, grid_c, Face(), Center(), Center()) || continue
        Q[1, density_layer(σθ_of(T[1, j, k], S[1, j, k]))] += u[1, j, k] * Axᶠᶜᶜ(1, j, k, grid_c)
    end
    @inbounds for k in 1:Nz, i in 1:Nx
        if wetcell(i, 1, k, grid_c, Center(), Face(), Center())
            Q[2, density_layer(σθ_of(T[i, 1, k], S[i, 1, k]))] += v[i, 1, k] * Ayᶜᶠᶜ(i, 1, k, grid_c)
        end
        if wetcell(i, Ny+1, k, grid_c, Center(), Face(), Center())
            Q[3, density_layer(σθ_of(T[i, Ny, k], S[i, Ny, k]))] -= v[i, Ny+1, k] * Ayᶜᶠᶜ(i, Ny+1, k, grid_c)
        end
    end
    return Q
end

function static_instability(T, S)
    n = 0; worst = 0.0
    @inbounds for j in 1:Ny, i in 1:Nx, k in 1:(Nz-1)
        (wetcell(i, j, k,   grid_c, Center(), Center(), Center()) &&
         wetcell(i, j, k+1, grid_c, Center(), Center(), Center())) || continue
        zr = 0.5 * (z_centers[k] + z_centers[k+1])
        Δρ = SWP.ρ(T[i, j, k+1], S[i, j, k+1], zr, eos) - SWP.ρ(T[i, j, k], S[i, j, k], zr, eos)
        if Δρ > 1e-4; n += 1; worst = max(worst, Δρ); end
    end
    return n, worst
end

function surface_freshwater_flux(m)
    use_atmos || return 0.0
    Jη = arr(m.interfaces.net_fluxes.ocean.η)
    Q = 0.0
    @inbounds for j in 1:Ny, i in 1:Nx
        immersed_inactive_node(i, j, Nz, grid_c, Center(), Center(), Center()) && continue
        Q += Jη[i, j, 1] * Azᶜᶜᶜ(i, j, Nz, grid_c)
    end
    return Q
end

function native_range(fts, floor_value)
    a = filter(x -> x > floor_value, vec(Array(interior(fts[1]))))
    isempty(a) ? (NaN, NaN) : (minimum(a), maximum(a))
end
T_ref = native_range(T_fts, 1.0); S_ref = native_range(S_fts, 5.0)
@info @sprintf("Reference (GLORYS native): T ∈ [%.3f, %.3f] | S ∈ [%.3f, %.3f]", T_ref[1], T_ref[2], S_ref[1], S_ref[2])
let (lo,_,hi,_) = wet_extrema(arr(ocean.model.tracers.T)), (slo,_,shi,_) = wet_extrema(arr(ocean.model.tracers.S))
    @info @sprintf("Initial condition (wet):    T ∈ [%.3f, %.3f] | S ∈ [%.3f, %.3f]", lo, hi, slo, shi)
end
let (W, S_, N) = glorys_residuals(Time(0.0))
    @info @sprintf("GLORYS residuals at t=0 (+ = inflow): west %+.3e  south %+.3e  north %+.3e  net %+.3e m³/s; NOT corrected",
                   W, S_, N, W + S_ + N)
end

# Calibrate the proxy: report where the σθ surfaces sit in depth, so the offset from the literature's
# γⁿ surfaces is visible rather than assumed. The papers put γⁿ=26.85 at roughly 80-320 m; if our
# σθ=27.38 surface sits much deeper than their ~700 m thermocline base, our "thermocline" layer is
# including intermediate water and will read high against theirs.
let T = arr(ocean.model.tracers.T), S = arr(ocean.model.tracers.S)
    for σ in (26.85, σ_THERMOCLINE, σ_INTERMEDIATE)
        d = Float64[]
        for j in 1:Ny, i in 1:Nx
            k = findlast(k -> wet3d[i, j, k] && σθ_of(T[i, j, k], S[i, j, k]) >= σ, 1:Nz)
            k === nothing && continue
            push!(d, -z_centers[k])
        end
        isempty(d) ? @info(@sprintf("σθ = %.2f surface: absent from this domain", σ)) :
            @info @sprintf("σθ = %.2f surface sits at %.0f m (median over wet columns; range %.0f–%.0f)",
                           σ, median(d), minimum(d), maximum(d))
    end
end

const salt₀ = Ref(0.0); const vol₀ = Ref(0.0)
let Ta0 = arr(ocean.model.tracers.T), Sa0 = arr(ocean.model.tracers.S)
    SV, V = salt_and_volume(Sa0)
    n, w = static_instability(Ta0, Sa0)
    Tc, To, ΔT, Ts = upwelling_index(Ta0)
    @info @sprintf("At rest: salt %.6e | volume %.6e | mean S %.4f | unstable pairs %d (worst Δρ %.3f)",
                   SV, V, SV/V, n, w)
    @info @sprintf("At rest: SST coastal %.3f, offshore %.3f, ΔT %+.3f °C | T at %.0f m over shelf %.3f °C",
                   Tc, To, ΔT, z_centers[k_src], Ts)
end

const wall_t0 = Ref(Base.time()); const wall_last = Ref(Base.time()); const sim_last = Ref(0.0)

function progress(sim)
    model = ocean_model_of(sim.model)
    sync_grid_c!()
    # ONE device->host copy per field per block, then everything below is plain array work. This is
    # what makes the script GPU-safe: the diagnostics all index cell by cell, which a device array
    # will not allow. ~17 MB per field at 0.02°, every 10 iterations -- negligible against the step.
    ua, va, wa = arr(model.velocities.u), arr(model.velocities.v), arr(model.velocities.w)
    Ta, Sa = arr(model.tracers.T), arr(model.tracers.S)
    ηa = arr(free_surface_field(model.free_surface))
    Tmin, _, Tmax, _ = wet_extrema(Ta)
    Smin, Imin, Smax, Imax = wet_extrema(Sa)
    _, Iw = findmax(abs, wa); _, Iu = findmax(abs, ua)

    @info @sprintf("iter %d | t = %s | Δt = %s | max|u,v,w| = (%.2e, %.2e, %.2e) | max|η| = %.3e | T ∈ [%.3f, %.3f] | S ∈ [%.3f, %.3f] (excursion %+.3f / %+.3f)",
                   iteration(sim), prettytime(time(sim)), prettytime(sim.Δt),
                   maximum(abs, ua), maximum(abs, va), maximum(abs, wa), maximum(abs, ηa),
                   Tmin, Tmax, Smin, Smax, Smin - S_ref[1], Smax - S_ref[2])
    @info @sprintf("        where: Smax %s | Smin %s | max|w| %s | max|u| %s",
                   whereis(Imax), whereis(Imin), whereis(Iw), whereis(Iu))

    W, So, N = boundary_transports(ua, va)
    Wr, Sr, Nr = glorys_residuals(Time(model.clock.time))
    @info @sprintf("        transport (+ = inflow): west %+.3e  south %+.3e  north %+.3e | net %+.3e m³/s -> %+.3e m/hour",
                   W, So, N, W+So+N, (W+So+N) / wet_surface_area * 3600)
    @info @sprintf("        GLORYS residual (uncorrected): west %+.3e  south %+.3e  north %+.3e | net %+.3e m³/s",
                   Wr, Sr, Nr, Wr+Sr+Nr)

    Q = layered_boundary_transports(Ta, Sa, ua, va)
    for (li, lab) in ((1, "Canary Cur σθ<27.10"), (2, "undercurrent  <27.70"), (3, "deep      σθ≥27.70"))
        @info @sprintf("        %s: west %+7.3f  south %+7.3f  north %+7.3f | net %+7.3f Sv",
                       lab, Q[1,li]/1e6, Q[2,li]/1e6, Q[3,li]/1e6, sum(Q[:,li])/1e6)
    end

    Tc, To, ΔT, Ts = upwelling_index(Ta)
    @info @sprintf("        UPWELLING: SST coastal %.3f  offshore %.3f  ΔT %+.3f °C | T at %.0f m over shelf %.3f °C | ΔT(fixed-λ, Jan-comparable) %+.3f",
                   Tc, To, ΔT, z_centers[k_src], Ts, upwelling_index_fixedλ(Ta))

    if use_tracer
        ca, aca = arr(model.tracers.c), arr(model.tracers.ac)
        C, AC, cmax, imax = tracer_budget(ca, aca)
        Wt, St, Nt = tracer_export(ca, ua, va)
        elapsed = max(time(sim), eps())
        @info @sprintf("        TRACER: inventory %.4e of %.4e released (%.2f%%) | mean age %s | max c %.3e %s",
                       C, elapsed, 100C/elapsed,
                       C > 1e-6*elapsed ? @sprintf("%.3f d", AC/C/86400) : "n/a (too little tracer)",
                       cmax, whereis(imax))
        @info @sprintf("        TRACER export (+ = in): west %+.3e  south %+.3e  north %+.3e | net %+.3e /s",
                       Wt, St, Nt, Wt+St+Nt)
    end

    SV, V = salt_and_volume(Sa)
    ninv, worstΔρ = static_instability(Ta, Sa)
    if vol₀[] == 0.0 && iteration(sim) > 0
        salt₀[], vol₀[] = SV, V
        @info @sprintf("        drift baseline set at t = %s (σ reconciled)", prettytime(time(sim)))
    end
    if use_atmos
        Qfw = surface_freshwater_flux(sim.model)
        @info @sprintf("        surface freshwater: %+.4e m³/s (%+.3f m/yr basin mean)",
                       Qfw, Qfw / wet_surface_area * 3.15576e7)
    end
    have = vol₀[] > 0
    @info @sprintf("        conservation: salt %s | volume %s | mean S %.4f | unstable pairs %d (worst Δρ %.3f)",
                   have ? @sprintf("%+.3e", (SV-salt₀[])/salt₀[]) : "(pending)",
                   have ? @sprintf("%+.3e", (V-vol₀[])/vol₀[])   : "(pending)",
                   SV/V, ninv, worstΔρ)

    now = Base.time(); wall = now - wall_t0[]
    Δwall, Δsim = now - wall_last[], time(sim) - sim_last[]
    wall_last[], sim_last[] = now, time(sim)
    rate = Δwall > 0 ? Δsim/Δwall : 0.0
    eta = (isfinite(sim.stop_time) && rate > 0) ? (sim.stop_time - time(sim))/rate : NaN
    @info @sprintf("        pacing: wall %s | %.0f× realtime | %.2f s/iter | ETA %s",
                   prettytime(wall), rate, Δwall/PROGRESS_ITERS, isnan(eta) ? "n/a" : prettytime(eta))
    flush(stderr)
    return nothing
end
# progress copies every 3-D field to the host and loops over all wet cells, so it is worth thinning
# for long runs (progress_iters).
add_callback!(simulation, progress, IterationInterval(PROGRESS_ITERS))

# --- Spin-up then the run -----------------------------------------------------
# Output filenames carry run_tag verbatim. It used to be ASSEMBLED from the parameters themselves,
# which made the filename a second, lossy record of the configuration -- and a silent one, since
# any parameter missing from the recipe (the release site, the sponge width) left no trace in the
# name. The config is the record now; the tag only has to be unique.
tag = get(ENV, "CANARY_TAG", "_ref90")

# Restart. Oceananigans' `set!(model, file)` cannot read these checkpoints: the Checkpointer serialises
# whatever model the simulation holds, which here is the coupled EarthSystemModel, and neither
# `checkpointer_address` nor `set!` has a method for it. The state is all there, so read it directly.
# Not restored: timestepper tendencies (the checkpoint stores none), so the first step after a restart
# is effectively first-order -- one step per restart, negligible over months.
function restore_ocean_state!(model, path)
    om = "simulation/model/ocean/model"
    JLD2.jldopen(path, "r") do file
        load(f, key) = copyto!(parent(f), on_architecture(architecture(model), file["$om/$key/data"]))
        for name in (:u, :v, :w)
            load(getproperty(model.velocities, name), "velocities/$name")
        end
        for name in keys(model.tracers)
            string(name) in keys(file["$om/tracers"]) ?
                load(model.tracers[name], "tracers/$name") :
                @warn "checkpoint has no tracer $name -- left at its initial condition"
        end
        fs = model.free_surface
        load(free_surface_field(fs), "free_surface/displacement")
        if hasproperty(fs, :barotropic_velocities)
            load(fs.barotropic_velocities.U, "free_surface/barotropic_velocities/U")
            load(fs.barotropic_velocities.V, "free_surface/barotropic_velocities/V")
        end
        ck = file["$om/clock"]
        model.clock.time, model.clock.iteration, model.clock.last_Δt = ck.time, ck.iteration, ck.last_Δt
    end
    update_state!(model)
    return nothing
end

# Restart. The checkpointer writes the OCEAN model, while `simulation` wraps the coupled model, so
# `run!(pickup=true)` would target the wrong object: restore the ocean model explicitly and carry its
# clock over to the coupled model. pickup: "auto" takes the newest checkpoint for this tag.
pickup_request = ""
pickup_file = if pickup_request == "auto"
    cands = filter(f -> startswith(basename(f), "canary_checkpoint$(tag)_iteration") && endswith(f, ".jld2"),
                   readdir(output_dir, join = true))
    isempty(cands) ? "" : sort(cands, by = mtime)[end]
else
    pickup_request
end

if isempty(pickup_file)
    simulation.Δt = SPINUP_DT
    simulation.stop_iteration = SPINUP_ITERS
    run!(simulation)
else
    isfile(pickup_file) || error("pickup file not found: $pickup_file")
    restore_ocean_state!(ocean.model, pickup_file)
    t0, it0 = ocean.model.clock.time, ocean.model.clock.iteration
    if use_atmos     # the coupled model keeps its own clock; without this the run restarts at t = 0
        simulation.model.clock.time = t0
        simulation.model.clock.iteration = it0
    end
    simulation.Δt = ocean.model.clock.last_Δt
    @info @sprintf("Picked up from %s: t = %s, iteration %d", basename(pickup_file), prettytime(t0), it0)
end

ocean_timescale(model) = cell_advection_timescale(ocean_model_of(model))
wizard = TimeStepWizard(; cfl = Float64(0.15), max_Δt = 15minutes, max_change = Float64(1.1),
                          cell_advection_timescale = ocean_timescale)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))
simulation.stop_iteration = Inf
simulation.stop_time = Float64(90) * days

output_fields   = merge(ocean.model.velocities, ocean.model.tracers)
# Output cadence. 6-hourly suits a 7-day run; a months-long tracer run wants daily or coarser, since
# volume scales with run length (output_hours, snapshot_hours for the 3-D snapshots).
output_hours   = Float64(6)
snapshot_hours = Float64(72)
output_schedule = TimeInterval(output_hours * hours)
@info @sprintf("Output: sections every %.1f h, 3-D snapshots every %.1f h", output_hours, snapshot_hours)
out(name) = joinpath(output_dir, "canary_glorys_$(name)$(tag)")

# Oceananigans renamed this writer keyword in 0.113 (overwrite_existing -> overwrite_files), for
# Checkpointer as well as JLD2Writer. Picking it from the method's own kwarg list keeps ONE code path
# across both versions; splat it into every writer below rather than naming the keyword directly.
const OVERWRITE = (:overwrite_files in Base.kwarg_decl(which(JLD2Writer, Tuple{Any, Any}))) ?
                  (; overwrite_files = true) : (; overwrite_existing = true)

simulation.output_writers[:surface] = JLD2Writer(ocean.model, output_fields;
    schedule = output_schedule, OVERWRITE..., indices = (:, :, Nz), filename = out("surface"))
simulation.output_writers[:west_section] = JLD2Writer(ocean.model, output_fields;
    schedule = output_schedule, OVERWRITE..., indices = (1, :, :), filename = out("west_section"))
# Zonal section at 26.5N, cutting the shelf, the upwelling front and the open ocean.
const j_shelf = argmin(abs.(φ_centers .- 26.5))
simulation.output_writers[:shelf_section] = JLD2Writer(ocean.model, output_fields;
    schedule = output_schedule, OVERWRITE..., indices = (:, j_shelf, :), filename = out("shelf_section"))
# All three open boundaries are saved, so the layered transports can be RE-BINNED offline: the layer
# definition (σθ as a γⁿ proxy, and where to put the thresholds) is the least settled part of this
# diagnostic, and re-deriving it from file is free where re-running is not. T, S and the grid all
# come back from these files, so σθ and area-weighted transport are fully recoverable.
#
# The north needs TWO writers. Face-in-y indexing puts the domain edge at Ny+1 for v while T and S
# live at Ny, and one slice cannot hold both -- `(:, Ny, :)` would silently give the last INTERIOR
# face, which is the off-by-one this script had until it was caught. West and south are exempt:
# index 1 is both the boundary face and the first interior cell.
simulation.output_writers[:south_section] = JLD2Writer(ocean.model, output_fields;
    schedule = output_schedule, OVERWRITE..., indices = (:, 1, :), filename = out("south_section"))
simulation.output_writers[:north_section] = JLD2Writer(ocean.model, output_fields;
    schedule = output_schedule, OVERWRITE..., indices = (:, Ny, :), filename = out("north_section"))
simulation.output_writers[:north_face_v] = JLD2Writer(ocean.model, (; v = ocean.model.velocities.v);
    schedule = output_schedule, OVERWRITE..., indices = (:, Ny+1, :), filename = out("north_face_v"))

# 3-D snapshot every snapshot_hours (~83 MB/frame at 0.02 deg, so ~2.5 GB for 90 days at 72 h).
# Cheap against the cost of re-running to answer a question this run could already have answered:
# interior sections, water-mass transformation, volume-integrated budgets.
simulation.output_writers[:snapshot3d] = JLD2Writer(ocean.model, output_fields;
    schedule = TimeInterval(snapshot_hours * hours), OVERWRITE..., filename = out("snapshot3d"))

simulation.output_writers[:free_surface] = JLD2Writer(ocean.model, (; η = free_surface_field(ocean.model.free_surface));
    schedule = output_schedule, OVERWRITE..., filename = out("free_surface"))
# One checkpoint on disk at a time (cleanup), so a months-long run can be resumed after a crash or a
# deliberate stop without accumulating ~560 MB files. checkpoint_days sets the cadence.
checkpoint_days = Float64(10)
simulation.output_writers[:checkpointer] = Checkpointer(ocean.model;
    schedule = TimeInterval(checkpoint_days * days), dir = output_dir, prefix = "canary_checkpoint$(tag)",
    OVERWRITE..., cleanup = true)
@info @sprintf("Checkpoint every %.1f days -> %s", checkpoint_days, "canary_checkpoint$(tag)_iteration*.jld2")

run!(simulation)
