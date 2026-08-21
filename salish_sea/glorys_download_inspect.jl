# glorys_download_inspect.jl
#
# Stage-2 step 0: download a single-day GLORYS snapshot over the Salish Sea
# open-boundary region and eyeball it, BEFORE wiring GLORYS into any boundary
# condition. Goal: confirm (a) Copernicus Marine credentials work, and (b) the
# regional subset covers the domain with sane values/coordinates.
#
# GLORYS = GLOBAL_MULTIYEAR_PHY 1/12° (~8 km), daily, 1993-01-01 .. 2021-06-30.
# Fields pulled match exactly what the OBCs + nudging band will need in Stage 2:
#   :temperature :salinity :u_velocity :v_velocity :free_surface
#
# NOTE on API: the shipped `examples/download_glorys_data.jl` references
# `DataWrangling.Copernicus.GLORYSDaily()`, but NumericalEarth v0.5.8 (this
# project's pin) has no `Copernicus` submodule -- GLORYS is `using`-ed straight
# into DataWrangling, so the correct name is `DataWrangling.GLORYSDaily()`.

using NumericalEarth
using Oceananigans
using CopernicusMarine
using CairoMakie
using Dates
using Printf
using Statistics
using Downloads   # NumericalEarth extends `Downloads.download(::Metadata)`, not `Base.download`

const DW = NumericalEarth.DataWrangling

# Node helpers. Two gotchas in this Oceananigans version:
#   1. `xnodes(::Field)` forwards an `indices` kwarg the grid method rejects.
#   2. On a LatitudeLongitudeGrid, `xnodes`/`ynodes` return *physical metres* (x is even
#      2-D, since Δx = R·cosφ·Δλ). For lon/lat degrees use `λnodes`/`φnodes` (1-D).
using Oceananigans.Fields: location
hx(f) = λnodes(f.grid, map(L -> L(), location(f))...)
hy(f) = φnodes(f.grid, map(L -> L(), location(f))...)

# --- 0. Credential check -----------------------------------------------------
# CopernicusMarine.jl reads COPERNICUSMARINE_SERVICE_USERNAME / _PASSWORD, or a
# stored login from `CopernicusMarine.login(username=..., password=...)`.
have_user = haskey(ENV, "COPERNICUSMARINE_SERVICE_USERNAME")
have_pass = haskey(ENV, "COPERNICUSMARINE_SERVICE_PASSWORD")
if have_user && have_pass
    @info "Copernicus credentials found in ENV (username=$(ENV["COPERNICUSMARINE_SERVICE_USERNAME"]))."
else
    @warn """Copernicus credentials NOT found in ENV.
    If you have not run a persistent `CopernicusMarine.login(...)`, the downloads below will fail with an auth error.
    Fix (either one):
      1. In this shell before launching Julia:
           export COPERNICUSMARINE_SERVICE_USERNAME=your_user
           export COPERNICUSMARINE_SERVICE_PASSWORD=your_pass
      2. Once, from Julia (stores credentials persistently, no ENV needed after):
           using CopernicusMarine; CopernicusMarine.login(username="your_user", password="your_pass")
    """
end

# --- 1. Region + date --------------------------------------------------------
# Domain from salish_simple_obc.jl: λ=(235.5, 238), φ=(47, 49.75)  [0-360° convention].
# Pad outward ~0.25° (~3 GLORYS cells) so boundary interpolation has exterior data on all sides.
pad = 0.25
λ₁, λ₂ = 235.5, 238.0
φ₁, φ₂ = 47.0, 49.75
region = DW.BoundingBox(longitude = (λ₁ - pad, λ₂ + pad),
                        latitude  = (φ₁ - pad, φ₂ + pad))

dataset = DW.GLORYSDaily()
date    = DateTime(1993, 1, 1)   # matches start_date in the main run

variables = (:temperature, :salinity, :u_velocity, :v_velocity, :free_surface)

@info "Region" longitude=(λ₁ - pad, λ₂ + pad) latitude=(φ₁ - pad, φ₂ + pad) date

# --- 2. Download + load each field, print diagnostics ------------------------
fields = Dict{Symbol,Any}()
for var in variables
    try
        meta = DW.Metadatum(var; dataset, region, date)
        path = Downloads.download(meta)
        @info "Downloaded $var" path
        f = Field(meta, inpainting = nothing)   # land stays NaN/missing, not filled
        fields[var] = f

        vals   = filter(isfinite, vec(Array(interior(f))))
        λs, φs = hx(f), hy(f)
        nwet   = length(vals)
        ntot   = length(interior(f))
        @info @sprintf("  %-13s size=%s  wet=%d/%d (%.0f%%)  range=[%.4g, %.4g]  mean=%.4g",
                       var, string(size(f)), nwet, ntot, 100nwet/ntot,
                       isempty(vals) ? NaN : minimum(vals),
                       isempty(vals) ? NaN : maximum(vals),
                       isempty(vals) ? NaN : mean(vals))
        @info @sprintf("  %-13s λ ∈ [%.3f, %.3f] (N=%d)   φ ∈ [%.3f, %.3f] (N=%d)",
                       "", first(λs), last(λs), length(λs),
                       first(φs), last(φs), length(φs))
    catch err
        @error "Failed on $var" exception=(err, catch_backtrace())
    end
end

# --- 3. Quick surface maps to eyeball the subset -----------------------------
# Surface is the top vertical level (Oceananigans z increases upward -> k=end).
surf(f) = ndims(interior(f)) == 3 && size(interior(f), 3) > 1 ?
          Array(interior(f))[:, :, end] : dropdims(Array(interior(f)), dims=3)

try
    fig = Figure(size = (1100, 800))
    panels = [(:temperature, "surface T (°C)"),
              (:salinity,    "surface S (psu)"),
              (:free_surface,"η (m)")]
    for (n, (var, label)) in enumerate(panels)
        haskey(fields, var) || continue
        f = fields[var]
        ax = Axis(fig[fld1(n,2), mod1(n,2)]; title = label, xlabel = "λ", ylabel = "φ")
        hm = heatmap!(ax, hx(f), hy(f), surf(f); nan_color = :gray70)
        Colorbar(fig[fld1(n,2), mod1(n,2), Right()], hm)
    end
    # surface current speed, if both velocity components loaded
    if haskey(fields, :u_velocity) && haskey(fields, :v_velocity)
        u, v = fields[:u_velocity], fields[:v_velocity]
        spd  = sqrt.(surf(u).^2 .+ surf(v).^2)
        ax = Axis(fig[2, 2]; title = "surface speed (m/s)", xlabel = "λ", ylabel = "φ")
        hm = heatmap!(ax, hx(u), hy(u), spd; nan_color = :gray70)
        Colorbar(fig[2, 2, Right()], hm)
    end
    out = "glorys_inspect_$(Dates.format(date, "yyyymmdd")).png"
    save(out, fig)
    @info "Saved surface maps -> $out"
catch err
    @error "Plotting failed (data may still have downloaded fine)" exception=(err, catch_backtrace())
end
