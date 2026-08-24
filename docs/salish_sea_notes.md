# Salish Sea / Juan de Fuca configuration — notes

Status as of 2026-08-21. Model: `salish_sea/salish_glorys_obc.jl`.
For the portable boundary-condition recipe see `obc-implementation-recipe.md`.

The purpose of this configuration was to **test the open-boundary implementation**, not to produce a
production Salish Sea hindcast. GLORYS was chosen deliberately as a simple parent product; LiveOcean
is the intended parent eventually. Read the limitations section before using any numbers.

## Configuration

| | |
|---|---|
| Domain | 235.5–238.0°E, 47.28–49.05°N (−124.5 to −122.0, Juan de Fuca focused) |
| Horizontal | 250 × 177 at 0.01° ≈ 0.74 × 1.11 km (`SALISH_RES_FACTOR` coarsens) |
| Vertical | 36 levels, 400 m, `ExponentialDiscretization(scale=150)` → 2.3 m surface cell, 13 levels in the top 50 m |
| Bathymetry | ETOPO via `regrid_bathymetry`, `minimum_depth=10`, `interpolation_passes=10`, `GridFittedBottom` |
| Open boundaries | west (Juan de Fuca) and north (Strait of Georgia at 49.05) |
| Parent | GLORYS 1/12° daily — T, S, u, v for boundaries and nudging; `zos` for the Flather exterior η |
| Nudging | `DatasetRestoring`, 12 h timescale, Gaussian mask, `sponge_width = 0.15°` |
| Advection | WENO(7) tracers, WENOVectorInvariant momentum (`ocean_simulation` defaults) |
| Timestepping | `TimeStepWizard(cfl=0.15, max_Δt=5min, max_change=1.1)`; Δt settles ≈ 3.4 min |
| Cost | ~12 s/iter, **~1.6 h wall per simulated day** on CPU |

Env flags: `SALISH_RES_FACTOR`, `SALISH_STOP_DAYS`, `SALISH_ZSTAR`, `SALISH_PARTIAL_CELLS`,
`SALISH_WENO_ORDER`. Outputs are tagged by non-default flags so baselines survive.

## Validation

A **7-day run at full resolution** (3460 iterations, 12.2 h wall) is stable and closes its budget:

| | value |
|---|---|
| S excursion beyond the GLORYS range | −0.166, flat for the whole run |
| max\|w\| | 1.35e-3 m/s, steady |
| max\|u\| | 0.39 m/s |
| net boundary transport | −420 m³/s |
| implied sea-level drift | −2.0e-4 m/hour (≈3 cm over 7 days) |
| boundary η vs GLORYS zos | −0.005 m |

Circulation is physically coherent: water enters through the northern boundary (+5.8e3 m³/s) and
leaves through Juan de Fuca (−6.8e3 m³/s), which is what should happen with the Fraser outside the
domain. Deep JdF inflow salinity `Sin` converged to **33.19** against an observed ~33.5.

## The instability, and its root cause

Early configurations blew up: S reached 36.3 and max|w| 5.5e-2 within an hour, growing exponentially
from a deep-basin cell. The **sole root cause was a density inversion in the initial condition**.

GLORYS has **no valid data below −320.6 m anywhere in this region**; deeper levels come back
zero-filled. Where the model bathymetry is deeper than that, interpolation blends a valid value with
zero and produces anything in between — observed S = 30.16 sitting under S = 33.92, and 15.24 once
`PartialCellBottom` activated deeper cells. That is 14 psu of fresh water under salty water,
convecting. The fix is `fill_down!` (vertical constant extrapolation below the deepest valid level)
applied to the restoring series, with initial conditions built by sampling that same repaired series
via `ic_from_fts!` rather than through an independent `set!` path.

Two things were ruled out by direct test along the way and are **not** contributing: tracer advection
order (WENO 5 vs 7 gave identical growth) and stair-step topography (`PartialCellBottom` made no
difference once the IC was repaired). `GridFittedBottom` is retained.

Other fixes made during the investigation, each independently worthwhile but none of them the
destabiliser: the vertical grid rebuild (the old 1000 m domain put half the levels below the deepest
seafloor and left a 14.8 m surface cell), GLORYS `zos` as the Flather exterior η, the net-transport
correction, and initializing η from `zos` instead of rest.

## Domain design

The northern boundary at **49.05°N** was chosen against measured alternatives:

- The original 49.75°N edge was **entirely land** — closed by `interpolation_passes=10`, not by
  geography. At 1–3 passes it opens (19–25 wet cells, to 287 m).
- 49.40–49.55 gives 2–3 disconnected segments (Texada Island) and touches the west boundary corner.
- **49.05**: single clean channel, 81 wet cells, 139 m deep at 10 passes, ~60 km clear of the west
  boundary. Critically it stays **above the 320 m GLORYS validity depth** — a 49.30 cut is 367 m and
  would prescribe extrapolated values across the deepest 13% of the section.
- Southern edge raised 47.0 → 47.28: the old rows 1–30 were 100% land.

`sponge_width` was cut 1.0° → 0.15°: at 1.0° the nudging mask was 0.80 in the deep Strait of Georgia
and never below 0.21 anywhere — near-global restoring rather than a boundary sponge. At 0.15° with a
49.05 boundary, Boundary Pass sits at mask 0.039 (306 h effective) and Juan de Fuca at zero.

**Cost of the choice**: the Fraser (mouth ≈49.12°N) is outside the domain, and the domain lost ~44%
of its water relative to the 49.75 version. Residence time and the freshwater budget of what remains
are correspondingly different objects.

## Known limitations

**GLORYS surface salinity is ~3–4 psu too fresh here.** Domain-wide surface S spans 26.95–28.87
against observed JdF values of ~31–32. This is **inherited from the parent, not generated by the
model** — at day 0.125 the model matches GLORYS to 0.08 psu at every wet point. At 1/12° GLORYS
cannot resolve 5–20 km channels with sills and smears Salish Sea freshwater across Juan de Fuca.
GLORYS is a poor parent for a surface-salinity study; LiveOcean or SalishSeaCast output would be
better. The model does add a further slow freshening of 0.26–0.79 psu over 7 days.

**No tides.** Qin ≈ 5.0e4 m³/s against MacCready et al. (2021, doi:10.1029/2020JC016738) reporting
~1.4e5 m³/s at the system mouth, while ΔS = 3.41 against an observed ~1.5–2. Those are the same
error: Knudsen requires Qin·Sin ≈ |Qout·Sout|, so a weak exchange must carry a large property
contrast. The signature is **insufficient vertical mixing**, and the cause is the absence of tidal
mixing at Haro Strait and Boundary Pass. Note their Qin is TEF (salinity coordinates, tidally
averaged) while ours is Eulerian, so some of the gap is method rather than physics.

**Not equilibrated at 7 days.** Qin had plateaued but Qout, ΔS, the Knudsen residual (−93% → −18%)
and the net through-flow were all still trending. 14–30 days would be needed for equilibrium.

**No rivers.** `net_transport_target = 0` is correct only because of this. NumericalEarth ≥ 0.6.1
provides `FreshwaterExchange` for adding them; the target must then become a net outflow equal to
the discharge.

## Next steps

1. **Tides (TPXO/FES).** The Flather BC already takes the right form — `η_ext += Σ A cos(ωt−φ)`,
   `U_ext += Σ U cos(...)` — so the boundary machinery needs no change. But four things cascade:
   Δt drops ~3.7×, the nudging would damp the tides, the transport correction must become tidally
   averaged, and the exchange diagnostic needs TEF or tidal averaging. Data needs a licence request
   (OSU for TPXO, AVISO for FES). See the vertical-coordinate note below for z\*.
2. **Back to LiveOcean** for boundary and initial conditions — hourly output is available but
   expensive for extended runs. This is what SalishSeaCast does for its own JdF boundary.
3. **z\*** is implemented behind `SALISH_ZSTAR=1` and mechanically works (~6% overhead), but was
   validated against the older Oceananigans and shows a small unexplained S divergence over an hour.
   Re-check before relying on it.

## Vertical coordinate: when z\* is worth it

Correcting an earlier claim in these notes: a static grid in Oceananigans is a **linear free
surface** — cell thicknesses are fixed and η never changes them (`reconcile_vertical_coordinate!`
is a no-op for `ZCoordinate`, `∂t_σ = 0`). It does **not** fail when η exceeds the top cell, and
there is no vanishing-cell instability. The cost is accuracy, not viability.

The quantity that decides it is **max|η| / Δz_top**, because with fixed volumes the surface tracer
concentration is wrong by O(η/Δz_top) — a rising free surface should dilute the surface layer and
does not — and surface fluxes are spread over Δz_top rather than the true Δz_top + η. Below ~0.1 the
error is negligible; 0.1–0.3 is a modest bias; above ~0.5 the surface budget is materially wrong;
above 1 it is meaningless at the surface, though it still runs.

So **thin surface cells select z\*, not tides**. Tides are just the extreme case: η of 0.3–0.5 m is
routine from steric, mesoscale and boundary-imposed signal alone. The Salish configuration runs
static z at Δz_top = 2.3 m with η of 0.2–0.5 m, a ratio of 0.09–0.22, which is fine. A ~1 m surface
cell anywhere would not be.

Two supporting signals: `default_vertical_coordinate` returns `ZStarCoordinate` whenever the grid is
mutable, and NumericalEarth's `FreshwaterExchange` exists specifically to cancel the z\* ambient-carry
term for surface tracer fluxes — the ecosystem treats z\* as the expected mode for realistic surface
forcing, with static z as the fallback.

## Reference models

- **LiveOcean** (UW, ROMS): covers OR/WA/Vancouver Island coastal waters plus the Salish Sea, open
  boundaries in the **open Pacific** forced by HYCOM; ~500 m in the Salish Sea, 3 km at the
  boundaries; 30 terrain-following levels; 45 rivers, 8 tidal constituents.
- **SalishSeaCast** (UBC, NEMO): the closer analogue — two open boundaries, west through Juan de Fuca
  (forced by LiveOcean) and north through **Johnstone Strait** (a monthly T/S climatology). Grid is
  curvilinear/rotated, 398 × 898 at ~500 m, aligned to the Salish Sea axis. Bathymetry is
  deliberately smoothed at the open boundaries.

Neither cuts mid-Strait-of-Georgia, so the 49.05 boundary has no published precedent to validate
against. SalishSeaCast's northern extent is plausibly driven by its operational remit for BC coastal
waters rather than by numerics — the tell being that it forces that boundary with a climatology.
