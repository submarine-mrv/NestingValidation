# Red Sea + Gulf of Aden configuration

Second OBC test domain, adapted from `salish_sea/`. **Status: written and smoke-tested at coarse
resolution; not yet validated.** Model: `redsea_glorys_obc.jl`.

## Why this domain

The Salish Sea tested the boundary machinery but could never test the thing nesting is *for*:
channels there were 5–20 km against a 9 km parent, so GLORYS could not represent the geometry it was
forcing. Here the open boundary sits in deep open ocean the parent resolves well, while
**Bab el-Mandeb (~36 km, ~4 GLORYS cells) is inside the domain** for the model to resolve.

Selected over Gulf of Mexico (Mississippi freshwater, chaotic Loop Current makes attribution hard),
North Sea (tide-dominated shallow shelf — every Salish difficulty at larger scale) and Canary
upwelling (cleanest boundaries but needs wind forcing this config has none of).

## Verified geometry

All measured from ETOPO via `regrid_bathymetry`, not assumed:

| | |
|---|---|
| Domain | 32–50°E, 10–30°N |
| Open boundaries | **east edge only** — south (10.0N), north (30.0N) and west (32.0E) are all land |
| East section at 50°E | single segment, 334 cells over 11.52–14.85°N, 2848 m deep, **719 km²**, ~37 GLORYS cells across |
| Hanish sill | **136 m at 13.68°N** (literature ~137 m) |
| Narrowest main channel | **35.8 km at 12.69°N**, ~4 GLORYS cells |
| Max depth | 3176 m, at 49.58°E / 13.18°N — in the Gulf of Aden beside the boundary |
| Red Sea proper | max 2478 m, **median 425 m** |
| Wet fraction | 16.1% (Salish Sea was 6.4%) |

For scale: the Salish Sea west boundary was **1.6 km²** in 1–2 parent cells. This one is 450× larger
and in water GLORYS resolves.

Suez and Aqaba terminate naturally inside the domain, so no artificial closure is needed.

## GLORYS validity — the check that mattered

The Salish Sea instability came from GLORYS having no valid data below 320 m there while the
bathymetry reached 353 m. Same check here, over a box straddling the boundary (49–51°E, 10.5–16°N):

```
DEEPEST GLORYS LEVEL WITH VALID DATA: -3226.1 m     (100% valid at every level above)
domain max depth                    : -3176 m
```

**50 m of margin**, and the profile is physically right — a salinity maximum of 35.50–36.37 at
460–770 m, which is the Red Sea Outflow Water signature spreading into the Gulf of Aden. So the
parent supplies the right water masses at the boundary, not merely numbers.

**Consequence: `depth` must stay ≤ ~3200 m.** Rounding it up would put model cells below the parent's
valid data and recreate the Salish failure exactly.

## Resolution and cost (CPU)

| Δ | km at 20°N | columns | 7-day run |
|---|---|---|---|
| 0.05° | 5.23 | 144,000 | ~28 h |
| 0.025° | 2.62 | 576,000 | **~8–9 days** |

2.5 km is the right target but needs GPU or HPC. Develop at 0.05° (`REDSEA_DELTA`).

## Vertical grid

`Nz = 32`, `depth = 3200`, `z_scale = 480` → **0.94 m surface cell, 17 levels in the top 150 m**.

The top 150 m is what matters: the sill is 136 m, so the whole Bab el-Mandeb exchange lives above it.
Depth stays 3200 despite the Red Sea reaching only 2478 m, because the deep part is beside the
boundary; only 0.7% of columns exceed 2500 m, so those levels are cheap to leave coarse (GLORYS
salinity is 34.776–34.802 over 2871–3226 m — homogeneous).

**z\* is on by default**, because of the 0.94 m surface cell rather than tides. See the
vertical-coordinate section of `../docs/salish_sea_notes.md`: a static grid is a *linear* free
surface that never fails but carries an O(max|η| / Δz_top) error, and here that ratio is O(1).
Evaporation being this basin's driver sharpens it — a freshwater flux over a fixed-volume thin cell
cannot concentrate salt in a thinning layer. `REDSEA_ZSTAR=0` falls back to static z, in which case
raise `z_scale` to ~600 for a 2.8 m cell.

## Differences from the Salish Sea script

- **One open boundary, on the east.** The Flather sign flips: `U = Uᵉˣᵗ + √(gH)(ηᵇ − ηᵉˣᵗ)`, and
  inflow is `−u`, so every transport sign is inverted.
- **`interpolation_passes = 5`, not 10**, with a startup check that reports the Bab el-Mandeb
  cross-section and warns if smoothing closed it. Ten passes shut a comparable Salish channel.
- **Sponge 0.5° on the east edge only**, ~7° from the strait, so the passage runs free.
- **`net_transport_target` sign is inverted in principle**: evaporation-driven, so once surface
  forcing is added the target becomes net *inflow*, where a river-forced domain needs net outflow.
- **Exchange diagnostic is at the strait, not the boundary.** At 50°E nothing has been transformed
  yet — the boundary ΔS came out at 0.03 psu. The Bab el-Mandeb section is the number to compare
  with observations (~0.3–0.4 Sv each way).

## Status and what is untested

Smoke tests at 0.1–0.2° run clean end to end: BCs, transport correction, diagnostics, output writers,
and z\* all engage without error. The IC is faithful (`T ∈ [2.012, 28.203]`, `S ∈ [34.779, 41.869]`,
inside the parent range) and the strait diagnostic gives **ΔS ≈ −0.3**, correctly negative for an
inverse estuary — fresher in at the surface, saltier out at depth.

Not yet established:

- **No validation run.** Coarse smoke tests only. Strait transports are ~20× below literature, but at
  11 km the strait is 8 badly-resolved cells and the runs are 0.02 days from rest, so that number
  means nothing yet.
- **z\* is unexercised on the upgraded packages** (Oceananigans 0.110.17 / NumericalEarth 0.6.1) and
  showed a small unexplained salinity drift in the Salish hour-test against the older version.
- **No surface forcing.** The Red Sea is evaporation-driven; without it the interior will not develop
  a realistic overturning. Fine for a boundary test, not for physics.
- **GLORYS validity was checked in the boundary box only.** The Red Sea interior deserves the same
  check — in the Salish Sea it was the interior nudging, not the boundary, that did the damage.
- **Seasonal timescales.** Meaningful Red Sea circulation runs are months; days to weeks suffices for
  boundary testing.

![candidate domain](redsea_candidate_domain.png)
