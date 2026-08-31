# Canary upwelling configuration (NW Africa)

Third OBC test domain. **Status: written, geometry and parent data verified, smoke-tested end to end
at 0.1°. No validation run yet.** Model: `canary_glorys_obc.jl`.

## Why this domain

Each domain so far has tested one thing the previous could not.

| | Salish Sea | Red Sea | **Canary** |
|---|---|---|---|
| open boundaries | 2, in 1–2 parent cells | 1, well resolved | **3, all in deep open ocean** |
| unresolved feature | none — parent could not resolve the channels either | Bab el-Mandeb, inside | **island wakes + upwelling filaments** |
| interior driver | rivers (outside the domain) | evaporation | **wind** |

This is the first domain where **the atmosphere is the point rather than a correction**. It was
considered and rejected when the Red Sea was chosen — "cleanest boundaries but needs wind forcing this
config has none of" — and became viable only once JRA55 was wired and validated in `red_sea/`.

## Verified geometry

Measured from ETOPO via `regrid_bathymetry`, not assumed.

| | |
|---|---|
| Domain | 21.5–10.0°W, 25.0–29.5°N |
| Open boundaries | **west, south, north** — the east edge at 10°W is entirely land (African coast) |
| West (−21.5) | 100% wet, 1 segment, 4578–4885 m, **2388 km²** |
| North (29.5) | 99% wet, 1 segment, 30–4885 m, **3504 km²** |
| South (25.0) | 58% wet, 1 segment, 14–4584 m, **1712 km²** |
| Wet fraction | **72.7%** (Red Sea 16.1%, Salish Sea 6.4%) |
| Max depth | 4885 m at −21.48°E, 29.48°N — the NW boundary corner |
| Median wet depth | 3504 m |

Combined open-boundary area is **7.6e3 km²**, ten times the Red Sea's 719 km².

All six major islands resolve, even at the 0.1° smoke resolution (La Palma 6/16 land cells,
El Hierro 2/16, Tenerife 12/16, Gran Canaria 14/16, Fuerteventura 10/16, Lanzarote 6/16). The startup
check reports these counts — over-smoothing dissolves the wake-generating topography silently, so
check before trusting a run.

## The depth constraint — opposite in sign to the Red Sea

GLORYS's **deepest valid level over this region is −4836.6 m**, and only 2.9% of columns have data
there. Raw ETOPO reaches 4885 m, so **the model floor would sit below the parent** — the Salish Sea
failure condition. The Red Sea had +50 m of headroom; this has −48 m.

It is cheap to fix because the deep water is rare:

| deeper than | columns | |
|---|---|---|
| 4700 m | 803 | 5.3% |
| **4800 m** | **181** | **1.2%** |
| 4836 m | 48 | 0.3% |

**`depth = 4800` is a hard cap, not a convenience.** It clips 1.2% of columns by at most 85 m of
abyssal plain and keeps every model cell inside the parent's valid data. Do not raise it.

## Parent data — `glorys_interior_check.jl`

Adapted from `red_sea/`. Numbers at **0.02°, 1993-01-01**: 93,973 wet columns, 2,626,186 wet cells.

**Coverage.** 88.65% of wet columns sit below the parent's deepest valid *centre*, but only **1.68% of
wet cells** do — the deficit is a thin sliver at the bottom of many columns, not a systematic hole.
Where it matters most it barely exists:

| boundary | wet columns | cells below parent |
|---|---|---|
| west | 225 | **0 / 7200 (0.00%)** |
| south | 332 | 73 / 8663 (0.84%) |
| north | 569 | 178 / 16479 (1.08%) |

**Provenance.** No non-physical values. Water masses are a textbook North Atlantic column:

| band | S (median) | T (median) |
|---|---|---|
| surface / Ekman | 36.933 | 20.55 |
| upwelling source (100–300 m) | 36.489 | 17.57 |
| NACW/AAIW (300–1000 m) | 35.604 | 11.12 |
| NADW (1000–3000 m) | 35.148 | 4.67 |
| abyssal | 34.918 | 2.21 |

**Stability — much cleaner than the Red Sea.** 40,321 unstable pairs (1.59%), but the tail that
matters is **5 pairs at Δρ ≥ 0.1, worst 0.133 kg/m³**, against 700 pairs at worst 0.92 for the Red Sea.
All five sit in the top 12 m at La Palma and the shelf, in columns where GLORYS has land and the model
has water. This domain has no hypersaline coastal water for horizontal inpainting to import, which was
the Red Sea's actual failure mechanism.

![interior check](canary_glorys_interior_check.png)

## The GLORYS residual is 35× the Red Sea's — watch this

There is **no net-transport correction**, following the measurement in `red_sea/` that a Flather
boundary absorbs the parent's spurious net transport unaided. But that was measured at 1.44e5 m³/s,
and here the residual is far larger:

| | Canary | Red Sea |
|---|---|---|
| net residual at t=0 | **−5.009e6 m³/s** (−5 Sv) | −1.44e5 |
| per boundary | west −6.02e6, south +1.72e6, north −7.10e5 | — |
| uncompensated sea level | **−1.05 m/day** | −0.02 m/day |
| η offset Flather needs | **1.3 cm** | 0.4 cm |

The west edge dominates, which is the Canary Current — a few Sv through a 2388 km² section is
physically ordinary; it is the *net* that has to close. The scaling still says the boundary absorbs
it, because the compensating velocity is only 6.6e-4 m/s over an area ten times larger. **But this
extrapolates the no-correction conclusion 35× beyond where it was tested**, and uncorrected it would
be a metre a day. `progress` reports actual transport per boundary alongside the uncorrected residual
every block, so a divergence is visible immediately. This is the sharpest available test of that
decision.

## Resolution and cost

| Δ | km at 27.25°N | cells | 7-day run |
|---|---|---|---|
| 0.1° | 9.9 | 0.17 M | smoke tests only |
| 0.05° | 5.0 | 0.66 M | ~2 h |
| **0.02°** | **2.06** | **4.1 M** | **~21 h** (extrapolated from `red_sea/`) |
| 0.01° | 1.03 | 16.6 M | ~7 days — GPU/HPC only |

The domain is 7× smaller in area than the Red Sea, which is why 2 km is affordable here when it was
not there. **0.02° is the CPU validation resolution**; 1 km is the target for island wakes and
upwelling filaments and needs a GPU.

## Vertical grid

`Nz = 32`, `depth = 4800`, `z_scale = 700` → **1.21 m surface cell, 11 levels in the top 50 m, 16 in
the top 150 m**, at the cost of ~930 m bottom cells.

That trade is deliberate: upwelling lives in the Ekman layer and the top few hundred metres, and the
source water for coastal upwelling here is around 100 m. The abyss is quiescent and can carry coarse
cells. 32 levels over 4800 m is 50% more depth than the Red Sea on the same budget, so the deep
resolution is genuinely thin — a limitation for anything involving deep water masses.

z\* is on by default, as in `red_sea/`: with a 1.2 m surface cell the O(max|η|/Δz_top) error of a
static (linear free surface) grid is not small. `CANARY_ZSTAR=0` falls back.

## Differences from the other scripts

- **Three open boundaries.** `salish_sea/` is the template (west + north); south is new. Inflow is
  `+u` at west, `+v` at south, `−v` at north. Oceananigans applies the side-dependent Flather sign
  itself, so the closures return the exterior value only — **do not pre-flip signs**, which is exactly
  the double inversion `red_sea/` carried through its whole validation run.
- **Sponge on all three edges**, combined with `max()` so corners do not double-count. 0.25°
  (`CANARY_SPONGE`) — the upwelling front sits against the eastern coast, far from every sponge.
- **An upwelling index in the progress output**: coastal minus offshore SST at matched latitudes, plus
  mean temperature at the ~100 m source depth. This is the signal the domain exists to produce, and
  the other two configurations needed no equivalent. At initialisation it reads **ΔT = −2.25 °C** —
  the upwelling is already in GLORYS; the question is whether the model holds or sharpens it.
- **No net-transport correction** (see above).

## Env flags

| flag | default | |
|---|---|---|
| `CANARY_DELTA` | `0.02` | horizontal resolution, degrees |
| `CANARY_START` | `1993-01-01` | run start |
| `CANARY_DEPTH` | `4800` | **hard cap** — see the depth section |
| `CANARY_ZSCALE` | `700` | vertical stretching |
| `CANARY_SPONGE` | `0.25` | nudging band width, degrees |
| `CANARY_ATMOS` | `JRA` | `JRA` or `none` (boundary-only control) |
| `CANARY_JRA_DATASET` | `repeat` | `repeat` or `multiyear` |
| `CANARY_ZSTAR` / `CANARY_PASSES` / `CANARY_WENO_ORDER` / `CANARY_STOP_DAYS` / `CANARY_TAG` | | as in `red_sea/` |

## Status and what is untested

Smoke-tested end to end at 0.1°: all three boundaries carry transport, z\*, JRA55, the conservation
and stability monitors, the upwelling index and the output writers all engage. Drift ~1e-6 over the
test. Not yet established:

- **No validation run at any resolution.** Everything above is a 0.1° mechanical smoke test.
- **The 5 Sv GLORYS residual is untested against the no-correction decision** (see above). This is the
  first thing a real run should be judged on.
- **January is the weak-upwelling season.** The window matches `salish_sea/` and `red_sea/`, but
  Canary upwelling peaks in boreal summer at these latitudes. `CANARY_START=1993-07-01` gives the
  strong case, and the domain's own signal will be clearer there.
- **The south boundary is only 58% wet** — the African coast cuts through its eastern end, where a
  sponge meets a coastal boundary. Structurally the least clean part of the domain, and the most
  likely place to need the boundary moved rather than the sponge tuned.
- **Deep resolution is thin**: ~930 m bottom cells. Fine for upwelling, not for deep water masses.
- **1 km is untested** and is where island wakes and filaments actually live.
