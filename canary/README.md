# Canary upwelling configuration (NW Africa)

Third OBC test domain. **Status: validated — stable 7-day JRA55-forced run at 0.02° (2 km),
38.8 h wall. January only; July not yet run.** Model: `canary_glorys_obc.jl`.

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

## The GLORYS residual, and the no-correction decision

There is **no net-transport correction**, following the measurement in `red_sea/` that a Flather
boundary absorbs the parent's spurious net transport unaided. That was measured at 1.44e5 m³/s; here
the residual is ~18× larger, so this domain is the sharpest available test of that decision.

| | Canary (0.02°) | Red Sea |
|---|---|---|
| net residual at t = 0 | **−2.033e6 m³/s** | −1.44e5 |
| per boundary at t = 0 | west −5.31e6, south +2.10e6, north +1.18e6 | — |
| typical net once running | −2.6e6 | — |
| uncompensated sea level | **−0.55 m/day** | −0.03 m/day |
| η offset Flather needs | 0.69 cm | 0.04 cm |

**The decision holds.** At t = 2 days the residual reads **−2.605e6 m³/s** while the model's actual net
transport through the same three sections is **+1.345e5 m³/s** — the boundary absorbs ~95% of a 2.6 Sv
imbalance without help, and volume drift sits at −3.0e-6 with mean S flat. The Red Sea conclusion
extrapolates to 18× its measured magnitude.

The west edge dominates the residual, which is the Canary Current: a few Sv through a 2388 km² section
is physically ordinary, and it is the *net* that has to close.

> **Earlier versions of this file claimed −5.009e6 m³/s and "35× the Red Sea".** That was wrong on two
> counts: it was measured at 0.1° rather than the run resolution, and it included an off-by-one in the
> north-boundary diagnostic (sampling y-face `Ny`, an interior face, instead of `Ny+1`, the domain
> edge). Correcting the index moved the north residual by 3.5 Sv *and flipped its sign*. The boundary
> *conditions* were never affected — `north_V_obc` and `v_north_exterior` always used `Ny+1` — only
> the diagnostics were.

## Layered transports: comparing to the Canary Current literature

Full-depth net transport is the volume budget, not a quantity any paper reports. The inverse box
models (Hernández-Guerra 2005; Pérez-Hernández 2013; Casanova-Masjoan 2020) integrate over **neutral
density (γⁿ)** layers — thermocline γⁿ < 27.38, intermediate 27.38–27.922, deep > 27.922.

γⁿ needs the Jackett & McDougall lookup, which this stack does not have, so **σθ = ρ(T,S,0) − 1000
stands in**. Two consequences, both measured rather than assumed:

- **The deep threshold had to move.** σθ in this domain maxes at **27.798**, so a γⁿ = 27.922 bound is
  unreachable and the deep bin reads zero. It is set to **27.70**, which is NADW here.
- **The layers sit too deep.** The run prints where each σθ surface actually lies:

  | surface | our σθ | literature γⁿ |
  |---|---|---|
  | 26.85 | **406 m** (median) | ~80–320 m |
  | 27.38 | **1195 m** | ~700 m |

  So our "thermocline" swallows the AAIW the literature counts as intermediate and **reads high**
  against their numbers. Quote it with that caveat.

**Fixed-depth bins would be worse, not simpler.** Isopycnals heave seasonally, so a fixed split aliases
that heave into apparent transport change — corrupting exactly the January-versus-July comparison the
diagnostic exists to make.

Early result at t = 2 days, thermocline layer: **north +5.7 Sv (inflow), south +2.0, west −3.2**.
The north inflow is the Canary Current entering the domain and compares well with Hernández-Guerra's
**4.7 ± 0.8 Sv southward** north of the islands, allowing for the deep-layer bias. The **south reading
+2.0 Sv inflow is not yet explained** — the Canary Current picture says the flow should leave through
the south — and it has been persistently positive, so it is not noise.

All three boundaries are also **saved as sections** (plus a separate `north_face_v` at `Ny+1`, since
face-in-y indexing cannot put `v` and `T`/`S` in one slice), so the layer definition can be re-binned
offline — including with a proper γⁿ implementation later — without re-running.

## Resolution and cost

| Δ | km at 27.25°N | cells | 7-day run |
|---|---|---|---|
| 0.1° | 9.9 | 0.17 M | smoke tests only |
| 0.05° | 5.0 | 0.66 M | ~2 h |
| **0.02°** | **2.06** | **4.1 M** | **~73 h measured** (not the ~21 h extrapolated) |
| 0.01° | 1.03 | 16.6 M | ~7 days — GPU/HPC only |

The domain is 7× smaller in area than the Red Sea, which is why 2 km is affordable here when it was
not there. **0.02° is the CPU validation resolution**; 1 km is the target for island wakes and
upwelling filaments and needs a GPU.

**The 21 h figure was wrong and the reason is worth recording.** It scaled Δt from the Red Sea by grid
spacing, but the binding constraint here is **vertical, not horizontal**: with a 1.4 m surface cell and
wind-driven convection, Δt settles at 1.4–3.3 min where the horizontal CFL alone (max|u| ≈ 0.5 m/s over
2 km) would allow ~7 min. At ~40 s/iter that is **~73 h for 7 days**. Raising `CANARY_ZSCALE` would buy
timestep at the cost of the surface resolution that motivated `scale = 700`.

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

## Validation — 7-day run at 0.02°, JRA55-forced, January

3620 iterations, **38.8 h wall** at ~40 s/iter, Δt settling at 1.4–4.5 min.

### The no-correction decision holds at 21× the magnitude it was measured

This is what the domain was chosen to test. Mean over t ≥ 4 d:

| | |
|---|---|
| GLORYS residual through the three sections | **−3.033e6 m³/s** |
| model's actual net transport | **+2.5e4 m³/s** |
| **absorbed by the Flather boundary** | **99.2%** |
| volume drift | −1.2e-6 |
| salt drift | +2.7e-6 |
| mean S | 35.3068 → 35.3069 |

`red_sea/` measured the no-correction conclusion at 1.44e5 m³/s. It extrapolates to **21× that**
without a correction term and without drift.

### The Canary Current is there — but the first diagnostic hid it

The on-the-fly thermocline transport drifted from **+8.5 Sv (north inflow) at day 0.5 to −0.5 Sv
(outflow) at day 7**, which read as the model losing the Canary Current. It was not. Re-binning the
saved boundary sections offline into finer density bands (mean over t ≥ 4 d, Sv, + = into the domain):

| band | mean depth | west | south | north |
|---|---|---|---|---|
| σθ<26.5 | ~70–110 m | −0.81 | +0.21 | −0.09 |
| 26.5–26.85 | ~210–330 m | −1.01 | +0.51 | **+1.28** |
| 26.85–27.1 | ~420–520 m | −1.11 | +0.52 | **+4.58** |
| 27.1–27.38 | ~780–840 m | −2.10 | **+4.38** | **−4.78** |
| 27.38–27.70 | ~1250–1470 m | +0.28 | **+5.27** | **−2.60** |
| σθ>27.70 | ~2800–3200 m | −1.29 | −2.78 | −0.41 |

Two opposing currents sat inside the single σθ<27.38 bin:

- **Above ~500 m: north inflow totalling +5.8 Sv** — the Canary Current, southward through the
  domain. Compares well with Hernández-Guerra's **4.7 ± 0.8 Sv** north of the islands, allowing for
  the layer sitting deeper than theirs.
- **700–1500 m: north outflow, south inflow, ~4–5 Sv** — a **poleward undercurrent**, flowing the
  opposite way. A real feature of this eastern boundary (*The Canary Intermediate Poleward
  Undercurrent*, JPO 2021), and visible in the cross-shore section as a northward band hugging the
  African slope.

As the undercurrent spun up over the run it overwhelmed the surface signal and flipped the sum's sign.
**The bin has since been retuned to σθ < 27.10 (~500 m)** — matching the *depth* the literature's
γⁿ < 27.38 layer occupies rather than its numerical label, since σθ ≠ γⁿ.

This is what the saved boundary sections bought: the diagnosis needed re-binning, and re-binning from
file cost minutes where re-running would have cost 39 hours.

### Upwelling: present, unchanged, and weakly tested

`ΔT` (coastal minus offshore SST) went **−2.260 → −2.289 °C**, mean −2.268 over t ≥ 4 d. The model
neither sharpened nor destroyed the gradient GLORYS supplied. Temperature at the 109 m source depth
over the shelf held at 17.53 °C. **January is the weak season here** — this is a null result on the
domain's headline signal, and the July run is the real test.

![cross-shore section](canary_section_day7.png)

### Open after this run

- **The undercurrent at 4–5 Sv looks large** for an eastern boundary undercurrent (typically O(1 Sv)).
  No literature value was found to check it against.
- **~4 Sv of apparent transfer from upper to deeper layers** in the mean, too large to be real
  diapycnal flux. Most likely isopycnal heave changing layer volumes, since the binning is
  instantaneous — not verified.
- **The thermocline circulation was still drifting at day 7**, so 7 days does not settle this domain's
  transports any more than it settled the Red Sea's strait.

## Status and what is untested

Smoke-tested end to end at 0.1°: all three boundaries carry transport, z\*, JRA55, the conservation
and stability monitors, the upwelling index and the output writers all engage. Drift ~1e-6 over the
test. Not yet established:

- ~~No validation run.~~ **Done** — 7 days at 0.02°, stable, 38.8 h.
- ~~The GLORYS residual is untested against the no-correction decision.~~ **Done** — 99.2% absorbed.
- **July is not run.** January is the weak-upwelling season and the ΔT result is null; the domain's
  headline signal is untested.
- **January is the weak-upwelling season.** The window matches `salish_sea/` and `red_sea/`, but
  Canary upwelling peaks in boreal summer at these latitudes. `CANARY_START=1993-07-01` gives the
  strong case, and the domain's own signal will be clearer there.
- **The south boundary is only 58% wet** — the African coast cuts through its eastern end, where a
  sponge meets a coastal boundary. Structurally the least clean part of the domain, and the most
  likely place to need the boundary moved rather than the sponge tuned.
- **Deep resolution is thin**: ~930 m bottom cells. Fine for upwelling, not for deep water masses.
- **1 km is untested** and is where island wakes and filaments actually live.
