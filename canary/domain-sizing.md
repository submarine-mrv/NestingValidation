# How far south must the Canary domain reach?

Status as of 2026-09-25. Model: `canary/canary_glorys_obc.jl`. Analysis: `canary/virtual_boundary.jl`.

The Canary configuration exists to answer a residence-time question: water released at a point on
Gran Canaria (−15.38°E, 28.16°N) — a CDR-style continuous surface source — drifts southwest with the
Canary Current, and we want to know **how much of it the domain still holds after months**, and
therefore **how big the domain has to be** for that number to mean anything.

The current domain (φ ∈ [25.0, 29.5]) puts its southern open boundary **352 km** from the release.
A 90-day run at 0.02° ends with **88.8% of released tracer still inside**, and the leak is almost
entirely southward: at day 90 the boundary export splits **south −9.8e−2, west −1.1e−2, north
−1.1e−2 /s**, so south carries ~81% of it. The domain is losing tracer through exactly the boundary
that a larger domain would move.

## The method: virtual boundaries

A run only holds data inside its own boundaries, so a *larger* domain cannot be measured from a
smaller one's output. What can be measured is the family of **nested sub-domains**: pick a virtual
southern boundary, integrate the tracer inventory north of it, and repeat for several latitudes.
That gives retention as a function of distance-to-boundary `L` over a measured range, from one run.

`canary/virtual_boundary.jl` does this, writing one CSV row per 3-D snapshot frame. Virtual
boundaries are kept at or north of 25.75°N so none sits inside the real sponge (a 0.25° Gaussian at
25.0, felt to ~0.75°).

**On the bias — the obvious argument is wrong, and a real run settled it.** Tracer crossing a virtual
boundary southward can drift back north and be counted as retained again, which suggests these are
upper bounds. This note asserted exactly that for some time. Measurement says otherwise: the effect
is real but small, and it is outweighed by one the argument omits — a *real* boundary further south
permanently flushes tracer that, in a larger domain, would have returned. Net, the method
**under**-estimates retention by up to ~4 points. See *The check: a real φ_S = 22.0 run* below, which
compares every number in this section against a domain that actually extends to 22°N.

## What the 90-day run measured

90 days at 0.02° (575 × 225 × 32), GPU, JRA55-forced, 13.6 h wall. Release normalised to
∫S dV = 1 unit/s, so tracer released by time *t* is exactly *t* in seconds — which makes every
retention figure below a clean fraction with no calibration step.

Retention north of each virtual boundary, as a percentage of **released** (not of what is still in
the domain):

| φ_b | L (km) | day 30 | day 60 | day 90 | τ (days) |
|---|---|---|---|---|---|
| 27.50 | 73 | 74.3 | 67.3 | 57.5 | 68 ± 2 |
| 27.25 | 101 | 89.7 | 75.0 | 60.1 | 87 ± 2 |
| 27.00 | 129 | 94.6 | 82.6 | 63.9 | 114 ± 4 |
| 26.50 | 185 | 97.5 | 92.0 | 75.6 | 201 ± 10 |
| 26.00 | 241 | 97.5 | 92.6 | 84.6 | 298 ± 9 |
| 25.75 | 268 | 97.5 | 93.1 | 87.0 | 332 ± 8 |
| *25.00 (real)* | *352* | — | — | *88.8* | — |

τ is fitted to a well-mixed flush model, `N(t) = Sτ(1 − e^(−t/τ))` with S = 1 unit/s. The model is
wrong at early times — the plume needs 16–40 days just to reach these latitudes, so nothing leaks
before then — but the fit is insensitive to that: refitting over t ≥ 30, 45 and 60 days moves every
τ by less than 5%. The quoted ± is the fit's formal 1σ and is much smaller than the structural
uncertainty below.

## τ(L), and how far it can honestly be pushed

τ rises steeply and near-linearly with L over the measured range:

    τ ≈ 1.26 (±0.14) d/km × L − 33 (±14) d          r = 0.994

A linear τ(L) is what **advective** escape looks like (τ = L/U, constant U); diffusive escape would
give τ ∝ L². The implied bulk escape speed is U ≈ 0.8 km/d ≈ 0.9 cm/s — much slower than the
first-arrival speed, which is 7.4 km/d (8.6 cm/s) measured from when each boundary first passes 5%
escape. Those are consistent: a fast filament reaches the boundary early and a slow bulk drift sets
the flush rate.

**The linear form is not established, only preferred-ish.** A quadratic fits the same six points
better (rms residual 6.7 d against 17.5 d). With six correlated points there is no basis to choose,
and the two forms diverge fast outside the measured range. Extrapolation beyond L ≈ 270 km should be
read as order-of-magnitude:

| φ_S | L (km) | τ from linear fit |
|---|---|---|
| 25.0 (current) | 352 | ~410 d |
| 24.0 | 463 | ~550 d |
| 23.0 | 574 | ~690 d |
| 22.0 | 686 | ~830 d |

Even taking these as upper bounds and halving them, the current domain already has τ ≫ the 90-day
run length. **The domain is not the limiting factor on a 90-day experiment; it would be on a
multi-year one.**

## Tracer advection: bound-preserving WENO is not usable here

Plain WENO undershoots into c < 0 in the dilute plume edge — 6.8% of surface cells at day 90, and
enough to make the pointwise age `ac/c` meaningless there (ages of 29000 days in a 90-day run), which
is why `animate_tracer.jl` masks age below 1% of the instantaneous max. The obvious fix is WENO's
`bounds` argument, and `tracer_bounds: ["true"]` wires it up. **It must stay off.**

Two 0.1° runs differing only in the tracer scheme (`_d01_bwenotest` vs `_d01_plainweno`, same grid,
same times, identical T ranges) measure the release against a source normalised to exactly 1 unit/s:

| scheme | d(inventory)/dt, day 1 | day 2 | min c |
|---|---|---|---|
| plain WENO | 0.9996 | 0.9993 | −2.7e−6 |
| bound-preserving | **1.2357** | **1.2514** | −1.1e−6 |

Bound-preserving WENO **creates tracer at ~25%/day**, and does not even eliminate the negatives it
was meant to fix. The cause is visible in Oceananigans 0.110.17: `BoundsPreservingWENO` defines a
single `div_Uc` (`src/Advection/bounds_preserving_tracer_advection_operators.jl`) that passes `ρ = 1`
and divides by `Vᶜᶜᶜ`, with no specialisation for a mutable vertical coordinate — and the file's own
comment hedges on whether it is immersed-boundary safe. Under z\* with `GridFittedBottom`, both apply.

This matters more than a tracer-scheme detail: inventory is the quantity this entire study measures,
so enabling it would have inflated every retention number by ~25% while looking perfectly stable.

## The finding that actually constrains the domain

Extending south stops helping at **φ_S ≈ 22.8°N**, and not for any oceanographic reason.

First, a correction to an earlier version of this note, because it changes how the number should be
read. Distances from the release (−15.38°E, 28.16°N) to each open boundary:

| boundary | distance | share of export, 25°N run at day 90 |
|---|---|---|
| **north** (29.5°N) | **149 km** | 10% |
| west (−21.5°E) | 600 km | 9% |
| south (25.0°N) | 352 km | **81%** |

**North is the nearest open boundary, and always has been** — in the 22°N domain too, since only the
southern edge moved. An earlier draft claimed west becomes "the nearest open boundary" past φ_S =
22.8; that is simply false, and the framing was wrong. Proximity does not set the leak: the Canary
Current advects the plume southwest, so the southern boundary took 81% of the export despite sitting
2.4× further away than the northern one.

What the crossover actually says is narrower, and still useful: **west and south are equidistant from
the release at φ_S = 22.76**, so beyond that the western boundary is the nearer obstacle *along the
direction the tracer actually travels*. Extending south alone then stops buying distance in the
direction that matters.

The 22°N run supports this. At day 32 its export is **north −4.3e−2 /s, south and west numerically
zero** (1e−60): with the southern boundary moved to 686 km the plume has not reached it, and what
leaves does so through the near northern edge. So moving the southern boundary works as intended —
and it promotes the *northern* boundary, not the western one, to dominant leak in the interim.

This is why the `lat_south` knob is only half a solution, though the other half is not the one the
earlier draft named: a domain-sizing study that moves the southern boundary alone leaves both the
northern boundary at 149 km and the western at 600 km fixed.

Cost, for reference — extending south is cheap in cells but not free in data:

| φ_S | Ny at 0.02° | cost vs current |
|---|---|---|
| 25.0 | 225 | 1.00× |
| 24.0 | 275 | 1.22× |
| 23.0 | 325 | 1.44× |
| 22.0 | 375 | 1.67× |

The 90-day run would go from 13.6 h to ~23 h at φ_S = 22.0. **The GLORYS staging is the real cost**:
the parent bounding box follows the domain with a 0.5° margin and its extent is part of every cached
filename, so changing φ_S invalidates the entire cache and forces a fresh download. Nothing already
staged is reusable.

## What this does and does not settle

Settled:

- The current domain holds 88.8% of a 90-day release; the southern boundary carries ~81% of the
  loss.
- Retention timescales are advective, not diffusive, over 73–268 km, with a bulk escape speed of
  ~0.9 cm/s and a first-arrival speed of ~8.6 cm/s.
- τ at the current domain size is ~400 days (upper bound), so 90-day and even 6-month experiments
  are not domain-limited.
- At 90 days only. Moving the southern boundary past ~22.8°N stops buying distance along the direction the tracer
  travels, because west and south are equidistant from the release there. Note this is *not* a
  statement about the nearest boundary: north, at 149 km, is nearest in both domains.

Not settled:

- **Every number here is an upper bound** on retention, by the re-entry bias above. The size of that
  bias is unmeasured. The way to measure it is a real run at φ_S = 22.0 compared against this run's
  virtual boundary at the same latitude — that is what `lat_south` is for, and it has not been
  run.
- **The extrapolation past L ≈ 270 km is not validated**, and linear vs quadratic τ(L) is unresolved.
  The φ_S = 22.0 run would test that too, since its virtual boundaries would extend the measured
  range to ~690 km.
- **One season, one release site, one year.** The Canary Current is strongly seasonal and this is a
  single 90-day window.
- **Reproducibility gap**: `canary/out3mo/vb.csv` was computed on the GCP VM from
  `canary_glorys_snapshot3d_zstar_jraryf_tracer3mo.jld2`, which stayed on that machine
  (`/home/vtamsitt/nesting/canary/`) and is not in the repo. Only the derived CSV came back. The run
  **started 1993-01-01** — not recorded in `tracer3mo.log`, but established by its t = 0 state being
  bit-identical to the January 7-day run. The log header should print the start date.

## The check: a real φ_S = 22.0 run, measured against the extrapolation

90 days at 0.02°, 1.29 days wall on an L4, identical to the 25°N run in every respect but the
southern boundary — same start (1993-01-01), same resolution, same forcing, same plain-WENO tracer
scheme. `canary/out3mo/vb_phis22.csv`. Domain inventory at day 90: **91.36% of released**, against
88.82% for the 25°N domain.

### The virtual-boundary method is accurate to a few percent — and biased the other way

Retention north of each shared latitude at day 90, as a percentage of released:

| φ_b | L (km) | 25°N run (virtual) | 22°N run (measured) | error |
|---|---|---|---|---|
| 25.75 | 268 | 86.96% | 89.74% | −2.78 |
| 26.00 | 241 | 84.64% | 87.73% | −3.09 |
| 26.50 | 185 | 75.59% | 79.50% | −3.91 |
| 27.00 | 129 | 63.92% | 66.12% | −2.19 |
| 27.25 | 101 | 60.09% | 60.23% | −0.14 |
| 27.50 | 74 | 57.46% | 54.21% | +3.25 |

**This note previously claimed every virtual-boundary retention number was an UPPER bound. That was
wrong.** The stated reasoning — tracer crossing a virtual boundary southward can return and be
counted as retained again — is real but small, and it is outweighed by an effect the argument missed:
in the 25°N run the *real* boundary at 25.0 permanently flushed tracer that, in the larger domain,
stays and circulates back north. Depletion beats double-counting. The method therefore slightly
UNDER-estimates retention, by ≤4 points at every latitude tested, with the sign flipping only at
27.5 where the virtual boundary is far from the real one.

Either way the error is a few percent, which is the useful finding: **virtual boundaries inside one
run predict a real larger domain to within ~4 percentage points**, so the technique is sound and the
90-day retention numbers can be trusted at that level.

### The extrapolated thresholds held up

Predicted from the 25°N τ(L) fit against measured in the 22°N run, at 90 days:

| target retention | predicted L | measured L | error |
|---|---|---|---|
| 60% | 94 km | ~101 km | 7% low |
| 70% | 121 km | ~140 km | 16% low |
| 80% | 174 km | ~187 km | 7% low |
| 90% | 331 km | ~300 km | 10% high |

Accurate to 7–16% in distance, and the direction of the error at 60–80% matches the caveat already
recorded here: the well-mixed model over-predicts retention mid-range, so its latitudes are biased
north and the true boundary sits slightly further south.

### Retention saturates, so the southward extension was not needed

| φ_S | L (km) | retention at day 90 |
|---|---|---|
| 23.0 | 574 | 91.34% |
| 24.0 | 463 | 91.26% |
| 25.0 | 352 | 90.90% |
| 25.5 | 296 | 90.35% |
| whole 22°N domain | — | 91.36% |

Going from 352 km to 574 km buys **0.44 percentage points**. Essentially all the tracer still in the
domain at day 90 is north of 23°N, and almost all of it north of 25°N. **At 90 days the 25°N domain
was already large enough**; the extension to 22°N cost 1.67× the cells and ~29 h of GPU to confirm
that the boundary was not the binding constraint. That is a real answer, and the cheaper way to have
reached it was the virtual-boundary analysis, which said so to within a few percent from a run that
already existed.

The caveat that survives: this is a 90-day statement. τ(L) at the current domain size is several
hundred days, so a multi-year experiment would see the boundary, and the saturation above would not
hold.

## What this means for the CDR question — the 90-day answer is not the answer

**"The 25°N domain was already large enough" is a statement about 90 days only, and it is easy to
misread as settling the CDR question. It does not.** Air-sea CO₂ equilibration runs months to about a
year, so the timescale that matters for whether uptake can be accounted for inside the domain is
~1 year, not 90 days. Over 90 days the plume has barely begun to spread and no domain size matters.

The φ_S = 22.0 run does sharpen the answer, by extending the measured range of τ(L) from 73–268 km to
73–574 km. Refitted on its 13 virtual boundaries, every point now constrained to better than 6%:

    τ = 0.854 d/km × L + 106 d

| φ_S | L (km) | τ (d) | 90 d | 180 d | **1 yr** | 2 yr |
|---|---|---|---|---|---|---|
| **25.0 (current)** | 352 | 407 | 89.7% | 80.8% | **66.0%** | 46.4% |
| 24.0 | 463 | 502 | 91.5% | 84.0% | 71.0% | 52.7% |
| 23.0 | 574 | 597 | 92.8% | 86.3% | 74.8% | 57.7% |
| 22.0 | 686 | 692 | 93.8% | 88.0% | 77.7% | 61.8% |
| 21.0 | 797 | 787 | 94.5% | 89.4% | 80.0% | 65.2% |

The 90-day column reproduces the measurements (89.7% predicted against 88.8% measured at 25°N), which
is the fit's own consistency check.

**At one year the current domain holds ~66%, so about a third of the treated water is outside it; by
two years more than half is.** That is the concern this note opened with, and it survives.

Southern boundary needed for a target retention *at one year*:

| target | boundary |
|---|---|
| 50% | 26.9°N |
| 60% | 25.9°N |
| 70% | 24.2°N |
| 80% | **21.0°N** |

**This is more pessimistic than the earlier extrapolation, not less.** The old fit put 80%-at-one-year
at φ_S ≈ 22.9°N; the better-constrained one puts it at 21.0°N, because the old fit over-estimated τ at
long range by roughly 20%.

**The weak link has moved from space to time.** Extrapolating across domain size is now validated to
7–16%. Extrapolating from a 90-day record to one- and two-year retention is not: those τ values are
inferred from curvature in a 90-day series, and the well-mixed flush model is carrying the
extrapolation. Treat the one-year row as the best available estimate, not a measurement. Settling it
needs a run of order a year — at the 22°N domain's ~29 h per 90 days, about 5 days of GPU.

## The coastal release is a different experiment, and its retention is an upper bound

The 90-day reference run (`_ref90`, Oceananigans 0.113.3) releases from the **south coast of Gran
Canaria** at 27.8078°N, 15.3971°W in 60 m of water. Every earlier number in this note comes from the
**north tip** at 28.16°N in 535 m — exposed, in the Canary Current. The two are not comparable, and
the difference is large: at day 56 the coastal run held **99.3%** of what it had released, against
~94% for the offshore release at the same age, with a peak concentration half again as high.

The animation makes the reason plain. The coastal plume stays a compact patch against the coast; the
offshore one spread across the domain. Boundary export follows: the coastal run exports almost
entirely through the NORTH boundary at ~2e-2 /s, with south six orders of magnitude smaller, where
the offshore run was south-dominated.

**Do not read this as a reservoir result.** Five things changed at once — Oceananigans 0.110.17 →
0.113.3, tangential velocity BCs added, the release site moved, oblique radiation on the velocities,
and TracerReservoir on the tracers. Geography is the most plausible driver, since a release in 60 m
of sheltered water on the island's lee ventilates far more slowly than one off the exposed tip, and
barely any tracer reaches a boundary for the reservoir to act on.

**And treat 99.3% as an upper bound.** At 0.02° the model does not resolve the processes that
disperse a nearshore plume: submesoscale fronts and filaments live at O(1 km), the bottom boundary
layer in 60 m of water spans only a few of the 12 levels there, and the configuration has no tides,
no surface waves and no river fluxes. Every one of those omissions removes mixing, so retention is
biased HIGH and the coastal-versus-offshore gap is likely narrower than measured here. The direction
is certain even though the magnitude is not.

Settling it needs a controlled comparison — the same configuration at both release sites, varying
only the release position — and ideally a resolution check. 0.01° fits the current hardware: the
0.02° run uses 3.4 GB of the L4's 23 GB, so 4x the cells is ~14 GB, and the cost is wall-clock
(roughly 8x) rather than capability.

## The coastal plume recirculates around the islands, which may matter more than retention

Rescaling the tracer panel so the dilute plume is visible (saturating at 1% of the peak — the peak is
7.1e-4 and its own 99.995th percentile is 7.1e-5, so nearly all the mass sits in a few cells at the
source) changes the reading of the coastal run.

**The plume travels north-west, wrapping around Gran Canaria and Tenerife — not south-west with the
Canary Current.** That is consistent with boundary export being almost entirely northward, where the
offshore release was south-dominated. The island wake appears to recirculate the release rather than
the mean flow carrying it downstream.

If that holds, the implication runs past domain sizing. Air-sea CO₂ equilibration takes months, so
tracer held in an island recirculation for weeks may be **substantially equilibrated before it ever
reaches the Canary Current or the gyre**. Where uptake happens, and how OAE efficiency should be
counted, would then look quite different from a picture in which treated water is advected promptly
downstream.

That makes the nearshore resolution question load-bearing rather than a refinement. A recirculation
set by island wakes at 2 km is exactly the structure most likely to change at higher resolution, and
it now bears on the uptake calculation and not merely on how much tracer stays inside the domain.
The 1/60° comparison is the cheap test; `w` at ~109 m, which shows filamentary submesoscale motion
strongest in the wakes and at the shelf break, is the field to compare.
