# How far south must the Canary domain reach?

Status as of 2026-09-21. Model: `canary/canary_glorys_obc.jl`. Analysis: `canary/virtual_boundary.jl`.

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

**The bias runs one way, and it matters.** Tracer that crosses a virtual boundary southward can
drift back north, and this method counts it as retained again. A *real* boundary there would have
flushed it, because the prescribed exterior carries no tracer and nothing comes back in. So every
retention number below is an **upper bound**, and every timescale derived from it is an upper bound
too. A real domain with its southern boundary at 27.0°N would retain *less* than the 63.9% quoted
for the virtual boundary there.

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
`bounds` argument, and `CANARY_TRACER_BOUNDS=1` wires it up. **It must stay off.**

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

This is why the `CANARY_PHI_S` knob is only half a solution, though the other half is not the one the
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
- Moving the southern boundary past ~22.8°N stops buying distance along the direction the tracer
  travels, because west and south are equidistant from the release there. Note this is *not* a
  statement about the nearest boundary: north, at 149 km, is nearest in both domains.

Not settled:

- **Every number here is an upper bound** on retention, by the re-entry bias above. The size of that
  bias is unmeasured. The way to measure it is a real run at φ_S = 22.0 compared against this run's
  virtual boundary at the same latitude — that is what `CANARY_PHI_S` is for, and it has not been
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
