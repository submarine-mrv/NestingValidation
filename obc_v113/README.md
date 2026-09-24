# Testing Oceananigans v0.113 oblique radiation

Isolated environment for evaluating `ObliqueRadiation` (CliMA/Oceananigans.jl#5962) before deciding
whether to adopt it in `canary/`. **Its `Project.toml` is deliberately separate**: the repo's main
environment stays on Oceananigans 0.110.17, which is what analyses the runs currently in flight.

    julia --project=obc_v113 -e 'using Pkg; Pkg.instantiate()'
    julia --project=obc_v113 obc_v113/oblique_tracer_outflow.jl

## The upgrade is two jobs, not one

`ObliqueRadiation` is pure Oceananigans and can be tested with no NumericalEarth at all. Moving
`canary/` to v0.113 is a different matter:

| | requires | |
|---|---|---|
| test `ObliqueRadiation` | Oceananigans alone | done, below |
| move `canary/` to v0.113 | **NumericalEarth 0.6.1 → 0.8.0** | not started |

The registry pins it: NumericalEarth 0.8 is the first version admitting Oceananigans 0.111-0.113,
and 0.6.1 (ours) caps at 0.111. That crosses two minor versions of the package supplying
`DatasetRestoring`, `Metadata`, `interpolate_to_grid`, `ocean_simulation` and `EarthSystemModel` --
a migration, not a bump. Keep it off the critical path.

## What the swap would actually touch

`canary_glorys_obc.jl` builds one `radiation` object and feeds it to both `NormalFlowBoundaryCondition`
(normal velocity) and `ValueBoundaryCondition` (T, S, and the release tracers). The keyword signatures
are identical, so those call sites change by one word:

    ObliqueRadiation(; inflow_timescale = 0, outflow_timescale = Inf, use_boundary_velocity = false)

Both schemes carry state (`φᵇ, φ₁, φ₁ˡ`, plus `previous_boundary`/`previous_interior` for oblique) but
every field is `Nothing` at construction and materialised per boundary, so sharing one prototype
across u/v, T/S and three boundaries is safe -- as `NormalRadiation` already does.

**But the PR's headline gain is not available to us as a one-word change.** Its enstrophy result
(0.923 against 1.611 for Dirichlet) comes from radiating the TANGENTIAL velocity, and this
configuration sets no open boundary condition on the tangential velocity at all: `u_bcs` is west only,
`v_bcs` south and north, i.e. the boundary-normal component in each case. Capturing that part means
ADDING boundary conditions, not swapping a scheme.

## Measured: tracer leaving an open boundary at an angle

`oblique_tracer_outflow.jl`. Prescribed uniform flow at angle θ to the boundary normal, so truth is
analytic (rigid translation of a Gaussian, periodic in the tangential direction) and every departure
is a boundary artefact. `rms` is the error against analytic while the blob straddles the boundary;
`resid` is the tracer mass left behind after it should have fully exited (ideal 0).

| θ | resid Dirichlet | resid Normal | resid Oblique | Oblique vs Normal |
|---|---|---|---|---|
| 0° | 2.47e+00 | 1.09e-04 | 7.15e-04 | +558% |
| 15° | 2.47e+00 | 3.65e-04 | 7.28e-04 | +99% |
| 30° | 2.47e+00 | 8.91e-04 | 7.91e-04 | −11% |
| 45° | 2.48e+00 | 1.62e-03 | 1.13e-03 | −31% |
| 60° | 2.54e+00 | 2.76e-03 | 1.35e-03 | **−51%** |

Three readings, in order of how much they should change our behaviour:

1. **Radiation versus Dirichlet is the giant effect** -- three to four orders of magnitude. We already
   have radiation. Everything below is a refinement on top of a decision already made correctly.
2. **The schemes cross over near 30°.** `NormalRadiation`'s residual grows steeply with angle
   (1.1e-4 → 2.8e-3, a factor 25); `ObliqueRadiation`'s is nearly flat (7.2e-4 → 1.4e-3, a factor 2).
   Oblique buys a **flatter angle dependence at the cost of a higher floor** -- worse than normal for
   near-normal incidence, half the error by 60°.
3. **Our case sits on the favourable side.** The release plume leaves the southern boundary heading
   southwest, i.e. tens of degrees off the boundary normal, which is where oblique wins.

θ = 0 is NOT a control here: the schemes coincide when the tangential GRADIENT OF THE FIELD vanishes,
not the tangential velocity, and a blob has dc/dy != 0 at the boundary even in purely normal flow.
The exact-agreement control is Oceananigans' own unit test, which passes (11/11 in
`test/test_open_boundary_conditions_hydrostatic.jl`, including "ObliqueRadiation reduces to
NormalRadiation at normal incidence" and "is mirror-symmetric along the boundary").

## What this does not establish

- **Tracer only.** The velocity boundary conditions are untested here, and the tangential-velocity
  question above is untouched.
- **Prescribed uniform flow, `Centered()` advection, one resolution, 2-D.** Production is WENO(7) on a
  real flow with an immersed boundary, where the numerics may not carry over.
- **No nesting.** The decisive test is a truncation experiment -- a large truth run, a child sub-box
  driven from it, interior error compared against truth. `validation/open_boundaries/oscillating_flow.jl`
  in the package is the natural base: a cylinder wake whose flow direction sweeps a full 360°, so every
  open boundary sees every angle.

## Migration scope: NumericalEarth 0.6.1 → 0.8.0

Staged in `v113_staging/` (a third environment, so the main one keeps reading 0.110.17-written output
while the runs in flight are analysed). Resolves clean: NumericalEarth 0.8.0, Oceananigans 0.113.1,
CopernicusMarine unchanged at 0.2.0.

**The API surface is intact.** Every symbol the canary scripts use exists in 0.8.0:
`MultiYearJRA55`, `RepeatYearJRA55`, `JRA55PrescribedAtmosphere`, `JRA55PrescribedRadiation`,
`OceanOnlyModel`, `DataWrangling`, `ocean_simulation`, `EarthSystemModel`, plus `DW.Metadata`,
`DW.Metadatum`, `DW.DatasetRestoring`, `DW.BoundingBox`, `DW.GLORYSDaily`,
`DW.NearestNeighborInpainting`, `DW.interpolate_to_grid`. The two signatures most likely to drift
match our call sites unchanged:

    Metadata(variable_name; dataset, dates, dir, region, filename, start_date, end_date)
    DatasetRestoring(metadata, arch_or_grid; rate, mask, time_indices_in_memory, time_indexing,
                     inpainting, cache_inpainted_data, field_name)

Symbols existing is not the same as behaviour being preserved, so the actual gate is a regression
pair: the same code and inputs run at 0.1° on both environments, diffed field by field. Nothing about
the port is interesting if that diff is clean.

## NestedModels: revisit when the ocean version lands

NumericalEarth 0.8 adds a `NestedModels` submodule -- `NestedModel`, `NestedSimulation`,
`parent_boundary_conditions`, `parent_forcings`, `exchange_state!`, `reconstruct_parent_state`,
`blend_parent_terrain!`. `parent_boundary_conditions` builds `FieldBoundaryConditions` for a child's
lateral boundaries straight from a parent's `FieldTimeSeries`, taking `variables`, `sides` and
per-field `schemes` -- which is what `canary_glorys_obc.jl` hand-rolls.

**It is atmosphere-only today** (the submodule's sources mention `atmosphere` and never `ocean`), so
it does not apply to this repo yet. The ocean abstraction is planned upstream; when it lands it should
replace our hand-written boundary wiring with a function call. Keep the hand-written version until
then, and re-check this on each NumericalEarth bump.

## Truncation experiment: the tangential boundary condition is the whole effect

`truncation_truth.jl` then `truncation_child.jl`. Truth is a doubly-periodic 2-D domain -- periodic
means it has **no open boundaries anywhere**, so every point in it, including the ring where the
child's boundaries sit, is uncontaminated. The child is a sub-box on the *same* spacing, so boundary
values transfer by index with no interpolation error of their own. Any difference the child shows is
caused by its open boundaries and nothing else.

Two design mistakes were needed to get a meaningful measurement, and both are worth recording because
they are easy to repeat:

**1. A perfect parent makes the experiment vacuous.** Given exact truth on the child's own grid at
output cadence, a hard Dirichlet boundary is optimal by construction, and every radiation scheme
loses to it (dirichlet 8.8e-3 against normal 2.4e-2). Radiation exists to let features leave when the
exterior is WRONG. The parent must therefore be degraded -- `PARENT_R` coarsens it in space,
`PARENT_K` in time -- mimicking GLORYS at 1/12° sampled daily against a 0.02° child.

**2. Degrading the parent is not enough if the child has no dynamics of its own.** Initialised from
truth and driven by truth, the child has nothing to radiate, so 8x/10x degradation barely moved the
result (9.4e-3 against 2.5e-2). The fix is an obstacle inside the child's footprint plus a mean flow:
truth resolves the wake, the coarsened parent cannot, and the wake has to leave through the boundary.

**3. Pointwise error is the wrong metric for a chaotic wake.** Once the wake goes chaotic the child
decorrelates from truth whatever the boundary does, and that swamps the boundary signal -- every arm
landed at RMS ≈ 0.25, a quarter of the mean flow. The enstrophy RATIO against truth is statistical
and survives decorrelation: above 1 the boundary manufactures vorticity, below 1 it damps it.

### Result

Enstrophy ratio against truth over the child interior, ideal 1.0, mean over the run:

| arm | 8x coarser / 10x sparser | 4x / 5x |
|---|---|---|
| dirichlet | 0.183 | 0.197 |
| normal (normal components only -- what `canary/` does) | 0.187 | 0.203 |
| **normal + tangential** | **0.246** | **0.281** |
| oblique (normal components only) | 0.184 | 0.201 |
| **oblique + tangential** | **0.240** | **0.277** |

**Adding a tangential-velocity open boundary condition is worth +32% to +39% of retained enstrophy.
Switching NormalRadiation to ObliqueRadiation is worth nothing measurable** -- 0.184 against 0.187
without the tangential BC, 0.240 against 0.246 with it, i.e. within noise and if anything marginally
favouring `NormalRadiation`. The ordering is identical at both degradation levels.

This inverts the priority. `canary/` sets no tangential open boundary condition at all (`u_bcs` is
west only, `v_bcs` south and north -- the boundary-normal component in each case), and **a tangential
radiation BC already works on Oceananigans 0.110.17**, verified directly. So the change that appears
to matter needs no migration, while the migration buys a scheme change that did not register here.

### What this does not establish

- **2-D, idealised, one flow.** A cylinder wake in decaying turbulence is not the Canary upwelling.
- **Every arm is heavily damped** (ratio 0.18-0.28, so 70-80% of truth's enstrophy is lost whatever
  the boundary does). That is a stressed regime, and the ranking might not survive in a gentler one.
- **Velocities only.** Tracers -- the thing the domain-sizing study actually measures -- are not in
  this metric, and the separate `oblique_tracer_outflow.jl` test did favour `ObliqueRadiation` at
  incidence angles beyond ~30°. The two results are not in conflict: different fields, different
  regimes, and the tracer test had an analytic truth this one cannot have.

## Migration diff, and what the time convention does and does not explain

Same code, same inputs, 0.1° for 2 days. `v113t0` is 0.113 with `CANARY_PARENT_T0=0`, which restores
0.6.1's indexing and so isolates the time-convention change from everything else.

| field | 0.110.17 | 0.113 | 0.113, old origin | explained by the 12 h shift? |
|---|---|---|---|---|
| v (min) | −0.53849 | −0.45778 | **−0.53094** | **yes** — 15% gap closes to 1.4% |
| v (abs mean) | 2.9134e−2 | 2.8652e−2 | 2.9001e−2 | mostly — 1.65% to 0.46% |
| u (abs mean) | 3.0205e−2 | 2.9743e−2 | 2.9960e−2 | partly — 1.53% to 0.81% |
| T (mean) | 9.75984 | 9.72780 | 9.72794 | **no** — unchanged by the origin |
| T (max) | 21.638 | 21.518 | 21.518 | **no** |
| c (abs mean) | 3.1498e−9 | 3.2497e−9 | 3.2606e−9 | **no** |

So the split is clean: **the dynamics differences are the time convention, the thermodynamic ones are
not.** Restoring the old origin recovers the velocity extremes almost exactly, and does nothing at
all for temperature or the tracer. T's maximum sits 0.12 °C lower on 0.113 whatever the parent
timestamps say, which means a second, unrelated change somewhere in Oceananigans 0.110.17 → 0.113.1
or NumericalEarth 0.6.1 → 0.8.0 — surface fluxes and the JRA55 coupling are the obvious suspects.

That second change is now localised: **it is in the initial condition, not the physics.** Peak
temperature at *iteration 0* is 21.691 on 0.110.17 and 21.570 on 0.113, and at t = 2 days it is
21.638 against 21.518 — the same 0.12 °C offset, constant rather than growing. So nothing accumulates
over the run; the parent is simply regridded onto the child grid slightly differently. Surface fluxes
and the JRA55 coupling, the obvious first suspects, are not involved. Both 0.113 runs agree exactly
here, which also confirms the offset is independent of the timestamp convention — the initial
condition is built from the first parent frame either way.

A constant initial-condition offset of this size is benign for the runs this repo does, but it should
be confirmed as a regridding change rather than a data-selection change before `v113_staging`
replaces the main environment.
