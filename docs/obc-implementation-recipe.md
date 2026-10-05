# Open boundary conditions: implementation recipe

Distilled from the Salish Sea / GLORYS build (Jul–Aug 2026, `salish_sea/salish_glorys_obc.jl`).
This is the **transferable** part: what to do in ANY regional config driven by a parent product.
Names below are for **Oceananigans ≥ 0.110.17**, i.e. after `ss/open-boundary-conditions` merged.

## The wiring that works

- **u/v (3-D normal flow):** `NormalFlowBoundaryCondition(closure; scheme = NormalRadiation(inflow_timescale=6hours, outflow_timescale=Inf), discrete_form = true)`
- **U/V (barotropic transport, Flather):** `GravityWaveRadiationBoundaryCondition(closure; discrete_form = true)` — the closure returns a **2-tuple `(U_ext, η_ext)`**
- **η (Chapman):** `SurfaceWaveRadiationBoundaryCondition()`
- **T/S:** `ValueBoundaryCondition(closure; scheme = radiation, discrete_form = true)`
- `ocean_simulation` passes `u, U, v, V, η, T, S` boundary_conditions straight through.

**Renames at the merge** (pre-merge branch → 0.110.17), in case you meet older code:
`Radiation` → `NormalRadiation`; `ImplicitGravityWaveRadiation(BoundaryCondition)` →
`SurfaceWaveRadiation(BoundaryCondition)`; `radiation_scheme.jl` → `normal_radiation.jl`.
`GravityWaveRadiation` and `NormalFlowBoundaryCondition` were unchanged, the Flather value is still
a 2-tuple, and the Chapman halo fill is byte-identical. Upgrading reproduced a 1-hour res-3 run to
3–4 significant figures with identical salinity, at the same cost per iteration.

**Flather is exactly:** west `U = Uᵉˣᵗ − √(gH)(ηᵇ − ηᵉˣᵗ)`, north `V = Vᵉˣᵗ + √(gH)(ηᵇ − ηᵉˣᵗ)`, with ηᵇ the face-average of adjacent interior η. Sign errors here are easy and expensive — derive before trusting.

## Traps that cost real time

1. **The Chapman BC never calls `getbc`.** Its halo fill is `η[0,j,k] = (ηᵇ + C·η₁)/(1+C)`, pure radiation of interior η outward. The `0` in `ImplicitGravityWaveRadiationBoundaryCondition()` is ignored, not applied. The **Flather 2-tuple is the only route an exterior η reaches the solution.**
2. **Reuse the `DatasetRestoring` FieldTimeSeries for the BC closures.** Those live in `model.forcing`, so `update_model_field_time_series!` advances them automatically. A standalone FTS held only by a closure is **opaque to that traversal and will never advance in time.** (Requires defining the restorings BEFORE the BCs.) For a field with no restoring (e.g. sea surface height), hold **all** time indices in memory instead.
3. **`DatasetRestoring` defaults to `time_indices_in_memory = 2`** and silently reloads slices from the cached files as the clock advances — which undoes any in-memory repair you made. Pass the full length if you patched the data.
4. **`set!(model.tracers.T, …)` fails once BCs are clock-dependent** (`DiscreteBoundaryFunction … not callable`): setting triggers a halo fill needing `clock.time`, which plain `set!` has none. Fix: set into scratch `CenterField(grid)` with default BCs, then `interior(model.tracers.T) .= interior(T_ic)`. Same for η and velocities.
5. **`set!(field, ::Metadatum)` re-reads the cached inpainted file** and does NOT go through any FieldTimeSeries you repaired. If the parent data needs fixing, **repair in ONE place and have every consumer read from there** — build ICs by sampling the repaired FTS via `interpolate_to_grid`, not by a second independent path. This one produced a multi-day wild goose chase.
6. **Parent data is zero-filled below its deepest valid level**, and those zeros ARE sampled at wet model points wherever your bathymetry is deeper than the parent's. Interpolation then blends a valid value with zero and yields *anything in between* — so an absolute floor cannot clean it. Repair by **vertical fill-down** (replace sub-floor values with the deepest valid value above) and check the parent's actual validity depth per level before trusting a boundary that reaches below it.
7. **Initialize η from the parent's sea surface height, not from rest.** Starting at η=0 while the boundary says +0.2 m manufactures a basin-filling transient lasting hours (L/√(gH) sets the timescale) that swamps every bounded test. This was worth ~100× in the volume budget.
8. **`Cyclical()` time indexing wraps the start of the run to the END of the parent record.** NumericalEarth 0.8 stamps `GLORYSDaily` at the **centre** of its daily-mean window (12:00), so a series built from `run_start : Day(1) : …` has its first sample 12 h *after* t = 0. `Cyclical()` does not persist across that gap — it wraps, so the **last day of the record drives the first half-day of the run**, at the open boundary *and* in the sponge. Pass `time_indexing = Clamp()` instead (persistence), to `DatasetRestoring` as well as to any standalone FTS: **`DatasetRestoring` defaults to `Cyclical()` internally**, and under trap 2 its series *is* the boundary exterior, so clamping only the standalone one fixes nothing visible. The lead-in day the upstream warning suggests is not an option for GLORYS — the record begins 1993-01-01, which is the usual `run_start`. Invisible before the migration: under NumericalEarth 0.6.1 the stamp was midnight, t = 0 sat exactly on the first sample, and there was no gap to fill.

## Net-transport correction — measured, and removed

Parent velocities interpolated onto a different (usually better-resolved) bathymetry do **not**
conserve volume. The standard fix is to subtract a spatially uniform normal-velocity
`Δu = (Q_raw − Q_target)/A` from the boundary exteriors, applied to **both** the 3-D normal flow and
the depth-integrated Flather transport, recomputed from the **raw** parent data in a callback.

**In the Red Sea domain this was measured and found not to matter, and has been removed.** Two
2-day runs at 0.05° whose prescribed exterior nets differed by **3.0e5 m³/s** were indistinguishable:
strait exchange within 0.6%, `max|η|` 0.3911 vs 0.3903 m, `max|u|` 0.6620 vs 0.6621 m/s, volume drift
within 1.5%. The no-correction case lies *between* those two, so it is bracketed without a third run.

The mechanism is that **Flather is a restoring feedback, not an integrator**:
`U = Uᵉˣᵗ + √(gH)(ηᵇ − ηᵉˣᵗ)`. To source a 1.5e5 m³/s residual, η settles ~4 mm below the exterior and
supplies it indefinitely — no secular accumulation, which is what the volume drift shows (it
oscillates and returns rather than ramping). The same argument covers evaporation: a Red Sea
freshwater demand of ~1.3e4 m³/s is an *order of magnitude smaller* than the residual and needs
~0.2 mm of offset, so the boundary supplies it unaided, seasonally too.

There is a second, sharper piece of evidence. The Red Sea script carried a **sign error** in this
correction through its entire 7-day validation run — the formula was negated for the east boundary
*and* the correction still added, inverting twice, so the exterior was driven to `2Q − Q_target`,
doubling the imbalance it was meant to remove. Nothing in the solution revealed it. A component whose
sign can be inverted for a full validation run with no detectable symptom is not doing much work.

**Keep `Q_raw` as a diagnostic.** It is a parent/child mismatch measure, not a physical signal, and it
**flips sign with resolution** — `+3.15e5 m³/s` at 0.2° against `−1.49e5` at 0.05° for the same domain
— because parent and child disagree about bathymetry and cell geometry. It also varies day to day with
the parent fields. Reporting it tells you how badly the two grids disagree; applying it does not
appear to buy anything.

**Caveat on generality.** This was measured on a domain with an unusually forgiving boundary: a
719 km² section carrying 24–45 Sv of gross exchange, so the residual is ~0.4% of throughflow. The
Salish Sea west boundary was 1.6 km². The restoring argument does not depend on section size, so the
same conclusion is expected there, but it is **untested** — `salish_sea/` still applies the
correction. Do not assume a rigid-lid or prescribed-transport boundary behaves this way; the whole
argument rests on having a radiating barotropic mode and a free surface.

## Diagnostics that earned their place

- **Argmax localisation**: report `(i,j,k)` + depth + a WEST-bnd/NORTH-bnd/interior tag for the extrema. Turns "something is growing" into "it is at the bottom of the deepest basin, interior" in one line.
- **Wet-only extrema** (land cells hold zeros that corrupt minima).
- **Excursion referenced to the FORCING DATA range, not the IC** — the boundary can legitimately import water outside the initial range; referencing the IC reports physical adjustment as spurious.
- **Section transports** per boundary, plus the implied domain-mean sea-level tendency in m/hour.
- **Exchange flow**: accumulate positive and negative `u·A` separately, with flux-weighted salinity of each limb → Qin, Qout, Sin, Sout, ΔS and a Knudsen residual (as % of Qin·Sin). The **net** tells you almost nothing about whether the exchange is right.
- **Boundary η vs parent η mismatch** — the direct read on the Flather term.
- **Wall-clock pacing with a windowed rate**; use s/iter (stable from the first callback) for planning, not the ETA (pessimistic until Δt caps).

## Adding freshwater sources (rivers)

NumericalEarth ≥ 0.6.1 provides `FreshwaterExchange` in `ocean_simulation`: each source carries a
prescribed tracer concentration (zero salinity for pure water, its own temperature for heat), and
only the live surface value is read, which cancels the z\* "ambient carry" to machine precision.
Adding a source is a two-line change there (`carrying += Jʷ`, `content += c·Jʷ`). Prefer this over
hand-rolling a virtual salt flux through `additional_surface_fluxes`. If you reinstate a
net-transport correction, its target must then become a net **outflow** equal to the discharge — but
see above for whether it is worth having at all.

## Method lessons

- **Validate at production resolution.** A coarse smoke test inverted the sign of the boundary through-flow and showed no instability where the fine grid had one. Coarse runs are for mechanical smoke tests only.
- **Before trusting a negative result, verify the fix actually reached the thing you tested.** A "not the cause" conclusion was wrong for days because the repair never touched the initial condition.
- **Change one thing per run**, and tag outputs by configuration so baselines survive.
- Bounded 1-hour tests are enough for stability; **transport and water-mass questions need days** (advective transit from boundary to interior).

## Cross-check: the Mid-Atlantic Bight configuration

`NumericalEarthSandbox/scripts/02_mab_glorys_obc.jl` (enri66) is an independent GLORYS-forced open
boundary configuration. Comparing it against `canary/` is the only external check this recipe has
had, and it is worth keeping current as both evolve.

### Four choices agree, independently

| | MAB | canary/ |
|---|---|---|
| barotropic U/V | `GravityWaveRadiationBoundaryCondition` (Flather) | same |
| free surface η | `SurfaceWaveRadiationBoundaryCondition` (Chapman) | same |
| tracers T/S | `ValueBoundaryCondition` + `NormalRadiation`, outflow `Inf` | same |
| Flather exterior | depth integral through the **model's** wet column at model z-levels | same |

The last one is the subtle one and the agreement is the most reassuring part: MAB offers three ways
to build `Uᵉˣᵗ` (`native` over GLORYS's own vertical grid, `masked` through the model wet column,
`legacy` unmasked) and defaults to the masked form, which is what `west_U_obc` and friends here have
always computed. Two implementations converging on that independently is decent evidence it is right.

### Four differences, in the order they matter

1. **MAB sets TANGENTIAL velocity boundary conditions; this configuration never did.** Theirs is
   `ValueBoundaryCondition(Interpolated(fts); scheme = NormalRadiation(...))` on the tangential
   component at every boundary, with a `prescribed` mode that swaps it for stiff Dirichlet. That is
   independent corroboration of the truncation experiment in `obc_v113/README.md`, which put the
   tangential condition at +32-39% of retained enstrophy, and of Oceananigans PR #5962, which
   attributes its headline result to the same thing. Three lines of evidence, three different routes.
   `canary/` now sets it unconditionally, in exactly that form.

2. **Normal velocity: MAB uses `PerturbationAdvection`, this uses the radiation scheme.** They apply
   different schemes to the normal component (`PerturbationAdvection`) and to the tangential
   component and tracers (`NormalRadiation`); `canary/` used one scheme throughout.
   `normal_scheme: ["perturbation"]` reproduces their split. Available on 0.110.17 as well as
   0.113, so this needs no migration.

3. **MAB has no sponge at all** — interior restoring was removed and all forcing enters through the
   boundaries. `canary/` runs a 0.25° Gaussian at a 12 h timescale. The narrow width here is
   deliberate (the Salish lesson: a wide sponge pins the features the domain exists to resolve), but
   MAB went further and dropped it entirely. A sponge-off arm is cheap and worth measuring.

4. **Inflow timescale: MAB 1 day, `canary/` 6 hours.** Ours is roughly 4x stiffer on inflow.

### Read their tuning with care

MAB appears to run its child at **1/12°, the same resolution as GLORYS**, so it is not really
downscaling: the child resolves little the parent cannot supply, and the open boundaries carry
correspondingly less stress. `canary/` runs 0.02° against a 1/12° parent, a ~4.5x refinement, which
is the harder regime and precisely the one where the truncation experiment found the tangential
condition mattering. Treat the agreement above as a strong check on the **wiring**, and their
timescales and sponge choice as a weaker guide to **tuning** here.

## Shelved: upstreaming parent-bathymetry matching into NumericalEarth

Trap 6 above — child bathymetry deeper than the parent's valid data — is fixed here by a *global*
cap (`depth: [4800]`), which is the wrong generalisation twice over. The constraint binds on the
open-boundary rim, where the deficit actually bites (the west/north corner, per
`canary/glorys_interior_check.jl`), yet the cap clips 1.2% of all columns by up to 85 m. Worse, it
couples: `z_reference` spans `0..-depth`, so lowering `depth` restretches every vertical level in
the domain to satisfy a constraint on a thin band.

`enri66/NumericalEarthSandbox`'s `scripts/glorys_bathymetry.jl` has the right shape instead —
`match_boundary_bathymetry!(bottom_height, grid, data_dir; region, sides, n_match, minimum_depth)`,
a raised-cosine blend of the child's bottom toward the parent's over `n_match` cells inward from the
named sides. Three things it gets right: it reads GLORYS's **static `deptho`** rather than the
valid/NaN mask of an inpainted tracer (the mask says where values survive inpainting, not where the
seafloor is); it blends either direction, where a cap can only shallow the child; and it preserves
both land masks rather than blending them into a meaningless fractional mask.

If this is picked up, the gaps to close, in order:

1. `n_match` is a fixed cell count. The atmosphere side already hit this and documented it:
   `default_terrain_blend_width` exists because a fixed count steepens the parent→child transition
   as ~1/Δx and destabilised high-resolution runs at the boundary corner. Bottom slope has the same
   geometry, and 0.02° to 1/60° spans 5× here. Take a physical blend length instead.
2. It is purely 2-D, so it does **not** settle trap 6. Matching `deptho` makes the two seafloors
   agree; feasibility is whether every child cell needing a parent value has one, which depends on
   the child's z levels against the parent's valid ones. Those can disagree with `deptho` matched.
3. `deptho` is a CMEMS convention — the function's best property is its least portable one. A
   parent-agnostic version needs a per-dataset `parent_bottom_height(metadata)` trait whose fallback
   errors honestly rather than silently reaching for an inpainted field's mask.

Where it would belong: the repair as `NestedModels/parent_bathymetry.jl`, beside the existing
`blend_parent_terrain!` (same three-layer split — grid-level kernel, model-specific orchestration,
resolution-aware default); the vertical check as a separate `validate_parent_depth_coverage` in the
`DataWrangling` `validate_dataset_coverage` family, fatal at open boundaries and warn-and-inpaint in
the interior. Not inside `regrid_bathymetry`: it caches on `bathymetry_regridding_key`, which knows
nothing about a parent dataset or which sides are open, so a cached field could silently be one
matched to a different parent. Two functions, not one — and only the second would have flagged
Canary before the 39-hour run. Coastline/mask mismatch is a third, topological problem that shares
the diagnosis layer and nothing else.

## ObliqueRadiation breaks on a lateral open boundary that touches land

`ObliqueRadiation` (Oceananigans 0.113, CliMA/Oceananigans.jl#5962) diverges in this configuration at
0.02°, and the cause is specific: **the scheme cannot be used on a lateral open boundary that
contains any land cells.** It is not a resolution problem and not a strong-gradient problem in
general.

### What happened

A 90-day run with `radiation: ["oblique"]` on all three boundaries reached **S = 45.7 psu** (an
excursion of +8.4 beyond the GLORYS range) within **36 simulated minutes**, with `max|η|` at 2.0 m
against a normal 0.35 m and worst Δρ oscillating to 5.2 kg/m³. The reported location was cell
**(332, 1)** — j = 1 is the south boundary, and i = 332 is its **last wet cell**, with land at i = 333.

### The controlled test

Three one-day runs at 0.02°, identical but for which boundary carries the oblique scheme
(`oblique_sides`):

| boundary | land cells | max S | worst Δρ | result |
|---|---|---|---|---|
| west | **0** / 225 | 37.30 | 0.20 | **stable**, full day, exit 0 |
| north | **6** / 575 (1%) | **40.63** | 0.77 | diverged, +3.3 psu |
| south | **243** / 575 (42%) | **47.06** | 0.91 | diverged, +9.8 psu |

Three things follow. **Any land contact is enough** — six cells out of 575 breaks it. **Severity
scales with the land fraction**, monotonically. And **resolution is not the cause**: the west boundary
runs the same 2 km grid and carries the strongest transport in the domain (the Canary Current,
−6.9e6 m³/s) with no trouble at all.

### Why

The oblique scheme adds a term `NormalRadiation` does not have:

    cₜ = −(∂φ/∂t)(∂φ/∂t̂) / |∇φ|²

`t̂` is the boundary-TANGENTIAL direction, which on the south boundary runs along the coast. At a wet
cell whose along-boundary neighbour is land, that derivative is evaluated across the land/water
interface, `cₜ` becomes meaningless, and the tangential advection term injects tracer. The docstring
says the tangential term "is applied on the lateral boundaries only" and says nothing about a lateral
boundary that terminates in land.

Note the pattern: the bounds-preserving WENO problem recorded in `canary-domain-sizing.md` was also
an immersed-boundary interaction, and that file's own source comment hedges on immersed-boundary
safety. **New boundary/advection schemes in Oceananigans are worth testing against land before use.**

### Correction: oblique DOES work with land elsewhere — on velocities, not tracers

An earlier version of this section concluded that ObliqueRadiation is simply unusable on a boundary
touching land. That is too strong, and the counter-example is the Mid-Atlantic Bight configuration
itself. `NumericalEarthSandbox/scripts/04_mab_glorys_tides_reservoirs.jl` — not `02_…`, which is the
older script — runs oblique on a domain whose open boundaries are cut by the US east coast, and works.

Three differences matter, and they point somewhere more specific than "land":

1. **MAB applies oblique to VELOCITIES, not tracers.** Its default tracer scheme is `TracerReservoir`;
   oblique is the default for the boundary-normal velocity and optional for the tangential one. Our
   configuration fed one `radiation` object to everything including T and S — and **what blew up was
   salinity**. So the failure implicates oblique on TRACERS at a land-terminated boundary, which MAB
   never exercises.
2. **MAB runs at 1/12°, the same resolution as GLORYS** — no refinement. `canary/` runs 0.02° against
   a 1/12° parent, ~4.5x. The interior carries structure at the boundary that the exterior cannot
   supply, which is exactly what makes a tangential derivative there ill-posed.
3. **MAB passes `phase_speed_weight`, which does not exist in released 0.113.3.** They track an
   unreleased branch, so the scheme has damping we do not have. (0.113.3 does have `target_transport`,
   absent from 0.113.1, so it is still moving.)

### The stencil, from the source

Released v0.113.4 computes the along-boundary difference with index clamping and no land check:

    @inline function tangential_differences(φ, t, k, b)
        T = size(φ, 1)
        @inbounds begin
            φ₀ = φ[t, k, b]
            φ₋ = φ[max(t - 1, 1), k, b]
            φ₊ = φ[min(t + 1, T), k, b]
        end
        return φ₀ - φ₋, φ₊ - φ₀
    end

Nothing in the file references immersed boundaries, masking or inactive cells. At a wet boundary cell
whose along-boundary neighbour is land, `φ₊` or `φ₋` is whatever occupies the dry slot, so the
difference is meaningless — which is the mechanism behind the dose-response above.

Not raised upstream, deliberately: the constraint is likely already understood by the handful of
people exercising the PR, and the `NestedOceanModel` abstraction coming to NumericalEarth is expected
to set these schemes automatically with fuller documentation. Recorded here because it constrains
what `canary/` can use today.

`phase_speed_weight`, which the MAB scripts pass, exists in no release (0.113.4 included) and would
damp the symptom rather than mask the stencil, so it is not the lever.

### Consequence here

`oblique_sides: ["west"]` is the only safe setting *with oblique applied to everything*. The
untested configuration, and the one MAB's success points to, is **oblique on the velocities with
`TracerReservoir` or `NormalRadiation` on T/S**. That is the next thing to try before concluding
anything about oblique in this domain.

`radiation: ["normal"]` is unaffected, and is what every run before 0.113 used: the 7-day 0.02° shakedown with
everything else identical ran stably (max S 37.0, worst Δρ 0.27).

**A run can diverge and still exit 0.** The south arm completed its full day and returned success
with S at 44.2 psu. Check the salinity excursion, not the exit code.
