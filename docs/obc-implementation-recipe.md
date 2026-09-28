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
   `CANARY_NORMAL_SCHEME=perturbation` reproduces their split. Available on 0.110.17 as well as
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
