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

## Net-transport correction

Parent velocities interpolated onto a different (usually better-resolved) bathymetry do **not** conserve volume. Standard fix: subtract a spatially uniform normal-velocity `Δu = (Q_raw − Q_target)/A` from the boundary exteriors, applied to **both** the 3-D normal flow and the depth-integrated Flather transport so baroclinic and barotropic data stay consistent. Recompute from the **raw** parent data in a callback (idempotent, not a feedback loop); `Q_target = 0` with no rivers, negative (net outflow) once rivers are added. Δu can look tiny and still be doing enormous work — multiply by the boundary area before dismissing it.

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
hand-rolling a virtual salt flux through `additional_surface_fluxes`. Note the sign convention on any
net-transport correction must then become a net **outflow** equal to the discharge.

## Method lessons

- **Validate at production resolution.** A coarse smoke test inverted the sign of the boundary through-flow and showed no instability where the fine grid had one. Coarse runs are for mechanical smoke tests only.
- **Before trusting a negative result, verify the fix actually reached the thing you tested.** A "not the cause" conclusion was wrong for days because the repair never touched the initial condition.
- **Change one thing per run**, and tag outputs by configuration so baselines survive.
- Bounded 1-hour tests are enough for stability; **transport and water-mass questions need days** (advective transit from boundary to interior).
