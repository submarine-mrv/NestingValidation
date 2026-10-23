# Nested regional ocean examples — open boundaries driven by GLORYS

Two self-contained scripts showing how to drive a regional Oceananigans/NumericalEarth ocean from a
coarser parent product through open lateral boundaries. They are distilled from three validated
configurations in this repository; the working versions carry diagnostics, parameter sweeps and
investigation history that these deliberately omit.

| file | what it adds |
|---|---|
| `one_open_boundary_redsea.jl` | The minimal case. ONE open boundary, parent sampled live inside the condition closures. **Read this first.** |
| `three_open_boundaries_canary.jl` | THREE open boundaries, the crossed tangential wiring that implies, per-side scheme selection, and the parent precomputed onto the faces so the kernels are GPU-compilable. |

Both run as-is given Oceananigans ≥ 0.113.2 and NumericalEarth 0.8, plus Copernicus credentials for
the GLORYS download.

## The pattern: five conditions, not one

An open boundary is not a single boundary condition. Each open side needs all five of these, and
leaving any one out fails quietly rather than loudly:

| field | condition | carries |
|---|---|---|
| boundary-NORMAL velocity | `NormalFlowBoundaryCondition` + a radiation scheme | the baroclinic inflow/outflow |
| boundary-TANGENTIAL velocity | `ValueBoundaryCondition` + a radiation scheme | the along-boundary flow |
| barotropic transport `U`/`V` | `GravityWaveRadiationBoundaryCondition` (Flather) | the depth-integrated transport **and the exterior η** |
| free surface `η` | `SurfaceWaveRadiationBoundaryCondition` (Chapman) | nothing — see below |
| tracers `T`, `S` | `ValueBoundaryCondition` + `TracerReservoir` | the parent water mass |

Three choices in that table are not guessable from the API:

- **The Chapman condition never calls `getbc`.** Its halo fill is pure radiation of interior η
  outward, so whatever you pass it is ignored, not applied. The **Flather 2-tuple is the only route
  an exterior η reaches the solution**. Pass it, but do not expect it to impose anything.
- **The tangential component needs a condition too.** It is easy to omit, and the default is a
  no-flux wall, which is wrong on an open edge. Measured at +32–39% of retained enstrophy in an
  idealised truncation test.
- **Return the exterior value only.** Oceananigans applies the side-dependent Flather sign itself
  (west/south negative, north/east positive). Pre-flipping double-inverts it.

## Traps that cost real time

Each is marked `TRAP:` at the point it bites in the code.

1. **Reuse the `DatasetRestoring` `FieldTimeSeries` in the boundary closures.** Those live in
   `model.forcing`, so `update_model_field_time_series!` advances them. A standalone FTS held only
   by a closure is opaque to that traversal and **never advances in time**. This requires defining
   the restorings before the boundary conditions.
2. **`time_indices_in_memory` defaults to 2** and silently reloads slices from disk as the clock
   advances — undoing any in-memory repair. Hold them all if you repaired the data.
3. **`time_indexing` defaults to `Cyclical()`.** NumericalEarth 0.8 stamps `GLORYSDaily` at the
   *centre* of its daily-mean window, so the first sample lands 12 h after `t = 0` and `Cyclical`
   fills that gap by wrapping to the **far end of the record** — the last day of the series driving
   the first half-day of the run, at the boundary and in the sponge. Use `Clamp()`.
4. **Parent data is zero-filled below its deepest valid level**, and those zeros are sampled at wet
   model points wherever the child is deeper. Interpolation then blends a valid value with zero and
   yields anything between, so an absolute floor cannot clean it. Repair by vertical fill-down.
5. **`set!(field, ::Metadatum)` re-reads the cached file** and never sees that repair. Repair in ONE
   place and have every consumer read from there, initial conditions included.
6. **Assigning initial conditions directly to model tracers** triggers a halo fill that needs
   `clock.time`, which plain `set!` cannot supply once the BCs are clock-dependent. Write to a
   scratch field and copy.
7. **Initialise η from the parent, not from rest.** Starting at 0 while the boundary says +0.2 m
   manufactures a basin-filling transient lasting hours.
8. **`ocean_simulation` merges boundary conditions side by side with its defaults**, so side-only
   conditions keep the default top flux condition and surface forcing still reaches the ocean. Worth
   re-checking on a version bump: a replace-rather-than-merge would silently disable it.

## Choosing the schemes

- **Oblique radiation** (Raymond & Kuo 1984) diagnoses a tangential phase speed as well as a normal
  one, so its stencil runs *along* the boundary. On a side that terminates in coast that derivative
  straddles a land/water interface. Use it on clean sides; keep `NormalRadiation` where a boundary
  ends in land. **Never give oblique to tracers** — that combination drove salinity to 45.7 psu.
- **`TracerReservoir`** gives the boundary a memory: a plain `Value` condition hands re-entering
  water the exterior value regardless of what just left. Set the inflow length scale to ~0.3× the
  volume-weighted median parcel excursion across the section, measured from the domain's own output.
  Three different geometries measured 4.9–5.3 km, so 5 km is a reasonable first guess — but measure.
- **Site the boundary where the parent resolves well**, far from any constriction you care about.
  A boundary inside a narrow channel is prescribed from very few parent cells: in one of our domains
  the parent had only 7–9 wet cells across sections the child resolved with 176–250.

## Where the long version lives

`docs/obc-implementation-recipe.md` in this repository carries the full reasoning, the controlled
experiments behind each choice, and the net-transport correction analysis (measured, then removed).
