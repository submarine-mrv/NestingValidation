# nesting_implementation

Regional ocean configurations in Oceananigans, testing open boundary conditions driven by a
parent product.

- `salish_sea/` — Salish Sea / Juan de Fuca, GLORYS-forced. Validated: stable 7-day run at 0.74 km.
- `red_sea/` — Red Sea + Gulf of Aden, GLORYS-forced, optional JRA55 atmosphere. Validated both ways:
  stable 7-day runs at 5.2 km (17.5 h unforced, 22.9 h forced). Strait transports not yet converged.
- `examples/` — minimal, self-contained open-boundary examples to share: one open boundary
  (Red Sea) and three (Canary), with the pattern and its traps in `examples/README.md`.
- `canary/` — Canary upwelling, NW Africa. Three open boundaries, wind-driven interior, islands
  inside the domain. Validated in both seasons at 2 km (January 38.8 h, July 36.3 h). Runs on GCP GPU
  and matches CPU (see `canary/README.md` and `docs/environment.md`); 0.02° not yet run on GPU.

## Working agreements

- **Write notes and findings to this repo, never to Claude's memory store.** Durable notes belong
  in `docs/`; instructions that should persist across sessions belong in this file.
- **Keep code comments lean.** Annotate the non-obvious decision in 1-3 lines and stop; investigation
  history goes in `docs/`, not in source. Delete comments a change has falsified.

## Environment

- **Julia 1.12** — now the juliaup *default* channel, so any directory works. (It used to be a
  per-directory override on the old `salish_obcs` path; renaming the repo orphaned that and
  silently dropped `julia` to 1.10, where precompilation dies with a `StaticData` error.)
- **Oceananigans ≥ 0.110.17** (after `ss/open-boundary-conditions` merged) and **NumericalEarth 0.6.1**.
  No git pin — do not reinstate a `[sources]` entry for Oceananigans.
- A transient `signal 11` segfault on the first `using NumericalEarth` just needs a retry.
- `CopernicusMarine.jl`'s bundled binary does not run on this macOS (13.2); a pip-venv shim is
  installed at the package's scratch path — see `docs/environment.md`, which also covers the data
  caches, the GLORYS download gotchas, and run cost.

## Docs

- `docs/obc-implementation-recipe.md` — portable recipe for parent-driven open boundaries: the BC
  wiring that works, the traps, the net-transport correction, the diagnostics worth having, and a
  cross-check against the Mid-Atlantic Bight configuration in NumericalEarthSandbox.
- `docs/salish_sea_notes.md` — this configuration: settings, validation, the instability root cause,
  why the domain boundaries sit where they do, and known limitations.
- `docs/environment.md` — Julia/juliaup, the CopernicusMarine macOS shim, caches, run cost.
- `docs/canary-domain-sizing.md` — how far south the Canary domain must reach: the virtual-boundary
  method, retention over 73-268 km, the tau(L) fit, and why extending south past 22.8N is wasted
  unless the west boundary moves too.
- `red_sea/README.md` — that domain’s verified geometry, the GLORYS depth checks (boundary *and*
  interior, the latter via `red_sea/glorys_interior_check.jl`), and what is untested.
