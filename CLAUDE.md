# nesting_implementation

Regional ocean configurations in Oceananigans, testing open boundary conditions driven by a
parent product. Currently: `salish_sea/` (Salish Sea / Juan de Fuca, GLORYS-forced).

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
  installed at the package's scratch path. See `docs/`.
- Long runs: ~12 s/iter at 250x177x36 on CPU, roughly 1.6 h per simulated day.

## Docs

- `docs/obc-implementation-recipe.md` — portable recipe for parent-driven open boundaries: the BC
  wiring that works, the traps, the net-transport correction, the diagnostics worth having.
