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
  and matches CPU (see `canary/README.md` and `ENVIRONMENT.md`); 0.02° not yet run on GPU.

## Working agreements

- **Write notes and findings to this repo, never to Claude's memory store.** Instructions that
  should persist across sessions belong in this file.
- **`docs/` is the external user guide only.** It is read by people who did not run these
  experiments, so it holds portable, durable guidance and nothing else. Do **not** add
  investigation notes, measurement studies, or machine-specific traps to it.
- **Findings go where they can't drift from the code that made them.** A measurement belongs in
  the header of the script that produces it; the reasoning behind it belongs in the commit message
  that lands it; a correction to a wrong claim belongs at the point of the claim. Domain-specific
  notes belong in that domain's `README.md`, operational traps in `ENVIRONMENT.md`.
- **Keep code comments lean.** Annotate the non-obvious decision in 1-3 lines and stop; keep
  investigation history out of source. Delete comments a change has falsified.

## Environment

- **Julia 1.12** — now the juliaup *default* channel, so any directory works. (It used to be a
  per-directory override on the old `salish_obcs` path; renaming the repo orphaned that and
  silently dropped `julia` to 1.10, where precompilation dies with a `StaticData` error.)
- **Oceananigans ≥ 0.110.17** (after `ss/open-boundary-conditions` merged) and **NumericalEarth 0.6.1**.
  No git pin — do not reinstate a `[sources]` entry for Oceananigans.
- A transient `signal 11` segfault on the first `using NumericalEarth` just needs a retry.
- `CopernicusMarine.jl`'s bundled binary does not run on this macOS (13.2); a pip-venv shim is
  installed at the package's scratch path — see `ENVIRONMENT.md`, which also covers the data
  caches, the GLORYS download gotchas, and run cost.

## Docs

The user guide, for external readers:

- `docs/obc-implementation-recipe.md` — portable recipe for parent-driven open boundaries: the BC
  wiring that works, the traps, the net-transport correction, the diagnostics worth having, and a
  cross-check against the Mid-Atlantic Bight configuration in NumericalEarthSandbox.

Everything else is internal and lives next to what it describes:

- `ENVIRONMENT.md` — Julia/juliaup, the CopernicusMarine macOS shim, caches, output buffering,
  GCP GPU runs, run cost.
- `salish_sea/README.md` — that configuration: settings, validation, the instability root cause,
  why the domain boundaries sit where they do, and known limitations.
- `canary/domain-sizing.md` — how far south the Canary domain must reach: the virtual-boundary
  method, retention over 73-268 km, the tau(L) fit, and why extending south past 22.8N is wasted
  unless the west boundary moves too.
- `red_sea/README.md` — that domain’s verified geometry, the GLORYS depth checks (boundary *and*
  interior, the latter via `red_sea/glorys_interior_check.jl`), and what is untested.

## Where this stands (2026-10-28)

- **Red Sea** — the three-arm scheme comparison is done (`e4eafc5`); the recommended production
  configuration is in `red_sea/README.md`. **Unmeasured:** the solution's sensitivity to the
  `TracerReservoir` inflow length `L`. The test is two 2-day 0.05° arms at `L` = 20 km and 60 km
  against the recorded `L` = 5 km arm. Background and the cross-domain spread:
  `parent_excursion_survey.jl` and commit `fd77716`.
- **Salish** — the tangential BC and oblique/reservoir scheme selection were added in `a50ac90`
  but have **never been run**. The baseline is still the original validated 7-day run (its output
  holds 56 snapshots ending at exactly 7.0 days). Worth doing: at 8.3× refinement Salish is the
  decisive test of the hypothesis in `red_sea/README.md` that oblique radiation's failure tracks
  refinement rather than land geometry — Canary failed at 4.2×, Red Sea survived at 1.7×.
- **PR 676 (`nested_ocean_model`)** — `pr676_staging/redsea_nested.jl` is corrected but not yet
  re-run. The first attempt was **invalid**: `Nz`/`scale` of 40/1800 against the baseline's 32/480
  gave 4 levels in the top 150 m instead of 17, which halves the strait transport, because the
  Hanish sill is 136 m. Separately, PR 676's defaults differ from the baseline on four axes —
  inflow timescale 3days vs 6hours, outflow 360days vs `Inf`, reservoir `L` 20km vs 5km,
  relaxation 1/6h vs 1/12h — so the schemes must be pinned explicitly to compare implementations
  rather than configurations. They now are. What the invalid run did establish: a clean 2-day exit
  and 25.7 s/iter against the baseline's 32.3, so the exchanger adds no meaningful overhead.
