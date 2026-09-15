# Environment notes

Machine-specific setup and traps. Verified 2026-08-21 on macOS 13.2.1 (Darwin 22.3.0, arm64).

## Julia

**1.12 is required.** It is now the juliaup *default* channel, so any directory works:

```
juliaup default 1.12
```

Historically this was a per-directory override on the old `salish_obcs` path. Renaming the repo
orphaned that override and silently dropped `julia` to 1.10, where precompilation dies with an
unhelpful `UndefVarError: StaticData` inside `PrecompileTools`. If you ever see that error, check
`julia --version` first.

A transient `signal 11` segfault on the first `using NumericalEarth` after a while (often after a
machine restart, mid-precompile) just needs a retry — nothing is corrupted.

## CopernicusMarine: the macOS shim

**Problem.** `CopernicusMarine.jl` ships a PyInstaller binary of the `copernicusmarine` CLI that
cannot run on this machine. The root cause is the OS, not Julia: the bundled NumPy 2.1.3 links
Accelerate's ILP64 LAPACK (symbol `_cblas_caxpy$NEWLAPACK$ILP64`), which only exists on macOS 13.3+.
The failure is

```
Symbol not found: _cblas_caxpy$NEWLAPACK$ILP64 ... built for macOS 14.0 which is newer than running OS
```

Both `login` and `download`/`subset` shell out to that same binary, so credentials do not help.

**The workaround in place.** `CopernicusMarine.jl` invokes a hardcoded path
`joinpath(SCRATCH[], "copernicusmarine")` with no override hook, but it only re-downloads the binary
if that file is *missing*. So the binary was replaced with a shim pointing at a working
OpenBLAS-wheel install of the same version:

```sh
python3 -m venv ~/.venvs/copernicusmarine
~/.venvs/copernicusmarine/bin/pip install "copernicusmarine==2.4.1"

cd ~/.julia/scratchspaces/cd43e856-93a3-40c8-bc9e-6146cdce14fa/copernicusmarine-v2.4.1/
mv copernicusmarine copernicusmarine.pyinstaller.bak
cat > copernicusmarine <<'EOF'
#!/bin/sh
exec "$HOME/.venvs/copernicusmarine/bin/copernicusmarine" "$@"
EOF
chmod +x copernicusmarine
```

**It breaks if** the package bumps `TOOLBOX_VERSION` (the scratch directory name changes, a fresh
broken binary is downloaded, and you must re-shim with the venv pinned to the new version), or if
anything calls `CopernicusMarine.executable(force=true)`. The permanent fix is macOS ≥ 13.3.

**Credentials.** `COPERNICUSMARINE_SERVICE_USERNAME` / `COPERNICUSMARINE_SERVICE_PASSWORD` in the
environment, or a stored login via `CopernicusMarine.login(username=..., password=...)` (which works
through the shim). Only needed when data is actually downloaded — the caches below usually serve.

## GLORYS download gotchas

- **Any script that downloads GLORYS must `using CopernicusMarine`.** The backend lives in
  `NumericalEarthCopernicusMarineExt`, which only loads when that package is. Without it you get
  *"No download method for Metadata{GLORYSDaily...} is available (is the backend package loaded?)"*.
- The download method is `Downloads.download(::Metadata)`, **not** `Base.download` — `using Downloads`
  and qualify. Only matters when calling download directly; `DatasetRestoring` and `set!` do it
  internally.
- GLORYS's native grid is **−180..180** longitude while the model grid is 0–360. The regrid
  reconciles this internally; no manual fix needed. But if you sample `λnodes` of a GLORYS grid
  yourself, do **not** subtract 360.
- On a `LatitudeLongitudeGrid`, `xnodes`/`ynodes` return **metres** (x is 2-D: Δx = R·cosφ·Δλ). For
  degrees use `λnodes`/`φnodes`, which are 1-D.

## Caches

Both live under `~/.julia/scratchspaces/904d977b-046a-4731-8b86-9235c0d1ef02/`:

- `bathymetry_cache/bathymetry_{Nx}x{Ny}_{lon}_{lat}_{hash}.jld2` — keyed on grid size and extent, so
  changing the domain or `res_factor` regenerates it once (a few minutes, no network).
- `GLORYS/{var}_GLORYSDaily_{date}_{region}.nc` and `..._inpainted.jld2` — keyed on the **region**, so
  `glorys_region` in the model script is deliberately fixed at (235.25–238.25, 46.75–50.00) rather
  than derived from the domain edges. Tying it to the domain would invalidate the cache on every
  boundary tweak and force a re-download.

## Run cost

~12 s/iter at 250 × 177 × 36 on CPU, i.e. **~1.6 h wall per simulated day** once Δt settles near its
CFL limit of ~3.4 min. A 7-day run is ~2580 iterations and ~12 h. The `pacing:` line in the progress
output reports live s/iter and an ETA — use s/iter for planning, since the ETA reads pessimistic
until Δt caps.

## GCP GPU runs

Verified 2026-09-15 on project `<gcp-project>`. No Docker: the VM runs Julia directly.

**Image.** Family `canary-julia112-cuda` (`canary-julia112-cuda-v1`), built on the ModelBenchmarking
`climaocean-env` image (NVIDIA driver 535, Container Toolkit, gcsfuse) plus Julia 1.12.2 and this
repo's depot, precompiled on an L4 so the CUDA-aware cache is real. CUDA.jl fetches its own runtime
(13.4), so the old driver is fine. v1 accidentally includes ~23 GB of local data copies and gcsfuse
cache (a cleanup step failed before the snapshot); rebuild from a cleaned disk if storage matters.

**Machine.** `g2-standard-8` (L4 24 GB, 8 vCPU, 31 GB RAM) in `us-central1-a`. Quota at time of writing:
20 L4, 16 A100, 8 T4 in us-central1, none in use.

**Data.** Mounted with gcsfuse rather than copied; both buckets are US multi-region.

    gcsfuse --implicit-dirs <obs-bucket> /mnt/obs
    gcsfuse --implicit-dirs --cache-dir /var/gcsfuse-cache --file-cache-max-size-mb -1 \
            --file-cache-cache-file-for-range-read <forcing-bucket> /mnt/forcing-cache
    CANARY_GLORYS_DIR=/mnt/obs/glorys  CANARY_JRA55_DIR=/mnt/forcing-cache/jra55-ryf

Canary GLORYS (both windows, 288 MB) was uploaded to `<obs-bucket>/glorys/`; JRA55 RYF was
already in `<forcing-bucket>/jra55-ryf/`. Measured: on first load mount vs local copy
differs only by noise (compile time dominates); on JRA55 window reloads a plain mount is 17% slower
(1.17 vs 1.00 s per reload) and a mount with file cache matches a copy (0.97 s, 9 GB of local cache).

**Traps, all hit:**

- **Memory during GPU compilation can exhaust the VM**, which then thrashes rather than crashes — SSH
  times out and the kernel logs nothing. Diagnose with Cloud Monitoring (CPU flat, disk reads in the
  hundreds of thousands per minute), reset the VM, and run with a watchdog that kills Julia below 2 GB
  free. `kill -USR1 <julia pid>` prints a stack profile, which is how the compile blowup was located.
- **Julia's stderr is block-buffered when redirected**, so a log can stay empty for many minutes.
- **An illegal CUDA address is reported at the next synchronisation, not where it happens.**
  `CUDA_LAUNCH_BLOCKING=1` makes the error name the kernel that caused it.
- **Frequent SSH polling gets throttled** ("timed out during banner exchange"); poll every few minutes.
- `pgrep -f <pattern>` inside a `bash -c "... <pattern> ..."` matches itself.
- A full GPU launch spends ~8 minutes compiling before the first iteration at 0.1°.
