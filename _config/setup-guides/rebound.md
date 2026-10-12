# Setup: REBOUND ecosystem (N-body, P6 only)

REBOUND is a C library with Python bindings for N-body integration. Only study/03 runs for N-body studies need it. Facts below come from the tooling check of 2026-10-12 (workflow map Section 8).

## Install (venv A: REBOUND, REBOUNDx, celmech)

1. `python -m venv .venv-nbody` in `python/nbody/`, then activate it.
2. `pip install rebound` (binary wheels). REBOUNDx ships as source only: on Windows it needs the MSVC Build Tools.
3. celmech is a C extension; its Windows build and REBOUND 5 compatibility are unverified. If it fails, use WSL2.

## Install (venv B: ASSIST)

ASSIST pins `rebound<5` and supports Linux and macOS; on Windows use WSL2. It also needs JPL ephemeris files (about 1 GB for DE440 plus the small-body file; larger kernels 2-14 GB) under the data root.

## Verify

- `python -c "import rebound; print(rebound.__version__)"`
- For bit-reproducible runs, build REBOUND from source with `FFP_CONTRACT_OFF`; record version and build flags in the manifest.

## How this repo uses it

Adapters in `python/nbody/` write HDF5 per `/schema/run-store.md`; Julia reads those files. No in-process bridge.
