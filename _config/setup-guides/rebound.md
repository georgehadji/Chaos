# Setup: REBOUND ecosystem (N-body, P6 only)

REBOUND is a C library with Python bindings for N-body integration. Only study/03 runs for N-body studies need it. Facts below come from the tooling check of 2026-10-12 (workflow map Section 8).

## Install (venv A: REBOUND, REBOUNDx, celmech)

1. `uv venv -p 3.12 .venv-nbody` in `python/nbody/` (see `python-and-git.md`).
2. `uv pip install -p .venv-nbody/Scripts/python.exe rebound` (binary wheels).
3. REBOUNDx and celmech: on Linux and macOS, `pip install reboundx celmech`. On native Windows neither builds from PyPI; run `python python/nbody/windows-build/build_windows.py --venv python/nbody/.venv-nbody` instead (section "Native Windows build").

## Native Windows build

Verified 2026-10-12 on Windows 11 x64, CPython 3.12, REBOUND 5.2.2: REBOUNDx 5.1.0 and celmech 1.5.8 build, install and pass `smoke_tests.py` (12 checks) and celmech's 24 upstream tests. Evidence: `python/nbody/windows-build/verification.txt`.

- Needs: MSVC Build Tools 2022, `uv` and `git` on PATH. The script downloads the sdists from PyPI (sha256 pinned), applies the patches in `patches/`, builds both wheels and installs them; nothing outside the venv and its work folder.
- REBOUNDx builds with MSVC (VLA and `rand_r` fixes, import library generated from REBOUND's DLL). celmech's C library builds with `zig cc` from the `ziglang` PyPI package, because MSVC's C compiler lacks C99 complex arithmetic.
- `celmech-1.5.8-modern-deps-compat.diff` fixes celmech 1.5.8 against REBOUND 5, SciPy 1.18 and Python 3 on every platform (`N_real`, `xrange`, `vode` shape); without it 15 of 24 upstream tests fail.
- Limits: the REBOUNDx wheel is tied to the REBOUND version it was linked against (rebuild after upgrading REBOUND); `stochastic_forces` uses a portable `rand_r`, so stochastic runs are reproducible on Windows but not bit-identical to Linux; only the REBOUNDx effects in `smoke_tests.py` were run.

## Install (venv B: ASSIST)

ASSIST pins `rebound<5` and supports Linux and macOS; on Windows use WSL2. It also needs JPL ephemeris files (about 1 GB for DE440 plus the small-body file; larger kernels 2-14 GB) under the data root.

## Verify

- `python -c "import rebound; print(rebound.__version__)"`
- For bit-reproducible runs, build REBOUND from source with `FFP_CONTRACT_OFF`; record version and build flags in the manifest.

## How this repo uses it

Adapters in `python/nbody/` write HDF5 per `/schema/run-store.md`; Julia reads those files. No in-process bridge.
