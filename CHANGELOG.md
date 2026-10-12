# Changelog

All notable changes to `Chaoslab` and this repository. Versions follow Semantic Versioning (`_config/build-order.md`, section "Releases").

## Unreleased

### Added (M0 bootstrap)

- ICM workspace: four workspaces, 19 stage contracts, `_config/` references, run-store schema v1, setup questionnaire.
- Julia package skeleton `Chaoslab` (empty module), test entry with Aqua, JET and a layering test that also tests itself.
- ICM scripts (Python standard library): `icm_validate.py`, `icm_status.py`, `bibcheck.py`, `promote.py`, with offline unit tests.
- CI on Windows and Ubuntu, Julia LTS and current release, two threads; a separate job runs the ICM checks and script tests.
- Compute policy: GPU first with automatic CPU fallback, one kernel for both, allocation-free hot loops, thread-count invariance, backend recorded per run (decision D17).
- Native Windows build of REBOUNDx 5.1.0 and celmech 1.5.8 (`python/nbody/windows-build/`): patches, build script, smoke tests; closes risk K2 (decision D18).

### Added (M1 literature)

- Four literature runs (lyapunov-algorithms, standard-map-stickiness, chaos-indicators, symplectic-integrators): bibliography 180 to 353 verified entries, 40 located benchmark values, sticky-orbit candidates for the P1 panel (decision D19).
- `bibcheck.py`: DataCite fallback, old-style arXiv IDs, duplicate records under a second DOI form, key suffixes per `conventions.md`, proceedings years, em dashes as `---`. `promote.py`: blocks duplicate titles; flips only keys the topic adds.

### Deferred, with reason

- `Manifest.toml`: nothing to pin until the first real dependency (M2); then one manifest per Julia version in CI (decision D13).
- `bin/chaoslab.jl`, `bench/`, `python/nbody/`: no caller before M2 or M3.
- Coverage reporting: nothing to cover before the first module (M2).
- `CITATION.cff` and Zenodo: after `setup` fills the author data, before the first tag `v0.1.0` (M3).
