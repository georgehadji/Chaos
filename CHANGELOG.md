# Changelog

All notable changes to `Chaoslab` and this repository. Versions follow Semantic Versioning (`_config/build-order.md`, section "Releases").

## Unreleased

### Added (M0 bootstrap)

- ICM workspace: four workspaces, 19 stage contracts, `_config/` references, run-store schema v1, setup questionnaire.
- Julia package skeleton `Chaoslab` (empty module), test entry with Aqua, JET and a layering test that also tests itself.
- ICM scripts (Python standard library): `icm_validate.py`, `icm_status.py`, `bibcheck.py`, `promote.py`, with offline unit tests.
- CI on Windows and Ubuntu, Julia LTS and current release, two threads; a separate job runs the ICM checks and script tests.
- Compute policy: GPU first with automatic CPU fallback, one kernel for both, allocation-free hot loops, thread-count invariance, backend recorded per run (decision D17).

### Deferred, with reason

- `Manifest.toml`: nothing to pin until the first real dependency (M2); then one manifest per Julia version in CI (decision D13).
- `bin/chaoslab.jl`, `bench/`, `python/nbody/`: no caller before M2 or M3.
- Coverage reporting: nothing to cover before the first module (M2).
- `CITATION.cff` and Zenodo: after `setup` fills the author data, before the first tag `v0.1.0` (M3).
