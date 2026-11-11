# Integration log: henon-heiles

Milestone M5 (catalog row F3, D116), on `main`. Stage 05 of `workspaces/module-build`.

## Catalog rows flipped

| File | Row or section | Change |
|------|----------------|--------|
| `_config/catalogs/systems.md` | F3 | planned -> verified (GPU part of criterion 16 not run; portrait at 200 x 200, the paper's 500 x 500 not run) |
| `_config/catalogs/systems.md` | F3 source cell | HenonHeiles1964 -> AntonopoulosManosSkokos2010 Sect. 3.1, Hairer2008 Sect. 2 (HenonHeiles1964 not read; spec Inconsistencies (e)) |
| `_config/catalogs/systems.md` | F3 known-answer cell | the portraits stated as data (fractions, SALI at T_f = 2000; spec Inconsistencies (j)) |
| `_config/bibliography/indicators.md` | AntonopoulosManosSkokos2010 | "not yet read" -> Sect. 3.1 read (spec Inconsistencies (d)) |

Evidence: `stages/04-verify/output/henon-heiles-verification.md`, accepted at the stage 04 checkpoint (owner: "Accept, commit", 2026-11-11): 17 of 17 criteria pass after Amendment D117; oracles O1 to O6 (256-bit closed forms, ForwardDiff and complex step, TaylorIntegration BigFloat references, ChaosTools spectrum, AntonopoulosManosSkokos2010 fraction 0.546 at E = 1/8, escape continuation); bitwise 1 against N CPU threads; independent review without CRITICAL or HIGH finding.

Open: review findings M2 (section Jacobian restarted on the sub-step path, within the order-8 tolerance) and L3, L5 to L10 (LOW). Not run: the GPU part of criterion 16 (`verify/gpu.jl` writes `data/gpu-benchmark.toml`; owner decision), the 500 x 500 portrait (preregistered study size).

## Stage 05 changes

- `_config/catalogs/systems.md` row F3, `_config/bibliography/indicators.md` entry AntonopoulosManosSkokos2010 (above).
- `CHANGELOG.md`: Unreleased, Added, the Henon-Heiles entry.
- `src/Systems/CONTEXT.md`, `src/Integrators/CONTEXT.md`: already name F3 and the section and escape tools (stage 03); unchanged here.
- Local checks: `julia --project=docs docs/make.jl` exit 0, no error or warning; quality and layering 14 of 14 pass after D117; henon-heiles suite 6904 pass; whole package suite 1426407 pass, 0 fail, 3 broken at stage 03 (not rerun after the D117 fixes, CI runs it on push).

## Commits

- d011607 module-build henon-heiles stage 01: spec approved (D116)
- d4d37b4 module-build henon-heiles stage 02: RED tests and test plan
- 4c83095 module-build henon-heiles stage 03: Systems energy-surface tools, Integrators section and escape (GREEN)
- 8bd0625 module-build henon-heiles stage 04: verified after Amendment D117
- (this commit) module-build henon-heiles stage 05: catalog row F3 verified, changelog, integration log

## Files changed by the module

Short hashes of the henon-heiles commits that touched each file. Stage paths are relative to `workspaces/module-build/`.

- `CHANGELOG.md` (stage 05)
- `_config/catalogs/systems.md`, `_config/bibliography/indicators.md` (stage 05)
- `_config/bibliography/benchmarks.md` (d011607)
- `_config/decision-log.md` (d011607, 8bd0625)
- `bench/henon-heiles.jl` (4c83095)
- `docs/src/api/systems.md`, `docs/src/api/integrators.md` (4c83095)
- `src/Systems/` `henon_heiles.jl` (4c83095, 8bd0625); `flow_models.jl`, `Systems.jl`, `CONTEXT.md` (4c83095)
- `src/Integrators/` `section.jl`, `escape.jl` (4c83095, 8bd0625); `Integrators.jl`, `CONTEXT.md` (4c83095)
- `test/henon-heiles/` `runtests.jl`, `hh_helpers.jl`, `hh_ac01_05.jl`, `hh_ac13_17.jl` (d4d37b4); `hh_ac06_09.jl`, `hh_ac10_12.jl` (d4d37b4, 8bd0625)
- `verify/henon-heiles.jl` (8bd0625)
- `01-spec/output/henon-heiles-spec.md` (d011607, d4d37b4, 8bd0625)
- `02-tests/output/henon-heiles-test-plan.md` (d4d37b4, 8bd0625)
- `03-implement/output/henon-heiles-impl-notes.md` (4c83095)
- `04-verify/output/` `henon-heiles-verification.md`, `henon-heiles-verify-log.txt`, `henon-heiles-d117-verify-log.txt`, `henon-heiles-energy.png`, `henon-heiles-portrait.png` (8bd0625)
- `05-integrate/output/henon-heiles-integration-log.md` (this commit)

## Audit

| Check | Result |
|---|---|
| Status matches evidence | F3 verified on 17 of 17 criteria; the GPU part of criterion 16 and the 500 x 500 portrait named in the row's qualifier |
| Log | commits and files listed above; the stage 05 SHA is the commit that adds this file |
| Changelog | `CHANGELOG.md` Unreleased, Added: Henon-Heiles entry |

Not changed: `_config/architecture.md` (Module map: systems with `step`, `tangent`, `invariants`; spec Inconsistencies (f)); `src/precompile.jl` has no henon-heiles workload.
