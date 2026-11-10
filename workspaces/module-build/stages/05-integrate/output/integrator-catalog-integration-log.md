# Integration log: integrator-catalog

Milestone M5 (third M5 module, D111), on `main`. Stage 05 of `workspaces/module-build`.

## Catalog rows flipped

| File | Row or section | Change |
|------|----------------|--------|
| `_config/catalogs/integrators.md` | S1, S3, S4, S6, S7, S9, S13, P2, P3, R1 | planned -> verified |
| `_config/catalogs/integrators.md` | S2, S8, S12, P1, G1 | planned -> verified (GPU runs of criterion 18 not run; P1: for SABA3) |
| `_config/catalogs/integrators.md` | S5 | planned -> verified (independent read of ForestRuth1990 Eqs. (4.3), (4.4) not run; GPU runs of criterion 18 not run) |
| `_config/catalogs/integrators.md` | S11 | planned -> verified (criterion 3(b2), McLachlan1995 Table 2 against Yoshida1990 Table 1, not run) |
| `_config/catalogs/integrators.md` | G2 | planned -> verified (Newton option, closed-form coefficients and quadratic invariants deferred; GPU runs of criterion 18 not run) |
| `_config/catalogs/integrators.md` | G3 | planned -> verified (tangent in the extended space not scored, open point P11) |
| `_config/catalogs/integrators.md` | source cells | S1, S2: Yoshida1990, Ruth1983 (HairerLubichWanner2006 not read); S5: ForestRuth1990, Yoshida1990, Omelyan2002 (CandyRozmus1991 not read); S7: KahanLi1997, McLachlan1995 (Suzuki1990 not read, spec Inconsistencies (c)); S9: CalvoSanzSerna1993a Eq. (2.8) (spec P5); G1, G2: Hairer2008 (HairerLubichWanner2006 not read) |
| `_config/catalogs/integrators.md` | notes cells | S7 orders and piece counts; S11 solutions built; S12 s9odr6a, s17odr8a; P1 n = 1 to 10; P2 also SBABC; P3 the tables built and not built (P8); G2 compensated sums (D113); G3 the exponent 1/(l - 1) (spec Inconsistencies (d)); R1 tangents by finite differences, package extension |
| `_config/catalogs/integrators.md` | S10, S14, G4, R2, R3, N1 to N10 | unchanged (planned; S10 and R2 deferred for want of read sources, the rest out of scope) |
| `_config/catalogs/systems.md` | F1, F4 | planned -> verified (`HarmonicOscillator`; `Pendulum` with the near-integrable form, near-separatrix FTLE not built) |
| `_config/catalogs/systems.md` | F2 | planned -> verified (`Kepler`, Hamiltonian and first integrals; the exact solution by Kepler's equation deferred) |

Evidence: `stages/04-verify/output/integrator-catalog-verification.md`, accepted at the stage 04 checkpoint (owner: "Accept, commit", 2026-11-09): full-size suite 248 241 assertions, 0 fail, 0 error, 0 broken; oracles O6 (OrdinaryDiffEqSymplecticRK 2.2.3), O7 (TaylorIntegration), O8 (REBOUND arrays, ABA864 b4 to 2.2e-39), O9 (OrdinaryDiffEqSDIRK ImplicitMidpoint) pass; criterion 15(d) passes under Amendment D114; independent reviews without CRITICAL or HIGH finding.

Open: 9 LOW findings of the first review (R4 to R7, R9 to R13) and of the D113 review 1 MEDIUM and 3 LOW (D1 resolved by D114; D2 `isolated_force_evaluations` has no method for `CorrectedMethod`, D3 the cyclic count checks only stage 1, D4 the Double64 compensated sum is not error-free). Not run: the GPU device runs of criterion 18 and O11 (`verify/gpu.jl` writes `data/gpu-benchmark.toml`; owner decision), criterion 3(b2) (McLachlan1995 text not in the repository), the independent read of ForestRuth1990 for criterion 5.

## Stage 05 changes

- Docs: new page `docs/src/api/integrators.md`, `Chaoslab.Integrators` added to the layers of `docs/make.jl` and to `docs/src/api.md`; `docs/src/index.md` layer table (Integrators row, Systems row names the flows, the Indicators row no longer says only SALI); `docs/src/api/systems.md` names the flows. Closes the CI documentation failure on the integrator-catalog commits (34 docstrings not in the manual).
- `src/Integrators/stepper.jl`: the `tangent_advance` docstring was attached to the private `_tangent_drift` (a comment and two definitions sat between it and the function), so the exported `tangent_advance` had no docstring and the reference in the `IntegratorMap` docstring did not resolve; moved above the function. No code change.
- `src/Integrators/CONTEXT.md`: the GPU sentence states the device runs were not run at stage 04.
- `CHANGELOG.md`: Unreleased, Added, `Chaoslab.Integrators` entry.
- Local checks: `julia --project=docs docs/make.jl` exit 0, no error; quality and layering 14 of 14 pass (Aqua, JET); the integrator-catalog folder suite 7317 pass, 0 fail, 2 broken (the Taylor extension skips, TaylorIntegration not loaded). The whole package suite was not rerun locally (CI runs it on push).

## Commits

- 06f8e6f module-build integrator-catalog stage 01: spec approved (D111)
- 0cca869 module-build integrator-catalog stage 02: RED tests and test plan
- e2b278c module-build integrator-catalog Amendment D112: test-side assertions corrected
- 2e35b3f module-build integrator-catalog stage 03: Integrators module, flows, Taylor extension
- 5f8f211 module-build integrator-catalog fix round D113
- 9d58dfe module-build integrator-catalog stage 04: verified after D113 and Amendment D114
- (this commit) module-build integrator-catalog stage 05: catalog rows verified, docs, changelog, integration log

## Files changed by the module

Short hashes of the integrator-catalog commits that touched each file. Stage paths are relative to `workspaces/module-build/`.

- `CHANGELOG.md` (stage 05)
- `Project.toml` (2e35b3f)
- `_config/catalogs/integrators.md`, `_config/catalogs/systems.md` (stage 05)
- `_config/decision-log.md` (06f8e6f, e2b278c, 5f8f211, 9d58dfe)
- `_config/risks.md` (5f8f211)
- `bench/integrator-catalog.jl` (2e35b3f)
- `docs/make.jl`, `docs/src/api.md`, `docs/src/index.md`, `docs/src/api/systems.md`, `docs/src/api/integrators.md` (stage 05)
- `ext/ChaoslabTaylorExt.jl` (2e35b3f); `ext/ChaoslabCUDAExt.jl` (5f8f211)
- `src/Chaoslab.jl` (2e35b3f)
- `src/Integrators/` `Integrators.jl`, `ensemble.jl`, `gauss.jl`, `splitting.jl`, `tao.jl` (2e35b3f, 5f8f211); `stepper.jl` (2e35b3f, 5f8f211, stage 05); `CONTEXT.md` (2e35b3f, 5f8f211, stage 05); `integrator_map.jl`, `tables.jl`, `tables_near.jl` (2e35b3f)
- `src/Systems/` `Systems.jl`, `flow_models.jl`, `CONTEXT.md` (2e35b3f); `flow_core.jl` (2e35b3f, 5f8f211)
- `test/integrator-catalog/` `runtests.jl`, `ic_ac17_23.jl` (0cca869, 5f8f211 for `ic_ac17_23.jl`); `ic_ac01_03.jl`, `ic_ac04_08.jl`, `ic_flows.jl` (0cca869, e2b278c); `ic_ac09_12.jl`, `ic_ac13_16.jl`, `ic_helpers.jl` (0cca869, e2b278c, 5f8f211)
- `test/analysis/an_p2_figures.jl`, `test/analysis/an_p3_dispatch.jl` (e2b278c: the TaylorIntegration weak dependency admitted)
- `verify/integrator-catalog.jl` (5f8f211, 9d58dfe); `verify/Project.toml`, `verify/Manifest.toml` (5f8f211)
- `01-spec/output/integrator-catalog-spec.md` (06f8e6f, e2b278c, 5f8f211, 9d58dfe)
- `02-tests/output/integrator-catalog-test-plan.md` (0cca869, e2b278c, 5f8f211, 9d58dfe)
- `03-implement/output/integrator-catalog-impl-notes.md` (2e35b3f, 5f8f211)
- `04-verify/output/` `integrator-catalog-verification.md`, `integrator-catalog-verify-log.txt`, `integrator-catalog-order.png`, `integrator-catalog-brouwer.png` (9d58dfe)
- `05-integrate/output/integrator-catalog-integration-log.md` (this commit)

## Audit

| Check | Result |
|---|---|
| Status matches evidence | every flipped row is backed by criteria that pass in the report; each part not run or deferred (GPU device runs, 3(b2), ForestRuth1990 read, G2 deferred parts, G3 tangent, F2 exact solution) is named in the row's qualifier |
| Log | commits and files listed above; the stage 05 SHA is the commit that adds this file |
| Changelog | `CHANGELOG.md` Unreleased, Added: `Chaoslab.Integrators` entry |

Not changed: `_config/architecture.md` (still lists `AbstractFlow` with dimension <= 8 and puts Integrators below Tangent; spec Inconsistencies (a)); `src/precompile.jl` has no integrator workload.
