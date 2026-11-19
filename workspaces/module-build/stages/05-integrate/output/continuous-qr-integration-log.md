# Integration log: continuous-qr

Milestone M5 (catalog row L8, D120, D121), on `main`. Stage 05 of `workspaces/module-build`.

## Catalog rows flipped

| File | Row or section | Change |
|------|----------------|--------|
| `_config/catalogs/lyapunov.md` | L8 | planned -> verified (CPU only: no device kernel, criterion 22 device part not applicable); method cell restated with A for the Jacobian, the k < d form and the normalised rates (spec Inconsistencies (a), D121 M2); "Gives" per unit time and the first k exponents; main-issue cell restated from the evidence (error set by the Q-scheme order, penalty term, R1 time ratio against Geist1990's ~40x) |
| `_config/catalogs/lyapunov.md` | L8 source cell | Goldhirsch1987, Dieci1997, ChristiansenRugh1997 -> Skokos2010 (Sect. 6.4), the three named as not on disk (spec Sources) |
| `_config/bibliography/lyapunov.md` | Goldhirsch1987, Dieci1997, ChristiansenRugh1997 | "not on disk" recorded (spec Sources) |
| `_config/bibliography/benchmarks.md` | row "Mean relative CPU time of Lyapunov-exponent methods" | R1 (0.36 of the discrete RGS time on F3 orbit C) noted as not a reproduction |

Evidence: `stages/04-verify/output/continuous-qr-verification.md`, accepted at the stage 04 checkpoint (owner: "Accept", 2026-11-19): 22 of 22 criteria pass after Amendment D121; oracles BigFloat(256) closed forms, TaylorIntegration BigFloat(256) variational references, the discrete QR spectrum, ChaosTools 3.5.4; bitwise 1 against 8 CPU threads; independent review without CRITICAL or HIGH finding; M1 to M4 handled under D121 (M2 kept over the exact trace cancellation, owner: "Keep M2, amend").

Open: review findings L1 to L6 (LOW). Not applicable: the device part of criterion 22 (no kernel in this release).

## Stage 05 changes

- `_config/catalogs/lyapunov.md` row L8, `_config/bibliography/lyapunov.md`, `_config/bibliography/benchmarks.md` (above).
- `CHANGELOG.md`: Unreleased, Added, the continuous QR entry.
- `src/Lyapunov/CONTEXT.md`: the held-items row for continuous-qr now names the missing device kernel and `verify/continuous-qr.jl` (75 lines).
- Local checks: `julia --project=docs docs/make.jl` (result in the stage 05 commit message); after D121 the package continuous-qr suite 1625 pass, Aqua, JET and layering 14 pass, FULL1 to FULL4 680, 88, 67536, 107 pass; whole package suite recorded at stage 03 (not rerun after D121; CI runs it on push).

## Commits

- 5b6a036 module-build continuous-qr stage 01: spec approved (D120)
- 165c3ed module-build continuous-qr stage 02: RED tests and test plan
- 8850000 module-build continuous-qr stage 03: continuous QR spectrum in Lyapunov (GREEN)
- 352943d module-build continuous-qr stage 04: accepted after Amendment D121
- (this commit) module-build continuous-qr stage 05: catalog row L8 verified, changelog, integration log

## Files changed by the module

Short hashes of the continuous-qr commits that touched each file. Stage paths are relative to `workspaces/module-build/`.

- `CHANGELOG.md` (stage 05)
- `_config/catalogs/lyapunov.md`, `_config/bibliography/lyapunov.md`, `_config/bibliography/benchmarks.md` (stage 05)
- `_config/decision-log.md` (5b6a036, 352943d)
- `bench/continuous-qr.jl` (8850000)
- `docs/src/api/lyapunov.md` (8850000)
- `src/Lyapunov/` `continuous_qr.jl` (8850000, 352943d), `Lyapunov.jl` (8850000), `CONTEXT.md` (8850000, 352943d, stage 05)
- `test/continuous-qr/` `runtests.jl` (165c3ed); `cq_helpers.jl`, `cq_ac01_05.jl`, `cq_ac06_10.jl`, `cq_ac11_15.jl`, `cq_ac16_22.jl` (165c3ed, 352943d)
- `verify/continuous-qr.jl` (352943d)
- `01-spec/output/continuous-qr-spec.md` (5b6a036, 352943d)
- `02-tests/output/continuous-qr-test-plan.md` (165c3ed, 352943d)
- `03-implement/output/continuous-qr-impl-notes.md` (8850000)
- `04-verify/output/` `continuous-qr-verification.md`, `continuous-qr-verify-log.txt`, `continuous-qr-d121-verify-log.txt`, `continuous-qr-convergence.png`, `continuous-qr-orthogonality.png` (352943d)
- `05-integrate/output/continuous-qr-integration-log.md` (this commit)

## Audit

| Check | Result |
|---|---|
| Status matches evidence | L8 verified on 22 of 22 criteria; the CPU-only qualifier names the absent device part |
| Log | commits and files listed above; the stage 05 SHA is the commit that adds this file |
| Changelog | `CHANGELOG.md` Unreleased, Added: continuous QR entry |

Not changed: `_config/build-order.md` row M5 P3 (status by milestone); `src/precompile.jl` has no continuous-qr workload; `_config/risks.md` (reviewed at M5 exit, after clv-ginelli).
