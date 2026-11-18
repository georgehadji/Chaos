# Integration log: ofli

Milestone M5 (catalog row I6, D118, D119), on `main`. Stage 05 of `workspaces/module-build`.

## Catalog rows flipped

| File | Row or section | Change |
|------|----------------|--------|
| `_config/catalogs/indicators.md` | I6 | planned -> verified (GPU part of criterion 18 not run); vectors "1" -> "1 or full basis" (spec P8, Inconsistencies (a)); regular and chaotic cells and the main-issue cell restated from the evidence (bounded value on a periodic orbit lost after a time of order h^-p, R1; flow held by the reducer, P10; `ofli_total`, D119) |
| `_config/catalogs/indicators.md` | I6 source cell | Fouchard2002 -> Darriba2012, Maffione2013 (Fouchard2002 abstract only; spec Sources) |
| `_config/bibliography/indicators.md` | Fouchard2002, Darriba2012, Maffione2013 | read status and sections recorded (spec Sources) |
| `_config/bibliography/benchmarks.md` | new row "OFLI behaviour by orbit type (flows)" | Darriba2012 qualitative statements, no numbers (spec P8, F7) |

Evidence: `stages/04-verify/output/ofli-verification.md`, accepted at the stage 04 checkpoint (owner: "Accept", 2026-11-18): 18 of 18 criteria pass after Amendment D119; oracles BigFloat(256) closed forms, same-orbit BigFloat references with singular values of P Y, TaylorIntegration BigFloat references of F3, ChaosTools sums cross-check; bitwise 1 against 8 CPU threads; independent review without CRITICAL or HIGH finding; M1 to M3 fixed under D119.

Open: review findings L1 to L4 (LOW). Not run: the GPU part of criterion 18 (`verify/gpu.jl` writes `data/gpu-benchmark.toml`; owner decision).

## Stage 05 changes

- `_config/catalogs/indicators.md` row I6, `_config/bibliography/indicators.md`, `_config/bibliography/benchmarks.md` (above).
- `CHANGELOG.md`: Unreleased, Added, the OFLI entry.
- `src/Indicators/CONTEXT.md`: already names `OFLI`, `ofli_initial`, `ofli_total` (stages 03, 04); unchanged here.
- Local checks: `julia --project=docs docs/make.jl` exit 0, no error or warning; after D119 the package ofli suite 349 pass, Aqua, JET and layering 14 pass, FULL1 to FULL4 78, 841, 380, 22 pass; whole package suite 1426757 pass, 0 fail, 3 broken at stage 03 (not rerun after D119; CI runs it on push).

## Commits

- e870cdb module-build ofli stage 01: spec approved (D118)
- 73ddda1 module-build ofli stage 02: RED tests and test plan
- 7eaf9d4 module-build ofli stage 03: OFLI reducer and ofli_initial (GREEN)
- edcea52 module-build ofli stage 04: verification, Amendment D119 and fixes M1-M3
- 02e6b20 module-build ofli stage 04: FULL3 and FULL4 reruns after D119 pass
- 784db9f module-build ofli stage 04: accepted after Amendment D119
- (this commit) module-build ofli stage 05: catalog row I6 verified, changelog, integration log

## Files changed by the module

Short hashes of the ofli commits that touched each file. Stage paths are relative to `workspaces/module-build/`.

- `CHANGELOG.md` (stage 05)
- `_config/catalogs/indicators.md`, `_config/bibliography/indicators.md`, `_config/bibliography/benchmarks.md` (stage 05)
- `_config/decision-log.md` (e870cdb, edcea52)
- `bench/ofli.jl` (7eaf9d4)
- `docs/src/api/indicators.md` (7eaf9d4, edcea52)
- `src/Indicators/` `ofli.jl`, `Indicators.jl`, `CONTEXT.md` (7eaf9d4, edcea52)
- `test/ofli/` `runtests.jl`, `of_helpers.jl` (73ddda1); `of_ac01_05.jl`, `of_ac06_10.jl`, `of_ac11_15.jl`, `of_ac16_18.jl` (73ddda1, edcea52)
- `verify/ofli.jl` (edcea52)
- `01-spec/output/ofli-spec.md` (e870cdb, edcea52)
- `02-tests/output/ofli-test-plan.md` (73ddda1, edcea52)
- `03-implement/output/ofli-impl-notes.md` (7eaf9d4)
- `04-verify/output/` `ofli-verification.md` (edcea52, 02e6b20, 784db9f), `ofli-verify-log.txt`, `ofli-plateau.png`, `ofli-classes.png` (edcea52), `ofli-d119-verify-log.txt` (edcea52, 02e6b20)
- `05-integrate/output/ofli-integration-log.md` (this commit)

## Audit

| Check | Result |
|---|---|
| Status matches evidence | I6 verified on 18 of 18 criteria; the GPU part of criterion 18 named in the row's qualifier |
| Log | commits and files listed above; the stage 05 SHA is the commit that adds this file |
| Changelog | `CHANGELOG.md` Unreleased, Added: OFLI entry |

Not changed: `_config/build-order.md` row M5 P3 (names `ofli`, status by milestone); row I7 (OFLI2) stays catalog only (spec P7); `src/precompile.jl` has no ofli workload.
