# Verification: analysis

Module-build stage 04 for the slug `analysis` (decision D60). Branch `module/analysis`, working tree (nothing committed before stage 05). Verify scripts: `verify/analysis.jl` (checks, run under `julia --project=verify -t 23`) and `verify/analysis-independent.jl` (independent recomputation, 444 lines, written by a fresh Sonnet subagent from the spec and the preregistration without access to `src/Analysis`, the tests or the notes). Real inputs: manifest of commit e35ee5a of branch `study/p1-le-benchmark-run` (1363 of 1363 chunks ok), the store `data/p1-le-benchmark/` (1363 files, 475 MB), panel and preregistration of `main`, the new parameter file `workspaces/study/stages/04-analyze/references/p1-le-benchmark-analysis-parameters.toml` (P1).

Raw output: `analysis-verify-log.txt` (final run, 2026-10-20 15:50, 23 threads: PASS 40, FAIL 0, INFO 2), `analysis-independent-results.txt` (the independent script's results, 1021 lines, read by the comparison), plot `analysis-k0-closed-form.png`.

Run history, for the record. Run 1 stopped on a defect of the module that only the real store shows (items of the tau sweep with a breakdown hold NaN rows and were refused; fixed, see finding 2). Run 2 passed every check of the run itself and showed two defects of the script (field indices of the independent file, integer parsing) and, in the comparison, 15 orbits whose ladder test differed from the independent value (finding 4). Run 3, after the fixes of the review and the common half point in both implementations, is the log. Between run 3 and the end of the stage the only change to `src/Analysis` was a type assertion for JET (`pkgdir(...)::String`).

Discipline. The hypothesis outcomes were visible to the agent from run 2 on (they are printed in the comparison). No threshold, definition or reading was changed in the light of a result: every change is the repair of a defect (findings 1 to 3, 5, review) or the choice of one common half point for the ladder, argued by the reviewer from the grids of the items before any comparison had run. This report states agreement and measurements only and does not interpret the results; that is study stage 04.

## Evidence per criterion

| Criterion | Evidence | Result |
|-----------|----------|--------|
| 1 Refusals (E17) | 17 malformed inputs of the approved suite plus failed chunk, missing chunk list, missing control scope, parameter ranges, each refused before any output with the cause named, output directory not created; `Pkg.test` 54916 of 54916 | PASS |
| 2 Reader | fixture round trip of all 55 chunks (every series, point, flag, attribute); on the real store the module's reader and the independent reader (HDF5.jl only) give equal exponents, flags and horizons for 883 of 883 orbits | PASS |
| 3 Real-store closed form (O1) | 40 class T orbits at K = 0, every recorded n of the 8 methods: largest deviation 1.67e-16 (L1, L3, L4, L4b), 2.22e-16 (L5, L5b, L6), 5.55e-17 (L9 against its singular-value form); bound 1e-9. The independent script: 2.22e-16 and 5.55e-17. Plot `analysis-k0-closed-form.png` | PASS (O1 corrected for L9, finding 1) |
| 4 H1 metrics, 5 H1 decision, 6 controls, 7 R and H2, 8 n_conv, 10 H3 ratio | approved unit tests (boundaries at and one float below) 434 of 434; on the real store Delta_QB, abs(L9 - L5), abs(L2 - L5), R(n_H2), n_conv and its defined flag equal the independent value in 883 of 883 orbits (exact or 1e-12 relative) | PASS |
| 7 sets, 12 labels | E_ch (209), W (355), S (80), dropped (1 orbit, drift) equal; raw labels, final labels, forced tests, ladder agreement, flags equal in 883 of 883 orbits | PASS |
| 9 Bootstrap | exact enumeration of 6^6 resamples, stream rule against SHA-256 of the independent tool, coverage 900 to 1000 of 1000 (unit tests); H2 and H3 interval ends of the module against the independent script's own generator: within 0.05 of the interval width (the check passes) | PASS |
| 11 H4, 13 H5 | unit tests; real store: the eight values of H4 and the counts of H5a, H5b, H5c equal the independent ones, verdict words equal | PASS |
| 14 Order and parameters | unit tests (six perturbations change only what depends on them; decisions in order H1 to H5) | PASS |
| 15 Output content | unit tests; the real tables file has the T1 to T7 and F1 to F7 sections, readings, deviations, header with hashes, module version and commit | PASS |
| 16 Determinism | real store: tables file at 23 threads and at 1 thread byte-identical (runs of 54.0 s and 39.8 s); fixture runs at 1 and 4 threads | PASS |
| 17 Read-only | SHA-256 of the 1367 input files (1363 chunks, manifest, panel, preregistration, parameters) equal before the run, after the run and after the second run | PASS |
| 18 Figures | real store: seven files `p1-le-benchmark-fig-F1.pdf` to `F7.pdf`, sizes 106.0, 64.8, 9.5, 7.7, 11.2, 6.8 and 13.9 kB, each starting `%PDF-`; the digest lines of the tables file equal the digests of the figure data | PASS |
| 19 Independent recomputation (O3) | `verify/analysis-independent.jl` against the module: 883 orbits; 13 per-orbit quantities 0 differences; sets, dropped, S without positive growth equal; H1 counts, H2 and H3 (sizes, medians, ratio exact; interval ends within 0.05 of width), H4, H5 equal; 47 label-count entries 0 differing | PASS |
| 20 Layering and dependencies | `Pkg.test`: layering 2 of 2, Aqua 11 of 11, JET clean (after one type assertion); no plotting package in `src/Analysis`; Figures is a package extension with the weak dependency CairoMakie | PASS |
| 21 Performance | `analyze` with figures on the real store: 54.0 s wall at 23 threads (first call, compilation included), peak resident memory 2.45 GiB, 17.1 GiB allocated, 2.3 s in the garbage collector; limits 600 s and 8 GiB | PASS |
| 22 Command | `bin/chaoslab.jl analyze` on the real inputs under the verify environment: exit code 0, tables file and seven figure files written; unknown option: exit code 1 with the usage text; without CairoMakie exit code 2 (unit test) | PASS |

Package suite at the end of the stage: `Pkg.test` 54916 pass, 0 fail, 0 error in 6 min 43 s (analysis 434, experiments 4861, the earlier modules unchanged); coverage of `src/Analysis` 98.3 percent (698 of 710 executable lines).

## Findings of the verification

1. **Spec error, oracle O1 for L9.** The closed form of the spec holds for L1 and the QR family. L9 stores ln(sigma_max)/n of the product, which for the integrable map is ln[(n + sqrt(n^2 + 4))/2]/n; against the form of the spec it differs by 0.0416 (n = 3). Found by the independent author, confirmed on the real store (log line "finding"). Spec corrected ("Corrections found in stage 04", item 1); no change of the module.
2. **Defect, breakdown items refused.** Items of the tau sweep with the flag `breakdown_tau` store `t = 0` and NaN after the breakdown. The finiteness check refused the first one (chunk 502, orbit 791). The reviewer found the same independently. Fixed (items skipped by the check, listed in F7 with step and kind), regression test `an_breakdown.jl`.
3. **Defect, Y5 control without `K`.** The real panel's control 883 has no `K` and no `K_text`; the module raised a `KeyError`. Found by the reviewer. Fixed, the fixture now omits them.
4. **Reading, ladder half point.** The Float64 and Double64 items (horizon N_1) have no step N_B/2, the BigFloat item (horizon N_B) has; comparing windows of different length made 15 S_L orbits differ in the ladder test between the two implementations. Decided: one common half point for the three precisions. Both implementations changed; then 0 differences.
5. **Typing (JET).** An untyped `*`, `join` and a generic function without methods were reported by JET; replaced by annotations and a function reference for the extension.
6. **Data note, orbit 881.** The independent script reports that its drift flag is not reproduced from the recorded points of its lifted L5 items (largest |p| 4.8e4 against the bound 1e6). The module's own check of every item it reads (the drift flag against the recomputation from the stored points, the runner's rule with windings) passed on the real store, so the flag stands on other items of the orbit. No action.

## Independent review (fresh Sonnet subagent, spec and code, not the notes)

Read-only review of `src/Analysis/*.jl` and the extension against the spec. 2 CRITICAL, 2 HIGH, 6 MEDIUM, 12 LOW. Dispositions:

| Severity | Finding | Disposition |
|----------|---------|-------------|
| CRITICAL | `K` read for the Y5 control that has none | fixed (finding 3) |
| CRITICAL | `breakdown_tau` items refuse the run | fixed (finding 2) |
| HIGH | a chunk that failed twice crashes the run in `orbit_record` and enters F1, F5, F6 | closed by a named refusal of any manifest chunk with a status other than ok (the study has none: 1363 of 1363 ok); full handling of failed chunks is not built, recorded in the spec corrections, item 3 |
| HIGH | ladder compares labels of different half points | fixed (finding 4) |
| MEDIUM | diagnostics flag read only on the N_1 item | now read on the N_1 and the label-horizon L5 item |
| MEDIUM | `control_stop` chunk accepted | covered by the refusal of non-ok chunks |
| MEDIUM | silent fallbacks (no long item, fewer than two directions, no BigFloat item) | refused with the orbit named |
| MEDIUM | control with a positive exponent and no scope skipped | refused |
| MEDIUM | no final re-hash of the inputs | added before the tables are written |
| MEDIUM | breakdown step and kind not tabulated, missing keys with only breakdowns | tabulated under F7, keys kept |
| LOW | lambda_ref <= 0 in F2; F2 not pooled | orbits left out; pooled series over W added |
| LOW | header lacks module version and commit | added |
| LOW | panel start hash and chunk attributes not compared | added |
| LOW | parameters not range checked | added |
| LOW | per-orbit rebuild of two maps | hoisted |
| LOW | F6 mean is an unstated choice | added to the readings |
| LOW | figure drawer exception, empty F7 | status `:figures_failed`, exit code 2; guard added |
| LOW | malformed manifest raises KeyError | named refusal |
| LOW | strong-chaos K constant outside E0 | kept by decision of the owner (spec corrections, item 7) |
| LOW | tie at the H4 tolerance not exact in binary64 | not changed: a tie has probability zero on measured values; the approved test uses the computed distance |
| LOW | percentile level (1 - 0.95)/2 is 0.025000000000000022 | not changed: harmless at 1e4 resamples |

Regression tests for the fixed findings: `an_breakdown.jl`, `an_review.jl` (listed in the test plan, "Tests added in stage 04").

## Audit

| Check | Result |
|-------|--------|
| Evidence: every acceptance criterion has a number or a plot | pass |
| Review closed: no open CRITICAL or HIGH finding | pass (the HIGH on failed chunks is closed by a refusal, a limit stated in the spec corrections; the owner may ask for full handling) |

## Open for the owner and for stage 05

- The K values of S are a constant of the module (kept by decision). A generic form needs a spec amendment and stage 02 changes.
- Chunks that failed twice are refused, not analysed. The study has none.
- The experiments suite took 3 min 53 s in the last full run, over its 3 minute ceiling (D35), as in stage 05 of `experiments`; the analysis suite adds 48 s and is under it.
- The final manifest `p1-le-benchmark-manifest.json` is on branch `study/p1-le-benchmark-run`, not on `main`; study stage 04 needs it on the branch it runs from. The parameter file was added to the study stage 04 folder on this branch and belongs to that stage.
- `verify/Manifest.toml` gained the extension entry; `Manifest.toml` for Julia 1.10 is still unresolved (earlier open item).
- Study stage 04 (running `analyze`, drafting the findings) has not been started.
