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


## Amendment D89 (stage 04)

Module-build stage 04 for Amendment D89 (study `p2-indicators`; decision D89). Module at commit b792724 (`Chaoslab` 0.1.2), working tree clean at the start of the stage; this stage adds only `verify/analysis-p2.jl` and the output files named here. Verify scripts: `verify/analysis-p2.jl` (module side: TRANSCRIPTION, RUN, REPORT, COMPARE, P1UNCH, CLI; run under `julia --project=verify -t 23`), `verify/analysis-p2-independent.jl` (independent reader and store-identity checker, 617 lines, HDF5, SHA and TOML only, no Chaoslab code) and `verify/analysis-p2-decision.py` (independent decision layer, 1081 lines, numpy and scipy, own bootstrap generator). Per the headers of the two independent files they were written from the spec and the frozen preregistration without reading `src/Analysis`.

Raw output: `analysis-p2-verify-log.txt` (final run, 2026-10-29 23:37, 23 threads: PASS 133, FAIL 4, INFO 88; the four FAIL lines are the comparison differences of finding 1 and finding 2 below, left in the log as they came out), `analysis-p2-independent-results.txt` (14180 lines, SHA-256 ef2f1714be4b7614ddd07635037052c2718162b4ccc664a8df8f1a793af513a7, read by the comparison). Nothing was written under `data/` (listing of the data tree, 3374 files: SHA-256 fb9f38c1...5e20a before and after). No export file was copied into the repository.

Run history, for the record.
- Independent side (run before this stage's script, results file as above): the reader read 695 chunk files with 0 anomalies and the decision layer passed its 16 self checks. It reported no store defect.
- Module side, TRANSCRIPTION: first run, 4 FAIL, all defects of my script (the list of machine keys was sorted as a tuple; the four label keys looked up the P1 sections by their short name while the P1 headings carry a suffix, "H1 (claim C1, ...)"). Script repaired, table untouched; second run 0 failures.
- Module side, development runs (not the evidence): RUN alone (53.5 s, 1.89 GiB, all checks PASS) and REPORT plus COMPARE (the same 21 differences as the final run). Final run: all sections in one process, 23 threads, nothing else running on the machine.
- Discipline. No threshold, definition or reading was changed in the light of a result. The parameter table was not edited (0 typo edits). Nothing in `src/`, `test/`, `data/`, the spec, the preregistration or `_config/` was touched. The module and the independent side disagree in 21 of 14,146 compared quantities; none was resolved by changing either side (findings 1 and 2).

## Evidence per criterion (Amendment D89)

(S) criteria rest on the folder suite `test/analysis/` (re-run in this stage with `julia --project -t 4 --check-bounds=yes`: 1213 of 1213 pass in 2 min 57.7 s); the file names are those of the test plan. (V) criteria quote the number and the line of `analysis-p2-verify-log.txt`.

| Criterion | Evidence | Result |
|-----------|----------|--------|
| D89a Dispatch, P1 unchanged | (S) `an_p2_dispatch.jl` "D89a P1 path unchanged" (golden values taken before the amendment). (V) log P1UNCH: `analyze` on the real P1 inputs, status ok in 26.9 s; the new tables file has 18103 lines, the accepted `p1-le-benchmark-tables.md` 18103; the only two lines that differ are "module version" (0.0.1 to 0.1.2) and "git commit" (d31b655 to b792724); the seven figure data digests are equal; the P1 table has 37 keys and is read by the P1 reader. CLI: a P1 table with the P2 preregistration, and `--p1-store` with the P1 preregistration, exit code 1, nothing written | PASS |
| D89b Parameter table | (S) `an_p2_params.jl`. (V) TRANSCRIPTION: the module reads the real table (58 keys, 11 machine keys, the quote and text checks of E22 and A2, `bootstrap_seed`, label keys against the P1 table); my line-by-line transcription below: 47 keys with a quote and 11 machine keys, 0 failures, 0 edits, 4 weak quotes (INFO) | PASS |
| D89c Input refusals | (S) `an_p2_refusals.jl` (54 malformed inputs, every rule a to p of E36). The real inputs raised no refusal in any of the runs below (E36 a to p all passed on 695 chunks) | PASS |
| D89d Reader and identity | (S) `an_p2_refusals.jl` "D89d reader and identity". (V) independent reader: 883 panel and 2304 grid orbits satisfy E24 (stored `params` and `ic` texts against the panel file and the `[[run.h4]]` tables, K index and label recovered); the module's sets and statistics read from the same items equal the independent ones (D89w) | PASS |
| D89e Budget rule | (S) `an_p2_verdicts.jl` "D89e budget rule" (the two real grids). (V) step used per indicator and budget: 60 of 60 equal to the independent value; 200, 447, 891, 1778, 3548, 7079, 14125, 28184, 56234, 112202 for LE, SALI, FLI, MEGNO and 200 ... 28184, 51200, 102400 for RLI and STR (E26, P15) | PASS |
| D89f Verdicts R1, D89g Verdicts R2 | (S) `an_p2_verdicts.jl` (boundary values, undecided causes, flags). (V) verdict counts chaotic, regular, undecided, not applicable for W2, S2, T_s and the regular set, 8 indicators, 10 budgets, both rules: 2560 of 2560 equal | PASS |
| D89h Detection | (S) `an_p2_detection.jl` (118,098 sequences). (V) outcome counts, M_X, number of orbits and the distribution of t_det for 5 sets, 8 indicators, both rules: det 400 of 400 and tdet 426 of 426 equal | PASS |
| D89i Agreement | (S) `an_p2_agreement.jl`. (V) F1: A, kappa, N and the four counts of 28 pairs, 2 sets, 2 budgets (784 of 784), Spearman (112 of 112), mean kappa (4 of 4; the H1 statistic differs from the independent mean by 5.6e-17) | PASS |
| D89j AUC and regular range | (S) `an_p2_agreement.jl`, `an_p2_pipeline.jl`. (V) AUC with class counts, pooled and by source, per K value (4680 of 4680 equal), regular range lo, hi, median, n (2880 of 2880), leaving budget counts (112 of 112) and the T_s summary per K (192 of 192) | PASS |
| D89k Sets | (S) `an_p2_sets.jl`, `an_p2_variant.jl`. (V) set sizes total, panel, grid and per K for 17 sets (150 of 150), drops (5 of 5), W2 and S2 overlap 0; conservation of counts: T1 holds 3187 orbits once and rho counts sum to 3187 | PASS |
| D89l Reference class | (S) `an_p2_sets.jl` "D89l". (V) rho counts by source and K (48 of 48), per cell (87 of 87), P1 label pairs per cell (47 of 47); see O10 below | PASS |
| D89m Statistics, D89n Decisions | (S) `an_p2_stats.jl`. (V) four versions (R1 I_7, R1 I_8 on E8, R2 I_7, R2 I_8): 162 + 181 + 162 + 181 values compared; all 28 decision words (six hypotheses and H4 on T_s^design, four versions) are equal on both sides (the H4 design word of the extension is NOT DECIDABLE on both); ROBUST flags 12 of 12 equal; after the stage 04 rework all 181 + 181 + 162 + 162 equal (the extension design statistic is on T_s^design intersected with E8; the R2 H3 `npass` differs by 1, a tolerated count, finding 2) | PASS (rework) |
| D89o Bootstrap | (S) `an_p2_bootstrap.jl`. (V) interval ends of H1 to H5 against the independent generator: 100 compared, all within max(0.05 of the interval width, one lattice step) after the rework; with the 0.05 rule alone 93 of 100 (the 7 that were not: 4 belong to the design statistic of the extension, 3 are one lattice step) | PASS under the tolerance settled at the checkpoint (finding 2) |
| D89p Output content | (S) `an_p2_output.jl`. (V) the tables file has the sections Header, Decisions, T1 to T17, F1 to F6 with 6 `digest:` lines, readings, deviations, exploratory items; no clock and no absolute path in it (0 matches of a grep); the six digests equal the digests of the data tables; T2: 320 pooled rows with chaotic + regular + undecided + not applicable equal to the set size, by-source and by-K groups sum to the set sizes (I17) | PASS |
| D89q Determinism | (S) `an_p2_output.jl`. (V) tables file at 23, 1 and 4 threads byte-identical (SHA-256 25d5d244f074...ad0e6e, 1.2 MB, 14113 lines), also the file written by the command into another folder | PASS |
| D89r Read-only | (S) `an_p2_output.jl`, `an_p2_refusals.jl`. (V) SHA-256 of 2066 input files (695 P2 chunks, 1363 P1 chunks, both manifests, panel, both preregistrations, both parameter tables, `numerics-defaults.toml`) equal before the run, after run 1 (23 threads), after run 2 (1 thread), after run 3 (4 threads) and after a refused run (`--out` inside `data/`, folder not created); listing of the data tree (path, size, mtime) has the same SHA-256 after run 1, after run 3 and after the refusal | PASS |
| D89s Command | (S) `an_p2_output.jl` "D89s command". (V) `bin/chaoslab.jl analyze` on the real P2 inputs with the P1 reference inputs by the study layout: exit code 2 and the line "figures skipped"; its tables file is byte-identical to the file of the module call with explicit `p1_*` keywords; only the tables file in `--out`; unknown option: exit code 1 | PASS |
| D89t Layering | (S) `an_p2_dispatch.jl` "D89t". `src/Analysis` has four `using` lines (Experiments, SHA, StableRNGs, TOML), no plotting import | PASS |
| D89u Figure data | (S) `an_p2_output.jl` "D89u". (V) the six digest lines of the tables file equal the digests of the data tables; the command exits 2 with "figures skipped". The PDF half is not testable: no drawer exists (P11) | PASS (PDF half not applicable) |
| D89v Store identities | prefix consistency L1 and tangent chain, horizon 10^6 against 10^7, 172 orbits of S_D: 340,560 values compared, 0 differences (independent reader); flag counts recomputed from the stored words equal the manifest's `indicator_flags` for all 213 window chunks (136 with listed flags), 0 differing (reader) and the module run did not refuse (E36 p); E24 for all 3187 orbits; omissions fq1 and bk1 = {883} equal to `not_applicable`; stored `control_deviation` equals its recomputation (10 controls, largest difference 0); drift flag: only orbit 881; K = 0 line printed in the log for 10 windows (fq1: 8.8e-06 at W = 100 down to 5.6e-17 at W >= 12800; pa1: 2.2e-16 or 1.1e-16 at every window; 5 orbits of the gate set have no finite fq1 distance, the five with p_0 = 0), equal to the independent values (20 of 20) | PASS |
| D89w Independent recomputation | see "Comparison" below: 14,146 quantities compared, 14,146 equal after the stage 04 rework (21 differences before it; 1 of the 14,146 is a tolerated count) | PASS (rework) |
| D89x Real-run report | see "Real-run report" below: exit code 2 (figures skipped), 48.8 s at 23 threads, 1.88 GiB peak | PASS |
| D89y Pilot rule check | (S) `an_p2_verdicts.jl` "D89y" (192 verdicts from the 24 tabulated rows). The independent decision layer reproduces the 192 tabulated verdicts from the tabulated scores (192 of 192). The (V) half needs a pilot store, which is not in the repository (P30) | PASS (V half not applicable) |

## Transcription of the parameter table (E22, line by line)

For each of the 58 keys of `workspaces/study/stages/04-analyze/references/p2-indicators-analysis-parameters.toml` the script checks: the section exists in the frozen preregistration (headings outside code fences, text up to the next heading of the same or a higher level); the quote occurs verbatim in that section; `text` occurs in the quote; the value equals the number I read in the cited sentence (the expected values are typed in the script from reading, not copied from the table), with the type (integer or real) the sentence implies. The 11 machine keys equal the values parsed from the two TOML blocks (arms `lyap_p`, `win_p`, `dir_tg`, indicator `fq1`, `run.health`, the 2304 `[[run.h4]]` tables, `panel_prereg_sha256_accepted`); `budgets[k]` = 2 `windows[k]`, `N_c` = twice the largest window = last budget, `b_primary` = `budgets[7]` (b_6), `grid_orbit_count` = 8 x 24 x 12 = 2304 = number of tables; `bootstrap_seed` = seed of `numerics-defaults.toml`; the four label keys equal the P1 table.

Result: 58 of 58 keys pass; 0 values wrong; 0 edits made to the table (no typo found). Four keys have a weak quote, reported as INFO: `label_lambda_min` and `label_growth_ratio` quote "the same two label rules (P1 "Label rules"", `label_rotation_tol` quotes "structure axis from the rotation numbers of the recorded orbit points", `label_estimator_bound` quotes "above the P1 H1a bound at N_1 gives U". None of the four quotes contains the number; the values (1e-4, 5, 1e-7, 1e-8) are held by the P1 table and by the P1 text ("Label rules": 1e-4, r(N) >= 5, 1e-7; "H1": <= 1e-8), which the script checked, and the R1 row of LE in the same preregistration states "the same two numbers" (1e-4, 5) on purpose. The quote check cannot protect these four against a wrong number; the equality with the P1 table can.

## Real-run report (D89x)

The numbers below are the log lines of the REPORT section, recorded as they came out; no interpretation (that is the study's stage 04).

Run: `Chaoslab.Analysis.analyze` for study `p2-indicators` through the module API, `figures = true`, P1 reference inputs given by keyword, store `data/p2-indicators` (695 chunks, 26,964 items) and `data/p1-le-benchmark` (1363 chunks) read-only.

| Quantity | Value |
|----------|-------|
| Status, exit code | `figures_skipped`, exit code 2 (the drawers of F1 to F6 do not exist, P11) |
| Wall time, 23 threads, first call (compilation included) | 48.8 s (limit 900 s in the script; spec limit 15 min). Later calls in the same process, without compilation: 38.2 s with `threads = 1`, 29.2 s with `threads = 4`; the P1 analysis alone (P1UNCH) 26.9 s |
| Peak resident memory after the 23-thread run | 1.88 GiB (limit 8 GB); after all three runs 2.55 GiB; allocated 34.7 GiB, garbage collector 2.7 s |
| Tables file | 1.2 MB, 14,113 lines, 6 digests |

Set sizes (total; panel; grid; by K value):
- E_all 2974; 789; 2185; K = 0: 80, 0.5: 425, 0.75: 258, K_c: 442, 1.25: 265, pi/2: 454, 2: 279, 3: 286, 2 pi: 405, 10: 80.
- E8 2236; 51; 2185; 0: 8, 0.5: 275, 0.75: 258, K_c: 276, 1.25: 265, pi/2: 283, 2: 279, 3: 286, 2 pi: 298, 10: 8.
- W2 214; 125; 89; 0.5: 99, 0.75: 23, K_c: 83, 1.25: 1, pi/2: 5, 2: 2, 2 pi: 1. W2_8 90 (1; 89).
- S2 448; 161; 287; 2 pi: 368, 10: 80. S2_8 303 (16; 287). W2 and S2 overlap 0.
- T_s 30; 15; 15; 0.5: 8, 0.75: 2, K_c: 12, 1.25: 1, pi/2: 5, 2: 2. T_s_8 15 (0; 15).
- T_s^design 96; 96; 0; 0.5: 21, K_c: 37, pi/2: 37, 2 pi: 1.
- Regular set 1159; 315; 844. Gate set 80 (K = 0 panel cells, classes T and R). S_D* 109, S_L* 51.
- Reference class C, R, U: panel 475, 315, 93; grid 1341, 844, 119. Dropped from E_all: 213 (13 controls, 200 with rho = U); no non-control orbit carries the drift or excluded flag (orbit 881 is a control).
- Decidability against `min_set` = 20: T_s decidable (30 orbits); T_s_8 not decidable (15); W2 214, S2 448, T_s^design 96, S_D* 109, S_L* 51 decidable.

Decisions of the primary version (R1, I_7, E_all), statistics and 95 percent intervals:

| Hypothesis | Statistics | Verdict |
|------------|-----------|---------|
| H1 on W2 (214 orbits) | U 0.30374 [0.24299, 0.36449]; A_min 0.75701 (MEGNO-E1); U on S2 0; mean kappa 0.48043; 21 A_XY in the log | SUPPORTED |
| H2 on S2 (448) | M_X = 200 for all seven indicators; Q 1 [1, 1] | NOT SUPPORTED |
| H3 (W2 214, S2 448) | M_W2: LE 1600, SALI 1600, FLI 3200, MEGNO 3200, E1 2400, E2 400, WBA 200; M_S2 = 200 for all; R_X [lo, hi]: LE 8 [8, 16], SALI 8 [4, 16], FLI 16 [8, 16], MEGNO 16 [8, 16], E1 12 [4, 32], E2 2 [1, 4], WBA 1 [1, 1]; indicators passing 3 | NOT SUPPORTED |
| H4 on T_s (30) | phi_X: SALI 0.2, FLI 0.26667, MEGNO 0.23333, E1 0.7, E2 0.63333, WBA 0.06667 (30 decided each); Delta_phi 0.63333 [0.46667, 0.8]. On T_s^design (96): Delta_phi 0.21875 [0.13542, 0.30208] | SUPPORTED; design: NOT SUPPORTED (statistic of record) |
| H5 on W2 (214) | F_X: LE 0.17757, SALI 0.18224, FLI 0.20561, MEGNO 0.23944, E1 0.20093, E2 0.08411, WBA 0.02804; F_med 0.18224 [0.13084, 0.22897] | SUPPORTED |
| H6 (S_D* 109, S_L* 51) | delta_dir: LE 0 (7 dropped), SALI 0.04587, FLI 0, MEGNO 0 (9 dropped); delta_prec: LE 0 (3), SALI 0, FLI 0, MEGNO 0 (4), E1 0.01961, E2 0, WBA 0 | NOT SUPPORTED |

ROBUST flags (primary against the I_8 extension on E8): H1 ROBUST, H2 ROBUST, H3 ROBUST, H4 NOT COMPARABLE (T_s_8 has 15 orbits), H5 ROBUST, H6 NOT ROBUST. Extension verdicts: H1 SUPPORTED (U 0.27778 [0.18889, 0.37778], A_min 0.76667, 90 orbits), H2 NOT SUPPORTED, H3 NOT SUPPORTED (1 indicator passes), H4 NOT DECIDABLE, H5 SUPPORTED (F_med 0.18333 [0.11111, 0.25]), H6 SUPPORTED (delta_dir SALI 0.07843).

R2 version (reported, not decisive): H1 SUPPORTED (U 0.22897 [0.17290, 0.28505], A_min 0.79907), H2 NOT SUPPORTED, H3 NOT SUPPORTED (1 indicator passes), H4 SUPPORTED (Delta_phi 0.7), H5 NOT SUPPORTED, H6 SUPPORTED. The full statistics of the extension and of R2 are in the log.

Rule-validity gate (share of the 80 gate orbits decided and chaotic at the last budget; limit 0.10): LE, SALI, FLI, MEGNO, E2: 0 of 80 decided; E1 and WBA: 0 of 75 decided (5 undecided each, the orbits with p_0 = 0 of spec P17; their ids were not re-derived here); RLI: 0 of 8 decided (72 not applicable). No indicator violates the gate; no second version was produced.

## Comparison with the independent implementation (D89w, O6, O10)

Method. The comparison script builds, from the result of the module call, a table under the names of the independent results file (documented at the top of `verify/analysis-p2-decision.py`) and compares every independent line that has a counterpart: integers, texts, verdicts, flags and sets exactly; reals to 1e-12 relative (the K = 0 line to 1e-12 relative or 1e-15 absolute: values of order 1e-16 differ by one rounding); bootstrap interval ends (U, Q, R_X, Delta_phi, F_med, lower and upper) within 0.05 of the width of the module's interval. The R2 detection outcomes are computed from the module's R2 verdict series by the module's detection function (the result does not expose them); the extension and R2 rows come from the stored decisions.

Result after the stage 04 rework (re-run 2026-10-30, T_s^design intersected with E8 in the extension, tolerance max(0.05 of the width, one lattice step) for U, Delta_phi and F_med ends): 14,146 quantities compared, 14,146 equal; no independent name lacks a module counterpart. Of the 14,146, 100 are interval ends, all within the tolerance (93 within 0.05 of the width; 3 differ by one lattice step; 4 that belonged to the design statistic of the extension now equal). One count, R2 I_7 H3 `npass` (independent 0, module 1), differs by the one lower end that moved and is a tolerated count (logged as `DIFF | tolerated count`). Before the rework: 14,125 equal, 21 differ (R1 I_7 162 of 162, R1 I_8 170 of 181, R2 I_7 160 of 162, R2 I_8 173 of 181), recorded below.

The 21 differences before the rework, with both values (now settled; the log keeps only the tolerated count):
1. T_s^design in the extension, 17 values (finding 1): R1 I_8 H4 design: n 2 (independent) against 96 (module); phi SALI, FLI, MEGNO, E1, E2, WBA 0, 0, 0, 0, 0, 0 against 0.05208, 0.05208, 0.0625, 0.22917, 0.07292, 0.01042; Delta_phi 0 against 0.22917; interval ends 0, 0 against undefined, undefined. R2 I_8 H4 design: n 2 against 96; phi E1, E2, WBA 0, 0, 0 against 0.22917, 0.07292, 0.01042; Delta_phi 0 against 0.22917; interval ends 0, 0 against undefined, undefined.
2. One lattice step of a bootstrap percentile, 3 interval ends and 1 count (finding 2): R1 I_8 H1 U upper end 0.36667 against 0.37778 (set of 90 orbits, step 1/90); R2 I_8 H1 U lower end 0.13333 against 0.14444 (step 1/90); R2 I_7 H4 Delta_phi lower end 0.53333 against 0.56667 (30 orbits, step 1/30); R2 I_7 H3 `npass` 0 against 1 (R_E1 = 12, lower end 3 against 4, the rule is "> 3").

O10. The panel orbits per (K index, class, instability label, structure label): the labels computed by the module inside the P2 run (the base labelling on the P1 inputs) give 47 cells and 883 orbits, 0 cells differ from table T6 of the accepted `p1-le-benchmark-tables.md`; the independent reader got the same 47 cells and its final labels equal those of the D60 independent recomputation for 883 of 883 orbits (reader lines in the log).

O9 identities from the independent reader (copied from the log): prefix consistency for 172 S_D orbits, 340,560 values, 0 differences; flag counts of 213 window chunks (136 with listed flags), 0 differing; omissions fq1 and bk1 {883} equal to the manifest; control_deviation of 10 controls, largest difference 0.0; drift flag recomputation: orbit 881 only; E24: 2304 grid and 883 panel orbits, 0 bad; 695 files, 16,165 series items, 10,111 window items, 0 anomalies. O8: the pilot rule check, 192 of 192 verdicts. O7: kappa and AUC closed forms, the exact enumeration of the bootstrap of U and Q over 6^6 resamples (generator within 4 standard errors) and Spearman against `scipy.stats.spearmanr`, 16 self checks, 16 PASS.

Not compared (the independent side did not compute them, or the module does not expose them): the second versions of T3 to T8 (no indicator violates the gate, none exists); the number of dropped bootstrap replicates (`U_repdrop`, 0 on the independent side; the module reports no such count); per-orbit verdicts and scores (only aggregates are compared, through the counts, medians, AUC and ranges); the by-source and by-K point values of T2 and T3 to T8 other than the set sizes, rho counts, AUC and range; the F5 lattice rows, the F6 percentiles of the stretching distance, T17 cost and the layout and digests of the written tables file; the exact bootstrap draws (the generators differ by design); R2 detection outcomes beyond what the stored decisions give.

## Findings of the verification

1. **Reading difference, T_s^design in the I_8 extension (LOW; no decision depends on it).** The independent side intersects the design ensemble with E8 (2 orbits) because E29 writes "a set restricted to E8 ... with the suffix 8"; the module uses the whole ensemble (96 orbits) with the label set name `T_s_design_8`. E33 names the sets restricted to E8 (W2_8, S2_8, T_s_8, S_D* and S_L* intersected with E8) and T_s^design is not among them; the module follows the list and records the reading (`P2_READINGS`, module reading 4 of the implementation notes, review reading 11). The spec is silent, so neither side is wrong; SETTLED at the checkpoint (2026-10-30): T_s^design is intersected with E8 in the extension (module, tests, reading 4 and 11 changed). 17 of the 21 differences, all gone after the rework. The statistic is "of record", not a decision, and the extension of H4 is NOT COMPARABLE on the primary sets anyway. 
2. **Tolerance of D89w and O6 for small sets (LOW).** For W2_8 (90 orbits) and T_s (30 orbits) the bootstrap statistics live on a lattice of step 1/90 and 1/30, wider than 0.05 of the interval width (0.008 to 0.015). Three interval ends differ by exactly one lattice step, which two correct generators can produce when the 2.5 or 97.5 percent point sits at a lattice point (the standard error of the empirical distribution function at B = 10^4 is 0.0016); one count follows (R2 H3, `npass` 0 or 1 of 7, needs 5, no effect on the verdict). SETTLED at the checkpoint (2026-10-30): an interval end of U, Delta_phi and F_med agrees within max(0.05 of the width, one lattice step); the ratios Q and R_X keep 0.05 of the width; a count that follows from a moved end may differ by the number of moved ends (one count). The tolerance was settled after the data were seen, which the owner did knowingly; I did not verify the mechanism by re-running either generator with other seeds.
3. **Reading difference, flip counts (INFO).** The independent decision layer counts an orbit without an item of an indicator (RLI outside E8) at flip count 0; the module counts only orbits with an item. Compared after adding the number of such orbits at flip count 0; 116 of 116 equal.
4. **Result not exposed (INFO).** `analyze` returns the R1 detection table only; the R2 detection outcomes (needed for the R2 rows of T10 and for a check of F2 under R2) are recomputed by calling the module's detection function on the R2 verdict series. A result field would make the R2 version checkable without that call.
5. **Weak quotes (INFO).** Four label keys of the parameter table carry quotes without a number (see the transcription); their values are protected by the equality with the P1 table, not by the quote check. A2 of the owner decisions already states the four label keys get a quote from the section that cites them.
6. **Module and stores (positive).** The 695 P2 chunks and the P1 labelling pass all refusal rules on the first call; 213 chunks with window items satisfy the flag and omission identities; the stored control deviations equal their recomputation; the module and an independent HDF5 reader agree on 2304 + 883 orbit identities.

## Independent review of the code (fresh Sonnet subagent, spec and code)

Read-only review of `src/Analysis/p2_*.jl` against the spec: 0 CRITICAL, 1 HIGH, 7 MEDIUM, 4 LOW. Fixed in b792724 (tests in `test/analysis/an_p2_review.jl`, written from the spec before the fixes: 79 tests, 46 failing and 2 erroring before, all passing after) except the three findings that became readings.

| Finding | Disposition |
|---------|-------------|
| 1 (HIGH) K pooling by K index instead of K value (AUC, regular range, leaving budget, share below the median, by-K rows of T2 to T14) | fixed: pooling by the binary64 K value, panel and grid together; keys of `res.auc`, `res.range`, `res.below_median` are K values |
| 2 type 7 percentile with an infinite neighbour gave NaN | fixed |
| 3 a chunk with a stripped attribute (`plan_sha256`, `prereg_sha256`, `panel_sha256`, `status`, `study`) was accepted | fixed: named refusal in P2 |
| 4 a numeric `text` could differ from the value | fixed: must equal the value |
| 5 leaving budget when the regular range is missing at a budget | fixed: such a budget counts as "not above" |
| 6 K = 0 line with a non-finite distance gave NaN | fixed: filtered and counted |
| 7 unknown indicator kind and wrong bit names were a `KeyError` or accepted | fixed: named refusals |
| 9 a control with a positive exponent and no `control_deviation` was accepted | fixed: refused |
| 12 a growth ratio overflowing to Inf gave an Inf score | fixed: chaotic, no score |
| 8, 10, 11 the 0.05 share of dropped replicates and the design class S are module constants; the 5 percent rule is joint per label; T_s^design is intersected with E8 in the extension (changed at the checkpoint) | kept as readings, listed in `P2_READINGS` (finding 1 above is the third) |

Review closed: no open CRITICAL or HIGH.

## Ambiguities listed by the independent author

The independent author's own list of 12 items is not in the repository (the results file has no header). The reading differences that this comparison made visible, which are the ones that matter for the module: (a) T_s^design in the extension (finding 1); (b) orbits without an item in the flip-count distribution (finding 3); (c) the dropped-replicate rule (more than 5 percent of the replicates undefined makes the interval undefined), applied per statistic by the independent side and jointly per label by the module (P20); it matters for the H4 design statistic of the extension, where RLI is decided on 2 of 96 orbits and many replicates have no RLI orbit (the module interval is undefined, the independent one, on its set of 2, is defined); (d) the leaving budget when the regular range is missing at some budget: the independent side skips the K value, the module counts the budget as "not above" (review fix 5); no K value that has T_s orbits has an incomplete range on the real store, so no difference arises. Items (a) to (c) explain every difference found; (d) has no effect on a number.

## Audit

| Check | Result |
|-------|--------|
| Evidence: every acceptance criterion D89a to D89y has a number or a test file | pass |
| Review closed: no open CRITICAL or HIGH finding | pass |
| Real-run limits (D89x): time at most 15 min, memory at most 8 GB | pass (48.8 s, 1.88 GiB) |
| Independent recomputation (D89w): all exact quantities equal, interval ends within the tolerance | pass after the rework: 14,146 of 14,146 (one count tolerated); the tolerance and the T_s^design reading were settled by the owner after the first comparison, so this pass is not a blind test of them |

## Not verified

- E36 o (an input changed while the run goes): the module compares the hashes taken at the start with those before the tables are written; there is no test hook, and the real run was not disturbed on purpose. My own SHA-256 of all inputs before and after each run is the only check.
- The figure drawers and the six PDF files of F1 to F6 (P11): do not exist; the PDF half of D89u is not applicable. The data and digests are verified.
- Timing on cold caches: the P2 and P1 stores had been read by the SHA-256 pass of the script just before the first run, so the 51 s is a warm-cache time. The 1-thread and 4-thread runs are keyword settings (`threads = 1`, `threads = 4`) inside a 23-thread session, not separate `julia -t 1` processes; the reading of the files is serial in every case.
- The bootstrap generators are not compared draw by draw (O6 compares the interval ends).
- The independent side did not compute the items listed under "Not compared" above.
- The pilot store (P30): the (V) half of D89y.
- Coverage of the new files under the suite was not re-measured; the GPU does not apply (no device path in the module).

## Amendment D91: figure drawers F1 to F6 of study p2-indicators (stage 04)

Verify scripts: `verify/analysis-p2-figures.jl` and `verify/analysis-p2-figures.py` (independent author, written from the spec only; no file of `src/`, `ext/` or `test/` was read). Raw output: `analysis-d91-verify-log.txt` (final run 2026-10-30 17:36, 23 threads, after the review fixes below: PASS 143, FAIL 6, INFO 18). Contact sheets for the owner visual check: `analysis-d91-F1.png` to `analysis-d91-F6.png` (72 dpi). Real P2 command under the verify project: exit 0, six PDF files, tables equal the accepted stage 04 file except the git commit line, 94.7 s whole command, 3.21 GB peak, figure step 24.0 s cold (0.9 s warm).

### Evidence per criterion

| Criterion | Evidence | Result |
|-----------|----------|--------|
| D91a command, six files | exit 0, six PDFs with `%PDF-` signature and one page, no skipped or failed line | pass |
| D91b figure is the block, digests | digests recomputed by own code from the accepted file, the rows the drawer received and the digest lines agree for F1 to F6 | pass |
| D91c tables unchanged, digests | tables file equals the accepted file except the git commit line (with and without figures); six digest lines equal the D91c values | pass |
| D91d failing drawer | non-positive and zero D_str: exception names figure and row ("F6 row 100: ..."); after the fix no figure file is written or replaced | pass |
| D91e P1 unchanged | P1 digests F1 to F7 equal the accepted ones; P1 tables differ only in the version and git commit lines; a copy of the tree before D91 (fc4cace) gives all 11 P1 files pixel-identical to the current module; the paper-folder copies of F3 to F7 and F1a to F2b are pixel-identical. The six tracked study-folder files `p1-le-benchmark-fig-F1` to `F6.pdf` are not (see finding V1) | 6 FAIL against the stale tracked files; pass against the pre-D91 tree |
| D91f layout | widths 991 to 1000 pt, heights at most 1.3 times the width, smallest text scaled to 504 pt 8.06 to 8.14 pt | pass |
| D91g text | all names, titles, N ranges, notes and legends present once as specified | pass |
| D91h plotted numbers | all F5 cells, 2445 F6 points, 320 F3 points, 700 F4 points, 200 F4 gaps empty, F1 cell colours decoded within 0.0049 of the colour bar span (limit 0.005; thin margin) | pass |
| D91i drawing counts | panels, lines and markers per figure as specified; the F4 chance-level line is absent in the two empty panels (E42: such a panel shows no lines) | pass, 1 INFO |
| D91j time and memory | 24.0 s cold figure step (bound 120 s), 3.21 GB (bound 8 GB) | pass |
| D91k determinism | two runs: tables byte-identical, six digests equal, rendered pixels of the six PDFs equal (raw bytes differ in embedded metadata) | pass |
| D91l owner visual check | contact sheets written | open: the owner decides |
| D91m edge fixtures | seven fixtures (no RLI in S2, no RLI anywhere, three K with all AUC nothing, 30-orbit set, one K index, class count changing with the step, one budget in F1) draw six files each | pass |
| D91n layering (S) | folder suite 1287 of 1287 after the fixes (3 min); dependency files unchanged | pass |
| D91o read-only | SHA-256 of 2436 data files, 28 other inputs and the tree listing equal before and after | pass |
| I22 to I29 | figure is the block, one encoding (indicator colours and markers equal in F2 to F4), no silent alteration, dependency files unchanged | pass |

### Independent code review (fresh Sonnet subagent, spec and code)

0 CRITICAL, 0 HIGH, 2 MEDIUM, 10 LOW.

| Finding | Disposition |
|---------|-------------|
| M1 a failure in figure k left F1 to F(k-1) of the new run beside F(k+1) to F6 of an earlier run | fixed: the six pages are drawn to temporary files and moved into place only when all six are drawn; a failure leaves the folder as it was; re-verified (failed-figure run writes only the tables file) |
| M2 F4: RLI has no AUC at K = 0.5 and K = 6.28 (and every AUC is nothing at K = 0 and 10); the missing line was not stated on the figure | fixed: a fourth title line in every panel ("AUC of every indicator", "no AUC: RLI", "no AUC at all"); F4 page 1000 x 870 pt |
| L1 every error names its figure | fixed: any error of a drawer is rethrown as "F<n>: ..." |
| L2 duplicate rows were summed (F2) or the first taken (F3, F4) | fixed: refused with a named error |
| L3 axis labels missing (F1 matrices, F2 bottom row y axis) | fixed ("indicator"); F1 matrix size reduced to keep the page under 1000 pt |
| L4 a `nothing` cell in F1, F3, F6 or "n.a." in F5 raises | kept: the real block has none; the error names figure and row (E38 b) |
| L5 E42 text is set with hard line breaks | kept: a text check must normalise white space (done in the verify script) |
| L6 F2 "N differs from the set size" test uses the maximum N | kept: the label compares with the largest N of the set, equal to the set size in the real block |
| L7 grey-scale: T_s (grey) against R (blue) are close, `:vik` has two dark ends | open for the owner visual check (D91l) |
| L8 `p2order` appends unknown sets in dictionary order | kept: the real block holds only the four known sets |

No test changed. The fixes live in `ext/p2_figures.jl` and `ext/p2_figures_b.jl` only; the package suite does not load them except for the layering checks.

### Findings of the verification

V1. **Tracked P1 study-folder figures are stale (LOW, no D91 cause).** The tracked `p1-le-benchmark-fig-F1` to `F6.pdf` in `workspaces/study/stages/04-analyze/output` date from 91abcfa (2026-10-20), before the figure redraw (caaf95a, 812d484); F3 has another page size (480 x 315 against 615 x 360 pt). D91 did not change P1 (pre-D91 tree gives identical pictures). The six FAIL lines of D91e are this and nothing else. Owner decision: refresh the tracked files from a P1 run (a P1 change, not D91), or accept; the criterion D91e is met in substance by the comparison with the tree before D91.
V2. **F4 holds nine lines in six panels only (INFO).** The spec's "nine lines each" holds in six of ten panels; two panels have eight (RLI has no AUC), two none; each panel now says so in its title.
V3. **Kappa decode margin (INFO).** The largest F1 colour error is 0.0049 of the span against the limit 0.005; the render at 4 px/pt, not the figure, limits it.

## Not verified (D91)

- D91l: the look of the figures, colour-vision and grey-scale legibility (owner).
- The full package suite was not rerun after the review fixes (17 min); only the analysis folder suite (1287 of 1287). The fixes are in the extension, which the package suite does not load; the base run at b2d7ebe passed.
- Aqua, JET, coverage and bench for the extension files.
- The P11 reading text in `P2_READINGS` still says that the figure drawers do not exist; it sits in the accepted stage 04 tables file and is left unchanged (owner decision).
