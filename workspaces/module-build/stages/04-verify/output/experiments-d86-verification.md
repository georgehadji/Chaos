# Verification: experiments, Amendment D86 (recording of the chaos indicators)

Module-build stage 04, steps 1 and 3, for the third scope of `experiments` (decisions D86 and D87). Code at HEAD 5b1f91d (stage 03 checkpoint); the verify scripts are new and uncommitted. Nothing under `src/`, `test/` or `bench/` was edited. The independent review (step 2) is a separate run, merged at the end.

Raw output: `experiments-d86-verify-log.txt` (CPU script) and `gpu-rtx5060-d86-verify-log.txt` (device script); figure `experiments-d86-device.png`. Scripts in `verify/`: `experiments-indicators.jl` (main) with `-prep.jl`, `-campaign.jl` (child process), `-plan.jl` (child process, D86b), `-reader.jl` (independent reader, HDF5.jl and SHA only, no Chaoslab), `-series.jl`, `-closed.jl`, `-hash.jl`, `-memory.jl`, `-source.jl`; `experiments-indicators-device.jl`; the reduced block `experiments-indicators-block.toml` (a verification input, values chosen for the check, not study values; nothing like a p2 value).

## Runs

| Run | What | Result | Time |
|-----|------|--------|------|
| `experiments-indicators.jl`, 23 threads | F1, D86b [full], 4 campaigns, reader, manifest, files, O18, O19, rerun, pilot, memory, source search | PASS 113, FAIL 0, INFO 63 (the first full run raised one script error in the pilot section, a sort of mixed keys; section PILOT was rerun alone and spliced into the log, noted in it) | 813 s plus 77 s for the pilot rerun |
| `experiments-indicators-device.jl`, 23 threads, RTX 5060 Laptop | table, per-spec records, parity, O21, Double64, device campaigns | PASS 32, FAIL 0, FALLBACK 4, INFO 14; a second run of sections RECORD, PARITY and PLOT (PASS 12, INFO 6, appended to the log) made the figure | about 480 s plus 90 s |
| ICM validator | `python scripts/icm_validate.py` | ALL CHECKS PASS | - |

Reduced campaign (D86y): the real preregistration with its first `[run]` block replaced by the reduced block, and a reduced copy of the real panel file (60 drawn orbits per K for the 6 K values, all 13 controls; one orbit per cell in the subset S_L, 23 cells, because the grammar has no per-cell selection). 2793 items, 352 chunks, 9 arms. One spec of each kind: tangent chain FLI, MEGNO, SALI at N = 4000; pair RLI then MEGNO (N = 300, finding F2); stretching 100 bins at N = 4000 with the counts of S_L; frequency, phase advance and birkhoff with windows 250 and 1000 (horizons 2000, 2001, 2000); Float64 on all panel orbits and Double64 on S_L; one injected chunk failure (NaN in a series column, attempts 1 and 2); one injected kill after 70 chunks and one during the data of chunk 71.

Quoted from the orchestrator (not rerun): experiments folder 6830 pass, 0 fail; full `Pkg.test()` 1,400,142 pass, 0 fail, 1 broken (outside experiments); coverage of `src/Experiments` 96.1 percent of 3440 lines (new `run_indicator*.jl` files 94 to 100 percent).

## Findings in the code or the spec

| # | Severity | File:line | Evidence | Note |
|---|----------|-----------|----------|------|
| F1 | HIGH (defect) | `src/Experiments/run_indicator_store.jl:86` (`_h5val`; raised through `_put_provenance!`, line 31, and `_put_attrs!`, `run_store.jl:75`) | Section F1 of the CPU log: the unpatched campaign (own process) writes 237 chunk files, then aborts with `MethodError: no method matching unsafe_convert(::Type{Ptr{Nothing}}, ::Double64)` at the first Double64 `stretching` chunk. The stretching reducer publishes `spectrum_lo`, `spectrum_hi`, `spectrum_cell_width` as Double64 and `_h5val` passes any `AbstractFloat` through. | Every Double64 stretching chunk is unwritable and the campaign dies with an exception (not a failed chunk). The suite did not catch it (no Double64 stretching chunk is written in the D86 tests). A one-line fix (convert a non-Float64 float to its decimal text, as the method for other types does) is for stage 03. All campaign scripts here define the missing method in the script only (`verify/experiments-indicators-prep.jl`, `EXI_PATCH=1`); every campaign below ran with that workaround. |
| F2 | note (design limit, not a defect) | R18 item 5, rli R8 | The pair on the linear control Y5 at N = 4000: the shadow orbit overflows at step 770, the chunk of the controls cell fails twice with `state = -Inf ... step 770, orbit 883` (section O18 of the log). The first campaign with the pair at N = 4000 failed that chunk and, with 352 chunks, was still inside the stop fraction. | The P2 preregistration must give the pair a horizon below about 700 for Y5, or no Y5 in the pair arm. The reduced block uses N = 300 for the pair. |
| F3 | LOW (spec wording) | E31 against `src/Experiments/run_indicator_health.jl:_ind_nan_ok` | For S = 0 the reducer sets `zero` only and `reldig` is NaN (28 rows, the fixed-point and Y5 controls). E31 says `reldig` is finite unless `reldig_undefined`; the code accepts NaN in `reldig` under `zero` or `reldig_undefined`. The independent reader implements the code's rule and lists each case as a NOTE. | The flag set by the code is `reldig_undefined` only when the value is NaN and `zero` is clear. The spec text should say "unless zero or reldig_undefined". |

Owner notes from the numbers:

- `frequency` with `mean_removal = true`: on the K_std = 0 line the recovered frequency differs from p/(2 pi) by up to about 3e-9 (22 of 48 values outside 16 eps in a first run of section O18). With `mean_removal = false` all 48 are within 16 eps (worst ratio 0.012). The closed form of naff-fma criterion 7(b) is for the second case; the preregistration should state the choice and know the cost.
- Cost: the frequency kind costs 4.1 to 6.6 microseconds of thread time per orbit-step at windows 250 and 1000 (the transform dominates; the orbit pass alone is 0.005 s per 12 pilot orbits), against 0.4 to 0.8 for the tangent chain and about 0.5 for the window reducers; the budget lines of the reduced block (mine, arbitrary) were exceeded by the frequency, pair and phase-advance specs (lines `COST` of the log). The pair arm at N = 300 shows the per-orbit setup cost (12 900 ns per orbit-step). The pilot of the real plan gives the numbers for the budgets.
- Device: records per spec on 1024 orbits (device wall ns per orbit-step against CPU 23 threads, last run): tangent 69 against 37 (device 1.9 times slower), pair 297 against 696 (device 0.43), phase advance 30.6 against 50.7 (0.60, but 1.12 in the first run: within noise), birkhoff 62 against 44 (1.4 times slower; 2.2 in the first run). `auto` follows the records; only the pair (and, in one of two takes, phase advance) resolves to the device. The record file `data/gpu-benchmark.toml` is unchanged (SHA-256 equal before and after); the per-spec records were written to the scratch folder.
- The criterion D86y reads "the study's data folder holds 0 files". The P1 campaign has since finished and `data/p1-le-benchmark` holds 1741 files (1363 chunks, the manifest rows and the `stale/` folder). The check made is that no entry under `data/` (path, size, modification time) changed: 2679 entries before and after; every store of the scripts is under the scratch folder.
- Reader independence: h5py is not installed and was not installed. The reader of O20 is Julia with HDF5.jl in a separate process without any Chaoslab code (as O13); the one tool outside HDF5.jl is `h5dump` of HDF5_jll (49 indicator datasets, 0 differences).

## Fixes after the review (2026-10-29)

Orchestrator, after the verifier and the reviewer reported. `src/` changed in three places; the verifier's first run (logs without the suffix) was on 5b1f91d with the F1 workaround `EXI_PATCH=1`; the rerun (logs `*-fixed.txt`) is on the fixed code with `EXI_PATCH=0`.

| Finding | Fix | Test | Evidence |
|---------|-----|------|----------|
| F1 (HIGH) | `_h5val` in `run_indicator_store.jl` passes only strings, `Bool`, `Base.BitInteger`, `Float32`, `Float64`; any other value (Double64, BigFloat, tuple) is stored as its exact text | `D86 review F1` (Double64 stretching campaign on h4, provenance `spectrum_lo` is text); fails with the MethodError when the fix is reverted | CPU rerun: the unpatched campaign now completes (352 chunks ran, 1 failed = the injected one); the only FAIL line of the rerun is the F1 line that expected the abort |
| R1 (MEDIUM) | frequency flag bits `refined1`, `refined2` are `!window.refined` (set only for a window that was not refined) in `run_indicator_record.jl` | the reference flags of D86i(d) now include both bits (`D86_FLAGS_CHECKED`) | folder suite; no unrefined window occurs in the fixture, so the bit is checked only as 0 |
| R3 (MEDIUM) | `plan_run` raises `ArgumentError` naming the spec when no orbit at all is selected for it | `D86 review R3` (panel without S_L orbits, arm with `orbits = "S_L"`) | passes alone and in the folder suite |

Owner decisions of the checkpoint (Amendment D88, spec section and decision-log row D88): R2 fixed (the group `/pooled/<id>` carries the attribute `range` and has no `counts` when the orbits of the chunk do not share one range; test `D88 pooled counts carry their range and are not written across ranges`, both cases occur in the fixture; the tests D86i and D86s that assumed `counts` in every chunk were adapted; the independent reader checks the rule), R4 settled on E25 (no `undersampled` or SALI-floor flag, no SALI-floor root attribute; R18 item 2 and E24 edited), F3 settled (E31 edited to "unless `zero` or `reldig_undefined`").

Runs on the final code: experiments folder 6844 pass, 0 fail, 0 error (6830 plus 14 new; 6 min 43 s); full `Pkg.test()` 1,400,156 pass, 0 fail, 0 error, 1 broken (outside experiments), 15 min; CPU verify rerun (`EXI_PATCH=0`) PASS 113, FAIL 1 (the F1 expectation above), INFO 62, 343 s, peak RSS 1.21 GiB, the reader accepts the D88 layout; device rerun PASS 32, FAIL 0, FALLBACK 4, INFO 14, peak RSS 2.62 GiB (`experiments-d86-verify-log-fixed.txt`, `gpu-rtx5060-d86-verify-log-fixed.txt`; the device figure was redrawn). Coverage was not rerun after the edits.

## Evidence per criterion

(S) rows quote the package suite; (V) and (G) rows are the scripts of this stage. Run: C = CPU log, G = device log.

| Criterion | Tag | Evidence | Carried by | Result |
|-----------|-----|----------|------------|--------|
| D86a | S | suite (69 malformed blocks, canonical text, merge); V adds: a spec with `threshold`, `class`, `detect` or `first_crossing`, a spec id in `control_methods`, a spec without `bins` are each refused naming the key | `ex_d86_grammar.jl`; C section D86W | PASS |
| D86b | S, [full] | real P1 preregistration and panel, code_version of the running campaign read from `chunk-00001.h5` (read only): plan hash 81b04cdf5fd6cfd2 equal for the code of ba21bcc (git worktree, removed afterwards), the current code and the manifest of the finished campaign; plan text hash and block hash equal; 15569 of 15569 item hashes equal; ordinals, keys and hashes of 1363 of 1363 chunks equal, old against new and new against the manifest; no `kind=` or `spec=` in an old chunk key. The `*.jl`, `Project.toml` and lock of ba21bcc equal those of 1c3de19 (one documentation file differs) | `verify/experiments-indicators.jl` (D86B), `-plan.jl`; C | PASS, 0 differences |
| D86c | S | suite; V: omitted orbits `fq: [883]`, `bk: [883]` in the reduced plan, plan hash equal on a second planning; 2793 items, 352 chunks | C (SETUP) | PASS |
| D86d | S | suite; V (O19): the text and hash of all 2793 items of the reduced plan written here from E27 and hashed with `sha256sum`: 0 text and 0 hash differences; the 6 spec texts (written from the TOML) and spec hashes equal; 352 of 352 chunk keys rebuilt from item fields equal; 352 chunk hashes (`sha256sum` of the chunk text) and the plan hash equal; the files of campaign B carry the same chunk_key, chunk_hash, plan_sha256, item hashes and spec_sha256 | `-hash.jl`; C | PASS |
| D86e | S | suite; V: the chain record of the 187 files of tg, pr, st equals the forms the module must set (all but the last non-replacing), spec text and hash equal the block's, one provenance group per member, three distinct `indicator` values in a chain of three | C (FILES) | PASS |
| D86f | S | suite only (injected frames and starts) | `ex_d86_chains.jl` | (S) |
| D86g | S | suite; V on 389 ok tangent orbits (Float64, Double64, wrapped, lifted): `t` equals the grid, all columns finite, FLI nondecreasing, SALI in [0, sqrt 2 (1 + 2 eps)] (0 violations); against the stand-alone replacing FLI, MEGNO, SALI run to the step on the identical orbit: 1740 values (31 orbits incl. 4 controls and 3 Double64, 5 seeded grid steps and the last), 0 differences, bitwise (O16, I20, I22); orbit point and monitor columns equal those of the L1 item for 716 orbits (I21); closed forms through the recorder (K_std = 0 shear, direction 0, wrapped and lifted): FLI_1 and MEGNO column 1 exactly 0 (10848 values), FLI_B within 8 eps (n + 3) of (1/2) ln(n^2 + 1) (worst ratio 0.021), SALI within the A15 form (0.0048), MEGNO column 2 within 6 B_S (0.0060) over 2712 rows each | C (FILES, O18) | PASS |
| D86h | S | suite; V: pair columns equal the stand-alone RLI and MEGNO on the PairMap with the stored xi_0 and the block frame, 768 values (Float64 and Double64), 0 differences; stored Dx_r equals the realized separation of `pair_start`; xi_0 unit vectors; Dx_r / Dx in [0.9994, 1] for all 396 pair orbits; RLI and DL exactly 0 for the shear and Y5 (1469 rows) | C | PASS |
| D86i | S | suite; V: D at 3 grid steps and the last, 17 scalars and the range equal the stand-alone final with workspace (155 values, 0 differences); counts stored exactly for the S_L orbits (419 orbits, 0 disagreements) and equal the stand-alone counts; pooled counts of 3 chunks equal the exact integer sum of the stand-alone counts of every orbit | C | PASS |
| D86j | S | suite; V: (a) capacity floor(f free / bytes) with B = 32 W_max (Float64) and 64 W_max (Double64) and halving (58823 to 29411); (b) 96 bytes per orbit and window; (d) frequency tables equal `frequency_map` on the plain-iteration series bitwise (30 rows Float64, 8 rows Double64); (f) the largest dataset of a frequency orbit has 254 elements, none holds the 2000-sample series. The orbit points equal the L1 item's (I21). (c) resident bytes after a chunk: suite only | `-memory.jl`, C (FILES) | PASS |
| D86k | S | suite; V: K_std = 0 line within 16 eps of p/(2 pi) mod 1 (no mean removal): 48 values, worst ratio 0.012; fixed-point controls give an ok chunk with the `degenerate` bit counted; flag rows over the campaign listed (below_floor 120, degenerate 26, mirror 258 and 268, near_zero, tie) | C | PASS |
| D86l | S | suite; V: phase-advance tables equal the pure functions on the plain-iteration series (30 rows, bitwise); lifted observable within 16 eps + max abs q eps of p/(2 pi) in both modes (96 values, worst ratio 0.035; the D87 bound); `undefined window` bit set on 12 rows, chunks ok | C | PASS |
| D86m | S | suite; V: weighted Birkhoff tables (4 observables, 2 windows) equal `weighted_birkhoff` on the plain-iteration series, 120 rows bitwise; `zero` set on 9 rows with `below_floor`; the K_std = 0 criterion 3 and 31 against 32 observables: suite only | C | PASS |
| D86n | S | suite; V: bytes(orbit) of every (spec, precision, mode) of the reduced plan equals E28 written from the reducers (14 of 14); sizes 72, 96, 200, 88 and stretching 104 + 8016; all worked numbers of R17 item 4 reproduced (18864, 424088 and 106022 at 4 GB, 136000, 58823, 16896, 473484, 8224, 972762, 8544, 936329); Double64 B = 256000; one-orbit refusal: suite only | C (MEM) | PASS |
| D86o | S | suite; G: backend table at 23 threads with the record in force (no per-spec record): every Float64 spec resolves to the CPU, seventh reason for tg, pr, pa, bk and sixth for st and fq; with the per-spec records the rows follow them (gpu where `gpu_faster`); a record with another `spec_sha256` gives the CPU, the seventh reason and the note | G (TABLE, RECORD) | PASS |
| D86p | S, G | suite on the KA CPU backend; G: tg, pr, pa, bk (373 or 372 Float64 items each) compile and run on the device with provenance gpu; fq and st give the sixth reason; Double64 of all four falls back with reason 3 ("kernel does not compile for the element type"); finite records and the health rule hold. O21 below. Injected out-of-memory and compile faults: suite only | G | PASS (4 FALLBACK rows are the Double64 passes) |
| D86q | S | suite; V (I24): campaign B (23 threads, capacity 1000) against T4 (4 threads, capacity 7, own process) and against C (own process, capacity from the free memory): 0 differing datasets or attributes outside times, host fields, `chunk_size`, `launches_per_chunk` | C (READER) | PASS |
| D86r | S | suite only; V: the injected failure fails its chunk twice and no other chunk fails in 4 campaigns | C (CAMPAIGN) | (S), PASS |
| D86s | S, V | package read-back equals the written file (suite); V: the independent reader read 352 of 352 files in each of 4 stores with 0 failed checks (R15 root attributes, dtypes, shapes, attribute lists of every dataset, provenance groups, pooled group, no `lyap`, E31); the module's `read_chunk` equals it on 47812 attributes and 33977 datasets (sha-256 of the raw bytes, dtype, shape), 0 differ, in each store; `h5dump` agrees on 49 datasets | `-reader.jl`; C | PASS |
| D86t | S | suite (20 kill points, 10 resumes); V: killed after 70 chunks and during the data of chunk 71 (one `.tmp`, no final file), resumed on 23 threads: 0 differences against the uninterrupted campaign (I14); a rerun on the finished store runs 0 chunks and modifies 0 files | C | PASS |
| D86u | S | suite; V: R7 keys, chunk and arm entries; `indicators` table (kind, text, hash, chain, budget, bytes terms, windows, run length, frame source, direction rule, observable names); `indicator_flags` of 327 chunks equal the recount from the files; `not_applicable` lists 883 for fq and bk; `transform_seconds` on 47 frequency chunks; per-arm per-spec cost lines; same in manifest A | C (MANIFEST) | PASS |
| D86v | S | suite; V: CPU rerun of 20 items per arm (9 arms, 180 items): 180 bitwise equal; pilot report: window kinds at 2000, 2001, others at 300, bytes(orbit), capacity, flag shares, orbit pass and transform seconds for frequency, "not study data" | C (RERUN, PILOT) | PASS |
| D86w | S, V | suite; V: no p2 preregistration exists, so the search is for the reduced block's values (250, 1000, 2001, 4000, 8.0, 1e-12, the budget lines) as numeric tokens in `src/Experiments`: 4 hits, all explained (`PILOT_MAX_DRAWS`, `WARM_STEPS`, `PAR_STEPS_PER_TASK`, the placeholder tolerance of the allocation probe); 0 defaults of a spec key; the only code line in `run_indicator*.jl` and `run_block.jl` with `threshold` is the monitor threshold of R5; of 101 distinct attribute and dataset names written, the only matches for threshold, detect, cross, class, sticky are `diagnostics@threshold` (R5) and the panel labels `class_label`, `design_class` | `-source.jl`; C | PASS |
| D86x | S | suite only | `ex_d86_guards.jl` | (S) |
| D86y | V, G, [full] | below | C, G | PASS (with F1 worked around) |
| D86z | stage 05 | not run here | - | open |

D86y in detail:

| Pass condition | Number |
|----------------|--------|
| Campaign completes or stops only by the injected causes | 4 campaigns (A resumed, B, C, T4): status `complete with failed chunks`, failed chunks exactly [25] (the injected one); the device campaigns likewise |
| Manifest validates against R7 and R19 | PASS for A and B (above) |
| Every file passes O20 | 352 of 352 files in each of A, B, C, T4 and in the device stores GC and G2: 0 failed checks |
| Resume equals an uninterrupted second campaign (I14) | 0 differences |
| CPU rerun of D86v on 20 items | 180 of 180 bitwise |
| O21 on the device chains | below |
| Backend table follows the records at the run's thread count; CPU-only kinds record the sixth reason | table and campaigns G1 (auto with the per-spec records: pair on the device, the others on the CPU) and G2 (gpu request: tg, pr, pa, bk on the device; st, fq on the CPU with "reducer is CPU-only"; Double64 on the CPU); one tg chunk of G2 shows the reason "NumericalHealthError raised only on the GPU" (the injected failure, rerun on the CPU as designed) |
| Wall time below 15 minutes on 23 threads | 108.2 s for `run_campaign` in its own process, 127.1 s including Julia start, package load and compilation (limit 900 s); warm campaign in the main process 32.7 s |
| Peak host memory below 8 GB | 1.01 GiB (Sys.maxrss of that process, 0.55 GiB before the campaign) |
| Study data folder | no entry changed (see the note above) |

O21 on the device (log G, Float64, same orbits as the campaign):

| Chain | Form | Result |
|-------|------|--------|
| phase advance (372 orbits, 2 windows) | same samples, algebra only (naff-fma D82 form): device reducer against `weighted_phase_frequency` on the samples of the device orbit (trace kernel), max / (2 B_E2) | 0.0023 (1468 values, 20 undefined on both sides); the trace's final state equals the recorder's orbit point bitwise for 372 of 372 |
| birkhoff (372 orbits, 4 observables, 2 windows) | weighted-birkhoff 13(b): device reducer against `weighted_birkhoff` on the samples the device computed (b1) and on the CPU evaluation at the device states (b2) | b1: 0.25 (avg, in 4 eps S) and 0.14 (delta, in 8 eps S); b2: 0.25 and 0.14; gates 1 |
| tangent chain FLI, MEGNO, SALI (372 orbits) | recorder device columns at N = 4000 against the verified device path (`propagate`, same chain, same frames, backend gpu) | 372 of 372 bitwise equal, largest difference 0 |
| pair RLI, MEGNO (372 orbits, N = 300) | same, on the PairMap with the stored xi_0 and block frame | 372 of 372 bitwise equal |
| SALI of the chain, regular class (K = 0.5, 24 orbits, N = 4000) | sali criterion 18(a): |L - L_ref| / (8 eps (Psi1 + Psi2)) and |SALI - SALI_ref| / B_SALI against BigFloat(256) at the device's own orbit points | 0.00069 and 0.00086; gates 1 |

The FLI, MEGNO and RLI columns are gated through the bitwise equality with the `propagate` device path, whose own gates (fli-ofli 16, megno 15, rli 15(b)) are in `gpu-fli-verify-log.txt` (23 PASS, 0 FAIL), `gpu-megno-verify-log.txt` (32, 0), `gpu-rli-verify-log.txt` (23, 0); the BigFloat reference gate was rerun for SALI only (`gpu-sali-d68b-verify-log.txt`, 15 PASS, 0 FAIL). The orbit-level gap between device and CPU recordings is reported and not gated (R20 item 2): the device and the CPU records are bitwise equal for 60 of 373 tangent, 41 of 373 pair, 13 of 373 phase advance and 12 of 372 birkhoff orbits (the rest differ with the chaotic divergence of the two orbits; largest stored difference 5546, 80, 2.9 and 865 in tangent, pair, phase advance, birkhoff; the figure shows the per-orbit gap); stretching and frequency (CPU only) are bitwise equal between the device campaign and the CPU campaign for every orbit.

## Deviations accepted at the D86 stage 03 checkpoint (D87): did each hold?

| # | Deviation | Evidence that it held |
|---|-----------|-----------------------|
| 1 | window kinds step one handoff at a time and reset the deviation column | the phase-advance, birkhoff and frequency tables equal the pure functions on the plain-iteration series bitwise (158 rows); the same loop is in the device kernel and the device tables pass O21 |
| 2 | E1 buffer per orbit | the campaign with 1.01 GiB peak memory; `bytes(orbit)` still charges B (39216 for fq in the reduced block: B = 32000; 136000 in the worked example); capacity equals E28 exactly |
| 3 | pilot horizon ignores window-only arms | pilot rows at 2000, 2001 and 300; the block's pilot horizon (300) is below every arm, so the exemption itself is not exercised here (suite D86m) |
| 4 | `launches_per_chunk` of a CPU indicator chunk is 1 | all 328 indicator chunks of the capacity-7 campaign T4 carry 1; the L1 chunks carry ceil(orbits / 7) |
| 5 | capacity not fitting one orbit creates the empty store and raises | suite only |
| 6 | `first_full_step` is 0 when never reached | not exercised (the reduced block reaches it at 200 handoffs); suite |
| 7 | window kinds record the kind as members and chain; provenance of the first item | the reader check `members == [kind]` and the provenance group named by the kind pass on 352 files; the range is in `scalars@range` for every stretching orbit |
| 8 | count matrices of a chunk stay in memory | not measured separately; peak memory 1.01 GiB for 70 stretching chunks of 60 orbits at 100 bins |
| 9 | `_same` compares attributes with `isequal` | the failed pair chunk (O18, y5long) is read back; not compared against a rerun here |
| 10 | frame preparation in `_open_session` | the campaigns start without refusal; dependent-frame refusal: suite only |
| 11 | cost cap reads the spec budget; `specs` only for arms with indicators | R19 check: every indicator arm has the per-spec lines; the cap is disabled in the reduced block (`first_chunk_cost_factor = 1e6`), the 2.01 and 1.99 cases are suite |
| lifted | lifted-mode frequency bound 16 eps + max abs q eps | O18, 96 phase-advance values within it, worst ratio 0.035 |

## Not verified

- Injected device faults (out of memory, compile), the 20-point kill matrix, the atomic-write hooks, the dependent frames, the one-orbit refusal and the cost-cap thresholds: suite only.
- The independent reader is HDF5.jl in a separate process, not h5py (not installed); the second tool is `h5dump`.
- FLI, MEGNO and RLI device columns against a BigFloat reference at the device orbit were not rerun (bitwise equality with the verified `propagate` device path and the logs of those criteria).
- Stretching and birkhoff closed forms (Y5, K_std = 0 weighted Birkhoff criterion 3 at window 1000) and the p2 values (no preregistration yet) are not covered by the scripts. Detection times and sticky definitions are out of scope.
- The command line entry (`bin/`, M3) was not run: the campaigns call `run_campaign` and `run_pilot` directly.
- Device cost records are of 1024 repeated orbits and vary between takes (tangent 1.4 to 3.1 times, phase advance 0.6 to 1.1); `auto` for phase advance is not stable and the real study should take its own record.
- D86z (documents) is stage 05.

## Independent review

Reviewer: separate Sonnet agent, spec and code only (`experiments-d86-review.md`): 0 CRITICAL, 0 HIGH, 3 MEDIUM, 9 LOW, all by reading; nothing was run. I confirmed R1 and R2 in the code myself.

| ID | Sev | Status |
|----|-----|--------|
| R1 | MEDIUM | fixed (above) |
| R2 | MEDIUM | fixed as D88 (above); the finding was: the pooled stretching counts of a chunk are summed over orbits whose derived range can differ (chunks "controls" and "h4" hold different K and the Y5 matrix) and the pool stores no range. Recommended: store the range as an attribute of `pooled/<id>` and write no `counts` when the orbits of the chunk do not share it (schema change to E29 and `schema/run-store.md`) |
| R3 | MEDIUM | fixed (above) |
| R4 | LOW | settled on E25 (D88); the finding was: E25 gives the tangent family and stretching no flag word but R18 item 2 lists `undersampled` and "SALI floor reached" as counted flags, and E24 asks for the SALI floor as a root attribute; the code follows E25. Owner chooses the reading |
| R5 | LOW | open: the plan-time check of the derived range tests the system name only; a linear control with no stretching (largest singular value within 8 eps of 1) fails at run time inside a worker and aborts the campaign |
| R6, R7, R10 | LOW | open: per-spec cost report and cost cap judge the first group of an arm only, and pool precisions; the pilot report keeps the last group when a spec has several |
| R8 | LOW | open: the write-time read-back does not read `/pooled/<id>` back |
| R9 | LOW | open: `indicator_manifest` rescans the items at every manifest write (about 0.1 s per write at 2 million items, not measured) |
| R11 | LOW | open: `_stretch_D` restates the D formula of the reducer; guarded by the bitwise test D86i |
| R12 | LOW | open: the stretching workspace is allocated per orbit, not per slot (garbage, not wrong data) |

Not reviewed by the reviewer: the reducers, the tests as evidence, the `verify/` scripts, the schema text, an independent HDF5 reader, performance numbers. Everything the reviewer checked and found correct is listed in `experiments-d86-review.md`. The review file has LF line endings (git normalises them).

Audit: no open CRITICAL or HIGH finding and no open MEDIUM. Open: the LOW findings R5 to R12 and the design note F2 (the p2 preregistration must keep the pair horizon below about 700 for Y5).
