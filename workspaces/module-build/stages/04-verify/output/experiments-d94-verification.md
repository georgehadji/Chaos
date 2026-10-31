# Verification: experiments, Amendment D94 (stage 04)

Date: 2026-10-31 to 2026-11-01, Julia 1.13.2, code at ab61826 plus the review fixes below (committed with this report). Spec: `experiments-spec.md`, "Amendment D94" (E32 to E43, R22 to R27, I28 to I37, O23 to O26, criteria D94a to D94v). Test plan: `experiments-test-plan.md`, "Amendment D94". Notes: `experiments-d94-impl-notes.md`.

Verify scripts:
- `verify/experiments-d94.jl` with `experiments-d94-panels.jl`, `experiments-d94-screen.jl`, `experiments-d94-review.jl` and the input `experiments-d94-reduced.toml`: D94b [full], the reduced real panel, O24, O25, D94t source search, process-level thread counts, the O26 store, review probes REV, O26 runner side.
- `verify/experiments-d94-oracle.py` (O23) and `verify/experiments-d94-log.py` (O26 writer): independent author. It read only the spec, `schema/run-store.md` and `verify/*.py`; it read no file of `src/`, `test/`, `bench/` or the stage 02 and 03 folders.

Raw output:
- `experiments-d94-verify-log.txt`: Julia side before the fixes, 49 PASS, 0 FAIL.
- `experiments-d94-verify-log-fixed.txt`: after the fixes; REV, O26R, THREADS, P1PLAN, D94 script and suite summary.
- `experiments-d94-oracle-log.txt`: O23 on fixture A, the reduced panel and the three thread builds; every run 0 FAIL.

All run artifacts are under the session scratchpad (`d94v/`); nothing under `data/` was read or written.

## Evidence per criterion

| Criterion | Evidence | Result |
|-----------|----------|--------|
| D94a, D94c to D94p, D94t (S), D94u | package suite after the fixes: whole package 1403999 pass, 0 fail, 1 broken (existing, weighted-birkhoff); D94 files 2390 pass; quality 12 (Aqua, JET clean) | pass |
| D94b (S) | golden canonical SHA-256 of four lattice panels and the D86b block, plan and chunk hashes (suite) | pass |
| D94b [full] panel | P1 panel rebuilt under the preregistration of commit 3991155 (header hash 650b3558, the one the panel was built under) equals `p1-le-benchmark-panel.toml` byte for byte except `git_sha`, `package_version` and `prereg.path`; under the current amended preregistration the extra differences are only the prereg header, computed hash and amendment list | pass |
| D94b [full] plan | P1 plan hash `6408ce04...` equal at 7f70aa7 and at the fixed code; 1363 chunk ordinals, keys and hashes, 15569 items, 0 differing; keys equal the accepted run's manifest (its recorded plan hash cannot be recomputed: the manifest lacks the `code_version` of that run) | pass |
| D94q | O23 (independent Python) on the reduced real panel (real K list, seed 20261030, p3 quota and fraction tables read from the preregistration prose, `db_size` 100, `screen_points` 4096, `r_max` 3, no controls; 2193 orbits, 215 in D, 22 of 90 round cells empty at this size): 52 PASS, 0 differences in points (bitwise, 1193 pool orbits), point ids, ids, round cells, D flags and sizes, cell, run and orbit seeds, R and X starts, D_B for R = 1, 2, 3 (99, 100, 100); fixture A 39 PASS; generator: 1000 of 1000 words and floats of StableRNGs 1.1.0 equal (uniform = low 52 bits times 2^-52) | pass |
| D94r O24 | 200 random points per K, 1200 in all: lambda_L5 equals lambda_L1 bitwise (ratio to the O5 bound 0); class bits from L5 and from L1 equal the sidecar on all points; 0 near-threshold, 0 failures | pass |
| D94r O25 | 0 starts equal to a P1 start or one of the 2304 p2 grid starts (the p2 panel is the P1 panel file); smallest torus distance to P1 starts 1.27e-6 (orbit 1598, class R, K = 0.5, near P1 orbit 136, both drawn from the same periodic set), to the p2 grid 1.29e-3 | pass |
| D94p, process level | `julia -t 1, 4, 23`: fixture A panel, pool sidecar and plan text byte-identical (before and after the fixes); the O23 runs of the three builds pass | pass |
| O26 | independent writer's log byte-identical to the runner's log (3 entries) and its ref digests reproduced; runner side: the valid log accepted by `run_campaign` and by `cli --round-log` (exit 0); 22 of 24 one-byte alterations refused by the library and the entry point naming the entry, 0 files; 2 accepted (finding V1) | pass, 1 finding |
| D94t (V) | 0 hits of the distinctive p3 values in `src/Experiments` (small numbers in unrelated code and one docstring example listed, not scored) | pass |
| D94s | not run: the real test panel needs the amended study preregistration (P69), which does not exist. Cost estimate only: one K at `screen_points` 262144 in 110 s on 23 threads, six K about 0.18 h (bound 2 h), peak about 1.6 GiB (bound 8 GB) | not run |

## Independent code review (fresh Sonnet subagent, spec and diff, not the notes)

One CRITICAL, two HIGH, three MEDIUM, nine LOW. Owner decisions 2026-10-31: fix in code with verify probes (no new package test, as V3 of D93); H2 add a store check with a recorded test rework; M3 clarify the spec text at stage 05.

Fixed (each with a REV probe, 16 PASS):
- C1 (`run_plan.jl`): an indicator spec held only by a closure arm made every plan before closure, and the dry check plan, refuse "selected by no orbit". The study block has this shape. Now only the specs of arms planned in this phase need an item; others are checked structurally; the full check at closure.
- H1 (`run_rounds.jl`): a `D_B` arm without `release` was planned as each_round, moving ref chunk ordinals after closure (I33, E41). Now refused unless `release` is `closure`.
- H2 (`run_rounds.jl`, `run_session.jl`): the digest was checked against the plan only. `run_campaign` now also refuses an entry when a ref chunk file of its rounds is missing, not complete, or has another `chunk_hash` (root attributes only, E43), naming the entry, 0 files.
- M1: a panel with `r_max` above 1 and a block without release keys, or a `round_log` with such a block, is refused (R25 item 2).
- M2: one rule: an absent `release` means each_round in planning, text and digest (comment on `_release`).
- L4 duplicate K-index keys refused; L5 a missing `d_fraction` entry refused like `quota`; L6 bounds (`r_max` at most 1000, `screen_points` at most 2^30, fraction exponent at most 1000 before `10^ex`); L7 named refusal for a panel without rounds in `_subsets_at`; L9 invalid UTF-8 in the log named.

Recorded, not changed: M3 (pool size recorded is the count at the draw of the class, pinned reading 1; E35 says "at the start of the round"; spec note at stage 05); L1 (R24 asserts only partly: no size, D count, id gap or round range assert); L2 (manifest has no per-chunk round or "rounds run" field); L3 (`run_campaign` keywords `arms`, `up_to` beyond R26 item 3, library only); L8 (an item identical in a closure and an each_round arm is merged; cannot arise in the study design).

Test rework (owner decision H2): `test/experiments/ex_d94_run.jl`, D94o, two lines (the brief said one): the closed-plan run copies the ref chunks of the resumed store into its store before running, and the ran count is the closure chunks. No assertion weakened in kind; no other test changed.

## Findings of the verification

- V1 (O26). A count altered in the last entry of a log is accepted: no later entry chains it, and the runner does not read counts (they belong to Analysis, P81). A count altered in an earlier entry is refused at the next entry through the chain. Proposed: the test-side Analysis amendment, which writes the log, recomputes the counts of every entry from the ref store and refuses a mismatch; the Experiments guard stays as is.
- V2 (oracle). The per-orbit seed rule (E16) points to the tangent-core spec, which the oracle author was not given; it inferred the rule from the P1 panel (0 differences there and on every D94 panel). This part of O23 is not independent.
- V3 (spec text). The spec leaves open: binary64 `(2 pi) u` in `phi` and the K = 0 R start (a 256-bit evaluation differs by 1 ulp in some cases); the order of the K = 0 R fractions (denominator, then numerator); the byte encoding of the class-bit vector (T 0x01, C 0x02, S 0x04). The oracle matched the code's choices; proposed: write them into the spec at stage 05 with M3.
- V4. The REV probes were run on the fixed code only; the pre-fix failures of C1 and H1 were reproduced with the reviewer's probes before the fixes.

### Audit

| Check | Pass condition | Result |
|-------|----------------|--------|
| Evidence | Every acceptance criterion has a number or a plot | pass, except D94s (not run, needs P69; cost estimated) and D94v (stage 05) |
| Review closed | No open CRITICAL or HIGH review finding | pass (C1, H1, H2 fixed and probed) |

## Not verified (D94)

- D94s and O23 to O26 at full size on the real test panel (needs the study amendment P69).
- Byte identity across machines.
- The independent check of the per-orbit seed rule (V2).
- A count altered in the last log entry (V1).

## Owner decision at the checkpoint

Accepted by the owner on 2026-11-01, with V1 to V4 recorded as above (V1 goes to the test-side Analysis amendment; V3 and M3 go into the spec at stage 05).
