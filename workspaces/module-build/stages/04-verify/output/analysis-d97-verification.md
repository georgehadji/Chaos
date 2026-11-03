# Verification: analysis, Amendment D97 (stage 04)

Date: 2026-11-03, Julia 1.13.2, 4 threads, code at cc38b10. Spec: `analysis-spec.md`, "Amendment D97" (E74 to E82, I47 to I51, O28 to O30, criteria D97a to D97k, P117 to P121). Test plan: `analysis-test-plan.md`, "Amendment D97". Notes: `analysis-d97-impl-notes.md`.

Verify scripts (all in `verify/`, default project unless stated; none writes under `data/`, each refuses paths inside `/data/`):

- `verify/experiments-d96.jl` (changed, D96 behaviour kept): the patch include is conditional on `<fixture>/patches/<file>` existing; new mode `analyze` (count_round on the closed log, then `Analysis.analyze`); environment `FAULTS=1` loads a `mutate` hook that fails chunks twice (seven chunks, six arms). Modes `lib` and `count` unchanged.
- `verify/experiments-d96-cli.sh` (unchanged): CLI `run --round-log` by rounds with the library `count_round` between runs.
- `verify/analysis-d97-cli.sh` (new): `bin/chaoslab.jl run --round-log`, `count-round` between the runs, closure, resume, `analyze`; round log and outputs against the library ones.
- `verify/analysis-d97-oracle.py`, `verify/analysis-d97-compare.py` (new): the independent reader O20 adapted to the store of the runner. They supersede `analysis-p3t-oracle.py` and `analysis-p3t-compare.py` for a runner store; the D95 files are kept unchanged. The export script `verify/analysis-p3t-export.jl` is reused unchanged (verify environment).
- `verify/analysis-d97-runner.jl` (new): the runner's own log reader and checker on the log written by `count_round` and on the 39 one-byte alterations of `verify/experiments-d94-log.py`.

Raw output: `analysis-d97-verify-log.txt` (sections A to G).

Fixture: scratchpad `d97v/fixture`, the D96 fixture without `patches/` (13 files, SHA-256 list of `d96v/fixture-sha256.txt` checked, all OK). Study p3-abstention, panel seed 777, 46 orbits (2 gate, 44 with K > 0), 3 rounds, arms ref, ind, ind_pa1, d64, d64_pa1, chk, 350 chunks, closed at round 3. Nothing under `data/` read or written. No Analysis patch is loaded anywhere.

Two campaigns: clean (`d97v/out`), and with chunks failed twice (`d97v/outf3`): ref L1 `k3:T:r2`, ind tg `k2:C:r2`, ind_pa1 pa1 `k3:C:r2`, d64 tg `k3:C:r1`, d64_pa1 pa1 `k2:C:r1`, chk L5 Double64 `k2:C:r2`, chk L5 BigFloat `k3:C:r2`.

## Evidence per criterion

| Criterion | Evidence | Result |
|---|---|---|
| O28 library | clean campaign, no patch: rounds 1 to 3 (26, 22, 20 chunks), `count_round` continue, continue, close (15, 30, 44 orbits); closure 282 chunks; resume 0 run, 350 skipped; `count_round` on the closed log gives `closed`; `analyze` status ok | PASS |
| O28 CLI | `experiments-d96-cli.sh`: exit 0 at rounds 1 to 3, closure, resume; CLI log byte-identical to the library log. `analysis-d97-cli.sh`: `run`, `count-round` (exit 0 after each round and on the closed log), `analyze` exit 0; CLI `count-round` log equals the library log (60c708b5a15a...) | PASS |
| O28 log, independent writer | `experiments-d94-log.py` with the counts of the independent reader: bytes equal the `count_round` log in both campaigns (clean 60c708b5..., faults 93a95040...); the oracle also recomputes each entry equal to the file | PASS |
| O28 runner side | runner `read_round_log` and `check_round_log` accept the log; of 39 one-byte alterations 32 refused, 7 accepted (count lines of the last entry, as D94 O26, D95 O24) | PASS |
| O28 independent reader, clean | oracle checks 25 of 25; comparison 956 quantities, 956 equal, 0 differ | PASS |
| O29 | status `complete with failed chunks`; `analyze` ok; oracle 25 of 25; comparison 951 of 951 equal. Ref `k3:T:r2` excludes 24, 25 (log `excluded` 0, 2, 2); ind and ind_pa1 give "no verdict data" for 21, 22, 27, 28; d64 and d64_pa1 take 13 and 6 out of D (abs(D) 23 to 19); chk Double64 takes 22 out of D_B (abs(D_B) 4), chk BigFloat leaves 28 in D_B with `n.a.` in T3; T11 lists seven chunks with the orbits of their cells, equal to the oracle | PASS |
| D97a (S) | tests D97a; real store: every chunk key ends in `k<i>:<class>:r<round>`, counts and digest recomputed from these cells equal the file, I51 check passes | pass |
| D97b (S) | tests D97b; real: every verdict item read is at N_c (pa1 at N_c + 1), verdict outputs equal the module, abs(D) = 23 | pass |
| D97c (S) | tests D97c; real: `complete` and `complete with failed chunks` (7 failed entries) both accepted | pass |
| D97d (S, V) | tests D97d; real: failed chunks in each of the six arms give the exclusions of E77, equal to the oracle; refusal of a missing item of an unlisted chunk carried by tests only (see review items 1 and 2) | pass |
| D97e (S) | tests D97e; real: T10 measured costs equal the oracle to 12 digits for tg, fq1, bk1, pa1, L1 (clean: 160.811, 2458.74, 301.300, 189.728, 79.8726); sensitivity (iii) decidable | pass |
| D97f (S) | tests D97f; real: T3 header without "labels agree", omission text, five rows with two ladder flips (28, 31) equal to the oracle; one row with the missing BigFloat cell in the faults campaign; chk horizon 1000 equals N_B | pass |
| D97g (S) | tests only; the fixture closes at its last round, so no round above the closing round exists; the oracle's T2 and counts run to the closing round and equal the module | pass (S) |
| D97h (S) | tests D97h; real: the manifest names six arms, P114 and P115 texts in the tables file | pass |
| D97i (S) | `an_p3t_shape.jl` and fixture tests, carried by the suite | pass (S) |
| D97j, I47 (S) | whole package suite of stage 03: 1404780 pass, 0 fail, 1 broken (existing); not rerun here | pass (S) |
| D97k (V) | O28 and O29 pass; CLI log byte-identical to the library log; CLI `analyze` outputs byte-identical to the library output on the same store in both campaigns; independent implementation equals the module (956 and 951 quantities) | PASS |
| I48 | counts of each log entry equal the independent reader's, in both campaigns, and the D96 patched-run counts (only `ref_digest` and `prev_sha256` differ) | pass |
| I51 | oracle check: every orbit up to the closing round is in exactly one ref cell (the module does not check it, review item 7) | pass (oracle) |
| O30 (S), I49, I50 | tests only; I49 not measured here | not rerun |

## Independent review

Fresh Sonnet subagent, given the spec and the diff `f385423..cc38b10 -- src/Analysis`, not the notes or the test plan; read only. No CRITICAL or HIGH finding. E74, E75, E76, E78, E79, E80, E81 confirmed; E77 partly (items 1, 2).

| # | Finding | Severity | Disposition |
|---|---------|----------|-------------|
| 1 | `p3t_tables.jl` `_p3t_t3`, `p3t_data.jl` `p3t_dB`: D_B is formed from chk Double64 items only; a missing chk item (either precision) of a chunk that did not fail shows `n.a.` where E77 says a missing item of an unlisted chunk is refused | MEDIUM | fixed in the rework: refused when the orbit carries the other chk item and the chunk is not failed (`an_p3t_d97r.jl`) |
| 2 | `p3t_data.jl` `_p3t_double`: a D orbit missing a Double64 item leaves D silently whether or not its d64 or d64_pa1 chunk failed (E77: unlisted chunk, refusal); no separate count of orbits that left D by failure | MEDIUM | fixed in the rework: refused for all four Double64 items of an unlisted chunk; a failed chunk leaves D, counted in abs(D) (E70 and P104 b name no other place) |
| 3 | `p3t_store.jl:48`: the panel key `class` is not in the required keys; a panel without it gives a bare `KeyError`, not a named refusal | MEDIUM | fixed in the rework: `class` required, named refusal |
| 4 | `p3t_run.jl:49-55`: the runner writes `measured_ns_per_orbit_step = 0.0` for a spec with zero steps; Analysis accepts 0.0, so E78 "zero sum of steps, not decidable naming the spec" is missed for tg, fq1, bk1, pa1 (L1 handled) | MEDIUM | fixed in the rework: 0.0 is no value, (iii) not decidable naming the spec |
| 5 | T3 is skipped without the P1 parameters, so the N_B check is not evaluated | LOW | as D95; no change |
| 6 | `_p3t_chk_items`: duplicate (id, label) after the bit count is stripped overwrite silently; Dict order picks the named mismatch | LOW | owner decision |
| 7 | `p3t_inputs.jl:56-59`: the six arm names are required in the manifest, but a chunk of an unknown arm is not refused (I51 part) | LOW | owner decision |
| 8 | `p3t_data.jl:38`: a failed closure chunk with an off-form cell is not recognised; ends in an unnamed "item missing" refusal | LOW | owner decision |
| 9 | `p3t_data.jl:41`: the drift set is over all items of all arms, now including ind_pa1, d64, d64_pa1, chk (E64 b "any item read") | LOW | owner to confirm intended |
| 10 | `p3t_tables.jl:26`: the T2 header sentence text changed (I50 allows only named lines) | LOW | owner to accept |
| 11 | `p3t_verdict_tables` duplicates the body of `p3_verdicts` with a horizon added | LOW | optional |
| 12 | `p3t_run.jl:36-38`: chunk entry fields indexed directly; a missing field gives a bare `KeyError` | LOW | optional |
| 13 | size, I47, exact arithmetic: no function over 50 lines, no file over 400; no non-p3t file changed; ladder and lambda stay `Rational{BigInt}` | none | confirmed |

## Rework after the review (2026-11-03)

Owner decision: fix items 1 to 4. Tests `test/analysis/an_p3t_d97r.jl` (7 testsets, 21 assertions; RED 14 fail, 7 pass guards), code in `p3t_data.jl`, `p3t_store.jl`, `p3t_run.jl` (+33/-8). Analysis subset GREEN; whole package suite `Chaoslab | 1404801 pass, 1 broken, 1404802 total` (broken: existing weighted-birkhoff). Oracle rerun on both campaigns (scratchpad `d97v/rework`): 956 and 951 quantities equal, 0 differ; CLI outputs byte-identical to the library. Readings: the count of orbits that leave D by a failed chunk stays in abs(D); a missing chk item is detectable only when the orbit carries the other chk item (D_B is drawn at closure and not recorded elsewhere).

## Audit

| Check | Result |
|---|---|
| Evidence | every criterion has a number or run line; D97g, O30, I49 carried by tests only |
| Review closed | no open CRITICAL, HIGH or MEDIUM (items 1 to 4 fixed in the rework); LOW items 5 to 12 open for the owner |
| Analysis defect found by the oracle | none |

## Notes for the owner

- Independence: the adapted reader keeps the structure of the D95 independent file, written from the spec (E74 to E81) and the schema. After the first fault run it read one docstring of `src/Analysis/p3t_data.jl` (D_B); spec line 1042 says the same.
- Spec point, E77 and E79 for chk: D_B is the set of orbits with an L5 Double64 chk item. An orbit whose Double64 chk chunk failed leaves D_B and has no T3 row; "its T3 row shows the missing cell" holds when the BigFloat chunk fails (orbit 28, `n.a.`). The module follows the D_B definition (impl notes, deviation 4) and the stage 02 reading approved by the owner.
- E77 "the orbits leave D and are counted": the loss shows only through abs(D) and `hc.set_size`; no separate count (review item 2).
- The fixture gives HA, HB, HC "not decidable", T_s empty, W2 of 10 (7); HB coverage of an empty set is 0 (price 1), as in D95. Statistics are checked on small sets only.
- Not reached by the real fixture: a round above the closing round, T_s > 0, early close, class S, drift flags, paused arm, kill and resume. D97g stays with the tests.
- The first clean `analysis-d97-cli.sh` run printed two DIFFERENT lines comparing two different campaigns (not byte comparable); the final script analyses the library store through the CLI and compares those bytes (PASS, both campaigns).
- Chunk hashes change between reruns, so the clean log SHA-256 differs from the D96 report (4c96292f...); counts are equal.
- No plot: the amendment has no convergence or order quantity.
