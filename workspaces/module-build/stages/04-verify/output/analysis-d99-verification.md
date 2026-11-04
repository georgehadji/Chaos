# Verification: analysis, Amendment D99 (stage 04)

Date: 2026-11-04, Julia 1.13.2, code at 15f22d3 ("after"), pre-amendment code `git archive d5946d2` in scratch ("before"; its `src` and `bin` equal those of b2c1df0, the code of the committed p3 run). Spec: `analysis-spec.md`, "Amendment D99" (E83 to E86, I52, I53, O31, O32, D99a to D99i, P122 to P125). Test plan: `analysis-test-plan.md`, "D99". Notes: `analysis-d99-impl-notes.md`. Raw output: `analysis-d99-verify-log.txt` (sections A to G).

Nothing under `data/` written: `data/p3-abstention` read by `analyze` only (no file newer than the start of this stage; `git status` of `data/` and `workspaces/` clean). The committed study outputs are not written to: every run went to an empty scratch folder.

Scripts:

- `verify/analysis-d99-chain.jl` (new, default project): D99g. A study layout with the default file names; `calibrate`, `run`, `count-round` and `analyze` through their CLI parsers in process. Options only where the default is outside the layout: the stores (default under `data/`) and the manifest and round log of `run` (default the study folders of the repository).
- Scratch `d99v/build.jl` (test helpers, as the D96 builder): the toy p2 and P1 stores (`p3_cached("d96")`, n_s = 1000) and the layout `d99v/study/stages`. The five p2 and P1 input files it produces have the same bytes as the D96 fixture (log A).

## Evidence per criterion

| Criterion | Evidence | Result |
|---|---|---|
| D99a to D99f (S) | `an_p3t_d99.jl` 81 of 81 (stage 03); whole suite 1 404 893 pass, 0 fail, 1 existing broken | pass (S) |
| D99g (V) | D96 fixture layout, 3 rounds, closure: the default folder of `calibrate` equals the record folder of `count-round` and `analyze` (same normalized path); `run` exit 0 at rounds 1 to 3 and closure; `count-round` exit 0 after each round (continue, continue, close) and on the closed log; `analyze` exit 0, outputs in the default folder, record fingerprint bfd14662e98b equals the record in the default folder, unchanged after the chain (log B). Deviation: `calibrate` itself does not run on this fixture (notes 1) | PASS, with deviation |
| D99h, O31 (V) | real p3 inputs (store `data/p3-abstention`, read only). Before reproduces the committed tables except the line "results SHA-256" (git fingerprint `unknown` from the archive, log C). Before vs after: tables differ in exactly 7 lines, all named by I52 (T2 row D target `n.a.`, T3 introduction, HA decision line, added T5 sentence, T5 row a(T_s) - a(S2) `fails` to `ok`, added T8 line, P120) plus the line "results SHA-256", which follows from the machine-file changes E84 and E85 (log D). Machine file: 14 043 leaves, 4 differ: `git_commit`, `git_clean`, `hc.target_failures` (new, 45), `readings[21]` (P120) (log E). 0 numeric cells of T1 to T11 differ; decisions HA, HB, HC equal; data of F2 and F3 equal (no line of their sections differs); the round log is not written by `analyze` (unchanged in git) | PASS |
| D99c on the real run | after: "A2" 0, "study amendment" 0 in the tables file and machine file; "deviation of the analysis" 2 and 1 (log F) | pass |
| D99e on the real run | `HC: ... failures 76 ...`, then `Target failures n_f: 45; achieved failures: 76`; underpowered false (log F) | pass |
| D99i, I53 (S, V) | p1 and p2 code untouched (`git diff d5946d2 15f22d3` names 3 p3 files); p1, p2 and `calibrate` golden tests pass in the suite. Real `calibrate --out <scratch>` at 15f22d3: record 2 293 leaves, 2 differ (`git_commit`, `preregistration`: the committed record was written before amendment A1 was added to the preregistration); calibration tables differ only in the record digest, the preregistration hash and the git commit lines (log G) | pass |
| I52 | D99h above | pass |

## Independent review

Fresh Sonnet subagent, given the spec and `git diff d5946d2 15f22d3 -- src/Analysis`, not the notes or the test plan; read only. E84, E85, E86, I52 confirmed (E86: the three default paths normalize to the same string); E83 confirmed for a decided HA (the row label uses the same exact comparisons as `_p3t_ha_decide`).

| # | Finding | Severity | Disposition |
|---|---------|----------|-------------|
| 1 | `p3t_tables.jl:101`: when HA is not decidable, `failed` is empty, so the decision line prints "Status of item 3: ok" while a row with an undefined bound reads "fails". Confirmed on the D96 fixture: "HA: not decidable. Status of item 3: ok." with row a(T_s) - a(S2) "fails" (lower bound undefined). Not reached on the committed run (HA decided) | MEDIUM | fixed in the rework (below): the item status is "ok" when both row statuses are "ok", else "fails" |
| 2 | `_p3t_q` calls `p3_rat` on any non-Rational; a point given as text would throw. Points come from `_sub` as Rational or nothing | LOW | no change: unreachable |
| 3 | sizes: no changed function over 50 lines, files 197, 122, 101 lines | none | confirmed |

## Rework after the review (2026-11-04)

Owner decision: "Fix, tests first". Test `test/analysis/an_p3t_d99r.jl` (three hand fixtures with HA not decidable and an empty failed list; RED 10 pass, 2 fail, 0 error; test plan "Rework after stage 04 review (2026-11-04)"). Code: `_p3t_t5` in `p3t_tables.jl` (+2/-1) takes the status of item 3 in the decision line from its two rows (ok when both are ok), which equals the failed list whenever HA is decided. GREEN: p3 subset 956 of 956 (d99r 12 of 12); whole package suite 1 404 905 pass, 0 fail, 1 existing broken (weighted-birkhoff); quality 12 of 12, layering 2 of 2, analysis 2 763. Real p3 run (log H): the tables equal the first after run except the line "results SHA-256" (machine file: only `git_clean`, the tree being dirty during the run), so D99h stands. D96 fixture: "HA: not decidable. Status of item 3: fails." with rows ok, fails.

## Audit

| Check | Result |
|---|---|
| Evidence | every criterion has a number or a run line; D99a to D99f by tests |
| Review closed | no open finding: item 1 fixed in the rework, item 2 unreachable |

## Notes for the owner

1. D99g deviation. The record of the D96 fixture is built by the test helpers (`p3t_record_nt`), and its SHA-256 is pinned in the test-side table (`calibration_record_sha256`). `calibrate` on these toy inputs refuses (`recorded step 1000 is missing in chunk-00001.h5, orbit 1, spec L1 (lambda_1)`), and a record it wrote would not match the pinned hash. So the chain script resolves the default folder of `calibrate` (the function `calibrate` uses) and places the fixture record there; from there on every path is a default. The link "calibrate writes where count-round and analyze read" is shown by that path equality, by D99f, and by the real `calibrate` writing the same record bytes as committed apart from fingerprints (D99i). A full chain with `calibrate` writing would need a fixture whose test-side table is built from a `calibrate` record.
2. The log of section B was made before one label of the script was corrected ("the record calibrate wrote" became "the record in the default folder of calibrate"); the check is the same.
3. `run` takes `--manifest` because its default is the study folder of the repository (`_root()/workspaces/study/stages/03-run/output`), not relative to the preregistration as the analysis defaults are. Not in the scope of D99.
4. No plot: the amendment has no convergence or order quantity.
