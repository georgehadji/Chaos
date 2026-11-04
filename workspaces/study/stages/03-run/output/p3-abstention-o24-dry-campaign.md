# p3-abstention: O24 dry campaign (gate before the A1 freeze run)

Date 2026-11-04. Not test data: the panel is a small fixture panel, not the test panel, and nothing is written under `data/` (the gate of `../references/round-loop.md`, owner decision 2026-11-03, check O24 of D95). Code at commit e5d97f8 (Amendment A2 committed, preregistration hash f614a38449f611a42c56a0aa79bbc2d106dfaab679d1972a2b5349011df1a072). No Experiments or Analysis function is patched or redefined.

## Fixture

`verify/p3-abstention-o24-inputs.jl` reads the panel input block and the run block from section "Amendment A2" of the committed preregistration and shrinks them so that the toy p2 and P1 inputs of the Analysis test helpers can serve the counts. Unchanged from A2: the six arms (ref each_round; ind, ind_pa1, d64, d64_pa1, chk closure), their orbit sets, methods and precisions, the four indicator specs and their budget lines, the health table, the budgets (BigFloat 341296.4 ns, item 6 of A2), `controls = "none"`, the panel keys of the P1 rule and of D94. Shrunk: panel seed 777, 4 K (0, 0.5, K_c, 2 pi), 400 screening points, 3 rounds, tiny quotas, D fraction 0.5, `db_size` 5 (46 orbits); horizons 80000 (ref), 20480 and 20481 (indicator arms), 1000 (chk), windows 20 to 10240, `k19.threads` 4, pilot horizon 1000. Two-file seal (item 4 of A2): build copy X (header 18d89ac6...e735, list holds 64 zeros) for the panel, sealed copy Y (header d5efeec1...963d, list holds the X hash) for the runs.

## Result

| Check | Evidence | Result |
|-------|----------|--------|
| Library rounds | `verify/experiments-d96.jl lib`: rounds 1 to 3 status complete (26, 22, 20 chunks), indicator table empty before closure; `count_round` continue, continue, close (15, 30, 44 orbits) | pass |
| Closure and resume | closure run 282 chunks, indicator table bk1, fq1, pa1, tg; resume 0 run, 350 skipped | pass |
| CLI rounds | `verify/experiments-d96-cli.sh`: `bin/chaoslab.jl run --round-log` exit 0 at rounds 1 to 3, closure and resume; 350 chunk files in both stores | pass |
| Runner accepts the log | CLI round log byte-identical to the library log (SHA-256 b819407fdc4e4529b3088fd2b307256080b2efd3282c3508129630cb1d22c8f7) | pass |
| Analysis on the closed log | `verify/experiments-d96.jl analyze`: `count_round` action closed; `analyze` status ok, no refusal (HA, HB, HC not decidable on the toy) | pass |
| Reproducible | a second build from the script gives a byte-identical panel, the same X and Y hashes and the same round log hash | pass |

Wall time of the library campaign 53 s at 4 threads; panel build 167 s.

## Not covered

Real horizons (1e7, 102400, 1e6) and the real panel sizes; the first-chunk cost cap of chk at the real BigFloat cost; K19 at 23 threads; failed chunks (covered by D97 stage 04, O29); the real calibration record and test-side table (they come from the freeze run and A1).
