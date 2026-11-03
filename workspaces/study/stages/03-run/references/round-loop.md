# Run in rounds with a round log (Amendments D94 and D95)

For a study whose run block has arms with `release` (rounds and closure). The rules have one home: `workspaces/module-build/stages/01-spec/output/experiments-spec.md`, "Amendment D94" (R25 to R27, the runner side) and `analysis-spec.md`, "Amendment D95" (E62, E63, the counts and the log). This note only says what changes at this stage.

- Loop per round r: `bin/chaoslab.jl run ... --round-log <log>` runs the each_round arms up to the released round; then `bin/chaoslab.jl count-round <preregistration> <parameters> <test-parameters> --round-log <log> ...` verifies the whole log, counts round r from the ref arm only and appends one entry (`continue` or `close`).
- After `close`, the next `run` plans the closure arms; `analyze` runs at stage 04 on the closed log.
- Never edit the log by hand: every entry is chained and recounted from the store at each call; an altered log is refused.
- Gate before the A1 freeze run (owner decision 2026-11-03, check O24 of D95): after the study amendment P69, run a dry campaign with the real run block on a small panel that is not the test panel, alternating `run` and `count-round` by round to closure; the runner must accept the log (`run --round-log`, exit 0). No test data is produced.
