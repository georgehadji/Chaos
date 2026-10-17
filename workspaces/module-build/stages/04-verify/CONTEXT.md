# Stage 04: Verify

Prove the module against oracles and an independent review.

## Inputs

| Source | File/Location | Section/Scope | Why |
|--------|--------------|---------------|-----|
| Stage 01 | `../01-spec/output/[module]-spec.md` | "Acceptance criteria", "Oracle", and the criteria and corrections of each amendment section | what to prove |
| Stage 02 | `../02-tests/output/[module]-test-plan.md` | criterion-to-test table | which test files carry each criterion |
| Stage 03 | `../03-implement/output/[module]-impl-notes.md` | Full file | deviations to check |
| Code | `/src/<Module>/`, `/test/<module>/` | the files the test plan names | what is verified |
| Verify | `/verify/<module>.jl`, `/verify/gpu.jl` | the module's script; for [local] criteria the matching `gpu.jl` sections | oracle and device runs (D27, D29) |
| Config | `/_config/environment.md` | "Local GPU" | the device of every [local] criterion |

## Process

1. Run `verify/<module>.jl` (oracle environment `verify/`) at each precision the spec names, and the [local] criteria with `verify/gpu.jl` on the device of `_config/environment.md`; make the convergence or order plot.
2. Independent review (fresh Sonnet subagent, given the spec and the code, not the notes): math versus code, numerical traps.
3. Run the audit, then write the verification report: evidence for each acceptance criterion, review findings and their fixes.
4. **[Checkpoint]** Present the report.

## Checkpoints

| After Step | Agent Presents | Human Decides |
|------------|---------------|---------------|
| 3 | Evidence per criterion, plot, review findings | Accept, or send back to 01, 02 or 03 |

## Audit

| Check | Pass Condition |
|-------|---------------|
| Evidence | Every acceptance criterion has a number or a plot |
| Review closed | No open CRITICAL or HIGH review finding |

## Outputs

| Artifact | Location | Format |
|----------|----------|--------|
| Verification report | `output/[module]-verification.md` | Evidence table, review findings (read by 05) |
| Figures | `output/[module]-*.png` | Convergence and oracle plots |
| Verify log | `output/[module]-verify-log.txt` | Raw output of `verify/<module>.jl` (read by the report) |
| GPU logs | `output/gpu-*-verify-log.txt` | Raw output of `verify/gpu.jl`, one file per run, shared by the modules it covers |
