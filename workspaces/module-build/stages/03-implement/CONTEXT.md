# Stage 03: Implement

Write the Julia code until the tests pass (GREEN) and the quality gates are clean.

## Inputs

| Source | File/Location | Section/Scope | Why |
|--------|--------------|---------------|-----|
| Stage 01 | `../01-spec/output/[module]-spec.md` | Full file | what to build |
| Stage 02 | `../02-tests/output/[module]-test-plan.md`, `/test/<module>/` | Full | the target |
| Config | `/_config/architecture.md` | "Module map", "Layering" | structure |
| Config | `/_config/coding-conventions.md` | Full file | code rules |
| Code | `/src/` | exported names of lower layers only | what to reuse |

## Process

1. Implement in `/src/<Module>/` (Sonnet subagent) until the tests pass (GREEN). Tests are not edited; a needed test change goes back to 02.
2. Run the full suite, Aqua, JET, the layering test, coverage, and the `bench/` budget for this module.
3. Write implementation notes: deviations from the spec, performance numbers, known limits.
4. Run the audit, write the outputs. No checkpoint (audit only).

## Audit

| Check | Pass Condition |
|-------|---------------|
| GREEN | All tests pass, including every earlier module |
| Quality | Aqua and JET clean; layering test passes; coverage of this module at least 80% |
| Size | Functions under 50 lines, files under 400 |
| Tests untouched | `git diff` of `/test/<module>/` since 02 is empty |
| Numerics | No `@fastmath`; accumulated logs use compensated summation |

## Outputs

| Artifact | Location | Format |
|----------|----------|--------|
| Code | `/src/<Module>/` | Julia source |
| Benchmarks | `/bench/<module>.jl` | BenchmarkTools suite |
| Notes | `output/[module]-impl-notes.md` | Deviations, performance, limits (read by 04) |
