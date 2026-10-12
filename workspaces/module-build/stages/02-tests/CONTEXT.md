# Stage 02: Tests

Write failing tests (RED) that encode the spec.

## Inputs

| Source | File/Location | Section/Scope | Why |
|--------|--------------|---------------|-----|
| Previous stage | `../01-spec/output/[module]-spec.md` | Full file | what to test |
| Config | `/_config/numerics-policy.md`, `/_config/numerics-defaults.toml` | Full files | tolerances and seeds |
| Config | `/_config/coding-conventions.md` | "Tests" section | test style |
| Code | `/test/runtests.jl` | Full file | how suites are registered |
| Guide | `/_config/setup-guides/julia.md` | Full file, first run only | install and test commands |

## Process

1. Write test files in `/test/<module>/` (Sonnet subagent): invariants at seeded random points, order or convergence slopes, oracle comparisons, regression values from the bibliography benchmark section only.
2. Write the test plan: a table mapping each acceptance criterion to its tests.
3. Run the suite (`julia --project -e "using Pkg; Pkg.test()"`); the new tests must fail or error because the module does not exist yet, never because of a broken test (RED).
4. **[Checkpoint]** Present the plan and the RED run.
5. Run the audit, write the outputs.

## Checkpoints

| After Step | Agent Presents | Human Decides |
|------------|---------------|---------------|
| 3 | Test plan, RED output | Approve the tests or add cases |

## Audit

| Check | Pass Condition |
|-------|---------------|
| Coverage of spec | Every acceptance criterion maps to at least one test |
| RED | Every new test reports Fail or Error caused by missing module functionality (e.g. `UndefVarError`, `MethodError`); no parse errors, no errors inside test helpers |
| Reproducible | Every random draw uses a recorded StableRNGs seed |
| No memory values | Every regression value cites a bibliography benchmark row |

## Outputs

| Artifact | Location | Format |
|----------|----------|--------|
| Tests | `/test/<module>/` | Test.jl files |
| Test plan | `output/[module]-test-plan.md` | Criterion-to-test table, RED log |
