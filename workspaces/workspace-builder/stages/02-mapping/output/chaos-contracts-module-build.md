# Contracts: module-build workspace

DRAFT (Stage 02). One run per library module. `[module]` is a slug from `/_config/build-order.md` (e.g. `lyap-householder-qr`); `<Module>` is the Julia module folder it lands in (from the architecture module map). Paths starting with `/` are from the repo root.

## 01-spec

### Inputs

| Source | File/Location | Section/Scope | Why |
|--------|--------------|---------------|-----|
| User | conversation | module slug | per-run input |
| Config | `/_config/build-order.md` | the milestone row naming this module | scope and prerequisites |
| Config | `/_config/catalogs/` | this module's rows only | what to build, sources |
| Config | `/_config/bibliography/` | section entries for the rows' keys | sources |
| Config | `/_config/bibliography/benchmarks.md` | rows for these methods | reference numbers |
| Reference | `references/spec-format.md` | Full file | required spec sections |
| Config | `/_config/architecture.md` | "Module map" row, "Layering" | where it sits, what it may import |
| Config | `/_config/numerics-policy.md` | Full file | precision, coordinates, renormalization rules |
| Sources | `data/literature/` | papers cited by the rows | equations and coefficient tables |

### Process

1. Confirm the module and its catalog rows; check its prerequisites are `verified` in the catalogs.
2. Draft the spec (Sonnet subagent): equations with source key and equation number, algorithm in math or pseudocode, invariants, oracle (package and settings), acceptance criteria with tolerances, coefficient tables named by paper and table number.
3. **[Checkpoint]** Present the spec.
4. Run the audit, write the output.

### Checkpoints

| After Step | Agent Presents | Human Decides |
|------------|---------------|---------------|
| 2 | Equations, invariants, oracle, acceptance criteria | The math is right, or send back |

### Audit

| Check | Pass Condition |
|-------|---------------|
| Sourced | Every equation and coefficient table names a verified key (no `?`) and an equation or table number |
| WHAT not HOW | No Julia type, function or file names (Pattern 10) |
| Measurable | Every acceptance criterion has a pass/fail test and a tolerance |
| Oracle | At least one independent oracle or analytic result is named |

### Outputs

| Artifact | Location | Format |
|----------|----------|--------|
| Module spec | `output/[module]-spec.md` | Equations, invariants, oracle, acceptance criteria, sources |

## 02-tests

### Inputs

| Source | File/Location | Section/Scope | Why |
|--------|--------------|---------------|-----|
| Previous stage | `../01-spec/output/[module]-spec.md` | Full file | what to test |
| Config | `/_config/numerics-policy.md`, `/_config/numerics-defaults.toml` | Full files | tolerances and seeds |
| Config | `/_config/coding-conventions.md` | "Tests" section | test style |
| Code | `/test/runtests.jl` | Full file | how suites are registered |
| Guide | `/_config/setup-guides/julia.md` | Full file, first run only | install and test commands |

### Process

1. Write test files in `/test/<module>/` (Sonnet subagent): invariants at seeded random points, order or convergence slopes, oracle comparisons, regression values from the bibliography benchmark section only.
2. Write the test plan: a table mapping each acceptance criterion to its tests.
3. Run the suite (`julia --project -e "using Pkg; Pkg.test()"`); the new tests must fail or error because the module does not exist yet, never because of a broken test (RED).
4. **[Checkpoint]** Present the plan and the RED run.
5. Run the audit, write the outputs.

### Checkpoints

| After Step | Agent Presents | Human Decides |
|------------|---------------|---------------|
| 3 | Test plan, RED output | Approve the tests or add cases |

### Audit

| Check | Pass Condition |
|-------|---------------|
| Coverage of spec | Every acceptance criterion maps to at least one test |
| RED | Every new test reports Fail or Error caused by missing module functionality (e.g. `UndefVarError`, `MethodError`); no parse errors, no errors inside test helpers |
| Reproducible | Every random draw uses a recorded StableRNGs seed |
| No memory values | Every regression value cites a bibliography benchmark row |

### Outputs

| Artifact | Location | Format |
|----------|----------|--------|
| Tests | `/test/<module>/` | Test.jl files |
| Test plan | `output/[module]-test-plan.md` | Criterion-to-test table, RED log |

## 03-implement

### Inputs

| Source | File/Location | Section/Scope | Why |
|--------|--------------|---------------|-----|
| Stage 01 | `../01-spec/output/[module]-spec.md` | Full file | what to build |
| Stage 02 | `../02-tests/output/[module]-test-plan.md`, `/test/<module>/` | Full | the target |
| Config | `/_config/architecture.md` | "Module map", "Layering" | structure |
| Config | `/_config/coding-conventions.md` | Full file | code rules |
| Code | `/src/` | exported names of lower layers only | what to reuse |

### Process

1. Implement in `/src/<Module>/` (Sonnet subagent) until the tests pass (GREEN). Tests are not edited; a needed test change goes back to 02.
2. Run the full suite, Aqua, JET, the layering test, coverage, and the `bench/` budget for this module.
3. Write implementation notes: deviations from the spec, performance numbers, known limits.
4. Run the audit, write the outputs. No checkpoint (audit only).

### Audit

| Check | Pass Condition |
|-------|---------------|
| GREEN | All tests pass, including every earlier module |
| Quality | Aqua and JET clean; layering test passes; coverage of this module at least 80% |
| Size | Functions under 50 lines, files under 400 |
| Tests untouched | `git diff` of `/test/<module>/` since 02 is empty |
| Numerics | No `@fastmath`; accumulated logs use compensated summation |

### Outputs

| Artifact | Location | Format |
|----------|----------|--------|
| Code | `/src/<Module>/` | Julia source |
| Benchmarks | `/bench/<module>.jl` | BenchmarkTools suite |
| Notes | `output/[module]-impl-notes.md` | Deviations, performance, limits (read by 04) |

## 04-verify

### Inputs

| Source | File/Location | Section/Scope | Why |
|--------|--------------|---------------|-----|
| Stage 01 | `../01-spec/output/[module]-spec.md` | "Acceptance criteria", "Oracle" | what to prove |
| Stage 03 | `../03-implement/output/[module]-impl-notes.md` | Full file | deviations to check |
| Code | `/src/<Module>/`, `/test/<module>/` | Full | what is verified |

### Process

1. Run the oracle cross-check at each precision the spec names; make the convergence or order plot.
2. Independent review (fresh Sonnet subagent, given the spec and the code, not the notes): math versus code, numerical traps.
3. Write the verification report: evidence for each acceptance criterion, review findings and their fixes.
4. **[Checkpoint]** Present the report.

### Checkpoints

| After Step | Agent Presents | Human Decides |
|------------|---------------|---------------|
| 3 | Evidence per criterion, plot, review findings | Accept, or send back to 01, 02 or 03 |

### Audit

| Check | Pass Condition |
|-------|---------------|
| Evidence | Every acceptance criterion has a number or a plot |
| Review closed | No open CRITICAL or HIGH review finding |

### Outputs

| Artifact | Location | Format |
|----------|----------|--------|
| Verification report | `output/[module]-verification.md` | Evidence table, review findings (read by 05) |
| Figures | `output/[module]-*.png` | Convergence and oracle plots |

## 05-integrate

### Inputs

| Source | File/Location | Section/Scope | Why |
|--------|--------------|---------------|-----|
| Previous stage | `../04-verify/output/[module]-verification.md` | Full file | accepted evidence |
| Config | `/_config/catalogs/` | this module's rows | status to flip |
| Code | `/src/<Module>/CONTEXT.md`, `/CHANGELOG.md` | Full files | docs to update |

### Process

1. Flip the catalog `Status` of this module's rows to `verified` (promotion write into `/_config/`).
2. Update the folder `CONTEXT.md` and `CHANGELOG.md`.
3. **[Checkpoint]** Present the commit.
4. Commit on a branch. Push or PR only with explicit approval.

### Checkpoints

| After Step | Agent Presents | Human Decides |
|------------|---------------|---------------|
| 2 | Diff and commit message | Approve the commit; push yes or no |

### Outputs

| Artifact | Location | Format |
|----------|----------|--------|
| Status flip | `/_config/catalogs/` | promotion write (read by study/01 and later module runs) |
| Integration log | `output/[module]-integration-log.md` | Commit SHA, files changed (final deliverable) |
