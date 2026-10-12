# Stage 01: Spec

Write the math contract for one library module: equations, invariants, oracle, acceptance criteria.

## Inputs

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

## Process

1. Confirm the module and its catalog rows; check its prerequisites are `verified` in the catalogs.
2. Draft the spec (Sonnet subagent): equations with source key and equation number, algorithm in math or pseudocode, invariants, oracle (package and settings), acceptance criteria with tolerances, coefficient tables named by paper and table number.
3. **[Checkpoint]** Present the spec.
4. Run the audit, write the output.

## Checkpoints

| After Step | Agent Presents | Human Decides |
|------------|---------------|---------------|
| 2 | Equations, invariants, oracle, acceptance criteria | The math is right, or send back |

## Audit

| Check | Pass Condition |
|-------|---------------|
| Sourced | Every equation and coefficient table names a verified key (no `?`) and an equation or table number |
| WHAT not HOW | No Julia type, function or file names (Pattern 10) |
| Measurable | Every acceptance criterion has a pass/fail test and a tolerance |
| Oracle | At least one independent oracle or analytic result is named |

## Outputs

| Artifact | Location | Format |
|----------|----------|--------|
| Module spec | `output/[module]-spec.md` | Equations, invariants, oracle, acceptance criteria, sources |
