# Stage 04: Analyze

Produce the preregistered figures and tables, then the findings.

## Inputs

| Source | File/Location | Section/Scope | Why |
|--------|--------------|---------------|-----|
| Stage 01 | `../01-design/output/[study]-preregistration.md` | "Analysis plan", "Hypotheses" | what to compute |
| Stage 03 | `../03-run/output/[study]-manifest.json` | Full file | which data |

## Process

1. Run `julia --project bin/chaoslab.jl analyze` (no AI): figures and tables exactly as preregistered.
2. Draft findings (Sonnet subagent): one decision per hypothesis (confirmatory), then exploratory results labelled as such.
3. **[Checkpoint]** Present figures and findings.
4. Run the audit, write the outputs.

## Checkpoints

| After Step | Agent Presents | Human Decides |
|------------|---------------|---------------|
| 2 | Figures, tables, draft findings | Interpretation is right, or revise |

## Audit

| Check | Pass Condition |
|-------|---------------|
| Plan followed | Every confirmatory result uses the frozen analysis plan; deviations listed |
| Regenerable | Every figure and table rebuilds from the manifest with one command |

## Outputs

| Artifact | Location | Format |
|----------|----------|--------|
| Figures and tables | `output/[study]-fig-*.pdf`, `output/[study]-tables.md` | Read by 05 and paper |
| Findings | `output/[study]-findings.md` | Confirmatory decisions, exploratory notes |
