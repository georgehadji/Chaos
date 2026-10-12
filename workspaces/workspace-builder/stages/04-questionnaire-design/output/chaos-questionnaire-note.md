# Chaos: Questionnaire (ICM Stage 04)

The questionnaire lives where it runs: `/setup/questionnaire.md` (same reason as the Stage 03 manifest). Date: 2026-10-12.

## Placeholder coverage

| Placeholder | File | Question |
|-------------|------|----------|
| `{{AUTHOR_NAME}}` | `_config/publication.md` | Q1 |
| `{{AFFILIATION}}` | `_config/publication.md` | Q2 |
| `{{ORCID}}` | `_config/publication.md` | Q3 |
| `{{TARGET_JOURNALS}}` | `_config/publication.md` | Q4 |
| `{{CITATION_STYLE}}` | `_config/publication.md` | Q5 |
| `{{COMPUTE_TARGET}}` | `_config/environment.md` | Q6 |
| `{{MAX_CORES}}` | `_config/environment.md` | derived from Q6 |
| `{{DATA_ROOT}}` | `_config/environment.md` | Q7 |
| `{{HIGH_PRECISION_BACKEND}}` | `_config/environment.md` | Q8 |

Conditional: Q9 removes the paper workspace. No `{{?SECTION}}` blocks were needed.

## Per-run values (not in the questionnaire)

Collected by the entry stage of each run: topic slug (literature-review/01), module slug (module-build/01), study slug and paper (study/01), paper slug and target journal (paper/01).

## Not asked (decided, see `_config/decision-log.md`)

Repo, license, subagent models, language and Julia install route, default precision Float64.

## Audit

| Check | Result |
|-------|--------|
| Placeholder coverage | PASS: 9 of 9 placeholders covered (8 asked, 1 derived) |
| Per-run separation | PASS: no slugs or run parameters asked |
| Flat structure | PASS: Q1-Q9, one list |
| Defaults present | PASS: every question has a default |
