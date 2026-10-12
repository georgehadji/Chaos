# Stage 01: Scope

Turn a topic into an approved search plan.

## Inputs

| Source | File/Location | Section/Scope | Why |
|--------|--------------|---------------|-----|
| User | conversation | topic slug and the question it serves | per-run input |
| Config | `/_config/research-question.md` | "Primary question", "Three claims to test", "Paper sequence" | keeps the search tied to the claims |
| Config | `/_config/research-ideas.md` | "Prior art to verify before claiming any gap" | known gaps to check |
| Config | `/_config/bibliography/CONTEXT.md` | Full file | what is already covered, per section |

## Process

1. Ask for the topic slug and the question it answers.
2. Read the covered sections of the bibliography for this topic.
3. Draft the search plan: sub-questions, inclusion and exclusion rules, sources (Crossref, arXiv, ADS, INSPIRE), query strings, stop rule.
4. **[Checkpoint]** Present the plan.
5. Run the audit, revise, write the output.

## Checkpoints

| After Step | Agent Presents | Human Decides |
|------------|---------------|---------------|
| 3 | Sub-questions, rules, queries, stop rule | Approve or narrow the scope |

## Audit

| Check | Pass Condition |
|-------|---------------|
| Tied to claims | Every sub-question names a claim (C1-C3), a research line, or a catalog row ID |
| Testable rules | Each inclusion rule can be applied to a title and abstract with a yes/no answer |
| Stop rule | A numeric or source-based stop condition is stated |

## Outputs

| Artifact | Location | Format |
|----------|----------|--------|
| Search plan | `output/[topic]-search-plan.md` | Sub-questions, rules, queries, stop rule |
