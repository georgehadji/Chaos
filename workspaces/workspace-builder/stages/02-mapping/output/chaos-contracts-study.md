# Contracts: study workspace

DRAFT (Stage 02). One run per study. `[study]` is a slug such as `p1-le-benchmark`. A study reads library code only as the activated Julia package (`/src/`), never module-build `output/` folders (Pattern 14). Run data lives in `data/[study]/` (git-ignored); the manifest is tracked.

## 01-design

### Inputs

| Source | File/Location | Section/Scope | Why |
|--------|--------------|---------------|-----|
| User | conversation | study slug, paper it serves | per-run input |
| Config | `/_config/research-question.md` | "Three claims to test", "Test panel", "Outputs retained", the paper's row | what the study must decide |
| Config | `/_config/catalogs/` | rows in scope, `Status` column | only verified methods may run |
| Config | `/_config/value-types.md` | Full file | agree the value the study delivers (Pattern 13) |
| Config | `/_config/numerics-defaults.toml` | Full file | defaults the study may override |
| Config | `/_config/research-ideas.md` | the idea row this study serves | scope |
| Config | `/_config/environment.md` | Full file | compute limits for the budget |
| Reference | `references/preregistration-format.md` | Full file | required sections |

### Process

1. Agree the study scope and its value types with the human.
2. Draft the preregistration (Sonnet subagent): hypotheses, panel, methods, metrics, decision thresholds, compute budget, analysis plan, stopping rules.
3. **[Checkpoint]** Present the preregistration.
4. Run the audit; freeze: record the git SHA and the file's SHA-256 in its header and commit it. Later changes are amendments, logged before any data exists.

### Checkpoints

| After Step | Agent Presents | Human Decides |
|------------|---------------|---------------|
| 1 | Scope and value types | Confirm the study is worth running |
| 2 | Full preregistration | Approve and freeze, or revise |

### Audit

| Check | Pass Condition |
|-------|---------------|
| Verified methods | Every method in scope has catalog `Status` = `verified` |
| Decidable | Every hypothesis has a metric and a numeric decision rule |
| Single home | Thresholds appear only in this file |
| Budget | Estimated CPU hours fit the limits in `/_config/environment.md` |

### Outputs

| Artifact | Location | Format |
|----------|----------|--------|
| Preregistration | `output/[study]-preregistration.md` | Frozen; header holds git SHA and SHA-256 |

## 02-panel (skip for grid scans)

### Inputs

| Source | File/Location | Section/Scope | Why |
|--------|--------------|---------------|-----|
| Previous stage | `../01-design/output/[study]-preregistration.md` | "Panel" section | orbit classes to build |

### Process

1. Run `julia --project bin/chaoslab.jl panel` with the preregistration (no AI): picks initial conditions per class, writes the panel file and phase portraits.
2. **[Checkpoint]** Present the portraits.
3. If a class is wrong, the human amends `[study]-preregistration.md` in 01's output (logged amendment, allowed only before 03-run), then this stage reruns.

### Checkpoints

| After Step | Agent Presents | Human Decides |
|------------|---------------|---------------|
| 1 | Phase portraits with the chosen orbits marked | Orbits match their classes, or amend |

### Outputs

| Artifact | Location | Format |
|----------|----------|--------|
| Panel | `output/[study]-panel.toml` | Initial conditions, class labels, seeds (read by 03) |
| Portraits | `output/[study]-portraits.png` | Phase portraits (read by the human at the checkpoint) |

## 03-run

### Inputs

| Source | File/Location | Section/Scope | Why |
|--------|--------------|---------------|-----|
| Stage 01 | `../01-design/output/[study]-preregistration.md` | "Methods", "Budget" | what to run |
| Stage 02 | `../02-panel/output/[study]-panel.toml` | Full file (or the grid spec in 01) | where to run |
| Schema | `/schema/run-store.md` | Full file | HDF5 layout |
| Config | `/_config/numerics-defaults.toml` | Full file | health thresholds |
| Config | `/_config/environment.md` | Full file | data root, cores |
| Guide | `/_config/setup-guides/julia.md`, `/_config/setup-guides/rebound.md` | Full files, first run only | install and run commands (rebound only for N-body) |

### Process

1. Record the preregistration SHA-256 in the manifest, then run `julia --project bin/chaoslab.jl run` (no AI): chunked, resumable, health checks on every chunk. N-body studies first run the Python adapters in `/python/nbody/`, which write HDF5 per the schema.
2. Rerun failed chunks once; log any that fail again.
3. **[Checkpoint]** Present the health summary.

### Checkpoints

| After Step | Agent Presents | Human Decides |
|------------|---------------|---------------|
| 2 | Chunks done, failed, health-check hits | Accept, rerun, or stop |

### Audit

| Check | Pass Condition |
|-------|---------------|
| Complete | Every chunk is done or logged as failed with a reason |
| Provenance | Manifest git SHA matches a commit with a clean tree; the preregistration SHA-256 recorded at run start still matches at run end |

### Outputs

| Artifact | Location | Format |
|----------|----------|--------|
| Run data | `data/[study]/` | HDF5 per `/schema/run-store.md` (git-ignored) |
| Manifest | `output/[study]-manifest.json` | Provenance, chunk list, hashes (read by 04 and 05) |

## 04-analyze

### Inputs

| Source | File/Location | Section/Scope | Why |
|--------|--------------|---------------|-----|
| Stage 01 | `../01-design/output/[study]-preregistration.md` | "Analysis plan", "Hypotheses" | what to compute |
| Stage 03 | `../03-run/output/[study]-manifest.json` | Full file | which data |

### Process

1. Run `julia --project bin/chaoslab.jl analyze` (no AI): figures and tables exactly as preregistered.
2. Draft findings (Sonnet subagent): one decision per hypothesis (confirmatory), then exploratory results labelled as such.
3. **[Checkpoint]** Present figures and findings.
4. Run the audit, write the outputs.

### Checkpoints

| After Step | Agent Presents | Human Decides |
|------------|---------------|---------------|
| 2 | Figures, tables, draft findings | Interpretation is right, or revise |

### Audit

| Check | Pass Condition |
|-------|---------------|
| Plan followed | Every confirmatory result uses the frozen analysis plan; deviations listed |
| Regenerable | Every figure and table rebuilds from the manifest with one command |

### Outputs

| Artifact | Location | Format |
|----------|----------|--------|
| Figures and tables | `output/[study]-fig-*.pdf`, `output/[study]-tables.md` | Read by 05 and paper |
| Findings | `output/[study]-findings.md` | Confirmatory decisions, exploratory notes |

## 05-report

### Inputs

| Source | File/Location | Section/Scope | Why |
|--------|--------------|---------------|-----|
| Stage 04 | `../04-analyze/output/[study]-findings.md`, `[study]-tables.md` | Full files | results |
| Stage 03 | `../03-run/output/[study]-manifest.json` | Full file | traceability |
| Config | `/_config/bibliography/` | entries cited | citations |

### Process

1. Write the results document (Sonnet subagent): every number traced to a table or the manifest, every claim cited, limitations stated.
2. **[Checkpoint]** Present the results.
3. Run the audit, write the output, commit.

### Checkpoints

| After Step | Agent Presents | Human Decides |
|------------|---------------|---------------|
| 1 | Results document | Accept, or revise |

### Audit

| Check | Pass Condition |
|-------|---------------|
| Traced | Every number points to a table row or manifest field |
| Cited | Every citation key is verified (no `?`, no UNVERIFIED) |
| Style | No em dashes |

### Outputs

| Artifact | Location | Format |
|----------|----------|--------|
| Results | `output/[study]-results.md` | Final deliverable of the study (read by paper) |
