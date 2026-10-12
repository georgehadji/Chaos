# Stage 03: Run

Run the study mechanically with health checks and provenance.

## Inputs

| Source | File/Location | Section/Scope | Why |
|--------|--------------|---------------|-----|
| Stage 01 | `../01-design/output/[study]-preregistration.md` | "Methods", "Budget" | what to run |
| Stage 02 | `../02-panel/output/[study]-panel.toml` | Full file (or the grid spec in 01) | where to run |
| Schema | `/schema/run-store.md` | Full file | HDF5 layout |
| Config | `/_config/numerics-defaults.toml` | Full file | health thresholds |
| Config | `/_config/environment.md` | Full file | data root, cores |
| Guide | `/_config/setup-guides/julia.md`, `/_config/setup-guides/rebound.md` | Full files, first run only | install and run commands (rebound only for N-body) |

## Process

1. Record the preregistration SHA-256 in the manifest, then run `julia --project bin/chaoslab.jl run` (no AI): chunked, resumable, health checks on every chunk. N-body studies first run the Python adapters in `/python/nbody/`, which write HDF5 per the schema.
2. Rerun failed chunks once; log any that fail again.
3. **[Checkpoint]** Present the health summary.

## Checkpoints

| After Step | Agent Presents | Human Decides |
|------------|---------------|---------------|
| 2 | Chunks done, failed, health-check hits | Accept, rerun, or stop |

## Audit

| Check | Pass Condition |
|-------|---------------|
| Complete | Every chunk is done or logged as failed with a reason |
| Provenance | Manifest git SHA matches a commit with a clean tree; the preregistration SHA-256 recorded at run start still matches at run end |

## Outputs

| Artifact | Location | Format |
|----------|----------|--------|
| Run data | `data/[study]/` | HDF5 per `/schema/run-store.md` (git-ignored) |
| Manifest | `output/[study]-manifest.json` | Provenance, chunk list, hashes (read by 04 and 05) |
