# Decision Log

Layer 3, canonical. One row per structural decision: what was chosen, why, and when to look again. Append a row whenever a checkpoint decides something that changes structure, tools, scope or rules. Never rewrite a row; supersede it with a new one.

| ID | Date | Decision | Why | Revisit when |
|----|------|----------|-----|--------------|
| D1 | 2026-10-11 | Separate repo `georgehadji/Chaos`, built with ICM | Sequential, review-heavy research fits stage contracts and human gates | never |
| D2 | 2026-10-11 | Subagents on Sonnet; mechanical subagent work on Haiku; never Opus | owner's instruction | owner says so |
| D3 | 2026-10-12 | Julia core (`Chaoslab`); Python only for N-body adapters and ICM scripts; HDF5 files between them, no in-process bridge | precision is a run parameter (Float64 to BigFloat on one kernel); forward AD works on Windows; ChaosTools and TaylorIntegration as oracles | Julia blocks a milestone for more than two weeks |
| D4 | 2026-10-12 | License GPL-3.0-or-later; `_core/` stays MIT | matches REBOUND, REBOUNDx, celmech, ASSIST | never |
| D5 | 2026-10-12 | Never cite from memory; three verification classes; UNVERIFIED items never cited | citation errors found in the seed list (see `bibliography/corrections.md`) | never |
| D6 | 2026-10-12 | This folder holds only Chaos; the LLM orchestrator moved to its own folder and is never pushed from here | different project | never |
| D7 | 2026-10-12 | Studies run through one Julia entry point, `bin/chaoslab.jl`; `scripts/` is Python stdlib for ICM mechanics only | one language for all numerics | never |
| D8 | 2026-10-12 | Catalog `Status` is binary (`planned`, `verified`), written only by module-build/05 | one gated writer; progress comes from `status <slug>` | a third state is needed by a study |
| D9 | 2026-10-12 | ICM scaffold lives at the repo root; builder outputs keep manifests only | workspaces must live where they run; no duplicate copies | never |
| D10 | 2026-10-12 | No decision LLM and no reranker in any stage | decisions are human checkpoints or preregistered numeric rules; screening by written rules is reproducible | a literature topic returns more than about 300 candidates (reranker only) |
| D11 | 2026-10-12 | P6 (solar system) runs only if P5 (shadow-Hamiltonian closure) passes its falsification criteria | novelty audit: broad physical-vs-numerical chaos question is NO-GO | P5 results accepted |
| D12 | 2026-10-12 | Release policy in `build-order.md`, section "Releases" | the library is a product other researchers can cite | first tagged release |
