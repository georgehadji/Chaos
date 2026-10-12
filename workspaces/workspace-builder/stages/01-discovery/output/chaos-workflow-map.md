# Chaos: Workflow Map (ICM Stage 01 Discovery)

Output of `workspace-builder/stages/01-discovery`. Input to Stage 02 (mapping). Status: APPROVED at the Stage 01 checkpoint (2026-10-12): Julia core, GPL-3.0-or-later, initial commit pushed to georgehadji/Chaos. Local path: `E:\Documents\Vibe-Coding\Chaos Research Orchestrator`.

Companion artifacts in this folder:

| File | Holds | Becomes (Layer 3, at scaffolding) |
|------|-------|-----------------------------------|
| `chaos-research-question.md` | Question, claims C1-C3, test panel, outputs, paper sequence | `_config/research-question.md` |
| `chaos-architecture.md` | Languages, package layering, design patterns, scalability, robustness | `_config/architecture.md` + `_config/numerics-policy.md` |
| `chaos-catalog-lyapunov.md` | Every LE algorithm, status and phase | `_config/catalogs/lyapunov.md` |
| `chaos-catalog-indicators.md` | Every chaos indicator, status and phase | `_config/catalogs/indicators.md` |
| `chaos-catalog-integrators.md` | Every symplectic / reference integrator | `_config/catalogs/integrators.md` |
| `chaos-catalog-systems.md` | Maps, flows, N-body systems | `_config/catalogs/systems.md` |
| `chaos-bibliography.md`, `references.bib` | Verified references and benchmark values | `_config/bibliography/` (the .md, 255 lines, is split by its five sections at stage 03; the .bib is a data file, exempt from the 200-line rule) |

## 1. End to end

Start: a research question (in `chaos-research-question.md`). End: reproducible papers P1-P6, each backed by a tested Julia library (`Chaoslab.jl`, see `chaos-architecture.md`), versioned data with provenance, and verified citations.

## 2. Why ICM fits, and where it does not

The work is sequential and review-heavy: the math of every method must be checked by a human before code exists, and every study must be pre-registered before data exists. That is ICM's home ground. What ICM does not do (Section 5.2 of the paper) is handled by code, not by agents: parallel sweeps over 10^6 orbits, chunked storage and resumable runs live in `chaoslab.experiments` and in `scripts/`, which stages call.

## 3. Repo shape (one repo, four workspaces, one library)

```
Chaos/
  CLAUDE.md, CONTEXT.md          Layer 0, Layer 1
  _core/                         ICM conventions (vendored, MIT)
  _config/                       Layer 3 factory: research question, architecture,
                                 numerics policy, catalogs/, bibliography/
  workspaces/
    workspace-builder/           ICM builder (this plan is its output)
    literature-review/           run once per topic
    module-build/                run once per library module
    study/                       run once per study (P1, P2, ...)
    paper/                       run once per paper
  Project.toml, src/             Julia package Chaoslab (canonical code home)
  test/, bench/                  Test.jl suites, BenchmarkTools budgets
  python/nbody/                  optional N-body adapters (REBOUND ecosystem)
  scripts/                       mechanical work, no AI (Python stdlib)
  schema/                        HDF5 run-store layout shared by Julia and Python
  data/                          run outputs (git-ignored), manifests tracked
```

Workspaces are not numbered: they are independent pipelines. Dependency is one-way: literature-review promotes into `_config/bibliography/`; module-build reads `_config/` and writes `src/` + `test/`; study reads `src/` (as the activated Julia package) and `_config/`; paper reads study outputs. No workspace references a later one.

## 4. Workspaces and stages

### 4.1 literature-review (per topic, e.g. "lyapunov-algorithms")

| Stage | In | Out (`output/`) | Human gate |
|-------|----|-----------------|-----------|
| 01-scope | topic, `_config/research-question.md` | `[topic]-search-plan.md`: questions, inclusion rules, queries | approve scope |
| 02-search | search plan | `[topic]-candidates.md` with DOI / arXiv IDs (Sonnet subagents) | screen list |
| 03-verify | candidates | `[topic]-verified.bib`, `[topic]-verify-report.md` (`scripts/bibcheck.py` resolves every DOI on Crossref) | drop or fix failures |
| 04-synthesize | verified bib, papers | `[topic]-notes.md`: annotations, benchmark values with page cites, open questions | accept notes |
| 05-promote | notes, bib | diff against `_config/bibliography/` (mechanical merge by `scripts/promote.py`) | approve merge |

### 4.2 module-build (per module, e.g. "lyap-householder-qr")

| Stage | In | Out | Human gate |
|-------|----|-----|-----------|
| 01-spec | catalog row, bibliography, architecture | `[module]-spec.md`: equations, algorithm, invariants, oracles, acceptance criteria, sources. WHAT, not HOW (Pattern 10) | check the math |
| 02-tests | spec, numerics policy | test files in `test/`, `[module]-test-plan.md`; tests must fail (RED) | approve tests |
| 03-implement | spec, test plan, architecture | code in `src/<Module>/`, `[module]-impl-notes.md`; tests GREEN, Aqua, JET, layering test, perf budget | none (audit only) |
| 04-verify | code, spec | `[module]-verification.md`: oracle cross-check, convergence plot, independent review (Sonnet reviewer) | accept or send back |
| 05-integrate | verification | catalog status flip, folder CONTEXT/docs, changelog, commit (push needs approval) | approve commit |

### 4.3 study (per study, e.g. "p1-le-benchmark")

| Stage | In | Out | Human gate |
|-------|----|-----|-----------|
| 01-design | research question, catalogs | `[study]-preregistration.md`: hypotheses, panel, metrics, thresholds, budgets, analysis plan, frozen before data | approve and freeze |
| 02-panel | preregistration | `[study]-panel.toml` + phase portraits (`scripts/build_panel.py`) | inspect orbits |
| 03-run | panel, run config | `data/[study]/` (chunked store) + `[study]-manifest.json` (`scripts/run_study.py`, no AI); health checks | rerun failures |
| 04-analyze | manifest, data | figures, tables (`scripts/analyze.py`), `[study]-findings.md` | interpret |
| 05-report | findings, bibliography | `[study]-results.md`: every number traced to a manifest, every claim cited | accept |

### 4.4 paper (per paper)

01-outline, 02-draft, 03-review (adversarial Sonnet reviewer, ARS `/ars-reviewer`, `/ars-citation-check`), 04-finalize (LaTeX via ARS `/ars-format-convert`). Gates after outline and review.

## 5. Shared context (Layer 3, `_config/`)

research-question, architecture, numerics-policy (precision, lifted vs wrapped coordinates, renormalization intervals, tolerances), coding-conventions, the four catalogs, bibliography. Each fact has one home; stage CONTEXT files point to sections.

## 6. User-specific variables (setup questionnaire)

`{{AUTHOR_NAME}}`, `{{AFFILIATION}}`, `{{GITHUB_REPO}}` (default georgehadji/Chaos), `{{JULIA_VERSION}}` (current stable via juliaup), `{{PYTHON_VERSION}}` (3.12, scripts and N-body only), `{{DEFAULT_PRECISION}}` (Float64), `{{HIGH_PRECISION_BACKEND}}`, `{{COMPUTE_TARGET}}` (laptop / workstation / cluster), `{{MAX_CORES}}`, `{{DATA_ROOT}}`, `{{TARGET_JOURNALS}}`, `{{CITATION_STYLE}}`, `{{SUBAGENT_MODEL}}` (sonnet), `{{SIMPLE_TASK_MODEL}}` (haiku). Per-run values (topic, module name, study ID, K values) are asked by each entry stage, not in setup.

## 7. Optional stages

literature-review/04-synthesize can be skipped for small topics. study/02-panel is skipped for grid scans. The whole `paper/` workspace is optional. P5 solar-system work is optional (removes REBOUND prerequisites).

## 8. Tool prerequisites

| Tool | Stage | Required | Purpose |
|------|-------|----------|---------|
| Julia (juliaup, `winget install julia`); NOT installed on this machine yet | all code stages | yes | library, tests, figures; ChaosTools / TaylorIntegration as oracles |
| Python 3.12 (installed) | scripts/, N-body | yes | ICM mechanical scripts (stdlib only) |
| Git + GitHub CLI | module-build/05 | yes | commits, PRs |
| REBOUND (pip wheels since 3.28; source build with `FFP_CONTRACT_OFF` for bit-reproducible runs) | module-build, study (P5) | optional | N-body integrators, MEGNO. Native Windows: no OpenMP/MPI/OpenGL/AVX512, so WHFast512 and parallel REBOUND runs go to WSL2/Linux |
| WSL2 (Ubuntu) | study (P5) | optional | AVX512 / OpenMP REBOUND builds, cluster-like runs on Windows |
| REBOUNDx 5.1.0 (pip, sdist only: needs a C compiler, MSVC Build Tools on Windows) | module-build, study (P6) | optional | extra forces: GR, tides, J2, radiation, migration (only `gr_potential` keeps WHFast symplectic) |
| celmech 1.5.8 (pip, C extension with gcc-style flags; Windows build and REBOUND 5 compatibility unverified) | module-build, study (P6) | optional | resonant / secular models, canonical variables, FMFT |
| ASSIST 1.2.3 (pins `rebound<5`; Linux/macOS, Windows only via WSL) + DE440/441 + sb441-n16 (~1 GB; larger kernels 2-14 GB) | study (P6) | optional | ephemeris-quality test particles; separate venv or WSL2 |
| WHFast512 on native Windows | study (P6) | optional | conflicting sources (installation page: no AVX512 on native Windows; 5.1.0 changelog: Windows build supported); verify on this CPU |
| LaTeX (TeX Live / MiKTeX) | paper/04 | optional | typesetting |

## 9. Skills

Bundled: none (no installed skill covers numerical dynamics). Used from the installed ARS plugin: `/ars-lit-review`, `/ars-citation-check`, `/ars-reviewer`, `/ars-format-convert`. Global: `icm` (this method).

## 10. Agent model policy (Layer 0 rule)

The orchestrating agent runs stages. Delegated subagents: Sonnet for search, review, spec drafting, test writing; Haiku for mechanical tasks (file listing, formatting, status). Stage CONTEXT files name the model per delegated step.

## 11. Build order (runs of the workspaces)

| Milestone | Runs | Exit |
|-----------|------|------|
| M0 bootstrap | builder stages 02-05; install Julia; repo skeleton, `Project.toml`, CI, `_config/`, `schema/run-store.md` | `scripts/icm_validate.py` passes, CI green on Windows + Ubuntu |
| M1 literature | topics: lyapunov-algorithms, chaos-indicators, symplectic-integrators, standard-map-stickiness | bibliography promoted |
| M2 P1 modules | standard-map, tangent-core, precision, diagnostics, lyap-benettin, lyap-gram-schmidt, lyap-householder-qr, lyap-svd-stabilized, lyap-two-trajectory | all verified |
| M3 P1 study | p1-le-benchmark | results accepted |
| M4 P2 | modules: sali, gali, megno, fli-ofli, rli, stretching-numbers, naff-fma; study p2 | results accepted |
| M5 P3 | modules: coupled-maps, henon-heiles, integrator catalog, continuous-qr, clv-ginelli; study p3 | results accepted |
| M6 P4 | modules: embedding, lyap-wolf, lyap-rosenstein, lyap-kantz, lyap-eckmann-ruelle; study p4 | results accepted |
| M7 P5 (R2) | modules: forced-pendulum, shadow-hamiltonian (BCH), transport metrics; study p5 | falsification criteria evaluated |
| M8 P6 (conditional on P5) | modules: rebound-adapter, reboundx-effects, celmech-models, assist-adapter, horizons-fetch; study p6 | results accepted |

## 12. License (decided 2026-10-12: GPL-3.0-or-later)

`Chaos` is GPL-3.0, matching REBOUND, REBOUNDx, celmech and ASSIST. The vendored ICM files in `_core/` stay MIT with their notice kept in `_core/LICENSE-ICM` (MIT is GPL-compatible). `LICENSE` (GPL-3.0 text) and the SPDX id `GPL-3.0-or-later` (in `Project.toml` and `python/nbody/pyproject.toml`) are added at scaffolding (builder stage 03).
