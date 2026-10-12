# Chaos: Research Idea Backlog

Layer 3, canonical. Source: user-supplied scoping reviews (2026-10-12). Every cited claim here is UNVERIFIED until it passes literature-review/03-verify; keys with `?` must not be cited yet. Each idea maps to a research line in `chaos-research-question.md` (R1 LE vs indicators, R2 shadow-Hamiltonian transport, R3 scalable computation) and to a paper slot.

## Ranking rule

Rank = novelty x feasibility with the planned library x venue fit. Ideas that reuse the same `chaoslab` modules rank higher (one engine, many papers). Re-rank at each study/01-design checkpoint.

## Top 5 (adopted into the plan)

| Rank | Idea | Line | Paper | Needs (modules) |
|------|------|------|-------|-----------------|
| 1 | Unified, pre-registered benchmark of LE methods and chaos indicators on symplectic maps (standard, Henon, Froeschle, coupled maps, Henon-like cubic map), with runtime, detection latency, misclassification vs long-run reference; tests the Maffione2011a "FLI+RLI global, MEGNO+SALI local" package | R1 | P1, P2 | P1/P2 modules |
| 2 | Open benchmark suite + library (`chaoslab`) for LE / CLV / indicators on symplectic maps, N >= 4 | R3 | JOSS + P3 | all core modules |
| 3 | Block / sparse QR Benettin for banded tangent maps of coupled maps and lattices (N = 10..100) | R3 | P3 | L7c, Y4, F5 |
| 4 | Shadow-Hamiltonian predictive closure for excess transport (conditional GO) | R2 | P5 | BCH module, F7, FMA |
| 5 | Error analysis + mixed precision for QR-based LE on symplectic maps (where round-off, not truncation, limits lambda_i(n); pairing residual as error witness) | R1/R3 | P1 appendix or own paper | precision backends, diagnostics |

## Backlog (not scheduled; promote by editing this table)

| # | Idea | Cluster | Line | Notes |
|---|------|---------|------|-------|
| 6 | Hybrid QR/SVD for partial spectra (leading modes, avoid continuous-SVD degeneracy) | methods | R3 | after P3 |
| 7 | Adaptive reorthonormalization triggered by orthogonality loss / local stretching | methods | R3 | L7b; must not bias estimates |
| 8 | AD tangent maps for complex symplectic maps (generating functions, integrator maps) | methods | R3 | D2 in Lyapunov catalog |
| 9 | Efficient LE on quasi-periodic tori via torus embedding, adapted to maps | methods | R1 | Fiedler2024 |
| 10 | GPU LE / CLV / indicators for N-dim maps (Chaoticus-style, maps not flows) | hardware | R3 | Chaoticus2025 as prior art and oracle |
| 11 | Mixed precision FP16/32/64 for indicators on GPU | hardware | R3 | after 10 |
| 12 | Parallel block-QR CLVs for Hamiltonian lattices (distributed) | hardware | R3 | cluster only |
| 13 | LDI vs SALI vs GALI in high-dimensional maps | indicators | R1 | Sales2026 |
| 14 | Stability exponents and CLVs vs plain LE spectra (Cvitanovic-style) | indicators | R1 | ChaosBook lecture notes, grey lit |
| 15 | Entropy-based indicators and recurrence plots on symplectic maps | indicators | R1 | I16 |
| 16 | Performance of indicators on symplectic Henon-like polynomial maps in higher dimension (CERN line) | indicators | R1 | Bazzani2023, Montanari2025 |
| 17 | Indicator "packages" per task (global map vs single orbit), formalized | indicators | R1 | part of idea 1 |
| 18 | Surrogate ML models predicting LE / indicators from short orbits | ML | R1 | leakage risk; needs idea-1 data |
| 19 | Active-learning / Bayesian-optimization exploration of parameter space near regular/chaotic boundaries | automation | R3 | uses ICM study workspace |
| 20 | Agent-driven automated exploration pipeline proposing conjectures | automation | R3 | speculative |
| 21 | AI-assisted vs classical exploration: efficiency and discovery rate | automation | R3 | speculative |
| 22 | Discrete QR for delay / renewal equations | methods | none | out of scope (non-Hamiltonian) |

## Prior art to verify before claiming any gap

Geist1990 (QR/SVD comparison); Lu1997? (efficient discrete QR, Physica D pii S0167278996002163); simple discrete QR 2006? (Appl. Math. Comput. pii S0096300305005254); Dieci-Van Vleck toolbox (LESLIS/LESLIL, Hermite interpolation, Sacker-Sell spectrum); Strathclyde QR error analysis?; QR LE for DAEs?; Maffione2011a, Maffione2011b, Darriba2012 and arXiv:1108.2196; Skokos & Manos (arXiv:1412.7401); CERN performance analysis on symplectic Henon-like cubic map (cds.cern.ch/record/2856415); Chaoticus GPU package (arXiv:2507.00622); TNO stability preprint (arXiv:2607.13629); LDI paper (Chaos Solitons Fractals 205, 2026?); Nordita/ChaosBook Lyapunov lecture notes; DL classification of motion types (Sci. Rep. 2022, s41598-022-05696-9). Existing software to compare against: DynamicalSystems.jl / ChaosTools (lyapunov, lyapunovspectrum, gali, 0-1 test, expansion entropy), REBOUND (MEGNO, frequency analysis), gala, lyapynov, nolds, pyNAFF / nafflib.
