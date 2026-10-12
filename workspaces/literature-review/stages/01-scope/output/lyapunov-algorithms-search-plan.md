# Search plan: lyapunov-algorithms

Approved by the owner at the literature-review/01 checkpoint, 2026-10-12.

## Question

What does the literature already establish about agreement, finite-time behaviour, numerical stability and precision sensitivity of LE algorithms on area-preserving maps, and what benchmark values can serve as P1 regression targets?

## Sub-questions, rules and queries

| SQ | Sub-question | Ties to | Include if (yes/no on title + abstract) | Queries |
|----|--------------|---------|------------------------------------------|---------|
| L-SQ1 | Published lambda_1 values or curves lambda(K) for the standard map with stated method and length | C1, Y1, benchmarks | reports numerical LE values for the standard map (or Chirikov-Taylor map) | `"standard map" Lyapunov exponent numerical K`; `Chirikov standard map maximal Lyapunov exponent dependence on K`; `Lyapunov exponent standard map accelerator modes` |
| L-SQ2 | Distributions and convergence of finite-time LEs in mixed Hamiltonian phase space | C2, L1, I1 | studies finite-time LE distributions, convergence rates or time scales in Hamiltonian or area-preserving systems | `finite-time Lyapunov exponent distribution Hamiltonian`; `finite-time Lyapunov spectrum area-preserving map stickiness`; `convergence Lyapunov exponent estimate mixed phase space` |
| L-SQ3 | Effect of floating-point precision and round-off on computed LEs; multiple or extended precision | P1 precision arm, idea 5 | measures LE dependence on arithmetic precision, round-off or shadowing | `Lyapunov exponent round-off error floating point precision`; `multiple precision Lyapunov exponents chaotic map`; `double-double quadruple precision chaotic orbit` |
| L-SQ4 | Reorthonormalization interval, orthogonality loss and QR stability in practice | L3-L7 | analyses or measures Gram-Schmidt / Householder / Givens behaviour or the renormalization interval in LE computation | `reorthonormalization interval Lyapunov exponents`; `orthogonality loss Gram-Schmidt Lyapunov spectrum`; `Householder Givens QR Lyapunov exponents comparison` |
| L-SQ5 | Two-trajectory (shadow orbit) LE estimation: perturbation size and renormalization | L2 | studies the two-trajectory or separation method and its parameters | `two-trajectory method Lyapunov exponent perturbation size`; `separation of nearby orbits Lyapunov renormalization initial distance` |
| L-SQ6 | Prior art for a modern LE-algorithm benchmark (2010-2026) and software with LE routines | research-ideas prior art, R3 | compares two or more LE algorithms, or presents software computing LEs | `comparison Lyapunov exponent algorithms benchmark`; `software Lyapunov spectrum Julia Python package`; `DynamicalSystems.jl`; arXiv 2607.13629 (TNO preprint, from research-ideas) |

## Out of scope

Time-series LE estimators (topic for M6), CLVs beyond what P1 needs (M5).

## Sources

| Source | Access | Used for |
|--------|--------|----------|
| Crossref | REST API, no key (`CROSSREF_MAILTO` polite pool) | journal articles, DOIs |
| arXiv | export API, no key | preprints (nlin.CD, math.DS, astro-ph.EP, physics.comp-ph, physics.acc-ph) |
| OpenAlex | REST API, no key | abstracts for screening (Crossref often has none); every hit still needs a DOI or arXiv ID |
| INSPIRE | REST API, no key | accelerator-physics work only (chaos-indicators SQ5) |
| ADS | needs an API token; none set on this machine | skipped unless the owner sets `ADS_API_TOKEN` |

## Common inclusion and exclusion

- Include only if all hold: (a) peer-reviewed article, book, chapter, or arXiv preprint; (b) a DOI, arXiv ID or ADS bibcode found in this run; (c) the sub-question rule below answers yes on title plus abstract.
- Exclude: not about Hamiltonian, symplectic or area-preserving dynamics (unless the sub-question allows dissipative systems); applications that use an indicator without studying it; key or DOI already in `references.bib`; theses, slides, blogs.
- Every row names the sub-question that found it; every exclusion names the failed rule.

## Stop rule

Per query and source: screen the top 50 hits by relevance. A sub-question is saturated when two consecutive queries add no new included candidate. Hard stop at 60 included candidates per topic; above that the human narrows the scope.
