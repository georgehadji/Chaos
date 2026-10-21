# Search plan: weighted-birkhoff

Approved by the owner at the literature-review/01 checkpoint, 2026-10-21 (M4, decision D65).

## Question

How do weighted Birkhoff averages (WBA) classify orbits of area-preserving maps, how fast do they converge on tori, on chaotic orbits and on sticky orbits, how do they compare with Lyapunov exponents and frequency analysis on the same orbits, and does prior work already cover the multi-observable, multiscale convergence profile with coarse-graining of idea 25?

## Already in the bibliography (not searched again)

Das2017 (WBA, quantitative quasiperiodicity) and SanderMeiss2025 (exponents from WBA, compared with MEGNO); Laskar1990, LaskarFroeschleCelletti1992 and Laskar1993 (frequency analysis, catalog I12). They seed the citation chase in WB-SQ1 to WB-SQ3.

## Sub-questions, rules and queries

| SQ | Sub-question | Ties to | Include if | Queries |
|----|--------------|---------|------------|---------|
| WB-SQ1 | WBA as a regular versus chaotic classifier on maps; the owner's items arXiv:2001.00086 (Sander and Meiss, WBA classification of the standard map) and arXiv:1512.07286 (check whether it is Das2017) | C3, I13, idea 25 | uses a weighted Birkhoff average or its convergence to label orbits of a map or Hamiltonian system | `weighted Birkhoff average chaos detection`; `weighted Birkhoff average standard map`; `quantitative quasiperiodicity weighted Birkhoff`; the two arXiv IDs |
| WB-SQ2 | Convergence theory and practice: super-convergence on tori, behaviour on chaotic and sticky orbits, weight functions, digits thresholds, dependence on the observable | I13, idea 25 hypotheses H1 and H2 | states or measures the convergence rate of a WBA, or how it depends on the weight function or the observable | `weighted Birkhoff average convergence rate`; `Birkhoff average weight function quasiperiodic super-convergence`; `rotation number weighted Birkhoff`; forward citations of Das2017 |
| WB-SQ3 | WBA compared with Lyapunov exponents, frequency analysis (NAFF), the 0-1 test or variational indicators on the same orbits | C3, I1, I12, I15 | reports WBA and at least one other indicator on the same orbits or grid | `weighted Birkhoff average Lyapunov exponent comparison`; `weighted Birkhoff frequency analysis comparison`; forward citations of SanderMeiss2025 |
| WB-SQ4 | Frequency analysis background and the effect of filtering or coarse-graining on the classification; the owner's items Laskar arXiv:astro-ph/9906176 and the 2003 Chaos paper on filtering chaotic signals (author to identify) | I12, idea 25 hypothesis H3 | studies frequency analysis of maps, or how low-pass filtering or coarse-graining of an orbit changes a regular or chaotic label | `introduction to frequency map analysis Laskar`; `filtering chaotic signals` (Crossref, journal Chaos, 2003); `coarse graining chaos detection orbit filtering`; the arXiv ID |

## Out of scope

Ergodic-theory results on Birkhoff averages without a numerical or classification use; continuous-time flows unless the method is transferred to maps.

## Sources

As in `chaos-indicators-search-plan.md`: Crossref, arXiv (math.DS, nlin.CD, astro-ph.EP), OpenAlex for abstracts; ADS skipped without `ADS_API_TOKEN`; INSPIRE not used.

## Common inclusion and exclusion

- Include only if all hold: (a) peer-reviewed article, book, chapter, or arXiv preprint; (b) a DOI, arXiv ID or ADS bibcode found in this run; (c) the sub-question rule answers yes on title plus abstract.
- Exclude: not about Hamiltonian, symplectic or area-preserving dynamics, except WB-SQ4 frequency analysis of flows when it is applied to maps; key or DOI already in `references.bib`; theses, slides, blogs.
- Every row names the sub-question that found it; every exclusion names the failed rule.

## Stop rule

Per query and source: screen the top 50 hits by relevance. A sub-question is saturated when two consecutive queries add no new included candidate. The owner's four items (WB-SQ1, WB-SQ4) each end with a verified record or a stated verification failure. Hard stop at 40 included candidates; above that the owner narrows the scope.
