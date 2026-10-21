# Search plan: recurrence-indicators

Approved by the owner at the literature-review/01 checkpoint, 2026-10-21 (M4, decision D65).

## Question

Which recurrence-based and related phase-space indicators (recurrence-time entropy, recurrence-plot quantifiers, recurrence networks, Lagrangian descriptors, persistent homology) detect chaos and stickiness on area-preserving maps, how do they compare with the variational indicators of the P2 study on the same orbits, and how do their parameters (recurrence threshold, orbit length) change the label?

## Already in the bibliography (not searched again)

Sales2023 (doi:10.1063/5.0140613, one of the owner's identifiers), Marwan2007, Zou2016, DaquinKovacs2025, CaoChenZhang2026, TinocoErdelyi2026, JimenezLopezGarciaGarrido2026b, ChirikovShepelyansky1999 and the recurrence-time rows of `benchmarks.md`. They seed the citation chase in RI-SQ1 and RI-SQ2.

## Sub-questions, rules and queries

| SQ | Sub-question | Ties to | Include if | Queries |
|----|--------------|---------|------------|---------|
| RI-SQ1 | Recurrence-time entropy and recurrence-plot quantifiers as chaos or stickiness indicators on area-preserving maps; the owner's items doi:10.1016/j.physleta.2026.132224 (= arXiv:2606.23501), doi:10.1016/j.physleta.2009.12.057, doi:10.1007/978-3-031-91062-3_6 | C3, I16, idea 15 | defines or applies a recurrence-based quantifier to separate regular, chaotic or sticky orbits in an area-preserving map or Hamiltonian system | `recurrence time entropy stickiness`; `recurrence plot standard map chaos`; `recurrence quantification analysis Hamiltonian chaos detection`; the three identifiers |
| RI-SQ2 | Recurrence-based indicators compared with variational indicators or Lyapunov exponents on the same system | C3, I1, I2, I5, I8 | reports a recurrence-based quantifier and at least one of LE, SALI, GALI, MEGNO, FLI on the same orbits or grid | `recurrence plot Lyapunov exponent comparison map`; `recurrence quantification FLI MEGNO comparison`; forward citations of Sales2023 and Zou2016 |
| RI-SQ3 | Parameter and finite-time dependence: recurrence threshold epsilon, orbit length, sampling, embedding | C3, preregistration (P2 thresholds) | studies how a recurrence quantifier's value or label changes with threshold, length or sampling | `recurrence threshold selection recurrence plot`; `recurrence quantification orbit length dependence`; `recurrence time entropy convergence` |
| RI-SQ4 | Lagrangian descriptors and persistent-homology indicators on maps (owner's list) | I17, C3 | defines or tests a Lagrangian-descriptor or topological indicator on a map or Hamiltonian system and compares it with another indicator | `Lagrangian descriptors chaos indicator map`; `Lagrangian descriptor area-preserving map`; `persistent homology chaos indicator Hamiltonian` |
| RI-SQ5 | Resolve the owner's item arXiv:1008.2538 (SALI, 2 and 3 degrees of freedom) | I2 | the record is found and matches | arXiv ID lookup |

## Out of scope

Recurrence analysis of measured time series without a model map or Hamiltonian (that is P4); dissipative systems; applications that use a quantifier without studying it.

## Sources

As in `chaos-indicators-search-plan.md`: Crossref, arXiv (nlin.CD, math.DS, physics.data-an), OpenAlex for abstracts; ADS skipped without `ADS_API_TOKEN`; INSPIRE not used.

## Common inclusion and exclusion

- Include only if all hold: (a) peer-reviewed article, book, chapter, or arXiv preprint; (b) a DOI, arXiv ID or ADS bibcode found in this run; (c) the sub-question rule answers yes on title plus abstract.
- Exclude: not about Hamiltonian, symplectic or area-preserving dynamics; key or DOI already in `references.bib`; theses, slides, blogs.
- Every row names the sub-question that found it; every exclusion names the failed rule.

## Stop rule

Per query and source: screen the top 50 hits by relevance. A sub-question is saturated when two consecutive queries add no new included candidate. The owner's four identifiers (RI-SQ1, RI-SQ5) each end with a verified record or a stated verification failure. Hard stop at 40 included candidates; above that the owner narrows the scope.
