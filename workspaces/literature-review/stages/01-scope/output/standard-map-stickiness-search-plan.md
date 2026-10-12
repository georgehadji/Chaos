# Search plan: standard-map-stickiness

Approved by the owner at the literature-review/01 checkpoint, 2026-10-12.

## Question

What is known quantitatively about stickiness, trapping-time statistics, accelerator modes and their effect on finite-time measures in the standard map, so P1/P2 sticky orbits can be chosen and labelled from the literature?

## Sub-questions, rules and queries

| SQ | Sub-question | Ties to | Include if | Queries |
|----|--------------|---------|------------|---------|
| M-SQ1 | Trapping / recurrence-time statistics and power-law tails in area-preserving maps | C2, test panel "sticky" class | reports trapping or recurrence time distributions in an area-preserving map | `stickiness area-preserving map recurrence time distribution`; `power-law trapping times standard map`; `Poincare recurrences mixed phase space` |
| M-SQ2 | Accelerator modes of the standard map (lifted coordinates, anomalous transport) | Y1, lifted-coordinate rule | studies accelerator modes or anomalous diffusion in the standard map | `accelerator modes standard map`; `anomalous diffusion standard map ballistic`; `standard map diffusion coefficient K` |
| M-SQ3 | Sticky regions and finite-time LE or indicator signatures | C2, C3, panel orbit classes | links stickiness to finite-time LE distributions or indicator values | `stickiness finite-time Lyapunov exponents`; `sticky orbits Lyapunov exponent distribution`; `cantori stickiness finite-time` |
| M-SQ4 | Island-around-island hierarchy and cantori near the golden-mean torus | test panel near-critical K | describes island hierarchies, cantori or critical-torus structure in the standard map | `island hierarchy standard map stickiness`; `cantori standard map golden mean`; `critical invariant circle standard map` |
| M-SQ5 | Selecting representative orbits (periodic, KAM, chaotic, sticky) for the standard map | test panel design | gives a method or catalogue for locating periodic orbits or classifying orbits in the standard map | `periodic orbits standard map numerical location`; `standard map orbit classification phase space` |

## Out of scope

Quantum standard (kicked rotor) map, dissipative variants.

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
