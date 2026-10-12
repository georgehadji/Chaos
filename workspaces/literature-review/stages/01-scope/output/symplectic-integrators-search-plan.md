# Search plan: symplectic-integrators

Approved by the owner at the literature-review/01 checkpoint, 2026-10-12.

## Question

Which open source questions (pending catalog keys, tangent-map integration, round-off behaviour, modified-Hamiltonian computation) must be resolved before the integrator catalog (M5) and the shadow-Hamiltonian module (M7)?

## Sub-questions, rules and queries

| SQ | Sub-question | Ties to | Include if | Queries |
|----|--------------|---------|------------|---------|
| S-SQ1 | Resolve catalog S9: which Calvo and Sanz-Serna 1993 paper holds the fourth-order coefficients | S9 | is a 1993 Calvo / Sanz-Serna paper on symplectic Runge-Kutta-Nystrom or optimized fourth-order methods | Crossref author query `Calvo Sanz-Serna 1993`; `fourth-order symplectic Runge-Kutta-Nystrom methods` |
| S-SQ2 | Resolve catalog N7: source of Janus (bit-wise reversible integrator) | N7 (`ReinTamayo2018?`) | introduces or describes the Janus integrator | `Janus bit-wise reversible integrator`; `reversible integrator integer arithmetic N-body` |
| S-SQ3 | Integrating variational equations with symplectic schemes (tangent map method) since 2012 | Tangent, P3 | integrates deviation vectors or tangent maps with symplectic or splitting methods | `tangent map method symplectic integrator variational equations`; `symplectic integration variational equations chaos indicators` |
| S-SQ4 | Round-off error growth in long symplectic integrations (Brouwer's law, compensated summation) | precision arm, R2 controls | measures or reduces round-off error growth in symplectic integration | `round-off error symplectic integrator Brouwer law`; `compensated summation symplectic integration long-term`; `floating point error long-term integration` |
| S-SQ5 | Computing modified Hamiltonians of splitting methods (BCH series, symbolic or algorithmic) | Shadow module, R2, M7 | gives an algorithm or explicit terms for the modified Hamiltonian or BCH series of a splitting | `modified Hamiltonian splitting method explicit terms`; `Baker-Campbell-Hausdorff series algorithm Lie algebra`; `backward error analysis splitting methods modified equation computation` |
| S-SQ6 | Processed / corrector (modified-kernel) splittings | R2 third arm | introduces or analyses processed or kernel-modified symplectic methods | `processed symplectic integrator corrector`; `symplectic correctors kernel` |

## Out of scope

N-body-specific integrator engineering (covered; M8), non-symplectic adaptive schemes.

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
