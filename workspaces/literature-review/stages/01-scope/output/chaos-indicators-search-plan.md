# Search plan: chaos-indicators

Approved by the owner at the literature-review/01 checkpoint, 2026-10-12.

## Question

How do variational and non-variational indicators compare in classification, detection time and threshold choice on symplectic maps, and is the pending LGALI source resolvable?

## Sub-questions, rules and queries

| SQ | Sub-question | Ties to | Include if | Queries |
|----|--------------|---------|------------|---------|
| I-SQ1 | Comparative studies of indicators (SALI, GALI, MEGNO, FLI/OFLI, RLI, FMA, 0-1, LDs) published 2016-2026 | C3, I1-I12 | compares two or more chaos indicators on the same system | `comparison chaos indicators symplectic map`; `SALI GALI MEGNO FLI comparison`; `chaos indicators performance comparison Hamiltonian 2020` |
| I-SQ2 | Thresholds, detection time and classification stability | C3, preregistration | studies how thresholds or integration time affect regular/chaotic classification | `chaos indicator threshold classification time`; `detection time chaotic orbit SALI threshold`; `MEGNO threshold regular chaotic classification` |
| I-SQ3 | Indicators on sticky orbits and near separatrices | C2, C3 | reports indicator behaviour on sticky or weakly chaotic orbits | `sticky orbits chaos indicators`; `weak chaos detection indicator stickiness`; `SALI sticky orbit` |
| I-SQ4 | Resolve catalog I4: origin of LGALI / GALI for local dynamics | I4 (`Skokos2008?`) | defines or introduces LGALI, or GALI for local dynamics | `linear dependence GALI local dynamics`; `LGALI indicator`; `generalized alignment index local dynamics Skokos` |
| I-SQ5 | Indicators in accelerator physics (dynamic aperture), including the CERN CDS 2856415 item | prior art, Montanari2025 cluster | studies chaos indicators for accelerator lattices or Henon-like maps | INSPIRE: `chaos indicators dynamic aperture`; `reversibility error method accelerator`; `Henon map chaos indicators accelerator` |
| I-SQ6 | Machine-learning classification of orbits (prior art; DL Sci. Rep. 2022 from research-ideas) | prior art for P2 novelty | classifies regular vs chaotic orbits with a learned model | `deep learning classification regular chaotic orbits`; DOI 10.1038/s41598-022-05696-9 |

## Out of scope

Indicator applications to specific astronomical systems without method content.

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
