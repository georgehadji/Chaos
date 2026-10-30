# Search plan: indicators-addenda

Owner-supplied list, no open search. Scoped by the coordinator on 2026-10-30; the checkpoint is the owner's approval of this plan.

## Question

Which of the owner's newly supplied chaos-indicator sources can enter the bibliography, under which keys and verification class, and is the existing entry Simo2014 citable?

## Items

| # | Item (as held) | In the bibliography now | What the run must settle |
|---|----------------|-------------------------|--------------------------|
| 1 | Ramasubramanian, K. & Sriram, M. S., "Global geometric indicator of chaos and Lyapunov exponents in Hamiltonian systems", Phys. Rev. E 64, 046207 (2001), doi:10.1103/PhysRevE.64.046207 (owner-pasted record and abstract; no full text held) | no (only RamasubramanianSriram2000, a different paper) | new key RamasubramanianSriram2001; DOI record match |
| 2 | Froeschle, C. & Lega, E., "The fine structure of Hamiltonian systems revealed using the Fast Lyapunov Indicator", in Steves, Maciejewski & Hendry (eds.), Chaotic Worlds: From Order to Disorder in Gravitational N-Body Dynamical Systems, NATO Science Series II vol. 227, Springer 2006, pp. 131 to 165 (book PDF held; contents page read) | no | new key FroeschleLega2006; chapter DOI or book ISBN record |
| 3 | Lopez Pina, A., "A control approach to complexity: from molecular vibrational dynamics to competition games on networks", doctoral thesis, Universidad Politecnica de Madrid, February 2019 (PDF held; title and abstract pages read) | no | new key LopezPina2019; repository record (oa.upm.es) or registry record |
| 4 | Jimenez-Lopez, J. & Garcia-Garrido, V. J., "Learning the chaotic and regular nature of trajectories in Hamiltonian systems with Lagrangian descriptors" (PDF held is the SSRN preprint 4973458, marked not peer reviewed) | yes, as LopezGarrido2024 (arXiv 2407.18831) | the stored author names are inverted ("Lopez, Javier Jimenez"), so the key breaks the key format; check the arXiv record, whether a journal version exists, and propose the corrected key JimenezLopezGarciaGarrido2024 with a line in corrections.md (the key appears in 5 tracked files) |
| 5 | Simo, C., "Measuring the total amount of chaos in some Hamiltonian systems" (PDF held: AIMS manuscript template, preprint) | yes, as Simo2014 (note "verified via secondary sources only; no DOI"), and also listed in unverified.md (final venue not confirmed) | the two files disagree. Search for the published version (Crossref by title and author). Found and matched: class verified, entry updated, line removed from unverified.md. Not found: class stays verified-secondary only if an independent reference list or registry record confirms it, else the entry leaves references.bib and stays in unverified.md; the catalog row I12c that cites it is flagged |

Already in the bibliography, no action: Darriba2012 (the owner-pasted full text; the Froeschle-Lega 2006 chapter is in its reference list), Chaoticus2025 (arXiv 2507.00622; PDF held).

## Sub-questions, rules and queries

| SQ | Sub-question | Ties to | Include if (yes/no on title + abstract) | Queries |
|----|--------------|---------|------------------------------------------|---------|
| B-SQ1 | Is there a global (phase-space or energy-surface) measure of the amount of chaos to set beside the per-orbit indicators of C3? | claim C3; catalog row I12c (Simo2014); research ideas "Prior art to verify" | the item defines or measures a fraction or global average of chaos over a set of initial conditions or an energy surface | none: owner-supplied (items 1, 5) |
| B-SQ2 | Which primary source states the FLI behaviour and resolving power used for catalog rows I5 and I6? | claim C3; catalog rows I5, I6 | the item defines the FLI or OFLI and shows its behaviour on regular, resonant or chaotic orbits | none: owner-supplied (item 2) |
| B-SQ3 | Machine-learning and Lagrangian-descriptor classification of orbits | claim C3; catalog rows I17, I18 | the item classifies orbits as regular or chaotic with Lagrangian descriptors or a trained classifier | none: owner-supplied (item 4) |
| B-SQ4 | Dependence of the share of chaotic orbits on a physical parameter (energy, forcing frequency) in a Hamiltonian system | claim C3, research line R2 (physical vs numerical transport) | the item measures the share of chaotic trajectories as a function of a system parameter | none: owner-supplied (item 3; weak tie, the thesis is mostly control theory: the owner may drop it here) |

## Sources

Owner-supplied files outside the repository (the owner's Downloads folder) and one pasted record (item 1). Nothing is copied under `data/` in this run (owner instruction 2026-10-30); full texts stay where they are and are cited by their record, not by a local path. Identifiers are resolved at stage 03 on Crossref (items 1, 2, 5), arXiv (item 4), and the UPM repository or a thesis registry (item 3).

## Common inclusion and exclusion

- Include exactly items 1 to 5. No searching beyond the title and author queries of item 5 (published version) and item 4 (journal version); no snowballing.
- Every new key follows `_config/bibliography/conventions.md` (ASCII family names, up to three authors, year).
- A key change (item 4) and a class change (item 5) are recorded in `corrections.md` at stage 05.

## Stop rule

Source-based: the run stops when each of the five items has a verified record, a verified-secondary record with its confirming source named, or a stated verification failure.
