# Search plan: lyapunov-addenda

Owner-supplied list, no open search. Scoped by the coordinator on 2026-10-13; the checkpoint is the owner's own list.

## Question

Which sources back the Batch B Lyapunov module specs (decision D30: lyap-svd-stabilized, lyap-householder-qr, lyap-gram-schmidt, error control of QR methods), and which two Shevchenko-line items are wanted later for the study papers?

## Sub-questions, rules and queries

| SQ | Sub-question | Ties to | Include if (yes/no on title + abstract) | Queries |
|----|--------------|---------|------------------------------------------|---------|
| A-SQ1 | Is there a discrete SVD (or QR) algorithm for a long product of matrices usable at general d for lyap-svd-stabilized? | catalog rows L5, L9, L10, L10b in `_config/catalogs/lyapunov.md`; D30 | the item treats SVD or QR of a product or fundamental matrix, with an algorithm, for Lyapunov exponents | none: the owner supplied the items |
| A-SQ2 | What is known on the error and stability of QR-based LE computation? | L5, L8; D30 | the item gives error bounds or error control for QR methods | none: owner-supplied |
| A-SQ3 | Analytical Lyapunov time estimates and the relation between Lyapunov time and instability time (two items kept for later papers) | claim C1/C2 research line (Shevchenko line), paper sequence | the item is by Shevchenko or Cincotta et al. and states Lyapunov-time estimates for Hamiltonian resonance models | none: owner-supplied |

## Sources

Owner-supplied full texts in `data/literature/incoming/` (items 1-4, 7-9), two records without full text (items 5, 6). Identifiers are resolved at stage 03 on Crossref and arXiv (OpenAlex only as a registry record for the one item without a DOI).

## Common inclusion and exclusion

- Include exactly the owner's nine items. No searching, no snowballing. Item 9 is already in `references.bib` (DieciVanVleck2002): no new entry, full text only.
- Every row names its sub-question.

## Stop rule

Source-based: the run stops when each of the owner's nine items has a verified record or a stated verification failure.
