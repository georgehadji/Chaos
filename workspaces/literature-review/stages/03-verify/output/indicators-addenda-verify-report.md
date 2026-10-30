# Verification report: indicators-addenda

Stage 03. Input: included rows of `../02-search/output/indicators-addenda-candidates.md`, one identifier per paper (DOI for all five). Tool: `scripts/bibcheck.py` (Crossref, DataCite and arXiv records), then `scripts/verify_report.py`; both were run with `--root` on a scratch copy of the repository files they read, so nothing was written under `data/` (owner instruction). Every BibTeX field comes from the fetched record except the fixes listed below.

## Summary

- Fetch failures: 0 of 5.
- In `indicators-addenda-verified.bib`: 3 new entries, all class verified (DOI records): RamasubramanianSriram2001, FroeschleLega2006, LopezPina2019.
- 2 records update existing entries and are not in the .bib, because `verify_report.py` drops a record whose title equals an existing entry. Stage 05 applies them by hand (section "Updates to existing entries").
- Keys follow `conventions.md`; the three new keys do not clash with `references.bib` (checked by script). The proposed key JimenezLopezGarciaGarrido2025 does not clash either (existing: JimenezLopezGarciaGarrido2026a, 2026b).

## Fixes and notes

- RamasubramanianSriram2001: the Crossref record gives the second author as "Sriram, M."; set to "Sriram, M. S." as on the owner-pasted APS record and as in the authors' 2000 record (doi:10.1103/physreve.62.4850). Crossref gives no page or article number; the article number 046207 is the DOI suffix.
- FroeschleLega2006: the Crossref record is a book chapter (type book-chapter, container "Chaotic Worlds: From Order to Disorder in Gravitational N-Body Dynamical Systems", pages 131-165); the script wrote it as `@article` with the book title as journal. Set to `@incollection` with `booktitle`. Editors (Steves, Maciejewski, Hendry), series (NATO Science Series II, vol. 227) and publisher (Springer) are read from the held book PDF (cover, series page, copyright page: e-book ISBN 978-1-4020-4706-0, which the chapter DOI extends); its contents page gives the chapter at p. 131, the next chapter at p. 167.
- LopezPina2019: the record is a dissertation (type dissertation, DOI prefix of the UPM thesis repository); the script wrote `@misc`. Set to `@phdthesis` with school "Universidad Politecnica de Madrid" as printed on the title page of the held PDF (Tesis doctoral, February 2019, director J. C. Losada Gonzalez). The space before the colon in the record title is removed. Relevance is weak (B-SQ4): the chaos content is a laser-perturbed HCN molecule, where the share of chaotic trajectories rises when the laser frequency is resonant (abstract, p. iv).

## Updates to existing entries (stage 05, by hand)

1. Simo2014 (keep the key). The existing entry is `@misc`, grey literature, with the note "Verified via secondary sources only; no DOI", and the same item is in `unverified.md` ("final published venue not confirmed"). Found: the published version, Crossref record doi:10.3934/dcds.2014.34.5135. Class becomes verified. Stage 05: replace the entry with the record below and remove the Simo 2014 line from `unverified.md`; record the change in `corrections.md`; the annotation in `maps.md` and the catalog row I12c keep the key. The held PDF is the AIMS manuscript, not the journal typesetting: page numbers read in it are manuscript pages.

```
@article{Simo2014,
  author = {Simó, Carles},
  title = {Measuring the total amount of chaos in some Hamiltonian systems},
  year = {2014},
  journal = {Discrete and Continuous Dynamical Systems},
  volume = {34},
  number = {12},
  pages = {5135-5164},
  doi = {10.3934/dcds.2014.34.5135}
}
```

2. LopezGarrido2024 becomes JimenezLopezGarciaGarrido2025. The existing entry (arXiv 2407.18831) stores the authors as "Lopez, Javier Jimenez" and "Garrido, Victor Jose Garcia" (family names cut), so its key breaks the key format. Found: the journal version, Crossref record doi:10.1016/j.chaos.2024.115876 (Chaos, Solitons and Fractals 191, 115876, 2025), same title, family names Jimenez-Lopez and Garcia-Garrido. Class verified. Stage 05: replace the entry with the record below (the arXiv ID kept as eprint), update the annotation in `indicators.md`, and record the key change in `corrections.md`. The key LopezGarrido2024 appears elsewhere only in earlier literature-review outputs (02-search, 03-verify, 05-promote logs of chaos-indicators), which stay as they are; no study, spec or paper cites it. The held PDF is the SSRN preprint (doi:10.2139/ssrn.4973458, marked not peer reviewed).

```
@article{JimenezLopezGarciaGarrido2025,
  author = {Jiménez-López, Javier and García-Garrido, V. J.},
  title = {Learning the chaotic and regular nature of trajectories in Hamiltonian systems with Lagrangian descriptors},
  year = {2025},
  journal = {Chaos, Solitons \& Fractals},
  volume = {191},
  pages = {115876},
  doi = {10.1016/j.chaos.2024.115876},
  eprint = {2407.18831}
}
```

## Full-text checks

| Item | File (owner's Downloads, not filed) | Against the record |
|---|---|---|
| 1 RamasubramanianSriram2001 | none (pasted record and abstract) | title, authors, journal, volume, date (20 September 2001) and DOI equal the record |
| 2 FroeschleLega2006 | The_fine_structure_of_Hamiltonian_systems_revealed.pdf (the whole book, 344 pp.) | contents page: chapter title and authors, p. 131; copyright page: Springer 2006, e-book ISBN 978-1-4020-4706-0 |
| 3 LopezPina2019 | A CONTROL APPROACH TO COMPLEXITY.pdf | title page: title, author, UPM, Tesis doctoral, Febrero 2019 |
| 4 JimenezLopezGarciaGarrido2025 | Learning the Chaotic and Regular Nature Of Trajectories ... .pdf | title and authors match; SSRN preprint 4973458, not the journal typesetting |
| 5 Simo2014 | MEASURING THE TOTAL AMOUNT OF CHAOS IN SOME.pdf | title, author, abstract; AIMS manuscript template ("Manuscript submitted to AIMS' Journals"), not the journal typesetting |

## Dropped at verification

- JimenezLopezGarciaGarrido2025 and Simo2014 were dropped by the script as same-title records of existing entries; they are carried as updates above, not lost.

## UNVERIFIED

None: every included candidate resolved to a fetched record.
