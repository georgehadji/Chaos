# Build Order

Layer 3, canonical. Milestones and the module, topic and study slugs each run uses. Catalog `Phase` columns point here.

## Milestones

| Milestone | Runs | Exit |
|-----------|------|------|
| M0 bootstrap | builder stages 02-05; install Julia; repo skeleton, `Project.toml`, CI, `_config/`, `schema/run-store.md` | `scripts/icm_validate.py` passes, CI green on Windows + Ubuntu |
| M1 literature | topics: lyapunov-algorithms, chaos-indicators, symplectic-integrators, standard-map-stickiness | bibliography promoted |
| M2 P1 modules | standard-map, tangent-core, precision, diagnostics, lyap-benettin, lyap-gram-schmidt, lyap-householder-qr, lyap-svd-stabilized, lyap-two-trajectory | all verified |
| M3 P1 study | p1-le-benchmark | results accepted |
| M4 P2 | modules: sali, gali, megno, fli-ofli, rli, stretching-numbers, naff-fma; study p2 | results accepted |
| M5 P3 | modules: coupled-maps, henon-heiles, integrator catalog, continuous-qr, clv-ginelli; study p3 | results accepted |
| M6 P4 | modules: embedding, lyap-wolf, lyap-rosenstein, lyap-kantz, lyap-eckmann-ruelle; study p4 | results accepted |
| M7 P5 (R2) | modules: forced-pendulum, shadow-hamiltonian (BCH), transport metrics; study p5 | falsification criteria evaluated |
| M8 P6 (conditional on P5) | modules: rebound-adapter, reboundx-effects, celmech-models, assist-adapter, horizons-fetch; study p6 | results accepted |

Every milestone exit: review `risks.md`; record structural decisions in `decision-log.md`.

## Releases

- Version tags `vMAJOR.MINOR.PATCH` (Semantic Versioning) on the Julia package `Chaoslab`; `0.x` until the P1-P3 API settles.
- First tag `v0.1.0` at M3 exit (P1 modules verified, P1 study accepted); a minor version per later milestone.
- GitHub-Zenodo integration (set up in M0) gives every tag a citable DOI; each paper cites the tag and the study manifests it used.
- User documentation with Documenter.jl from M2; a JOSS paper on the library after M5.
- Each release: tests and CI green, `CHANGELOG.md` entry, catalog `Status` matches the code.
