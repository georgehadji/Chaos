# Catalogs

One row per method or system, with sources, build milestone (`Phase`, see `../build-order.md`) and `Status` (`planned` or `verified`; a qualifier in parentheses names a deferred part, for example `verified (d = 2)`; a row with an open criterion stays `planned`). Load only the rows a stage names.

| File | Rows |
|------|------|
| `lyapunov.md` | L1-L18: variational QR, SVD, CLV, FTLE, time-series LE methods; tangent sources, precision, built-in checks |
| `indicators.md` | I1-I18: SALI, GALI, FLI, OFLI, MEGNO, RLI, stretching numbers, FMA/NAFF, WBA, REM, 0-1 test, RQA |
| `integrators.md` | S1-S14, P1-P3, N1-N10, G1-G4, R1-R3: splitting, near-integrable, N-body, implicit, reference methods; required checks |
| `systems.md` | Y1-Y5 maps, F1-F7 flows, B1-B6 N-body, time-series sources |
