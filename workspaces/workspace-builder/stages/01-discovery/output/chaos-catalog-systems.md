# Catalog: Systems

Layer 3 draft. Every system ships: the map or Hamiltonian, its analytic tangent map / Hessian, its invariants, and a known-answer test case. Keys with `?` are pending verification.

## A. Area-preserving and symplectic maps (exact tangent maps, no integrator error)

| ID | System | Dim | Role | Known-answer tests | Source | Phase |
|----|--------|-----|------|--------------------|--------|-------|
| Y1 | Chirikov-Taylor standard map (lifted + wrapped coordinates) | 2 | P1/P2 workhorse | det DT = 1; K = 0 gives lambda = 0; large-K lambda_1 vs Chirikov estimate ln(K/2) | Chirikov1979, Greene1979 | M2 |
| Y2 | Henon area-preserving quadratic map | 2 | second 2D benchmark, escaping orbits | det = 1 | Henon1969 | M4 |
| Y3 | Coupled standard maps / Froeschle 4D map | 4 | GALI_k, full spectra, Arnold diffusion | pairing lambda_1 = -lambda_4, lambda_2 = -lambda_3 | Froeschle1972 | M5 |
| Y4 | N coupled standard maps (chain) | 2N | spectrum scaling, k-vector cost | pairing for all i | Skokos2004, Skokos2007 | M5 |
| Y5 | Linear symplectic maps (rotation, hyperbolic) | 2N | exact exponents for unit tests | lambda = ln|eigenvalue| exactly | none needed | M2 |

## B. Hamiltonian flows (need integrators)

| ID | System | DOF | Role | Known-answer tests | Source | Phase |
|----|--------|-----|------|--------------------|--------|-------|
| F1 | Harmonic oscillator | 1..N | exact solution, integrator order tests | energy, exact phase | none needed | M5 |
| F2 | Kepler problem | 2 | integrable, nonlinear, order tests | energy, angular momentum, Runge-Lenz | HairerLubichWanner2006 | M5 |
| F3 | Henon-Heiles | 2 | mixed phase space; escape above E = 1/6 | energy; E = 1/12, 1/8, 1/6 portraits | HenonHeiles1964 | M5 |
| F4 | Pendulum | 1 | separatrix, near-separatrix FTLE | energy | none needed | M5 |
| F5 | FPU-beta lattice / disordered Klein-Gordon chain | N | high-dimensional tangent dynamics, integrator efficiency | energy, pairing | Senyange2018, Danieli2019 | M5 (optional) |
| F6 | Circular restricted three-body problem (planar) | 2 | astro bridge, Jacobi constant; first R2 transfer target | Jacobi constant | MurrayDermott1999 | M8 |
| F7 | Forced pendulum, autonomous form H = p^2/2 - cos q + omega I + eps cos q cos phi | 2 | R2 workhorse: known separatrix, tunable stochastic layer, cheap | H conserved (autonomous form); eps = 0 integrable | Chirikov1979, LichtenbergLieberman1992 | M7 |

## C. N-body (REBOUND)

| ID | System | Role | Source | Phase |
|----|--------|------|--------|-------|
| B1 | Outer solar system (Sun + Jupiter..Neptune [+ Pluto]) | Lyapunov times, integrator study | SussmanWisdom1988, MurrayHolman1999 | M7 |
| B2 | Full solar system / inner planets | inner-planet chaos, Mercury | Laskar1989, LaskarGastineau2009 | M7 |
| B3 | Two-planet MEGNO cross-check configuration | cross-check our MEGNO vs REBOUND's | CincottaSimo2000, ReinTamayo2015 | M7 |

| B4 | Solar system with extra physics (GR, tides, J2, radiation, migration) | effect of non-Newtonian terms on Lyapunov times and on integrator choice | TamayoRein2020 (REBOUNDx) | M7 |
| B5 | Resonant / secular reduced models (Poincare variables, disturbing-function expansions, secular theory) | analytic skeleton of chaos: resonance overlap, secular resonances; compare reduced-model indicators with full N-body | HaddenTamayo2022 (celmech) | M7 |
| B6 | Ephemeris-quality test particles (asteroids, NEOs) in the DE440/441 field with variational equations | high-precision reference orbits; tangent vectors for MEGNO/LE of small bodies | Holman2023 (ASSIST) | M7 (optional) |

Initial conditions for B1-B2 come from JPL Horizons (DE440/DE441 ephemerides, Park2021) via `scripts/fetch_horizons.py`, cached with query date in the manifest. ASSIST needs the DE440/441 planet file plus the `sb441-n16` asteroid kernel (about 1 GB, downloaded once by `scripts/fetch_ephemerides.py` into `{{DATA_ROOT}}`, never committed).

Ecosystem packages (all GPL-3.0, adapters only, no re-implementation): REBOUND (integrators, MEGNO), REBOUNDx (extra forces), celmech (analytic and semi-analytic celestial mechanics), ASSIST (ephemeris-quality test particles).

## D. Time-series sources (P4)

Scalar observables x_n or p_n from Y1-Y3, optionally with additive noise and finite length; the tangent-map oracle from the same orbit is the ground truth.
