# Catalog: Systems

Layer 3, canonical. `Status` is `planned` or `verified`; only module-build/05 writes it. Every system ships: the map or Hamiltonian, its analytic tangent map / Hessian, its invariants, and a known-answer test case. Keys with `?` are pending verification.

## A. Area-preserving and symplectic maps (exact tangent maps, no integrator error)

| ID | System | Dim | Role | Known-answer tests | Source | Phase | Status |
|----|--------|-----|------|--------------------|--------|-------|--------|
| Y1 | Chirikov-Taylor standard map (lifted + wrapped coordinates) | 2 | P1/P2 workhorse | det DT = 1; K = 0 gives lambda = 0; large-K lambda_1 vs Chirikov estimate ln(K/2) | Chirikov1979, Greene1979 | M2 | verified |
| Y2 | Henon area-preserving quadratic map | 2 | second 2D benchmark, escaping orbits | det = 1 | Henon1969 | M4 | planned |
| Y3 | Coupled standard maps / Froeschle 4D map (variant of Skokos2010 Eq. (16), parameters nu, kappa, mu) | 4 | GALI_k, full spectra, Arnold diffusion | pairing lambda_1 = -lambda_4, lambda_2 = -lambda_3; mu = 0 is Y1 with K = -nu under the shear (x + p, p); orbits C and D of Skokos2003a (SALI of C about 0.28, lambda_1 of D about 8e-3) | Skokos2010, Skokos2003a (Froeschle1972 origin, not read) | M5 | verified |
| Y4 | N coupled standard maps (chain, 2 <= N <= 100, open ends or ring for N >= 3, nearest-neighbour sine coupling eps = 2 beta; general N derived from the printed N = 2 map, D103) | 2N | spectrum scaling, k-vector cost | pairing for all i; eps = 0 is N copies of Y1 bitwise; exact origin spectrum from the chain Laplacian | Maffione2011a (N = 2, Sect. 6) | M5 | verified |
| Y5 | Linear symplectic maps (rotation, hyperbolic) | 2N | exact exponents for unit tests | lambda = ln|eigenvalue| exactly | none needed | M2 | verified |

## B. Hamiltonian flows (need integrators)

| ID | System | DOF | Role | Known-answer tests | Source | Phase | Status |
|----|--------|-----|------|--------------------|--------|-------|--------|
| F1 | Harmonic oscillator | 1..N | exact solution, integrator order tests | energy, exact phase | none needed | M5 | verified (integrator-catalog: `HarmonicOscillator`) |
| F2 | Kepler problem | 2 | integrable, nonlinear, order tests | energy, angular momentum, Runge-Lenz | HairerLubichWanner2006 | M5 | verified (integrator-catalog: `Kepler`, Hamiltonian and first integrals; the exact solution by Kepler's equation deferred, HairerLubichWanner2006 not read) |
| F3 | Henon-Heiles | 2 | mixed phase space; escape above E = 1/6 | energy; critical points and escape energy; E = 1/12, 1/8, 1/6 portraits as data (chaotic fractions, SALI at T_f = 2000) | AntonopoulosManosSkokos2010 Sect. 3.1, Hairer2008 Sect. 2 (HenonHeiles1964 not read) | M5 | verified (henon-heiles: `HenonHeilesFlow`, energy-surface starts in Systems, `section_map` and `escape_orbit` in Integrators; GPU part of criterion 16 not run; portrait at 200 x 200, the paper's 500 x 500 not run) |
| F4 | Pendulum | 1 | separatrix, near-separatrix FTLE | energy | none needed | M5 | verified (integrator-catalog: `Pendulum`, also the near-integrable form p^2/2 - eps cos q; near-separatrix FTLE not built) |
| F5 | FPU-beta lattice / disordered Klein-Gordon chain | N | high-dimensional tangent dynamics, integrator efficiency | energy, pairing | Senyange2018, Danieli2019 | M5 (optional) | planned |
| F6 | Circular restricted three-body problem (planar) | 2 | astro bridge, Jacobi constant; first R2 transfer target | Jacobi constant | MurrayDermott1999 | M8 | planned |
| F7 | Forced pendulum, autonomous form H = p^2/2 - cos q + omega I + eps cos q cos phi | 2 | R2 workhorse: known separatrix, tunable stochastic layer, cheap | H conserved (autonomous form); eps = 0 integrable | Chirikov1979, LichtenbergLieberman1992 | M7 | planned |

## C. N-body (REBOUND)

| ID | System | Role | Source | Phase | Status |
|----|--------|------|--------|-------|--------|
| B1 | Outer solar system (Sun + Jupiter..Neptune [+ Pluto]) | Lyapunov times, integrator study | SussmanWisdom1988, MurrayHolman1999 | M7 | planned |
| B2 | Full solar system / inner planets | inner-planet chaos, Mercury | Laskar1989, LaskarGastineau2009 | M7 | planned |
| B3 | Two-planet MEGNO cross-check configuration | cross-check our MEGNO vs REBOUND's | CincottaSimo2000, ReinTamayo2015 | M7 | planned |

| B4 | Solar system with extra physics (GR, tides, J2, radiation, migration) | effect of non-Newtonian terms on Lyapunov times and on integrator choice | TamayoRein2020 (REBOUNDx) | M7 |
| B5 | Resonant / secular reduced models (Poincare variables, disturbing-function expansions, secular theory) | analytic skeleton of chaos: resonance overlap, secular resonances; compare reduced-model indicators with full N-body | HaddenTamayo2022 (celmech) | M7 |
| B6 | Ephemeris-quality test particles (asteroids, NEOs) in the DE440/441 field with variational equations | high-precision reference orbits; tangent vectors for MEGNO/LE of small bodies | Holman2023 (ASSIST) | M7 (optional) |

Initial conditions for B1-B2 come from JPL Horizons (DE440/DE441 ephemerides, Park2021) via `scripts/fetch_horizons.py`, cached with query date in the manifest. ASSIST needs the DE440/441 planet file plus the `sb441-n16` asteroid kernel (about 1 GB, downloaded once by `scripts/fetch_ephemerides.py` into the data root (`_config/environment.md`), never committed).

Ecosystem packages (all GPL-3.0, adapters only, no re-implementation): REBOUND (integrators, MEGNO), REBOUNDx (extra forces), celmech (analytic and semi-analytic celestial mechanics), ASSIST (ephemeris-quality test particles).

## D. Time-series sources (P4)

Scalar observables x_n or p_n from Y1-Y3, optionally with additive noise and finite length; the tangent-map oracle from the same orbit is the ground truth.
