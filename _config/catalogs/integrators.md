# Catalog: Integrators

Layer 3, canonical. `Status` is `planned` or `verified`; only module-build/05 writes it. Orders are the classical orders of each scheme; stage counts and coefficients are NOT given here and must be copied from the cited source during module-build/01-spec (Pattern 5: the paper is the canonical home of coefficients). Keys with `?` are pending verification.

Implementation rule (see architecture): every explicit splitting/composition method is a coefficient table run by ONE generic stepper. Adding a method = adding a cited table + an order test, no new loop. Where OrdinaryDiffEqSymplecticRK ships the method (SymplecticEuler, VelocityVerlet, VerletLeapfrog, PseudoVerletLeapfrog, McAte2-8, Ruth3, CandyRoz4, CalvoSanz4, Yoshida6, KahanLi6, KahanLi8, SofSpa10), its solver is an oracle for our table; it lacks Forest-Ruth, PEFRL, Blanes-Moan, SABA/SBAB, ABA864, Tao and WH. REBOUND ships leapfrog orders 2/4/6/8 and SABA2-4, SABA(10,4), SABA(8,6,4), SABA(10,6,4), SABAH variants.

## A. Explicit splitting for separable H = T(p) + V(q)

| ID | Method | Order | Notes | Source | Phase | Status |
|----|--------|-------|-------|--------|-------|--------|
| S1 | Symplectic Euler (A then B, and B then A) | 1 | baseline, composes to leapfrog | Yoshida1990, Ruth1983 (HairerLubichWanner2006 not read) | M5 | verified |
| S2 | Stormer-Verlet / leapfrog (DKD and KDK) | 2 | time-reversible; reference for all order tests | Yoshida1990, Ruth1983 (HairerLubichWanner2006 not read) | M5 | verified (GPU runs of criterion 18 not run) |
| S3 | Optimized 2-stage second-order (McLachlan) | 2 | smaller error constant than leapfrog | McLachlan1995 | M5 | verified |
| S4 | Ruth third-order | 3 | historical | Ruth1983 | M5 | verified |
| S5 | Forest-Ruth / Candy-Rozmus fourth-order | 4 | same as Yoshida triple jump of leapfrog | ForestRuth1990, Yoshida1990, Omelyan2002 (CandyRozmus1991 not read) | M5 | verified (independent read of ForestRuth1990 Eqs. (4.3), (4.4) not run; GPU runs of criterion 18 not run) |
| S6 | Yoshida triple-jump composition (4 from 2, 6 from 4, 8 from 6) | 4, 6, 8 | generic composition builder; negative steps | Yoshida1990 | M5 | verified |
| S7 | Suzuki fractal (5-stage) composition | 4+ | smaller negative steps; orders 4, 6, 8 (5, 25, 125 leapfrog pieces) | KahanLi1997, McLachlan1995 (Suzuki1990 not read; KahanLi1997 credits Suzuki's 1991 and 1992 papers, not in the bibliography) | M5 | verified |
| S8 | PEFRL (position-extended Forest-Ruth-like) | 4 | much smaller error constant than S5 | Omelyan2002 | M5 | verified (GPU runs of criterion 18 not run) |
| S9 | Calvo-Sanz-Serna fourth-order | 4 | optimized coefficients | CalvoSanzSerna1993a (Eq. (2.8)) | M5 | verified |
| S10 | Blanes-Moan SRKN (optimized RKN splitting) | 4, 6 | very efficient for separable H with T quadratic | BlanesMoan2002 | M5 | planned |
| S11 | Yoshida sixth/eighth order (solution A / D families) | 6, 8 | many stages; 6th order solution A, 8th order solutions A and D | Yoshida1990 | M5 | verified (criterion 3(b2), McLachlan1995 Table 2 against Yoshida1990 Table 1, not run) |
| S12 | Kahan-Li optimal compositions | 6, 8 | s9odr6a, s17odr8a | KahanLi1997 | M5 | verified (GPU runs of criterion 18 not run) |
| S13 | Sofroniou-Spaletta composition | 10 | high precision reference; many stages | SofroniouSpaletta2005 | M5 | verified |
| S14 | Processed (corrector) methods | var. | kernel + processor | Blanes1999 | catalog only | planned |

## B. Near-integrable H = A + eps B

| ID | Method | Order | Notes | Source | Phase | Status |
|----|--------|-------|-------|--------|-------|--------|
| P1 | SABA_n / SBAB_n | (2n, 2) in (tau, eps) | needs exact flow of A; n = 1 to 10 | LaskarRobutel2001 | M5 | verified (GPU runs of criterion 18 not run for SABA3) |
| P2 | SABAC_n (with corrector) | higher in eps | also SBABC_n | LaskarRobutel2001 | M5 | verified |
| P3 | Blanes et al. families, e.g. ABA864, ABAH864 | (8,6,4) | modified-potential variants; ABA864 and ABAH864 built, ABA104, ABA1064, ABAH844, ABAH1064 not | Blanes2013 | M5 | verified |

## C. N-body (via REBOUND, adapter only, no re-implementation)

| ID | Method | Kind | Source | Phase | Status |
|----|--------|------|--------|-------|--------|
| N1 | Wisdom-Holman (WHFast), Jacobi / democratic heliocentric, symplectic correctors | symplectic, 2nd order + correctors | WisdomHolman1991, WisdomHolmanTouma1996, ReinTamayo2015 | M7 | planned |
| N2 | High-order WH kernels (WHCKL, WHCKM, WHCKC) | symplectic | Rein2019 | M7 | planned |
| N3 | WHFast512 (AVX512) | symplectic; not available on native Windows builds (no AVX512 there), run under WSL2/Linux | Javaheri2023 | M7 | planned |
| N4 | SABA family in REBOUND | symplectic | LaskarRobutel2001, Rein2019 | M7 | planned |
| N5 | MERCURIUS, TRACE (hybrid, close encounters) | hybrid | REBOUND docs | M7 (catalog only unless needed) | planned |
| N6 | IAS15 (adaptive Gauss-Radau, 15th order) | non-symplectic, round-off limited | ReinSpiegel2015 | M7 | planned |
| N7 | Janus (bit-wise reversible) | reversible | ReinTamayo2018 | catalog only | planned |
| N8 | REBOUNDx effects inside WHFast / IAS15 (operator-split additional forces) | keeps symplectic splitting for conservative effects; dissipative ones break it | TamayoRein2020 | M7 | planned |
| N9 | ASSIST: IAS15 with DE440/441 perturbers, harmonics, GR, Marsden non-gravitational model, first-order variational equations | ephemeris-quality reference | Holman2023, Marsden1973 | M7 (optional) | planned |
| N10 | celmech: integration of reduced resonant / secular Hamiltonians (canonical variables, Lie-series transformations) | semi-analytic models | HaddenTamayo2022 | M7 | planned |

REBOUND usage rules (from rebound.hanno-rein.de quickstart, checked 2026-10-12): the adapter always sets `sim.integrator`, `sim.dt` and units explicitly (defaults are IAS15 and G = 1); the REBOUND version and build flags go into every run manifest; runs that must be bit-reproducible use a source build with `FFP_CONTRACT_OFF` set (no fused multiply-add contraction).

## D. Implicit and nonseparable

| ID | Method | Order | Notes | Source | Phase | Status |
|----|--------|-------|-------|--------|-------|--------|
| G1 | Implicit midpoint | 2 | symplectic for any H; Gauss s = 1 | Hairer2008 (HairerLubichWanner2006 not read) | M5 | verified (GPU runs of criterion 18 not run) |
| G2 | Gauss-Legendre collocation, s stages | 2s | symplectic for any H; fixed-point solve with compensated sums (D113) | Hairer2008 (HairerLubichWanner2006 not read) | M5 | verified (Newton option, closed-form coefficients and quadratic invariants deferred; GPU runs of criterion 18 not run) |
| G3 | Tao extended phase-space explicit method | 2 (composable) | nonseparable H, explicit; exponent 1/(l - 1) of Yoshida1990, not the printed 1/(l + 1) | Tao2016 | M5 | verified (tangent in the extended space not scored, open point P11) |
| G4 | Pihajoki extended phase space | 2 | leapfrog for nonseparable H | Pihajoki2015 | catalog only | planned |

## E. Non-symplectic references (truth generators for tests)

| ID | Method | Use | Source | Phase | Status |
|----|--------|-----|--------|-------|--------|
| R1 | Taylor series method, high order, arbitrary precision, jet transport for variational equations | reference orbits; reference tangents by finite differences (jet transport deferred) | JorbaZou2005 | M5 (TaylorIntegration.jl, BigFloat; package extension) | verified |
| R2 | High-order explicit RK (OrdinaryDiffEq DP8 / Vern9) | quick non-symplectic comparison | HairerNorsettWanner1993 | M5 | planned |
| R3 | Bulirsch-Stoer | catalog | HairerNorsettWanner1993 | catalog only | planned |

## F. Variational equations

Tangent map method: integrate the variational equations with the SAME splitting as the orbit, so tangent dynamics stays symplectic and consistent (SkokosGerlach2010). Efficiency comparisons: GerlachSkokos2011, Gerlach2012, Senyange2018, Danieli2019.

## F2. Modified (shadow) Hamiltonian

For every splitting in A and B, the leading terms of H_h come from the Baker-Campbell-Hausdorff expansion of the composed flows, i.e. nested Poisson brackets of T and V (Reich1999, HairerLubichWanner2006). Computed symbolically (Symbolics.jl) and code-generated to the numeric kernels (`build_function`); validated by checking that H_h is conserved to higher order than H along numerical orbits (slope test vs h). Needed by research line R2 (P5).

## G. Checks every integrator must pass (tests, module-build/02)

Measured order = nominal order (log-log slope of energy or state error vs step, with tolerance); time-reversibility (forward then backward) for symmetric methods; symplecticity of the one-step map Jacobian (J^T Omega J = Omega) at random points; bounded energy error over long runs (no secular drift beyond round-off); Brouwer-law round-off growth with compensated summation (Hairer2008).
