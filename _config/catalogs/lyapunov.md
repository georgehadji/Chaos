# Catalog: Lyapunov-Exponent Algorithms

Layer 3, canonical. `Status` is `planned` or `verified`; only module-build/05 writes it. One row per algorithm. `Phase` = milestone that builds it (see `_config/build-order.md`). `Source` keys resolve in `/_config/bibliography/references.bib`; a key marked `?` is pending verification and must not be cited until verified.

Every variational method below consumes the same tangent stream from `Tangent.propagate` (analytic tangent map for maps, variational equations integrated with the same symplectic scheme for flows). Methods differ only in how they turn that stream into exponents.

## A. Variational, maximal exponent

| ID | Method | Gives | Vectors | Main numerical issue | Source | Phase | Status |
|----|--------|-------|---------|----------------------|--------|-------|--------|
| L1 | Benettin single-vector renormalization (finite-time mLCE lambda_1(n)) | lambda_1(n) | 1 | slow O(ln n / n) convergence for regular orbits; threshold ambiguity in weak chaos | Benettin1980a, Skokos2010 | M2 | planned |
| L2 | Two-trajectory (shadow orbit) separation with periodic renormalization | approx. lambda_1 | 0 (second orbit) | perturbation too large: nonlinear saturation; too small: round-off; wrapped vs lifted coordinates | Benettin1976, Wolf1985 | M2 | planned |

## B. Variational, full spectrum by reorthonormalization (discrete QR family)

Shared step: propagate k tangent vectors, every `tau` steps factor W = Q R, accumulate ln|R_ii|, continue with Q. Strategies differ in the factorization.

| ID | Orthonormalizer | Cost (dim d, k vectors) | Stability notes | Source | Phase | Status |
|----|-----------------|-------------------------|-----------------|--------|-------|--------|
| L3 | Classical Gram-Schmidt (CGS) | O(d k^2) | loses orthogonality when vectors nearly align; kept as historical baseline | ShimadaNagashima1979, Benettin1980b | M2 | planned |
| L4 | Modified Gram-Schmidt (MGS) | O(d k^2) | orthogonality loss ~ eps * cond(W) | Skokos2010 | M2 | planned |
| L5 | Householder QR (LAPACK geqrf) | O(d k^2) | backward stable (assumption); default reference | Geist1990, Dieci1997 | M2 | planned |
| L6 | Givens rotations QR | O(d k^2) | stable; cheap for d = 2, 4 and banded tangent maps (assumption, counts derived in the lyap-householder-qr spec) | Householder1958, DanielGragg1976 | M2 | planned |
| L4b | Repeated (re-orthogonalized) Gram-Schmidt, RGS | O(d k^2), ~2x GS | restores orthogonality for nearly parallel columns | DanielGragg1976, Geist1990 | M2 | planned |
| L5b | Householder variants H1 (treppen iteration path) / H2 (diagonal iteration path); the augmented block [W, I] is bookkeeping that returns Q explicitly | O(d k^2) | as L5; H1 deferred, H2 built in M2 (D30) | Householder1958, Geist1990 | M2 | planned |
| L7 | Reorthonormalization interval study, fixed `tau` | n/a | `tau` too long: overflow / alignment; too short: cost | Skokos2010 | M2 (parameter, not module) | verified |
| L7b | Adaptive reorthonormalization (trigger on measured orthogonality loss or local stretching) | n/a | research line R3; must not bias estimates | open | M5 | planned |
| L7c | Block / sparse QR exploiting banded tangent maps of coupled maps and lattices | O(d k b) for bandwidth b | research line R3 | open | M5 | planned |

## C. Variational, full spectrum, continuous and alternative formulations

| ID | Method | Gives | Main issue | Source | Phase | Status |
|----|--------|-------|-----------|--------|-------|--------|
| L8 | Continuous QR (ODE for Q on the orthogonal group, flows only): dQ/dt = Q S, d(ln R_ii)/dt = (Q^T J Q)_ii | full spectrum | orthogonality drift (needs projection, Cayley or a nu (Q^T Q - I) correction); reported ~40x slower than discrete QR on a 29-dim Toda chain (Geist1990, Table III p.886) | Goldhirsch1987, Dieci1997, ChristiansenRugh1997 | M5 | planned |
| L9 | Discrete SVD of the finite-time product, stabilized (product SVD / QR-accumulated SVD, no naive multiplication) | finite-time singular values sigma_i(n); spectrum as n grows | finite-time singular values are not term-by-term equal to QR diagonals | GreeneKim1987, Dieci1997, Stewart1997, OliveiraStewart2000 | M2 (d = 2); general d M3 (D32) | planned |
| L10 | Continuous SVD (Greene-Kim): ODEs for U and ln sigma_i | spectrum + singular vectors | singular when exponents are (nearly) degenerate, i.e. on tori and periodic orbits, which Hamiltonian systems always have (lambda_i = 0 pairs); costliest | GreeneKim1987, Dieci1997, DieciVanVleck2007 | catalog only (degeneracy makes it unsuitable for symplectic systems) | planned |
| L10b | Single large-time SVD of the full product | spectrum | overflow / loss of small singular values; kept only as a negative control | Geist1990 | catalog only | planned |
| L11 | Rangarajan-Habib-Ryne parametrization (no rescaling, no reorthonormalization) | full spectrum | parametrization singularities; derived ODEs per dimension | Rangarajan1998 | M5 (catalog only) | planned |
| L12 | Covariant Lyapunov vectors (Ginelli forward QR + backward upper-triangular iteration) | CLVs, local exponents, tangencies | memory for stored R factors; checkpointing | Ginelli2007, WolfeSamelson2007, KuptsovParlitz2012 | M5 | planned |
| L12b | Oseledets vectors via wedge-product SVD / dichotomy projectors (comparative study of CLV algorithms) | CLVs | cost, stability ranking per Froyland2013 | Froyland2013 | M5 (one variant as cross-check of L12) | planned |
| L13 | Finite-time LE distributions / stretching numbers / local exponents | variability of local stretching | window length dependence | Voglis1994, PikovskyPoliti2016 | M4 | planned |
| L14 | FTLE fields (grid of finite-time lambda_1) | phase-space maps, transport barriers | horizon choice | Haller2001, Shadden2005 | M4 | planned |

## D. Data-driven (time series only, no tangent map)

| ID | Method | Gives | Main issue | Source | Phase | Status |
|----|--------|-------|-----------|--------|-------|--------|
| L15 | Wolf et al. (fiducial trajectory, replacement of neighbours) | lambda_1 (and more) | parameter sensitivity, neighbour replacement | Wolf1985 | M6 | planned |
| L16 | Eckmann-Ruelle / Eckmann-Kamphorst-Ruelle-Ciliberto (local Jacobian fits + QR) | spectrum | spurious exponents, embedding dimension | EckmannRuelle1985, Eckmann1986 | M6 | planned |
| L17 | Rosenstein (mean log divergence of nearest neighbours) | lambda_1 | Theiler window, linear-region choice | Rosenstein1993 | M6 | planned |
| L18 | Kantz (averaged neighbourhood divergence) | lambda_1 | neighbourhood size | Kantz1994 | M6 | planned |

Delay embedding (Takens, false nearest neighbours) is a supporting module for L15-L18. Known caveat to test: short series of conservative maps in weak chaos give unreliable estimates (Lega2000).

## D2. Tangent-map sources (apply to every variational row)

Analytic Jacobian (default, ground truth); automatic differentiation (forward mode, ForwardDiff.jl) for maps without hand-written Jacobians, validated against the analytic one; complex-step derivative for tests. AD tangents are part of research line R3. Oracles: ChaosTools `lyapunov` / `lyapunovspectrum` (no SALI, MEGNO, FLI there), TaylorIntegration `lyap_taylorinteg` at BigFloat.

Status: verified (M2 Batch A, tangent-core: analytic Jacobians checked against ForwardDiff and complex step; stage 04 reports in `workspaces/module-build/stages/04-verify/output/`, decision D22).

## D3. Benchmark values (pointer)

Geist1990, damped-driven Toda chain: relative CPU times and the finite-time lambda_1 and lambda_28 values were checked against the paper (Tables II and III, pp.885-886, decision D30, O12) and live in _config/bibliography/benchmarks.md, rows with source key Geist1990. A reproduction target only, never a test oracle.

## E. Precision variants (apply to every row)

Float64 (default); Double64 (DoubleFloats.jl) or MultiFloats; BigFloat for short oracle runs. Precision is a run parameter recorded in the manifest, never a code fork.

Status: verified (M2 Batch A, precision: Float64, Double64, BigFloat; stage 04 reports in `workspaces/module-build/stages/04-verify/output/`, decision D22).

## F. Built-in checks every estimator reports

symplectic pairing residual eps_pair(n); eps_det(n) for maps; orthogonality loss ||Q^T Q - I||; sum of exponents (zero for Hamiltonian); agreement of SALI decay with lambda_1 - lambda_2 (when both are run).

Status: verified (M2 Batch A, diagnostics: pairing, det, orthogonality and sum residuals; stage 04 reports in `workspaces/module-build/stages/04-verify/output/`, decision D22).
