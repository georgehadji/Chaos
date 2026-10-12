# Catalog: Lyapunov-Exponent Algorithms

Layer 3 draft. One row per algorithm. `Phase` = milestone that builds it (see `_config/build-order.md`). `Source` keys resolve in `references.bib`; a key marked `?` is pending verification and must not be cited until verified.

Every variational method below consumes the same tangent stream from `chaoslab.tangent.propagate` (analytic tangent map for maps, variational equations integrated with the same symplectic scheme for flows). Methods differ only in how they turn that stream into exponents.

## A. Variational, maximal exponent

| ID | Method | Gives | Vectors | Main numerical issue | Source | Phase |
|----|--------|-------|---------|----------------------|--------|-------|
| L1 | Benettin single-vector renormalization (finite-time mLCE lambda_1(n)) | lambda_1(n) | 1 | slow O(ln n / n) convergence for regular orbits; threshold ambiguity in weak chaos | Benettin1980a, Skokos2010 | M2 |
| L2 | Two-trajectory (shadow orbit) separation with periodic renormalization | approx. lambda_1 | 0 (second orbit) | perturbation too large: nonlinear saturation; too small: round-off; wrapped vs lifted coordinates | Benettin1976, Wolf1985 | M2 |

## B. Variational, full spectrum by reorthonormalization (discrete QR family)

Shared step: propagate k tangent vectors, every `tau` steps factor W = Q R, accumulate ln|R_ii|, continue with Q. Strategies differ in the factorization.

| ID | Orthonormalizer | Cost (dim d, k vectors) | Stability notes | Source | Phase |
|----|-----------------|-------------------------|-----------------|--------|-------|
| L3 | Classical Gram-Schmidt (CGS) | O(d k^2) | loses orthogonality when vectors nearly align; kept as historical baseline | ShimadaNagashima1979, Benettin1980b | M2 |
| L4 | Modified Gram-Schmidt (MGS) | O(d k^2) | orthogonality loss ~ eps * cond(W) | Skokos2010 | M2 |
| L5 | Householder QR (LAPACK geqrf) | O(d k^2) | backward stable; default reference | Geist1990, Dieci1997 | M2 |
| L6 | Givens rotations QR | O(d k^2) | stable; cheap for d = 2, 4 and banded tangent maps | Geist1990 | M2 |
| L4b | Repeated (re-orthogonalized) Gram-Schmidt, RGS | O(d k^2), ~2x GS | restores orthogonality for nearly parallel columns | DanielGragg1976, Geist1990 | M2 |
| L5b | Householder variants H1 (reflectors on W) / H2 (augmented block [W, I]) | O(d k^2) | as L5; H2 avoids forming Q explicitly | Householder1958, Geist1990 | M2 |
| L7 | Reorthonormalization interval study, fixed `tau` | n/a | `tau` too long: overflow / alignment; too short: cost | Skokos2010 | M2 (parameter, not module) |
| L7b | Adaptive reorthonormalization (trigger on measured orthogonality loss or local stretching) | n/a | research line R3; must not bias estimates | open | M5 |
| L7c | Block / sparse QR exploiting banded tangent maps of coupled maps and lattices | O(d k b) for bandwidth b | research line R3 | open | M5 |

## C. Variational, full spectrum, continuous and alternative formulations

| ID | Method | Gives | Main issue | Source | Phase |
|----|--------|-------|-----------|--------|-------|
| L8 | Continuous QR (ODE for Q on the orthogonal group, flows only): dQ/dt = Q S, d(ln R_ii)/dt = (Q^T J Q)_ii | full spectrum | orthogonality drift (needs projection, Cayley or a nu (Q^T Q - I) correction); reported ~40x slower than discrete QR on a 31-dim Toda chain (Geist1990, to verify) | Goldhirsch1987, Dieci1997, ChristiansenRugh1997 | M5 |
| L9 | Discrete SVD of the finite-time product, stabilized (product SVD / QR-accumulated SVD, no naive multiplication) | finite-time singular values sigma_i(n); spectrum as n grows | finite-time singular values are not term-by-term equal to QR diagonals | GreeneKim1987, Dieci1997 | M2 |
| L10 | Continuous SVD (Greene-Kim): ODEs for U and ln sigma_i | spectrum + singular vectors | singular when exponents are (nearly) degenerate, i.e. on tori and periodic orbits, which Hamiltonian systems always have (lambda_i = 0 pairs); costliest | GreeneKim1987, Dieci1997, DieciVanVleck2007 | catalog only (degeneracy makes it unsuitable for symplectic systems) |
| L10b | Single large-time SVD of the full product | spectrum | overflow / loss of small singular values; kept only as a negative control | Geist1990 | catalog only |
| L11 | Rangarajan-Habib-Ryne parametrization (no rescaling, no reorthonormalization) | full spectrum | parametrization singularities; derived ODEs per dimension | Rangarajan1998 | M5 (catalog only) |
| L12 | Covariant Lyapunov vectors (Ginelli forward QR + backward upper-triangular iteration) | CLVs, local exponents, tangencies | memory for stored R factors; checkpointing | Ginelli2007, WolfeSamelson2007, KuptsovParlitz2012 | M5 |
| L12b | Oseledets vectors via wedge-product SVD / dichotomy projectors (comparative study of CLV algorithms) | CLVs | cost, stability ranking per Froyland2013 | Froyland2013 | M5 (one variant as cross-check of L12) |
| L13 | Finite-time LE distributions / stretching numbers / local exponents | variability of local stretching | window length dependence | Voglis1994, PikovskyPoliti2016 | M4 |
| L14 | FTLE fields (grid of finite-time lambda_1) | phase-space maps, transport barriers | horizon choice | Haller2001, Shadden2005 | M4 |

## D. Data-driven (time series only, no tangent map)

| ID | Method | Gives | Main issue | Source | Phase |
|----|--------|-------|-----------|--------|-------|
| L15 | Wolf et al. (fiducial trajectory, replacement of neighbours) | lambda_1 (and more) | parameter sensitivity, neighbour replacement | Wolf1985 | M6 |
| L16 | Eckmann-Ruelle / Eckmann-Kamphorst-Ruelle-Ciliberto (local Jacobian fits + QR) | spectrum | spurious exponents, embedding dimension | EckmannRuelle1985, Eckmann1986 | M6 |
| L17 | Rosenstein (mean log divergence of nearest neighbours) | lambda_1 | Theiler window, linear-region choice | Rosenstein1993 | M6 |
| L18 | Kantz (averaged neighbourhood divergence) | lambda_1 | neighbourhood size | Kantz1994 | M6 |

Delay embedding (Takens, false nearest neighbours) is a supporting module for L15-L18. Known caveat to test: short series of conservative maps in weak chaos give unreliable estimates (Lega2000).

## D2. Tangent-map sources (apply to every variational row)

Analytic Jacobian (default, ground truth); automatic differentiation (forward mode, ForwardDiff.jl) for maps without hand-written Jacobians, validated against the analytic one; complex-step derivative for tests. AD tangents are part of research line R3. Oracles: ChaosTools `lyapunov` / `lyapunovspectrum` (no SALI, MEGNO, FLI there), TaylorIntegration `lyap_taylorinteg` at BigFloat.

## D3. Candidate benchmark values (user-supplied scoping review, verify before use)

Geist1990, damped-driven Toda chain: relative CPU time GS 1.00, RGS 1.02, H1 1.07, H2 1.17, continuous QR 40; discrete QR variants agree on lambda_1 and lambda_28 within finite-time fluctuation. Bibliographic record verified (Geist1990); the numbers themselves are UNVERIFIED (full text not accessible on 2026-10-12, see `chaos-bibliography.md`). Confirm against the paper in literature-review/03-verify before any use; even then a reproduction target only, never a test oracle.

## E. Precision variants (apply to every row)

Float64 (default); Double64 (DoubleFloats.jl) or MultiFloats; BigFloat for short oracle runs. Precision is a run parameter recorded in the manifest, never a code fork.

## F. Built-in checks every estimator reports

symplectic pairing residual eps_pair(n); eps_det(n) for maps; orthogonality loss ||Q^T Q - I||; sum of exponents (zero for Hamiltonian); agreement of SALI decay with lambda_1 - lambda_2 (when both are run).
