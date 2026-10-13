# Spec: standard-map

Owner checkpoint 2026-10-13 (decision D24): Float32 (binary32) is dropped from M2. Every Float32 item in this spec (element-type lists, the screening tier, the Float32 GPU gate and Float32 test cases) is inactive until revisited at M4; stage 02 writes no Float32 tests.

Library layer: Systems. Draft for the 01-spec checkpoint, branch m2-p1-modules. No Julia type, function or file names appear outside the Compute and Oracle sections (Pattern 10).

## Scope

- Catalog rows (`_config/catalogs/systems.md`, section A): Y1 (Chirikov-Taylor standard map, lifted and wrapped coordinates, dimension 2) and Y5 (linear symplectic maps: rotation and hyperbolic, dimension 2N).
- Prerequisites: none. This module imports nothing from Foundation (task brief; `_config/architecture.md` "Layering"). Its consumers are tangent-core, precision and the diagnostics and Lyapunov modules.
- The module supplies, for each system: one forward step, one analytic tangent step, one inverse step, the invariant residuals of the Invariants section, and the known answers of this spec. It performs no health checking (NaN/Inf, pairing) and no renormalization; those belong to tangent-core and Diagnostics.
- Domain: state any finite real pair; parameter K any finite real, including 0 and negative values. Non-finite input gives non-finite output; it is not an error here.
- Out of scope: Henon map (Y2), coupled maps (Y3, Y4), indicators, estimators, element-type policy (precision spec). The criteria here use the element types named by owner decision, M2 kickoff: binary64 (default), double-double (about 32 digits), 256-bit arbitrary precision (CPU oracle only), and binary32 only if the precision spec enables it.

## Equations

Notation. eps_T is the spacing of element type T at 1 (binary64: 2^-52 = 2.22e-16; double-double: nominal 2^-104 = 4.9e-32; binary32: 2^-23). Omega_2 = [[0,1],[-1,0]] in the ordering (x, p). The block form of Omega_2N and the symplectic condition M^T Omega M = Omega are those of Skokos2010, J_2N p.7 and Eq. (5) p.8.

1. Canonical map (`_config/research-question.md`, "Ground truth for variational methods (standard map)"), lifted coordinates, x and p real and unbounded:
   p' = p + K sin x,   x' = x + p'.
2. Wrapped view, display only: x_w = x - 2 pi floor(x / 2 pi) in [0, 2 pi), same for p. The map commutes with the shifts (derived): if (x, p) -> (x', p') then (x + 2 pi a, p + 2 pi b) -> (x' + 2 pi (a + b), p' + 2 pi b) for integers a, b, because sin has period 2 pi and p' - p does not depend on a, b. So the map is well defined on the torus (x mod 2 pi, p mod 2 pi) and a lifted orbit, once wrapped, is the torus orbit.
3. Analytic tangent map, ordering (x, p) (derived): the entries are dp'/dx = K cos x, dp'/dp = 1, dx'/dx = 1 + dp'/dx = 1 + K cos x, dx'/dp = dp'/dp = 1. With c = K cos x:
   DT = [[1 + c, 1], [c, 1]],   det DT = (1 + c) - c = 1.
   This equals the matrix in `_config/research-question.md`, "Ground truth". It depends on x only through cos x, so the lifted and wrapped tangent maps coincide. Because det = 1 and the dimension is 2, DT^T Omega_2 DT = det(DT) Omega_2 = Omega_2 (derived: the (1,2) entry of DT^T Omega_2 DT is the determinant, the diagonal entries vanish).
4. Inverse map (derived): from x' = x + p' get x = x' - p'; from p' = p + K sin x get p = p' - K sin x. So
   x = x' - p',   p = p' - K sin(x' - p').
   This is the time-reversed form printed in DatserisDiakonosSchmelcher2015 Eq. (9) p.15 with theta = x, p = p (checked: theta_{n+1} = theta_n - p_n, p_{n+1} = p_n - k sin theta_{n+1}). Inverse tangent at the image point, with c evaluated at the preimage x: DT^-1 = [[1, -1], [-c, 1 + c]] (derived: inverse of [[a,1],[c,1]] with a - c = 1 is [[1,-1],[-c,a]]).
5. Parameter symmetry (derived): sin(x + pi) = -sin x, so T_{-K}(x + pi, p) = (x' + pi, p') where (x', p') = T_K(x, p). Spectra and finite-time exponents of orbit (x, p) at K equal those of (x + pi, p) at -K.
6. Conversion of source forms to the canonical form (derived; K is unchanged in all of them):
   - Unit-torus form with x_u, y mod 1: y' = y + (K / 2 pi) sin(2 pi x_u), x_u' = x_u + y' (daSilvaMancheinBeims2018 Eq. (1) p.2 with the noise amplitude D = 0; ContopoulosHarsoula2008 Eq. (1) p.3). Set x = 2 pi x_u, p = 2 pi y: multiply both equations by 2 pi, giving p' = p + K sin x, x' = x + p'. "mod 1" becomes "mod 2 pi". Fixed point (1/2, 0) becomes (pi, 0); (0, 0) stays (0, 0). Jacobian traces are unchanged (a diagonal rescaling of both variables by the same factor is a similarity).
   - Form with theta, p mod 2 pi (DatserisDiakonosSchmelcher2015 Eq. (1) p.4): p' = p + k sin theta, theta' = theta + p'. Identical with x = theta, K = k.
   - Sign-flipped unit-torus form of Meiss2015 Eq. (4) p.3 (benchmarks row "Confinement and escape time of a sticky orbit"): y' = y - (k / 2 pi) sin(2 pi x_M), x_M' = x_M + y'. Set x = 2 pi (x_M + 1/2), p = 2 pi y, K = k. Then sin x = -sin(2 pi x_M), so p' = p + K sin x and x' = x + p'. The hyperbolic fixed point (0.5, 0) of that paper is (0, 0) here (trace 2 + K > 2, consistent with "hyperbolic").
   - McMillan form (SalaArtuso2015 Eq. 5, benchmarks row "Lyapunov exponent, single chaotic orbit"): conversion p = x - y is written in that row and in the precision spec, Equations item 2; not repeated here.
7. Fixed points on the torus (derived). p' = p mod 2 pi forces K sin x = 2 pi m for an integer m; x' = x mod 2 pi forces p' = 0 mod 2 pi, so p = 0 mod 2 pi. With m = 0: (0, 0) and (pi, 0), any K. With m != 0 (needs |K| > 2 pi |m|): x with sin x = 2 pi m / K, the accelerator-mode points; the lifted orbit then moves by 2 pi m per step in both variables.
   Trace of DT at a fixed point: 2 + K cos x. At (0, 0): 2 + K. At (pi, 0): 2 - K. At an m != 0 point: 2 +- sqrt(K^2 - (2 pi m)^2). Eigenvalues lambda_pm = (tr +- sqrt(tr^2 - 4)) / 2, product 1. Elliptic if |tr| < 2, hyperbolic if |tr| > 2, parabolic (Jordan block) if |tr| = 2.
   - (0, 0): elliptic for -4 < K < 0, hyperbolic for K > 0 or K < -4, parabolic at K = 0 and K = -4.
   - (pi, 0): elliptic for 0 < K < 4, hyperbolic for K < 0 or K > 4, parabolic at K = 0 and K = 4.
   - Exponent of a hyperbolic fixed point for a generic tangent vector: ln |lambda_+|, with |lambda_+| = (|tr| + sqrt(tr^2 - 4)) / 2. At a parabolic point the tangent power grows linearly in n and the finite-time exponent decays like ln n / n (see item 9).
   - Check against the benchmark row "Stability of the central fixed point at (1/2, 0)" (Quantity: "Elliptic for K < 4, hyperbolic for K > 4; trace of the Jacobian is 2 -/+ K"): with item 6, (1/2, 0) is (pi, 0), trace 2 - K; the row's reading holds for K > 0. At K = 4 exactly, tr = -2 and the point is linearly unstable by algebraic growth with exponent 0; this agrees with ContopoulosHarsoula2008 p.3 ("unstable for K >= 4") and is why the row says hyperbolic only for K > 4.
8. Eigenvectors of (0, 0), K > 0 (derived): the first row-independent equation of DT v = lambda v is K v_x + v_p = lambda v_p, so v_p / v_x = K / (lambda - 1). With s = sqrt(K (K + 4)), lambda_pm = (2 + K +- s) / 2: unstable direction (v_x, v_p) proportional to (1, (s - K) / 2); stable direction proportional to (1, -(K + s) / 2). Benchmarks row "Eigenvectors of the fixed point (0,0)" (DatserisDiakonosSchmelcher2015 Eqs. (7), (8) p.15) prints unstable ((-k + s)/2, 1) and stable (-1, (-k + s)/2), and at k about 1 (0.62, 1) and (-1, 0.62). Reading the printed pairs as (p, theta) (not (theta, p) as the row's note assumes) the unstable vector matches mine for every k; the stable vector matches only near k = 1 (at K = 2 the printed formula gives p/theta = -1.366, mine gives -2.732). So the printed stable formula is used only at k about 1, as a loose check (criterion 12).
9. Known answer K = 0 (derived): the map is p' = p, x' = x + p; DT = [[1,1],[0,1]] for every orbit, DT^n = [[1, n], [0, 1]]. For w(0) = (0, 1): ||w(n)|| = sqrt(n^2 + 1) and the finite-time exponent X(n) = (1/n) ln(||w(n)|| / ||w(0)||) = ln(sqrt(n^2 + 1)) / n, so n X(n) / ln n = 1 + O(1 / (n^2 ln n)). For w(0) = (1, 0) (eigenvector): X(n) = 0 exactly. The largest singular value of DT^n is (n + sqrt(n^2 + 4)) / 2. Definition of X: Skokos2010 Eq. (54) p.31. The law ||w|| proportional to t and X = O(ln t / t) for regular orbits: Skokos2010 Eq. (62) and the text after it, p.35 (derived there for flows in action-angle form; p.34 states a similar analysis holds for symplectic maps); the general form X(t) = chi + (b + z(t)) / t is Eq. (63) p.36.
10. Y5, linear symplectic maps, dimension 2N, ordering (q_1..q_N, p_1..p_N), Omega_2N = [[0, I], [-I, 0]]; state z, map z' = A z, tangent map equal to A for every state. Building blocks per degree of freedom j (derived):
   - Rotation by angle theta_j: R = [[cos th, sin th], [-sin th, cos th]], det 1, eigenvalues e^(+-i th), exponents 0 and 0.
   - Hyperbolic with parameter mu_j > 0: H = [[cosh mu, sinh mu], [sinh mu, cosh mu]], det = cosh^2 - sinh^2 = 1, eigenvalues e^(+-mu) (eigenvectors (1, +-1)), H symmetric so the singular values of H^n are e^(+-n mu).
   - Reflected hyperbolic -H: eigenvalues -e^(+-mu), same exponents (exponent is ln |eigenvalue|).
   A = direct sum of N such blocks, each acting on the pair (q_j, p_j) (so a coordinate permutation of the ordering above). Exact spectrum: {+mu_j and -mu_j for each hyperbolic block, 0 twice for each rotation block}, sorted descending; this is lambda = ln |eigenvalue| (catalog row Y5).
   - Non-normal variant: A_S = S^-1 A S with S = [[I, B], [0, I]], B symmetric N x N (symplectic: S^T Omega S = Omega, derived by multiplying out blocks, which uses B^T = B). A_S is symplectic, has the same eigenvalues and exponents; its matrix is not orthogonal-structured, so QR-based estimators are exercised nontrivially.
   - Inverse (derived from M^T Omega M = Omega): A^-1 = Omega^-1 A^T Omega.

## Algorithm

1. Forward step: item 1 of Equations, in the lifted variables. A wrapped-evolution mode applies the wrap of item 2 to x and p after every step; it is for the wrapped-vs-lifted study only and is not the reference.
2. Tangent step: item 3. The standard map's tangent step is evaluated at the pre-step x, the same x that the forward step uses.
3. Inverse step: item 4.
4. Y5: z' = A z; tangent A; inverse via Omega^-1 A^T Omega or the stored inverse blocks.
5. Parameters (K, mu_j, theta_j, B) are inputs, stored in the element type of the state. Constants that are not exactly representable (pi) are built in that type (precision spec, Algorithm item 2).
6. No coefficient tables. No numerics-defaults.toml key is needed by this module; no value is proposed.

## Invariants

Residuals the module exposes for any state (floating-point values in the element type T):
- I1: |det DT - 1|. Standard map and Y5.
- I2: max entry of |DT^T Omega DT - Omega|.
- I3: lifted vs wrapped: wrapped(lifted orbit) agrees with the torus orbit (item 2).
- I4: round trip: inverse(forward(z)) = z, and forward(inverse(z)) = z.
- I5: K -> -K, x -> x + pi symmetry (item 5).
- I6: the exact spectrum (item 10) and the fixed-point eigenvalues (item 7) are reproduced by the tangent products.
Error budget used by every orbit comparison below (derived, first-order). Let sigma(K) = (S + sqrt(S^2 - 4)) / 2 with S^2 = 2 K^2 + 2 |K| + 5. This bounds the largest singular value of DT for every x: the squared Frobenius norm of DT is 2 c^2 + 2 c + 3 with |c| <= |K|, maximal at c = |K|, and sigma_1 sigma_2 = 1 gives (sigma_1 + sigma_2)^2 = ||DT||_F^2 + 2. Let Gamma_n = sigma(K)^n, a bound on the growth of any perturbation over n steps (the inverse has the same singular values). A rounding error of relative size eps_T in one step, made at step k, is amplified by at most Gamma_(n-k) afterwards, so after n steps the orbit difference in the max norm is at most C eps_T M_n n Gamma_n with C = 16 (each step contributes about 6 rounding errors of the sin, product, two additions and the wrap; sin and cos from a GPU library may add 2 ulp), where M_n = max(1, |K|, max over k <= n of |x_k|, |p_k|) taken on the lifted orbit. Write tau(n) = 16 eps_T M_n n Gamma_n. The horizon n_H(K) is the largest n <= 20 with Gamma_n <= 10^6 (so tau stays far below 1 for binary64). The distance d used for torus comparisons is the maximum over x and p of the absolute difference reduced to (-pi, pi]; for lifted comparisons it is the plain absolute difference.

## Compute

Policy: `_config/architecture.md` "Compute" (decision D17); machine: GTX 1050 3 GB, compute capability 6.1, weak binary64 throughput (risk K15); CI has no GPU.
- Data parallel, GPU candidates: ensembles of independent orbits (grids of initial conditions, K sweeps, indicator maps). Orbits share no data; steps within an orbit are sequential. One KernelAbstractions kernel serves CPU and GPU. The step is light (one sin, one cos, about 6 additions and multiplications), so binary64 on cc 6.1 may be compute-bound; no GPU speedup is assumed here. The 04-verify report measures binary64 CPU threads vs GPU and the module ships GPU as default only if the cross-backend test passes.
- CPU only: a single long orbit, and arbitrary-precision (BigFloat, not isbits, so a GPU request falls back with the recorded reason "kernel does not compile for the element type"). Double64 (DoubleFloats, isbits) is attempted on the GPU; if the kernel does not compile for it, the run is on CPU with the same reason. Float32 is a GPU candidate for screening ensembles only if the precision spec enables it; its error budget is the Invariants error budget with eps = 2^-23.
- Memory per orbit (structure-of-arrays): 2 state numbers plus 4 tangent numbers (a 2 x 2 tangent matrix, or two deviation vectors of 2 numbers each) = 6 numbers: 48 bytes in binary64, 96 bytes in double-double, 24 bytes in binary32. For orientation 10^7 orbits is 480 MB in binary64; actual chunk sizes derive from free memory at run time (policy), not from this figure.
- Reproducibility: initial conditions of orbit j come from a random stream keyed by j, never by thread or device (policy); the CPU result for orbit j is bitwise identical with 1 and with N threads.
- CPU-vs-GPU agreement: not bitwise (chaotic orbits amplify FMA and library differences); the bound is tau(n) of the Invariants section for n <= n_H(K), max over the ensemble. Longer runs are compared at the estimator tolerances owned by the Lyapunov and Indicators specs.

## Oracle

In the package test suite (analytic, invariant, cross-method and small self-contained tests):
- O1 Analytic tangent vs derivative of the map at seeded random points: complex-step derivative (step 1e-30) and ForwardDiff, both applied to the forward map, compared with the analytic tangent. Seed from `_config/numerics-defaults.toml` (seed = 20261012), generator StableRNGs.
- O2 The same forward step evaluated in Float64, Double64 (DoubleFloats) and 256-bit BigFloat from identical exactly representable data, short horizon n_H(K).
- O3 Analytic known answers: K = 0 shear, fixed points, Y5 exponents, inverse round trip, symmetry, conversion checks from published fixed-point data (items 6 to 9).
- O4 CPU backend of the kernel vs a plain loop (stands in for the GPU in CI).
In the separate verification environment driven from module-build stage 04 (not in the package test suite):
- O5 ChaosTools standard map system as a cross-check of orbits: same seeded initial conditions, binary64, compared after wrapping on the torus over n <= n_H(K) (criterion 19). The convention of its parameter and angle ordering is confirmed against the canonical form at stage 04 before comparing; this spec does not assert it. Pass bound tau(n); a longer run (n = 100) is an unscored sanity line only.
- O6 ForwardDiff of the n-step map vs the product of analytic tangent matrices for n <= n_H(K) at 10^4 seeded points (heavy), and the GPU-vs-CPU comparison on the GTX 1050 (result recorded in the 04-verify report).
- O7 256-bit BigFloat reference orbits at the larger sample, and a measurement of the effective double-double accuracy of sin and cos (see criterion 14).
Published lambda(K) values are loose sanity checks only (risk K18); none is an acceptance criterion here. In particular the benchmark "lambda ~ ln(K/2)" (Chirikov law, K >> 1) is left to the Lyapunov modules.

## Acceptance criteria

Seeded random states: x and p uniform in [-4 pi, 4 pi], K from {0, +-0.5, +-1, +-2, +-4, +-5, +-6, +-10, +-20} and uniform in [-20, 20]; 1000 points per check unless stated. eps_T of the type under test; tau(n) and n_H(K) as defined in Invariants.
1. Tangent vs derivative (O1): every entry of the analytic DT differs from the derivative of oracle O1 (complex-step and automatic differentiation, both) by at most 8 eps_T (1 + |K|), at all points. Pass: all points.
2. Determinant and symplecticity: |det DT - 1| <= 8 eps_T (1 + |K|) and every entry of |DT^T Omega_2 DT - Omega_2| <= 8 eps_T (1 + |K|), at all points of criterion 1, in each of binary64 and double-double.
3. Inverse round trip: for n = 1 and n = n_H(K), n forward steps then n inverse steps return to the start with plain (lifted) difference <= tau(n); and inverse then forward likewise. Pass: all 200 orbits per K in the listed set.
4. Inverse tangent: max entry of |DT(z) DT_inv(T z) - I| <= 8 eps_T (1 + |K|)^2 at all points.
5. Lifted vs wrapped (I3): over n <= n_H(K), the lifted orbit wrapped at the end and the wrapped-evolution orbit differ in torus distance by <= tau(n), for 200 orbits per K in the listed set (M_n taken on the lifted orbit).
6. Symmetry (I5): the orbit of (x + pi, p) at -K, shifted back by pi, agrees with the orbit of (x, p) at K in torus distance <= tau(n), n <= n_H(K), with pi built in the element type.
7. K = 0 exactness: for x_0, p_0 multiples of 2^-10 with magnitude < 8, 1000 steps give x_n = x_0 + n p_0 and p_n = p_0 with zero difference (exact, 0 ulp); and the product of n one-step tangent matrices equals [[1, n], [0, 1]] exactly for every n in {1, 2, 10, 100, 1000, 10^4}.
8. K = 0 finite-time exponent, X(n) formed from the exact n-step product of criterion 7 applied to w(0) (one norm, one logarithm, no renormalize-and-accumulate path, which would carry about n eps): for w(0) = (0, 1), |X(n) - ln(sqrt(n^2 + 1)) / n| <= 8 eps_T for n in {1, 10, 100, 10^3, 10^4, 10^6}, and |n X(n) / ln n - 1| <= 10^-9 at n = 10^4 (expected 5.4e-10 from item 9); for w(0) = (1, 0), X(n) = 0 exactly for the same n. The decay is the O(ln t / t) law of Skokos2010 p.35, Eq. (62).
9. Fixed-point classification and eigenvalues (item 7): for K in {-6, -2, 1, 3, 4.5, 6} (all away from parabolic values) at (0, 0) and at (pi, 0) built in the element type, the class (elliptic or hyperbolic) is correct, tr(DT) = 2 + K cos x to 8 eps_T (1 + |K|), and |lambda_+| (hyperbolic) or the rotation angle (elliptic, cos phi = tr / 2) matches the closed form to 32 eps_T (1 + |K|). At K = 4, point (pi, 0), and at K = 0: DT is a Jordan block and DT^n has maximum entry 4n (K = 4, derived: DT = -I + N with N^2 = 0, so DT^n = (-1)^n (I - nN), entries 1 + 2n, -n, 4n, 1 - 2n) and n (K = 0); the ratio of the maximum entry to n is 4 and 1 within 1e-12 for n in {10, 100, 1000}.
10. Fixed point is exact: starting at (0, 0), K any listed value, 10^4 steps stay at (0, 0) with zero difference.
11. Benchmark trace row: for K in {1, 3, 3.9, 4.1, 5}, tr(DT) at (pi, 0) equals 2 - K to 8 eps_T (1 + |K|), is elliptic for K < 4 and hyperbolic for K > 4; tr(DT) at (0, 0) equals 2 + K (source forms converted by item 6).
12. Datseris eigenvectors at K = 1: derived unstable direction has v_p / v_x = (s - K) / 2 = 0.618 and the stable direction has v_p / v_x = -1.618, each to 10^-12; with the printed pairs read as (p, theta), the printed small-to-large component ratio 0.62 agrees with the derived magnitude 0.618 (unstable) and 1/1.618 = 0.618 (stable) to 0.005 (rounding of the printed digits). Loose check only.
13. Contopoulos orbit data (benchmarks row "Island centre and unstable 4/9 orbit point"), K = 5, unit-torus coordinates converted by item 6: the island centre O1 = (0.68002, 0.36002) returns to itself modulo 1 after 2 steps (ContopoulosHarsoula2008 Sec. 2 p.4: O1 and O1' form a period-2 family bifurcated from (1/2, 0) at K = 4); the point P1 = (0.642118, 0.369344) of the unstable 4/9 orbit returns after 18 steps and not after 9 (Sec. 3 p.5: nine points plus nine symmetric ones, multiplicity 18). Tolerance 1e-4 in each unit coordinate, i.e. 2 pi 1e-4 in canonical coordinates. Basis: the printed 6 and 5 digits (rounding 5e-7 and 5e-6) amplified by the monodromy (largest eigenvalue magnitude about 13 for P1 over 18 steps, 1 for the elliptic O1 over 2 steps; derived in a scratch run). Residuals observed while drafting (binary64): below 1e-4 for O1 at n = 2 and about 5e-6 for P1 at n = 18; P1 after 9 steps is off by 0.28 and 0.26. Loose sanity check only (risk K18 spirit).
14. Element-type agreement (O2): binary64 vs 256-bit, double-double vs 256-bit, and binary32 vs 256-bit if enabled, from identical exactly representable data, K in the listed set, 200 orbits each, n = n_H(K): torus distance and lifted difference <= tau(n) with eps_T of the lower precision. For double-double eps is the effective value measured in O7 (largest relative error of sin and cos on the seeded sample, not less than 2^-104 = 4.9e-32); pass requires the effective value <= 2^-100, else report to the owner. For binary32 use n = largest n <= 20 with Gamma_n <= 100.
15. Y5 symplecticity: for N in {1, 2, 3, 5}, random blocks (rotation, hyperbolic, reflected hyperbolic, mu_j in [0.05, 0.5], theta_j in [0, 2 pi)), A and A_S satisfy max |A^T Omega A - Omega| <= 16 eps_T ||A||^2 and |det A - 1| <= 16 eps_T ||A||^N.
16. Y5 exact exponents: for block-diagonal A, n in {10, 50}: |(1/n) ln sigma_max(A^n) - lambda_1| <= 16 eps_T and |(1/n) ln sigma_max(A^-n) + lambda_2N| <= 16 eps_T, with A^-1 from the inverse of item 10 and lambda_1, lambda_2N from the exact spectrum. For A_S the bound is 16 eps_T + ln(kappa(S)) / n with kappa(S) = ||S|| ||S^-1|| (derived: sigma_max(S^-1 D S) is within a factor kappa of sigma_max(D)). Pass: all N and parameter draws.
17. Y5 inverse: A^-1 A = I and A A^-1 = I to 16 eps_T ||A||^2 and the formula Omega^-1 A^T Omega agrees with the inverse of the block-diagonal form to the same bound.
18. Cross-backend (oracle O4, and the GPU comparison of O6): 10^5 seeded orbits at K in {1, 5}, n <= n_H(K): CPU-backend kernel vs plain loop, and GPU vs CPU, max torus distance <= tau(n). Result for orbit j is bitwise identical with 1 and with 8 CPU threads. The double-double kernel either meets the same bound on the GPU or the run records the reason "kernel does not compile for the element type" and runs on the CPU.
19. Oracle O5: the external package's standard map orbits vs this module for 100 seeded orbits, K in {0.5, 1, 2, 5}, n <= n_H(K): torus distance <= tau(n), binary64. If the convention check shows a different parameter or ordering, the comparison is made after the stated conversion and that conversion is written into the 04-verify report.

## Sources

Verified keys used (all in `_config/bibliography/references.bib`, none marked `?`):
- Skokos2010 (arXiv:0811.0882 pages): J_2N block form p.7; symplectic map Eq. (4) and condition Eq. (5) p.8; finite-time exponent Eq. (54) p.31; remark that the regular-orbit analysis carries to symplectic maps p.34; norm growth Eq. (62) and X = O(ln t / t) p.35; Eq. (63) p.36.
- daSilvaMancheinBeims2018: Eq. (1) p.2 and fixed points (0,0), (1/2,0) with trace |2 +- K| p.2 (read in a text extraction of the PDF).
- ContopoulosHarsoula2008: Eq. (1) p.3 and stability statement p.3; period-2 family (O1, O1') p.4 Sec. 2; orbit 4/9 points, multiplicity 18, island centre p.5, Sec. 3.
- DatserisDiakonosSchmelcher2015: Eq. (1) p.4; Eqs. (7), (8), (9) p.15.
- Meiss2015: Eq. (4) p.3 (form and caption only).
- Layer 3: `_config/research-question.md` "Ground truth"; `_config/bibliography/benchmarks.md` rows "Stability of the central fixed point at (1/2, 0)", "Eigenvectors of the fixed point (0,0)", "Island centre and unstable 4/9 orbit point", "Confinement and escape time of a sticky orbit" (form only), "Lyapunov exponent, single chaotic orbit" (form only); `_config/catalogs/systems.md` rows Y1, Y5; owner decisions, M2 kickoff (double-double extended precision, arbitrary precision for short CPU oracle runs, oracles in a separate environment); decision D17; risks K15, K18.
- Chirikov1979 and Greene1979 are listed for Y1 in the catalog but were not read here and are not cited for any equation or number.
Findings for the owner: (a) the Datseris stable-eigenvector formula (Eq. (7)) disagrees with the derived eigenvector except near k = 1, and the printed components follow the order (p, theta); (b) the 4/9 orbit point P1 has period 18, not 9, and O1 has period 2, which the benchmarks row does not say; (c) double-double sin and cos accuracy is unmeasured.

## Amendment D28 (2026-10-13)

Decision D28: long runs store the state wrapped (exact reduction modulo 2 pi each step) plus integer winding counts per coordinate. This amendment adds the wrapped system with counts as a second map next to the lifted one of Equations item 1. The lifted map is unchanged and stays the reference for the wrapped-vs-lifted sensitivity study. Pattern 10 applies: Julia names appear only in the Compute subsection of this amendment.

### Equations (wrapped system with counts)

Notation. tp is 2 pi built in the element type T (a rounded number, not the real 2 pi). State (x, p, m, n): x, p in [0, tp), m and n integers stored exactly in T. Lifted coordinates X = x + 2 pi m and P = p + 2 pi n. The phase-space (tangent) dimension stays 2; the state has 4 components.

W-E1. Reduction (derived). For a finite u, red(u) = (r, a) with a = floor(u / tp) and r = u - tp a, followed by two guards in this order: if r < 0 then r += tp and a -= 1; if r >= tp then r -= tp and a += 1. Why: floor(u / tp) can be one too large when u / tp rounds up to an integer (r is then a tiny negative number, the first guard fires, and r + tp may round to tp, which the second guard maps to 0 with the count restored); it can be one too small when tp a rounds low (r >= tp, second guard). Result: r in [0, tp) and u = r + tp a up to one rounding. For a non-finite u the result is non-finite and the count is meaningless (W-I5).

W-E2. Forward step (derived from Equations item 1):
- p~ = p + K sin x; (p', a) = red(p~); n' = n + a.
- x~ = x + p'; (x', b) = red(x~); m' = m + n' + b.
Check: P' = P + K sin X = p + 2 pi n + K sin x = p~ + 2 pi n = p' + 2 pi (n + a) = p' + 2 pi n'. X' = X + P' = x + p' + 2 pi (m + n') = x~ + 2 pi (m + n') = x' + 2 pi (m + n' + b) = x' + 2 pi m'. So (x', p', m', n') reproduces the lifted map exactly in real arithmetic. The arguments are bounded: |p~| <= tp + |K| and x~ < 2 tp, so each reduction costs one rounding of size O(eps_T), the same order as one rounding of the lifted step.

W-E3. Tangent step (derived): the matrix of Equations item 3 evaluated at x, DT = [[1 + c, 1], [c, 1]] with c = K cos x. The lifted point enters only through cos, which has period 2 pi, so the wrapped and lifted tangents coincide in real arithmetic; in T the cos argument is x (below tp) instead of X, so it is at least as accurate as in the lifted mode. det DT = 1. The counts carry no tangent component.

W-E4. Inverse step (derived from Equations item 4). From (x', p', m', n'): x~ = x' - p' (in (-tp, tp)); (x, b) = red(x~); m = m' - n' + b. Then p~ = p' - K sin x; (p, a) = red(p~); n = n' + a. Check: X = X' - P' = x~ + 2 pi (m' - n') = x + 2 pi (b + m' - n') = x + 2 pi m. P = P' - K sin X = p' - K sin x + 2 pi n' = p~ + 2 pi n' = p + 2 pi (a + n') = p + 2 pi n. Inverse tangent at the image point (derived): DT^-1 = [[1, -1], [-c, 1 + c]] with c = K cos(x' - p') (cos of the unreduced preimage coordinate, equal to cos x up to the period).

W-E5. Torus conversion and lift (derived). From a lifted point (X, P): (x, m) = red(X) and (p, n) = red(P). The lift of a state is (x + tp m, p + tp n) in T. In T the lift has absolute error about |m| tp eps_T, so for |m| beyond 1 / eps_T it is meaningless in T and is formed in a wider type from the exact counts.

W-E6. Circle length (derived). The wrapped system in T is the exact map on a circle of circumference tp, not 2 pi. Each wrap shifts the point by |2 pi - tp| <= pi eps_T, so after n steps with w wraps in total the cumulative shift is at most w pi eps_T. w is at most (M_n / 2 pi) + 2 n, so the shift is O(M_n eps_T) and lies inside the error budget tau(n) of Invariants (whose factor 16 and growth Gamma_n are far larger).

### Invariants (wrapped system)

- W-I1: after n <= n_H(K) steps, x and p are in [0, tp) and m, n are integers, and the torus distance of (x, p) to the lifted orbit wrapped at step n is at most tau(n) (M_n on the lifted orbit).
- W-I2: the lift of the state, computed in arbitrary precision from the exact counts, differs from the lifted orbit by at most tau(n) in the plain max norm.
- W-I3: counts are exact integers while |m| and |n| stay below the largest integer T stores exactly (2^53 for binary64, 2^106 for double-double, 2^precision for arbitrary precision); every update m' = m + n' + b and n' = n + a is then an exact addition.
- W-I4: inverse(forward(s)) = s: x, p to tau(n), counts exactly. A count can flip by one only when a coordinate lies within tau(n) of 0 or tp, a null set for continuous seeded data; a flip is reported with the orbit and step, not tolerated.
- W-I5 (health): if |m'| or |n'| reaches the limit of W-I3, forward and inverse return NaN in all four components, so the health check of tangent-core fires; a silently wrong count is not allowed. Non-finite input gives non-finite output as in Scope.
- W-I6: K = 0 is exact: p and n never change, and m_k = floor((x0 + k p0) / tp) in exact arithmetic when x0, p0 are dyadic and the sums are exact in T.

### Compute (amendment)

Per orbit: 4 state numbers plus 4 tangent numbers (a 2 x 2 matrix, k = 2 deviation vectors of 2 numbers) = 8 numbers: 64 bytes in binary64, 128 bytes in double-double, before reducer state and the 16-byte index. With tangent-core's formula and state length s = 4, d = 2, k = 2: sizeof(T) (4 + 4 + 4) + 16 = 112 bytes in binary64. Cost per step: one sin, two floors and a few additions, about the lifted cost.

Interface. System `WrappedStandardMap{T}(K)`, state `SVector{4,T}` = (x, p, m, n), exported with `torus_state(z)` (lifted pair to state) and `lift(s)` (state to lifted pair in T). `dimension` returns 2. `forward`, `inverse`, `tangent`, `inverse_tangent` use the lifted system's names; the tangents are 2 x 2. tp is `2 * pi` built in T (Systems imports nothing from Foundation, decision D23). The exactness limit is `maxintfloat(T)`.

### Acceptance criteria (amendment)

Seeded data: x0, p0 uniform in [0, tp) drawn in Float64 and cast to T (the lifted and wrapped orbits start from identical data), m = n = 0, K in {0.5, 2, 5, 10}.
- W1 Torus agreement with the lifted map: for every n <= n_H(K), 200 seeded orbits per K, Float64 and Double64 (effective eps as in criterion 14), the torus distance between the wrapped state (x, p) and the lifted orbit is at most tau(n).
- W2 Lifted reconstruction: the lift of the state, with the state widened exactly to 256-bit arithmetic and 2 pi taken at 256 bits, differs from the lifted-mode orbit by at most tau(n) (plain max norm) for n <= n_H(K), same orbits as W1.
- W3 Counts exact and bounded: after 10^6 steps at K = 5 (Float64), x and p are in [0, tp) and m, n are integer-valued. The counts are recomputed independently at 256 bits for the first 50 steps. Free-running a 256-bit orbit for 50 steps cannot be compared with the binary64 orbit (chaotic divergence: the exponent at K = 5 is about 0.9 per step, so a 1e-16 difference reaches order 1 after about 40 steps, see the impl notes for the measured step); the recomputation is therefore teacher-forced: at each of the first 50 steps the 256-bit increments a and b are computed from the exact binary64 state of that step and compared with the increments the binary64 map produced (n' - n and m' - m - n'). Additionally the free-running 256-bit orbit is compared for n <= n_H(K) in x, p (tau(n)) and counts (exact).
- W4 K = 0 known answer: p and n never change; for k <= 1000, m_k equals floor((x0 + k p0) / (2 pi)) computed at 256 bits, with x0 and p0 dyadic multiples of 2^-10 in [0, 2 pi). The binary64 state is exact in these sums, and the 256-bit floor differs from the binary64 floor only if x0 + k p0 lies within about 2 pi eps of a multiple of 2 pi; with dyadic data of 2^-10 spacing the distance of such a sum from a multiple of the irrational 2 pi is far above that (checked for the seeded cases in the test plan).
- W5 Tangent agreement: along the wrapped orbit for n <= n_H(K), the product of the tangent matrices satisfies the tangent-core criterion 4 bound (8 n eps_T times the product of Frobenius norms times ||W_0||_F against the 256-bit product at the same realized points), and differs from the lifted-mode tangent product by at most twice that bound plus the first-order orbit term n |K| tau(n) Gamma_n ||W_0||_F (derived: the orbits differ by at most tau(j) in x after j steps, so each tangent entry K cos x differs by at most |K| tau(j); a perturbation of one factor is amplified by at most Gamma over the other steps).
- W6 Inverse round trip: n forward then n inverse steps recover the state for n = 1 and n = n_H(K): x, p within tau(n) (torus distance) and counts exactly equal, 200 orbits per K, Float64 and Double64; likewise inverse then forward.
- W7 Zero allocation: forward, tangent, inverse and the per-orbit propagation of tangent-core with the wrapped system allocate 0 bytes after warm-up, Float64 and Double64.
- W8 Overflow guard: a state with m one below the binary64 exactness limit 2^53 and a step that increases m (the new x~ is at least tp) gives NaN in all four components; likewise for n; propagation of that state raises the numerical health error with kind "state"; a state with m two below the limit stays finite with the exact count.
