# Chaos: Research Question and Study Design

Layer 3, canonical. Canonical home for the research framing. Method details live in the catalogs; citations live in the bibliography.

## Primary question

Under controlled numerical conditions, which Lyapunov-exponent (LE) algorithms and variational chaos indicators most reliably and efficiently distinguish regular, weakly chaotic, strongly chaotic and sticky trajectories in the Chirikov-Taylor standard map?

Core distinction: (i) algorithms that compute the same asymptotic LE spectrum versus (ii) chaos indicators that classify finite-time dynamics and need not converge to an exponent.

## Three claims to test

| # | Claim | Testable form |
|---|-------|---------------|
| C1 | Asymptotic equivalence | Benettin, QR (Gram-Schmidt / Householder), stabilized SVD and well-conditioned tangent-map methods agree on lambda_1 when run long enough at adequate precision. |
| C2 | Finite-time inequivalence | Finite-time LE estimators differ substantially over realistic budgets, most of all in weak chaos, sticky layers and near separatrices. |
| C3 | Diagnostic inequivalence | SALI, GALI, MEGNO, FLI/OFLI, RLI, frequency and recurrence methods operationalize "chaos" differently; classifications and detection times diverge. |

Headline to defend: algorithms targeting the same asymptotic spectrum agree only after a regime-dependent convergence time; finite-time indicators can classify faster but do not answer the same dynamical question.

## Ground truth for variational methods (standard map)

Map, with x and p taken mod 2*pi for display only (lifted coordinates internally):

    p_{n+1} = p_n + K sin(x_n)
    x_{n+1} = x_n + p_{n+1}

Analytic tangent map in (x, p) ordering:

    DT_K(x_n, p_n) = [[1 + K cos x_n, 1], [K cos x_n, 1]],   det = 1

Symplectic, so the 2D spectrum is {lambda, -lambda}. Any persistent violation is a numerical diagnostic, not physics. Chaotic SALI decays like exp(-(lambda_1 - lambda_2) n) = exp(-2 lambda_1 n), which gives an internal cross-check against lambda_1.

## Test panel

K regimes: K = 0 (integrable); small K (KAM tori plus resonances); near-critical golden-mean breakup (K near Greene's K_c, value taken from bibliography); intermediate K (mixed phase space); large K (strong chaos).

Orbit classes per regime: KAM torus; resonant (near-)periodic orbit; chaotic-sea orbit; sticky orbit near island boundary or cantorus; orbit near a hyperbolic fixed point or separatrix.

The panel is a dynamical-regime benchmark set, not one lambda_1 per K.

## Outputs retained

Per orbit-method pair: full time series, not final scalars.

- Finite-time MLE lambda_1(n); full spectrum lambda_i(n) where relevant
- Symplectic pairing residual eps_pair(n) = max_i |lambda_i(n) + lambda_{2N+1-i}(n)|
- Area-preservation residual eps_det(n) = |det DT_K - 1|
- SALI, GALI_k, MEGNO, FLI/OFLI, RLI trajectories
- Detection time: first time a method reaches a pre-registered, stable classification
- Cost: wall time, tangent-vector count, memory, precision
- Sensitivity to deviation-vector orientation and magnitude, K, initial point, length, precision, wrapped vs lifted coordinates

## Expected sources of disagreement (hypotheses, not bugs)

Finite-time effects and stickiness; weak chaos needing long horizons; O(1/n)-type decay of regular-orbit FTLEs versus fixed thresholds; alignment (SALI/GALI) versus growth-rate time scales; two-trajectory perturbation size (nonlinear saturation vs round-off) and wrapping; QR diagonal order dependence vs finite-time singular values; data-driven reconstruction of conservative maps (short, noisy series); non-transferable thresholds across K, orbit type and budget.

Thresholds are pre-registered per study in the study workspace's design stage, before any data is looked at.

## Paper sequence

| Paper | Scope | Needs (modules) |
|-------|-------|-----------------|
| P1 | LE algorithms on the standard map: Benettin, QR (CGS, MGS, Householder), stabilized SVD, two-trajectory. Finite-time disagreement, symplectic consistency, cost, precision. Each orbit gets two preregistered labels (idea 23, decisions D36, D37): orbit structure (periodic, only where known to be periodic such as the panel's fixed point, O1 and P1; quasiperiodic; irregular; unresolved) and instability (positive growth supported; none detected on the tested horizon; unresolved). Unresolved means the precision ladder (Float64, Double64, BigFloat), the estimators, or a designated subset re-run with several initial deviation directions disagree. Finite-time rates use the stated canonical norm. | standard map, tangent core, LE estimators, diagnostics, precision backends |
| P2 | LE estimators vs indicators: lambda_1(n), SALI, GALI_2, MEGNO, FLI/OFLI, RLI on grids and a sticky-orbit ensemble. Detection latency, classification stability. | P1 + indicators, FMA/NAFF comparator |
| P3 | Higher dimensions: coupled / Froeschle 4D maps, Henon-Heiles and other flows. Full spectra, GALI_k (k>2), pairing, CLVs, integrator effects. | ND maps, flows, symplectic integrator catalog, continuous QR, CLVs |
| P4 | Observation-limited and noisy settings: tangent "oracle" vs Wolf, Rosenstein, Kantz, Eckmann-Ruelle estimators. | time-series estimators, embedding |
| P5 | Shadow-Hamiltonian predictive closure (research line R2, below). | splitting integrators, BCH / modified-Hamiltonian module, forced pendulum, transport metrics, FMA |
| P6 (conditional on P5) | Transfer of R2 to planar restricted three-body, then secular planetary model, then Solar-System ensembles (REBOUND, REBOUNDx, celmech, ASSIST). | N-body adapters, Horizons fetch |

## Research line R2: physical vs numerical transport (predictive closure)

Prior-art audit (user-supplied, 2026-10-12; every source must pass literature-review/03-verify before citing): a broad "distinguish physical from numerical chaos" project is NOT novel. Backward error analysis already gives a modified Hamiltonian H_h = H + h^p H_p + h^(p+2) H_(p+2) + ... for symmetric order-p symplectic methods; resonance shifts and separatrix changes follow from perturbation theory; Solar-System studies 2019-2026 already document integrator-dependent precession and instability statistics. Indicators (MEGNO, FLI, SALI, GALI, LE, FMA) measure the numerically realized map; alone they cannot attribute chaos to H or to H_h.

Surviving question (conditional GO): in a low-dimensional Hamiltonian with a known separatrix and a measurable transport coefficient, does a truncated modified Hamiltonian of a specified splitting quantitatively predict the step-size-dependent excess in action / frequency diffusion and the displacement of the stochastic-layer boundary, relative to a high-accuracy reference, over a pre-specified finite time?

Conservative claim: test whether an explicitly computed delta H_h predicts excess finite-time transport near a specified resonance or separatrix, rather than documenting disagreement between methods.

Minimum viable experiment:
- System: forced pendulum, autonomous form H = p^2/2 - cos q + omega I + eps cos q cos phi (systems catalog F7).
- Arms: reference (high-order adaptive at stringent tolerance plus an independent high-order symplectic run); leapfrog / Strang (leading O(h^2) shadow terms); fourth-order symmetric splitting; processed or modified-kernel scheme (changes the known shadow term).
- Controls: h over at least one decade; fixed output cadence, ensemble, total time and float format; precision control (float64 vs double-double / quad) to separate round-off from discretization.
- Pre-registered observables: D_J(T) = <[J(T) - J(0)]^2> / (2T) on an ensemble away from absorbing boundaries; resonance-centre shift (FMA); stochastic-layer area or section-based proxy; finite-time LE and MEGNO; escape / crossing probability through an action boundary; modified-energy deviation and the leading analytic delta H_h.
- Test: D(h) = D_ref + a h^p (generic) versus D(h) = D_ref + F[delta H_h; local resonance geometry] (mechanistic), judged on held-out step sizes and ensembles.
- Falsified if: delta H_h predicts neither resonance displacement nor transport bias better than generic scaling; no detectable excess over the feasible h range; splittings with different delta H_h give the same diffusion once order and h are controlled; results flip with diffusion window, ensemble boundary or reference solver.

## Research line R3: scalable LE / CLV / indicator computation for high-dimensional symplectic maps

Gap (user-supplied scoping review, to verify): no standard benchmark battery or maintained open library targets LE spectra, CLVs and chaos indicators on symplectic maps with N >= 4 at scale; Geist1990-style comparisons predate modern hardware and cover flows, not maps. Contributions, folded into P3 plus a software paper (e.g. JOSS) for `chaoslab`:
- Benchmark suite: standard, Froeschle, coupled standard maps (N = 2..100), fixed panels, published runtimes and accuracies.
- Block / sparse QR for banded tangent maps; adaptive reorthonormalization schedules; partial spectra.
- Automatic-differentiation tangent dynamics (validated against analytic Jacobians).
- GPU and mixed-precision backends behind the same interface (only after the CPU reference is verified).
- Later: active-learning exploration of parameter space focused on regular/chaotic boundaries.

Ranked alternatives if R2 fails or as follow-ups: observable-specific convergence certification for instability hazards; floating-point vs discretization crossover in chaotic transport; reversible / adaptive schemes and resonance-crossing fidelity; shadow-Hamiltonian-informed step selection from local resonance widths; calibrated uncertainty bounds for FMA transport diagnostics.
