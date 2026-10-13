# Diagnostics

Residuals every estimator reports, and the reducer that checks them inside the tangent stream. Imports Foundation, Systems and Tangent; never a Lyapunov module, so exponents arrive as data. Built by the diagnostics run (catalog section F).

Spec: `workspaces/module-build/stages/01-spec/output/diagnostics-spec.md` (DG-; B1 to B4 are its floating-point bounds, I1 to I5 its identities; Amendment D30). Reports: `stages/04-verify/output/diagnostics-verification.md`.

## Files

| File | Holds |
|------|-------|
| `Diagnostics.jl` | module and exports |
| `residuals.jl` | the five residual functions, the three thresholds, `check_residual` |
| `monitor.jl` | `ResidualMonitor` reducer, its isbits state, `monitor_exponents` |

## Exports

| Name | What it does |
|------|--------------|
| `eps_det(M)`, `eps_symp(M)` | `abs(det M - 1)`; max entry of `M' Omega M - Omega` |
| `eps_orth(Q)` | max entry of `Q'Q - I` |
| `eps_sum(lambda)`, `eps_pair(lambda)` | `abs(sum lambda)`; max over `i` of `abs(lambda_i + lambda_{2N+1-i})`, estimator order |
| `det_threshold`, `symp_threshold`, `sum_threshold` | thresholds as multipliers of `eps(T)` times a growth model (keys in `_config/numerics-defaults.toml`) |
| `check_residual` | throw `NumericalHealthError` for a non-finite value or one above its threshold |
| `ResidualMonitor(sys, det_tol, orth_tol)` | reducer: det, symplecticity and orthogonality at each handoff, running maxima in the state |
| `monitor_exponents` | feed exponents and window growth `G`: checks the sum, records the pairing residual |

## Invariants a caller can rely on

| Invariant | Guaranteed by |
|-----------|---------------|
| Exact integer maps give residuals exactly 0 | DG-AC1 |
| Det and symplecticity residuals stay under their thresholds on symplectic maps | DG B1, B3, DG-AC2, DG-AC6 |
| The sum residual has no part growing with `n` (compensated accumulation required of the estimator) | DG I2, B4, DG-AC3, DG-AC7 |
| For `2N = 2` pairing equals the sum residual; pairing indices are `i` and `2N+1-i`, never adjacent | DG I3, DG-AC4 |
| A fault just above a threshold raises with the right name, value, orbit and step; just below, nothing | DG-AC5 |
| Names are reported in the fixed order det, symplecticity, orthogonality, sum, pairing; a non-finite value always raises | DG Algorithm "Health error", DG-AC5 |
| Over an ensemble the lowest failing orbit id is reported, whatever the threads or backend | DG-AC10, DG-AC11 |
| Per-orbit results are bitwise equal for 1 and 4 threads | DG-AC10 |
| `reducer_update` and `monitor_exponents` allocate nothing after warm-up | DG-AC9 |
| The module never raises "orthonormalisation breakdown"; the orthonormalisers do | DG Amendment D30 |

## Tests and verification

- Package: `test/diagnostics/` (`runtests.jl`, `faultmaps.jl` with exact-fault maps).
- Oracle run: `verify/diagnostics.jl` (256-bit recomputation at the stored points, 10^6-step panel, thread child process); figures and log in `stages/04-verify/output/`.
- GPU line (monitor on the device, Float64): `verify/gpu.jl`, section GPU of the verification report.
- Notes: `stages/03-implement/output/diagnostics-impl-notes.md`; bench `bench/diagnostics.jl`.

## Deferred

| Item | Where |
|------|-------|
| Pairing health check for `N >= 2` at finalise (needs the allowance kappa, no default before the first `N >= 2` system) | DG open point O2, M3 |
| ChaosTools full-spectrum sum next to `eps_sum` (reported, not pass/fail) | lyap-benettin stage 04, D27 |
| Tighten or confirm `pairing_residual_tol` from the lyap-* maxima (measured far below the key) | DG open point O3, D30 |
| Energy and modified-energy errors; SALI versus `lambda_1 - lambda_2` | M5; M4 |
| `M_max` is not in `reducer_final` (the kernel-path check uses `det_max / eps`) | `diagnostics-verification.md`, review findings |
| Double64 monitor on the GPU (kernel does not compile; run falls back with reason 3) | tangent-core GPU section |
