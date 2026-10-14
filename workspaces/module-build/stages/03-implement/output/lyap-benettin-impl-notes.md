# Implementation notes: lyap-benettin

Code: `src/Lyapunov/` (`Lyapunov.jl`, `norm.jl` with `safe_norm`, `benettin.jl` with `Benettin` and `BenettinState`; 4 + 11 + 56 lines, no function over 50 lines), included and re-exported from `src/Chaoslab.jl` after Diagnostics. Imports Foundation, Systems and Tangent only (layering test passes). Bench: `bench/lyap-benettin.jl` (5 entries, not run: BenchmarkTools is not a dependency). The other Batch B names (QRSpectrum, orthonormalizers, SVDSpectrum, two-trajectory, `orth_loss`, `window_G` and so on) do not exist yet; their suites stay RED with UndefVarError only (checked in the full-suite logs: 120 errors, all UndefVarError).

## Status

`test/lyap-benettin/`: 15009 pass, 0 fail, 0 error, Julia 1.13.1, 2 threads. Whole package suite (`Pkg.test`, 2 threads): Aqua and JET clean (quality 12/12), layering 2/2, diagnostics 245, precision 150, standard-map 774, tangent-core 565, lyap-benettin 15009 on both Julia 1.13.1 and 1.10; the only non-pass results are the RED suites of the four later Batch B slugs (101 fail, 120 error, identical on both versions). Standard-map took 100 min on Julia 1.10 (9.5 min on 1.13), a Batch A runtime issue not touched here.

## Runtime against the 3 minute ceiling (D35)

Testset time of `lyap-benettin` on 2 threads, compile of the test code included (package already precompiled): 2 min 44 s (first run), 2 min 59.6 s (idle re-run, 3 min 53 s wall including process start and package load), 3 min 22 s inside the full suite on 1.13 and on 1.10 (machine shared with other runs). So the suite is at the ceiling and over it under load. The test plan lists no ordered cut list beyond the sample reductions already built in (its deviation 3), and tests are not edited here, so no cut was applied. Largest testsets: AC5 29 s, AC1 21 to 24 s, AC19 18 s, AC17 16 s, AC10 15 s, AC2 11 s. Owner decision needed: cut AC5 or AC1 sample sizes, or measure the ceiling on the idle machine.

## Deviations from the spec and contract

- `safe_norm` scales by the maximum absolute entry with a division (not by a power of two), as the contract says; meets the 2 eps bound of AC8.
- Zero or non-finite norm is recorded with code 1; the error has `value` = the norm, `threshold` 0.0 (Inf when the norm is non-finite), `column` 1, `what` "deviation matrix". A state that is already failed ignores later steps.
- `reducer_final` at n = 0 returns zeros (lambda 0); otherwise `lambda = total(sum) / n` with n the last successful handoff.
- `BenettinState{T}` fields: `sum`, `n`, `code`, `value`, `step` (contract).

## Performance (Julia 1.13.1, standard map K = 2, 1 thread, measured at top level, indicative only)

- `reducer_step`: about 0.2 us Float64, 1.5 us Double64.
- `propagate_orbit` 10^4 steps, chain (Benettin, ResidualMonitor): 0 bytes allocated in both types (AC10c). Time per step is dominated by `sin`/`cos` of the map on a chaotic lifted orbit: 0.5 us (Benettin alone) and 0.8 us (with monitor) Float64; 16 and 34 us Double64.

## Known limits

- GPU: the step is isbits and uses `log`, `sqrt`, division only; Double64 `log` compiling on the device is untested here (stage 04; a failure falls back per chunk with reason 3).
- Coverage not measured in this run (every line of `src/Lyapunov/` is exercised by the suite).
- Test-side `lb_Probe` states hold Vectors, so those tests use `propagate_orbit`, not `propagate`.
