# Implementation notes: ofli (catalog row I6; spec D118)

Code (new): `src/Indicators/ofli.jl` (163 lines), `bench/ofli.jl` (28). Edited: `src/Indicators/Indicators.jl` (exports `OFLI`, `ofli_initial`; include after `fli.jl`; `using LinearAlgebra: dot`), `src/Indicators/CONTEXT.md` (79 lines), `docs/src/api/indicators.md` (one sentence; the names appear through `@autodocs`). Nothing under `test/`, `_config/` or `data/` was touched. Tests: `test/ofli/`, plan `stages/02-tests/output/ofli-test-plan.md` (the API proposal is implemented as written, points U1 to U13 as the operative reading).

## Status

| Run | Result |
|-----|--------|
| New suite, first run on the code (temp driver, Julia 1.13.2, 1 thread) | the first run had 61 errors, all `UndefVarError: dot` (a missing `using LinearAlgebra: dot` in `Indicators`); after the one-line fix 345 pass, 0 fail, 0 error, 1 min 7 s |
| Aqua, JET (`test/quality.jl`), layering (`test/layering.jl`), `test/fli-ofli/`, `test/sali/`, `test/megno/`, `test/gali/`, one process | 1276546 pass, 0 fail, 0 error, 0 broken, 7 min 17 s |
| Coverage, `--code-coverage=@src/Indicators/ofli.jl`, new suite alone | 68 of 69 executable lines, 98.6 percent (the one zero is the header line of `reducer_error`, a multi-line definition whose body lines are covered); the `.cov` file was deleted |
| Whole package suite, `Pkg.test()` (main session) | 1426757 pass, 0 fail, 0 error, 3 broken (integrator-catalog 2, weighted-birkhoff 1, as before), 41 min 16 s |
| Not run | `verify/` scripts, `OF_FULL_RUN` sizes, device runs (stage 04), `--check-bounds=yes` |

## Design

- `OFLI{F}(flow, replace)`, constructor `OFLI(flow; replace = true)`. The state `OFLIState{T, k}` is isbits: an `SVector{k}` of `CompensatedSum`, the running maxima `g1` (first column) and `gB` (basis), the instantaneous value, the companion unprojected maxima `f1`, `fB` (initial value 0, the s = 0 term), `sin theta` at the two maximisers and its minimum, the count of minus infinity terms, the step counters and the maximisers, a failure code. About 13 numbers plus `k` compensated sums.
- `reducer_step`: one `vector_field(flow, x)` evaluation (candidate A, P1), `u = f / ||f||` (scaled norm, so scaling `f` changes nothing); per column `alpha = safe_norm(w)` and `beta = safe_norm(w - (u . w) u)`; all columns are checked before anything is updated, so a failure leaves sums and maxima unchanged. `g_i = total(S_i, before) + log(beta_i)` (minus infinity for `beta = 0`, counted), then `S_i = add(S_i, log(alpha_i))`: the operations of FLI and SALI, so the sums are bitwise equal to theirs (AC10a). Maxima by strict increase, ties keep the earlier maximiser. `replace = true` divides every entry by its column's `alpha` (the same division as FLI), `replace = false` returns `W` as it arrived.
- `ofli_initial(flow, x0, W0)` is the pure s = 0 term (U1): the reducer sees no `x0` or `W0` in `reducer_init`, so `ofli` and `ofli1` cover the handoffs only and the caller combines both, as the tests do.
- `check_frame(OFLI(flow), W0, x0; orbit)`: unit columns within `8 effective_eps` (kind "initial frame", column), and for a single column a sine against the field not above `8 effective_eps`; a basis keeps a column along the flow.
- Health: "vector field" (code 2) when `||f||` is not finite or not above `floatmin(T)`, checked before the columns; "deviation matrix" (code 1, column) as in FLI. `reducer_init` throws `ArgumentError` for `dimension(flow) != d` or `k` outside `1:d` (P11, mismatching flow). Provenance: `(indicator = "ofli", ofli_flow, ofli_replace)`.
- Layering: only Foundation and Systems are used (`vector_field`, `dimension`, `safe_norm`, `CompensatedSum`); no Integrators import, no new edge (`test/layering.jl` passes). `LinearAlgebra` (already a dependency) is imported for `dot` only.

## Deviations from the spec and the test plan

1. None of substance: every point of the API proposal and of U1 to U13 is implemented as written.
2. `reducer_init` also refuses `k > d` or `k < 1` with `ArgumentError` (the spec says `1 <= k <= d`, the plan names only the dimension guard).
3. `check_frame` for `k = 1` returns the unit-norm defect (as the other methods); the parallel test throws and returns nothing else. A start whose field is zero passes the frame check and is reported by the stream as "vector field" (nothing to project against).
4. The initial values of `sin1`, `sinB`, `sinmin` are `Inf` (no handoff seen yet); the first handoff sets `sinmin`.

No `@fastmath`. Accumulated logs: compensated (`add`, `total`) for every `S_i`; the maxima are plain comparisons of single terms.

## Calibration margins (the numbers the test plan asked to record at stage 03)

| Item | Realized | Bound |
|------|----------|-------|
| AC1 worst ratio (error / 16 eps max(1, abs g)), `p0 = (0.6, 0.8)`, columns 1 to 4 | Float64: 0.323 (tau 1), 0.539 (tau 10); Double64: 0.553 (tau 1), 0.046 (tau 10) | 1 |
| AC1 worst ratio, `p0 = (1/2, 0)`, columns 3, 4 | Float64: 0.750 (tau 1), 0.431 (tau 10); Double64: 0.509 (tau 1), 0.043 (tau 10) | 1 |
| AC3a identity (256 bits, 100 sample times) | 2.6e-77 | 1e-30 |
| AC4a identity, relative, amplitude 2 / 3 (Double64 order-10 reference) | 1.2e-21 / 1.1e-19 | 1e-18 |
| AC4b spread of ln min ||f|| over 10 periods, amplitude 2 / 3 | 0 / 0 (the golden-section minimum lands on the same value each period) | 1e-12 |
| AC4c error, order-4 table, h = 1/20, t <= 10 P, amplitude 2 / 3 | 0 / 0: the maximum of the series is the s = 0 term (the start is a turning point, where the field norm is smallest and `g` largest), so the series is constant and no handoff exceeds it | 1e-5 |
| AC17 order-8 error at t = 10, 30, 60, orbit R | 1.0e-14, 8.2e-15, 1.7e-14 | 1e-6 |
| AC17 order-8 error at t = 10, 30, 60, orbit C | 7.3e-15, 3.5e-14, 6.8e-13 | 1e-6 |
| AC17 order-4 error (reported), R / C | 1.3e-5 to 1.7e-5 / 5.2e-6, 7.2e-6, 5.8e-4 | reported |

The AC4c result means the test does not exercise a departure from the plateau at the chosen amplitudes: it passes at 0 error because the plateau is the s = 0 value on both sides. AC17 order 8 is at rounding level (1e-14 to 7e-13), several orders under the bound; the stage 04 Taylor form is the real discretization test.

## Performance (`bench/ofli.jl`, benchmark environment, 1 thread, medians; indicative, BenchmarkTools resolution 100 ns)

| Entry | Time | Allocations |
|-------|------|-------------|
| one reducer step, k = 4, Float64, replacing / reader | 0.20 / 0.20 us | 0 |
| one reducer step, k = 4, Double64 | 1.4 us | 0 |
| `ofli_initial`, k = 4 | 0.1 us | 0 |
| `check_frame` (k = 4) | constant folded by the benchmark (below resolution) | 0 |
| orbit C, 1000 time units (2 10^4 steps, KahanLi8, tau = 20), k = 4 | 14.9 ms | 0 |
| same, chain `(OFLI(replace = false), SALI())`, k = 2 | 15.3 ms | 0 |
| 1024 orbits, 200 time units, k = 4, CPU | 3.14 s | 3128 (the result vectors) |

The 1000 handoffs add about 0.2 ms to a run whose steps cost about 14.7 ms (the order-8 tangent step is 1.3 us in henon-heiles): the reducer is about 1.4 percent of the stream.

## Known limits

- The reducer holds the flow and cannot check it is the one the map integrates (spec P10): a mismatch gives wrong numbers silently; the run record carries the flow name for the experiment layer to compare with the map's.
- Candidate A: on a periodic orbit the plateau is lost after a time of order `h^-p` (O11, P6; measured at stage 04, R1). `ofli` and `ofli1` range over the handoffs only; the s = 0 term is `ofli_initial`.
- A state vector longer than the flow dimension (winding counts) is not handled: `vector_field` receives `x` as is; the flows of this layer have `s = d`.
- The companions `fli`, `fli1` are reported for any `k`; they equal the delivered FLI only at `k = 2` (AC10a).
- GPU: the state is isbits for Float64 and Double64 and the step allocates nothing, but the device kernel with the flow held in the reducer is not run here (stage 04, `verify/gpu.jl`); BigFloat is CPU only.
- `check_frame` for `k = 1` needs the field at `x0`; a zero field at `x0` is left to the stream's health error.

## Audit

| Check | Result |
|-------|--------|
| GREEN | new suite 345 pass, 0 fail; quality, layering and the fli-ofli, sali, megno, gali suites 1276546 pass, 0 fail; whole package suite 1426757 pass, 0 fail, 3 broken |
| Quality | Aqua and JET clean, layering passes (the combined run above); coverage of `ofli.jl` 68 of 69 executable lines (98.6 percent) by `--code-coverage=@<file>`, `.cov` deleted |
| Size | longest functions: `ofli_advance` 20 lines, `reducer_step` 14, `check_frame` 17; `ofli.jl` 163 lines; `src/Indicators/CONTEXT.md` 79 lines |
| Tests untouched | `git diff -- test/` empty |
| Numerics | no `@fastmath`; every accumulated log is a `CompensatedSum` |
| Em dashes | none in the new or edited files |
