# Implementation notes: coupled-maps (catalog Y3, Y4)

Code (new): `src/Systems/froeschle_map.jl` (126 lines), `src/Systems/standard_map_chain.jl` (166), `src/Systems/wrapped_standard_map_chain.jl` (94), `bench/coupled-maps.jl` (36). Edited: `src/Systems/Systems.jl`, `src/Systems/residuals.jl`, `src/Systems/CONTEXT.md`. Spec `coupled-maps-spec.md` with Amendment D104, test plan `coupled-maps-test-plan.md`. Tests untouched (`git diff --stat -- test` empty). No new dependency. Nothing committed.

## Status

| Run | Result |
|-----|--------|
| RED (stage 02) | 0 pass, 5 fail, 69 error (all `UndefVarError`) |
| New suite alone, `-t2`, testenv of stage 02 | 556 pass, 0 fail, 0 error, first run after the code was written, no iteration; 20.6 s of evaluation |
| Whole suite `julia --project -t 4 -e "using Pkg; Pkg.test()"`, Julia 1.13.2 | 1 405 464 pass, 0 fail, 0 error, 1 broken, 25 min 22.7 s |

Suites in the whole run: quality 12 (Aqua and JET), layering 2, analysis 2763, coupled-maps 556, diagnostics 245, experiments 9267, fli-ofli 263533, lyap-benettin 15009, lyap-gram-schmidt 18938, lyap-householder-qr 575, lyap-svd-stabilized 12915, lyap-two-trajectory 493, megno 1000962, naff-fma 900, precision 150, rli 64460, sali 9210, standard-map 474, stretching 3476, tangent-core 808, weighted-birkhoff 716 + 1 broken. The one broken is the skipped device testset AC13 of weighted-birkhoff (known since its stage 03 notes). No "Bad Image" error occurred.

## Change table

| File | Lines | Change |
|------|-------|--------|
| `src/Systems/froeschle_map.jl` | +126 new | `FroeschleMap`, `WrappedFroeschleMap`, `torus_state`, `lift` (two-argument) |
| `src/Systems/standard_map_chain.jl` | +166 new | `StandardMapChain`, in-place forms, banded frame actions |
| `src/Systems/wrapped_standard_map_chain.jl` | +94 new | `WrappedStandardMapChain`, `torus_state`, `lift` |
| `src/Systems/Systems.jl` | +33 -1 | exports (4 types, 6 in-place names), docstrings of the six new generic functions, three `include`s, module docstring |
| `src/Systems/residuals.jl` | +13 -0 | `symplectic_residual(::AbstractMatrix)` dense method |
| `src/Systems/CONTEXT.md` | +14 -4 | files, exports, invariants, tests, deferred row (72 lines) |
| `bench/coupled-maps.jl` | +36 new | forward, tangent, frame action for Y3 and Y4 at N = 2, 10, 100 |

## Design

- Y3 (`FroeschleMap`): state `SVector{4}`, operation order of the stub and of `StandardMap` (`z3 - nu * sin(x1) - f`), so at `mu = 0` a step is the Y1 step on the shear image with `K = -nu` to the bit (criterion 9 passes with `==`, no FMA contraction: Julia does not fuse without `@fastmath`). Tangent and inverse tangent as `SMatrix{4,4}` from the shared `_abc`. `WrappedFroeschleMap` state `SVector{8}`; counts per spec Algorithm item 3; all eight components NaN when a count reaches `maxintfloat`. `propagate` works with it unchanged (state longer than `dimension`).
- Y4 (`StandardMapChain`): fields `N::Int`, `K::T`, `eps::T`, `ring::Bool`; the constructor promotes `K` and `eps` and throws `ArgumentError` for `N` outside 2:100, an unknown boundary, a ring with `N < 3` (reading R5). State any `AbstractVector` (decision P8: plain `Vector`, no StaticArrays for the chain). Neighbours of site `j` come from `_left`, `_right` (0 = none, ring wraps), so open chain and ring share one code path. `_force` is `K sin q_j + eps * (sum of neighbour sines)`; at `eps = 0` the sum is 0 and the result equals Y1 bitwise.
- Zero allocation: the inverse tangent and its action need `q = q' - p'`; `_Preimage(z, N)` is an isbits indexable view, so no vector is built. `_grow` returns a row of G as a tuple (no allocation). The dense forms fill a caller matrix; the pure forms allocate only the result.
- Frame actions are row by row with G's band: `M W = [W_q + G W_q + W_p; G W_q + W_p]`, `M^-1 W = [W_q - W_p; W_p + G (W_p - W_q)]`; cost `O(N k)`.
- `WrappedStandardMapChain` holds a `StandardMapChain` in field `chain`; same arithmetic as W-E2 and W-E4 per site; `tangent`, `inverse_tangent` read the first `2N` entries; a shared `_guard_counts!` poisons all `4N` entries.
- `torus_state(sys, z)` and `lift(sys, s)` two-argument methods for the wrapped Y3 and the wrapped chain; the one-argument Y1 methods are unchanged.
- `symplectic_residual(::AbstractMatrix)` builds a dense `Omega`; the static `symplectic_form` at d = 200 would compile a 40 000 element tuple.

## Performance (`bench/coupled-maps.jl`, `verify` environment, `-t2`, BenchmarkTools medians, indicative, machine not idle)

| Entry | Time | Allocations |
|-------|------|-------------|
| Y3 `forward`, 1000 steps | 18.4 us (18 ns per step) | 0 |
| Y3 `forward`, `tangent`, `inverse`, 4 x 4 frame product (single call) | below the timer resolution (about 100 ns), 0 | 0 |
| Y3 wrapped `forward` | about 100 ns (resolution limit) | 0 |
| Y4 N = 2, 10 `forward!`, `tangent!`, frame actions | 0 to 200 ns (resolution limit) | 0 |
| Y4 N = 100 `forward!` / wrapped `forward!` | 0.7 us / 1.1 us | 0 |
| Y4 N = 100 `tangent!` (dense 200 x 200) | 15.2 us | 0 |
| Y4 N = 100 `tangent_action!` / `inverse_tangent_action!` on 4 vectors | 1.2 us / 1.2 us | 0 |
| Y4 N = 100 dense tangent times the same frame (`mul!`, tangent excluded) | 14.9 us | 0 |

The banded action is 12 times cheaper than the dense product at N = 100 and 100 times cheaper than building the dense tangent first. The N = 100 to N = 50 timing ratio of criterion 17 (V) is stage 04. No time budget for systems exists in `bench/` or `_config` (the spec states memory formulas only); none checked.

## Audit

| Check | Result |
|-------|--------|
| GREEN | new suite 556 of 556; whole suite 1 405 464 pass, 0 fail, 0 error, 1 broken (pre-existing) |
| Quality | Aqua and JET clean (quality 12 of 12, `Systems` is in the OWN_MODULES list); layering 2 of 2. Coverage by `--code-coverage=user` on the new suite, `.cov` files parsed and deleted: `froeschle_map.jl` 56 of 56, `standard_map_chain.jl` 110 of 110, `wrapped_standard_map_chain.jl` 52 of 52 executable lines (100 percent, 218 of 218); the new `symplectic_residual` method is covered (`residuals.jl` 21 of 21) |
| Size | longest function: `forward!` of the wrapped chain, `_f3_tangent`, `inverse_tangent_action!` (all 14 lines or less); files 126, 166, 94 lines |
| Tests untouched | `git diff --stat -- test` empty |
| Numerics | no `@fastmath`; no accumulated logs in this module (no estimator), so no compensated sum is needed |
| Em dashes | none in the new or edited files |

## Readings and deviations

- None in behaviour against the spec or the Assumed API. Names and signatures are the Assumed API of the test plan.
- Field names: `FroeschleMap.nu`, `kappa`, `mu`; `StandardMapChain.N`, `K`, `eps`, `ring` (extra, `Bool`); `WrappedStandardMapChain.chain`. The tests read `nu` and `N` only.
- No constructor that promotes mixed types for `FroeschleMap` (like `StandardMap`: build the parameters in the state's element type, use `exact`); the chain promotes `K` and `eps` as the Assumed API needs.
- The sum of neighbour sines for the ring runs left then right neighbour (the stub ran right then wrap at site 1); differs by rounding only, no test depends on the order.
- The wrapped chain with a `Dual` state is not supported (`maxintfloat`), as for `WrappedStandardMap`; the lifted chain accepts Dual and Complex (criteria 1 and 2 pass).

## Process deviations

- The first shell command that wrote the chain files failed to parse (a quote inside a heredoc), so the two files were written with the editor tool instead; no effect on content.
- The full suite was run once with `-t 4` as the instruction asks (earlier notes use `-t 2`); one run only, after the new suite was GREEN.
- Coverage was measured on the new suite alone with the stage 02 test environment, not inside `Pkg.test`.

## Known limits and not done here

- `StandardMapChain` states are plain vectors: `propagate`, the tangent-core stream and the GPU kernels need isbits `SVector` states, so the chain is not usable there yet (P8). Y3 is.
- Criteria 20 and 21, the N = 100 / N = 50 timing ratio, the GPU parts of 18 and the 1 versus 8 thread comparison are stage 04 (`verify/coupled-maps.jl`).
- Not edited (outside this stage): `src/precompile.jl`, `docs/src/api/*.md`, `CHANGELOG.md`, `_config/catalogs/systems.md` (catalog sources column and status, stage 05).

## Conflicts between spec, docs and code

- `_config/architecture.md` ("Module map") names `step`, `tangent`, `invariants` and StaticArrays for dimension up to 8; the code uses `forward`, `inverse`, `tangent`, `inverse_tangent`, `dimension` and residual functions, and the chain has dimension up to 200. Followed the code.
- `src/Systems/Systems.jl` docstrings said `forward` takes an `SVector`; the chain takes an `AbstractVector` (docstring of `forward` not changed, the chain docstring states it).
- `src/Systems/CONTEXT.md` says Systems imports nothing from Foundation; kept (the new files use only the Systems helpers `_sin`, `_cos`, `_pi`, `_reduce`).
- Spec criterion 12 off-block bound and criterion 18 thread wording: see Amendment D104 (applied by the tests).
- `Systems.jl` module docstring and the `forward` docstring still describe the Y1 `SVector` state; `docs/src/api/` was not checked for the new exports.

## Rework D105 (2026-11-05)

Fixed-size (SVector) Y4 states, lifted (2N) and wrapped (4N), and the banded `ChainTangent`, per spec Amendment D105. Tests untouched (`test/coupled-maps/cm_d105.jl` from 50fba26). Not committed.

### Results

| Check | Result |
|-------|--------|
| coupled-maps suite (`myred.jl`, 2 threads) | GREEN 630 of 630 (592 passing before, 16 fail + 22 error before; 556 existing + 74 D105: D105a 18, D105b 16, D105c 40) |
| Whole suite (`Pkg.test()`, `-t 4`, 25 min 28 s) | 1 405 539 pass, 0 fail, 0 error, 1 broken (the existing weighted-birkhoff one) |
| Quality | Aqua and JET clean (quality 12 of 12); layering 2 of 2 |
| Coverage (`--code-coverage=user`, coupled-maps suite) | `chain_tangent.jl` 51 of 51, `wrapped_standard_map_chain.jl` 91 of 91, `standard_map_chain.jl` 110 of 110 executable lines (100 percent); the `.cov` files this run created were deleted (432 older `.cov` files from 10-29 and 10-31 stay, git-ignored, not mine) |
| Size | functions at most 14 lines; files 96 (`chain_tangent.jl`), 153 (wrapped chain) |
| Numerics | no `@fastmath`; no dependency added; no em dashes |

### Change table

| File | Diff | Change |
|------|------|--------|
| `src/Systems/chain_tangent.jl` | new, 96 lines | `ChainTangent{T,N}` (SVectors `g`, `gl`, `gr`, flags `ring`, `inv`; `size` (2N, 2N); `getindex` with the dense expressions so `Matrix` is bitwise the dense tangent; `*` with an `SMatrix` via `_band`, `vcat`, `_rows` of pair.jl); `forward`, `inverse`, `tangent`, `inverse_tangent` for `SVector` of the lifted chain |
| `src/Systems/wrapped_standard_map_chain.jl` | +59, -0 | `forward`, `inverse`, `tangent`, `inverse_tangent`, `torus_state`, `lift` for `SVector`; helpers `_check_wrapped`, `_guard_counts`, `_reduced` |
| `src/Systems/Systems.jl` | +2, -1 | export `ChainTangent`, include `chain_tangent.jl` |
| `src/Systems/CONTEXT.md` | +3, -2 (74 lines) | files table, export row, deferred row |
| `bench/coupled-maps.jl` | +13 | SVector entries |

The Vector-state paths are not edited. The wrapped struct is unchanged (field `chain`, still isbits).

### Compile times (first call, one session, Float64)

| N | lifted forward | tangent * W | inverse_tangent * W | wrapped forward | wrapped tangent * W | wrapped inverse |
|---|---|---|---|---|---|---|
| 2 | 0.02 s | 0.10 s | 0.02 s | 0.02 s | 0.02 s | 0.03 s |
| 10 | 0.03 s | 0.32 s | 0.05 s | 0.07 s | 0.04 s | 0.04 s |
| 100 | 0.12 s | 3.08 s | 0.92 s | 0.76 s | 0.59 s | 0.60 s |

Each N compiles its own methods (the sizes are types). The D105a testset (N = 100 included, both types) took 13 to 15 s.

### Bench (BenchmarkTools, median; timer resolution on this Windows host is 100 ns, so N = 2 and N = 10 entries are at the floor), all 0 bytes

| Entry | N = 2 | N = 10 | N = 100 |
|-------|-------|--------|---------|
| vector `forward!` | <100 ns | 100 ns | 700 ns |
| SVector `forward` | <100 ns | 100 ns | 700 ns |
| vector wrapped `forward!` | 100 ns | 100 ns | 1100 ns |
| SVector wrapped `forward` | 100 ns | 100 ns | 1100 ns |
| `tangent_action!` (4 vectors) | <100 ns | 100 ns | 1200 ns |
| SVector `tangent(..) * W` (4 vectors, includes building the band) | <100 ns | 100 ns | 1400 ns |
| `inverse_tangent_action!` | <100 ns | 100 ns | 1200 ns |
| SVector `inverse_tangent(..) * W` | 100 ns | 100 ns | 1300 ns |
| dense `tangent!` / dense times frame | <100 / 100 ns | 200 / 200 ns | 15 600 / 15 000 ns |

The SVector forms cost the same as the in-place vector forms (within the timer floor and 15 percent at N = 100) and are 11 times cheaper than the dense tangent product at N = 100.

### Readings and deviations

- Bitwise equality is kept by the same expression order: the force and the band come from `_force` and the cosine expressions of `_grow` (written per kind: `K cos q_j`, `eps cos(q_j - q_left)` or `+0.0`, then `g0 + cl + cr`, `gl = -cl`, `gr = -cr`); `getindex` uses the dense `tangent!` / `inverse_tangent!` expressions (`G + e`, `-e`, `-G`), so `-0.0` appears where the dense matrix has it.
- Allocation: the first version (tuple concatenation with `...`, tuples of 5-tuples and of pairs) allocated for N >= 10 (inference limits for long tuples, up to 95 kB at N = 100). Final form: `vcat` of SVector blocks, three homogeneous `ntuple`s for the band, `ntuple(j -> R[j][1], Val(N))` instead of `map(first, R)`, and `::Type{T}` in `_chain_tangent`. Cost of the wrapped step: no duplicate reduction.
- A wrong state length throws `DimensionMismatch` (also `m.N` against the length), as `PairTangent` does for a frame.
- The `ChainTangent` type holds `inv` as a field, not a type parameter (one type for both operators, so one compile per N).
- The `StandardMapChain` and `Systems.jl` docstrings were not changed (`forward` still says SVector of Y1); `ChainTangent` has its docstring. `docs/src/api/` and `src/precompile.jl` not checked or edited.
- Line endings: the working-copy files of this repo are mixed; `Systems.jl` and the wrapped file were normalised to LF (the index form), so git diffs only the real lines; `chain_tangent.jl` is LF.

### Conflicts between spec, docs and code

- Spec criterion text for D105a says "banded frame action of criterion 17" while the product is compared with `tangent_action!` of the chain (RD6 of the test plan); no code conflict.
- `Systems.jl` generic docstrings (`forward`, `tangent`) still describe the Y1 `SVector` state; the chain now accepts both `Vector` and `SVector`, with the same results.
- `src/Systems/CONTEXT.md` header says "(Y3, Y4)" built by coupled-maps; unchanged.
- Carried from the first stage 03: `_config/architecture.md` Module map names `step`, `invariants`; code uses `forward`, `inverse`, `tangent`, `inverse_tangent`.
