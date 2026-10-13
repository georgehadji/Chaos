# Implementation notes: precision (Foundation part)

Code: `src/Foundation/Foundation.jl` (precision items: `precision_label`, `parse_precision`, `exact`, re-export of `Double64`). Bench: `bench/precision.jl`. Systems items of the precision suite use `StandardMap` from `src/Systems/` (another agent).

## Status

`test/precision/runtests.jl`: 127 pass, 0 fail, 0 error (30 s, `-t 2`), run after Systems compiled. AC6 realised ratio max 2.61 with 1 exceeding orbit (allowance 10); AC13 median n* = 16.

## Deviations from the spec

- Float32 is not supported (D24): `precision_label(Float32)` has no method (MethodError).
- `exact` takes `Union{Integer,AbstractFloat,Irrational}` through `T(x)`: BigFloat(pi) is correctly rounded at the precision in force, `Double64(pi)` is DoubleFloats' own constant (measured 0.02 eps_D from 512-bit pi). Rational: `T(num)/T(den)`. String: `parse(T, s)` (Double64 parse measured 0.06 eps_D for "0.1", 0.009 eps_D for "0.001"; BigFloat parse is correctly rounded at current precision).
- `parse_precision` throws `ArgumentError` on any string other than `Float64`, `Double64`, `BigFloat(<digits>)`.

## Performance (measured by hand, BenchmarkTools is not installed in the test env)

- `add` (Neumaier), 10^7 calls in a loop: Float64 7.2 ns, Double64 307 ns, 0 bytes allocated.
- `exact(Double64, pi)`: about 1.9 us per call (not a hot path).

## Known limits

- `exact(Double64, x::BigFloat)` / `exact(Double64, x::BigInt)` rely on the DoubleFloats constructor; not tested.
- `precision_label(BigFloat)` reads `precision(BigFloat)` at call time; callers must set the precision before building labels.

## Review fixes (stage 04)

- `parse_precision` accepts `AbstractString`, anchors with `\z` (no trailing newline accepted), parses bits with `tryparse` and throws `ArgumentError` for unparsable bits or bits < 64.
