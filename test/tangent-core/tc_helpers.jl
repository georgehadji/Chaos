# Helpers for the tangent-core suite (prefix tc_). Test maps and reducers are defined here, never in src/.
using Test, Chaoslab, StaticArrays, LinearAlgebra, StableRNGs, DoubleFloats

const TC_SEED = 20261012   # base seed, _config/numerics-defaults.toml; per-criterion offsets are added in the testsets

const tc_AbstractMap = isdefined(Chaoslab, :AbstractMap) ? Chaoslab.AbstractMap : Any

# Linear map z -> A z (Y5 maps); tangent is A everywhere.
struct tc_Lin{T,D,L} <: tc_AbstractMap
    A::SMatrix{D,D,T,L}
end
# Standard map p' = p + K sin x, x' = x + p'; tangent from spec E8.
struct tc_Std{T} <: tc_AbstractMap
    K::T
end
# Fault map on state (step counter, value, flag); only flagged orbits misbehave.
struct tc_Fault <: tc_AbstractMap
    kind::Symbol   # :inf_state (state becomes Inf when the counter reaches `at`) or :nan_W (tangent has a NaN at counter `at`)
    at::Int
end
# Map whose forward throws: proves that argument validation happens before any work.
struct tc_Boom <: tc_AbstractMap end

@static if isdefined(Chaoslab, :forward) && isdefined(Chaoslab, :tangent) && isdefined(Chaoslab, :dimension)
    Chaoslab.dimension(::tc_Lin{T,D}) where {T,D} = D
    Chaoslab.forward(m::tc_Lin, z) = m.A * z
    Chaoslab.tangent(m::tc_Lin, z) = m.A

    Chaoslab.dimension(::tc_Std) = 2
    function Chaoslab.forward(m::tc_Std, z)
        p2 = z[2] + m.K * sin(z[1])
        return SVector(z[1] + p2, p2)
    end
    function Chaoslab.tangent(m::tc_Std, z)
        c = m.K * cos(z[1])
        return SMatrix{2,2}(1 + c, c, one(c), one(c))
    end

    Chaoslab.dimension(::tc_Fault) = 3
    function Chaoslab.forward(m::tc_Fault, z)
        n = z[1] + 1
        y = (m.kind === :inf_state && z[3] == 1 && n == m.at) ? Inf : z[2] + 1
        return SVector(n, y, z[3])
    end
    function Chaoslab.tangent(m::tc_Fault, z)
        a = (m.kind === :nan_W && z[3] == 1 && z[1] == m.at) ? NaN : 1.0
        return SMatrix{3,3,Float64}(a, 0, 0, 0, 1, 0, 0, 0, 1)
    end

    Chaoslab.dimension(::tc_Boom) = 2
    Chaoslab.forward(::tc_Boom, z) = error("tc_Boom: forward was called")
    Chaoslab.tangent(::tc_Boom, z) = error("tc_Boom: tangent was called")
end

# Column-normalizing reducer: replaces W by W with unit columns and accumulates log column norms (CompensatedSum).
struct tc_LogNorm end
# Read-only reducer: counts handoffs, checks I2 (n == previous + tau), records max |W| seen.
struct tc_Reader
    tau::Int
end

@static if isdefined(Chaoslab, :reducer_init) && isdefined(Chaoslab, :reducer_update) && isdefined(Chaoslab, :reducer_final)
    function Chaoslab.reducer_init(::tc_LogNorm, ::Type{T}, ::Val{d}, ::Val{k}) where {T,d,k}
        return ntuple(_ -> CompensatedSum{T}(), Val(k))
    end
    function Chaoslab.reducer_update(::tc_LogNorm, st, n::Int, x, W::SMatrix{d,k,T}) where {d,k,T}
        s = ntuple(j -> sqrt(sum(abs2, ntuple(i -> W[(j - 1) * d + i], Val(d)))), Val(k))
        st2 = ntuple(j -> add(st[j], log(s[j])), Val(k))
        Wn = SMatrix{d,k,T}(ntuple(i -> W[i] / s[(i - 1) ÷ d + 1], Val(d * k)))
        return st2, Wn
    end
    Chaoslab.reducer_final(::tc_LogNorm, st) = (logs = map(total, st),)

    Chaoslab.reducer_init(::tc_Reader, ::Type{T}, ::Val, ::Val) where {T} = (0, 0, 0, zero(T))
    function Chaoslab.reducer_update(r::tc_Reader, st, n::Int, x, W)
        c, last, bad, m = st
        return (c + 1, n, bad + (n == last + r.tau ? 0 : 1), max(m, maximum(abs, W))), W
    end
    Chaoslab.reducer_final(::tc_Reader, st) = (count = st[1], last = st[2], bad = st[3], maxabs = st[4])
end

# Fibonacci numbers, F_0 = 0, F_1 = 1 (Int64 is exact up to F_92).
function tc_fib(n::Int)
    a, b = 0, 1
    for _ in 1:n
        a, b = b, a + b
    end
    return a
end
# A^n for A = [[2,1],[1,1]] (E9): [[F_{2n+1}, F_{2n}], [F_{2n}, F_{2n-1}]], as Float64 (exact below 2^53).
tc_Apow(n::Int) = SMatrix{2,2,Float64}(tc_fib(2n + 1), tc_fib(2n), tc_fib(2n), tc_fib(2n - 1))

# Exact conversion to BigFloat (Float64 and Double64).
tc_big(x::Double64) = BigFloat(x.hi) + BigFloat(x.lo)
tc_big(x::Real) = BigFloat(x)
tc_bigmat(M) = [tc_big(M[i, j]) for i in 1:size(M, 1), j in 1:size(M, 2)]

# Matrix power by repeated squaring at the current BigFloat precision.
function tc_matpow(A::AbstractMatrix, n::Int)
    R = Matrix{eltype(A)}(I, size(A)...)
    B = copy(A)
    while n > 0
        isodd(n) && (R = R * B)
        B = B * B
        n >>= 1
    end
    return R
end

# Kolmogorov-Smirnov statistic of a sample against the uniform distribution on [0, 1].
function tc_ks_uniform(u::AbstractVector)
    s = sort(u)
    n = length(s)
    return maximum(max(i / n - s[i], s[i] - (i - 1) / n) for i in 1:n)
end

# Run f(); return the NumericalHealthError it throws (nothing if it returns), rethrow anything else.
function tc_health_error(f)
    try
        f()
        return nothing
    catch e
        e isa NumericalHealthError ? e : rethrow()
    end
end

# Compensated and plain sums of n copies of x.
function tc_csum(n::Int, x::T) where {T}
    acc = CompensatedSum{T}()
    for _ in 1:n
        acc = add(acc, x)
    end
    return total(acc)
end
function tc_plainsum(n::Int, x::T) where {T}
    s = zero(T)
    for _ in 1:n
        s += x
    end
    return s
end

# Bytes allocated by one propagate_orbit call after a warm-up call (measured inside a function).
function tc_alloc_bytes(sys, x0, W0, N, tau, reducers)
    propagate_orbit(sys, x0, W0, N, tau, reducers)
    return @allocated propagate_orbit(sys, x0, W0, N, tau, reducers)
end

# AC4 residual: ||W_n - W_ref||_F and the bound 8 n eps(T) (prod ||M_i||_F) ||W_0||_F, both as Float64.
# W_ref = M_{n-1} ... M_0 W_0 at BigFloat(256) with M_i from E8 at the realized orbit points pts = [x_0, ..., x_{n-1}].
function tc_std_residual(Kt, pts, W0, Wn, epsT)
    return setprecision(BigFloat, 256) do
        Kb = tc_big(Kt)
        Wref = tc_bigmat(W0)
        prodM = one(BigFloat)
        for z in pts
            c = Kb * cos(tc_big(z[1]))
            M = [1+c one(c); c one(c)]
            prodM *= sqrt(sum(abs2, M))
            Wref = M * Wref
        end
        err = sqrt(sum(abs2, tc_bigmat(Wn) - Wref))
        bound = 8 * length(pts) * tc_big(epsT) * prodM * sqrt(sum(abs2, tc_bigmat(W0)))
        (Float64(err), Float64(bound))
    end
end

# Device plumbing shared by the AC10 and D29 testsets: the run-wide device is set and restored explicitly, so no test
# depends on whether CUDA happens to be loaded.
const tc_F = Chaoslab.Foundation
const tc_Tg = Chaoslab.Tangent
function tc_with_device(f, dev)
    old = tc_F.DEVICE[]
    tc_F.DEVICE[] = dev
    try
        return f()
    finally
        tc_F.DEVICE[] = old
    end
end
