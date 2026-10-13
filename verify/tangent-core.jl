# Stage 04 verification script for tangent-core (spec: workspaces/module-build/stages/01-spec/output/tangent-core-spec.md,
# including Amendment D28). Run from the repo root:  julia --project=verify -t 4 verify/tangent-core.jl
# One line per criterion: AC<n> | quantity | value | bound | verdict.  Child mode (-t 1 runs) is started by this script itself.
using Chaoslab, StaticArrays, LinearAlgebra, StableRNGs, DoubleFloats, Printf, Serialization

const CHILD = length(ARGS) >= 2 && ARGS[1] == "child"
CHILD || @eval using CairoMakie

const SEED = 20261012
const ROOT = normpath(joinpath(@__DIR__, ".."))
const OUT = joinpath(ROOT, "workspaces", "module-build", "stages", "04-verify", "output")
const ROWS = Tuple{String,String}[]

fmt(x::AbstractFloat) = @sprintf("%.4g", x)
fmt(x) = string(x)
function row(id, qty, val, bound, verdict)
    push!(ROWS, (id, verdict))
    println(id, " | ", qty, " | ", fmt(val), " | ", fmt(bound), " | ", verdict)
    flush(stdout)
end
pf(b::Bool) = b ? "PASS" : "FAIL"
epsf(::Type{T}) where {T} = Float64(eps(T))

# ---------------------------------------------------------------------------------------------- test maps and reducers
struct Fault <: Chaoslab.AbstractMap
    kind::Symbol   # :inf_state (state Inf when the counter reaches `at`) or :nan_W (tangent has a NaN at counter `at`)
    at::Int
end
Chaoslab.dimension(::Fault) = 3
function Chaoslab.forward(m::Fault, z)
    n = z[1] + 1
    y = (m.kind === :inf_state && z[3] == 1 && n == m.at) ? Inf : z[2] + 1
    return SVector(n, y, z[3])
end
function Chaoslab.tangent(m::Fault, z)
    a = (m.kind === :nan_W && z[3] == 1 && z[1] == m.at) ? NaN : 1.0
    return SMatrix{3,3,Float64}(a, 0, 0, 0, 1, 0, 0, 0, 1)
end

struct Boom <: Chaoslab.AbstractMap end   # forward throws: argument validation must happen before any work
Chaoslab.dimension(::Boom) = 2
Chaoslab.forward(::Boom, z) = error("Boom: forward was called")
Chaoslab.tangent(::Boom, z) = error("Boom: tangent was called")

struct LogNorm end   # replaces W by unit columns, accumulates log column norms (compensated)
struct Reader; tau::Int; end   # read-only: handoff count, I2 check, max |W|
struct LenReader end   # records the length of the point it is handed (W9e)

function Chaoslab.reducer_init(::LogNorm, ::Type{T}, ::Val{d}, ::Val{k}) where {T,d,k}
    return ntuple(_ -> CompensatedSum{T}(), Val(k))
end
function Chaoslab.reducer_update(::LogNorm, st, n::Int, x, W::SMatrix{d,k,T}) where {d,k,T}
    s = ntuple(j -> sqrt(sum(abs2, ntuple(i -> W[(j - 1) * d + i], Val(d)))), Val(k))
    st2 = ntuple(j -> add(st[j], log(s[j])), Val(k))
    Wn = SMatrix{d,k,T}(ntuple(i -> W[i] / s[(i - 1) ÷ d + 1], Val(d * k)))
    return st2, Wn
end
Chaoslab.reducer_final(::LogNorm, st) = (logs = map(total, st),)

Chaoslab.reducer_init(::Reader, ::Type{T}, ::Val, ::Val) where {T} = (0, 0, 0, zero(T))
function Chaoslab.reducer_update(r::Reader, st, n::Int, x, W)
    c, last, bad, m = st
    return (c + 1, n, bad + (n == last + r.tau ? 0 : 1), max(m, maximum(abs, W))), W
end
Chaoslab.reducer_final(::Reader, st) = (count = st[1], last = st[2], bad = st[3], maxabs = st[4])

Chaoslab.reducer_init(::LenReader, ::Type, ::Val, ::Val) = 0
Chaoslab.reducer_update(::LenReader, st, n::Int, x, W) = (length(x), W)
Chaoslab.reducer_final(::LenReader, st) = (len = st,)

# state (x, p, counter, running sum of x), dimension 2; `fault` 3 or 4 sends that component to Inf when the counter reaches `at`
    struct V9 <: Chaoslab.AbstractMap
        K::Float64
        fault::Int
        at::Int
    end
    struct V9Boom <: Chaoslab.AbstractMap end
    Chaoslab.dimension(::V9) = 2
    Chaoslab.dimension(::V9Boom) = 2
    function Chaoslab.forward(m::V9, z)
        p = z[2] + m.K * sin(z[1])
        c = z[3] + 1
        s = z[4] + z[1]
        c = (m.fault == 3 && c == m.at) ? Inf : c
        s = (m.fault == 4 && c == m.at) ? Inf : s
        return SVector(z[1] + p, p, c, s)
    end
    function Chaoslab.tangent(m::V9, z)
        c = m.K * cos(z[1])
        return SMatrix{2,2}(1 + c, c, one(c), one(c))
    end
    Chaoslab.forward(::V9Boom, z) = error("V9Boom: forward was called")
    Chaoslab.tangent(::V9Boom, z) = error("V9Boom: tangent was called")

# ---------------------------------------------------------------------------------------------------------- helpers
bigf(x::Double64) = BigFloat(x.hi) + BigFloat(x.lo)
bigf(x::Real) = BigFloat(x)
bigmat(M) = [bigf(M[i, j]) for i in 1:size(M, 1), j in 1:size(M, 2)]
fro(M) = sqrt(sum(abs2, M))

function matpow(A::AbstractMatrix, n::Int)
    R = Matrix{eltype(A)}(I, size(A)...)
    B = copy(A)
    while n > 0
        isodd(n) && (R = R * B)
        B = B * B
        n >>= 1
    end
    return R
end

function ks_uniform(u)
    s = sort(u)
    n = length(s)
    return maximum(max(i / n - s[i], s[i] - (i - 1) / n) for i in 1:n)
end

function health_error(f)
    try
        f()
        return nothing
    catch e
        e isa NumericalHealthError ? e : rethrow()
    end
end

function allocbytes(sys, x0, W0, N, tau, red)
    propagate_orbit(sys, x0, W0, N, tau, red)   # warm-up
    return @allocated propagate_orbit(sys, x0, W0, N, tau, red)
end

function csum(n::Int, x::T) where {T}
    acc = CompensatedSum{T}()
    for _ in 1:n
        acc = add(acc, x)
    end
    return total(acc)
end
function plainsum(n::Int, x::T) where {T}
    s = zero(T)
    for _ in 1:n
        s += x
    end
    return s
end

# 20 seeded initial points on the torus [0, 2pi)^2 in element type T (as test AC4)
function sample_points(T, seed, nsamp = 20)
    rng = StableRNG(seed)
    return [SVector{2,T}(T(2pi * rand(rng)), T(2pi * rand(rng))) for _ in 1:nsamp]
end

# ------------------------------------------------------------------------------------------------ child mode helpers
digest(xs, Ws, fs) = hash((hash(collect(reinterpret(UInt64, collect(reinterpret(Float64, xs))))),
                           hash(collect(reinterpret(UInt64, collect(reinterpret(Float64, Ws))))),
                           hash(collect(reinterpret(UInt64, Float64[l for f in fs for l in f[1].logs]))),
                           hash([f[2] for f in fs])))

function ac6_digests()
    seed = SEED + 6
    n_orb, N = 4096, 1000
    sys = StandardMap(2pi)
    rng = StableRNG(seed)
    x0 = [SVector(2pi * rand(rng), 2pi * rand(rng)) for _ in 1:n_orb]
    reducers = (LogNorm(), Reader(1))
    out = Dict{Any,UInt}()
    for chunk in (64, 1000)
        res = propagate(sys, x0, N; k = 2, tau = 1, seed = seed, reducers = reducers, chunk = chunk)
        out[chunk] = digest(res.x, res.W, res.finals)
    end
    W0 = [initial_deviations(seed, j, Val(2), Val(2), Float64) for j in 1:n_orb]
    serial = [propagate_orbit(sys, x0[j], W0[j], N, 1, reducers; orbit = j) for j in 1:n_orb]
    out[:serial] = digest([s[1] for s in serial], [s[2] for s in serial], [map(reducer_final, reducers, s[3]) for s in serial])
    return out
end

function ac13_time(T, norb, N)
    seed = SEED + 13
    sys = StandardMap(T(2pi))
    rng = StableRNG(seed)
    x0 = [SVector{2,T}(T(2pi * rand(rng)), T(2pi * rand(rng))) for _ in 1:norb]
    kw = (k = 2, tau = 1, seed = seed, reducers = (LogNorm(),))
    propagate(sys, x0[1:64], 10; kw...)   # warm-up (compilation)
    return @elapsed propagate(sys, x0, N; kw...)
end

# sample sizes for AC13 (CPU side; the spec's 10^6 orbits x 10^4 steps is too long on this machine)
const AC13_SIZES = Dict(Float64 => (10^5, 10^3), Double64 => (2000, 10^3))

if CHILD
    BLAS.set_num_threads(1)
    res = Dict{Any,Any}(:nthreads => Threads.nthreads(), :ac6 => ac6_digests())
    for (T, (no, N)) in AC13_SIZES
        res[(:ac13, T)] = ac13_time(T, no, N)
    end
    serialize(ARGS[2], res)
    exit(0)
end

# ------------------------------------------------------------------------------------------------------------- header
println("# tangent-core verification, julia ", VERSION, ", threads ", Threads.nthreads(), ", cpu ", Sys.CPU_NAME,
        " (", Sys.CPU_THREADS, " hardware threads), base seed ", SEED)
println("# format: AC<n> | quantity | value | bound | verdict")
BLAS.set_num_threads(Threads.nthreads() > 1 ? 1 : BLAS.get_num_threads())

# ------------------------------------------------------------------------------------------------------------------ AC1
let
    seed = SEED + 1
    sys = LinearSymplecticMap(SMatrix{2,2,Float64}(2, 1, 1, 1))
    x0 = SVector(3.0, 5.0)
    W0 = SMatrix{2,2,Float64}(1, 0, 0, 1)
    mism = 0
    nbig = 0.0
    for n in 1:35
        res = propagate(sys, [x0], n; k = 2, tau = 1, seed = seed, W0 = [W0])
        a = 0   # A^n = [[F_{2n+1}, F_{2n}], [F_{2n}, F_{2n-1}]] computed in Int64
        An = Int.(matpow([2 1; 1 1], n))
        Wexp = SMatrix{2,2,Float64}(An[1, 1], An[2, 1], An[1, 2], An[2, 2])
        xexp = SVector(Float64(An[1, 1] * 3 + An[1, 2] * 5), Float64(An[2, 1] * 3 + An[2, 2] * 5))
        mism += count(res.W[1] .!= Wexp) + count(res.x[1] .!= xexp)
        nbig = max(nbig, maximum(abs, Wexp))
    end
    row("AC1", "mismatching entries of W_n and x_n vs integer A^n, n=1..35 (largest entry $(Int(nbig)) < 2^53)", mism, 0, pf(mism == 0))
end

# ------------------------------------------------------------------------------------------------------------------ AC2
let
    seed = SEED + 2
    N = 10_000
    sys = LinearSymplecticMap(SMatrix{2,2,Float64}(2, 1, 1, 1))
    x0 = [SVector(0.0, 0.0) for _ in 1:100]
    W0 = [initial_deviations(seed, j, Val(2), Val(1), Float64) for j in 1:100]
    ref = setprecision(BigFloat, 256) do
        AN = matpow(BigFloat[2 1; 1 1], N)
        [Float64(log(sqrt(sum(abs2, AN * bigmat(w))))) for w in W0]
    end
    lam = log((3 + sqrt(5)) / 2)
    for tau in (1, 5, 10)
        res = propagate(sys, x0, N; k = 1, tau = tau, seed = seed, W0 = W0, reducers = (LogNorm(), Reader(tau)))
        wa = wb = 0.0
        for j in 1:100
            S = res.finals[j][1].logs[1]
            wa = max(wa, abs(S + log(norm(res.W[j])) - ref[j]))
            wb = max(wb, abs(S / N - lam))
        end
        ok_i2 = all(f -> f[2].count == N ÷ tau && f[2].last == N && f[2].bad == 0, res.finals)
        row("AC2a", "tau=$tau max over 100 orbits |S_n + ln||w_n|| - ln||A^n w0||_BigFloat(256)|", wa, 8 * N * epsf(Float64), pf(wa <= 8 * N * epsf(Float64)))
        row("AC2b", "tau=$tau max |S_n/n - ln((3+sqrt5)/2)| (exact value $(fmt(lam)))", wb, 2e-3, pf(wb <= 2e-3))
        row("AC2c", "tau=$tau handoff count N/tau and I2 spacing on all 100 orbits", ok_i2, true, pf(ok_i2))
    end
end

# ------------------------------------------------------------------------------------------------------------------ AC3
let
    seed = SEED + 3
    sys = StandardMap(0.0)
    x0 = SVector(0.25, 0.5)
    W0 = SMatrix{2,2,Float64}(1, 0, 0, 1)
    for n in (1, 10, 1000, 10^6)
        res = propagate(sys, [x0], n; k = 2, tau = 1, seed = seed, W0 = [W0])
        bad = count(res.W[1] .!= SMatrix{2,2,Float64}(1, 0, n, 1)) + count(res.x[1] .!= SVector(0.25 + 0.5 * n, 0.5))
        row("AC3a", "n=$n mismatching entries of W_n = [[1,n],[0,1]] and x_n = x0 + n p0 (bitwise)", bad, 0, pf(bad == 0))
    end
    n = 10_000
    res = propagate(sys, [x0], n; k = 1, tau = 1, seed = seed, W0 = [SMatrix{2,1,Float64}(0, 1)], reducers = (LogNorm(),))
    S = res.finals[1][1].logs[1]
    e = abs(S - 0.5 * log(float(n)^2 + 1))
    row("AC3b", "n=1e4 |S_n - 0.5 ln(n^2+1)|", e, 8 * n * epsf(Float64), pf(e <= 8 * n * epsf(Float64)))
    m = S / n - log(n) / n
    row("AC3c", "n=1e4 S_n/n - ln(n)/n (must be <= 1e-12; negative means below ln n / n)", m, 1e-12, pf(m <= 1e-12))
end

# ------------------------------------------------------------------------------------------------------ AC4 and oracle
# Same realized orbit: reference product of E8 at BigFloat(256) along the points of the T orbit (AC4); full BigFloat orbit (oracle).
function ac4_curves(T, K, nmax, seed)
    sys = StandardMap(T(K))
    x0 = sample_points(T, seed)
    W0 = [initial_deviations(seed, j, Val(2), Val(2), T) for j in 1:length(x0)]
    nj = length(x0)
    ratio = zeros(nmax)      # max_j err / (8 n eps prod||M||_F ||W0||_F)   (<= 1 passes)
    normerr = zeros(nmax)    # max_j err / (prod||M||_F ||W0||_F)           (plotted against 8 n eps)
    relerr = zeros(nmax)     # max_j err / ||W_ref||_F
    badx = 0
    setprecision(BigFloat, 256) do
        Kb = bigf(sys.K)
        refs = [Vector{Matrix{BigFloat}}(undef, nmax) for _ in 1:nj]
        prods = [Vector{BigFloat}(undef, nmax) for _ in 1:nj]
        zs = [Vector{SVector{2,T}}(undef, nmax) for _ in 1:nj]
        for j in 1:nj
            z = x0[j]
            Wr = bigmat(W0[j])
            pr = one(BigFloat)
            for n in 1:nmax
                c = Kb * cos(bigf(z[1]))
                M = [1+c one(c); c one(c)]
                pr *= fro(M)
                Wr = M * Wr
                z = forward(sys, z)
                refs[j][n] = Wr
                prods[j][n] = pr
                zs[j][n] = z
            end
        end
        for n in 1:nmax
            res = propagate(sys, x0, n; k = 2, tau = 1, seed = seed, W0 = W0)
            for j in 1:nj
                badx += !isequal(res.x[j], zs[j][n])
                err = fro(bigmat(res.W[j]) - refs[j][n])
                den = prods[j][n] * fro(bigmat(W0[j]))
                normerr[n] = max(normerr[n], Float64(err / den))
                ratio[n] = max(ratio[n], Float64(err / (8 * n * bigf(eps(T)) * den)))
                relerr[n] = max(relerr[n], Float64(err / fro(refs[j][n])))
            end
        end
    end
    return (; ratio, normerr, relerr, badx)
end

function oracle_curve(T, K, nmax, seed)
    sys = StandardMap(T(K))
    x0 = sample_points(T, seed)
    W0 = [initial_deviations(seed, j, Val(2), Val(2), T) for j in 1:length(x0)]
    nj = length(x0)
    rel = zeros(nmax)
    setprecision(BigFloat, 256) do
        sb = StandardMap(bigf(sys.K))
        refs = [Vector{Matrix{BigFloat}}(undef, nmax) for _ in 1:nj]
        for j in 1:nj
            z = SVector{2,BigFloat}(bigf(x0[j][1]), bigf(x0[j][2]))
            Wr = bigmat(W0[j])
            for n in 1:nmax
                Wr = Matrix(tangent(sb, z)) * Wr
                z = forward(sb, z)
                refs[j][n] = Wr
            end
        end
        for n in 1:nmax
            res = propagate(sys, x0, n; k = 2, tau = 1, seed = seed, W0 = W0)
            for j in 1:nj
                rel[n] = max(rel[n], Float64(fro(bigmat(res.W[j]) - refs[j][n]) / fro(refs[j][n])))
            end
        end
    end
    return rel
end

const AC4_DATA = Dict{Tuple{Any,Float64},Any}()
let
    seed = SEED + 4
    for (K, nmax) in ((0.5, 100), (2pi, 40)), T in (Float64, Double64)
        AC4_DATA[(T, K)] = ac4_curves(T, K, nmax, seed + (T === Float64 ? 0 : 1000) + nmax)
    end
    for T in (Float64, Double64), (K, n) in ((0.5, 100), (2pi, 20))
        d = AC4_DATA[(T, K)]
        r = maximum(d.ratio[n])
        row("AC4", "T=$T K=$(K == 0.5 ? "0.5" : "2pi") n=$n max over 20 orbits of ||W_n-W_ref||_F / (8 n eps prod||M||_F ||W0||_F) [orbit mismatches vs forward: $(d.badx)]",
            r, 1.0, pf(r <= 1 && d.badx == 0))
    end
    for T in (Float64, Double64), K in (0.5, 2pi)
        d = AC4_DATA[(T, K)]
        r = maximum(d.ratio)
        row("AC4", "T=$T K=$(K == 0.5 ? "0.5" : "2pi") max of the same ratio over all n in 1..$(length(d.ratio)) (figure data)", r, 1.0, pf(r <= 1))
    end
end

const ORACLE = Dict{Any,Vector{Float64}}()
let
    seed = SEED + 4
    nmax = 100
    for T in (Float64, Double64)
        rel = oracle_curve(T, 0.5, nmax, seed + (T === Float64 ? 0 : 1000) + nmax)
        ORACLE[T] = rel
        bound = 16 * 0.5 * nmax^2 * epsf(T)
        row("ORACLE-BF", "T=$T K=0.5 n=100 max over 20 orbits of ||W_n(T orbit) - W_n(full BigFloat(256) orbit)||_F / ||W_n^BF||_F", rel[nmax], bound, pf(rel[nmax] <= bound))
        worst = maximum(rel[n] / (16 * 0.5 * n^2 * epsf(T)) for n in 1:nmax)
        row("ORACLE-BF", "T=$T K=0.5 max over n=1..100 of (relative error / (16 K n^2 eps))", worst, 1.0, pf(worst <= 1))
    end
    row("ORACLE-CT1", "ChaosTools regular orbit K=0.5 n=1e5 |lambda_1 Chaoslab - ChaosTools| <= 1e-4", "-", 1e-4, "DEFERRED lyap-benettin 04")
    row("ORACLE-CT2", "ChaosTools chaotic ensemble K=2pi 1000 points n=1e4 difference of means <= 4 SE", "-", "4 SE", "DEFERRED lyap-benettin 04")
end

# ------------------------------------------------------------------------------------------------------------------ AC5
let
    sys = StandardMap(2pi)
    x0 = SVector(1.0, 0.5)
    W0 = SMatrix{2,2,Float64}(1, 0, 0, 1)
    b1 = allocbytes(sys, x0, W0, 10_000, 1, (LogNorm(), Reader(1)))
    b4 = allocbytes(sys, x0, W0, 10_000, 4, (LogNorm(), Reader(4)))
    b0 = allocbytes(sys, x0, W0, 100, 1, ())
    row("AC5", "bytes allocated by propagate_orbit after warm-up, n=1e4, k=2, chain (LogNorm, Reader), tau=1 / tau=4 / empty chain n=100", "$b1 / $b4 / $b0", 0, pf(b1 == 0 && b4 == 0 && b0 == 0))
end

# --------------------------------------------------------------------------------------------- child (1 thread) for AC6, AC13
const CHILD_RES = let
    path = tempname()
    cmd = addenv(`$(Base.julia_cmd()) --project=$(dirname(Base.active_project())) -t 1 $(@__FILE__) child $path`, "JULIA_NUM_THREADS" => "1")
    println("# starting 1-thread child run (AC6 digests and AC13 timings)")
    flush(stdout)
    run(cmd)
    deserialize(path)
end

# ------------------------------------------------------------------------------------------------------------------ AC6
let
    mine = ac6_digests()
    nt = Threads.nthreads()
    c1, cN = CHILD_RES[:ac6], mine
    agree = [c1[64] == cN[64], c1[1000] == cN[1000], c1[64] == c1[1000], cN[64] == cN[1000],
             cN[64] == cN[:serial], c1[:serial] == cN[:serial], cN[64] == c1[64]]
    ok = all(agree)
    row("AC6", "bitwise-equal digests (x, W, logs, reader state) of 4096 orbits n=1e3: $(CHILD_RES[:nthreads]) thread vs $nt threads, chunk 64 and 1000, and serial propagate_orbit; agreeing comparisons of 7",
        count(agree), 7, pf(ok))
end

# ------------------------------------------------------------------------------------------------------------------ AC7
let
    seed = SEED + 7
    n = 10_000
    for (d, k) in ((2, 1), (2, 2), (4, 4), (8, 8))
        worst = 0.0
        for j in 1:n
            Q = initial_deviations(seed, j, Val(d), Val(k), Float64)
            worst = max(worst, maximum(abs, Q' * Q - I))
        end
        row("AC7a", "Float64 d=$d k=$k max over 1e4 orbits ||Q'Q - I||_max", worst, 8 * k * epsf(Float64), pf(worst <= 8 * k * epsf(Float64)))
    end
    for (d, k) in ((2, 1), (2, 2), (4, 4), (8, 8))
        worst = 0.0
        for j in 1:200
            Q = initial_deviations(seed, j, Val(d), Val(k), Double64)
            worst = max(worst, Float64(maximum(abs, Q' * Q - I)))
        end
        row("AC7a", "Double64 d=$d k=$k max over 200 orbits ||Q'Q - I||_max", worst, 8 * k * epsf(Double64), pf(worst <= 8 * k * epsf(Double64)))
    end
    ref = [initial_deviations(seed, j, Val(2), Val(2), Float64) for j in 1:n]
    thr = Vector{SMatrix{2,2,Float64,4}}(undef, n)
    Threads.@threads for j in 1:n
        thr[j] = initial_deviations(seed, j, Val(2), Val(2), Float64)
    end
    ident = LinearSymplecticMap(SMatrix{2,2,Float64}(1, 0, 0, 1))
    chunkok = true
    for chunk in (64, 1000)
        res = propagate(ident, [SVector(0.0, 0.0) for _ in 1:300], 4; k = 2, tau = 2, seed = seed, chunk = chunk)
        chunkok &= all(j -> isequal(res.W[j], ref[j]), 1:300)
    end
    checks = [isequal(ref, [initial_deviations(seed, j, Val(2), Val(2), Float64) for j in 1:n]),
              isequal(ref, reverse([initial_deviations(seed, j, Val(2), Val(2), Float64) for j in n:-1:1])),
              isequal(ref, thr), chunkok, length(unique(Tuple.(ref))) == n,
              initial_deviations(seed + 1, 1, Val(2), Val(2), Float64) != ref[1]]
    row("AC7b", "bitwise reproducibility (repeat, reversed order, $(Threads.nthreads()) threads, chunk 64 and 1000), 1e4 distinct matrices, seed sensitivity; passing of 6", count(checks), 6, pf(all(checks)))
    u = [mod2pi(atan(Q[2, 1], Q[1, 1])) / (2pi) for Q in (initial_deviations(seed, j, Val(2), Val(1), Float64) for j in 1:n)]
    ks = ks_uniform(u)
    row("AC7c", "KS statistic of the angle of the single deviation vector, d=2 k=1, 1e4 orbits", ks, 1.63 / sqrt(n), pf(ks <= 1.63 / sqrt(n)))
end

# ------------------------------------------------------------------------------------------------------------------ AC8
let
    seed = SEED + 8
    x0s(flagged) = [SVector(0.0, 0.0, j == flagged ? 1.0 : 0.0) for j in 1:6]
    ok = propagate(Fault(:inf_state, 7), x0s(0), 20; k = 1, tau = 5, seed = seed, chunk = 2)
    row("AC8", "control without a flagged orbit raises no error and reaches x[1] = 20", all(x -> x[1] == 20 && isfinite(x[2]), ok.x), true, pf(all(x -> x[1] == 20 && isfinite(x[2]), ok.x)))
    for (id, map, flagged, step, what) in (("a", Fault(:inf_state, 7), 3, 10, "state"), ("b", Fault(:nan_W, 3), 5, 5, "deviation matrix"))
        res = nothing
        err = health_error(() -> (res = propagate(map, x0s(flagged), 20; k = 1, tau = 5, seed = seed, chunk = 2)))
        got = err === nothing ? "no error" : (orbit = err.orbit, step = err.step, what = err.what, backend = err.backend)
        exp = (orbit = flagged, step = step, what = what, backend = :cpu)
        row("AC8$id", "$(map.kind) injected at step $(map.at), tau=5: (orbit, step, what, backend), nothing stored: $(res === nothing)", got, exp, pf(got == exp && res === nothing))
    end
    mu = exp(3.0)
    sys = LinearSymplecticMap(SMatrix{2,2,Float64}(mu, 0, 0, 1 / mu))
    z0 = [SVector(0.0, 0.0)]
    w0 = [SMatrix{2,1,Float64}(1, 0)]
    r200 = propagate(sys, z0, 200; k = 1, tau = 200, seed = seed, W0 = w0)
    rel = abs(r200.W[1][1] / exp(600.0) - 1)
    row("AC8c", "tau=200 n=200: finite W, relative deviation of W[1] from e^600", rel, 1e-10, pf(isfinite(r200.W[1][1]) && rel <= 1e-10))
    for tau in (300, 100)
        err = health_error(() -> propagate(sys, z0, 300; k = 1, tau = tau, seed = seed, W0 = w0))
        got = err === nothing ? "no error" : (step = err.step, what = err.what, orbit = err.orbit)
        row("AC8c", "N=300, tau=$tau overflow e^900: (step, what, orbit)", got, (step = 300, what = "deviation matrix", orbit = 1), pf(got == (step = 300, what = "deviation matrix", orbit = 1)))
    end
    nthrow = 0
    for (N, tau) in ((7, 5), (0, 1), (-5, 5))
        try
            propagate(Boom(), z0, N; k = 1, tau = tau, seed = seed)
        catch e
            nthrow += e isa ArgumentError
        end
    end
    row("AC8d", "N not a positive multiple of tau ((7,5), (0,1), (-5,5)): ArgumentError before any work, count of 3", nthrow, 3, pf(nthrow == 3))
end

# ------------------------------------------------------------------------------------------------------------------ AC9
let
    x = 0.1
    n = 10^8
    exact = n * Rational{BigInt}(x)
    relerr(v) = Float64(abs(Rational{BigInt}(v) - exact) / exact)
    ec = relerr(csum(n, x))
    ep = relerr(plainsum(n, x))
    row("AC9a", "compensated sum of 1e8 copies of fl(0.1): relative error", ec, epsf(Float64), pf(ec <= epsf(Float64)))
    row("AC9b", "plain running sum, same input (negative control; must exceed the bound, i.e. FAIL (a))", ep, epsf(Float64), ep > 100 * ec + 1e-10 && ep > epsf(Float64) ? "PASS" : "FAIL")
    n7 = 10^7
    exact7 = n7 * Rational{BigInt}(x)
    c7 = csum(n7, Double64(x))
    e7 = Float64(abs((Rational{BigInt}(c7.hi) + Rational{BigInt}(c7.lo)) - exact7) / exact7)
    row("AC9c", "Double64 compensated sum of 1e7 copies of fl(0.1): relative error", e7, epsf(Double64), pf(e7 <= epsf(Double64)))
    acc = CompensatedSum{Float64}()
    for v in (1.0, 1e100, 1.0, -1e100)
        acc = add(acc, v)
    end
    row("AC9d", "Neumaier branch check [1, 1e100, 1, -1e100] sums to 2 (extra)", total(acc), 2.0, pf(total(acc) == 2.0))
end

# ----------------------------------------------------------------------------------------------------------------- AC10
let
    seed = SEED + 10
    good = () -> (loadable = true, functional = true, device = "TestGPU", compiles = true, gpu_faster = true)
    withp(; kw...) = () -> merge(good(), NamedTuple(kw))
    nocall = () -> error("probe must not be called")
    allbad = () -> (loadable = false, functional = false, device = "", compiles = false, gpu_faster = false)
    checks = Pair{String,Bool}[]
    chk(name, b) = push!(checks, name => b)
    chk("reasons list", FALLBACK_REASONS == ["CUDA.jl not loadable", "no functional device", "kernel does not compile for the element type",
        "device out of memory", "NumericalHealthError raised only on the GPU", "reducer is CPU-only", "GPU not faster for the element type", "launch exceeds the display watchdog"])
    c = resolve_backend(:auto, Float64)
    chk("(a) auto, no CUDA", c.backend === :cpu && c.device == "cpu" && c.reason in ("CUDA.jl not loadable", "no functional device"))
    default_reason = c.reason
    c = resolve_backend(:cpu, Float64; probe = good)
    chk("(b) cpu request", c.backend === :cpu && c.device == "cpu" && c.reason == "")
    for req in (:auto, :gpu), T in (Float64, Double64)
        c = resolve_backend(req, T; probe = good)
        chk("working device $req $T", c.backend === :gpu && c.device == "TestGPU" && c.reason == "")
    end
    c = resolve_backend(:gpu, Float64; probe = withp(functional = false))
    chk("(c) gpu without device", c.backend === :cpu && c.device == "cpu" && c.reason == "no functional device")
    for req in (:auto, :gpu)
        c = resolve_backend(req, BigFloat; probe = nocall)
        chk("(d) BigFloat $req", c.backend === :cpu && c.reason == "kernel does not compile for the element type")
    end
    chk("(e) loadable", resolve_backend(:auto, Float64; probe = withp(loadable = false)).reason == "CUDA.jl not loadable")
    chk("(e) functional", resolve_backend(:gpu, Float64; probe = withp(functional = false)).reason == "no functional device")
    chk("(e) compiles", resolve_backend(:gpu, Float64; probe = withp(compiles = false)).reason == "kernel does not compile for the element type")
    chk("(e) cpu-only reducer", resolve_backend(:auto, Float64; probe = good, cpu_only_reducer = true).reason == "reducer is CPU-only")
    chk("(e) gpu not faster", resolve_backend(:auto, Float64; probe = withp(gpu_faster = false)).reason == "GPU not faster for the element type")
    chk("order allbad", resolve_backend(:auto, Float64; probe = allbad).reason == "CUDA.jl not loadable")
    chk("order allbad+cpu-only", resolve_backend(:auto, Float64; probe = allbad, cpu_only_reducer = true).reason == "reducer is CPU-only")
    chk("order BigFloat first", resolve_backend(:auto, BigFloat; probe = nocall, cpu_only_reducer = true).reason == "kernel does not compile for the element type")
    ident = LinearSymplecticMap(SMatrix{2,2,Float64}(1, 0, 0, 1))
    x0 = [SVector(0.0, 0.0) for _ in 1:200]
    prov = nothing
    for req in (:auto, :gpu, :cpu)
        p = propagate(ident, x0, 4; k = 2, tau = 2, seed = seed, backend = req, chunk = 64).provenance
        prov = p
        chk("(f) fields $req", all(nm -> hasproperty(p, nm), (:backend, :device, :backend_fallback_reason, :precision, :threads, :cpu, :seed, :tau, :k, :chunk_size)))
        chk("(f) values $req", string(p.backend) in ("cpu", "gpu") && p.device == "cpu" && p.precision == "Float64" &&
            p.seed == seed && p.tau == 2 && p.k == 2 && p.chunk_size == 64 && p.threads >= 1 && !isempty(p.cpu) &&
            (req === :cpu ? p.backend_fallback_reason == "" : p.backend_fallback_reason in ("CUDA.jl not loadable", "no functional device")))
    end
    failed = [n for (n, b) in checks if !b]
    row("AC10", "backend resolution and provenance checks passing of $(length(checks)); default-probe reason here: \"$default_reason\"; failed: $(isempty(failed) ? "none" : join(failed, ", "))",
        length(checks) - length(failed), length(checks), pf(isempty(failed)))
    println("# AC10 provenance record of a run: ", prov)
end

# ----------------------------------------------------------------------------------------------------------------- AC11
let
    free = 2^30
    vals = [bytes_per_orbit(Val(2), Val(2), Float64), bytes_per_orbit(Val(4), Val(4), Float64),
            chunk_size(free, Val(2), Val(2), Float64), chunk_size(free, Val(4), Val(4), Float64),
            chunk_size(free ÷ 2, Val(2), Val(2), Float64), chunk_size(free ÷ 2, Val(4), Val(4), Float64)]
    exp = [96, 240, 5592405, 2236962, fld(2^28, 96), fld(2^28, 240)]
    row("AC11", "M(2,2,Float64), M(4,4,Float64), chunk at 2^30 (d=2,k=2 / d=4,k=4), chunk at 2^29 (d=2,k=2 / d=4,k=4)", vals, exp, pf(vals == exp))
end

# ------------------------------------------------------------------------------------------------------- AC12 (GPU only)
let
    for (id, txt, bound) in (("AC12a", "regular class K=0.5, 1024 orbits, n=1e5: max |S_cpu/n - S_gpu/n|", "1e-6"),
                             ("AC12c", "chaotic ensemble K=10, 1024 orbits, n=1e4: per-orbit <= 6 sigma_FT, difference of means <= 4 SE", "6 sigma_FT, 4 SE"))
        row(id, txt, "-", bound, "STAGE04-GPU")
    end
    for K in (2pi, 10.0)
        G = sqrt((1 + K)^2 + K^2 + 2)
        row("AC12b", "K=$(K == 2pi ? "2pi" : "10") n=10 same-orbit |x_cpu - x_gpu| (bound 4 n eps G^n, G=$(fmt(G)))", "-", 4 * 10 * epsf(Float64) * G^10, "STAGE04-GPU")
    end
end

# ----------------------------------------------------------------------------------------------------------------- AC13
let
    BLAS.set_num_threads(1)
    for T in (Float64, Double64)
        no, N = AC13_SIZES[T]
        t1 = CHILD_RES[(:ac13, T)]
        tN = ac13_time(T, no, N)
        steps = no * N
        scale = (10^6 * 10^4) / steps
        row("AC13", "CPU wall s, T=$T, $no orbits x $N steps (d=2, k=2, K=2pi, LogNorm reducer, tau=1), 1 thread: $(fmt(1e9 * t1 / steps)) ns/orbit-step; extrapolated to 1e6 x 1e4: $(fmt(t1 * scale)) s",
            t1, "GPU<=CPU", "STAGE04-GPU")
        row("AC13", "CPU wall s, T=$T, $no orbits x $N steps, $(Threads.nthreads()) threads: $(fmt(1e9 * tN / steps)) ns/orbit-step; extrapolated to 1e6 x 1e4: $(fmt(tN * scale)) s; speedup vs 1 thread $(fmt(t1 / tN))",
            tN, "GPU<=CPU", "STAGE04-GPU")
    end
end

# ----------------------------------------------------------------------------------------------------------------- AC14
let
    tdir = joinpath(ROOT, "src", "Tangent")
    code(l) = first(split(l, '#'))
    hits(pat, dir) = [(basename(f), i, strip(l)) for f in filter(f -> endswith(f, ".jl"), readdir(dir; join = true))
                      for (i, l) in enumerate(eachline(f)) if occursin(pat, code(l))]
    fwd = hits("forward(sys", tdir)
    tan = hits("tangent(sys", tdir)
    body_calls = filter(h -> !startswith(h[3], "@inline function _orbit_body"), hits("_orbit_body(", tdir))
    outside = [(relpath(f, ROOT), i, strip(l)) for f in filter(f -> endswith(f, ".jl"), vcat(
                   [readdir(joinpath(ROOT, "src", d); join = true) for d in ("Foundation", "Diagnostics")]...))
               for (i, l) in enumerate(eachline(f)) if occursin(r"forward\(|propagate\(", code(l)) && !occursin("\"\"\"", l)]
    println("# AC14 sites: forward(sys: ", fwd, "; tangent(sys: ", tan, "; _orbit_body callers: ", body_calls)
    println("# AC14 outside src/Tangent (Foundation, Diagnostics) forward/propagate mentions: ", isempty(outside) ? "none" : outside)
    ok = length(fwd) == 1 && length(tan) == 1 && fwd[1][1] == tan[1][1] && length(body_calls) == 2 &&
         body_calls[1][1:2] != body_calls[2][1:2]
    loc(h) = "$(h[1]):$(h[2])"
    row("AC14", "per-step update sites in src/Tangent (forward(sys, x) at $(loc(fwd[1])), tangent(sys, x) at $(loc(tan[1])), both inside `_orbit_body`); callers of `_orbit_body`: $(loc(body_calls[1])) and $(loc(body_calls[2])) (`_orbit_kernel!` and `propagate_orbit`)",
        "$(length(fwd)) / $(length(tan)) sites", "1 / 1", pf(ok))
end

# ------------------------------------------------------------------------------------------------------ W9 (Amendment D28)
let
    seed = SEED + 9
    Kw = 2pi
    sys = Base.invokelatest(V9, Kw, 0, 0)
    x0 = SVector(1.0, 0.5, 0.0, 0.0)
    avail = try
        Base.invokelatest(propagate_orbit, sys, x0, initial_deviations(seed, 1, Val(2), Val(1), Float64), 3, 1, ())
        true
    catch e
        println("# W9 availability probe failed: ", sprint(showerror, e)[1:min(end, 200)])
        false
    end
    if !avail
        for id in ("W9a", "W9b", "W9c", "W9d", "W9e", "W9f", "W9g")
            row(id, "state longer than the tangent dimension (Amendment D28)", "-", "-", "PENDING D28")
        end
    else
        # (a)
        bad = 0
        for k in (1, 2), tau in (1, 3)
            W0 = initial_deviations(seed, 1, Val(2), Val(k), Float64)
            xr, Wr, _ = Base.invokelatest(propagate_orbit, sys, x0, W0, 12, tau, ())
            z, Wm = x0, W0
            for _ in 1:12
                Wm = tangent(sys, z) * Wm
                z = forward(sys, z)
            end
            bad += (xr != z) + (Wr != Wm) + (length(xr) != 4) + (size(Wr) != (2, k))
        end
        row("W9a", "k in {1,2} x tau in {1,3}: mismatches of the 4-component state and the 2 x k matrix vs manual iteration (bitwise), 4 cases", bad, 0, pf(bad == 0))
        # (b)
        rng = StableRNG(seed)
        xs = [SVector(2pi * rand(rng), 2pi * rand(rng), 0.0, 0.0) for _ in 1:257]
        bad = 0
        for k in (1, 2), tau in (1, 3)
            Ws = [initial_deviations(seed, j, Val(2), Val(k), Float64) for j in 1:257]
            ser = [Base.invokelatest(propagate_orbit, sys, xs[j], Ws[j], 12, tau, (); orbit = j) for j in 1:257]
            for chunk in (16, nothing)
                res = Base.invokelatest(propagate, sys, xs, 12; k = k, tau = tau, seed = seed, chunk = chunk)
                bad += !isequal(res.x, [s[1] for s in ser]) + !isequal(res.W, [s[2] for s in ser])
            end
        end
        row("W9b", "257 seeded orbits, chunk 16 and default, k x tau cases: chunked != serial (of 8 comparisons x 2 arrays)", bad, 0, pf(bad == 0))
        # (c)
        res = Any[]
        for fault in (3, 4)
            sf = Base.invokelatest(V9, Kw, fault, 7)
            err = health_error(() -> Base.invokelatest(propagate, sf, [x0], 12; k = 1, tau = 3, seed = seed))
            push!(res, err === nothing ? "no error" : (what = err.what, step = err.step))
        end
        want = [(what = "state", step = 9), (what = "state", step = 9)]
        row("W9c", "component 3 / component 4 set to Inf at step 7, tau=3: (what, step)", res, want, pf(res == want))
        # (d)
        W0 = initial_deviations(seed, 1, Val(2), Val(2), Float64)
        b = Base.invokelatest(allocbytes, sys, x0, W0, 12, 3, (LenReader(),))
        b0 = Base.invokelatest(allocbytes, sys, x0, W0, 100, 1, ())
        row("W9d", "bytes allocated by propagate_orbit with s=4 after warm-up (chain (LenReader), n=12 / empty chain, n=100)", "$b / $b0", 0, pf(b == 0 && b0 == 0))
        # (e)
        r = Base.invokelatest(propagate, sys, [x0], 12; k = 1, tau = 3, seed = seed, reducers = (LenReader(),))
        row("W9e", "length of the point handed to a reducer", r.finals[1][1].len, 4, pf(r.finals[1][1].len == 4))
        # (f)
        v = try
            [Base.invokelatest(bytes_per_orbit, Val(2), Val(2), Float64, 4), bytes_per_orbit(Val(2), Val(2), Float64),
             Base.invokelatest((; kw...) -> chunk_size(2^30, Val(2), Val(2), Float64; kw...), s = 4),
             chunk_size(2^30, Val(2), Val(2), Float64)]
        catch e
            sprint(showerror, e)
        end
        wantf = [112, 96, fld(2^29, 112), 5592405]
        row("W9f", "bytes per orbit s=4 d=2 k=2 / default, chunk at 2^30 for s=4 / default", v, wantf, pf(v == wantf))
        # (g)
        thrown = try
            Base.invokelatest(propagate, V9Boom(), [x0], 3; k = 3, tau = 1, seed = seed)
            "no error"
        catch e
            e isa ArgumentError ? "ArgumentError" : sprint(showerror, e)
        end
        row("W9g", "k=3 with d=2 (state length 4) rejected before any work", thrown, "ArgumentError", pf(thrown == "ArgumentError"))
    end
    nfail = count(r -> r[2] == "FAIL" && !startswith(r[1], "W9"), ROWS)
    row("W9h", "existing s = d suites unchanged: FAIL lines among AC1 to AC14 above (package suite is run by test/runtests.jl)", nfail, 0, pf(nfail == 0))
end

# ------------------------------------------------------------------------------------------------------------- figures
function make_figures()
    mkpath(OUT)
    cols = Dict(Float64 => :dodgerblue3, Double64 => :darkorange2)
    f1 = Figure(size = (1100, 460))
    for (i, (K, ttl)) in enumerate(((0.5, "K = 0.5"), (2pi, "K = 2π")))
        ax = Axis(f1[1, i]; xscale = log10, yscale = log10, xlabel = "n (map steps)",
                  ylabel = i == 1 ? "‖W_n − W_n^ref‖_F / (∏‖M_i‖_F ‖W_0‖_F)" : "", title = "Standard map, $ttl, 20 seeded orbits (max)")
        for T in (Float64, Double64)
            d = AC4_DATA[(T, K)]
            ns = collect(1:length(d.normerr))
            keep = d.normerr .> 0
            lines!(ax, ns[keep], d.normerr[keep]; color = cols[T], label = "$T error")
            kr = d.relerr .> 0
            lines!(ax, ns[kr], d.relerr[kr]; color = cols[T], linestyle = :dot, label = "$T relative to ‖W_ref‖")
            lines!(ax, ns, 8 .* ns .* epsf(T); color = cols[T], linestyle = :dash, label = "$T bound 8 n ε")
        end
        i == 2 && axislegend(ax; position = :lt, labelsize = 11)
    end
    save(joinpath(OUT, "tangent-core-ac4-bigfloat-product.png"), f1)
    f2 = Figure(size = (620, 460))
    ax = Axis(f2[1, 1]; xscale = log10, yscale = log10, xlabel = "n (map steps)", ylabel = "‖W_n − W_n^BF‖_F / ‖W_n^BF‖_F",
              title = "Standard map K = 0.5, T orbit vs full BigFloat(256) orbit, 20 orbits (max)")
    for T in (Float64, Double64)
        rel = ORACLE[T]
        ns = collect(1:length(rel))
        keep = rel .> 0
        lines!(ax, ns[keep], rel[keep]; color = cols[T], label = "$T")
        lines!(ax, ns, 16 .* 0.5 .* ns .^ 2 .* epsf(T); color = cols[T], linestyle = :dash, label = "$T bound 16 K n² ε")
    end
    axislegend(ax; position = :lt, labelsize = 11)
    save(joinpath(OUT, "tangent-core-oracle-fullbigfloat.png"), f2)
    return nothing
end
Base.invokelatest(make_figures)
println("# figures: tangent-core-ac4-bigfloat-product.png, tangent-core-oracle-fullbigfloat.png in ", OUT)

# ------------------------------------------------------------------------------------------------------------- summary
verdicts = Dict{String,Int}()
for (_, v) in ROWS
    verdicts[v] = get(verdicts, v, 0) + 1
end
println("# SUMMARY ", join(["$k=$v" for (k, v) in sort!(collect(verdicts))], ", "), " (", length(ROWS), " lines)")
