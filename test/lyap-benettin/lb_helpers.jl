# Helpers for the lyap-benettin suite (prefix lb_). Definitions only: nothing here calls a Lyapunov name at load time.
# Spec: workspaces/module-build/stages/01-spec/output/lyap-benettin-spec.md. API contract: M2 Batch B api-batch-b.md (D34).
using Test, Chaoslab, StaticArrays, LinearAlgebra, StableRNGs, DoubleFloats

const LB_SEED = 20261012          # base seed, _config/numerics-defaults.toml key `seed`
const LB_OFF = 21000              # this slug uses seed offsets 21000 to 21999 (base + LB_OFF + 10 * criterion + sub-index)
const LB_TYPES = (Float64, Double64)   # Float32 dropped (D24)
const LB_BITS = 256               # numerics-defaults.toml key `bigfloat_bits`
const LB_DET_TOL = 4.0            # numerics-defaults.toml key `det_residual_tol`
include(joinpath(@__DIR__, "lb_psi.jl"))   # Psi_n (S3, Amendment D35), shared with the householder and two-trajectory suites

lb_eps(::Type{T}) where {T} = Float64(effective_eps(T))     # effective eps: Float64 eps, Double64 2^-100
lb_BS(epsT, n, G, K = 0) = 8 * epsT * (n + G) + 16 * abs(K) * n^2 * epsT   # S3, Amendment D35: B_S; K = 0 where the tangent is constant (K_s)

# ------------------------------------------------------------------ BigFloat conversions and references

lb_big(x::Double64) = BigFloat(x.hi) + BigFloat(x.lo)       # exact; call inside setprecision(BigFloat, LB_BITS)
lb_big(x::AbstractFloat) = BigFloat(x)

# |value - ref| for value = S + C given in T (exact widening), evaluated at 256 bits, returned as Float64.
lb_err(S, C, ref) = setprecision(BigFloat, LB_BITS) do
    Float64(abs(lb_big(S) + lb_big(C) - ref))
end

# Orbit of the standard map in T: pts[i] = x_{i-1} (the point before step i, where M_i is evaluated), xN = x_n.
function lb_pts(sys, x0, n::Integer)
    pts = Vector{typeof(x0)}(undef, n)
    z = x0
    for i in 1:n
        pts[i] = z
        z = forward(sys, z)
    end
    return pts, z
end

# Reference ln(norm(Y_k w0) / norm(w0)) at k = every, 2 every, ... (E2, E8): Y_k = M_{k-1} ... M_0 with M_i the standard-map
# tangent [[1 + c, 1], [c, 1]], c = K cos x_i (E1), evaluated at BigFloat(256) on the realized points `pts`.
function lb_ref_std(Kt, pts, w0, every::Int)
    return setprecision(BigFloat, LB_BITS) do
        Kb = lb_big(Kt)
        w1, w2 = lb_big(w0[1]), lb_big(w0[2])
        n0 = log(w1^2 + w2^2) / 2
        out = BigFloat[]
        for (i, z) in enumerate(pts)
            c = Kb * cos(lb_big(z[1]))
            w1, w2 = (1 + c) * w1 + w2, c * w1 + w2
            i % every == 0 && push!(out, log(w1^2 + w2^2) / 2 - n0)
        end
        out
    end
end

# Worst of |S + C - ref[n / every]| / B_S(n, G, K) over the snapshots (n, S, C), in one 256-bit block.
function lb_worst_ratio(snaps, ref, every::Int, epsT, G, K)
    return setprecision(BigFloat, LB_BITS) do
        worst = 0.0
        for (n, S, C) in snaps
            worst = max(worst, Float64(abs(lb_big(S) + lb_big(C) - ref[div(n, every)])) / lb_BS(epsT, n, G, K))
        end
        worst
    end
end

# Full BigFloat(256) orbit and tangent product for K = Kf (Float64, exactly representable), start x0 widened exactly.
# Returns (|S + C - ref|, rho_n, M_state) with ref = ln(norm(Y w0)/norm(w0)), rho_n = sigma_max(Y) / norm(Y w0) and M_state = max(1, |K|,
# largest coordinate magnitude on the lifted orbit), both coordinates (S1; the state-magnitude factor of criteria 7 and 11, Amendment D35).
function lb_ref_full(Kf::Float64, x0, w0, n::Int, S, C)
    return setprecision(BigFloat, LB_BITS) do
        K = BigFloat(Kf)
        x, p = lb_big(x0[1]), lb_big(x0[2])
        a, b, c_, d = one(K), zero(K), zero(K), one(K)
        M = max(1.0, abs(Kf))
        for _ in 1:n
            M = max(M, Float64(abs(x)), Float64(abs(p)))
            c = K * cos(x)
            a, b, c_, d = (1 + c) * a + c_, (1 + c) * b + d, c * a + c_, c * b + d
            p = p + K * sin(x)
            x = x + p
        end
        w1, w2 = lb_big(w0[1]), lb_big(w0[2])
        y1, y2 = a * w1 + b * w2, c_ * w1 + d * w2
        ref = (log(y1^2 + y2^2) - log(w1^2 + w2^2)) / 2
        F = a^2 + b^2 + c_^2 + d^2
        dt = a * d - b * c_
        smax = sqrt((F + sqrt(F^2 - 4 * dt^2)) / 2)
        M = max(M, Float64(abs(x)), Float64(abs(p)))
        (Float64(abs(lb_big(S) + lb_big(C) - ref)), Float64(smax / sqrt(y1^2 + y2^2)), M)
    end
end

# Squared Frobenius norm of the first-window product (S3), from the Float64 tangents Ms[i] at pts[i].
function lb_Gfirst(Ms, tau::Int)
    P = one(SMatrix{2, 2, Float64})
    for i in 1:tau
        P = Ms[i] * P
    end
    return sum(abs2, P)
end
# Phi_j: squared Frobenius norm of the product over window j (steps (j-1) tau + 1 .. j tau).
function lb_Phis(Ms, tau::Int)
    out = Float64[]
    for j in 1:div(length(Ms), tau)
        P = one(SMatrix{2, 2, Float64})
        for i in ((j - 1) * tau + 1):(j * tau)
            P = Ms[i] * P
        end
        push!(out, sum(abs2, P))
    end
    return out
end

# E14: sigma(K) = (S + sqrt(S^2 - 4)) / 2, S^2 = 2 K^2 + 2 |K| + 5, evaluated in the element type of K.
function lb_sigma(K::T) where {T}
    S2 = 2 * K^2 + 2 * abs(K) + 5
    return (sqrt(S2) + sqrt(S2 - 4)) / 2
end

# Range invariant (criterion 18): worst of alpha / (sigma^tau (1 + delta_j)) and sigma^-tau (1 - delta_j) / alpha over all windows,
# delta_j = 16 eps tau Phi_j. Pass if <= 1. alphas are in T; the comparison is made in T.
function lb_range_violation(alphas, Phis, K::T, tau::Int, epsT::Float64) where {T}
    s = lb_sigma(K)
    hi0 = s^tau
    lo0 = inv(hi0)
    worst = 0.0
    for (a, Phi) in zip(alphas, Phis)
        delta = T(16 * epsT * tau * Phi)
        worst = max(worst, Float64(a / (hi0 * (1 + delta))), Float64(lo0 * (1 - delta) / a))
    end
    return worst
end

# ------------------------------------------------------------------ test maps (never in src/)

# Fault map on the state (n, y, flag), dimension 2 (one type, so the stream kernel compiles once): the counter n advances by one per
# step; only flagged orbits misbehave. kind :zero: the tangent is the zero matrix (the deviation vector collapses at the first
# handoff, criterion 9a); kind :nan: the tangent has a NaN at counter value `at` (criterion 9b).
struct lb_Fault <: Chaoslab.AbstractMap
    kind::Symbol
    at::Int
end
Chaoslab.dimension(::lb_Fault) = 2
Chaoslab.forward(::lb_Fault, z) = SVector(z[1] + 1, z[2], z[3])
function Chaoslab.tangent(m::lb_Fault, z::SVector{3, T}) where {T}
    z[3] == 1 || return one(SMatrix{2, 2, T})
    m.kind === :zero && return zero(SMatrix{2, 2, T})
    return SMatrix{2, 2, T}(z[1] == m.at ? T(NaN) : one(T), 0, 0, 1)
end

# Map whose forward and tangent throw: proves argument validation happens before any work.
struct lb_Boom <: Chaoslab.AbstractMap end
Chaoslab.dimension(::lb_Boom) = 2
Chaoslab.forward(::lb_Boom, z) = error("lb_Boom: forward was called")
Chaoslab.tangent(::lb_Boom, z) = error("lb_Boom: tangent was called")

# ------------------------------------------------------------------ test reducers (protocol of tangent-core / Batch B contract)

# Wraps a reducer: delegates init, step, failed, error to `inner`; before the inner step it records alpha_j = norm of the
# (unnormalized, k = 1) deviation matrix, and after it the (n, sum, comp) of the inner final at the handoff steps in `steps`
# (or at every handoff when `all`); `alphas = false` skips the alpha record. CPU only: the state holds Vectors, so use it with
# propagate_orbit. Wrapping Benettin in it everywhere keeps the number of distinct compiled argument-type combinations small.
struct lb_Probe{R}
    inner::R
    steps::Vector{Int}
    all::Bool
    alphas::Bool
end
lb_Probe(inner, steps::AbstractVector{<:Integer} = Int[]; all::Bool = false, alphas::Bool = true) =
    lb_Probe(inner, collect(Int, steps), all, alphas)

function Chaoslab.reducer_init(r::lb_Probe, ::Type{T}, v::Val, k::Val) where {T}
    return (Chaoslab.reducer_init(r.inner, T, v, k), 1, Tuple{Int, T, T}[], T[])
end
function Chaoslab.reducer_step(r::lb_Probe, st, n::Int, x, W)
    inner, pos, snaps, alphas = st
    r.alphas && push!(alphas, sqrt(sum(abs2, W)))
    inner2, W2 = Chaoslab.reducer_step(r.inner, inner, n, x, W)
    if r.all || (pos <= length(r.steps) && n == r.steps[pos])
        f = Chaoslab.reducer_final(r.inner, inner2)
        push!(snaps, (n, f.sum[1], f.comp[1]))
        pos += 1
    end
    return (inner2, pos, snaps, alphas), W2
end
Chaoslab.reducer_failed(r::lb_Probe, st) = Chaoslab.reducer_failed(r.inner, st[1])
Chaoslab.reducer_error(r::lb_Probe, st, orbit, backend) = Chaoslab.reducer_error(r.inner, st[1], orbit, backend)
Chaoslab.reducer_final(r::lb_Probe, st) = (inner = Chaoslab.reducer_final(r.inner, st[1]), snaps = st[3], alphas = st[4])

# Read-only reducer recording the orbit point at every handoff (criterion 10a).
struct lb_Trace end
Chaoslab.reducer_init(::lb_Trace, ::Type, ::Val, ::Val) = Any[]
function Chaoslab.reducer_update(::lb_Trace, st, n::Int, x, W)
    push!(st, x)
    return st, W
end
Chaoslab.reducer_final(::lb_Trace, st) = (xs = st,)

# ------------------------------------------------------------------ small utilities

# Run f(); return the NumericalHealthError it throws (nothing if it returns), rethrow anything else.
function lb_health_error(f)
    try
        f()
        return nothing
    catch e
        e isa NumericalHealthError ? e : rethrow()
    end
end

# K values of the panel, built in T without a Float64 detour (precision spec Algorithm 2).
lb_K(::Type{T}, label::Symbol) where {T} =
    label === :two_pi ? exact(T, pi) * 2 : label === :half ? exact(T, 0.5) : label === :ten ? exact(T, 10) : error("unknown K $label")

# Float64 starts uniform in [0, 2 pi)^2 from a StableRNG, widened exactly into T by the caller.
lb_starts(rng, m::Integer) = [SVector(2pi * rand(rng), 2pi * rand(rng)) for _ in 1:m]
lb_cast(::Type{T}, z::SVector{2}) where {T} = SVector{2, T}(T(z[1]), T(z[2]))

# Seeded unit deviation vector of orbit j as an SVector{2,T} (the column the stream starts from, k = 1).
lb_w0(seed, j, ::Type{T}) where {T} = initial_deviations(seed, j, Val(2), Val(1), T)[:, 1]

# Bytes allocated by one propagate_orbit call after a warm-up call, measured inside a function.
function lb_alloc_bytes(sys, x0, W0, N, tau, reducers)
    propagate_orbit(sys, x0, W0, N, tau, reducers)
    return @allocated propagate_orbit(sys, x0, W0, N, tau, reducers)
end

# Serial stand-in for `propagate` on the CPU: the same per-orbit body (propagate_orbit), orbit ids 1:length(x0), the first failing
# orbit in id order throws. Used where the kernel path adds nothing to the criterion (it costs one KernelAbstractions compile per
# argument-type combination); the kernel path itself is tested in AC9, AC10 and the provenance test. With `probe = true` and
# Benettin first in the chain, Benettin runs inside an lb_Probe (forwarding every call) and finals[j][1] is Benettin's own final.
function lb_run(sys, x0::AbstractVector{SVector{s, T}}, N::Integer; k::Integer = 1, tau::Integer = 1, seed::Integer,
                reducers::Tuple = (Benettin(),), W0 = nothing, probe::Bool = true) where {s, T}
    d = dimension(sys)
    wrap = probe && !isempty(reducers) && first(reducers) isa Benettin
    chain = wrap ? (lb_Probe(first(reducers); alphas = false), Base.tail(reducers)...) : reducers
    out = map(eachindex(x0)) do j
        W = W0 === nothing ? initial_deviations(seed, j, Val(d), Val(k), T) : W0[j]
        x, W1, st = propagate_orbit(sys, x0[j], W, N, tau, chain; orbit = j)
        f = map(reducer_final, chain, st)
        (x, W1, wrap ? (first(f).inner, Base.tail(f)...) : f)
    end
    return (x = [o[1] for o in out], W = [o[2] for o in out], finals = [o[3] for o in out])
end

# Benettin's final for one orbit (probe-wrapped chain, see lb_run).
lb_final(sys, x0, W0, N, tau; orbit::Integer = 1) = lb_run(sys, [x0], N; tau = tau, seed = 0, W0 = [W0]).finals[1][1]

# Regular-class starts at K = 0.5 (S7), shared by AC7, AC11 and AC17: candidates c = 1, 2, ... take their (x, p) from the c-th draw of
# a StableRNG (seed base + LB_OFF + 70) and the orbit id c for the seeded deviation vector; a candidate is regular when
# X(1e5) <= 5e-4 (Float64, tau = 1, from the estimator under test). Returns (ids, starts, rejected, seed) for the first 50 regular
# candidates; rejected counts the replaced ones (reported, not hidden). Candidates are screened in batches of 60.
function lb_regular()
    return get!(LB_CACHE, :regular) do
        seed = LB_SEED + LB_OFF + 70
        sys = StandardMap(0.5)
        rng = StableRNG(seed)
        cands = SVector{2, Float64}[]
        ok = Int[]
        while length(ok) < 50
            lo = length(cands) + 1
            append!(cands, lb_starts(rng, 60))
            W0 = [initial_deviations(seed, c, Val(2), Val(1), Float64) for c in lo:length(cands)]   # orbit id = candidate number
            res = lb_run(sys, cands[lo:end], 100_000; k = 1, tau = 1, seed = seed, W0 = W0)
            append!(ok, (lo - 1) .+ findall(<=(5.0e-4), [f[1].lambda[1] for f in res.finals]))
        end
        keep = ok[1:50]
        (ids = keep, starts = cands[keep], rejected = last(keep) - 50, seed = seed)
    end
end

# Memoized computations shared by several criteria (criteria 5, 6, 18 use the same runs).
const LB_CACHE = Dict{Symbol, Any}()
