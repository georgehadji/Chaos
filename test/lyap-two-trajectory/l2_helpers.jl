# Helpers for the lyap-two-trajectory suite (prefix l2_). Definitions only: nothing here calls a Lyapunov name at load time,
# so a missing Chaoslab.Lyapunov errors the testsets that need it and nothing else.
using Test, Chaoslab, StaticArrays, LinearAlgebra, StableRNGs

# Amplification sum Psi_n (lyap-benettin S3, Amendment D35): one definition, in the lyap-benettin folder
@isdefined(LB_PSI_LOADED) || include(joinpath(@__DIR__, "..", "lyap-benettin", "lb_psi.jl"))

const L2_SEED = 20261012                  # base seed, _config/numerics-defaults.toml
const L2_BASE = L2_SEED + 25000           # this slug uses 20261012 + 25000 .. 25999; testset AC<n> uses 50 n + sub (sub < 50)
l2_seed(ac::Integer, sub::Integer = 0) = L2_BASE + 50 * ac + sub
const L2_TWO_PI = 2 * pi

# ---- BigFloat helpers (call inside setprecision(BigFloat, b) do ... end)
l2_big(x::Double64) = BigFloat(x.hi) + BigFloat(x.lo)
l2_big(x::Real) = BigFloat(x)
l2_fibbig(m) = (a = big(0); b = big(1); for _ in 1:m; a, b = b, a + b; end; a)
# ln ||A^n e_1|| for A = [[2,1],[1,1]] (A^n e_1 = (F_{2n+1}, F_{2n}))
l2_ln_cat(n) = setprecision(BigFloat, 256) do
    log(sqrt(BigFloat(l2_fibbig(2n + 1))^2 + BigFloat(l2_fibbig(2n))^2))
end

# ---- tolerance formulas of the spec (T5, T6, T8, L1 S3); eps is the effective eps of the element type
l2_BS(::Type{T}, n, G_first, K = 0) where {T} = 8 * Float64(effective_eps(T)) * (n + G_first) + 16 * abs(K) * n^2 * Float64(effective_eps(T))   # L1 S3 as amended (K_s = 0 for constant tangents)
function l2_sigma(K)
    S2 = 2 * K^2 + 2 * abs(K) + 5
    S = sqrt(S2)
    return (S + sqrt(S2 - 4)) / 2
end
l2_Enl(K, d0) = abs(K) * l2_sigma(K) * d0 / sqrt(2)
l2_Ero(K, d0, Mstate, epsT) = (abs(K) + 4) * Mstate * epsT * l2_sigma(K) / d0
l2_D2T(K, d0, Mstate, epsT) = 2 * (l2_Enl(K, d0) + l2_Ero(K, d0, Mstate, epsT))
l2_d0star(K, Mstate, epsT) = sqrt(sqrt(2) * (abs(K) + 4) * Mstate * epsT / abs(K))
l2_valid(K, d0, tau, Mstate, epsT) = d0 * l2_sigma(K)^tau <= 0.1 && l2_Ero(K, d0, Mstate, epsT) <= 0.5       # T5, T6 validity

# ---- orbit side (iterating `forward` in T is bitwise the stream's orbit)
l2_sys(mode::Symbol, K) = mode === :wrapped ? WrappedStandardMap(K) : StandardMap(K)
l2_state(mode::Symbol, ::Type{T}, z) where {T} = mode === :wrapped ? torus_state(SVector{2, T}(T(z[1]), T(z[2]))) : SVector{2, T}(T(z[1]), T(z[2]))
l2_starts(rng, m) = [SVector(L2_TWO_PI * rand(rng), L2_TWO_PI * rand(rng)) for _ in 1:m]
function l2_Gfirst(sys, x0, tau)          # squared Frobenius norm of the first-window tangent product (L1 S3)
    x, P = x0, one(tangent(sys, x0))
    for _ in 1:tau
        P = tangent(sys, x) * P
        x = forward(sys, x)
    end
    return sum(abs2, Float64(P[i]) for i in 1:length(P))
end
function l2_Mstate(sys, x0, N)             # max(1, |K|, largest lifted coordinate magnitude) along the orbit (lifted mode)
    x, m = x0, max(1.0, abs(Float64(sys.K)))
    for _ in 1:N
        m = max(m, Float64(abs(x[1])), Float64(abs(x[2])))
        x = forward(sys, x)
    end
    return max(m, Float64(abs(x[1])), Float64(abs(x[2])))
end

# L1 series X(n) of the identical orbit and direction: Benettin with k = 1 (first column of the seeded frame)
l2_X(sys, x0s, n, seed; tau = 1) = [f[1].lambda[1] for f in propagate(sys, x0s, n; k = 1, tau = tau, seed = seed, reducers = (Benettin(),)).finals]
l2_X2(sys, x0s, n, seed, d0; tau = 1, U0 = nothing) = [f[1].X2 for f in propagate_shadow(sys, x0s, n; d0 = d0, tau = tau, seed = seed, U0 = U0).finals]

# S7 regular class: X(1e5) <= 5e-4, from the L1 estimator itself; returns the first m passing starts
function l2_regular_starts(K, m, rng)
    cand = l2_starts(rng, 3m)
    X = l2_X(StandardMap(K), cand, 100_000, L2_SEED)
    ok = X .<= 5e-4
    @assert count(ok) >= m "fewer than $m regular starts among $(3m) candidates"
    return cand[ok][1:m]
end

# median and least-squares slope
function l2_median(v)
    s = sort(v)
    n = length(s)
    return isodd(n) ? s[(n + 1) ÷ 2] : (s[n ÷ 2] + s[n ÷ 2 + 1]) / 2
end
function l2_slope(x, y)
    xm, ym = sum(x) / length(x), sum(y) / length(y)
    return sum((x .- xm) .* (y .- ym)) / sum((x .- xm) .^ 2)
end
l2_mean(a) = sum(a) / length(a)

# ---- test maps (Batch A names only)
# Fixed-state map for the health cases, state (a, b), d = 2: a advances by 1; on an orbit whose b exceeds 1/2 and once a reaches `thr` the
# `kind` happens: :collapse sends the point to the constant (thr, 0) (both orbits of the pair land on one point), :inf_ref / :inf_shadow /
# :inf_both set a to Inf for the reference (b == 1 exactly), the shadow (b != 1) or both.
struct l2_Fault <: AbstractMap
    kind::Symbol
    thr::Float64
end
Chaoslab.dimension(::l2_Fault) = 2
function Chaoslab.forward(m::l2_Fault, z)
    a = z[1] + 1
    (z[2] > 0.5 && a >= m.thr) || return SVector(a, z[2])
    m.kind === :collapse && return SVector(m.thr, 0.0)
    hit = m.kind === :inf_both || (m.kind === :inf_ref && z[2] == 1.0) || (m.kind === :inf_shadow && z[2] != 1.0)
    return hit ? SVector(Inf, z[2]) : SVector(a, z[2])
end

# Map whose forward throws: proves that a rejected d0 stops before any work.
struct l2_Boom <: AbstractMap end
Chaoslab.dimension(::l2_Boom) = 2
Chaoslab.forward(::l2_Boom, z) = error("l2_Boom: forward was called")

# Test double with the PLAIN wrapped difference (the artifact of T3): the wrapped standard map, but `separation` subtracts the wrapped coordinates.
struct l2_Plain{T} <: AbstractMap
    inner::WrappedStandardMap{T}
end
Chaoslab.dimension(::l2_Plain) = 2
Chaoslab.forward(m::l2_Plain, z) = forward(m.inner, z)
Chaoslab.tangent(m::l2_Plain, z) = tangent(m.inner, z)
@static if isdefined(Chaoslab, :Lyapunov) && isdefined(Chaoslab.Lyapunov, :separation)
    Chaoslab.Lyapunov.separation(::l2_Plain, ref, shadow) = SVector(shadow[1] - ref[1], shadow[2] - ref[2])
    Chaoslab.Lyapunov.displaced(m::l2_Plain, ref, delta) = Chaoslab.Lyapunov.displaced(m.inner, ref, delta)
end

function l2_health_error(f)
    try
        f()
        return nothing
    catch e
        e isa NumericalHealthError ? e : rethrow()
    end
end

# Psi_n of every start (T5a, L1 S3 as amended): the reference orbit of x0s[j] in the element type and the seeded direction of orbit j
# (the first column of the seeded frame, the direction `propagate_shadow` uses by default); the orbit ids are the positions in x0s.
l2_Psis(sys, x0s, seed, ::Type{T}, n) where {T} =
    [lb_Psi(sys, x0s[j], initial_deviations(seed, j, Val(2), Val(1), T)[:, 1], n) for j in eachindex(x0s)]
# Pathwise bound of T5a / T6a: factor * Delta_2T * Psi_n / n + B_S^L1(n) / n, factor = sigma^(2 (tau - 1)) (1 for tau = 1)
l2_pathbound(D2T, Psi, n, BS; factor = 1.0) = factor * D2T * Psi / n + BS / n
