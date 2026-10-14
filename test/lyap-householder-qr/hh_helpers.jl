# Helpers for the lyap-householder-qr suite (prefix hh_). Definitions only: nothing here calls a Lyapunov name at load time,
# so a missing Chaoslab.Lyapunov errors the testsets that need it and nothing else.
using Test, Chaoslab, StaticArrays, LinearAlgebra, StableRNGs

# Amplification sum Psi_n (lyap-benettin S3, Amendment D35): one definition, in the lyap-benettin folder
@isdefined(LB_PSI_LOADED) || include(joinpath(@__DIR__, "..", "lyap-benettin", "lb_psi.jl"))

const HH_SEED = 20261012                  # base seed, _config/numerics-defaults.toml
const HH_BASE = HH_SEED + 24000           # this slug uses 20261012 + 24000 .. 24999; testset AC<n> uses 50 n + sub (sub < 50)
hh_seed(ac::Integer, sub::Integer = 0) = HH_BASE + 50 * ac + sub

const HH_TWO_PI = 2 * pi

# (d, k) panel of AC1
const HH_DK = ((Val(2), Val(1)), (Val(2), Val(2)), (Val(3), Val(2)), (Val(4), Val(2)), (Val(4), Val(4)), (Val(8), Val(3)), (Val(8), Val(8)))

hh_algs(::Type{T}) where {T} = T === Float64 ? (HouseholderQR(), AugmentedQR(), GivensQR(), LapackQR()) :
                                               (HouseholderQR(), AugmentedQR(), GivensQR())

# Decision D38 (package suite under 3 minutes; the suite is compile-bound, cost = distinct (T, d, k, variant) specializations):
# the package suite keeps one representative per criterion; verify/lyap-householder-qr.jl defines HH_FULL_RUN = true before
# including runtests.jl and so runs the full cross-product of the same testsets (hh_sel_algs, hh_pick, hh_ac1_cases, ...).
const HH_FULL = @isdefined(HH_FULL_RUN) ? HH_FULL_RUN : false
hh_sel_algs(::Type{T}) where {T} = (HH_FULL || T === Float64) ? hh_algs(T) : (HouseholderQR(),)   # Double64: V-H kernel only
hh_one_alg(::Type{T}) where {T} = HH_FULL ? hh_algs(T) : (HouseholderQR(),)                        # V-H kernel only (propagate-level cases)
hh_pick(full, reduced) = HH_FULL ? full : reduced
# (index, element) pairs of `xs` at the indices `idx` (all of them in the full run); the index keeps the seed offsets of the full run
hh_sel(xs, idx) = [(i, xs[i]) for i in (HH_FULL ? eachindex(xs) : idx)]
# AC1 data cases (T index, T, shape index, Val d, Val k) shared by AC1 and AC3
function hh_ac1_cases()
    all = [(it, T, idk, vd, vk) for (it, T) in enumerate((Float64, Double64)) for (idk, (vd, vk)) in enumerate(HH_DK)]
    return HH_FULL ? all : filter(c -> (c[1], c[3]) in ((1, 2), (1, 4), (2, 4)), all)   # Float64 (2,2), (4,2); Double64 (4,2)
end
hh_kernel_algs() = (HouseholderQR(), AugmentedQR(), GivensQR())
hh_name(alg) = string(nameof(typeof(alg)))

# E9 counts: reflections r(d,k), rotations g(d,k)
hh_r(d, k) = k < d ? k : k - 1
hh_g(d, k) = k * d - k * (k + 1) ÷ 2
hh_count(alg, d, k) = hh_name(alg) == "GivensQR" ? hh_g(d, k) : hh_r(d, k)

# ---- number conversions and BigFloat reference algebra (call inside setprecision(BigFloat, b) do ... end)
hh_big(x::Double64) = BigFloat(x.hi) + BigFloat(x.lo)
hh_big(x::Real) = BigFloat(x)
hh_bigmat(M) = [hh_big(M[i, j]) for i in 1:size(M, 1), j in 1:size(M, 2)]
hh_cast(::Type{T}, z) where {T} = SVector{2, T}(T(z[1]), T(z[2]))

# QR of a stored d x k matrix at the current BigFloat precision with the E1 sign convention (positive diagonal).
function hh_ref_qr(W::SMatrix{d, k}) where {d, k}
    F = qr(hh_bigmat(W))
    Q = Matrix(F.Q)[:, 1:k]
    R = Matrix(F.R)[1:k, 1:k]
    for j in 1:k
        if R[j, j] < 0
            Q[:, j] .*= -1
            R[j, :] .*= -1
        end
    end
    return Q, R
end

# determinant by Gaussian elimination with partial pivoting at the current BigFloat precision
function hh_det(A::AbstractMatrix{T}) where {T}
    M = Matrix{T}(A)
    n = size(M, 1)
    d = one(T)
    for c in 1:n
        piv = c
        for i in (c + 1):n
            abs(M[i, c]) > abs(M[piv, c]) && (piv = i)
        end
        if piv != c
            for j in 1:n
                M[c, j], M[piv, j] = M[piv, j], M[c, j]
            end
            d = -d
        end
        d *= M[c, c]
        for i in (c + 1):n
            f = M[i, c] / M[c, c]
            for j in c:n
                M[i, j] -= f * M[c, j]
            end
        end
    end
    return d
end

# matrix power by repeated squaring at the current BigFloat precision
function hh_matpow(A::AbstractMatrix, n::Int)
    R = Matrix{eltype(A)}(I, size(A)...)
    B = copy(A)
    while n > 0
        isodd(n) && (R = R * B)
        B = B * B
        n >>= 1
    end
    return R
end

# ---- seeded synthetic matrices (AC1 to AC3): column j scaled by 10^u, u uniform in [-3, 3]
hh_draw(rng, ::Type{Float64}, s) = s * randn(rng)
hh_draw(rng, ::Type{T}, s) where {T} = T(s * randn(rng)) + T(s * randn(rng) * 2.0^-53)   # low word filled too (Double64)
function hh_randW(rng, ::Type{T}, ::Val{d}, ::Val{k}) where {T, d, k}
    M = Matrix{T}(undef, d, k)
    for j in 1:k
        s = 10.0^(6 * rand(rng) - 3)
        for i in 1:d
            M[i, j] = hh_draw(rng, T, s)
        end
    end
    return SMatrix{d, k, T}(M)
end

# 2-norm condition number of the first i columns of W (kappa_0 = 1), evaluated in Float64 (a scale for a tolerance only)
function hh_kappa(W::SMatrix{d, k}, i::Int) where {d, k}
    i == 0 && return 1.0
    sv = svdvals(Float64[Float64(W[a, b]) for a in 1:d, b in 1:i])
    return sv[1] / sv[end]
end
hh_colnorm(W, j) = sqrt(sum(abs2, Float64(W[a, j]) for a in 1:size(W, 1)))

# ---- tolerance formulas (all in units of the effective eps of the element type)
hh_theta(::Type{T}, G_max, tau, M_max) where {T} = 16 * Float64(effective_eps(T)) * (Float64(G_max) / tau + Float64(M_max))   # GS G13
hh_BS(::Type{T}, n, G_first, K = 0) where {T} = 8 * Float64(effective_eps(T)) * (n + G_first) + 16 * abs(K) * n^2 * Float64(effective_eps(T))   # L1 S3 as amended (K_s = 0 for constant tangents)
# theta_psi(n) = theta_sum + 8 eps Psi_n / n (Amendment D35): the unit of the individual exponents on the standard map; Psi_n of the orbit and the first column
hh_theta_psi(::Type{T}, G_max, tau, M_max, Psi, n) where {T} = hh_theta(T, G_max, tau, M_max) + 8 * Float64(effective_eps(T)) * Psi / n
# Psi_n of the orbit from the stored points zs (z_0 .. z_{n-1}) and the first column of the seeded frame
hh_Psi(sys, zs, W0) = lb_Psi([SMatrix{2, 2, Float64, 4}(map(Float64, tangent(sys, z))) for z in zs], Float64.(W0[:, 1]))

function hh_sigma(K)
    S2 = 2 * K^2 + 2 * abs(K) + 5
    S = sqrt(S2)
    return (S + sqrt(S2 - 4)) / 2
end
hh_Enl(K, d0) = abs(K) * hh_sigma(K) * d0 / sqrt(2)
hh_Ero(K, d0, Mstate, epsT) = (abs(K) + 4) * Mstate * epsT * hh_sigma(K) / d0
hh_D2T(K, d0, Mstate, epsT) = 2 * (hh_Enl(K, d0) + hh_Ero(K, d0, Mstate, epsT))
hh_d0star(K, Mstate, epsT) = sqrt(sqrt(2) * (abs(K) + 4) * Mstate * epsT / abs(K))

# ---- orbit-side quantities along the stored orbit (iterating `forward` in T is bitwise the stream's orbit)
function hh_Mmax(sys, x0, N)
    x, m = x0, 0.0
    for _ in 1:N
        m = max(m, Float64(maximum(abs, tangent(sys, x))))
        x = forward(sys, x)
    end
    return m
end
function hh_Gfirst(sys, x0, tau)          # squared Frobenius norm of the first-window product (L1 S3)
    x, P = x0, one(tangent(sys, x0))
    for _ in 1:tau
        P = tangent(sys, x) * P
        x = forward(sys, x)
    end
    return sum(abs2, Float64(P[i]) for i in 1:length(P))
end

# S7 starts and screening of the regular class at K = 0.5 (lambda_1(1e5) <= 5e-4, from this module's own k = 1 run)
hh_starts(rng, m) = [SVector(HH_TWO_PI * rand(rng), HH_TWO_PI * rand(rng)) for _ in 1:m]
function hh_regular_starts(K, m, rng)
    cand = hh_starts(rng, 3m)
    res = propagate(StandardMap(K), cand, 100_000; k = 1, tau = 1, seed = HH_SEED, reducers = (QRSpectrum(HouseholderQR()),))
    ok = [res.finals[j][1].lambda[1] <= 5e-4 for j in eachindex(cand)]
    @assert count(ok) >= m "fewer than $m regular starts among $(3m) candidates"
    return cand[ok][1:m], count(!, ok[1:(findall(ok)[m])])
end

# ---- test maps and reducers (Batch A names only)
# Standard map on (x, p) plus two extra components z[3] (step counter) and z[4] (sum of x): the W9 map of tangent-core, s = 4, d = 2.
struct hh_Long{T} <: AbstractMap
    K::T
end
struct hh_Short{T} <: AbstractMap      # the same map on the plain length-2 state
    K::T
end
Chaoslab.dimension(::Union{hh_Long, hh_Short}) = 2
function Chaoslab.forward(m::hh_Long, z)
    p2 = z[2] + m.K * sin(z[1])
    return SVector(z[1] + p2, p2, z[3] + 1, z[4] + z[1])
end
function Chaoslab.forward(m::hh_Short, z)
    p2 = z[2] + m.K * sin(z[1])
    return SVector(z[1] + p2, p2)
end
function Chaoslab.tangent(m::Union{hh_Long, hh_Short}, z)
    c = m.K * cos(z[1])
    return SMatrix{2, 2}(1 + c, c, one(c), one(c))
end

# Fault map on state (t, flag), d = 2: the tangent is the identity except `M` at the step whose pre-step counter is `at - 1`
# on a flagged orbit (flag > 0.5). Dyadic data, so the degenerate matrices below have an exactly zero pivot.
struct hh_Degen <: AbstractMap
    at::Int
    M::SMatrix{2, 2, Float64, 4}
end
Chaoslab.dimension(::hh_Degen) = 2
Chaoslab.forward(::hh_Degen, z) = SVector(z[1] + 1, z[2])
Chaoslab.tangent(m::hh_Degen, z) = (z[2] > 0.5 && z[1] == m.at - 1) ? m.M : SMatrix{2, 2, Float64}(1, 0, 0, 1)

# Read-only reducer: largest |W - I| seen and the number of handoffs (placed after the estimator it sees Q).
struct hh_Dev end
Chaoslab.reducer_init(::hh_Dev, ::Type{T}, ::Val, ::Val) where {T} = (zero(T), 0)
function Chaoslab.reducer_update(::hh_Dev, st, n::Int, x, W::SMatrix{d, k, T}) where {d, k, T}
    dev = zero(T)
    for j in 1:k, i in 1:d
        dev = max(dev, abs(W[i, j] - (i == j ? one(T) : zero(T))))
    end
    return (max(st[1], dev), st[2] + 1), W
end
Chaoslab.reducer_final(::hh_Dev, st) = (maxdev = st[1], count = st[2])

# ---- device plumbing (copy of the tangent-core pattern): the run-wide device is set and restored explicitly
const hh_F = Chaoslab.Foundation
const hh_Tg = Chaoslab.Tangent
struct hh_FakeDevice <: hh_F.AbstractDevice
    loadable::Bool
    functional::Bool
end
hh_F.probe(d::hh_FakeDevice, ::Type) = (loadable = d.loadable, functional = d.functional, device = "FakeGPU", compiles = true, gpu_faster = true)
hh_Tg.device_backend(::hh_FakeDevice) = hh_Tg.KernelAbstractions.CPU()
function hh_with_device(f, dev)
    old = hh_F.DEVICE[]
    hh_F.DEVICE[] = dev
    try
        return f()
    finally
        hh_F.DEVICE[] = old
    end
end

# ---- misc
function hh_health_error(f)
    try
        f()
        return nothing
    catch e
        e isa NumericalHealthError ? e : rethrow()
    end
end
# order-free equality of two NamedTuples (the contract fixes keys and values, not their order)
hh_nt_eq(a::NamedTuple, b::NamedTuple) = Set(keys(a)) == Set(keys(b)) && all(k -> getproperty(a, k) == getproperty(b, k), keys(a))
