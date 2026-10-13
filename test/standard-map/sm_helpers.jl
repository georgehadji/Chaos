# Helpers for the standard-map suite (prefix sm_). Definitions only: nothing here calls Chaoslab at load time.

const SM_SEED = 20261012                       # _config/numerics-defaults.toml, seed
const SM_KS = [0.0, 0.5, -0.5, 1.0, -1.0, 2.0, -2.0, 4.0, -4.0, 5.0, -5.0, 6.0, -6.0, 10.0, -10.0, 20.0, -20.0]
const SM_BITS = 256                            # numerics-defaults.toml, bigfloat_bits
const SM_EPS_DD_FLOOR = ldexp(1.0, -104)       # nominal double-double spacing, spec "Equations" notation
const SM_EPS_DD = ldexp(1.0, -100)             # accepted effective double-double eps (criterion 14 pass limit); DoubleFloats sin/cos measure 7.2e-31

sm_eps(::Type{Float64}) = ldexp(1.0, -52)
sm_eps(::Type{Double64}) = SM_EPS_DD    # nominal 2^-104 is below what DoubleFloats sin/cos and products reach; see test plan

# Invariants: sigma(K) bounds the largest singular value of DT; horizon n_H(K) = largest n <= 20 with sigma^n <= 1e6.
function sm_sigma(K::Real)
    S2 = 2 * K^2 + 2 * abs(K) + 5
    S = sqrt(S2)
    return (S + sqrt(S2 - 4)) / 2
end

function sm_nH(K::Real)
    s = sm_sigma(K)
    n = 0
    while n < 20 && s^(n + 1) <= 1.0e6
        n += 1
    end
    return max(n, 1)
end

# M_n = max(1, |K|, max_{k<=n} |x_k|, |p_k|) on the (lifted) orbit zs[1:n+1]
function sm_Mn(K::Real, zs, n::Integer)
    m = max(1.0, abs(Float64(K)))
    for k in 1:(n + 1), c in zs[k]
        m = max(m, abs(Float64(c)))
    end
    return m
end

sm_tau(epsT::Real, K::Real, zs, n::Integer) = 16 * epsT * sm_Mn(K, zs, n) * n * sm_sigma(K)^n

# Seeded cases (K, z0): m points per listed K (and m points with K uniform in [-20, 20]); z0 uniform in [-4 pi, 4 pi]^2.
function sm_cases(rng, m::Integer; listed::Bool = true, uniform::Bool = true)
    lim = 4 * Float64(pi)
    cases = Tuple{Float64, SVector{2, Float64}}[]
    if listed
        for K in SM_KS, _ in 1:m
            x = (2 * rand(rng) - 1) * lim
            p = (2 * rand(rng) - 1) * lim
            push!(cases, (K, SVector(x, p)))
        end
    end
    if uniform
        for _ in 1:m
            K = (2 * rand(rng) - 1) * 20.0
            x = (2 * rand(rng) - 1) * lim
            p = (2 * rand(rng) - 1) * lim
            push!(cases, (K, SVector(x, p)))
        end
    end
    return cases
end

# Exact widening of Float64 / Double64 / BigFloat to BigFloat (call inside setprecision(BigFloat, SM_BITS)).
sm_big(x::Double64) = BigFloat(x.hi) + BigFloat(x.lo)
sm_big(x::AbstractFloat) = BigFloat(x)

sm_cast(::Type{T}, z) where {T} = SVector{2, T}(T(z[1]), T(z[2]))

# Plain (lifted) max-norm difference, evaluated in 256-bit arithmetic and returned as Float64.
function sm_absdiff(a, b)
    return setprecision(BigFloat, SM_BITS) do
        Float64(maximum(abs(sm_big(a[i]) - sm_big(b[i])) for i in 1:length(a)))
    end
end

# Torus distance: max over components of |a - b| reduced to (-pi, pi], evaluated in 256-bit arithmetic.
function sm_tdist(a, b)
    return setprecision(BigFloat, SM_BITS) do
        tp = 2 * BigFloat(pi)
        m = 0.0
        for i in 1:length(a)
            d = sm_big(a[i]) - sm_big(b[i])
            d -= tp * round(d / tp)
            m = max(m, Float64(abs(d)))
        end
        m
    end
end

# n+1 states of the orbit under f (default forward; pass f = inverse for the time-reversed orbit)
function sm_orbit(sys, z, n::Integer; f = forward)
    zs = [z]
    for _ in 1:n
        push!(zs, f(sys, zs[end]))
    end
    return zs
end

# Independent plain-formula reference step (not Chaoslab): p' = p + K sin x, x' = x + p'
function sm_ref_orbit(K, z, n::Integer)
    zs = [z]
    for _ in 1:n
        x, p = zs[end]
        p2 = p + K * sin(x)
        push!(zs, SVector(x + p2, p2))
    end
    return zs
end

# Fixed-point class and magnitude from a trace (det = 1): hyperbolic -> |lambda_+|, elliptic -> rotation angle.
function sm_class(tr)
    if abs(tr) > 2
        return (:hyperbolic, (abs(tr) + sqrt(tr^2 - 4)) / 2)
    elseif abs(tr) < 2
        return (:elliptic, acos(tr / 2))
    else
        return (:parabolic, one(tr))
    end
end

sm_omega(::Type{T}) where {T} = SMatrix{2, 2, T}(zero(T), -one(T), one(T), zero(T))   # [[0, 1], [-1, 0]]

# Plain-loop dense linear algebra on Matrix{T} (compiles once per element type; the tests use it for Y5 instead of generic LinearAlgebra)
function sm_mm(A::AbstractMatrix{T}, B::AbstractMatrix{T}) where {T}
    C = zeros(T, size(A, 1), size(B, 2))
    for j in axes(B, 2), k in axes(A, 2), i in axes(A, 1)
        C[i, j] += A[i, k] * B[k, j]
    end
    return C
end

sm_tr(A::AbstractMatrix{T}) where {T} = permutedims(A)
sm_maxabs(A) = maximum(abs, A)
sm_fro(A) = sqrt(sum(abs2, A))

# determinant by Gaussian elimination with partial pivoting
function sm_det(A::AbstractMatrix{T}) where {T}
    M = Matrix{T}(A)
    n = size(M, 1)
    d = one(T)
    for k in 1:n
        piv = k
        for i in (k + 1):n
            abs(M[i, k]) > abs(M[piv, k]) && (piv = i)
        end
        if piv != k
            for j in 1:n
                M[k, j], M[piv, j] = M[piv, j], M[k, j]
            end
            d = -d
        end
        d *= M[k, k]
        for i in (k + 1):n
            f = M[i, k] / M[k, k]
            for j in k:n
                M[i, j] -= f * M[k, j]
            end
        end
    end
    return d
end
