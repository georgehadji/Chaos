# Shared helpers of the lyap-gram-schmidt suite (prefix gs_); lyap-svd-stabilized includes this file too.
# Lyapunov names (CGS, QRSpectrum, ...) appear only inside function bodies and testsets, so a missing
# Chaoslab.Lyapunov makes testsets error instead of breaking the include. Spec: lyap-gram-schmidt-spec.md.
using Test, Chaoslab, StaticArrays, LinearAlgebra, StableRNGs, DoubleFloats

const GS_HELPERS_LOADED = true            # include guard, tested by both runtests.jl files
const GS_SEED = 20261012                 # base seed, _config/numerics-defaults.toml
gs_seed(ac) = GS_SEED + 20000 + 100 * ac # lyap-gram-schmidt offsets 20100..21400
sv_seed(ac) = GS_SEED + 30000 + 100 * ac # lyap-svd-stabilized offsets 30100..31200

gs_algs() = (CGS(), MGS(), RGS())
gs_name(alg) = string(nameof(typeof(alg)))
gs_types() = (Float64, Double64)
# Double64 runs use RGS only (the compiled specialisations dominate the suite time; the variants are told apart in Float64).
gs_algsT(::Type{Float64}) = gs_algs()
gs_algsT(::Type{Double64}) = (RGS(),)
gs_eps(::Type{T}) where {T} = Float64(effective_eps(T))   # D34: the effective eps of the element type
# lyap-benettin S3 as amended (Amendment D35): B_S(n) = 8 eps (n + G_first) + 16 |K_s| n^2 eps, K_s = |K| for a standard-map orbit; written once in L1 S3
gs_BS(::Type{T}, n, G_first, K) where {T} = 8 * gs_eps(T) * (n + G_first) + 16 * abs(Float64(K isa Symbol ? 2pi : K)) * n^2 * gs_eps(T)

# Exact conversions.
gs_big(x::Double64) = BigFloat(x.hi) + BigFloat(x.lo)
gs_big(x::Real) = BigFloat(x)
gs_bigmat(M) = [gs_big(M[i, j]) for i in 1:size(M, 1), j in 1:size(M, 2)]
gs_rat(x::Double64) = Rational{BigInt}(x.hi) + Rational{BigInt}(x.lo)
gs_rat(x::Real) = Rational{BigInt}(x)
gs_absdiff(a, b::BigFloat) = Float64(abs(gs_big(a) - b))   # |a - b| as Float64, difference taken at 256 bits

# Standard-map panel value: a Float64 literal or :twopi (2 pi built in T, never through Float64).
gs_K(::Type{T}, K) where {T} = K === :twopi ? 2 * exact(T, pi) : T(K)
gs_KS() = (0.5, 2.0, 5.0, :twopi)
gs_stdmap(::Type{T}, K) where {T} = StandardMap(gs_K(T, K))

# Seeded starts: x, p uniform in [0, 2 pi) drawn in Float64, widened exactly (lyap-benettin S2, S7).
function gs_starts(::Type{T}, seed, m) where {T}
    rng = StableRNG(seed)
    return [SVector{2,T}(T(2 * Float64(pi) * rand(rng)), T(2 * Float64(pi) * rand(rng))) for _ in 1:m]
end

# Orbit points x_0 .. x_{N-1} (the points before each step) and the largest one-step tangent entry over them.
function gs_points(sys, x0, N)
    pts = Vector{typeof(x0)}(undef, N)
    x = x0
    for i in 1:N
        pts[i] = x
        x = forward(sys, x)
    end
    return pts
end
gs_Mmax(sys, pts) = maximum(Float64(maximum(abs, tangent(sys, x))) for x in pts)

# Tolerance unit theta_sum (G13): pairing_residual_tol (16) x eps_eff x (G_max / tau + M_max), as Float64.
gs_theta(::Type{T}, G_max, tau, M_max) where {T} = 16 * gs_eps(T) * (Float64(G_max) / tau + Float64(M_max))
# Linear maps: the tangent is the constant matrix A and every window has G = cond_2(A^tau) exactly (R_j has the
# singular values of A^tau Q_{j-1}, Q orthogonal); the estimator's G_max is only an upper bound for k >= 3 (API Q5).
gs_theta_lin(::Type{T}, A, tau) where {T} =
    gs_theta(T, cond(Matrix{Float64}(Float64.(A^tau))), tau, maximum(abs, Float64.(A)))

# Orthogonality bound of a recorded maximum: 8 k eps (RGS), 8 k G_w eps with G_w <= G + 1 (CGS, MGS at d = k = 2).
function gs_orth_bound(alg, ::Type{T}, k, G_max) where {T}
    g = alg isa RGS ? 1.0 : Float64(G_max) + 1
    return 8 * k * g * gs_eps(T)
end

# Linear maps of the Y5 family (z_0 = 0 is a fixed point, so x_n = 0 and the tangent is constant).
gs_cat(::Type{T}) where {T} = LinearSymplecticMap(SMatrix{2,2,T}(2, 1, 1, 1))
function gs_AS(::Type{T}) where {T}    # G12c: S^-1 diag(A_1, A_2) S, integer entries, det 1
    A1 = SMatrix{2,2,T}(2, 1, 1, 1)
    A2 = SMatrix{2,2,T}(3, 1, 2, 1)
    return conjugated(block_diagonal(A1, A2), SMatrix{2,2,T}(1, 2, 2, 1))
end
gs_zero(::Type{T}, d) where {T} = zero(SVector{d,T})

# Fibonacci numbers (Int64 exact up to F_92) and A^n = [[F_{2n+1}, F_{2n}], [F_{2n}, F_{2n-1}]] for the cat map.
function gs_fib(n::Int)
    a, b = 0, 1
    for _ in 1:n
        a, b = b, a + b
    end
    return a
end
# Integer entries above 2^53 need both parts: Double64(::Int64) rounds through Float64 (hi + lo split here keeps them exact).
gs_exactint(::Type{T}, n::Integer) where {T} = T(Float64(n)) + T(Float64(n - round(Int, Float64(n))))
gs_Apow(::Type{T}, n::Int) where {T} =
    SMatrix{2,2,T}(gs_exactint(T, gs_fib(2n + 1)), gs_exactint(T, gs_fib(2n)), gs_exactint(T, gs_fib(2n)), gs_exactint(T, gs_fib(2n - 1)))

# G12a: lambda_1(n) = ln mu + (1/(2n)) ln(c_+^2 + c_-^2 mu^-4n) for the first frame column q1; BigFloat(256).
function gs_cat_lambda1(q1, n)
    return setprecision(BigFloat, 256) do
        s5 = sqrt(BigFloat(5))
        mu = (3 + s5) / 2
        g = (s5 - 1) / 2                      # e_+ = (1, g) / nrm, e_- = (-g, 1) / nrm
        nrm = sqrt(1 + g^2)
        q = (gs_big(q1[1]), gs_big(q1[2]))
        cp = (q[1] + g * q[2]) / nrm
        cm = (-g * q[1] + q[2]) / nrm
        log(mu) + log(cp^2 + cm^2 * mu^(-4n)) / (2n)
    end
end
gs_lnmu() = setprecision(BigFloat, 256) do
    log((3 + sqrt(BigFloat(5))) / 2)
end
# Eigenbasis frame (e_+, e_-) of the cat map built in T.
function gs_cat_eigframe(::Type{T}) where {T}
    g = (sqrt(T(5)) - 1) / 2
    nrm = sqrt(1 + g^2)
    return SMatrix{2,2,T}(1 / nrm, g / nrm, -g / nrm, 1 / nrm)
end

# 256-bit chain on the stored orbit points of a standard map (AC4): tangent at 256 bits at the SAME points, the same
# widened frame, a reorthogonalised Gram-Schmidt (RGS) chain; returns S_i(n) as BigFloat. Never a separate 256-bit orbit.
# Scalar arithmetic on the entries of the 2 x 2 matrices (generic BigFloat matrix products are slow).
function gs_ref_std_S(Kt, pts, Q0, n)
    return setprecision(BigFloat, 256) do
        Kb = gs_big(Kt)
        a1, b1, a2, b2 = gs_big(Q0[1, 1]), gs_big(Q0[2, 1]), gs_big(Q0[1, 2]), gs_big(Q0[2, 2])   # columns (a1, b1), (a2, b2)
        P1, P2 = one(BigFloat), one(BigFloat)       # products of the pivots; one logarithm at the end (log is slow at 256 bits)
        for i in 1:n
            c = Kb * cos(gs_big(pts[i][1]))
            w1 = ((1 + c) * a1 + b1, c * a1 + b1)      # M = [[1 + c, 1], [c, 1]]
            v = ((1 + c) * a2 + b2, c * a2 + b2)
            r11 = sqrt(w1[1]^2 + w1[2]^2)
            q1 = (w1[1] / r11, w1[2] / r11)
            r12 = q1[1] * v[1] + q1[2] * v[2]
            v = (v[1] - r12 * q1[1], v[2] - r12 * q1[2])
            s = q1[1] * v[1] + q1[2] * v[2]            # second pass
            v = (v[1] - s * q1[1], v[2] - s * q1[2])
            r22 = sqrt(v[1]^2 + v[2]^2)
            a1, b1, a2, b2 = q1[1], q1[2], v[1] / r22, v[2] / r22
            P1 *= r11
            P2 *= r22
        end
        (log(P1), log(P2))
    end
end

# 256-bit product of the tangents of a standard map at the stored points; returns Y_n = [y11 y12; y21 y22] as a
# 2 x 2 Matrix{BigFloat} at the requested n values (a Dict n => Y).
function gs_ref_std_Y(Kt, pts, ns)
    return setprecision(BigFloat, 256) do
        Kb = gs_big(Kt)
        y11, y12, y21, y22 = one(BigFloat), zero(BigFloat), zero(BigFloat), one(BigFloat)
        out = Dict{Int,Matrix{BigFloat}}()
        for i in 1:maximum(ns)
            c = Kb * cos(gs_big(pts[i][1]))
            y11, y12, y21, y22 = (1 + c) * y11 + y21, (1 + c) * y12 + y22, c * y11 + y21, c * y12 + y22
            i in ns && (out[i] = [y11 y12; y21 y22])
        end
        out
    end
end

# Serial driver: the per-orbit body of propagate (propagate_orbit, bitwise the same result, AC9) over all orbits, with
# propagate's keywords and result fields (x, W, finals). Used by the sweeps to keep the number of compiled kernel
# specialisations low; tests of propagate itself (AC8, AC9, AC14) call propagate.
function gs_serial(sys, x0s, N; k, tau, seed, reducers, W0 = nothing)
    T = eltype(first(x0s))
    d = dimension(sys)
    Wi(j) = W0 === nothing ? initial_deviations(seed, j, Val(d), Val(k), T) : W0[j]
    out = [propagate_orbit(sys, x0s[j], Wi(j), N, tau, reducers; orbit = j) for j in eachindex(x0s)]
    return (x = [o[1] for o in out], W = [o[2] for o in out], finals = [map(reducer_final, reducers, o[3]) for o in out])
end

# Run the estimator on all orbits; `tol` is the orthogonality_tol multiplier (Inf measures the loss without raising).
function gs_run(sys, x0s, N; k, tau, seed, alg, tol = 8.0, W0 = nothing)
    return gs_serial(sys, x0s, N; k, tau, seed, reducers = (QRSpectrum(alg; orthogonality_tol = tol),), W0)
end
gs_fin(res, i) = res.finals[i][1]
# gs_run that returns (result, nothing) or (nothing, the NumericalHealthError it raised)
function gs_try_run(args...; kw...)
    try
        return gs_run(args...; kw...), nothing
    catch e
        e isa NumericalHealthError ? (nothing, e) : rethrow()
    end
end
gs_S(f) = f.sum .+ f.comp     # S_i(n) = sum + compensation

# Run f(); return the NumericalHealthError it throws (nothing if it returns), rethrow anything else.
function gs_health_error(f)
    try
        f()
        return nothing
    catch e
        e isa NumericalHealthError ? e : rethrow()
    end
end

# Bytes allocated by one propagate_orbit call after a warm-up call (measured inside a function).
function gs_alloc_bytes(sys, x0, W0, N, tau, reducers)
    propagate_orbit(sys, x0, W0, N, tau, reducers)
    return @allocated propagate_orbit(sys, x0, W0, N, tau, reducers)
end

# Bytes allocated by one orthonormalize call after a warm-up call (inside a function, not at testset scope).
function gs_alloc_orthonormalize(alg, W)
    orthonormalize(alg, W)
    return @allocated orthonormalize(alg, W)
end

# Plain per-orbit loop of a map and its tangent that drives the reducer by hand at each handoff (tangent-core W9 style);
# returns (x, W, state). Used where "bitwise equal to a plain per-orbit loop" is the criterion.
function gs_manual_loop(r, sys, x0::SVector{s,T}, W0::SMatrix{d,k,T}, N, tau) where {s,d,k,T}
    st = reducer_init(r, T, Val(d), Val(k))
    x, W = x0, W0
    for j in 1:(N ÷ tau)
        for _ in 1:tau
            W = tangent(sys, x) * W
            x = forward(sys, x)
        end
        st, W = reducer_step(r, st, j * tau, x, W)
    end
    return x, W, st
end

# Exact volume reference (AC11): ln vol_p(Y_n X) = (1/2) ln det(X' Y_n' Y_n X) with Y_n = A^n in BigInt, X the exact
# rationals of the frame entries, the Gram determinant exact, only the logarithm taken at 256 bits. Plain loops on
# BigInt and Rational{BigInt}: the generic LinearAlgebra paths for these types cost seconds of compilation.
function gs_bigint_matmul(A::Matrix{BigInt}, B::Matrix{BigInt})
    C = zeros(BigInt, size(A, 1), size(B, 2))
    for i in axes(A, 1), j in axes(B, 2), l in axes(A, 2)
        C[i, j] += A[i, l] * B[l, j]
    end
    return C
end
function gs_bigint_power(A::Matrix{BigInt}, n::Int)
    R = Matrix{BigInt}(I(size(A, 1)))
    B = copy(A)
    while n > 0
        isodd(n) && (R = gs_bigint_matmul(R, B))
        B = gs_bigint_matmul(B, B)
        n >>= 1
    end
    return R
end
function gs_exact_logvol(Aint::Matrix{BigInt}, n, Q0, p)
    Y = gs_bigint_power(Aint, n)
    d = size(Y, 1)
    V = [sum(Y[i, l] * gs_rat(Q0[l, j]) for l in 1:d) for i in 1:d, j in 1:p]          # Y_n X, exact
    G = Rational{BigInt}[sum(V[l, i] * V[l, j] for l in 1:d) for i in 1:p, j in 1:p]    # Gram matrix, exact, positive definite
    det = one(Rational{BigInt})
    for c in 1:p                                                                         # elimination without pivoting
        det *= G[c, c]
        for r in (c + 1):p
            f = G[r, c] / G[c, c]
            for k in c:p
                G[r, k] -= f * G[c, k]
            end
        end
    end
    return setprecision(BigFloat, 256) do
        log(BigFloat(det)) / 2
    end
end
gs_intmat(A) = BigInt.(round.(Int, Float64.(Matrix(A))))
