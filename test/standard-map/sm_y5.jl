# Y5 linear symplectic maps: criteria 15 to 17 (spec item 10). Blocks: rotation, hyperbolic H, reflected -H.

# Random block parameters: kind 1 rotation (theta in [0, 2 pi)), 2 hyperbolic, 3 reflected hyperbolic (mu in [0.05, 0.5]).
function sm_y5_params(rng, N::Integer)
    return [(rand(rng, 1:3), 2 * Float64(pi) * rand(rng), 0.05 + 0.45 * rand(rng)) for _ in 1:N]
end

function sm_y5_block(::Type{T}, (kind, theta, mu)) where {T}
    kind == 1 && return rotation_block(T, T(theta))
    H = hyperbolic_block(T, T(mu))
    return kind == 2 ? H : -H
end

# block of the inverse: R(-theta) = R'; H(-mu) = H^-1; (-H)^-1 = -H(-mu)
function sm_y5_inv_block(::Type{T}, (kind, theta, mu)) where {T}
    kind == 1 && return rotation_block(T, T(-theta))
    H = hyperbolic_block(T, T(-mu))
    return kind == 2 ? H : -H
end

# exact spectrum lambda = ln|eigenvalue|, sorted descending: +-mu per hyperbolic block, 0 twice per rotation block
function sm_y5_spectrum(params)
    lam = Float64[]
    for (kind, _, mu) in params
        append!(lam, kind == 1 ? (0.0, 0.0) : (mu, -mu))
    end
    return sort!(lam; rev = true)
end

sm_y5_A(::Type{T}, params, ::Val{N}) where {T, N} =
    block_diagonal(ntuple(j -> sm_y5_block(T, params[j]), Val(N))...)
sm_y5_Ainv(::Type{T}, params, ::Val{N}) where {T, N} =
    block_diagonal(ntuple(j -> sm_y5_inv_block(T, params[j]), Val(N))...)

function sm_y5_B(::Type{T}, rng, ::Val{N}) where {T, N}
    M = [2 * rand(rng) - 1 for _ in 1:N, _ in 1:N]
    return SMatrix{N, N, T}(T.((M + M') / 2))
end

# S = [I B; 0 I] and its inverse [I -B; 0 I] as 2N x 2N static matrices
function sm_y5_S(::Type{T}, B, ::Val{N}) where {T, N}
    Z = zeros(T, N, N)
    Id = Matrix{T}(I, N, N)
    return SMatrix{2N, 2N, T}([Id Matrix(B); Z Id]), SMatrix{2N, 2N, T}([Id -Matrix(B); Z Id])
end

# The tests do their own matrix arithmetic on dense Matrix copies (cheap to compile); Chaoslab functions get the SMatrix.
# Frobenius norm is used for ||A|| (spec leaves the norm open; it dominates the 2-norm).
function sm_ac15_ratios(::Type{T}, rng, ::Val{N}) where {T, N}
    params = sm_y5_params(rng, N)
    A = sm_y5_A(T, params, Val(N))
    AS = conjugated(A, sm_y5_B(T, rng, Val(N)))
    Om = symplectic_form(T, 2N)
    Omm = Matrix(Om)
    worst = 0.0
    for M in (A, AS)
        Mm = Matrix(M)
        nrm = Float64(sm_fro(Mm))
        worst = max(worst, Float64(sm_maxabs(sm_mm(sm_mm(sm_tr(Mm), Omm), Mm) - Omm)) / (16 * sm_eps(T) * nrm^2))
        worst = max(worst, Float64(symplectic_residual(M)) / (16 * sm_eps(T) * nrm^2))
        worst = max(worst, Float64(abs(sm_det(Mm) - 1)) / (16 * sm_eps(T) * nrm^N))
        worst = max(worst, Float64(det_residual(M)) / (16 * sm_eps(T) * nrm^N))
    end
    # map-level contract of item 10: dimension 2N, tangent equal to A at every state, forward z -> A z
    sys = LinearSymplecticMap(A)
    z = SVector{2N, T}(T.([2 * rand(rng) - 1 for _ in 1:(2N)]))
    Am = Matrix(A)
    ok = (dimension(sys) == 2N) && (Matrix(tangent(sys, z)) == Am) &&
        (sm_maxabs(Vector(forward(sys, z)) - vec(sm_mm(Am, reshape(Vector(z), :, 1)))) <= 16 * sm_eps(T) * sm_fro(Am) * sm_fro(Vector(z)))
    omega_ok = Omm == [zeros(T, N, N) Matrix{T}(I, N, N); -Matrix{T}(I, N, N) zeros(T, N, N)]
    return (worst, ok, omega_ok)
end

# sigma_max(M) in Float64 from the SVD
sm_smax(M) = svdvals(Matrix(M))[1]

function sm_ac16_ratios(rng, ::Val{N}) where {N}
    params = sm_y5_params(rng, N)
    lam = sm_y5_spectrum(params)
    A = sm_y5_A(Float64, params, Val(N))
    B = sm_y5_B(Float64, rng, Val(N))
    AS = conjugated(A, B)
    S, Sinv = sm_y5_S(Float64, B, Val(N))
    kappa = opnorm(Matrix(S)) * opnorm(Matrix(Sinv))
    worst_n = 0.0      # block-diagonal A, bound 16 eps
    worst_s = 0.0      # A_S, bound 16 eps + ln(kappa)/n
    for (M, bound_extra) in ((A, 0.0), (AS, log(kappa))), n in (10, 50)
        Mm = Matrix(M)
        Mi = Matrix(LinearSymplecticMap(M).Ainv)
        tol = 16 * sm_eps(Float64) + bound_extra / n
        r1 = abs(log(sm_smax(Mm^n)) / n - lam[1]) / tol
        r2 = abs(log(sm_smax(Mi^n)) / n + lam[end]) / tol
        if bound_extra == 0.0
            worst_n = max(worst_n, r1, r2)
        else
            worst_s = max(worst_s, r1, r2)
        end
    end
    return (worst_n, worst_s)
end

function sm_ac17_ratios(::Type{T}, rng, ::Val{N}) where {T, N}
    params = sm_y5_params(rng, N)
    A = sm_y5_A(T, params, Val(N))
    Ai = LinearSymplecticMap(A).Ainv              # Omega^-1 A' Omega
    Ab = sm_y5_Ainv(T, params, Val(N))            # inverse of the block-diagonal form
    Am, Aim = Matrix(A), Matrix(Ai)
    tol = 16 * sm_eps(T) * sm_fro(Am)^2
    Id = Matrix{T}(I, 2N, 2N)
    z = SVector{2N, T}(T.([2 * rand(rng) - 1 for _ in 1:(2N)]))
    sys = LinearSymplecticMap(A)
    r = max(
        Float64(sm_maxabs(sm_mm(Aim, Am) - Id)),
        Float64(sm_maxabs(sm_mm(Am, Aim) - Id)),
        Float64(sm_maxabs(Aim - Matrix(Ab))),
        Float64(sm_maxabs(Vector(inverse(sys, forward(sys, z))) - Vector(z))),
        Float64(sm_maxabs(Matrix(inverse_tangent(sys, z)) - Aim)),
    )
    return r / Float64(tol)
end

@testset "AC15 Y5 symplecticity" begin
    # max|A' Omega A - Omega| <= 16 eps ||A||^2 and |det A - 1| <= 16 eps ||A||^N for A and A_S; 10 draws per N in {1,2,3,5}
    rng = StableRNG(SM_SEED + 15001)
    for T in (Float64, Double64), N in (1, 2, 3, 5), _ in 1:10
        worst, ok, omega_ok = sm_ac15_ratios(T, rng, Val(N))
        @test worst <= 1
        @test ok
        @test omega_ok
    end
end

@testset "AC16 Y5 exact exponents" begin
    # |(1/n) ln smax(A^n) - lambda_1| <= 16 eps and |(1/n) ln smax(A^-n) + lambda_2N| <= 16 eps, n in {10, 50}; A_S adds ln(kappa(S))/n
    # Float64 only: sigma_max needs an SVD, which is not available for Double64 in the test environment
    rng = StableRNG(SM_SEED + 16001)
    for N in (1, 2, 3, 5), _ in 1:10
        worst_n, worst_s = sm_ac16_ratios(rng, Val(N))
        @test worst_n <= 1
        @test worst_s <= 1
    end
end

@testset "AC17 Y5 inverse" begin
    # A^-1 A = A A^-1 = I and Omega^-1 A' Omega = block inverse, max entry <= 16 eps ||A||^2; 10 draws per N
    rng = StableRNG(SM_SEED + 17001)
    for T in (Float64, Double64), N in (1, 2, 3, 5), _ in 1:10
        @test sm_ac17_ratios(T, rng, Val(N)) <= 1
    end
end
