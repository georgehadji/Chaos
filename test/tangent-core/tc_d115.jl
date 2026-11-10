# Amendment D115 (frame product with explicit fma): criteria D115a to D115e. D115f is the whole suite, D115g the CI matrix, D115h stage 04.
# Spec: tangent-core-spec.md, "Amendment D115". Interface: Chaoslab.Foundation.fused_product(A, B) (not exported; does not exist before stage 03).
# RED: UndefVarError for fused_product. Seeds: TC_SEED + 11500 + sub.

const TC115_SEED = TC_SEED + 11500
tc_fp(A, B) = Chaoslab.Foundation.fused_product(A, B)

# Scalar reference of item 1: first product, then fma in ascending inner index.
function tc115_ref(A::AbstractMatrix{T}, B::AbstractMatrix{T}) where {T}
    m, p = size(A)
    n = size(B, 2)
    C = Matrix{T}(undef, m, n)
    for j in 1:n, i in 1:m
        acc = A[i, 1] * B[1, j]
        for q in 2:p
            acc = fma(A[i, q], B[q, j], acc)
        end
        C[i, j] = acc
    end
    return C
end

@testset "D115a definition" begin
    rng = StableRNG(TC115_SEED + 1)
    for T in (Float64, Float32), m in (2, 4, 6, 8), n in (1, 2, 3, 4, 8), _ in 1:20
        A = SMatrix{m, m, T}(randn(rng, T, m, m))
        B = SMatrix{m, n, T}(randn(rng, T, m, n) .* T(10) .^ rand(rng, -3:3, m, n))
        C = tc_fp(A, B)
        @test C isa SMatrix{m, n, T}
        @test isequal(Matrix(C), tc115_ref(Matrix(A), Matrix(B)))
    end
end

@testset "D115b column independence" begin
    rng = StableRNG(TC115_SEED + 2)
    for T in (Float64, Float32), m in (2, 4, 6), n in (2, 3, 4), _ in 1:20
        A = SMatrix{m, m, T}(randn(rng, T, m, m))
        B = SMatrix{m, n, T}(randn(rng, T, m, n))
        C = tc_fp(A, B)
        for c in 1:n
            @test isequal(C[:, SOneTo(c)], tc_fp(A, B[:, SOneTo(c)]))
            @test isequal(C[:, SVector(c)], tc_fp(A, B[:, SVector(c)]))
        end
    end
end

@testset "D115c other element types" begin
    rng = StableRNG(TC115_SEED + 3)
    for T in (Double64, BigFloat), _ in 1:5
        A = SMatrix{4, 4, T}(T.(randn(rng, 4, 4)))
        B = SMatrix{4, 3, T}(T.(randn(rng, 4, 3)))
        @test isequal(tc_fp(A, B), A * B)
    end
end

# Reference stream of item 2: no reducer, the frame updated by fused_product before the point advances.
function tc115_loop(sys, x, W, N)
    for _ in 1:N
        W = tc_fp(tangent(sys, x), W)
        x = forward(sys, x)
    end
    return x, W
end

@testset "D115d stream" begin
    rng = StableRNG(TC115_SEED + 4)
    N = 500
    std = StandardMap(6.0)
    fro = FroeschleMap(0.5, 0.1, 0.001)
    cases = [(std, SVector(0.1, 0.2), (1, 2)),
             (fro, SVector(0.5, 1.0, 1.5, 2.0), (2, 3, 4)),
             (PairMap(std), SVector(0.1, 0.2, 0.1 + 1e-8, 0.2), (1, 2))]
    for (sys, x0, ks) in cases
        d = dimension(sys)
        W0 = SMatrix{d, maximum(ks)}(randn(rng, d, maximum(ks)))
        full = nothing
        for k in ks
            Wk = W0[:, SOneTo(k)]
            x, W, _ = propagate_orbit(sys, x0, Wk, N, 1, ())
            xr, Wr = tc115_loop(sys, x0, Wk, N)
            @test isequal(x, xr)
            @test isequal(W, Wr)
            k == maximum(ks) && (full = W)
        end
        for k in ks
            x, W, _ = propagate_orbit(sys, x0, W0[:, SOneTo(k)], N, 1, ())
            @test isequal(W, full[:, SOneTo(k)])
        end
    end
end

@testset "D115e no allocation" begin
    forced = Base.JLOptions().check_bounds == 1
    A, B = @SMatrix(rand(4, 4)), @SMatrix(rand(4, 4))
    f(A, B) = @allocated tc_fp(A, B)
    f(A, B)
    if forced
        @test_skip f(A, B) == 0
    else
        @test f(A, B) == 0
    end
end
