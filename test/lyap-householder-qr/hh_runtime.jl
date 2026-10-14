# AC9 to AC12: invariance, zero allocation, health and scaling, state independence.

function hh_ac9(alg, sys, x0, W0, seed, n, chunks)
    reducers = (QRSpectrum(alg),)
    serial = [propagate_orbit(sys, x0[j], W0[j], n, 1, reducers; orbit = j) for j in eachindex(x0)]
    sx, sW = [s[1] for s in serial], [s[2] for s in serial]
    sf = [map(reducer_final, reducers, s[3]) for s in serial]
    for chunk in chunks
        res = propagate(sys, x0, n; k = 2, tau = 1, seed = seed, reducers = reducers, chunk = chunk, W0 = W0)
        @test isequal(res.x, sx)
        @test isequal(res.W, sW)          # final Q
        @test isequal(res.finals, sf)     # lambda, G_max, orth_max, count_max
    end
end

@testset "AC9 thread and chunk invariance" begin
    # K = 2 pi, 4096 seeded orbits, n = 1e3, k = 2, tau = 1; propagate (any thread count of this process) equals the serial map of
    # propagate_orbit bitwise for chunk sizes 64, 1000 and the default; kernel path here, library path in the nested testset
    seed = hh_seed(9)
    m, n = 4096, 1000
    sys = StandardMap(HH_TWO_PI)
    x0 = hh_starts(StableRNG(seed), m)
    W0 = [initial_deviations(seed, j, Val(2), Val(2), Float64) for j in 1:m]
    for alg in hh_kernel_algs()
        hh_ac9(alg, sys, x0, W0, seed, n, alg isa HouseholderQR ? (64, 1000, nothing) : (64, 1000))
    end
    @testset "AC9 library path (alignment-dependent library kernels may differ: record and ask the owner)" begin
        hh_ac9(LapackQR(), sys, x0, W0, seed, n, (64, 1000))
    end
end

function hh_step_bytes(r, st, x, W)
    reducer_step(r, st, 10, x, W)
    return @allocated reducer_step(r, st, 10, x, W)
end
function hh_orbit_bytes(sys, x0, W0, N, tau, reducers)
    propagate_orbit(sys, x0, W0, N, tau, reducers)
    return @allocated propagate_orbit(sys, x0, W0, N, tau, reducers)
end
function hh_ac10(alg, ::Type{T}, vd::Val{d}, vk::Val{k}, seed) where {T, d, k}
    W = hh_randW(StableRNG(seed), T, vd, vk)
    st = reducer_init(QRSpectrum(alg), T, vd, vk)
    return hh_step_bytes(QRSpectrum(alg), st, zero(SVector{d, T}), W)
end
const HH_AC10_DK = ((Val(2), Val(1)), (Val(2), Val(2)), (Val(4), Val(4)), (Val(8), Val(8)))

@testset "AC10 zero allocation" begin
    # one handoff update (reducer_step) allocates 0 bytes after warm-up: kernel path, T in {Float64, Double64}, (d,k) in {(2,1), (2,2), (4,4), (8,8)};
    # plus propagate_orbit with the standard map (kernel path, d = k = 2); BigFloat exempt
    seed = hh_seed(10)
    for T in (Float64, Double64), alg in hh_kernel_algs(), (i, (vd, vk)) in enumerate(HH_AC10_DK)
        @test hh_ac10(alg, T, vd, vk, seed + i) == 0
    end
    for T in (Float64, Double64), alg in hh_kernel_algs()
        sys = StandardMap(T(HH_TWO_PI))
        @test hh_orbit_bytes(sys, SVector{2, T}(T(1.0), T(2.0)), SMatrix{2, 2, T}(1, 0, 0, 1), 100, 1, (QRSpectrum(alg),)) == 0
    end
    @testset "AC10 library path (reducer_step with the reused workspace; report to the owner if it cannot reach 0)" begin
        for (i, (vd, vk)) in enumerate(HH_AC10_DK)
            @test hh_ac10(LapackQR(), Float64, vd, vk, seed + 10 + i) == 0
        end
    end
end

# one degenerate handoff through the stream: orbit 3 of 5 (flagged) has W = M at the handoff n = 15 (tau = 5, fault at step 12)
function hh_ac11_stream(alg, M, col)
    seed = hh_seed(11)
    I2 = SMatrix{2, 2, Float64}(1, 0, 0, 1)
    sys = hh_Degen(12, M)
    x0 = [SVector(0.0, i == 3 ? 1.0 : 0.0) for i in 1:5]
    r = (QRSpectrum(alg),)
    for err in (hh_health_error(() -> propagate(sys, x0, 20; k = 2, tau = 5, seed = seed, W0 = fill(I2, 5), reducers = r, chunk = 2)),
                hh_health_error(() -> propagate_orbit(sys, x0[3], I2, 20, 5, r; orbit = 3)))
        @test err isa NumericalHealthError
        err === nothing && continue
        @test err.what == "orthonormalisation breakdown"
        @test err.orbit == 3 && err.step == 15 && err.column == col && err.backend === :cpu
        @test err.value == 0.0                    # the pivot (contract table); the threshold of the zero case is not asserted (Q11)
    end
    ok = propagate(sys, x0[[1, 2, 4, 5]], 20; k = 2, tau = 5, seed = seed, W0 = fill(I2, 4), reducers = r)
    @test all(f -> all(iszero, f[1].sum) && f[1].n == 20, ok.finals)
end

# direct reducer_step on a good state, then a degenerate W: recorded failure, accumulators and W unchanged
function hh_ac11_direct(alg, Wbad::SMatrix{d, k, Float64}, Wgood, col) where {d, k}
    r = QRSpectrum(alg)
    x = zero(SVector{d, Float64})
    st0 = reducer_init(r, Float64, Val(d), Val(k))
    st1, _ = reducer_step(r, st0, 1, x, Wgood)
    @test !reducer_failed(r, st1)
    st2, W2 = reducer_step(r, st1, 2, x, Wbad)
    @test reducer_failed(r, st2)
    @test W2 === Wbad
    f1, f2 = reducer_final(r, st1), reducer_final(r, st2)
    @test f2.sum == f1.sum && f2.comp == f1.comp
    e = reducer_error(r, st2, 7, :cpu)
    @test e isa NumericalHealthError && e.what == "orthonormalisation breakdown"
    @test e.orbit == 7 && e.column == col && e.step == 2 && e.backend === :cpu && e.value == 0.0
    @test_throws NumericalHealthError reducer_update(r, st1, 2, x, Wbad)
end

function hh_ac11_scale(::Type{T}, vd::Val{d}, vk::Val{k}, s, alg, seed) where {T, d, k}
    eps_T = Float64(effective_eps(T))
    Qt = initial_deviations(seed, 1, vd, vk, T)
    W = T(s) * Qt
    Q, R, cnt, status, col = orthonormalize(alg, W)
    @test status == QR_OK
    @test all(isfinite, Q) && all(isfinite, R)
    @test all(i -> abs(R[i, i] / T(s) - 1) <= 8 * d * eps_T, 1:k)
    @test Float64(maximum(abs, Q - Qt)) <= 16 * d * eps_T
    r = QRSpectrum(alg)
    st, _ = reducer_step(r, reducer_init(r, T, vd, vk), 1, zero(SVector{d, T}), W)
    @test !reducer_failed(r, st)
    lam = reducer_final(r, st).lambda
    @test all(i -> abs(Float64(lam[i]) - log(s)) <= 8 * d * eps_T + 4 * eps_T * abs(log(s)), 1:k)   # no overflow or underflow of the norm
end

@testset "AC11 health and scaling" begin
    seed = hh_seed(11)
    # (a) zero column and (b) identical columns (dyadic data: exact zero pivot in every variant), kind "orthonormalisation breakdown" (D34, Q3)
    cases = (("zero column 2", SMatrix{2, 2, Float64}(1, 2, 0, 0), 2),
             ("identical columns", SMatrix{2, 2, Float64}(1, 0, 1, 0), 2),
             ("zero column 1", SMatrix{2, 2, Float64}(0, 0, 1, 2), 1))
    Wgood = SMatrix{2, 2, Float64}(2, 1, 0.5, 3)
    for alg in hh_algs(Float64), (name, M, col) in cases
        @testset "AC11 $(hh_name(alg)) $name through the stream" begin
            hh_ac11_stream(alg, M, col)
        end
        @testset "AC11 $(hh_name(alg)) $name direct" begin
            hh_ac11_direct(alg, M, Wgood, col)
        end
    end
    for alg in hh_algs(Float64)
        @testset "AC11 $(hh_name(alg)) identical columns, d = 4" begin
            hh_ac11_direct(alg, SMatrix{4, 2, Float64}(1, 0, 0, 0, 1, 0, 0, 0), SMatrix{4, 2, Float64}(1, 2, 0, 1, 0, 1, 3, 0.5), 2)
        end
    end
    # (c) W = s Q_test, s in {1e-200, 1e200}
    for T in (Float64, Double64), (i, (vd, vk)) in enumerate(((Val(2), Val(2)), (Val(4), Val(4)), (Val(4), Val(2)))), s in (1e-200, 1e200), alg in hh_algs(T)
        hh_ac11_scale(T, vd, vk, s, alg, seed + 20 + i)
    end
end

@testset "AC12 state independence" begin
    # tangent-core W9 map (state length 4, dimension 2) against the plain length-2 system: lambda, G, Q bitwise identical, k in {1, 2}
    seed = hh_seed(12)
    N, tau = 500, 2
    for T in (Float64, Double64), k in (1, 2), alg in hh_algs(T)
        K = T(HH_TWO_PI)
        z2 = [hh_cast(T, z) for z in hh_starts(StableRNG(seed), 20)]
        z4 = [SVector(z[1], z[2], zero(T), zero(T)) for z in z2]
        long = propagate(hh_Long(K), z4, N; k = k, tau = tau, seed = seed, reducers = (QRSpectrum(alg),))
        short = propagate(hh_Short(K), z2, N; k = k, tau = tau, seed = seed, reducers = (QRSpectrum(alg),))
        @test isequal(long.W, short.W)
        @test isequal(long.finals, short.finals)
    end
end
