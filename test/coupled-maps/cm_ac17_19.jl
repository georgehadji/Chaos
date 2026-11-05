# Criteria 17 to 19 (package-test parts): banded tangent action, map-level cross-backend checks, zero allocation.
# The timing part of 17, the GPU parts of 18 and the 1 vs 8 thread comparison are stage 04 (see the test plan).

# ---- AC17: action of M and M^-1 on a frame of 4 vectors, banded structure vs the dense product ----
function cm_ac17_ratios(c::CMCase, rng)
    sys = cm_sys(Float64, c)
    d = cm_d(c)
    z = cm_state(Float64, c, cm_draw(rng, d, 4 * Float64(pi)))
    W = reshape(cm_draw(rng, d * 4, 1.0), d, 4)
    tol = 16 * cm_eps(Float64) * (1 + cm_amax(c))
    out = zeros(d, 4)
    tangent_action!(out, sys, z, W)
    r1 = maximum(abs, out - cm_ref_tangent(c, c.p, z) * W) / tol
    r2 = maximum(abs, out - Matrix(tangent(sys, z)) * W) / tol
    zi = forward(sys, z)
    inverse_tangent_action!(out, sys, zi, W)
    r3 = maximum(abs, out - cm_ref_inverse_tangent(c, c.p, zi) * W) / tol
    r4 = maximum(abs, out - Matrix(inverse_tangent(sys, zi)) * W) / tol
    return max(r1, r2, r3, r4)
end

@testset "AC17 large N and bandwidth" begin
    # N = 100 (d = 200), open chain and ring, (K, eps) in {(3, 0.6102), (5, 1)}: M W and M^-1 W for 4 seeded vectors, banded vs dense, 16 eps (1 + a_max) per entry;
    # also N = 10 and 50; the dense references are the module tangent and the plain-formula tangent (the time scaling is stage 04)
    rng = StableRNG(CM_SEED + 17001)
    worst = 0.0
    for N in (10, 50, 100), ring in (false, true), (K, e) in ((3.0, 0.6102), (5.0, 1.0))
        worst = max(worst, cm_ac17_ratios(cm_y4(N, K, e; ring), rng))
    end
    @test worst <= 1
end

# ---- AC18 (map level): serial vs threaded evaluation, CPU kernel vs plain loop ----
function cm_run_serial(sys, Z, n)
    out = similar(Z)
    for j in eachindex(Z)
        z = Z[j]
        for _ in 1:n
            z = forward(sys, z)
        end
        out[j] = z
    end
    return out
end

function cm_run_threaded(sys, Z, n)
    out = similar(Z)
    Threads.@threads for j in eachindex(Z)
        z = Z[j]
        for _ in 1:n
            z = forward(sys, z)
        end
        out[j] = z
    end
    return out
end

@testset "AC18 map-level reproducibility" begin
    # Y3 at the published parameters and (2, 0.5, 0.1): 10^5 seeded orbits, n = n_H, serial loop vs Threads.@threads loop bitwise;
    # Y4 N in {2, 5}: 10^4 orbits (Vector states); wrapped Y3 likewise (10^4)
    rng = StableRNG(CM_SEED + 18001)
    lim = 4 * Float64(pi)
    for p in ((0.5, 0.1, 1.0e-3), (2.0, 0.5, 0.1))
        c = cm_y3(p...)
        Z = [SVector{4, Float64}(ntuple(_ -> (2 * rand(rng) - 1) * lim, 4)) for _ in 1:100_000]
        sys = cm_sys(Float64, c)
        @test isequal(cm_run_serial(sys, Z, cm_nH(c)), cm_run_threaded(sys, Z, cm_nH(c)))
        Zw = [SVector{8, Float64}(ntuple(i -> i <= 4 ? 2 * Float64(pi) * rand(rng) : 0.0, 8)) for _ in 1:10_000]
        wsys = cm_wsys(Float64, c)
        @test isequal(cm_run_serial(wsys, Zw, cm_nH(c)), cm_run_threaded(wsys, Zw, cm_nH(c)))
    end
    for N in (2, 5)
        c = cm_y4(N, 3.0, 0.2)
        Z = [cm_draw(rng, 2N, lim) for _ in 1:10_000]
        sys = cm_sys(Float64, c)
        @test isequal(cm_run_serial(sys, Z, cm_nH(c)), cm_run_threaded(sys, Z, cm_nH(c)))
    end
end

@testset "AC18 CPU kernel vs plain loop" begin
    # propagate (CPU backend, KernelAbstractions) realizes the orbit of forward bitwise: Y3 lifted and wrapped, 2000 orbits, n = n_H, k = 4
    rng = StableRNG(CM_SEED + 18002)
    for p in ((0.5, 0.1, 1.0e-3), (2.0, 0.5, 0.1))
        c = cm_y3(p...)
        n = cm_nH(c)
        X0 = [SVector{4, Float64}(ntuple(_ -> (2 * rand(rng) - 1) * 4 * Float64(pi), 4)) for _ in 1:2000]
        res = propagate(cm_sys(Float64, c), X0, n; k = 4, seed = CM_SEED)
        ref = cm_run_serial(cm_sys(Float64, c), X0, n)
        @test all(j -> isequal(res.x[j], ref[j]), eachindex(X0))
        Xw = [SVector{8, Float64}(ntuple(i -> i <= 4 ? 2 * Float64(pi) * rand(rng) : 0.0, 8)) for _ in 1:2000]
        resw = propagate(cm_wsys(Float64, c), Xw, n; k = 4, seed = CM_SEED)
        refw = cm_run_serial(cm_wsys(Float64, c), Xw, n)
        @test all(j -> isequal(resw.x[j], refw[j]), eachindex(Xw))
    end
end

# ---- AC19: zero allocation ----
# SVector systems (Y3, lifted and wrapped): the pure functions; chain: the in-place functions with a preallocated workspace
function cm_alloc_pure(sys, z)
    forward(sys, z); inverse(sys, z); tangent(sys, z); inverse_tangent(sys, z)
    return (@allocated(forward(sys, z)), @allocated(inverse(sys, z)), @allocated(tangent(sys, z)), @allocated(inverse_tangent(sys, z)))
end

function cm_alloc_inplace(sys, z, out, M)
    forward!(out, sys, z); inverse!(out, sys, z); tangent!(M, sys, z); inverse_tangent!(M, sys, z)
    return (@allocated(forward!(out, sys, z)), @allocated(inverse!(out, sys, z)), @allocated(tangent!(M, sys, z)), @allocated(inverse_tangent!(M, sys, z)))
end

function cm_alloc_action(sys, z, WO, W)
    tangent_action!(WO, sys, z, W); inverse_tangent_action!(WO, sys, z, W)
    return (@allocated(tangent_action!(WO, sys, z, W)), @allocated(inverse_tangent_action!(WO, sys, z, W)))
end

@testset "AC19 zero allocation" begin
    # Y3 (lifted and wrapped): 0 bytes after warm-up, Float64 and Double64; Y4 N in {2, 10, 100}: in-place forms 0 bytes, lifted and wrapped, plus the frame actions;
    # the in-place forms agree with the pure forms bitwise
    for T in (Float64, Double64)
        c = cm_y3(0.5, 0.1, 1.0e-3)
        z = cm_state(T, c, [1.25, -0.5, 0.75, 2.0])
        @test all(==(0), cm_alloc_pure(cm_sys(T, c), z))
        @test all(==(0), cm_alloc_pure(cm_wsys(T, c), cm_wstate(T, c, [1.25, 0.5, 0.75, 2.0])))
    end
    for T in (Float64, Double64), N in (2, 10, 100)
        c = cm_y4(N, 3.0, 0.2)
        d = 2N
        rng = StableRNG(CM_SEED + 19001 + N)
        z = cm_state(T, c, cm_draw(rng, d, 4.0))
        sys = cm_sys(T, c)
        out, M = similar(z), zeros(T, d, d)
        @test all(==(0), cm_alloc_inplace(sys, z, out, M))
        forward!(out, sys, z)
        @test out == forward(sys, z)
        tangent!(M, sys, z)
        @test M == Matrix(tangent(sys, z))
        W, WO = reshape([T(v) for v in cm_draw(rng, d * 4, 1.0)], d, 4), zeros(T, d, 4)
        @test all(==(0), cm_alloc_action(sys, z, WO, W))
        wsys = cm_wsys(T, c)
        s = cm_wstate(T, c, cm_draw01(rng, d))
        sout, Mw = similar(s), zeros(T, d, d)
        @test all(==(0), cm_alloc_inplace(wsys, s, sout, Mw))
        forward!(sout, wsys, s)
        @test sout == forward(wsys, s)
    end
end
