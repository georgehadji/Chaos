# AC6 to AC10: tangencies and angles, flow geometry on F3, checkpointing and storage, k < d and column independence, execution invariance and allocation.
# Every testset starts by naming the new symbols it uses, so a missing module errors that testset (UndefVarError).

# ------------------------------------------------------------------ AC6: tangencies and angles
@testset "AC6a closed-form angles: the cat map (pi/2) and the non-normal symplectic map (arctan(1/3), pi/4), Float64 and Double64" begin
    _ = (ClvGinelli, propagate_clv_orbit)
    for T in (Float64, Double64)
        setprecision(BigFloat, CLV_BITS) do
            fc = clv_fixture(:c, T); W0 = clv_frame(clv_seed(6, 1), 1, 2, 2, T); C0 = clv_rand_C(StableRNG(clv_seed(6, 2)), 2, T)
            rc = clv_go(fc.map, zero(SVector{2, T}), 100, 20, 100; W0 = W0, C0 = C0, angle_sets = (1:1, 2:2))
            for n in 0:20
                @test abs(rc.cosines[1, 2, n + 1]) <= 1.0e3 * clv_eps(T)
                @test Float64(abs(clv_big(rc.angles[1, n + 1]) - BigFloat(pi) / 2)) <= 1.0e3 * clv_eps(T)
            end
            fb = clv_fixture(:b, T); W4 = clv_frame(clv_seed(6, 3), 1, 4, 4, T); C4 = clv_rand_C(StableRNG(clv_seed(6, 4)), 4, T)
            rb = clv_go(fb.map, zero(SVector{4, T}), 100, 20, 100; W0 = W4, C0 = C4, angle_sets = (1:2, 3:4))
            B = clv_B_fix(T, fb.gamma, 100); exact = (atan(BigFloat(1) / 3), BigFloat(pi) / 4)
            for n in 0:20, q in 1:2
                @test Float64(abs(clv_big(rb.angles[q, n + 1]) - exact[q])) <= B                              # ascending principal angles
            end
        end
    end
end

@testset "AC6b pairwise cosines from C agree with those of V = G C within 16 eps (maps and F3)" begin
    _ = (ClvGinelli, propagate_clv_orbit)
    T = Float64; rng = StableRNG(clv_seed(6, 5))
    systems = Any[(clv_cm(T), clv_starts(rng, 4, 1, T)[1], 1, 1.0), (clv_std(T), clv_starts(rng, 2, 1, T)[1], 1, 1.0), (clv_im(T), clv_orbit_C(), 20, CLV_H)]
    worst = 0.0
    for (is, (map, x0, tau, h)) in enumerate(systems)
        d = dimension(map); W0 = clv_frame(clv_seed(6, 10 + is), 1, d, d, T); C0 = clv_rand_C(StableRNG(clv_seed(6, 20 + is)), d, T)
        res = clv_go(map, x0, 40, clv_size(40, 200), 40; tau = tau, h = h, W0 = W0, C0 = C0)
        for n in 0:(size(res.C, 3) - 1), i in 1:d, j in 1:d
            V = res.V[:, :, n + 1]; C = res.C[:, :, n + 1]
            e1 = abs(res.cosines[i, j, n + 1] - abs(dot(V[:, i], V[:, j]))); e2 = abs(res.cosines[i, j, n + 1] - abs(dot(C[:, i], C[:, j])))
            worst = max(worst, max(e1, e2) / (16 * clv_eps(T))); @test max(e1, e2) <= 16 * clv_eps(T)
        end
    end
    @info "AC6b worst ratio error / (16 eps)" worst
end

@testset "AC6c planted tangency: the angle equals arctan(delta) to relative 1e-2 for delta in {1e-2, 1e-4, 1e-6}, Float64 and Double64" begin
    _ = (ClvGinelli, propagate_clv_orbit)
    for T in (Float64, Double64), delta in (1.0e-2, 1.0e-4, 1.0e-6)
        map = clv_planted(delta, T); W0 = clv_frame(clv_seed(6, 30), 1, 2, 2, T); C0 = clv_rand_C(StableRNG(clv_seed(6, 31)), 2, T)
        res = clv_go(map, zero(SVector{2, T}), 100, 5, 100; W0 = W0, C0 = C0, angle_sets = (1:1, 2:2))
        for n in 0:5
            ang = Float64(res.angles[1, n + 1]); @test abs(ang / atan(delta) - 1) <= 1.0e-2
            @test ang <= clv_pair_angle(res.V[:, 1, n + 1], res.V[:, 2, n + 1]) * (1 + 1.0e-2)             # the smallest principal angle never exceeds the pairwise angle of its two vectors
        end
    end
end

# ------------------------------------------------------------------ AC7: flow geometry on F3 (chaotic starts; orbit C has gap 0.007 and needs n* above 3000, see the plan, U5)
@testset "AC7 flow geometry on F3: f(z) lies in span(v_2, v_3), grad H is orthogonal to v_1 and v_4, the forward exponents pair (V15)" begin
    _ = (ClvGinelli, propagate_clv_orbit)
    T = Float64; flow = HenonHeilesFlow(); im = clv_im(T); starts = clv_chaotic_starts(clv_size(2, 20), 1); worst = 0.0
    for (is, z0) in enumerate(starts)
        lam = clv_pilot(im, z0, 20, CLV_H); Nn = clv_nstar(lam, 1.0; skip23 = true); ns = clv_size(30, 500)
        W0 = clv_frame(clv_seed(7, is), 1, 4, 4, T); C0 = clv_rand_C(StableRNG(clv_seed(7, 20 + is)), 4, T)
        res = clv_go(im, z0, Nn, ns, Nn; tau = 20, h = CLV_H, W0 = W0, C0 = C0, outputs = (:V,))
        pts = clv_orbit(im, res.x0, ns * 20)
        for n in 0:ns
            z = pts[20 * n + 1]; f = vector_field(flow, z); g = hamiltonian_gradient(flow, z); V = res.V[:, :, n + 1]
            a = clv_sin_sub(reshape(f, 4, 1), V[:, 2:3]); b = max(abs(dot(g, V[:, 1])), abs(dot(g, V[:, 4]))) / norm(g)
            worst = max(worst, a, b); @test a <= 1.0e-6 && b <= 1.0e-6
        end
    end
    # (c) the forward exponents of the window, at least 95 percent of the chaotic starts (the delivered convention of CQ criterion 11; package: shorter windows)
    m = clv_size(10, 200); cs = clv_chaotic_starts(m, 2); nok = 0
    for (is, z0) in enumerate(cs)
        res = clv_go(im, z0, 0, clv_size(300, 2000), 1; tau = 20, h = CLV_H, W0 = clv_frame(clv_seed(7, 100 + is), 1, 4, 4, T), C0 = clv_eye(4, 4, T), outputs = ())
        l = res.lambda_fwd; abs(l[2]) <= 2.0e-2 && abs(l[3]) <= 2.0e-2 && abs(l[1] + l[4]) <= 2.0e-2 && (nok += 1)
    end
    @test nok >= 0.95 * length(cs)
    @info "AC7 worst sine/cosine and forward-exponent pass fraction" worst nok = nok / length(cs)
end

# ------------------------------------------------------------------ AC8: checkpointing and storage
@testset "AC8 bitwise equality of the checkpointed and the store-all runs, h_b in {1, 7, 64, 500}, a shorter last block" begin
    _ = (ClvGinelli, propagate_clv_orbit)
    T = Float64; map = clv_cm(T); x0 = clv_starts(StableRNG(clv_seed(8, 1)), 4, 1, T)[1]
    W0 = clv_frame(clv_seed(8, 2), 1, 4, 4, T); C0 = clv_rand_C(StableRNG(clv_seed(8, 3)), 4, T); nf, ns, nb = 30, 100, 100
    ref = clv_go(map, x0, nf, ns, nb; W0 = W0, C0 = C0, angle_sets = (1:2, 3:4))
    for hb in (1, 7, 64, 500)
        ck = clv_go(map, x0, nf, ns, nb; W0 = W0, C0 = C0, angle_sets = (1:2, 3:4), cfg = (storage = :checkpoint, block = hb))
        @test ck.C == ref.C && ck.V == ref.V && ck.lambda == ref.lambda && ck.cosines == ref.cosines && ck.R == ref.R && ck.angles == ref.angles && ck.x == ref.x
    end
    im = clv_im(T); z = clv_orbit_C(); Wf = clv_frame(clv_seed(8, 4), 1, 4, 4, T)
    a = clv_go(im, z, 10, 30, 30; tau = 20, h = CLV_H, W0 = Wf, C0 = C0); b = clv_go(im, z, 10, 30, 30; tau = 20, h = CLV_H, W0 = Wf, C0 = C0, cfg = (storage = :checkpoint, block = 7))
    @test a.C == b.C && a.V == b.V && a.lambda == b.lambda
end

@testset "AC8b number of interval-map evaluations: n_f + N (store-all) and n_f + (N - h_last) + N (checkpointed)" begin
    _ = (ClvGinelli, propagate_clv_orbit)
    T = Float64; x0 = clv_starts(StableRNG(clv_seed(8, 5)), 4, 1, T)[1]; W0 = clv_frame(clv_seed(8, 6), 1, 4, 4, T); C0 = clv_rand_C(StableRNG(clv_seed(8, 7)), 4, T)
    nf, ns, nb = 30, 100, 100; N = ns + nb
    cnt = clv_CountMap(clv_cm(T)); clv_go(cnt, x0, nf, ns, nb; W0 = W0, C0 = C0); @test cnt.nf[] == nf + N && cnt.nt[] == nf + N
    for hb in (1, 7, 64, 500)
        cnt = clv_CountMap(clv_cm(T)); clv_go(cnt, x0, nf, ns, nb; W0 = W0, C0 = C0, cfg = (storage = :checkpoint, block = hb))
        nc = cld(N, hb); hlast = N - (nc - 1) * hb; expected = nf + (N - hlast) + N
        @test cnt.nf[] == expected && cnt.nt[] == expected
    end
end

@testset "AC8c memory formula: bytes of the stored R factors, frames and working buffers equal clv_bytes of V17 (288576 and 24528 for the examples)" begin
    _ = (ClvGinelli, propagate_clv_orbit, clv_bytes)
    @test clv_bytes(Val(4), Val(4), Float64; n_s = 1000, n_b = 1000) == 288576
    @test clv_bytes(Val(4), Val(4), Float64; n_s = 1000, n_b = 1000, storage = :checkpoint, block = 100) == 24528
    T = Float64; map = clv_cm(T); x0 = clv_starts(StableRNG(clv_seed(8, 8)), 4, 1, T)[1]; W0 = clv_frame(clv_seed(8, 9), 1, 4, 4, T); C0 = clv_rand_C(StableRNG(clv_seed(8, 10)), 4, T)
    for (ns, nb, st, hb) in ((1000, 1000, :store, 0), (1000, 1000, :checkpoint, 100), (100, 100, :checkpoint, 7), (50, 80, :store, 0))
        cfg = st === :store ? (;) : (storage = :checkpoint, block = hb)
        res = clv_go(map, x0, 10, ns, nb; W0 = W0, C0 = C0, outputs = (:V,), cfg = cfg)
        @test res.bytes == clv_bytes(Val(4), Val(4), T; n_s = ns, n_b = nb, storage = st === :store ? :store : :checkpoint, block = hb)
    end
    ra = clv_go(map, x0, 10, 100, 100; W0 = W0, C0 = C0, outputs = ())                                    # angles and exponents only: no frame term (V17)
    @test ra.bytes == clv_bytes(Val(4), Val(4), T; n_s = 100, n_b = 100, frames = false)
    @test clv_bytes(Val(4), Val(4), Double64; n_s = 100, n_b = 100) == 2 * clv_bytes(Val(4), Val(4), Float64; n_s = 100, n_b = 100) - 16   # sizeof scales, the 16 bytes of counters do not
end

@testset "AC8d streaming sink: samples delivered one at a time in decreasing n, nothing stored, peak bytes within bytes_ckpt plus one block of output" begin
    _ = (ClvGinelli, propagate_clv_orbit, clv_bytes)
    T = Float64; map = clv_cm(T); x0 = clv_starts(StableRNG(clv_seed(8, 11)), 4, 1, T)[1]; W0 = clv_frame(clv_seed(8, 12), 1, 4, 4, T); C0 = clv_rand_C(StableRNG(clv_seed(8, 13)), 4, T)
    nf, ns, nb, hb = 10, 100, 100, 16
    ref = clv_go(map, x0, nf, ns, nb; W0 = W0, C0 = C0)
    got = Dict{Int, Any}(); order = Int[]
    sink = (n, C, lam) -> (push!(order, n); got[n] = (copy(C), copy(lam)); nothing)
    res = clv_go(map, x0, nf, ns, nb; W0 = W0, C0 = C0, cfg = (storage = :checkpoint, block = hb), sink = sink, outputs = ())
    @test res.C === nothing && order == collect(ns:-1:0)
    @test all(n -> got[n][1] == ref.C[:, :, n + 1] && got[n][2] == ref.lambda[:, n + 1], 0:ns)
    @test res.bytes <= clv_bytes(Val(4), Val(4), T; n_s = ns, n_b = nb, storage = :checkpoint, block = hb) + sizeof(T) * hb * (2 * 16 + 4)
end

# ------------------------------------------------------------------ AC9: k < d and column independence
@testset "AC9 first k columns of C, V and Lambda are those of the k = d run, bitwise, for Householder and repeated Gram-Schmidt; k = 0 or k > d is rejected" begin
    _ = (ClvGinelli, propagate_clv_orbit)
    T = Float64; rng = StableRNG(clv_seed(9, 1))
    systems = Any[(clv_cm(T), clv_starts(rng, 4, 1, T)[1], 1, 1.0), (clv_im(T), clv_orbit_C(), 20, CLV_H)]
    for (is, (map, x0, tau, h)) in enumerate(systems), alg in (HouseholderQR(), RGS()), k in (1, 2, 4)
        W0 = clv_frame(clv_seed(9, 2), 1, 4, 4, T); C0 = clv_rand_C(StableRNG(clv_seed(9, 3)), 4, T)
        full = clv_go(map, x0, 20, 40, 20; tau = tau, h = h, W0 = W0, C0 = C0, alg = alg)
        part = clv_go(map, x0, 20, 40, 20; tau = tau, h = h, W0 = SMatrix{4, k, T}(W0[:, 1:k]), C0 = SMatrix{k, k, T}(C0[1:k, 1:k]), alg = alg)
        @test part.C == full.C[1:k, 1:k, :] && part.V == full.V[:, 1:k, :] && part.lambda == full.lambda[1:k, :] && part.lambda_fwd == full.lambda_fwd[1:k]
    end
    map = clv_cm(T); x0 = clv_starts(StableRNG(clv_seed(9, 4)), 4, 1, T)[1]
    @test_throws ArgumentError clv_go(map, x0, 5, 5, 5; k = 0)
    @test_throws ArgumentError clv_go(map, x0, 5, 5, 5; k = 5)
end

@testset "AC9c starting coefficient matrix: upper triangular, positive diagonal, unit columns, column j independent of k, keyed by (seed, orbit) (P13)" begin
    _ = (clv_initial_coefficients,)
    for T in (Float64, Double64)
        C4 = clv_initial_coefficients(clv_seed(9, 10), 3, Val(4), T); C2 = clv_initial_coefficients(clv_seed(9, 10), 3, Val(2), T); C1 = clv_initial_coefficients(clv_seed(9, 10), 3, Val(1), T)
        @test all(i -> all(j -> i <= j || iszero(C4[i, j]), 1:4), 1:4) && all(i -> C4[i, i] > 0, 1:4)
        @test all(j -> abs(norm(C4[:, j]) - 1) <= 8 * clv_eps(T), 1:4)
        @test C4[1:2, 1:2] == C2 && C4[1, 1] == C1[1, 1] && C4[1:3, 3] == clv_initial_coefficients(clv_seed(9, 10), 3, Val(3), T)[1:3, 3]   # column j does not depend on k
        @test clv_initial_coefficients(clv_seed(9, 10), 3, Val(4), T) == C4                                                  # same seed and orbit: same matrix
        @test clv_initial_coefficients(clv_seed(9, 10), 4, Val(4), T) != C4 && clv_initial_coefficients(clv_seed(9, 11), 3, Val(4), T) != C4
    end
end

# ------------------------------------------------------------------ AC10: execution invariance and allocation
@testset "AC10a execution invariance: 1 against 4 tasks, chunk sizes 64 and 1000, ensemble position; orbit points bitwise those of the integrator and of the delivered spectrum" begin
    _ = (ClvGinelli, propagate_clv_orbit, propagate_clv)
    T = Float64; M = clv_size(96, 1024); nf, ns, nb = clv_size(40, 150), clv_size(40, 200), clv_size(40, 200)
    xs = clv_starts(StableRNG(clv_seed(10, 1)), 2, M, T); map = clv_std(T); clv = ClvGinelli(map)
    r1 = propagate_clv(clv, xs, nf, ns, nb; k = 2, seed = 1, chunk = 64, ntasks = 1, outputs = (:V,))
    r4 = propagate_clv(clv, xs, nf, ns, nb; k = 2, seed = 1, chunk = 1000, ntasks = 4, outputs = (:V,))
    @test length(r1.results) == M && length(r4.results) == M
    @test all(j -> r1.results[j].C == r4.results[j].C && r1.results[j].V == r4.results[j].V && r1.results[j].lambda == r4.results[j].lambda && r1.results[j].x == r4.results[j].x, 1:M)
    for j in unique((1, 17, M))
        single = propagate_clv_orbit(clv, xs[j], nf, ns, nb; k = 2, seed = 1, orbit = j, outputs = (:V,))
        @test single.C == r1.results[j].C && single.V == r1.results[j].V && single.lambda == r1.results[j].lambda
    end
    for k in (1, 2), j in (1, 5)
        res = propagate_clv_orbit(clv, xs[j], nf, ns, nb; k = k, seed = 1, orbit = j)
        @test res.x == last(clv_orbit(map, xs[j], nf + ns + nb))
        sp = propagate(map, [xs[j]], nf + ns + nb; k = k, tau = 1, seed = 1, reducers = (QRSpectrum(HouseholderQR()),))
        @test res.x == sp.x[1]
    end
end

@testset "AC10b forward exponents are bitwise those of the delivered discrete spectrum (O-c): standard map and F3" begin
    _ = (ClvGinelli, propagate_clv_orbit)
    T = Float64
    for (map, x0, tau, h, d) in ((clv_std(T), SVector{2, T}(1.3, 0.7), 1, 1.0, 2), (clv_im(T), clv_orbit_C(), 20, CLV_H, 4))
        W0 = clv_frame(clv_seed(10, 5), 1, d, d, T); ns = clv_size(60, 200)
        for alg in (HouseholderQR(), RGS())
            res = clv_go(map, x0, 0, ns, 10; tau = tau, h = h, W0 = W0, C0 = clv_eye(d, d, T), alg = alg)
            sp = propagate(map, [x0], ns * tau; k = d, tau = tau, seed = 1, reducers = (QRSpectrum(alg),), W0 = [W0]).finals[1][1].lambda ./ h
            @test collect(res.lambda_fwd) == collect(sp)
        end
    end
end

@testset "AC10c no allocation per handoff in the forward and backward loops, Float64 and Double64 (200 extra handoffs add nothing beyond one checkpoint)" begin
    _ = (ClvGinelli, propagate_clv_orbit)
    for T in (Float64, Double64)
        map = clv_std(T); x0 = SVector{2, T}(T(1.3), T(0.7)); W0 = clv_frame(clv_seed(10, 8), 1, 2, 2, T); C0 = clv_rand_C(StableRNG(clv_seed(10, 9)), 2, T)
        go(clvobj, nf, nb) = propagate_clv_orbit(clvobj, x0, nf, 100, nb; k = 2, seed = 1, W0 = W0, C0 = C0, outputs = ())
        c0 = ClvGinelli(map); a = clv_alloc_bytes(go, c0, 100, 100); b = clv_alloc_bytes(go, c0, 300, 100)
        @test abs(b - a) <= 256                                                                            # the forward transient stores nothing; 256 bytes absorb @allocated noise of the result arrays (D123), 200 handoffs at 16 bytes each would add 3200
        cc = ClvGinelli(map; storage = :checkpoint, block = 200); a2 = clv_alloc_bytes(go, cc, 100, 100); b2 = clv_alloc_bytes(go, cc, 100, 300)
        @test -256 <= b2 - a2 <= sizeof(T) * (2 + 2 * 2) + 256                                             # one more checkpoint (state and frame) at most, +-256 bytes of @allocated noise (D123); per-handoff bytes would be far above
    end
end
