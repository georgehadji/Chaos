# AC14 X starts (E12), AC15 R placement (E13), AC16 draw (E14), AC17 streams and seeds (E15, E16), AC18 Y5 control (E19b).
# Placement and draw functions take the stream as an argument; the references below are written from the spec.

@testset "AC14 X starts" begin
    rng = StableRNG(ex_seed(14))
    l_lo, l_hi = -4.0, -2.0
    # consumption order (E12): u1 then u2; r = 10^l, l = l_lo + (l_hi - l_lo) u1, phi = 2 pi u2
    probe = StableRNG(ex_seed(14, 1))
    s = EX.place_X(probe, l_lo, l_hi)
    ref = StableRNG(ex_seed(14, 1))
    u1, u2 = rand(ref), rand(ref)
    @test abs(s.r - 10^(l_lo + (l_hi - l_lo) * u1)) <= 1e-12 * s.r
    @test abs(s.phi - 2pi * u2) <= 1e-12
    @test abs(s.x - s.r * cos(s.phi)) <= 1e-12 * s.r && abs(s.p - s.r * sin(s.phi)) <= 1e-12 * s.r
    draws = [EX.place_X(rng, l_lo, l_hi) for _ in 1:10_000]
    @test all(d -> l_lo - 1e-12 <= log10(d.r) <= l_hi + 1e-12, draws)
    @test all(d -> abs(EX.torus_distance((d.x, d.p), (0.0, 0.0)) - d.r) <= 1e-12 * d.r, draws)      # r << pi: the torus distance is r
    @test ex_ks_uniform([(log10(d.r) - l_lo) / (l_hi - l_lo) for d in draws]) <= 1.95 / sqrt(10_000)
    @test ex_ks_uniform([d.phi / 2pi for d in draws]) <= 1.95 / sqrt(10_000)
end

@testset "AC15 R placement" begin
    # K > 0 (E13): index 1 + floor(u1 L), angle 2 pi u2, the sum formed at 256 bits and rounded once
    setprecision(BigFloat, 256) do
        pts = [SVector{2, BigFloat}(BigFloat(pi) * k / 3, BigFloat(pi) * (7 - k) / 5) for k in 1:5]
        dR = 1e-3
        for sub in 1:20
            r = EX.place_R(StableRNG(ex_seed(15, sub)), pts, dR)
            ref = StableRNG(ex_seed(15, sub))
            u1, u2 = rand(ref), rand(ref)
            @test r.index == 1 + floor(Int, u1 * length(pts))
            @test abs(r.phi - 2pi * u2) <= 1e-12
            zb = (pts[r.index][1] + BigFloat(dR) * cos(BigFloat(r.phi)), pts[r.index][2] + BigFloat(dR) * sin(BigFloat(r.phi)))
            @test r.x == Float64(zb[1]) && r.p == Float64(zb[2])             # bit for bit from (point, phi)
            @test abs(EX.torus_distance((r.x, r.p), (Float64(pts[r.index][1]), Float64(pts[r.index][2]))) - dR) <= 1e-12 * dR + 4 * EX_EPS * 8
        end
    end
    # K = 0: reduced fractions a/b, 1 <= b <= q_max, in the order (b, a); p0 = 2 pi a / b, x0 = 2 pi u2, no displacement
    q_max = 6
    fr = [(a, b) for b in 1:q_max for a in 0:(b - 1) if gcd(a, b) == 1]
    for sub in 1:30
        r = EX.place_R0(StableRNG(ex_seed(15, 20 + sub)), q_max)
        ref = StableRNG(ex_seed(15, 20 + sub))
        u1, u2 = rand(ref), rand(ref)
        a, b = fr[1 + floor(Int, u1 * length(fr))]
        @test (r.a, r.b) == (a, b) && gcd(r.a, r.b) == 1 && r.b <= q_max
        @test r.p == setprecision(() -> Float64(2 * BigFloat(pi) * a / b), BigFloat, EX_BITS)
        @test abs(r.x - 2pi * u2) <= 1e-12
    end
end

@testset "AC16 draw" begin
    @testset "distinct qualifying ids; equals the partial Fisher-Yates reference" begin
        cands = collect(3:3:300)                               # ascending ids of qualifying cells
        for sub in 1:20
            got = EX.draw_candidates(StableRNG(ex_seed(16, sub)), cands, 10)
            @test got == ex_ref_draw(StableRNG(ex_seed(16, sub)), cands, 10)
            @test length(unique(got)) == 10 && all(in(cands), got)
        end
    end
    @testset "inclusion frequency over 2000 stream seeds" begin
        M, n_c, S = 1000, 10, 2000
        cnt = zeros(Int, M)
        for s in rand(StableRNG(ex_seed(16, 30)), UInt64, S)                # unrelated stream seeds
            foreach(c -> cnt[c] += 1, EX.draw_candidates(StableRNG(s), collect(1:M), n_c))
        end
        mean_c = S * n_c / M
        sd = sqrt(S * (n_c / M) * (1 - n_c / M))
        @test all(c -> abs(c - mean_c) <= 5 * sd, cnt)
        @test sum(cnt) == S * n_c
    end
    @testset "short and empty cells" begin
        @test sort(EX.draw_candidates(StableRNG(ex_seed(16, 31)), [4, 8, 15], 10)) == [4, 8, 15]   # M < n_c: M orbits
        @test isempty(EX.draw_candidates(StableRNG(ex_seed(16, 32)), Int[], 10))                    # M = 0: nothing, never filled
    end
end

@testset "AC17 streams and seeds" begin
    seen = Set{UInt64}()
    for i in 1:10, k in ("T", "R", "C", "S", "X"), s in 0:9
        push!(seen, EX.cell_stream_seed(EX_SEED + s, i, k))
    end
    @test length(seen) == 500
    # the rule of E15: first 8 bytes, big endian, of SHA-256 over "chaoslab-panel-1|<seed>|<i>|<k>" (independent SHA-256)
    for (i, k) in ((1, "T"), (3, "S"), (6, "X"), (10, "C"))
        @test EX.cell_stream_seed(EX_SEED, i, k) == ex_cell_seed(EX_SEED, i, k)
        @test EX.cell_stream_seed(EX_SEED, i, Symbol(k)) == ex_cell_seed(EX_SEED, i, k)
    end
    # the first 20 uniforms of one fixed cell at the default seed: the project's stable generator keyed by that seed.
    # No stored list: a list of generator output is not a bibliography benchmark row, so it is re-derived here (plan, ambiguity A5).
    @test rand(EX.cell_rng(EX_SEED, 2, "C"), 20) == rand(StableRNG(ex_cell_seed(EX_SEED, 2, "C")), 20)
    # I9: a cell's stream does not depend on what another cell consumed
    a = rand(EX.cell_rng(EX_SEED, 4, "R"), 5)
    rand(EX.cell_rng(EX_SEED, 1, "S"), 100)
    @test rand(EX.cell_rng(EX_SEED, 4, "R"), 5) == a
    # E16: the per-orbit seed is tangent-core's rule: it keys the generator that makes initial_deviations(seed, id)
    for (s, id) in ((EX_SEED, 1), (EX_SEED, 77), (EX_SEED + 5, 1093))
        os = EX.orbit_seed(s, id)
        g = randn(StableRNG(os), 4)
        Q, R = qr(SMatrix{2, 2, Float64}(g...))
        Qs = Q * Diagonal(sign.(diag(R)))
        @test Qs[:, 1] ≈ initial_deviations(s, id, Val(2), Val(1), Float64)[:, 1] atol = 1e-15
    end
end

@testset "AC18 Y5 control" begin
    c = EX.y5_control((3 + sqrt(5)) / 2)
    @test collect(c.x0) == [0.0, 0.0]
    setprecision(BigFloat, EX_BITS) do
        mu = (3 + sqrt(BigFloat(5))) / 2
        ev = sort(eigvals(Matrix{Float64}(c.sys.A)))
        @test abs(ev[1] - Float64(1 / mu)) <= 16 * EX_EPS * Float64(1 / mu)
        @test abs(ev[2] - Float64(mu)) <= 16 * EX_EPS * Float64(mu)
        @test abs(c.exponent - Float64(log(mu))) <= 16 * EX_EPS * Float64(log(mu))
    end
    @test c.sys isa LinearSymplecticMap
    @test abs(det(Matrix{Float64}(c.sys.A)) - 1) <= 16 * EX_EPS
    # the exponent is realised by the tangent stream along the expanding eigenvector, from the fixed point (0, 0)
    dir = SMatrix{2, 1, Float64}(1 / sqrt(2), 1 / sqrt(2))                  # the eigenvector (1, 1) of the cosh/sinh block for mu
    r = propagate(c.sys, [SVector{2, Float64}(c.x0...)], 200; k = 1, tau = 1, seed = ex_seed(18), reducers = (Benettin(),), W0 = [dir])
    @test abs(r.finals[1][1].lambda[1] - c.exponent) <= 1e-12
end
