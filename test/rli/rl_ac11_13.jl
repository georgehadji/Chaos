# AC11 to AC13: health and frames, chain / pair / execution invariance, wrapped against lifted state.
# AC14 (ChaosTools two one-column tangent systems, verification environment) and AC15 (CPU vs GPU, D68 form) are stage 04; the test plan maps them.

# Pair states of the fault map: (n, y, flag) for block A then block B.
rl_fault_x0s(m, flagA, flagB) = [SVector{6, Float64}(0, 0, j in flagA ? 1 : 0, 0, 0, j in flagB ? 1 : 0) for j in 1:m]
rl_unit_frames(m, seed) = [rl_frame(rl_xi0(seed, j, Float64)) for j in 1:m]

@testset "AC11a zero tangent, column and pair id of the health error" begin
    _ = (RLI, PairMap)
    # A zero norm is not caught by the finiteness check of the stream: the reducer records it. Zero tangent in block A of a flagged pair: error at the first handoff, kind
    # "deviation matrix", the pair id, the step, backend :cpu, column 1 (the orbit A), value 0, nothing stored; in block B only: column 2; in both: column 1; two flagged pairs:
    # the lowest id first. Float64 through propagate (kernel path, chunk 2), the chains {RLI} and {RLI (non-replacing), MEGNO}; Double64 through the per-orbit body. RLI needs
    # k = 2 and an even dimension: ArgumentError before any work.
    seed = rl_seed(11, 1)
    P = PairMap(rl_Fault(:zero, 0))
    Ws = rl_unit_frames(6, seed)
    for chain in ((RLI(),), (RLI(replace = false), MEGNO()))
        ok = propagate(P, rl_fault_x0s(6, (), ()), 12; k = 2, tau = 3, seed = seed, reducers = chain, chunk = 2, W0 = Ws)
        @test length(ok.finals) == 6
        for (fa, fb, orbit, col) in (((4,), (), 4, 1), ((), (5,), 5, 2), ((3,), (3,), 3, 1))
            res = nothing
            err = rl_health_error(() -> (res = propagate(P, rl_fault_x0s(6, fa, fb), 12; k = 2, tau = 3, seed = seed, reducers = chain, chunk = 2, W0 = Ws)))
            @test err isa NumericalHealthError
            @test res === nothing
            @test err.what == "deviation matrix" && err.orbit == orbit && err.step == 3 && err.backend === :cpu && err.column == col && err.value == 0.0
        end
        res = nothing
        err = rl_health_error(() -> (res = propagate(P, rl_fault_x0s(6, (5,), (2,)), 12; k = 2, tau = 3, seed = seed, reducers = chain, chunk = 2, W0 = Ws)))
        @test err isa NumericalHealthError && err.orbit == 2 && err.column == 2 && err.step == 3 && res === nothing        # lowest pair id first
    end
    x = SVector{6, Double64}(0, 0, 1, 0, 0, 0)
    W0 = rl_frame(rl_xi0(seed, 1, Double64))
    err = rl_health_error(() -> propagate_orbit(P, x, W0, 12, 3, (RLI(),); orbit = 4))
    @test err isa NumericalHealthError && err.what == "deviation matrix" && err.orbit == 4 && err.step == 3 && err.backend === :cpu && err.column == 1 && err.value == 0.0
    err = rl_health_error(() -> propagate_orbit(P, SVector{6, Double64}(0, 0, 0, 0, 0, 1), W0, 12, 3, (RLI(),); orbit = 4))
    @test err isa NumericalHealthError && err.column == 2
    # RLI needs k = 2 (and an even dimension): ArgumentError before any work (rl_Boom throws if the map is called).
    e = try
        propagate(PairMap(rl_Boom()), [SVector(0.0, 0.0, 0.0, 0.0)], 4; k = 1, tau = 2, seed = seed, reducers = (RLI(),))
        nothing
    catch ex
        ex
    end
    @test e isa ArgumentError && occursin("RLI needs k = 2, got k = 1", e.msg)
    e = try
        reducer_init(RLI(), Float64, Val(3), Val(2))
        nothing
    catch ex
        ex
    end
    @test e isa ArgumentError && occursin("RLI needs an even dimension, got d = 3", e.msg)
end

@testset "AC11b non-finite state or deviation matrix is not masked" begin
    _ = (RLI, PairMap)
    # Non-finite entries in either block are tangent-core AC8; only that the RLI chains do not change that verdict: NaN in the tangent of block B: kind "deviation matrix",
    # column 0, threshold Inf; NaN in the state of block B: kind "state".
    seed = rl_seed(11, 2)
    Ws = rl_unit_frames(6, seed)
    Pn = PairMap(rl_Fault(:nan, 3))
    for chain in ((RLI(),), (RLI(replace = false), MEGNO()))
        res = nothing
        err = rl_health_error(() -> (res = propagate(Pn, rl_fault_x0s(6, (), (3,)), 20; k = 2, tau = 5, seed = seed, reducers = chain, chunk = 2, W0 = Ws)))
        @test err isa NumericalHealthError && res === nothing
        @test err.what == "deviation matrix" && err.orbit == 3 && err.step == 5 && err.backend === :cpu && err.column == 0 && err.threshold == Inf && !isfinite(err.value)
        x0s = rl_fault_x0s(6, (), ())
        x0s[2] = SVector{6, Float64}(0, 0, 0, NaN, 0, 0)
        res = nothing
        err = rl_health_error(() -> (res = propagate(Pn, x0s, 20; k = 2, tau = 5, seed = seed, reducers = chain, chunk = 2, W0 = Ws)))
        @test err isa NumericalHealthError && res === nothing && err.what == "state" && err.orbit == 2 && err.step == 5
    end
    x = SVector{6, Double64}(0, 0, 0, 0, 0, 1)
    err = rl_health_error(() -> propagate_orbit(Pn, x, rl_frame(rl_xi0(seed, 1, Double64)), 20, 5, (RLI(),); orbit = 3))
    @test err isa NumericalHealthError && err.what == "deviation matrix" && err.orbit == 3 && err.step == 5 && err.column == 0
end

@testset "AC11c supplied frame check" begin
    _ = (RLI, check_frame)
    # A supplied block frame is accepted only if both columns are unit vectors to roundoff, each column is supported in its own block, and the two block vectors are bitwise equal
    # (RL5, I4; sali Amendment D67 and fli D72 item 7 form). Rejected with a numerical health error of kind "initial frame" (step 0, the orbit id, the offending column, value above
    # the threshold, backend :cpu): a column of norm 2, 3, 0 or 1 + 1000 eps; a column with support in both blocks; two block vectors that differ (column 2). Accepted: the block
    # frames of seeded xi_0 for ids 1..1e4 (Double64: 1..1000), a frame with columns of norm 1 + eps (check_frame returns the larger defect, at most 2 eps), dependent or exchanged
    # unit vectors. ArgumentError: a frame without 2 columns, or with an odd number of rows.
    r = RLI()
    seed = rl_seed(11, 3)
    for T in RL_TYPES
        epsT = rl_eps(T)
        big1 = T(1) + 1000 * T(epsT)
        h6, h8 = T(0.6), T(0.8)
        F(a...) = SMatrix{4, 2, T}(a...)
        cases = ((F(2, 0, 0, 0, 0, 0, 1, 0), 1), (F(1, 0, 0, 0, 0, 0, 3, 0), 2), (F(1, 0, 0, 0, 0, 0, 0, 0), 2), (F(0, 0, 0, 0, 0, 0, 1, 0), 1),
                 (F(big1, 0, 0, 0, 0, 0, big1, 0), 1), (F(1, 0, 0, 0, 0, 0, -big1, 0), 2),
                 (F(h6, 0, h8, 0, 0, 0, 1, 0), 1),                    # column 1 supported in both blocks
                 (F(1, 0, 0, 0, h6, 0, h8, 0), 2),                    # column 2 supported in both blocks
                 (F(1, 0, 0, 0, 0, 0, 0, 1), 2),                      # unit columns in their own blocks, block vectors (1, 0) and (0, 1) differ
                 (F(0, 1, 0, 0, 0, 0, 1, 0), 2))                      # block vectors (0, 1) and (1, 0) differ
        for (W, col) in cases
            e = rl_health_error(() -> check_frame(r, W; orbit = 7))
            @test e isa NumericalHealthError && e.what == "initial frame" && e.orbit == 7 && e.step == 0 && e.backend === :cpu && e.column == col
            @test e isa NumericalHealthError && e.value > e.threshold >= 0
        end
        for W in (F(1, 0, 0, 0, 0, 0, 1, 0), F(0, 1, 0, 0, 0, 0, 0, 1), F(0, -1, 0, 0, 0, 0, 0, -1), F(1 + T(epsT), 0, 0, 0, 0, 0, 1 + T(epsT), 0))
            delta = nothing
            @test rl_health_error(() -> (delta = check_frame(r, W; orbit = 3))) === nothing
            @test delta !== nothing && 0 <= delta <= 2 * epsT
        end
        @test_throws ArgumentError check_frame(r, SMatrix{4, 3, T}(1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0))
        @test_throws ArgumentError check_frame(r, SMatrix{3, 2, T}(1, 0, 0, 0, 0, 1))
        nid = T === Float64 ? 10_000 : 1_000
        rejected = count(j -> rl_health_error(() -> check_frame(r, rl_frame(rl_xi0(seed, j, T)); orbit = j)) !== nothing, 1:nid)
        @test rejected == 0
    end
end

@testset "AC11c rejection before any work" begin
    _ = (RLI, check_frame, PairMap)
    # Sali Amendment D67 form: the frame check runs before the stream, the tangent-core stream is unchanged. The caller checks every supplied frame first, so a rejected frame never
    # reaches the map (rl_Boom throws if forward or tangent is called); an accepted frame proceeds to the map.
    seed = rl_seed(11, 4)
    good = rl_frame(rl_xi0(seed, 1, Float64))
    bad = SMatrix{4, 2, Float64}(2, 0, 0, 0, 0, 0, 1, 0)
    W0 = [good, bad]
    z0 = SVector(0.0, 0.0, 0.0, 0.0)
    checked_run = () -> begin
        foreach(j -> check_frame(RLI(), W0[j]; orbit = j), eachindex(W0))
        propagate_orbit(PairMap(rl_Boom()), z0, W0[2], 4, 2, (RLI(),))
    end
    err = rl_health_error(checked_run)
    @test err isa NumericalHealthError && err.what == "initial frame" && err.orbit == 2 && err.step == 0
    err = rl_health_error(() -> check_frame(RLI(), bad; orbit = 5))
    @test err isa NumericalHealthError && err.what == "initial frame" && err.orbit == 5 && err.step == 0
    @test rl_health_error(() -> check_frame(RLI(), good; orbit = 1)) === nothing
    @test_throws ErrorException propagate_orbit(PairMap(rl_Boom()), z0, good, 4, 2, (RLI(),))
end

@testset "AC11d overflow" begin
    _ = (RLI, PairMap)
    # A = diag(mu, 1/mu), mu = e^3, Float64, chains {RLI} and {RLI (non-replacing), MEGNO}: tau = 200 raises no error at n = 200 (the columns arrive at e^600 and RLI reads them before
    # they are normalised; one handoff, the two blocks run the same operations so DL = RLI = 0 and s_A = s_B); tau = 300 raises the stream's error at handoff 300 (e^900 overflows),
    # as sali 13(d).
    mu = exp(3.0)
    sys = LinearSymplecticMap(SMatrix{2, 2, Float64}(mu, 0, 0, 1 / mu))
    xi = rl_xi0(rl_seed(11, 5), 1, Float64)
    W0 = rl_frame(xi)
    z0 = SVector(0.0, 0.0, 0.0, 0.0)                                   # both orbits at the fixed point: the state stays finite, only the deviation matrix overflows
    for chain in ((RLI(),), (RLI(replace = false), MEGNO()))
        @test rl_health_error(() -> propagate_orbit(PairMap(sys), z0, W0, 200, 200, chain)) === nothing
        _, _, st = propagate_orbit(PairMap(sys), z0, W0, 200, 200, chain)
        f = reducer_final(RLI(), st[1])
        @test f.n == 200 && f.handoffs == 1
        @test f.DL == 0 && f.RLI == 0 && f.sum[1] == f.sum[2] && f.comp[1] == f.comp[2]
        @test abs(abs(f.sum[1] + f.comp[1]) - 600) < 20
        err = rl_health_error(() -> propagate_orbit(PairMap(sys), z0, W0, 300, 300, chain))
        @test err isa NumericalHealthError && err.what == "deviation matrix" && err.step == 300 && err.backend === :cpu
    end
end

@testset "AC11e run record" begin
    _ = (RLI, PairMap, pair_start)
    # The run record carries tau, precision and k (stream), the indicator name and the convention tag (reducer_provenance). In a chain the later reducer wins a shared key (D34 A2), so
    # the indicator name is tested on {RLI} alone; Dx, Dx_r, the state mode and the pair id come from pair_start (the record of the start); the direction rule, the seed, the frame
    # source and the chain position are recorded by Experiments, not by this module (the convention tag names the direction rule).
    @test reducer_provenance(RLI()).indicator == "rli"
    @test reducer_provenance(RLI(replace = false)).indicator == "rli"
    @test reducer_provenance(RLI()).convention == RL_CONVENTION
    for wrapped in (true, false)
        sys = rl_sys(2.0, wrapped)
        xi = rl_xi0(rl_seed(11, 6), 1, Float64)
        rec = pair_start(sys, SVector(1.0, 2.0), xi, 1e-9; orbit = 1)
        @test all(k -> haskey(rec, k), (:z, :Dx, :Dx_r, :mode))
        @test rec.Dx == 1e-9 && 0.5e-9 <= rec.Dx_r <= 2e-9 && rec.mode == (wrapped ? "wrapped" : "lifted")
        res = propagate(PairMap(sys), [rec.z], 8; k = 2, tau = 2, seed = rl_seed(11, 6), reducers = (RLI(),), W0 = [rl_frame(xi)])
        @test res.provenance.indicator == "rli" && res.provenance.convention == RL_CONVENTION
        @test res.provenance.k == 2 && res.provenance.tau == 2 && res.provenance.precision == "Float64"
        f = res.finals[1][1]
        @test f.n == 8 && f.handoffs == 4 && f.sum isa SVector{2, Float64} && f.comp isa SVector{2, Float64}
        @test f.DL isa Float64 && f.RLI isa Float64 && f.Q isa Float64
    end
end

# ------------------------------------------------------------------ AC12

@testset "AC12a orbit and chain invariance" begin
    _ = (RLI, MEGNO, PairMap)
    # (a) The pair orbit at every handoff is bitwise equal for the chains {}, {RLI alone, the replacer}, {RLI (non-replacing), MEGNO}, {MEGNO (non-replacing), RLI} (the one-replacer
    # rule of D71, D72); (b) the RLI output of each chain is bitwise equal to its stand-alone value (the non-replacing position changes nothing), the final deviation matrix is the
    # same in all chains with a replacer, and the sums of RLI and MEGNO are bitwise equal per column; K = 2 pi, 5 starts, N = 400, tau in {1, 4}, Float64.
    T = Float64
    c = rl_case(:two_pi, T, 5)
    P = PairMap(c.sys)
    prs = rl_pairs(c, 1e-9)
    for tau in (1, 4), j in 1:5
        z, W = prs[j].z, rl_frame(c.xi0[j])
        x0_, _, st0 = propagate_orbit(P, z, W, 400, tau, (rl_Trace(),); orbit = j)
        x1, W1, st1 = propagate_orbit(P, z, W, 400, tau, (RLI(), rl_Trace()); orbit = j)
        x2, W2, st2 = propagate_orbit(P, z, W, 400, tau, (RLI(replace = false), MEGNO(), rl_Trace()); orbit = j)
        x3, W3, st3 = propagate_orbit(P, z, W, 400, tau, (MEGNO(replace = false), RLI(), rl_Trace()); orbit = j)
        xm, Wm, stm = propagate_orbit(P, z, W, 400, tau, (MEGNO(), rl_Trace()); orbit = j)
        tr = (reducer_final(rl_Trace(), st0[1]).xs, reducer_final(rl_Trace(), st1[2]).xs, reducer_final(rl_Trace(), st2[3]).xs, reducer_final(rl_Trace(), st3[3]).xs,
              reducer_final(rl_Trace(), stm[2]).xs)
        @test length(tr[1]) == 400 ÷ tau
        @test rl_nbad(t -> isequal(tr[1], t), tr) == 0
        @test x0_ === x1 && x0_ === x2 && x0_ === x3 && x0_ === xm
        f0 = reducer_final(RLI(), st1[1])
        @test isequal(f0, reducer_final(RLI(replace = false), st2[1])) && isequal(f0, reducer_final(RLI(), st3[2]))
        m0 = reducer_final(MEGNO(), stm[1])
        @test isequal(m0, reducer_final(MEGNO(), st2[2])) && isequal(m0, reducer_final(MEGNO(replace = false), st3[1]))
        @test isequal(W1, W2) && isequal(W1, W3) && isequal(W1, Wm)
        @test f0.sum == m0.sum && f0.comp == m0.comp && f0.n == 400 && f0.handoffs == 400 ÷ tau
    end
end

@testset "AC12b propagate equals the serial per-pair body" begin
    _ = (RLI, MEGNO, PairMap)
    # propagate equals the serial map of propagate_orbit bitwise (x, W, finals) for chunk 64 and 1000: 4096 pairs, K = 2 pi, wrapped, n = 1e3, Float64, chains {RLI} (tau in {1, 4})
    # and {RLI (non-replacing), MEGNO} (tau = 4). The thread count is the run's (-t 2 under Pkg.test), never changed here; 1 versus N threads in child processes is a stage 04
    # check (test plan, criterion 12(b)).
    seed = rl_seed(12, 10)
    T = Float64
    m, N = 4096, 1000
    sys = rl_sys(rl_K(T, :two_pi), true)
    P = PairMap(sys)
    starts = [rl_cast(T, s) for s in rl_starts(StableRNG(seed + 2), m)]
    xis = [rl_xi0(seed, j, T) for j in 1:m]
    z0 = [rl_pair(sys, starts[j], xis[j], 1e-9).z for j in 1:m]
    W0 = [rl_frame(x) for x in xis]
    for (reducers, taus) in (((RLI(),), (1, 4)), ((RLI(replace = false), MEGNO()), (4,)))
        for tau in taus
            serial = [propagate_orbit(P, z0[j], W0[j], N, tau, reducers; orbit = j) for j in 1:m]
            xs, Ws = [s[1] for s in serial], [s[2] for s in serial]
            fs = [map(reducer_final, reducers, s[3]) for s in serial]
            for chunk in (64, 1000)
                res = propagate(P, z0, N; k = 2, tau = tau, seed = seed, reducers = reducers, chunk = chunk, W0 = W0)
                same_x, same_W, same_f = isequal(res.x, xs), isequal(res.W, Ws), isequal(res.finals, fs)      # booleans: a failure report must not dump 4096 pairs
                @test same_x
                @test same_W
                @test same_f
            end
        end
    end
end

@testset "AC12c allocation" begin
    _ = (RLI, MEGNO, PairMap)
    # 0 bytes per pair after warm-up, n = 1e4, both types, tau in {1, 4}, chains {RLI} and {RLI (non-replacing), MEGNO}, wrapped state.
    seed = rl_seed(12, 20)
    for T in RL_TYPES, tau in (1, 4), reducers in ((RLI(),), (RLI(replace = false), MEGNO()))
        sys = rl_sys(rl_K(T, :two_pi), true)
        xi = rl_xi0(seed, 1, T)
        pr = rl_pair(sys, rl_cast(T, SVector(1.0, 2.0)), xi, 1e-9)
        @test rl_alloc_bytes(PairMap(sys), pr.z, rl_frame(xi), 10_000, tau, reducers) == 0
    end
end

@testset "AC12d position in the ensemble" begin
    _ = (RLI, PairMap)
    # I8: the series of pair j depends only on (system, parameters, initial point, Dx, direction rule, seed, j, tau, N, state mode, frame), not on its position in the ensemble: the
    # ensemble run with permuted positions (the supplied frames follow their pair ids) gives bitwise the same finals, points and matrices per pair; K = 2 pi, 256 pairs, n = 400,
    # chunk 64, tau in {1, 4}.
    seed = rl_seed(12, 30)
    c = rl_case(:two_pi, Float64, 256)
    P = PairMap(c.sys)
    prs = rl_pairs(c, 1e-9)
    z0, W0 = [p.z for p in prs], rl_frames(c)
    perm = sortperm(rand(StableRNG(seed), 256))              # a seeded permutation
    for tau in (1, 4)
        a = propagate(P, z0, 400; k = 2, tau = tau, seed = seed, reducers = (RLI(),), chunk = 64, W0 = W0)
        b = propagate(P, z0[perm], 400; k = 2, tau = tau, seed = seed, reducers = (RLI(),), chunk = 64, W0 = W0[perm])
        @test isequal(a.finals[perm], b.finals) && isequal(a.x[perm], b.x) && isequal(a.W[perm], b.W)
    end
end

# ------------------------------------------------------------------ AC13

@testset "AC13 wrapped versus lifted state" begin
    _ = (RLI, PairMap)
    # (a) K = 0, dyadic starts: RLI = 0 in both modes and the wrapped (state length 8, D28) and the lifted (length 4) outputs bitwise identical (all handoffs, the whole final), tolerance 0,
    # tau in {1, 3}, n = 300, 50 starts. (b) regular class K = 0.5 (S7), n in {10, 100}, Dx in {1e-8, 1e-12}, Float64 (50 starts): |RLI_wrapped - RLI_lifted| <=
    # 2 * 64 K_std M_state n eps rho_max, rho_max the largest over the handoffs <= n and both orbits of the rho of L1 criterion 7 (the L1 criterion 11(b) form per vector with the
    # factor 2 of RL11 for the difference of two orbits; borrowed, not derived), M_state the larger of the two orbits.
    seed = rl_seed(13)
    rng = StableRNG(seed)
    dy() = rand(rng, 0:2048) / 1024
    T = Float64
    for tau in (1, 3)
        xs = [SVector{2, T}(T(dy()), T(dy())) for _ in 1:50]
        xis = [rl_xi0(seed, j, T) for j in 1:50]
        res = map((true, false)) do wrapped
            sys = rl_sys(exact(T, 0), wrapped)
            rl_run((RLI(),), PairMap(sys), [rl_pair(sys, xs[j], xis[j], 2.0^-30).z for j in 1:50], [rl_frame(x) for x in xis], 300, tau; all = true)
        end
        @test rl_nbad(j -> isequal(res[1][j][1].snaps, res[2][j][1].snaps), 1:50) == 0
        @test rl_nbad(j -> all(f -> f.DL == 0 && f.RLI == 0, res[1][j][1].snaps), 1:50) == 0
    end
    m = 50
    cw, cl = rl_case(:half, T, m; wrapped = true), rl_case(:half, T, m; wrapped = false)
    for Dx in (1e-8, 1e-12), n in (10, 100)
        rw = rl_run((RLI(),), PairMap(cw.sys), [p.z for p in rl_pairs(cw, Dx)], rl_frames(cw), n, 1; steps = [n])
        prl = rl_pairs(cl, Dx)
        rl_ = rl_run((RLI(),), PairMap(cl.sys), [p.z for p in prl], rl_frames(cl), n, 1; steps = [n])
        bad = 0
        for j in 1:m
            a, b = rw[j][1].snaps[1], rl_[j][1].snaps[1]
            ra, Ma = rl_rho_max(0.5, cl.x0[j], cl.xi0[j], n)
            rb, Mb = rl_rho_max(0.5, prl[j].zB, cl.xi0[j], n)
            B = 2 * 64 * 0.5 * max(Ma, Mb) * n * eps(Float64) * max(ra, rb)
            bad += !(abs(a.RLI - b.RLI) <= B)
        end
        @test bad == 0
    end
end
