# AC11 to AC13: health and frames, chain and execution invariance (the one-replacer rule of D71, the replace option of D72), wrapped versus lifted state.
# AC14 (ChaosTools one-column tangent oracle line) and AC15 (local GPU) are stage 04; the test plan maps them.

# The chains with the one replacer last (D71): MEGNO alone replaces; MEGNO reads before SALI; FLI and MEGNO read before SALI.
mg_chains() = ((MEGNO(),), (MEGNO(replace = false), SALI()), (FLI(replace = false), MEGNO(replace = false), SALI()))

@testset "AC11a zero tangent, column and orbit of the health error" begin
    _ = (MEGNO, FLI, SALI)
    # A zero norm is not caught by the finiteness check of the stream: the reducer records it. Tangent zero for flagged orbits: error at the first handoff, kind
    # "deviation matrix", orbit id, step, backend, column 1, value 0, nothing stored. Column 2 collapsing alone (tangent diag(1, 0), frame I): column 2. Float64 through
    # propagate (kernel path, chunk 2), the three chains; Double64 through the per-orbit body. MEGNO needs k = 2.
    seed = mg_seed(11, 1)
    x0s(flagged) = [SVector{3, Float64}(0, 0, j in flagged ? 1 : 0) for j in 1:6]
    zero_map, col2_map = mg_Fault(:zero, 0), mg_Fault(:col2, 0)
    for chain in mg_chains()
        ok = propagate(zero_map, x0s(()), 12; k = 2, tau = 3, seed = seed, reducers = chain, chunk = 2)
        @test length(ok.finals) == 6
        res = nothing
        err = mg_health_error(() -> (res = propagate(zero_map, x0s((4,)), 12; k = 2, tau = 3, seed = seed, reducers = chain, chunk = 2)))
        @test err isa NumericalHealthError
        @test res === nothing
        @test err.what == "deviation matrix" && err.orbit == 4 && err.step == 3 && err.backend === :cpu
        @test err.column == 1 && err.value == 0.0
        res = nothing
        err = mg_health_error(() -> (res = propagate(zero_map, x0s((5, 2)), 12; k = 2, tau = 3, seed = seed, reducers = chain, chunk = 2)))
        @test err isa NumericalHealthError && err.orbit == 2 && err.step == 3 && res === nothing        # lowest orbit id first
        eye = [SMatrix{2, 2, Float64}(1, 0, 0, 1) for _ in 1:6]
        err = mg_health_error(() -> propagate(col2_map, x0s((5,)), 12; k = 2, tau = 3, seed = seed, reducers = chain, chunk = 2, W0 = eye))
        @test err isa NumericalHealthError && err.what == "deviation matrix" && err.orbit == 5 && err.step == 3
        @test err.column == 2 && err.value == 0.0
    end
    x = SVector{3, Double64}(0, 0, 1)
    W0 = mg_frame(seed, 1, Double64)
    err = mg_health_error(() -> propagate_orbit(zero_map, x, W0, 12, 3, (MEGNO(),); orbit = 4))
    @test err isa NumericalHealthError && err.what == "deviation matrix" && err.orbit == 4 && err.step == 3
    @test err.backend === :cpu && err.column == 1 && err.value == 0.0
    err = mg_health_error(() -> propagate_orbit(col2_map, x, SMatrix{2, 2, Double64}(1, 0, 0, 1), 12, 3, (MEGNO(),); orbit = 4))
    @test err isa NumericalHealthError && err.column == 2 && err.step == 3
    # MEGNO needs k = 2: ArgumentError before any work (mg_Boom throws if the map is called).
    e = try
        propagate(mg_Boom(), [SVector(0.0, 0.0)], 4; k = 1, tau = 2, seed = seed, reducers = (MEGNO(),))
        nothing
    catch ex
        ex
    end
    @test e isa ArgumentError && occursin("MEGNO needs k = 2, got k = 1", e.msg)
end

@testset "AC11b non-finite deviation matrix is not masked" begin
    _ = (MEGNO, FLI, SALI)
    # Non-finite entries are tangent-core AC8; only that the MEGNO chains do not change that verdict (kind "deviation matrix", column 0, threshold Inf).
    seed = mg_seed(11, 2)
    nan_map = mg_Fault(:nan, 3)
    x0s(flagged) = [SVector{3, Float64}(0, 0, j in flagged ? 1 : 0) for j in 1:6]
    for chain in mg_chains()
        res = nothing
        err = mg_health_error(() -> (res = propagate(nan_map, x0s((3,)), 20; k = 2, tau = 5, seed = seed, reducers = chain, chunk = 2)))
        @test err isa NumericalHealthError && res === nothing
        @test err.what == "deviation matrix" && err.orbit == 3 && err.step == 5 && err.backend === :cpu
        @test err.column == 0 && err.threshold == Inf && !isfinite(err.value)
    end
    x = SVector{3, Double64}(0, 0, 1)
    err = mg_health_error(() -> propagate_orbit(nan_map, x, mg_frame(seed, 1, Double64), 20, 5, (MEGNO(),); orbit = 3))
    @test err isa NumericalHealthError && err.what == "deviation matrix" && err.orbit == 3 && err.step == 5 && err.column == 0
end

@testset "AC11c supplied frame check" begin
    _ = (MEGNO, check_frame)
    # A supplied frame is accepted only if both columns are unit vectors to roundoff (the convention s(0) = 0 needs it; independence and orthogonality are not needed by
    # MEGNO and not checked here, the caller reports orthonormality): a column of norm 2, 3, 0 or 1 + 1000 eps is rejected with a numerical health error of kind
    # "initial frame" (step 0, the orbit id, the offending column, a value above a threshold that lies in (0, 1e-6)); unit columns, also dependent ones, and a column of
    # norm 1 + eps are accepted and check_frame returns the larger norm defect (at most 2 eps); the seeded frame is accepted for every orbit id 1..1e4 (Double64: 1..1000).
    r = MEGNO()
    seed = mg_seed(11, 3)
    for T in MG_TYPES
        epsT = mg_eps(T)
        big1 = T(1) + 1000 * T(epsT)
        for (W, col) in ((SMatrix{2, 2, T}(2, 0, 0, 1), 1), (SMatrix{2, 2, T}(1, 0, 0, 3), 2), (SMatrix{2, 2, T}(1, 0, 0, 0), 2),
                         (SMatrix{2, 2, T}(0, 0, 0, 1), 1), (SMatrix{2, 2, T}(big1, 0, 0, 1), 1), (SMatrix{2, 2, T}(1, 0, 0, -big1), 2))
            e = mg_health_error(() -> check_frame(r, W; orbit = 7))
            @test e isa NumericalHealthError && e.what == "initial frame" && e.orbit == 7 && e.step == 0 && e.backend === :cpu && e.column == col
            @test e isa NumericalHealthError && 0 < e.threshold < 1e-6 && e.value > e.threshold
        end
        for W in (SMatrix{2, 2, T}(1, 0, 0, 1), SMatrix{2, 2, T}(1, 0, 1, 0), SMatrix{2, 2, T}(0, 1, 0, -1), SMatrix{2, 2, T}(1 + T(epsT), 0, 0, 1))
            delta = nothing
            @test mg_health_error(() -> (delta = check_frame(r, W; orbit = 3))) === nothing
            @test delta !== nothing && 0 <= delta <= 2 * epsT
        end
        @test_throws ArgumentError check_frame(r, SMatrix{2, 3, T}(1, 0, 0, 1, 0, 0))
        nid = T === Float64 ? 10_000 : 1_000
        rejected = count(j -> mg_health_error(() -> check_frame(r, mg_frame(seed, j, T); orbit = j)) !== nothing, 1:nid)
        @test rejected == 0
    end
end

@testset "AC11c rejection before any work" begin
    _ = (MEGNO, check_frame)
    # Sali Amendment D67 form: the frame check runs before the stream, the tangent-core stream is unchanged. The caller checks every supplied frame first, so a rejected
    # frame never reaches the map (mg_Boom throws if forward or tangent is called); an accepted frame proceeds to the map.
    seed = mg_seed(11, 4)
    bad = SMatrix{2, 2, Float64}(2, 0, 0, 1)
    W0 = [mg_frame(seed, 1, Float64), bad]
    checked_run = () -> begin
        foreach(j -> check_frame(MEGNO(), W0[j]; orbit = j), eachindex(W0))
        propagate_orbit(mg_Boom(), SVector(0.0, 0.0), W0[2], 4, 2, (MEGNO(),))
    end
    err = mg_health_error(checked_run)
    @test err isa NumericalHealthError && err.what == "initial frame" && err.orbit == 2 && err.step == 0
    err = mg_health_error(() -> check_frame(MEGNO(), bad; orbit = 5))
    @test err isa NumericalHealthError && err.what == "initial frame" && err.orbit == 5 && err.step == 0
    @test mg_health_error(() -> check_frame(MEGNO(), W0[1]; orbit = 1)) === nothing
    @test_throws ErrorException propagate_orbit(mg_Boom(), SVector(0.0, 0.0), W0[1], 4, 2, (MEGNO(),))
end

@testset "AC11d overflow" begin
    _ = (MEGNO, SALI)
    # A = diag(mu, 1/mu), mu = e^3, Float64, chains {MEGNO} and {MEGNO (not the replacer), SALI}: tau = 200 raises no error at n = 200 (the columns arrive at e^600 and MEGNO reads
    # them before they are normalised; one handoff, so Y = Ybar = 2 s and Ybar20 = 3 s); tau = 300 raises the stream's error at handoff 300 (e^900 overflows), as sali 13(d).
    mu = exp(3.0)
    sys = LinearSymplecticMap(SMatrix{2, 2, Float64}(mu, 0, 0, 1 / mu))
    W0 = mg_frame(mg_seed(11, 5), 1, Float64)
    x0 = SVector(0.0, 0.0)
    for chain in ((MEGNO(),), (MEGNO(replace = false), SALI()))
        @test mg_health_error(() -> propagate_orbit(sys, x0, W0, 200, 200, chain)) === nothing
        _, _, st = propagate_orbit(sys, x0, W0, 200, 200, chain)
        f = reducer_final(MEGNO(), st[1])
        @test f.n == 200 && f.handoffs == 1
        for i in 1:2
            S = f.sum[i] + f.comp[i]
            @test abs(abs(S) - 600) < 20 && isapprox(f.Y[i], 2 * S; rtol = 1e-13) && isapprox(f.Ybar[i], 2 * S; rtol = 1e-13) && isapprox(f.Ybar20[i], 3 * S; rtol = 1e-13)
        end
        err = mg_health_error(() -> propagate_orbit(sys, x0, W0, 300, 300, chain))
        @test err isa NumericalHealthError && err.what == "deviation matrix" && err.step == 300 && err.backend === :cpu
    end
end

@testset "AC11e run record" begin
    _ = (MEGNO,)
    # The run record carries tau, precision and k (stream), the indicator name and the convention tag (reducer_provenance). In a chain the later reducer wins a shared key (D34 A2),
    # so the indicator name is tested on {MEGNO} alone; the frame source and the chain position are recorded by Experiments, not by this module.
    @test reducer_provenance(MEGNO()).indicator == "megno"
    @test reducer_provenance(MEGNO(replace = false)).indicator == "megno"
    @test reducer_provenance(MEGNO()).convention == MG_CONVENTION
    sys = StandardMap(2.0)
    res = propagate(sys, [SVector(1.0, 2.0)], 8; k = 2, tau = 2, seed = mg_seed(11, 6), reducers = (MEGNO(),))
    @test res.provenance.indicator == "megno" && res.provenance.convention == MG_CONVENTION
    @test res.provenance.k == 2 && res.provenance.tau == 2 && res.provenance.precision == "Float64"
    f = res.finals[1][1]
    @test f.n == 8 && f.handoffs == 4 && f.sum isa SVector{2, Float64} && f.comp isa SVector{2, Float64}
    @test f.Y isa SVector{2, Float64} && f.Ybar isa SVector{2, Float64} && f.Ybar20 isa SVector{2, Float64}
    @test f.R1 isa SVector{2, Float64} && f.R2 isa SVector{2, Float64} && f.Q isa SVector{2, Float64} && f.Q2 isa SVector{2, Float64}
end

# ------------------------------------------------------------------ AC12

@testset "AC12a orbit and chain invariance" begin
    _ = (MEGNO, FLI, SALI, Benettin)
    # (a) the orbit at every handoff is bitwise equal for the chains {}, {L1} on k = 1, {MEGNO alone, the replacer} on k = 2, {MEGNO (not the replacer), SALI},
    # {FLI (not the replacer), MEGNO (not the replacer), SALI}, {MEGNO (not the replacer), FLI (not the replacer), SALI}, {SALI} and {} on k = 2 (L1 S2); (b) the MEGNO, FLI and
    # SALI outputs of each chain are bitwise equal to their stand-alone values (the non-replacing position changes nothing), the final deviation matrix is the same in all chains
    # with a replacer, and the sums of the three reducers (and of the L1 runs of each column) are bitwise equal on the identical orbit; K = 2 pi, 5 starts, N = 400, tau in {1, 4},
    # Float64.
    seed = mg_seed(12)
    T = Float64
    for tau in (1, 4)
        sys = StandardMap(mg_K(T, :two_pi))
        for (j, s) in enumerate(mg_starts(StableRNG(seed + 1), 5))
            x0 = mg_cast(T, s)
            W2 = mg_frame(seed, j, T)
            W1 = SMatrix{2, 1, T}(W2[:, 1])
            xa, _, sta = propagate_orbit(sys, x0, W1, 400, tau, (mg_Trace(),); orbit = j)
            xb, _, stb = propagate_orbit(sys, x0, W1, 400, tau, (Benettin(), mg_Trace()); orbit = j)
            xc, Wc, stc = propagate_orbit(sys, x0, W2, 400, tau, (MEGNO(), mg_Trace()); orbit = j)
            xd, Wd, std = propagate_orbit(sys, x0, W2, 400, tau, (MEGNO(replace = false), SALI(), mg_Trace()); orbit = j)
            xe, We, ste = propagate_orbit(sys, x0, W2, 400, tau, (FLI(replace = false), MEGNO(replace = false), SALI(), mg_Trace()); orbit = j)
            xf, Wf, stf = propagate_orbit(sys, x0, W2, 400, tau, (MEGNO(replace = false), FLI(replace = false), SALI(), mg_Trace()); orbit = j)
            xg, Wg, stg = propagate_orbit(sys, x0, W2, 400, tau, (SALI(), mg_Trace()); orbit = j)
            xh, _, sth = propagate_orbit(sys, x0, W2, 400, tau, (mg_Trace(),); orbit = j)
            xi, Wi, sti = propagate_orbit(sys, x0, W2, 400, tau, (FLI(), mg_Trace()); orbit = j)
            tra = reducer_final(mg_Trace(), sta[1]).xs
            trs = (reducer_final(mg_Trace(), stb[2]).xs, reducer_final(mg_Trace(), stc[2]).xs, reducer_final(mg_Trace(), std[3]).xs, reducer_final(mg_Trace(), ste[4]).xs,
                   reducer_final(mg_Trace(), stf[4]).xs, reducer_final(mg_Trace(), stg[2]).xs, reducer_final(mg_Trace(), sth[1]).xs)
            @test length(tra) == 400 ÷ tau
            @test mg_nbad(tr -> isequal(tra, tr), trs) == 0
            @test xa === xb && xa === xc && xa === xd && xa === xe && xa === xf && xa === xg && xa === xh
            # (b)
            m0 = reducer_final(MEGNO(), stc[1])
            @test isequal(m0, reducer_final(MEGNO(replace = false), std[1])) && isequal(m0, reducer_final(MEGNO(replace = false), ste[2])) &&
                  isequal(m0, reducer_final(MEGNO(replace = false), stf[1]))                           # MEGNO: stand-alone (replacer) equals the non-replacing positions
            f0 = reducer_final(FLI(), sti[1])
            @test isequal(f0, reducer_final(FLI(replace = false), ste[1])) && isequal(f0, reducer_final(FLI(replace = false), stf[2]))
            s0 = reducer_final(SALI(), stg[1])
            @test isequal(s0, reducer_final(SALI(), std[2])) && isequal(s0, reducer_final(SALI(), ste[3])) && isequal(s0, reducer_final(SALI(), stf[3]))
            @test isequal(Wc, Wg) && isequal(Wd, Wg) && isequal(We, Wg) && isequal(Wf, Wg) && isequal(Wi, Wg)    # the same unit columns in all chains with a replacer
            # the sums: bitwise equal to SALI, FLI and the L1 run of each column
            @test m0.sum == s0.sum && m0.comp == s0.comp && m0.sum == f0.sum && m0.comp == f0.comp
            @test m0.n == 400 && m0.handoffs == 400 ÷ tau
            for i in 1:2
                _, _, stl = propagate_orbit(sys, x0, SMatrix{2, 1, T}(W2[:, i]), 400, tau, (Benettin(),); orbit = j)
                l = reducer_final(Benettin(), stl[1])
                @test m0.sum[i] + m0.comp[i] == l.sum[1] + l.comp[1]
            end
        end
    end
end

@testset "AC12b propagate equals the serial per-orbit body" begin
    _ = (MEGNO, FLI, SALI)
    # propagate equals the serial map of propagate_orbit bitwise (x, W, finals) for chunk 64 and 1000: 4096 orbits, K = 2 pi, n = 1e3, Float64, chains {MEGNO (not the replacer),
    # SALI} (tau in {1, 4}), {FLI (not the replacer), MEGNO (not the replacer), SALI} and {MEGNO} (tau = 4). The thread count is the run's (-t 2 under Pkg.test), never changed
    # here; 1 versus N threads in child processes is a stage 04 check (test plan, criterion 12(c)).
    seed = mg_seed(12, 10)
    T = Float64
    for (reducers, taus) in (((MEGNO(replace = false), SALI()), (1, 4)), ((FLI(replace = false), MEGNO(replace = false), SALI()), (4,)), ((MEGNO(),), (4,)))
        for tau in taus
            sys = StandardMap(mg_K(T, :two_pi))
            m, N = 4096, 1000
            x0 = [mg_cast(T, s) for s in mg_starts(StableRNG(seed + 2), m)]
            W0 = [mg_frame(seed, j, T) for j in 1:m]
            serial = [propagate_orbit(sys, x0[j], W0[j], N, tau, reducers; orbit = j) for j in 1:m]
            xs, Ws = [s[1] for s in serial], [s[2] for s in serial]
            fs = [map(reducer_final, reducers, s[3]) for s in serial]
            for chunk in (64, 1000)
                res = propagate(sys, x0, N; k = 2, tau = tau, seed = seed, reducers = reducers, chunk = chunk, W0 = W0)
                same_x, same_W, same_f = isequal(res.x, xs), isequal(res.W, Ws), isequal(res.finals, fs)      # booleans: a failure report must not dump 4096 orbits
                @test same_x
                @test same_W
                @test same_f
            end
        end
    end
end

@testset "AC12c allocation" begin
    _ = (MEGNO, SALI)
    # 0 bytes per orbit after warm-up, n = 1e4, both types, tau in {1, 4}, chains {MEGNO} and {MEGNO (not the replacer), SALI}.
    seed = mg_seed(12, 20)
    for T in MG_TYPES, tau in (1, 4), reducers in ((MEGNO(),), (MEGNO(replace = false), SALI()))
        sys = StandardMap(mg_K(T, :two_pi))
        @test mg_alloc_bytes(sys, mg_cast(T, SVector(1.0, 2.0)), mg_frame(seed, 1, T), 10_000, tau, reducers) == 0
    end
end

@testset "AC12d column exchange, sign change" begin
    _ = (MEGNO,)
    # I4: exchanging the two columns of a supplied frame exchanges every per-column output of the final (Y, Ybar, Ybar20, sums, compensation, the weighted sums) bitwise; a sign
    # change of either column leaves every output bitwise unchanged. K = 2 pi and the K = 0.5 regular class, 5 starts, N = 400, tau in {1, 4}, every handoff, Float64 and Double64.
    seed = mg_seed(12, 30)
    flds = (:Y, :Ybar, :Ybar20, :sum, :comp, :R1, :R2, :Q, :Q2)
    for T in MG_TYPES, label in (:two_pi, :half), tau in (1, 4)
        c = mg_case(label, T, 5)
        W0 = c.W0
        swap(W) = SMatrix{2, 2, T}(W[1, 2], W[2, 2], W[1, 1], W[2, 1])
        neg1(W) = SMatrix{2, 2, T}(-W[1, 1], -W[2, 1], W[1, 2], W[2, 2])
        neg2(W) = SMatrix{2, 2, T}(W[1, 1], W[2, 1], -W[1, 2], -W[2, 2])
        base = mg_run((MEGNO(),), c.sys, c.x0, W0, 400, tau; all = true)
        sw = mg_run((MEGNO(),), c.sys, c.x0, swap.(W0), 400, tau; all = true)
        n1 = mg_run((MEGNO(),), c.sys, c.x0, neg1.(W0), 400, tau; all = true)
        n2 = mg_run((MEGNO(),), c.sys, c.x0, neg2.(W0), 400, tau; all = true)
        for j in 1:5
            a, b = base[j][1].snaps, sw[j][1].snaps
            @test length(a) == 400 ÷ tau
            for f in flds
                ok = isequal([reverse(getproperty(s, f)) for s in a], [getproperty(s, f) for s in b])
                @test ok
            end
            for other in (n1, n2)
                ok = isequal(a, other[j][1].snaps)
                @test ok
            end
        end
    end
end

# ------------------------------------------------------------------ AC13

@testset "AC13 wrapped versus lifted state" begin
    _ = (MEGNO,)
    # (a) K = 0, dyadic starts: the wrapped (state length 4, D28) and the lifted series bitwise identical (all handoffs, the whole final), tolerance 0, tau in {1, 3}.
    # (b) regular class K = 0.5 (S7), n in {10, 100}, Float64 (50 starts), both columns, the three series: |M_wrapped - M_lifted| <= 6 * 64 K_std M_state n^2 eps rho_max, rho_max
    # the largest over the handoffs <= n of the rho of L1 criterion 7 for the column (the L1 criterion 11(b) form per vector with the coefficient sum of M11; borrowed, not derived
    # for weighted sums); comparisons with rho_max > 50 are excluded (at most 10 percent of the 100 comparisons per n).
    seed = mg_seed(13)
    rng = StableRNG(seed)
    dy() = rand(rng, 0:2048) / 1024
    T = Float64
    for tau in (1, 3)
        zs = [SVector{2, T}(T(dy()), T(dy())) for _ in 1:50]
        W0 = [mg_frame(seed, j, T) for j in 1:50]
        rl = mg_run((MEGNO(),), StandardMap(exact(T, 0)), zs, W0, 300, tau; all = true)
        rw = mg_run((MEGNO(),), WrappedStandardMap(exact(T, 0)), torus_state.(zs), W0, 300, tau; all = true)
        @test mg_nbad(j -> isequal(rl[j][1].snaps, rw[j][1].snaps), 1:50) == 0
    end
    m = 50
    c = mg_case(:half, T, m)
    for n in (10, 100)
        rl = mg_run((MEGNO(),), c.sys, c.x0, c.W0, n, 1; steps = [n])
        rw = mg_run((MEGNO(),), WrappedStandardMap(0.5), torus_state.(c.x0), c.W0, n, 1; steps = [n])
        excluded = 0
        for j in 1:m, i in 1:2
            a, b = rl[j][1].snaps[1], rw[j][1].snaps[1]
            rmax, Mst = mg_rho_max(0.5, c.x0[j], c.W0[j][:, i], n)
            if rmax > 50
                excluded += 1
                continue
            end
            B = 6 * 64 * 0.5 * Mst * n^2 * eps(Float64) * rmax
            @test abs(a.Y[i] - b.Y[i]) <= B && abs(a.Ybar[i] - b.Ybar[i]) <= B && abs(a.Ybar20[i] - b.Ybar20[i]) <= B
        end
        @test excluded <= m * 2 ÷ 10
    end
end
