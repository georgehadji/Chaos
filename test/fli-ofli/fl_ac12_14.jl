# AC12 to AC14: health and frames, chain and execution invariance (the one-replacer rule of D71), wrapped versus lifted state.
# AC15 (ChaosTools gali oracle line) and AC16 (local GPU) are stage 04; the test plan maps them.

@testset "AC12a zero tangent, column and orbit of the health error" begin
    _ = (FLI, SALI)
    # A zero norm is not caught by the finiteness check of the stream: the reducer records it. Tangent zero for flagged orbits: error at the first handoff, kind
    # "deviation matrix", orbit id, step, backend, column 1, value 0, nothing stored. Column 2 collapsing alone (tangent diag(1, 0), frame I): column 2. Float64 through
    # propagate (kernel path, chunk 2), chains {FLI} and {FLI (not the replacer), SALI}; Double64 through the per-orbit body. FLI needs k = 2.
    seed = fl_seed(12, 1)
    x0s(flagged) = [SVector{3, Float64}(0, 0, j in flagged ? 1 : 0) for j in 1:6]
    zero_map, col2_map = fl_Fault(:zero, 0), fl_Fault(:col2, 0)
    for chain in ((FLI(),), (FLI(replace = false), SALI()))
        ok = propagate(zero_map, x0s(()), 12; k = 2, tau = 3, seed = seed, reducers = chain, chunk = 2)
        @test length(ok.finals) == 6
        res = nothing
        err = fl_health_error(() -> (res = propagate(zero_map, x0s((4,)), 12; k = 2, tau = 3, seed = seed, reducers = chain, chunk = 2)))
        @test err isa NumericalHealthError
        @test res === nothing
        @test err.what == "deviation matrix" && err.orbit == 4 && err.step == 3 && err.backend === :cpu
        @test err.column == 1 && err.value == 0.0
        res = nothing
        err = fl_health_error(() -> (res = propagate(zero_map, x0s((5, 2)), 12; k = 2, tau = 3, seed = seed, reducers = chain, chunk = 2)))
        @test err isa NumericalHealthError && err.orbit == 2 && err.step == 3 && res === nothing        # lowest orbit id first
        eye = [SMatrix{2, 2, Float64}(1, 0, 0, 1) for _ in 1:6]
        err = fl_health_error(() -> propagate(col2_map, x0s((5,)), 12; k = 2, tau = 3, seed = seed, reducers = chain, chunk = 2, W0 = eye))
        @test err isa NumericalHealthError && err.what == "deviation matrix" && err.orbit == 5 && err.step == 3
        @test err.column == 2 && err.value == 0.0
    end
    x = SVector{3, Double64}(0, 0, 1)
    W0 = fl_frame(seed, 1, Double64)
    err = fl_health_error(() -> propagate_orbit(zero_map, x, W0, 12, 3, (FLI(),); orbit = 4))
    @test err isa NumericalHealthError && err.what == "deviation matrix" && err.orbit == 4 && err.step == 3
    @test err.backend === :cpu && err.column == 1 && err.value == 0.0
    err = fl_health_error(() -> propagate_orbit(col2_map, x, SMatrix{2, 2, Double64}(1, 0, 0, 1), 12, 3, (FLI(),); orbit = 4))
    @test err isa NumericalHealthError && err.column == 2 && err.step == 3
    # FLI needs k = 2: ArgumentError before any work (fl_Boom throws if the map is called).
    e = try
        propagate(fl_Boom(), [SVector(0.0, 0.0)], 4; k = 1, tau = 2, seed = seed, reducers = (FLI(),))
        nothing
    catch ex
        ex
    end
    @test e isa ArgumentError && occursin("FLI needs k = 2, got k = 1", e.msg)
end

@testset "AC12b non-finite deviation matrix is not masked" begin
    _ = (FLI, SALI)
    # Non-finite entries are tangent-core AC8; only that the FLI chains do not change that verdict (kind "deviation matrix", column 0, threshold Inf).
    seed = fl_seed(12, 2)
    nan_map = fl_Fault(:nan, 3)
    x0s(flagged) = [SVector{3, Float64}(0, 0, j in flagged ? 1 : 0) for j in 1:6]
    for chain in ((FLI(),), (FLI(replace = false), SALI()))
        res = nothing
        err = fl_health_error(() -> (res = propagate(nan_map, x0s((3,)), 20; k = 2, tau = 5, seed = seed, reducers = chain, chunk = 2)))
        @test err isa NumericalHealthError && res === nothing
        @test err.what == "deviation matrix" && err.orbit == 3 && err.step == 5 && err.backend === :cpu
        @test err.column == 0 && err.threshold == Inf && !isfinite(err.value)
    end
    x = SVector{3, Double64}(0, 0, 1)
    err = fl_health_error(() -> propagate_orbit(nan_map, x, fl_frame(seed, 1, Double64), 20, 5, (FLI(),); orbit = 3))
    @test err isa NumericalHealthError && err.what == "deviation matrix" && err.orbit == 3 && err.step == 5 && err.column == 0
end

@testset "AC12c supplied frame check" begin
    _ = (FLI, check_frame)
    # A supplied frame is accepted only if both columns are unit vectors to roundoff (the convention S(0) = 0 needs it; independence and orthogonality are not needed
    # by FLI and not checked here, the caller reports orthonormality): a column of norm 2, 3, 0 or 1 + 1000 eps is rejected with a numerical health error of kind
    # "initial frame" (step 0, the orbit id, the offending column, a value above a threshold that lies in (0, 1e-6)); unit columns, also dependent ones, and a column
    # of norm 1 + eps are accepted and check_frame returns the larger norm defect (at most 2 eps); the seeded frame is accepted for every orbit id 1..1e4 (Double64: 1..1000).
    r = FLI()
    seed = fl_seed(12, 3)
    for T in FL_TYPES
        epsT = fl_eps(T)
        big1 = T(1) + 1000 * T(epsT)
        for (W, col) in ((SMatrix{2, 2, T}(2, 0, 0, 1), 1), (SMatrix{2, 2, T}(1, 0, 0, 3), 2), (SMatrix{2, 2, T}(1, 0, 0, 0), 2),
                         (SMatrix{2, 2, T}(0, 0, 0, 1), 1), (SMatrix{2, 2, T}(big1, 0, 0, 1), 1), (SMatrix{2, 2, T}(1, 0, 0, -big1), 2))
            e = fl_health_error(() -> check_frame(r, W; orbit = 7))
            @test e isa NumericalHealthError && e.what == "initial frame" && e.orbit == 7 && e.step == 0 && e.backend === :cpu && e.column == col
            @test e isa NumericalHealthError && 0 < e.threshold < 1e-6 && e.value > e.threshold
        end
        for W in (SMatrix{2, 2, T}(1, 0, 0, 1), SMatrix{2, 2, T}(1, 0, 1, 0), SMatrix{2, 2, T}(0, 1, 0, -1), SMatrix{2, 2, T}(1 + T(epsT), 0, 0, 1))
            delta = nothing
            @test fl_health_error(() -> (delta = check_frame(r, W; orbit = 3))) === nothing
            @test delta !== nothing && 0 <= delta <= 2 * epsT
        end
        @test_throws ArgumentError check_frame(r, SMatrix{2, 3, T}(1, 0, 0, 1, 0, 0))
        nid = T === Float64 ? 10_000 : 1_000
        rejected = count(j -> fl_health_error(() -> check_frame(r, fl_frame(seed, j, T); orbit = j)) !== nothing, 1:nid)
        @test rejected == 0
    end
end

@testset "AC12c rejection before any work" begin
    _ = (FLI, check_frame)
    # Sali Amendment D67 form: the frame check runs before the stream, the tangent-core stream is unchanged. The caller checks every supplied frame first, so a rejected
    # frame never reaches the map (fl_Boom throws if forward or tangent is called); an accepted frame proceeds to the map.
    seed = fl_seed(12, 4)
    bad = SMatrix{2, 2, Float64}(2, 0, 0, 1)
    W0 = [fl_frame(seed, 1, Float64), bad]
    checked_run = () -> begin
        foreach(j -> check_frame(FLI(), W0[j]; orbit = j), eachindex(W0))
        propagate_orbit(fl_Boom(), SVector(0.0, 0.0), W0[2], 4, 2, (FLI(),))
    end
    err = fl_health_error(checked_run)
    @test err isa NumericalHealthError && err.what == "initial frame" && err.orbit == 2 && err.step == 0
    err = fl_health_error(() -> check_frame(FLI(), bad; orbit = 5))
    @test err isa NumericalHealthError && err.what == "initial frame" && err.orbit == 5 && err.step == 0
    @test fl_health_error(() -> check_frame(FLI(), W0[1]; orbit = 1)) === nothing
    @test_throws ErrorException propagate_orbit(fl_Boom(), SVector(0.0, 0.0), W0[1], 4, 2, (FLI(),))
end

@testset "AC12d overflow" begin
    _ = (FLI, SALI)
    # A = diag(mu, 1/mu), mu = e^3, Float64, chains {FLI (not the replacer), SALI} and {FLI}: tau = 200 raises no error at n = 200 (the columns arrive at e^600 and FLI reads
    # them before they are normalised) and FLI_B is 600 within the log of the first-row entries of the frame; tau = 300 raises the stream's error at handoff 300 (e^900
    # overflows), as sali 13(d).
    mu = exp(3.0)
    sys = LinearSymplecticMap(SMatrix{2, 2, Float64}(mu, 0, 0, 1 / mu))
    W0 = fl_frame(fl_seed(12, 5), 1, Float64)
    x0 = SVector(0.0, 0.0)
    for chain in ((FLI(replace = false), SALI()), (FLI(),))
        @test fl_health_error(() -> propagate_orbit(sys, x0, W0, 200, 200, chain)) === nothing
        _, _, st = propagate_orbit(sys, x0, W0, 200, 200, chain)
        f = reducer_final(FLI(), st[1])
        @test f.n == 200 && 599 < f.fli < 601
        err = fl_health_error(() -> propagate_orbit(sys, x0, W0, 300, 300, chain))
        @test err isa NumericalHealthError && err.what == "deviation matrix" && err.step == 300 && err.backend === :cpu
    end
end

@testset "AC12e run record" begin
    _ = (FLI,)
    # The run record carries tau, precision and k (stream) and the indicator name (reducer_provenance). In a chain the later reducer wins a shared key (D34 A2), so the
    # indicator name is tested on {FLI} alone; the frame source and the chain position are recorded by Experiments, not by this module.
    @test reducer_provenance(FLI()).indicator == "fli"
    @test reducer_provenance(FLI(replace = false)).indicator == "fli"
    sys = StandardMap(2.0)
    res = propagate(sys, [SVector(1.0, 2.0)], 8; k = 2, tau = 2, seed = fl_seed(12, 6), reducers = (FLI(),))
    @test res.provenance.indicator == "fli"
    @test res.provenance.k == 2 && res.provenance.tau == 2 && res.provenance.precision == "Float64"
    f = res.finals[1][1]
    @test f.n == 8 && f.sum isa SVector{2, Float64} && f.comp isa SVector{2, Float64}
    @test f.fli isa Float64 && f.fli1 isa Float64 && f.stepB isa Integer && f.colB isa Integer && f.step1 isa Integer
end

# ------------------------------------------------------------------ AC13

@testset "AC13a orbit and chain invariance" begin
    _ = (FLI, SALI, Benettin)
    # (a) the orbit at every handoff is bitwise equal for the chains {}, {L1} on k = 1, {FLI alone, the replacer} on k = 2, {FLI (not the replacer), SALI} on k = 2, {SALI}
    # on k = 2 and {} on k = 2 (L1 S2); (b) in {FLI, SALI}: the sali outputs (SALI, L, S_i, compensation) are bitwise equal to its stand-alone values, the FLI outputs are
    # bitwise equal to the stand-alone FLI outputs (the non-replacing position changes nothing), and the final deviation matrix is the same in all three runs;
    # K = 2 pi, 5 starts, N = 400, tau in {1, 4}, Float64.
    seed = fl_seed(13)
    T = Float64
    for tau in (1, 4)
        sys = StandardMap(fl_K(T, :two_pi))
        for (j, s) in enumerate(fl_starts(StableRNG(seed + 1), 5))
            x0 = fl_cast(T, s)
            W2 = fl_frame(seed, j, T)
            W1 = SMatrix{2, 1, T}(W2[:, 1])
            xa, _, sta = propagate_orbit(sys, x0, W1, 400, tau, (fl_Trace(),); orbit = j)
            xb, _, stb = propagate_orbit(sys, x0, W1, 400, tau, (Benettin(), fl_Trace()); orbit = j)
            xc, Wc, stc = propagate_orbit(sys, x0, W2, 400, tau, (FLI(), fl_Trace()); orbit = j)
            xd, Wd, std = propagate_orbit(sys, x0, W2, 400, tau, (FLI(replace = false), SALI(), fl_Trace()); orbit = j)
            xe, _, ste = propagate_orbit(sys, x0, W2, 400, tau, (fl_Trace(),); orbit = j)
            xf, Wf, stf = propagate_orbit(sys, x0, W2, 400, tau, (SALI(), fl_Trace()); orbit = j)
            tra, trb = reducer_final(fl_Trace(), sta[1]).xs, reducer_final(fl_Trace(), stb[2]).xs
            trc, trd = reducer_final(fl_Trace(), stc[2]).xs, reducer_final(fl_Trace(), std[3]).xs
            tre, trf = reducer_final(fl_Trace(), ste[1]).xs, reducer_final(fl_Trace(), stf[2]).xs
            @test length(tra) == 400 ÷ tau
            @test isequal(tra, trb) && isequal(tra, trc) && isequal(tra, trd) && isequal(tra, tre) && isequal(tra, trf)
            @test xa === xb && xa === xc && xa === xd && xa === xe && xa === xf
            # (b)
            fa, fd1 = reducer_final(FLI(), stc[1]), reducer_final(FLI(replace = false), std[1])
            @test isequal(fa, fd1)                                                  # FLI: stand-alone (replacer) equals the non-replacing position
            @test isequal(reducer_final(SALI(), std[2]), reducer_final(SALI(), stf[1]))   # SALI: unchanged by adding FLI before it
            @test isequal(Wc, Wf) && isequal(Wd, Wf)                                # the same unit columns in all chains with a replacer
            @test fa.n == 400 && fa.fli >= fa.fli1 >= 0
        end
    end
end

@testset "AC13b propagate equals the serial per-orbit body" begin
    _ = (FLI, SALI)
    # propagate equals the serial map of propagate_orbit bitwise (x, W, finals) for chunk 64 and 1000: 4096 orbits, K = 2 pi, n = 1e3, tau in {1, 4}, Float64, chains
    # {FLI (not the replacer), SALI} and {FLI}. The thread count is the run's (-t 2 under Pkg.test), never changed here; 1 versus N threads in child processes is a
    # stage 04 check (test plan, criterion 13(c)).
    seed = fl_seed(13, 10)
    T = Float64
    for (reducers, taus) in (((FLI(replace = false), SALI()), (1, 4)), ((FLI(),), (4,)))
        for tau in taus
            sys = StandardMap(fl_K(T, :two_pi))
            m, N = 4096, 1000
            x0 = [fl_cast(T, s) for s in fl_starts(StableRNG(seed + 2), m)]
            W0 = [fl_frame(seed, j, T) for j in 1:m]
            serial = [propagate_orbit(sys, x0[j], W0[j], N, tau, reducers; orbit = j) for j in 1:m]
            xs, Ws = [s[1] for s in serial], [s[2] for s in serial]
            fs = [map(reducer_final, reducers, s[3]) for s in serial]
            for chunk in (64, 1000)
                res = propagate(sys, x0, N; k = 2, tau = tau, seed = seed, reducers = reducers, chunk = chunk, W0 = W0)
                @test isequal(res.x, xs)
                @test isequal(res.W, Ws)
                @test isequal(res.finals, fs)
            end
        end
    end
end

@testset "AC13c allocation" begin
    _ = (FLI, SALI)
    # 0 bytes per orbit after warm-up, n = 1e4, both types, tau in {1, 4}, chains {FLI} and {FLI (not the replacer), SALI}.
    seed = fl_seed(13, 20)
    for T in FL_TYPES, tau in (1, 4), reducers in ((FLI(),), (FLI(replace = false), SALI()))
        sys = StandardMap(fl_K(T, :two_pi))
        @test fl_alloc_bytes(sys, fl_cast(T, SVector(1.0, 2.0)), fl_frame(seed, 1, T), 10_000, tau, reducers) == 0
    end
end

@testset "AC13d column exchange, sign change" begin
    _ = (FLI,)
    # I4: exchanging the two columns of a supplied frame leaves FLI_B (and the pair of sums, exchanged) bitwise unchanged and FLI_1 becomes the running maximum of the
    # other column's sum with the 0 (bitwise); a sign change of either column leaves both forms and the sums bitwise unchanged. K = 2 pi and the K = 0.5 regular class,
    # 5 starts, N = 400, tau in {1, 4}, every handoff, Float64 and Double64.
    seed = fl_seed(13, 30)
    for T in FL_TYPES, label in (:two_pi, :half), tau in (1, 4)
        c = fl_case(label, T, 5)
        W0 = c.W0
        swap(W) = SMatrix{2, 2, T}(W[1, 2], W[2, 2], W[1, 1], W[2, 1])
        neg1(W) = SMatrix{2, 2, T}(-W[1, 1], -W[2, 1], W[1, 2], W[2, 2])
        neg2(W) = SMatrix{2, 2, T}(W[1, 1], W[2, 1], -W[1, 2], -W[2, 2])
        base = fl_run((FLI(),), c.sys, c.x0, W0, 400, tau; all = true)
        sw = fl_run((FLI(),), c.sys, c.x0, swap.(W0), 400, tau; all = true)
        n1 = fl_run((FLI(),), c.sys, c.x0, neg1.(W0), 400, tau; all = true)
        n2 = fl_run((FLI(),), c.sys, c.x0, neg2.(W0), 400, tau; all = true)
        for j in 1:5
            a, b = base[j][1].snaps, sw[j][1].snaps
            @test length(a) == 400 ÷ tau
            @test isequal([s.fli for s in a], [s.fli for s in b])
            @test isequal([s.sum[1] + s.comp[1] for s in a], [s.sum[2] + s.comp[2] for s in b])
            @test isequal([s.fli1 for s in b], accumulate(max, [s.sum[2] + s.comp[2] for s in a]; init = zero(T)))
            for other in (n1, n2)
                o = other[j][1].snaps
                @test isequal([s.fli for s in a], [s.fli for s in o]) && isequal([s.fli1 for s in a], [s.fli1 for s in o])
                @test isequal([s.sum for s in a], [s.sum for s in o]) && isequal([s.comp for s in a], [s.comp for s in o])
            end
        end
    end
end

# ------------------------------------------------------------------ AC14

# rho_j (sigma_max(Y_j) / norm(Y_j w), L1 criterion 7) and M_state (max(1, |K|, largest |coordinate| on the lifted orbit), S1) for the standard map at K = Kf
# (Float64), start x0 and vector w0, after n steps, in BigFloat(256) (copied from lb_ref_full of the lyap-benettin helpers, which cannot be loaded here).
function fl_rho_M(Kf::Float64, x0, w0, n::Int)
    return setprecision(BigFloat, FL_BITS) do
        K = BigFloat(Kf)
        x, p = fl_big(x0[1]), fl_big(x0[2])
        a, b, c_, d = one(K), zero(K), zero(K), one(K)
        M = max(1.0, abs(Kf))
        for _ in 1:n
            M = max(M, Float64(abs(x)), Float64(abs(p)))
            c = K * cos(x)
            a, b, c_, d = (1 + c) * a + c_, (1 + c) * b + d, c * a + c_, c * b + d
            p = p + K * sin(x)
            x = x + p
        end
        w1, w2 = fl_big(w0[1]), fl_big(w0[2])
        y1, y2 = a * w1 + b * w2, c_ * w1 + d * w2
        F = a^2 + b^2 + c_^2 + d^2
        dt = a * d - b * c_
        smax = sqrt((F + sqrt(F^2 - 4 * dt^2)) / 2)
        M = max(M, Float64(abs(x)), Float64(abs(p)))
        (Float64(smax / sqrt(y1^2 + y2^2)), M)
    end
end

@testset "AC14 wrapped versus lifted state" begin
    _ = (FLI,)
    # (a) K = 0, dyadic starts: the wrapped (state length 4, D28) and the lifted series bitwise identical (all handoffs, both forms, the sums and the maximisers), tolerance 0.
    # (b) regular class K = 0.5 (S7), n in {10, 100}, Float64 (50 starts): |FLI_wrapped - FLI_lifted| <= 64 K M_state j^2 eps (rho_j^(1) + rho_j^(2)) at the larger of the
    # bounds evaluated at the two maximiser handoffs (j = handoff of the wrapped run, then of the lifted run; the L1 criterion 11(b) form with the two-point rule of Q11,
    # borrowed, not derived for the running maximum); starts with max rho > 50 at a maximiser are excluded (at most 10 percent).
    seed = fl_seed(14)
    rng = StableRNG(seed)
    dy() = rand(rng, 0:2048) / 1024
    T = Float64
    for tau in (1, 3)
        zs = [SVector{2, T}(T(dy()), T(dy())) for _ in 1:50]
        W0 = [fl_frame(seed, j, T) for j in 1:50]
        rl = fl_run((FLI(),), StandardMap(exact(T, 0)), zs, W0, 300, tau; all = true)
        rw = fl_run((FLI(),), WrappedStandardMap(exact(T, 0)), torus_state.(zs), W0, 300, tau; all = true)
        @test all(j -> isequal(rl[j][1].snaps, rw[j][1].snaps), 1:50)
    end
    m = 50
    c = fl_case(:half, T, m)
    for n in (10, 100)
        rl = fl_run((FLI(),), c.sys, c.x0, c.W0, n, 1; steps = [n])
        rw = fl_run((FLI(),), WrappedStandardMap(0.5), torus_state.(c.x0), c.W0, n, 1; steps = [n])
        excluded = 0
        for j in 1:m
            a, b = rl[j][1].snaps[1], rw[j][1].snaps[1]
            rho(i, step) = fl_rho_M(0.5, c.x0[j], c.W0[j][:, i], step)
            bound(step) = step == 0 ? (0.0, 0.0) : (r1 = rho(1, step); r2 = rho(2, step); (64 * 0.5 * r1[2] * step^2 * eps(Float64) * (r1[1] + r2[1]), max(r1[1], r2[1])))
            pairs = ((a.stepB, b.stepB), (a.step1, b.step1))
            forms = ((a.fli, b.fli), (a.fli1, b.fli1))
            for (k, (sa, sb)) in enumerate(pairs)
                ba, bb = bound(sa), bound(sb)
                if max(ba[2], bb[2]) > 50
                    excluded += 1
                    continue
                end
                @test abs(forms[k][1] - forms[k][2]) <= max(ba[1], bb[1])
            end
        end
        @test excluded <= m ÷ 10 * 2
    end
end
