# AC11 to AC13: health and frames, chain / execution / memory invariance, precision ladder.
# AC14 (ChaosTools independent tangent line), AC15 (ChaosTools local growth line) and AC16 (CPU vs GPU, D68 form) are stage 04; the test plan maps them.
# Every testset starts with a line naming the new symbols it uses (UndefVarError when missing).

# Fault-map states (n, y, flag): flagged orbits misbehave.
sn_fault_x0s(m, flagged) = [SVector{3, Float64}(0, 0, j in flagged ? 1 : 0) for j in 1:m]
sn_identity_frames(::Type{T}, m) where {T} = [SMatrix{2, 2, T}(1, 0, 0, 1) for _ in 1:m]

# One orbit of a chain through the per-orbit body with the caller's workspace: (x, W, finals), finals one entry per reducer (the three-argument reducer_final).
function sn_orbit(chain::Tuple, sys, x0, W0, N::Int, tau::Int; orbit::Int = 1)
    T = eltype(x0)
    ws = zeros(UInt32, sn_wslen(chain, T, size(W0, 1)))
    x, W, st = propagate_orbit(sys, x0, W0, N, tau, chain; orbit = orbit, workspace = ws)
    return x, W, map((r, s) -> Chaoslab.reducer_final(r, s, ws), chain, st)
end

@testset "AC11a zero tangent: column and orbit id of the health error" begin
    _ = (StretchingNumbers,)
    # A zero norm is not caught by the finiteness check of the stream: the reducer records it. The zero-tangent fault map on flagged orbits, k = 2: an error at the first handoff, kind "deviation
    # matrix", the orbit id, the step, backend :cpu, column 1 (both columns zero, or column 1), value 0, nothing stored; a tangent that zeroes the second column only (supplied identity frame):
    # column 2; two flagged orbits: the lowest id first. Float64 through propagate (kernel path with the workspace, chunk 2), the chains {SN} and {SN (reader), MEGNO}; Double64 through the
    # per-orbit body. k = 1 is an ArgumentError before any work.
    seed = sn_seed(11, 1)
    Ws = sn_identity_frames(Float64, 6)
    rng_ = (-1.0, 1.0)
    for chain in ((StretchingNumbers(; range = rng_, bins = 100, interval = 102),), (StretchingNumbers(; range = rng_, bins = 100, interval = 102, replace = false), MEGNO()))
        ok = propagate(sn_Fault(:zero), sn_fault_x0s(6, ()), 12; k = 2, tau = 3, seed = seed, reducers = chain, chunk = 2, W0 = Ws)
        @test length(ok.finals) == 6
        for (kind, flagged, orbit, col) in ((:zero, (4,), 4, 1), (:zero2, (5,), 5, 2), (:zero, (3, 6), 3, 1), (:zero2, (6, 2), 2, 2))
            res = nothing
            err = sn_health_error(() -> (res = propagate(sn_Fault(kind), sn_fault_x0s(6, flagged), 12; k = 2, tau = 3, seed = seed, reducers = chain, chunk = 2, W0 = Ws)))
            @test err isa NumericalHealthError
            @test res === nothing
            @test err.what == "deviation matrix" && err.orbit == orbit && err.step == 3 && err.backend === :cpu && err.column == col && err.value == 0.0 && err.threshold == 0.0
        end
    end
    r = StretchingNumbers(; range = rng_, bins = 100, interval = 102)
    x = SVector{3, Double64}(0, 0, 1)
    err = sn_health_error(() -> propagate_orbit(sn_Fault(:zero), x, SMatrix{2, 2, Double64}(1, 0, 0, 1), 12, 3, (r,); orbit = 4))
    @test err isa NumericalHealthError && err.what == "deviation matrix" && err.orbit == 4 && err.step == 3 && err.backend === :cpu && err.column == 1 && err.value == 0.0
    err = sn_health_error(() -> propagate_orbit(sn_Fault(:zero2), x, SMatrix{2, 2, Double64}(1, 0, 0, 1), 12, 3, (r,); orbit = 4))
    @test err isa NumericalHealthError && err.column == 2
    e = sn_error(() -> propagate(sn_Fault(:zero), sn_fault_x0s(2, ()), 6; k = 1, tau = 3, seed = seed, reducers = (r,)))
    @test e isa ArgumentError && occursin("StretchingNumbers needs k = 2, got k = 1", e.msg)
end

@testset "AC11c frame check by the caller" begin
    _ = (StretchingNumbers, check_frame)
    # sali Amendment D67 form, fli Amendment D72 item 7. A supplied frame with a non-unit column (defect above 8 eps) or two parallel columns (the smaller of |w1 - w2| and |w1 + w2| at most
    # 8 eps) is rejected before any work with kind "initial frame", step 0, the orbit id, backend :cpu, the column (the non-unit one; column 2 for parallel columns), threshold 8 eps;
    # accepted: the seeded frames of ids 1..10^4 (Double64: 1..1000), unit pairs at angles 0.1 and 1e-3, a nearly parallel pair 40 eps apart. The return value is the larger norm defect.
    # A frame without two columns is an ArgumentError.
    r = StretchingNumbers(; range = (-1.0, 1.0), bins = 100, interval = 100)
    seed = sn_seed(11, 3)
    for T in SN_TYPES
        epsT = sn_eps(T)
        thr = 8 * epsT
        nj = T === Float64 ? 10_000 : 1000
        bad = 0
        for j in 1:nj
            d = check_frame(r, sn_frame(seed, j, T); orbit = j)
            bad += !(0 <= d <= 8 * effective_eps(T))
        end
        @test bad == 0
        w = SVector{2, T}(exact(T, 3 // 5), exact(T, 4 // 5))               # unit to roundoff in T (the Float64 literals are not unit in Double64)
        perp = SVector{2, T}(-w[2], w[1])
        for ang in (0.1, 1e-3)
            W = hcat(w, cos(T(ang)) * w + sin(T(ang)) * perp)
            @test check_frame(r, SMatrix{2, 2, T}(W); orbit = 3) <= 8 * effective_eps(T)
        end
        @test check_frame(r, SMatrix{2, 2, T}(hcat(w, w + 40 * effective_eps(T) * perp)); orbit = 3) <= 8 * effective_eps(T)
        two = T(2)
        # (frame, expected column, norm defect (true) or parallel (false))
        badframes = ((SMatrix{2, 2, T}(two * w[1], two * w[2], perp[1], perp[2]), 1, true), (SMatrix{2, 2, T}(w[1], w[2], 3 * perp[1], 3 * perp[2]), 2, true),
                     (SMatrix{2, 2, T}(w[1], w[2], perp[1] / 2, perp[2] / 2), 2, true), (SMatrix{2, 2, T}((1 + 1000 * effective_eps(T)) * w[1], (1 + 1000 * effective_eps(T)) * w[2], perp[1], perp[2]), 1, true),
                     (SMatrix{2, 2, T}(0, 0, perp[1], perp[2]), 1, true), (SMatrix{2, 2, T}(w[1], w[2], 0, 0), 2, true),
                     (SMatrix{2, 2, T}(w[1], w[2], w[1], w[2]), 2, false), (SMatrix{2, 2, T}(w[1], w[2], -w[1], -w[2]), 2, false),
                     (SMatrix{2, 2, T}(hcat(w, w + 4 * effective_eps(T) * perp)), 2, false), (SMatrix{2, 2, T}(hcat(w, -w + 4 * effective_eps(T) * perp)), 2, false))
        for (W, col, isnorm) in badframes
            e = sn_health_error(() -> check_frame(r, W; orbit = 7))
            @test e isa NumericalHealthError && e.what == "initial frame" && e.step == 0 && e.orbit == 7 && e.backend === :cpu && e.column == col && e.threshold == Float64(thr)
            @test e isa NumericalHealthError && (isnorm ? e.value > e.threshold : e.value <= e.threshold)
        end
        nanframe = SMatrix{2, 2, T}(w[1], w[2], T(NaN), perp[2])
        e = sn_health_error(() -> check_frame(r, nanframe; orbit = 7))
        @test e isa NumericalHealthError && e.what == "initial frame" && e.column == 2
        @test sn_error(() -> check_frame(r, SMatrix{2, 3, T}(1, 0, 0, 1, 0, 0))) isa ArgumentError
        @test sn_error(() -> check_frame(r, SMatrix{2, 1, T}(1, 0))) isa ArgumentError
    end
end

@testset "AC11d overflow of the deviation matrix" begin
    _ = (StretchingNumbers,)
    # As sali 13(d). Test map diag(mu, 1/mu), mu = e^3, supplied unit columns: tau = 200 raises no error at n = 200 (one handoff; z = +3 and -3), tau = 300 raises the stream's health error
    # at handoff 300 (the first column overflows, kind "deviation matrix", step 300), nothing stored. Range [-4, 4], B = 1001, the interval equal to tau.
    mu = exp(3.0)
    sys = sn_Diag(mu, 1 / mu)
    r(tau) = StretchingNumbers(; range = (-4.0, 4.0), bins = 1001, interval = tau)
    x, W, fin = sn_orbit((r(200),), sys, SVector(0.0, 0.0), SMatrix{2, 2, Float64}(1, 0, 0, 1), 200, 200)
    f = fin[1]
    @test f.handoffs == 1 && f.n == 200 && abs(f.mean[1] - 3) <= 1e-12 && abs(f.mean[2] + 3) <= 1e-12 && sum(sn_cnt(f, 1)) == 1 && sum(sn_cnt(f, 2)) == 1
    err = sn_health_error(() -> sn_orbit((r(300),), sys, SVector(0.0, 0.0), SMatrix{2, 2, Float64}(1, 0, 0, 1), 300, 300; orbit = 5))
    @test err isa NumericalHealthError && err.what == "deviation matrix" && err.step == 300 && err.orbit == 5 && err.backend === :cpu
    res = nothing
    err = sn_health_error(() -> (res = propagate(sys, [SVector(0.0, 0.0)], 300; k = 2, tau = 300, seed = 1, reducers = (r(300),), W0 = [SMatrix{2, 2, Float64}(1, 0, 0, 1)])))
    @test err isa NumericalHealthError && err.step == 300 && res === nothing
end

@testset "AC11e run record" begin
    _ = (StretchingNumbers, stretching_range)
    # The run record of propagate carries, from the stream, tau, k, precision, seed and the backend, and, from the reducer (reducer_provenance, merged over the chain): the indicator name, the
    # convention tag, B, s_lo, s_hi, h, the interval, the minimum points per cell and the replace role. The state mode, the frame source and the chain are recorded by Experiments (as
    # sali D67 item 2, fli D72, rli D76 item 6). The final fields of one orbit have the documented types.
    sys = WrappedStandardMap(2pi * 1.0)
    lo, hi = stretching_range(sys)
    B = 200
    r = StretchingNumbers(; range = (lo, hi), bins = B, interval = 200, min_per_cell = 50, replace = false)
    x0 = [torus_state(SVector(0.5 + j, 1.0)) for j in 1:3]
    res = propagate(sys, x0, 400; k = 2, tau = 2, seed = 11, reducers = (r, MEGNO()))
    pv = res.provenance
    @test pv.indicator == "megno"                                    # later reducers win on a shared key (D34 A2); the stretching keys below are unshared
    pr = reducer_provenance(r)
    @test pr.indicator == "stretching" && pr.convention == SN_CONVENTION
    @test pr.spectrum_bins == B && pr.spectrum_interval == 200 && pr.spectrum_min_per_cell == 50 && pr.replace === false
    @test pr.spectrum_lo == Float64(lo) && pr.spectrum_hi == Float64(hi) && pr.spectrum_cell_width == (Float64(hi) - Float64(lo)) / B
    @test all(k -> getproperty(pv, k) == getproperty(pr, k), (:spectrum_bins, :spectrum_lo, :spectrum_hi, :spectrum_cell_width, :spectrum_interval, :spectrum_min_per_cell, :replace))
    @test pv.tau == 2 && pv.k == 2 && pv.precision == "Float64" && pv.seed == 11
    res1 = propagate(sys, x0, 400; k = 2, tau = 2, seed = 11, reducers = (StretchingNumbers(; range = (lo, hi), bins = B, interval = 200),))
    @test res1.provenance.indicator == "stretching" && res1.provenance.convention == SN_CONVENTION && res1.provenance.replace === true
    f = res1.finals[1][1]
    @test f.n == 400 && f.handoffs == 200
    @test f.D isa Float64 && f.undersampled isa Bool && size(f.counts) == (B + 2, 2) && eltype(f.counts) <: Integer
    @test f.sum isa SVector{2, Float64} && f.comp isa SVector{2, Float64} && f.sum2 isa SVector{2, Float64} && f.comp2 isa SVector{2, Float64} && f.mean isa SVector{2, Float64} && f.std isa SVector{2, Float64}
    @test f.below isa SVector{2, <:Integer} && f.above isa SVector{2, <:Integer}
end

@testset "AC12a orbit and chain invariance" begin
    _ = (StretchingNumbers, stretching_range)
    # I8. K = 2 pi, 5 starts, tau in {1, 4}, both types: (a) the orbit points at every handoff are bitwise equal for the chains {empty}, {SN alone, replacer}, {FLI (reader), SN}, {SN (reader),
    # MEGNO}, {MEGNO} (the chains without a normalising reducer overflow after about 600 steps on a chaotic orbit, so they run 300 steps and are compared with the first 300 steps of the
    # normalising chains, which run to 1000); the final points are equal; (b) the other reducers' outputs are bitwise equal to their stand-alone values (FLI alone, MEGNO alone), and the
    # stretching output is bitwise equal to its stand-alone value in every chain.
    for T in SN_TYPES, tau in (1, 4)
        c = sn_case(:two_pi, T, 5)
        N0 = 300
        for j in 1:5
            sys, x0, W0 = c.sys, c.x0[j], c.W0[j]
            mk(; replace = true) = sn_red(sys, 1000, tau; replace = replace)
            _, _, f_empty = sn_orbit((sn_Trace(),), sys, x0, W0, N0, tau; orbit = j)
            chains = ((mk(), sn_Trace()), (FLI(replace = false), mk(), sn_Trace()), (mk(replace = false), MEGNO(), sn_Trace()), (MEGNO(), sn_Trace()))
            outs = [sn_orbit(ch, sys, x0, W0, 1000, tau; orbit = j) for ch in chains]
            trace(o) = o[3][end].xs
            for o in outs
                @test isequal(trace(o)[1:(N0 ÷ tau)], f_empty[1].xs[1:(N0 ÷ tau)]) && length(trace(o)) == 1000 ÷ tau
                @test o[1] === outs[1][1]
            end
            solo = sn_orbit((mk(),), sys, x0, W0, 1000, tau; orbit = j)
            fli_solo = sn_orbit((FLI(),), sys, x0, W0, 1000, tau; orbit = j)
            megno_solo = sn_orbit((MEGNO(),), sys, x0, W0, 1000, tau; orbit = j)
            @test isequal(outs[1][3][1], solo[3][1]) && isequal(outs[2][3][2], solo[3][1]) && isequal(outs[3][3][1], solo[3][1])
            @test isequal(outs[2][3][1], fli_solo[3][1]) && isequal(outs[3][3][2], megno_solo[3][1]) && isequal(outs[4][3][1], megno_solo[3][1])
            @test isequal(outs[1][2], solo[2])                         # the final W of the replacing chain is the stand-alone one
        end
    end
end

@testset "AC12b propagate equals the serial per-orbit body" begin
    _ = (StretchingNumbers, stretching_range)
    # I8. 4096 orbits (Float64; Double64 256), K = 2 pi, n = 1e3, tau in {1, 4}, B = 1000: the points, the deviation matrices, the counts, the sums and D of propagate (kernel path with the workspace,
    # the run's threads) are bitwise equal to the serial per-orbit body for chunk sizes 64 and 1000. The explicit 1, 2, 4 thread counts in child processes are stage 04 (V).
    seed = sn_seed(12, 2)
    for (T, m) in ((Float64, 4096), (Double64, 256)), tau in (1, 4)
        sys = WrappedStandardMap(exact(T, pi) * 2)
        rng = StableRNG(seed + 1)
        x0 = [torus_state(sn_cast(T, s)) for s in sn_starts(rng, m)]
        r = sn_red(sys, 1000, tau)
        serial = [sn_orbit((r,), sys, x0[j], initial_deviations(seed, j, Val(2), Val(2), T), 1000, tau; orbit = j) for j in 1:m]
        for chunk in (64, 1000)
            res = propagate(sys, x0, 1000; k = 2, tau = tau, seed = seed, reducers = (r,), chunk = chunk)
            @test sn_nbad(j -> res.x[j] === serial[j][1] && res.W[j] === serial[j][2] && isequal(res.finals[j][1], serial[j][3][1]), 1:m) == 0
        end
    end
end

@testset "AC12c no allocation after warm-up" begin
    _ = (StretchingNumbers, stretching_range)
    # I8. 0 bytes allocated per orbit after warm-up (the workspace is the caller's), n = 1e4, both types, tau in {1, 4}, B = 1000, the chains {SN} and {SN (reader), MEGNO}, wrapped state.
    for T in SN_TYPES, tau in (1, 4)
        c = sn_case(:two_pi, T, 1)
        for chain in ((sn_red(c.sys, 1000, tau),), (sn_red(c.sys, 1000, tau; replace = false), MEGNO()))
            ws = zeros(UInt32, sn_wslen(chain, T, 2))
            @test sn_alloc_bytes(c.sys, c.x0[1], c.W0[1], 10_000, tau, chain, ws) == 0
        end
    end
end

@testset "AC12d position in the ensemble" begin
    _ = (StretchingNumbers, stretching_range)
    # I8. The series of orbit j is the same whatever its position in the ensemble: 256 orbits, K = 2 pi, n = 500, tau = 4, supplied frames that follow the ids, positions permuted (seeded), finals,
    # points and matrices bitwise equal after un-permuting; Float64 and Double64.
    seed = sn_seed(12, 4)
    for T in SN_TYPES
        c = sn_case(:two_pi, T, 256)
        r = sn_red(c.sys, 1000, 4)
        perm = sortperm(rand(StableRNG(seed), 256))
        a = propagate(c.sys, c.x0, 500; k = 2, tau = 4, seed = seed, reducers = (r,), chunk = 64, W0 = c.W0)
        b = propagate(c.sys, c.x0[perm], 500; k = 2, tau = 4, seed = seed, reducers = (r,), chunk = 100, W0 = c.W0[perm])
        @test sn_nbad(i -> a.x[perm[i]] === b.x[i] && a.W[perm[i]] === b.W[i] && isequal(a.finals[perm[i]], b.finals[i]), 1:256) == 0
    end
end

@testset "AC12e memory formula and launch capacity" begin
    _ = (StretchingNumbers, reducer_bytes, reducer_workspace)
    # Compute, experiments D53g. The reducer bytes per orbit are the state bytes plus 4 (B + 2) per column (two columns): reducer_bytes(r, T, Val(2), Val(2)) = sizeof(state) + 8 (B + 2) exactly for
    # B in {100, 1000}; the state is 80 bytes in Float64 and 144 in Double64 in the spec, plus the record of a failure (value, step, column) in this plan: 80 to 112 and 144 to 192 are accepted and the
    # actual size is reported. The stream bytes are the tangent-core formula (96 lifted and 112 wrapped in Float64). The chunk size (tangent-core) and the launch capacity (experiments) with
    # these bytes are floor(f free / (M + reducer bytes [+ series bytes])) exactly for three free-memory values.
    for T in SN_TYPES, B in (100, 1000)
        r = StretchingNumbers(; range = (-1.0, 1.0), bins = B, interval = 100)
        st = reducer_init(r, T, Val(2), Val(2))
        rb = reducer_bytes(r, T, Val(2), Val(2))
        @test reducer_workspace(r, T, Val(2), Val(2)) == 2 * (B + 2)
        @test rb == sizeof(st) + 8 * (B + 2)
        @test (T === Float64 ? 80 <= sizeof(st) <= 112 : 144 <= sizeof(st) <= 192) && isbitstype(typeof(st))
        @info "AC12e reducer scalar bytes" T B sizeof_state = sizeof(st) reducer_bytes = rb
        for (s, M) in ((2, sizeof(T) * (2 + 4 + 4) + 16), (4, sizeof(T) * (4 + 4 + 4) + 16))
            @test bytes_per_orbit(Val(2), Val(2), T, s) == M
            for free in (2^28, 3 * 10^9, 2^33 + 12345)
                @test chunk_size(free, Val(2), Val(2), T; f = 0.5, reducer_bytes = rb, s = s) == max(1, fld(free, 2 * (M + rb)))
                @test launch_capacity(free; f = 0.5, M = M, reducer = rb, columns = 8, grid_steps = 400) == fld(free, 2 * (M + rb + 8 * 8 * 400))
            end
        end
    end
    @test bytes_per_orbit(Val(2), Val(2), Float64, 4) == 112 && bytes_per_orbit(Val(2), Val(2), Float64, 2) == 96
    # the existing reducers request no workspace and their bytes are the state's
    for r in (MEGNO(), FLI(), SALI())
        @test reducer_workspace(r, Float64, Val(2), Val(2)) == 0 && reducer_bytes(r, Float64, Val(2), Val(2)) == sizeof(reducer_init(r, Float64, Val(2), Val(2)))
    end
end

@testset "AC12f step-chunked relaunch carries the workspace" begin
    _ = (StretchingNumbers, stretching_range)
    # tangent-core Amendment D44 form. A run split into 1, 4 and 16 launches (forced by the test hook Tangent._LAUNCH_STEPS, the CPU backend; N = 1600, tau = 4, launches of 1600, 400 and 100
    # steps) gives bitwise equal points, matrices, counts, sums and D, and the provenance reports the launches; K = 2 pi, 64 orbits, B = 1000, chunk 16, Float64 and Double64.
    seed = sn_seed(12, 6)
    for T in SN_TYPES
        c = sn_case(:two_pi, T, 64)
        r = sn_red(c.sys, 1000, 4)
        run_with(L) = begin
            old = Chaoslab.Tangent._LAUNCH_STEPS[]
            Chaoslab.Tangent._LAUNCH_STEPS[] = L
            try
                propagate(c.sys, c.x0, 1600; k = 2, tau = 4, seed = seed, reducers = (r,), chunk = 16, W0 = c.W0)
            finally
                Chaoslab.Tangent._LAUNCH_STEPS[] = old
            end
        end
        base = run_with(1600)
        @test base.provenance.launches_per_chunk == 1
        for (L, launches) in ((400, 4), (100, 16))
            res = run_with(L)
            @test res.provenance.launches_per_chunk == launches && res.provenance.launch_steps == L
            @test isequal(res.x, base.x) && isequal(res.W, base.W) && isequal(res.finals, base.finals)
        end
    end
end

@testset "AC13a precision ladder on the regular class" begin
    _ = (StretchingNumbers, stretching_range)
    # K = 0.5 regular class (sn_panel-sized: 4 starts, SN_FULL: 20), n = 1e4, tau = 1, B = 1000, Float64 against Double64 on the same initial points: each precision's counts lie in the envelope of its
    # own BigFloat reference (SN14), and for every cell |N_F64,c - N_D64,c| <= the larger of the two envelope widths N+_c - N-_c (the two orbits differ at rounding level; the envelope
    # of SN14 is the only width the spec defines, Deviation: it does not contain the orbit divergence term, which is negligible here).
    N, B = 10_000, 1000
    m = SN_FULL ? 20 : 4
    widths, cnts = Dict{Any, Any}(), Dict{Any, Any}()
    nbad_env = 0
    for T in SN_TYPES
        epsT = sn_eps(T)
        c = sn_case(:half, T, m)
        refs = sn_refs(:half, T, m, N)
        lo, hi = stretching_range(c.sys)
        res = sn_run((sn_red(c.sys, B, 1),), c.sys, c.x0, c.W0, N, 1; steps = [N])
        for j in 1:m, i in 1:2
            zs = sn_z_from_s(i == 1 ? refs[j].s1 : refs[j].s2, 1)
            eps = sn_halfwidths(refs[j], i, epsT, 1, N, lo, hi)
            e = sn_envelope(zs, eps, lo, hi, B, [N])[N]
            cnt = sn_cnt(res[j][1].snaps[1], i)
            nbad_env += !(all(e.Nm .<= cnt) && all(cnt .<= e.Np))
            widths[(T, j, i)] = e.Np .- e.Nm
            cnts[(T, j, i)] = cnt
        end
    end
    @test nbad_env == 0
    @test sn_nbad(ji -> all(abs.(cnts[(Float64, ji...)] .- cnts[(Double64, ji...)]) .<= max.(widths[(Float64, ji...)], widths[(Double64, ji...)])), [(j, i) for j in 1:m for i in 1:2]) == 0
end

@testset "AC13b precision ladder, ensemble means" begin
    _ = (StretchingNumbers, stretching_range, spectral_distance)
    # K in {2 pi, 5}, 1000 starts (the same initial points widened), n = 1e3, tau = 1, B = 1000, Float64 against Double64 through propagate: the ensemble means of the final mean z (both columns)
    # within 4 combined standard errors (L1 S5). Reported (R8): the distance D between the pooled spectra of the two precisions against the split-half distance of one precision.
    N, B = 1000, 1000
    for label in (:two_pi, :five)
        means, pooled = Dict{Any, Any}(), Dict{Any, Any}()
        for T in SN_TYPES
            K = sn_K(T, label)
            sys = WrappedStandardMap(K)
            rng = StableRNG(sn_seed(13, 2))
            x0 = [torus_state(sn_cast(T, s)) for s in sn_starts(rng, 1000)]
            r = sn_red(sys, B, 1)
            res = propagate(sys, x0, N; k = 2, tau = 1, seed = sn_seed(13, 3), reducers = (r,))
            fin = [f[1] for f in res.finals]
            means[T] = [Float64[f.mean[col] for f in fin] for col in 1:2]
            pooled[T] = [reduce(+, (sn_cnt(f, col) for f in fin)) for col in 1:2]
            pooled[(T, :half)] = [(reduce(+, (sn_cnt(f, col) for f in fin[1:500])), reduce(+, (sn_cnt(f, col) for f in fin[501:end]))) for col in 1:2]
        end
        for col in 1:2
            @test sn_agree(means[Float64][col], means[Double64][col])
        end
        h = 2 * Float64(stretching_range(WrappedStandardMap(sn_K(Float64, label)))[2]) / B
        dprec = Float64(spectral_distance(pooled[Float64][1], pooled[Double64][1], h))
        dsplit = Float64(spectral_distance(pooled[(Float64, :half)][1][1], pooled[(Float64, :half)][1][2], h))
        @info "AC13b (R8) pooled spectrum distance between precisions against the split-half distance of Float64" label dprec dsplit
    end
end

@testset "AC13c BigFloat against the 512-bit reference" begin
    _ = (StretchingNumbers, stretching_range)
    # K = 2 pi, 50 starts, n = 60, tau = 1, B = 100: BigFloat(256) on its own orbits against the 512-bit reference at its realized points: the counts of both columns lie in the envelope of SN14 with the
    # BigFloat epsilon (effective_eps at 256 bits), so equal to the reference counts unless a sample is within 1e-70 of an edge.
    seed = sn_seed(13, 4)
    setprecision(BigFloat, 256) do
        K = exact(BigFloat, pi) * 2
        sys = WrappedStandardMap(K)
        rng = StableRNG(seed)
        m, N, B = 50, 60, 100
        lo, hi = stretching_range(sys)
        epsB = Float64(effective_eps(BigFloat))
        r = StretchingNumbers(; range = (lo, hi), bins = B, interval = 60)
        nbad = 0
        for (j, s) in enumerate(sn_starts(rng, m))
            x0 = torus_state(SVector{2, BigFloat}(BigFloat(s[1]), BigFloat(s[2])))
            W0 = initial_deviations(seed, j, Val(2), Val(2), BigFloat)
            ref = sn_ref_orbit(sys, K, x0, W0, N; bits = 512)
            f = sn_run((r,), sys, [x0], [W0], N, 1; steps = [N])[1][1].snaps[1]
            for i in 1:2
                zs = sn_z_from_s(i == 1 ? ref.s1 : ref.s2, 1; bits = 512)
                eps = sn_halfwidths(ref, i, epsB, 1, N, lo, hi)
                e = sn_envelope(zs, eps, lo, hi, B, [N]; bits = 512)[N]
                cnt = sn_cnt(f, i)
                nbad += !(all(e.Nm .<= cnt) && all(cnt .<= e.Np))
            end
        end
        @test nbad == 0
    end
end

@testset "AC13d representation at K = 10" begin
    _ = (StretchingNumbers, stretching_range)
    # K = 10, n = 1e4, Float64 (5 starts, tau in {1, 4}): no overflow sample under the default range, D finite and inside [0, sqrt(2 / h)] (within 4 eps), the standard deviation of z positive and
    # finite, the counts sum to k: nothing is clipped, floored or capped.
    T = Float64
    epsT = sn_eps(T)
    for tau in (1, 4)
        c = sn_case(:ten, T, 5)
        lo, hi = stretching_range(c.sys)
        B = 1000
        h = (hi - lo) / B
        res = sn_run((sn_red(c.sys, B, tau),), c.sys, c.x0, c.W0, 10_000, tau; steps = [10_000])
        for j in 1:5
            f = res[j][1].snaps[1]
            @test f.below == [0, 0] && f.above == [0, 0] && sum(sn_cnt(f, 1)) == 10_000 ÷ tau && sum(sn_cnt(f, 2)) == 10_000 ÷ tau
            @test isfinite(f.D) && f.D > 0 && f.D <= sqrt(2 / h) * (1 + 4 * epsT)
            @test all(isfinite, f.std) && all(>(0), f.std) && all(isfinite, f.mean) && all(>(0), f.mean)
        end
    end
end
