# API tests: invariants I1 and the maximisers on real orbits, the reducer protocol (arithmetic of one handoff, failures), pure functions, prefix invariance.
# AC15 and AC16 are stage 04 (see the test plan).

@testset "API invariants and maximisers" begin
    _ = (FLI,)
    # I1 (Q4, Q6) on real orbits, K in {0.5 (regular class), 2 pi}, N = 400, tau in {1, 4}, every handoff, both types (Float64 10 starts, Double64 4): FLI_1 and FLI_B are
    # nondecreasing, FLI_B >= FLI_1 >= max(0, S_1) and FLI_B >= max(S_1, S_2). The recorded maximisers: the value of S_colB at the handoff stepB is FLI_B bitwise (0 and
    # (stepB, colB) = (0, 0) for the s = 0 term), the same for FLI_1 with column 1, stepB is a handoff not beyond n, and at every handoff where a form strictly increases its
    # maximiser is this handoff. (The convention for ties is not tested.)
    for T in FL_TYPES, label in (:half, :two_pi), tau in (1, 4)
        m = T === Float64 ? 10 : 4
        c = fl_case(label, T, m)
        run = fl_run((FLI(),), c.sys, c.x0, c.W0, 400, tau; all = true)
        for r in run
            sn = r[1].snaps
            @test length(sn) == 400 ÷ tau
            S1 = [s.sum[1] + s.comp[1] for s in sn]
            S2 = [s.sum[2] + s.comp[2] for s in sn]
            f1, fB = [s.fli1 for s in sn], [s.fli for s in sn]
            @test issorted(f1) && issorted(fB)
            @test all(i -> fB[i] >= f1[i] >= max(zero(T), S1[i]) && fB[i] >= max(S1[i], S2[i]), eachindex(sn))
            Sc = (S1, S2)
            for (i, s) in enumerate(sn)
                @test s.n == i * tau
                if s.stepB == 0
                    @test s.colB == 0 && s.fli == 0
                else
                    @test s.stepB % tau == 0 && 0 < s.stepB <= s.n && s.colB in (1, 2)
                    @test s.fli == Sc[s.colB][s.stepB ÷ tau]
                end
                if s.step1 == 0
                    @test s.fli1 == 0
                else
                    @test s.step1 % tau == 0 && 0 < s.step1 <= s.n
                    @test s.fli1 == S1[s.step1 ÷ tau]
                end
                prevB, prev1 = i == 1 ? zero(T) : fB[i - 1], i == 1 ? zero(T) : f1[i - 1]
                fB[i] > prevB && @test s.stepB == s.n
                f1[i] > prev1 && @test s.step1 == s.n
            end
        end
    end
end

@testset "API pure functions" begin
    _ = (fli_unlogged,)
    # fli_unlogged(F) = exp(F) in the type of F when it is finite, `nothing` when the exponential overflows (above 709.78 in Float64): the display conversion of the log
    # form (spec Representation); it never returns an infinity.
    for T in FL_TYPES
        epsT = fl_eps(T)
        @test fli_unlogged(T(0)) isa T && fli_unlogged(T(0)) == 1
        for F in (T(1), T(0.25), T(10), T(300), T(700))
            v = fli_unlogged(F)
            @test v isa T
            @test fl_err(v, setprecision(() -> exp(fl_big(F)), BigFloat, FL_BITS)) <= 8 * epsT * Float64(v)
        end
    end
    @test fli_unlogged(709.78) !== nothing && isfinite(fli_unlogged(709.78))
    @test fli_unlogged(709.79) === nothing && fli_unlogged(1.0e5) === nothing
end

@testset "API reducer protocol" begin
    _ = (FLI,)
    # init is total except k != 2 (ArgumentError "FLI needs k = 2, got k = $k"); the n = 0 final; the arithmetic of one handoff on the columns (3, 4) and (0, 5): alpha = 5, 5,
    # S_i = ln 5, FLI_1 = FLI_B = ln 5 attained at this handoff; the replacer returns the column-normalised matrix, the non-replacing position returns W unchanged; a later
    # smaller column does not lower the running maxima (memory); column 2 larger sets colB = 2; contracting columns leave the s = 0 term (0, step 0, column 0); a failure
    # (zero or non-finite norm) is recorded without throwing in reducer_step, sums and W unchanged, thrown by reducer_update (orbit 0, backend :cpu).
    @test_throws ArgumentError reducer_init(FLI(), Float64, Val(2), Val(1))
    @test_throws ArgumentError reducer_init(FLI(replace = false), Float64, Val(2), Val(3))
    e1 = try
        reducer_init(FLI(), Float64, Val(2), Val(1))
    catch ex
        ex
    end
    @test e1 isa ArgumentError && occursin("FLI needs k = 2, got k = 1", e1.msg)
    for T in FL_TYPES, replace in (true, false)
        epsT = fl_eps(T)
        lg(x) = setprecision(() -> log(BigFloat(x)), BigFloat, FL_BITS)
        r = FLI(replace = replace)
        st0 = reducer_init(r, T, Val(2), Val(2))
        f0 = reducer_final(r, st0)
        @test f0.n == 0 && all(iszero, f0.sum) && all(iszero, f0.comp) && f0.fli == 0 && f0.fli1 == 0
        @test f0.stepB == 0 && f0.colB == 0 && f0.step1 == 0
        @test f0.sum isa SVector{2, T} && f0.comp isa SVector{2, T} && f0.fli isa T && f0.fli1 isa T
        @test !reducer_failed(r, st0)
        x = SVector{2, T}(0, 0)
        W = SMatrix{2, 2, T}(3, 4, 0, 5)
        st1, W1 = reducer_step(r, st0, 7, x, W)
        f1 = reducer_final(r, st1)
        @test f1.n == 7
        @test fl_err(f1.sum[1] + f1.comp[1], lg(5)) <= 2 * epsT && fl_err(f1.sum[2] + f1.comp[2], lg(5)) <= 2 * epsT
        @test fl_err(f1.fli, lg(5)) <= 2 * epsT && fl_err(f1.fli1, lg(5)) <= 2 * epsT
        @test f1.stepB == 7 && f1.step1 == 7 && f1.colB in (1, 2)
        if replace
            @test all(abs.(W1 - SMatrix{2, 2, T}(T(3) / 5, T(4) / 5, 0, 1)) .<= 2 * epsT)
        else
            @test W1 == W
        end
        stu, Wu = reducer_update(r, st0, 7, x, W)
        @test isequal(stu, st1) && isequal(Wu, W1)
        st2, _ = reducer_step(r, st1, 14, x, SMatrix{2, 2, T}(6, 8, 0, 10))               # alpha = 10, 10
        f2 = reducer_final(r, st2)
        @test f2.n == 14 && fl_err(f2.fli, lg(50)) <= 4 * epsT && fl_err(f2.fli1, lg(50)) <= 4 * epsT
        @test f2.stepB == 14 && f2.step1 == 14
        st3, _ = reducer_step(r, st2, 21, x, SMatrix{2, 2, T}(0.375, 0.5, 0, 0.5))        # alpha = 5/8, 1/2: both sums fall, the maxima stay
        f3 = reducer_final(r, st3)
        @test f3.n == 21 && f3.fli == f2.fli && f3.fli1 == f2.fli1 && f3.stepB == 14 && f3.step1 == 14 && f3.colB == f2.colB
        @test fl_err(f3.sum[1] + f3.comp[1], lg(50 * 5 / 8)) <= 6 * epsT && fl_err(f3.sum[2] + f3.comp[2], lg(25)) <= 6 * epsT
        # column 2 larger
        stb, _ = reducer_step(r, st0, 7, x, SMatrix{2, 2, T}(3, 4, 0, 10))
        fb = reducer_final(r, stb)
        @test fl_err(fb.fli, lg(10)) <= 2 * epsT && fb.colB == 2 && fb.stepB == 7
        @test fl_err(fb.fli1, lg(5)) <= 2 * epsT && fb.step1 == 7
        # contracting columns: the s = 0 term
        stc, _ = reducer_step(r, st0, 7, x, SMatrix{2, 2, T}(0.375, 0.5, 0, 0.5))
        fc = reducer_final(r, stc)
        @test fc.fli == 0 && fc.fli1 == 0 && fc.stepB == 0 && fc.colB == 0 && fc.step1 == 0
        @test fl_err(fc.sum[1] + fc.comp[1], lg(5 / 8)) <= 2 * epsT && fl_err(fc.sum[2] + fc.comp[2], lg(1 / 2)) <= 2 * epsT
        # failure: column 2 zero, column 1 zero, a non-finite column
        Wz = SMatrix{2, 2, T}(3, 4, 0, 0)
        stf, Wf = reducer_step(r, st2, 28, x, Wz)
        @test reducer_failed(r, stf) && !reducer_failed(r, st2)
        @test Wf === Wz
        ff = reducer_final(r, stf)
        @test ff.sum == f2.sum && ff.comp == f2.comp && ff.n == f2.n && ff.fli == f2.fli && ff.fli1 == f2.fli1
        e = reducer_error(r, stf, 9, :cpu)
        @test e isa NumericalHealthError && e.what == "deviation matrix" && e.orbit == 9 && e.step == 28
        @test e.backend === :cpu && e.column == 2 && e.value == 0.0
        e2 = fl_health_error(() -> reducer_update(r, st2, 28, x, Wz))
        @test e2 isa NumericalHealthError && e2.what == "deviation matrix" && e2.orbit == 0 && e2.backend === :cpu && e2.step == 28 && e2.column == 2
        st1z, _ = reducer_step(r, st2, 28, x, SMatrix{2, 2, T}(0, 0, 0, 5))
        @test reducer_failed(r, st1z) && reducer_error(r, st1z, 1, :cpu).column == 1
        stn, _ = reducer_step(r, st2, 28, x, SMatrix{2, 2, T}(3, 4, T(NaN), 5))
        @test reducer_failed(r, stn)
        sti, _ = reducer_step(r, st2, 28, x, SMatrix{2, 2, T}(3, 4, T(Inf), 5))
        @test reducer_failed(r, sti)
    end
    @test reducer_provenance(FLI()).indicator == "fli"
end

@testset "API prefix invariance" begin
    _ = (FLI,)
    # The state after handoff j is bitwise identical for every N >= j tau (same seed, orbit, start, tau): a run of length j tau returns exactly the series value at j tau.
    # All-handoff snapshots of one run of 400 steps against separate runs, tau in {1, 4}, both types; the maximisers are part of the final and are compared too.
    seed = fl_seed(15, 20)
    for T in FL_TYPES, tau in (1, 4)
        sys = StandardMap(fl_K(T, :two_pi))
        x0 = fl_cast(T, SVector(1.0, 2.0))
        W0 = fl_frame(seed, 3, T)
        probe = fl_Probe(FLI(); all = true)
        _, _, st = propagate_orbit(sys, x0, W0, 400, tau, (probe,); orbit = 3)
        snaps = reducer_final(probe, st[1]).snaps
        @test length(snaps) == 400 ÷ tau
        for j in (1, 5, 25, 100 ÷ tau)
            N = j * tau
            _, _, s2 = propagate_orbit(sys, x0, W0, N, tau, (FLI(),); orbit = 3)
            f = reducer_final(FLI(), s2[1])
            @test snaps[j].n == N && f.n == N
            @test isequal(snaps[j], f)
        end
    end
end
