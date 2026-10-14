# Amendment D34 (M2 Batch B API contract, A1b and A3), included after the diagnostics testset in runtests.jl.
# Existing testsets are untouched and keep passing (none compares a ResidualMonitor final NamedTuple as a whole).
#   A1b  det_threshold, symp_threshold and sum_threshold use Foundation.effective_eps (Double64: 2^-100, not eps = 2^-104), so a
#        Double64 chain no longer raises false alarms where the rounding model (precision spec Amendments) is met.
#   A3   reducer_final(::ResidualMonitor, st) gains M_max and G_max so the caller can form sum_threshold at the end of a run.
# The helpers dg_catch, dg_DET_TOL, dg_TYPES come from runtests.jl. Probe maps are LinearSymplecticMap (constant tangent), so no
# test-defined map type is needed.

# Renormalising reducer for k = 1 (this file's own copy: the one in faultmaps.jl is loaded lazily by other testsets).
struct dg34_Renorm end
Chaoslab.reducer_init(::dg34_Renorm, ::Type, ::Val, ::Val) = nothing
Chaoslab.reducer_update(::dg34_Renorm, st, n::Int, x::SVector, W::SMatrix) = (st, W / sqrt(sum(abs2, W)))
Chaoslab.reducer_final(::dg34_Renorm, st) = (;)

@testset "D34 thresholds use effective_eps" begin
    # det: tol * epsE * (1 + max|M|); symp: tol * epsE * N^2 * (1 + max|M|)^2; sum: tol * epsE * (G_max / tau + M_max).
    # Powers of two make the products exact up to the one rounding of (1 + max|M|), hence the 4 eps relative slack.
    for T in dg_TYPES
        epsE = T(Float64(effective_eps(T)))
        slack = 4 * Float64(eps(T))
        M2 = SMatrix{2, 2, T}(2, 1, 1, 1)
        M4 = SMatrix{4, 4, T}(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1)
        @test isapprox(Float64(det_threshold(4.0, T, M2)), Float64(T(4) * epsE * 3); rtol = slack)
        @test isapprox(Float64(symp_threshold(4.0, T, M2)), Float64(T(4) * epsE * 1 * 9); rtol = slack)
        @test isapprox(Float64(symp_threshold(4.0, T, M4)), Float64(T(4) * epsE * 4 * 4); rtol = slack)
        @test isapprox(Float64(sum_threshold(16.0, T, T(6.5), 2, T(2))), Float64(T(16) * epsE * (T(6.5) / 2 + 2)); rtol = slack)
    end
    # Double64 is exactly sixteen times the thresholds built with eps(Double64); Float64 is unchanged
    M2 = SMatrix{2, 2, Double64}(2, 1, 1, 1)
    @test isapprox(Float64(det_threshold(4.0, Double64, M2) / (4 * eps(Double64) * 3)), 16.0; rtol = 1e-12)
    @test isapprox(Float64(sum_threshold(16.0, Double64, 5, 1, 2) / (16 * eps(Double64) * 7)), 16.0; rtol = 1e-12)
    M2f = SMatrix{2, 2, Float64}(2, 1, 1, 1)
    @test det_threshold(4.0, Float64, M2f) == 4.0 * eps(Float64) * 3.0
    @test symp_threshold(4.0, Float64, M2f) == 4.0 * eps(Float64) * 9.0
    @test sum_threshold(16.0, Float64, 5.0, 1, 2.0) == 16.0 * eps(Float64) * 7.0
end

@testset "D34 Double64 no false alarm" begin
    # Tangent diag(1 + delta, 1) has det residual exactly delta. With the effective eps the quiet side (theta/4, theta/2] must not
    # raise over 100 handoffs and the raise side [2 theta, 4 theta) raises with value delta and the effective threshold; with the
    # nominal eps(Double64) the quiet delta (about 2^-98) is sixteen times the threshold and the monitor raised a false alarm.
    for T in dg_TYPES
        epsE = Float64(effective_eps(T))
        theta(delta) = 4 * epsE * (2 + delta)                   # det_threshold(4, T, M) with max|M| = 1 + delta
        d0 = 8 * epsE
        raise = [exp2(j) for j in -130:20 if 2 * theta(d0) <= exp2(j) < 4 * theta(d0)]
        quiet = [exp2(j) for j in -130:20 if theta(d0) / 4 < exp2(j) <= theta(d0) / 2]
        @test !isempty(raise) && !isempty(quiet)
        run(delta, N) = begin
            sys = LinearSymplecticMap(SMatrix{2, 2, T}(1 + T(delta), 0, 0, 1))
            mon = ResidualMonitor(sys, dg_DET_TOL, Inf)
            dg_catch(() -> propagate_orbit(sys, SVector{2, T}(0, 0), SMatrix{2, 1, T}(1, 0), N, 1, (dg34_Renorm(), mon)))
        end
        for delta in quiet
            @test run(delta, 100) === nothing
        end
        for delta in raise
            e = run(delta, 1)
            @test e isa NumericalHealthError && e.what == "det" && e.value == delta && e.step == 1 && e.column == 0
            @test isapprox(e.threshold, theta(delta); rtol = 1e-9)   # the effective threshold, not the nominal one
        end
    end
    # a realistic Double64 orbit stays quiet (standard map K = 5, 2000 steps, seeded start)
    rng = StableRNG(dg_SEED + 2900)
    sys = StandardMap(Double64(5))
    mon = ResidualMonitor(sys, dg_DET_TOL, Inf)
    x0 = SVector{2, Double64}(2pi * rand(rng), 2pi * rand(rng))
    @test dg_catch(() -> propagate_orbit(sys, x0, SMatrix{2, 1, Double64}(1, 0), 2000, 1, (dg34_Renorm(), mon))) === nothing
end

@testset "D34 ResidualMonitor final M_max and G_max" begin
    # M_max is the largest max|DT| over the handoffs of the run (recomputed here by stepping the orbit), G_max is 0 in a plain run and
    # the running maximum of the G values passed to monitor_exponents afterwards; both are T-valued and sum_threshold accepts them.
    for T in dg_TYPES
        sys = StandardMap(T(2))
        mon = ResidualMonitor(sys, dg_DET_TOL, Inf)
        x0 = SVector{2, T}(1, 2)
        N = 100
        x, W, st = propagate_orbit(sys, x0, SMatrix{2, 1, T}(1, 0), N, 1, (dg34_Renorm(), mon))
        fin = reducer_final(mon, st[2])
        for name in (:det_max, :symp_max, :orth_max, :sum_max, :pair_max, :M_max, :G_max)
            @test hasproperty(fin, name)
        end
        z = x0
        mref = zero(T)
        for _ in 1:N
            z = forward(sys, z)                               # the monitor evaluates DT at the point after the window
            mref = max(mref, maximum(abs, tangent(sys, z)))
        end
        @test fin.M_max isa T && fin.G_max isa T
        @test fin.M_max == mref
        @test fin.G_max == 0
        @test fin.det_max >= 0 && fin.orth_max == 0           # existing fields keep their meaning (k = 1: no orthogonality)
        # monitor_exponents folds G and M_max in; the final shows the running maxima
        lam = SVector{2, T}(0.5, -0.5)
        st2 = monitor_exponents(st[2], lam, T(6.8), N; pairing_tol = 16.0, tau = 1, M_max = T(10), step = N, orbit = 1)
        fin2 = reducer_final(mon, st2)
        @test fin2.G_max == T(6.8) && fin2.M_max == T(10)
        st3 = monitor_exponents(st2, lam, T(1), N; pairing_tol = 16.0, tau = 1, M_max = T(1), step = N, orbit = 1)
        fin3 = reducer_final(mon, st3)
        @test fin3.G_max == T(6.8) && fin3.M_max == T(10)
        # the caller's end-of-run check from the final (Q1): sum_threshold(pairing_tol, T, G_max, tau, M_max)
        @test check_residual("sum", eps_sum(lam), sum_threshold(16.0, T, fin3.G_max, 1, fin3.M_max); step = N, orbit = 1) === nothing
    end
end
