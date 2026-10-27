# Self-test of the reference code of wb_helpers.jl against hand-computed values. It uses no name of the module under test (only verified modules: Systems, Tangent, Indicators' existing reducers), so it
# must pass before and after stage 03 (GREEN in the RED run).

@testset "SELF reference code against hand values" begin
    # weights from the definition: w_0 = 0, w(1/2) = exp(-4) (p = 1) and exp(-16) (p = 2), symmetric under jj <-> N - jj (jj = 1..N-1)
    for p in (1, 2)
        w = wb_wref(16, p)
        @test w[1] == 0 && abs(w[9] - exp(-BigFloat(4)^p)) < 1e-70
        @test all(j -> abs(w[j] - w[18 - j]) < 1e-70, 2:16) && all(>(0), w[2:end])
    end
    # the module-independent weight vector in T equals the BigFloat definition to a few eps (it is built from the existing verified birkhoff_weight)
    for T in WB_TYPES
        wv = wb_wvec(T, 16, 1)
        ref = wb_wref(16, 1)
        @test maximum(Float64(abs(wb_big(wv[j]) - ref[j]) / max(ref[j], 1e-300)) for j in 2:16) <= 8 * wb_eps(T)
        @test wv[1] == 0 && wv[2] == wv[16] && wv[5] == wv[13]
    end
    # window sums: a constant gives the constant in both windows; f_j = j gives N/2 and 3N/2 (the weighted mean of jj is N/2 by symmetry); S is the largest |f| of the 2N samples
    N = 16
    c = wb_dec("0.3")
    r = wb_ref_avgs(fill(c, 2N), N, 1)
    @test abs(r.avg1 - c) < 1e-70 && abs(r.avg2 - c) < 1e-70 && r.S == c && abs(r.B - sum(wb_wref(N, 1))) < 1e-70
    fj = BigFloat.(0:(2N - 1))
    r = wb_ref_avgs(fj, N, 1)
    @test abs(r.avg1 - N / 2) < 1e-70 && abs(r.avg2 - 3N / 2) < 1e-70 && r.S == 2N - 1
    r2 = wb_ref_avgs(fj, N, 2)
    @test abs(r2.avg1 - N / 2) < 1e-70 && abs(r2.avg2 - 3N / 2) < 1e-70
    rw = wb_ref_avgs(fj, N, 1; w = ones(N))                       # equal weights: the plain mean of 0..N-1 is (N - 1)/2
    @test abs(rw.avg1 - (N - 1) / 2) < 1e-70 && abs(rw.avg2 - (3N - 1) / 2) < 1e-70
    # the index of WB15: hand values (eps = 2^-52, S = 1, avg1 = 1, avg2 = 1 - 1/1000)
    e52 = 2.0^-52
    one_ = BigFloat(1)
    ix = wb_ref_index(one_, one_ - wb_dec("0.001"), one_, e52)
    @test abs(ix.delta - wb_dec("0.001")) < 1e-70 && abs(ix.absdig - 3) < 1e-70
    m = (2 - wb_dec("0.001")) / 2
    @test abs(ix.reldig - (3 + log10(m))) < 1e-70 && abs(ix.dig - ix.absdig) < 1e-70 && !ix.below_floor && !ix.zero
    ixf = wb_ref_index(one_, one_, one_, e52)                      # Delta = 0: absdig = -log10(32 eps S) = 52 log10 2 - log10 32 = 14.1479...
    @test abs(ixf.absdig - (52 * log10(BigFloat(2)) - log10(BigFloat(32)))) < 1e-70 && ixf.below_floor && !ixf.zero
    @test abs(ixf.reldig - ixf.absdig) < 1e-70                     # m = S = 1: reldig = absdig
    ixs = wb_ref_index(wb_dec("1e-16"), wb_dec("2e-16"), one_, e52) # m = 1.5e-16 <= phi = 7.1e-15: reldig undefined, absdig at the floor
    @test isnan(ixs.reldig) && abs(ixs.dig - ixs.absdig) < 1e-70 && ixs.below_floor
    ixz = wb_ref_index(zero(BigFloat), zero(BigFloat), zero(BigFloat), e52)
    @test ixz.zero && isnan(ixz.absdig) && isnan(ixz.reldig) && isnan(ixz.dig) && ixz.delta == 0
    # exact rotation signals: the recurrence equals the direct evaluation, cos 4 pi theta is the k = 2 signal
    alpha = sqrt(BigFloat(2)) - 1
    th0 = BigFloat(1) / 7
    cs, sn = wb_rot_big(th0, alpha, 1, 2000)
    c2, s2 = wb_rot_big(th0, alpha, 2, 2000)
    tp = 2 * wb_pi()
    for j in (0, 1, 777, 1999)
        @test abs(cs[j + 1] - cos(tp * (th0 + j * alpha))) < 1e-70 && abs(sn[j + 1] - sin(tp * (th0 + j * alpha))) < 1e-70 && abs(c2[j + 1] - cos(2 * tp * (th0 + j * alpha))) < 1e-70
    end
    @test maximum(abs.(cs .^ 2 .+ sn .^ 2 .- 1)) < 1e-70
    # the Poisson form WB17 (criterion 12(a)) against the direct weighted sums of the reference (N = 50, 100): they agree to 1e-30 (the spec's dry run, mpmath: 1.9e-21 at N = 50, 2.1e-30 at N = 100)
    for Nn in (50, 100)
        ser, _ = wb_rot_big(th0, alpha, 1, 2Nn)
        direct = wb_ref_avgs(ser, Nn, 1).avg1
        pois = wb_poisson_avg(Nn, th0, alpha)
        @test abs(direct - pois) < 1e-25
        @test abs(direct) > 1e-30 && abs(wb_poisson_avg(Nn, th0, alpha; nq = 1400) - pois) < 1e-30      # the quadrature is converged in the number of nodes
    end
    # the BigFloat orbit of the standard map against the wrapped Float64 orbit (K = 0.5, 200 steps): the lifted points agree to 1e-9 (the rounding perturbation of the Float64 orbit grows)
    z0 = SVector(1.1, 0.7)
    big_ = wb_orbit_big(0.5, z0, 200)
    sys = WrappedStandardMap(0.5)
    sts = wb_states(sys, torus_state(z0), 200)
    lifts = [lift(s) for s in sts]
    @test maximum(Float64(abs(big_[j][1] - lifts[j][1])) + Float64(abs(big_[j][2] - lifts[j][2])) for j in 1:200) < 1e-9
    # the five standard-map observables in BigFloat at a state: cos q, sin q, cos p, sin p, lifted p = p + 2 pi n (n = 0 here)
    s0 = SVector(0.3, 0.5, 0.0, 0.0)
    fr = wb_fref_std([s0])
    @test abs(fr[1, 1] - cos(BigFloat(0.3))) < 1e-70 && abs(fr[1, 4] - sin(BigFloat(0.5))) < 1e-70 && abs(fr[1, 5] - BigFloat(0.5)) < 1e-70
    s1 = SVector(0.3, 0.5, 2.0, 3.0)
    @test abs(wb_fref_std([s1])[1, 5] - (BigFloat(0.5) + 3 * wb_tp(Float64))) < 1e-70
    # criterion 7 pairing (derived, F3): the lifted increment of q between the samples j and j + 1 equals the kick form p + K sin q (lifted p) of the sample j; here the BigFloat orbit and K = 0.9
    big9 = wb_orbit_big(0.9, SVector(2.0, 3.0), 60)
    setprecision(BigFloat, WB_BITS) do
        for j in 1:59
            q, p = big9[j]
            @test abs((big9[j + 1][1] - q) - (p + BigFloat(0.9) * sin(q))) < 1e-70
        end
    end
    # the K = 0 map (WB16): p is constant exactly, q advances by p0 per step modulo 2 pi_T; the lifted p of the wrapped state is p0 exactly (n = 0)
    for T in WB_TYPES
        sys0 = WrappedStandardMap(wb_K(T, :zero))
        s0 = torus_state(SVector{2, T}(T(1.7), T(3.3)))
        sts = wb_states(sys0, s0, 400)
        @test all(s -> s[2] == s0[2] && s[4] == 0, sts) && length(wb_states_tau(sys0, s0, 10, 4)) == 10
    end
    # the callables of the helpers (any isbits callable of the orbit point is an observable): constant, affine, coordinate, kick, bad-at, tuple
    x = SVector(0.25, 0.5, 1.0, 2.0)
    @test wb_Const(0.3)(x) === 0.3 && wb_Aff(wb_Coord(1), 4.0, 2.0)(x) === 3.0 && wb_Coord(2)(x) === 0.5
    @test abs(wb_StdKick(0.5)(x) - (0.5 + 2 * pi * 2.0 + 0.5 * sin(0.25))) < 1e-14 && isnan(wb_BadAt(0.1)(x)) && wb_BadAt(0.3)(x) === 0.25 && wb_Tuple()(x) == (0.25, 0.25)
    @test isbitstype(typeof(wb_Aff(wb_Coord(1), 1.0, 0.0))) && isbitstype(typeof(wb_StdKick(0.5)))
    # the cases: the regular class and the panel starts are widened exactly; the stream helper with existing reducers gives the trace of the orbit
    c = wb_case(:two_pi, Float64, 3, 1)
    @test length(c.x0) == 3 && c.fseed == wb_seed(0, 1) && all(s -> 0 <= s[1] < 2pi && s[3] == 0, c.x0)
    x0 = c.x0[1]
    x, W, fin = wb_run((wb_Trace(), MEGNO()), c.sys, x0, wb_frame(wb_seed(0, 2), 1, Float64, 2, 2), 50, 1)
    @test length(fin[1].xs) == 50 && fin[1].xs[end] == x && fin[2].handoffs == 50 && wb_states(c.sys, x0, 50)[end] == x
    @test isnan(wb_quant([NaN, NaN], 0.5)) && wb_quant([1.0, 2.0, 3.0], 0.5) == 2.0
end
