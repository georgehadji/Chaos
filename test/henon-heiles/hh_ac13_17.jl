# AC13 to AC17: full spectrum, known tangent at the fixed points, GALI_k laws, determinism across chunks and threads, no allocation.
# Every testset starts by naming the new symbols it uses, so a missing module errors that testset (UndefVarError).
# Exponents of the estimators are per step; hh_lyap divides by the step size, so every exponent below is per time unit unless named "step".
# Stage 04 only: the Taylor and ChaosTools spectra (O3, O4) of AC13(d) at the spec's reference, the GPU part of AC16 (see the test plan).

# the first seeded candidate at energy E (window as the section domain) whose finite-time exponent over 500 time units exceeds 0.02: a chaotic start
function hh_chaotic_start(flow, m, E, rng; tau = 0.05)
    a, b = section_interval(flow, E); pm = sqrt(2E)
    for _ in 1:200
        z, st = random_energy_start(flow, E, (a, b), (-pm, pm), rng)
        st === :ok && hh_lambda1(m, flow, z, tau, 500.0) > 0.02 && return z
    end
    return nothing
end

@testset "AC13a to c full spectrum of chaotic orbits: pairing, zero exponents, sum" begin
    _ = (energy_start, random_energy_start, section_interval)
    flow = HenonHeilesFlow(); m = KahanLi(Float64, 8); tau = 0.05; T = hh_size(20000, 100000)
    pt = hh_defaults()["pairing_residual_tol"]
    zC, _ = energy_start(flow, 1 / 8, -0.015, 0.0)
    z6 = hh_chaotic_start(flow, m, 1 / 6, StableRNG(hh_seed(13, 1)))
    @test z6 isa SVector{4, Float64}
    for (name, z0) in (("C", zC), ("E1/6", z6)), alg in (RGS(), MGS())
        r = hh_lyap(m, flow, z0, tau, T; alg = alg); lam = r.lambda
        thr_step = sum_threshold(pt, Float64, r.G_max, 20, r.M)                                  # diagnostics-spec Tolerances, per step; the renormalization interval is 20 steps
        allow = 10 * log(T) / T                                                                  # finite-time allowance of the pairing (diagnostics-spec O2: kappa has no default), per time
        @test lam[1] > 0.005                                                                     # chaotic
        @test eps_pair(lam) <= thr_step / tau + allow                                            # (a) |l1 + l4|, |l2 + l3|
        @test abs(lam[2]) <= 10 * log(T) / T && abs(lam[3]) <= 10 * log(T) / T                   # (b) the two exponents that tend to 0
        @test eps_sum(r.lstep) <= thr_step                                                       # (c) sum of the exponents
        @info "AC13 reported" name alg lambda = lam
    end
end

@testset "AC13d spectrum at T = 100 against the Double64 order-10 reference (package oracle; Taylor and ChaosTools are stage 04)" begin
    _ = (energy_start,)
    flow = HenonHeilesFlow(); m = KahanLi(Float64, 8); tau = 0.05
    zC, _ = energy_start(flow, 1 / 8, -0.015, 0.0)
    W0 = initial_deviations(hh_seed(13, 2), 1, Val(4), Val(4), Float64)
    r = hh_lyap(m, flow, zC, tau, 100.0; W0 = W0)
    D = Double64
    rd = hh_lyap(hh_refmethod(), flow, SVector{4, D}(D.(zC)), D(1) / 80, 100.0; W0 = SMatrix{4, 4, D}(D.(W0)), renorm = 80)
    @test maximum(abs.(r.lambda .- Float64.(rd.lambda))) <= 1.0e-6
end

@testset "AC13e integrator effect on the ensemble mean of lambda_1" begin
    _ = (random_energy_start, section_interval)
    flow = HenonHeilesFlow(); E = 1 / 8; a, b = section_interval(flow, E); pm = sqrt(2E)
    rng = StableRNG(hh_seed(13, 3)); nens = hh_size(30, 100); Te = hh_size(2000, 10000)
    starts = [random_energy_start(flow, E, (a, b), (-pm, pm), rng)[1] for _ in 1:nens]
    W0s = [initial_deviations(hh_seed(13, 4), i, Val(4), Val(1), Float64) for i in 1:nens]
    lam1(m, tau) = begin
        im = IntegratorMap(m, flow, tau); rn = round(Int, 1 / tau)
        res = propagate(im, starts, round(Int, Te / tau); k = 1, tau = rn, seed = 1, reducers = (Benettin(),), W0 = W0s)
        [f[1].lambda[1] / tau for f in res.finals]
    end
    l8 = lam1(KahanLi(Float64, 8), 1 / 20); l8h = lam1(KahanLi(Float64, 8), 1 / 40); l4 = lam1(ForestRuth(Float64), 1 / 20)
    se = hh_std(l8) / sqrt(nens)
    @test abs(hh_mean(l8) - hh_mean(l8h)) <= 3 * se                                              # tau = 1/20 against tau = 1/40
    @test abs(hh_mean(l8) - hh_mean(l4)) <= 3 * se                                               # order-8 against order-4 table
    @info "AC13e reported" mean8 = hh_mean(l8) mean8_half = hh_mean(l8h) mean4 = hh_mean(l4) se
end

@testset "AC14 known tangent at the fixed points S0 and O" begin
    _ = (critical_points,)
    flow = HenonHeilesFlow(); m = KahanLi(Float64, 8); tau = 0.05; im = IntegratorMap(m, flow, tau)
    cp = critical_points(flow, Float64); I4 = one(SMatrix{4, 4, Float64})
    zS = SVector(cp[2][1], cp[2][2], 0.0, 0.0); zO = SVector(cp[1][1], cp[1][2], 0.0, 0.0)
    # S0 = (0, 1): the tangent flow has constant matrix diag-block Hess = diag(3, -1) (Equation 16)
    t = 5.0; s3 = sqrt(3.0); c = cos(s3 * t); s = sin(s3 * t)
    Mex = zeros(4, 4)
    Mex[1, 1] = c; Mex[1, 3] = s / s3; Mex[3, 1] = -s3 * s; Mex[3, 3] = c                      # (q1, p1) block
    Mex[2, 2] = cosh(t); Mex[2, 4] = sinh(t); Mex[4, 2] = sinh(t); Mex[4, 4] = cosh(t)        # (q2, p2) block
    M = propagate(im, [zS], round(Int, t / tau); k = 4, tau = round(Int, t / tau), seed = 1, W0 = [I4]).W[1]
    @test maximum(abs, Matrix(M) .- Mex) <= 1.0e-9
    lam = sort(collect(hh_lyap(m, flow, zS, tau, 100.0).lambda); rev = true)                    # spectrum (1, 0, 0, -1) within 0.05 at T = 100
    @test maximum(abs.(lam .- (1.0, 0.0, 0.0, -1.0))) <= 0.05
    # O = (0, 0): rotation with unit frequencies (the harmonic oscillator of the catalog), error at most B(n) at T = 10, all exponents tiny at T = 1000
    n = round(Int, 10 / tau)
    Mo = propagate(im, [zO], n; k = 4, tau = n, seed = 1, W0 = [I4]).W[1]
    @test maximum(abs, Matrix(Mo) .- Matrix(exact_flow_matrix(HarmonicOscillator([1.0, 1.0]), 10.0))) <= hh_B(Float64, 1.0, n, tau)
    @test maximum(abs, hh_lyap(m, flow, zO, tau, 1000.0).lambda) <= 1.0e-12
end

@testset "AC15a GALI_k of the regular orbit R: slopes 0, -2, -4" begin
    _ = (energy_start,)
    flow = HenonHeilesFlow(); m = KahanLi(Float64, 8); tau = 0.05; im = IntegratorMap(m, flow, tau)
    zR, _ = energy_start(flow, 1 / 8, 0.1, 0.0)
    ts = [100.0, 141.0, 200.0, 283.0, 400.0, 566.0, 800.0, 1131.0, 1600.0, 2000.0]
    for sd in 1:2                                                                                # generic seeded frames: the law of B7 is for generic vectors (spec point U2)
        W0 = initial_deviations(hh_seed(15, sd), 1, Val(4), Val(4), Float64)
        G = [propagate(im, [zR], round(Int, T / tau); k = 4, tau = 20, seed = 1, reducers = (GALI(),), W0 = [W0]).finals[1][1].gali for T in ts]
        for (k, expected) in ((2, 0.0), (3, -2.0), (4, -4.0))
            vals = [g[k - 1] for g in G]; idx = findall(>=(10 * gali_floor(Float64, k)), vals)     # resolved values only (B9)
            @test length(idx) >= 6
            @test abs(hh_slope(log10.(ts[idx]), log10.(vals[idx])) - expected) <= 0.25
        end
    end
end

@testset "AC15b GALI_k of the chaotic orbit C: exponential rates E_k lambda_1" begin
    _ = (energy_start,)
    flow = HenonHeilesFlow(); m = KahanLi(Float64, 8); tau = 0.05; im = IntegratorMap(m, flow, tau)
    zC, _ = energy_start(flow, 1 / 8, -0.015, 0.0)
    ts = collect(10.0:10.0:1000.0)
    for sd in 1:2
        W0 = initial_deviations(hh_seed(15, 10 + sd), 1, Val(4), Val(4), Float64)
        G = [propagate(im, [zC], round(Int, T / tau); k = 4, tau = 20, seed = 1, reducers = (GALI(),), W0 = [W0]).finals[1][1].gali for T in ts]
        S(T) = (r = hh_lyap(m, flow, zC, tau, T; W0 = W0); r.lstep .* r.n)                     # sums of the log growth (per step times steps)
        for k in 2:4
            vals = [g[k - 1] for g in G]
            idx = findall(v -> 1000 * gali_floor(Float64, k) <= v <= 1.0e-2, vals)             # the stretch above the measurement floor (B9)
            @test length(idx) >= 8
            t1, t2 = ts[idx[1]], ts[idx[end]]
            rate = -hh_slope(ts[idx], log.(vals[idx]))
            lam = sort(collect((S(t2) .- S(t1)) ./ (t2 - t1)); rev = true)               # exponents of the same stretch (finite-time, per time unit)
            Ek = sum(lam[1] - lam[j] for j in 2:k)
            @test 0.75 <= rate / Ek <= 1.25                                                    # within 25 percent of E_k lambda_1 (stretch exponents, spec point U3)
        end
    end
    # reported: the spec's literal reference, lambda_1 of the long run of AC13, against the same rates
    l1 = maximum(hh_lyap(m, flow, zC, tau, hh_size(20000.0, 100000.0)).lambda)
    @info "AC15b reported: lambda_1 of the long run (E_2 = 1, E_3 = 2, E_4 = 4 times this is the spec's literal reference)" l1
end

@testset "AC16 portrait labels and states are bitwise independent of chunking and threads" begin
    _ = (energy_grid,)
    m = KahanLi(Float64, 8); flow = HenonHeilesFlow(); n = hh_size(10, 40); Tf = hh_size(500.0, 2000.0)
    pa = hh_portrait(m, 1 / 8, n, Tf); pb = hh_portrait(m, 1 / 8, n, Tf; chunk = 7)
    @test pa.n > 0 && pa.x == pb.x && pa.W == pb.W && pa.chaotic == pb.chaotic
    im = IntegratorMap(m, flow, 1 / 20); N = round(Int, Tf / (1 / 20))
    serial = [propagate_orbit(im, z, HH_W2, N, 20, (SALI(),)) for z in pa.starts]               # the serial map of propagate_orbit
    @test [s[1] for s in serial] == pa.x
    @test [reducer_final(SALI(), s[3][1]).sali < 1.0e-8 for s in serial] == pa.chaotic
    @info "AC16 reported: threads of this run (the 1 and N thread comparison is the serial map above; the GPU comparison is stage 04)" Threads.nthreads()
end

@testset "AC17 no allocation after warm-up: step, tangent step, crossing solve, escape test" begin
    _ = (random_energy_start, section_interval, locate_crossing, locate_escape, escape_sector)
    flow = HenonHeilesFlow()
    for T in (Float64, Double64)
        m = KahanLi(T, 8); tau = T(1) / 20; one4 = one(SMatrix{4, 4, T}); rng = StableRNG(hh_seed(17, 1))
        z0, _ = random_energy_start(flow, T(1) / 8, (T(-1) / 10, T(1) / 10), (T(-1) / 5, T(1) / 5), rng)
        zb = hh_bracket(m, flow, z0, tau)
        @test zb !== nothing
        @test hh_alloc(advance, m, flow, zb, tau) == 0
        @test hh_alloc(tangent_advance, m, flow, zb, one4, tau) == 0
        @test hh_alloc(locate_crossing, m, flow, zb, tau) == 0
        # a state just inside r = R_esc whose next step is outside, from an escaping start at E = 1/4
        R = T(10); ze = nothing
        for _ in 1:100
            z, _ = random_energy_start(flow, T(1) / 4, (-T(1), T(1)), (-T(1), T(1)), rng)
            for _ in 1:4000
                zn = advance(m, flow, z, tau)
                (all(isfinite, zn) && hh_rad(zn) <= 2R) || break
                if hh_rad(z) <= R < hh_rad(zn)
                    ze = z; break
                end
                z = zn
            end
            ze === nothing || break
        end
        @test ze !== nothing
        @test hh_alloc(locate_escape, m, flow, ze, tau, R) == 0
        @test hh_alloc(escape_sector, flow, SVector(ze[1], ze[2])) == 0
    end
end
