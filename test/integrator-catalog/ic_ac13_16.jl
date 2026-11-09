# AC13 to AC16: long runs (bounded energy, angular momentum, Brouwer's law) and allocation / force-evaluation counts.
# Package sizes follow D38 (3 minute ceiling): runs of 1e4 steps instead of 1e5 (AC13, AC14), Brouwer at n up to 1e5 and 100 starts (spec: 3e6 and 200).
# verify/integrator-catalog.jl sets IC_FULL_RUN = true for the full sizes.

# energy samples of a run: H of the stored state in 256 bits at `nchk` equally spaced step counts
function ic_energy_series(m, kind, flow, z0, tau, nsteps, nchk; eps = 1.0)
    ns = [round(Int, j * nsteps / nchk) for j in 1:nchk]
    zs = advance_to(m, flow, ic_to(typeof(tau), z0), tau, ns)
    return setprecision(BigFloat, IC_BITS) do
        H0 = ic_H(kind, ic_to(BigFloat, z0); eps = eps)
        [Float64(abs(ic_H(kind, ic_big.(collect(z)); eps = eps) - H0)) for z in zs]
    end
end
ic_lsslope(y) = (n = length(y); x = collect(1.0:n); xm = sum(x) / n; sum((x .- xm) .* (y .- sum(y) / n)) / sum((x .- xm) .^ 2))

# D112: in binary64 the tables of order 6 or higher (and the generalized order (8, 6, 4) rows) are not truncation dominated at these steps: the energy error of 1e4 steps is the round-off random
# walk of Brouwer's law (1e-14 to 6e-14, growing like the square root of the step count), which the bounded-oscillation assertions are not built for. They are run in double-double (truncation
# 1e-17 to 1e-12, round-off 1e-29) and the others in binary64, the same assertions with eps_T of the element type used.
ic_ac13_type(M) = (M.p >= 6 || startswith(M.id, "ABA")) ? Double64 : Float64

@testset "AC13 energy bounded, no secular drift" begin
    _ = (advance_to, force_evaluations)
    nsteps = ic_size(10000, 100000); nchk = 1000
    cases = (("F1", :f1_3, IC_START[:f1_3], 1.0), ("F2e0.1", :f2, ic_kepler_state(0.1, 0.7), 1.0), ("F2e0.5", :f2, IC_START[:f2], 1.0), ("F4", :f4, IC_START[:f4], 1.0))
    for M in IC_S
        M.id in ("S1DK", "S1KD", "S4", "S3") && continue                                # criterion 13: S2 and S5 to S9, S11 to S13
        T = ic_ac13_type(M)
        tau = T(1) / ((startswith(M.id, "S2") || M.id == "S5" || M.id == "S8" || startswith(M.id, "S12")) ? 20 : 50)   # criterion 13: S2, S5, S8, S12 at 1/20
        for (name, kind, z0, eps) in cases
            E = ic_energy_series(M.mk(T), kind, ic_flow(T, kind; eps = eps), z0, tau, nsteps, nchk; eps = eps)
            half = nchk ÷ 2; Mn = ic_mmax(z0)
            @test maximum(E[(half + 1):end]) <= 2 * maximum(E[1:half]) + 16 * ic_eps(T) * Mn
            @test abs(ic_lsslope(E[(half + 1):end]) * half) <= 0.1 * maximum(E)
        end
    end
    for M in IC_P                                                                       # near-integrable rows on F4 with eps = 1e-2
        T = ic_ac13_type(M)
        E = ic_energy_series(M.mk(T), :f4, ic_flow(T, :f4; eps = 1e-2), [1.0, 0.3], T(1) / 50, nsteps, nchk; eps = 1e-2)
        half = nchk ÷ 2
        @test maximum(E[(half + 1):end]) <= 2 * maximum(E[1:half]) + 16 * ic_eps(T) * 1.0
        @test abs(ic_lsslope(E[(half + 1):end]) * half) <= 0.1 * maximum(E)
    end
end

@testset "AC14 angular momentum of F2 is preserved by every explicit method" begin
    _ = (advance_n, angular_momentum, runge_lenz)
    nsteps = ic_size(10000, 100000); tau = 1 / 50
    for e in (0.1, 0.5), M in IC_ALLX
        z0 = ic_kepler_state(e, 0.7); m = M.mk(Float64); flow = ic_flow(Float64, :f2)
        z = advance_n(m, flow, ic_to(Float64, z0), tau, nsteps)
        dL = setprecision(BigFloat, IC_BITS) do
            Float64(abs(ic_L(ic_to(BigFloat, collect(z))) - ic_L(ic_to(BigFloat, z0))))
        end
        @test dL <= ic_B(Float64, ic_mmax(z0), nsteps, tau)
        e == 0.5 && M.id in ("S2DKD", "S5", "S12-8") && @info "AC14 reported: Runge-Lenz drift (linear growth, not scored)" M.id maximum(abs, collect(runge_lenz(flow, z)) .- collect(runge_lenz(flow, ic_to(Float64, z0))))
    end
end

@testset "AC15 Brouwer's law on Henon-Heiles (T2), order 8 on the Stormer-Verlet base, package size" begin
    _ = (advance_to, KahanLi, HenonHeilesFlow)
    rng = StableRNG(ic_seed(15, 1))
    flow = HenonHeilesFlow(); m = KahanLi(Float64, 8); tau = 2 / 240
    ns = IC_FULL ? [1000, 3000, 10000, 30000, 100000, 300000, 1000000, 3000000] : [1000, 3000, 10000, 30000, 100000]
    nst = ic_size(100, 200)
    starts = Vector{Vector{Float64}}()
    for _ in 1:nst
        q1 = 0.0 + (2rand(rng) - 1) * 5e-3; q2 = 0.3 + (2rand(rng) - 1) * 5e-3; p2 = 0.2 + (2rand(rng) - 1) * 5e-3
        V = (q1^2 + q2^2) / 2 + q1^2 * q2 - q2^3 / 3
        push!(starts, [q1, q2, sqrt(2 * (1 / 8 - V - p2^2 / 2)), p2])                    # p1 > 0 from H0 = 1/8
    end
    function energy_errors(compensated)
        E = zeros(length(ns), nst)
        for (j, z0) in enumerate(starts)
            zs = advance_to(m, flow, ic_to(Float64, z0), tau, ns; compensated = compensated)
            setprecision(BigFloat, IC_BITS) do
                H0 = ic_H(:t2, ic_to(BigFloat, z0))
                for (i, z) in enumerate(zs)
                    E[i, j] = Float64(ic_H(:t2, ic_to(BigFloat, collect(z))) - H0)
                end
            end
        end
        return E
    end
    Ec = energy_errors(true); Eu = energy_errors(false)
    rms(E) = [sqrt(sum(abs2, E[i, :]) / nst) for i in axes(E, 1)]
    slope(r, lo) = (idx = findall(>=(lo), ns); (log(r[idx[end]]) - log(r[idx[1]])) / (log(ns[idx[end]]) - log(ns[idx[1]])))
    rc = rms(Ec); ru = rms(Eu)
    @info "AC15 reported: rms of the energy error with and without compensated summation" ns rc ru
    # (c) without compensation the law E ~ n^(1/2) holds in the explicit composition (Hairer2008 Fig. 2.1); package tolerance 0.2 (100 starts), full 0.1
    @test abs(slope(ru, 1000) - 0.5) <= (IC_FULL ? 0.1 : 0.2)
    # compensated summation reduces the end value by a large factor (spec: at least 10 at n = 3e6; drafting 300)
    @test ru[end] >= 10 * rc[end]
    # (b) the mean over the starts is zero within 3 standard errors, both forms
    for E in (Ec, Eu)
        se = sqrt(sum(abs2, E[end, :] .- sum(E[end, :]) / nst) / (nst - 1)) / sqrt(nst)
        @test abs(sum(E[end, :]) / nst) <= 3 * se
    end
    # (a) compensated slope within 0.1 over [1e4, 3e6]: needs the full size (the rms sits at the evaluation floor below n = 1e4 and near it up to 1e5)
    if IC_FULL
        idx = findall(>=(10000), ns)
        @test abs(slope(rc, 10000) - 0.5) <= 0.1
    end
end

# D113 (F1 of the stage 04 report). Criterion 15(d), Gauss s = 4 in the implementation of Hairer2008 Sect. 3, compensated summation: the mean of E(n) over the starts must be zero
# within 3 standard errors. Before D113 the stage sums and the update sum were plain sums and the mean was -99 standard errors at n = 3e6 (-14 at 1e5 with 100 starts, -5 at 1e4).
# Package size: 100 starts and n = 1e5 (about 5 s); the full size of the spec (200 starts, n up to 3e6, slope, negative control) stays in verify/integrator-catalog.jl, section AC15D.
@testset "AC15d Gauss s = 4 with compensated summation: no systematic energy drift (package size)" begin
    _ = (advance_to, Gauss, HenonHeilesFlow)
    rng = StableRNG(ic_seed(15, 3))
    flow = HenonHeilesFlow(); m = Gauss(Float64, 4); tau = 2 / 240
    ns = [10000, 30000, 100000]; nst = ic_size(100, 200)
    starts = Vector{Vector{Float64}}()
    for _ in 1:nst
        q1 = 0.0 + (2rand(rng) - 1) * 5e-3; q2 = 0.3 + (2rand(rng) - 1) * 5e-3; p2 = 0.2 + (2rand(rng) - 1) * 5e-3
        V = (q1^2 + q2^2) / 2 + q1^2 * q2 - q2^3 / 3
        push!(starts, [q1, q2, sqrt(2 * (1 / 8 - V - p2^2 / 2)), p2])
    end
    E = zeros(length(ns), nst)
    Threads.@threads for j in 1:nst
        zs = advance_to(m, flow, ic_to(Float64, starts[j]), tau, ns; compensated = true)
        setprecision(BigFloat, IC_BITS) do
            H0 = ic_H(:t2, ic_to(BigFloat, starts[j]))
            for (i, z) in enumerate(zs)
                E[i, j] = Float64(ic_H(:t2, ic_to(BigFloat, collect(z))) - H0)
            end
        end
    end
    mu = vec(sum(E, dims = 2)) ./ nst
    se = [sqrt(sum(abs2, E[i, :] .- mu[i]) / (nst - 1)) / sqrt(nst) for i in eachindex(ns)]
    rms = [sqrt(sum(abs2, E[i, :]) / nst) for i in eachindex(ns)]
    @info "AC15d reported: rms and mean over the starts in standard errors" ns rms mean_se = mu ./ se
    @test all(isfinite, E)
    @test abs(mu[end] / se[end]) <= 3
    @test abs(mu[end - 1] / se[end - 1]) <= 3
    # the cause (D113): the update sum and the stage sums are accurate to one rounding; the plain sums are off by several units in the last place
    F = ntuple(_ -> SVector(0.0, 0.0, 0.0, 0.0), 4); worst_u = 0.0; worst_s = 0.0
    rng2 = StableRNG(ic_seed(15, 4)); h = tau
    for _ in 1:ic_size(500, 2000)
        F = ntuple(_ -> SVector{4, Float64}(2 .* rand(rng2, 4) .- 1), 4)
        u = Chaoslab.Integrators._update_sum(m, F, h); st = Chaoslab.Integrators._stage_sum(m, F, h, 2)
        setprecision(BigFloat, IC_BITS) do
            ref_u = [BigFloat(h) * sum((ic_big(m.bhi[j]) + ic_big(m.blo[j])) * ic_big(F[j][c]) for j in 1:4) for c in 1:4]
            ref_s = [BigFloat(h) * sum((ic_big(m.Ahi[2, j]) + ic_big(m.Alo[2, j])) * ic_big(F[j][c]) for j in 1:4) for c in 1:4]
            for c in 1:4
                worst_u = max(worst_u, Float64(abs(ic_big(u[c]) - ref_u[c]) / abs(ref_u[c])))
                worst_s = max(worst_s, Float64(abs(ic_big(st[c]) - ref_s[c]) / abs(ref_s[c])))
            end
        end
    end
    @info "AC15d reported: worst relative error of the update sum and of a stage sum" worst_u worst_s
    @test worst_u <= 1.01 * ic_eps(Float64) && worst_s <= 1.01 * ic_eps(Float64)
end

@testset "AC16 no allocation after warm-up, force-evaluation counts" begin
    _ = (advance, tangent_advance, CountingFlow, kick_count, reset_counts!, force_evaluations, isolated_force_evaluations, advance_n)
    for T in (Float64, Double64), M in IC_ALLX, kind in (:f1_1, :f1_3, :f4)
        m = M.mk(T); flow = ic_flow(T, kind); z = ic_to(T, IC_START[kind]); d = ic_dim(kind)
        @test ic_alloc(advance, m, flow, z, T(0.05)) == 0
        W = SMatrix{d, d, T}(Matrix{T}(I, d, d))
        @test ic_alloc(tangent_advance, m, flow, z, W, T(0.05)) == 0
    end
    for M in IC_EXPL
        m = M.mk(Float64); c, d = coefficients(m)
        cf = CountingFlow(ic_flow(Float64, :f4)); z = ic_to(Float64, IC_START[:f4])
        advance(m, cf, z, 0.05)
        @test kick_count(cf) == isolated_force_evaluations(m) == ic_expected_evals(c, d)     # one step alone: measured = derived from the stored table
        # D113: the metadata is the cost per step of a long run, first step excluded (kicks merge across step boundaries)
        reset_counts!(cf); advance_n(m, cf, z, 0.05, 2); k2 = kick_count(cf)
        reset_counts!(cf); advance_n(m, cf, z, 0.05, 52); k52 = kick_count(cf)
        @test (k52 - k2) == 50 * force_evaluations(m) == 50 * ic_expected_cyclic(c, d)
        @test force_evaluations(m) <= isolated_force_evaluations(m) <= force_evaluations(m) + 1
    end
    # D113: corrected methods count the two corrector gradients per step on top of the isolated cost of the base
    for M in IC_PC
        m = M.mk(Float64)
        cf = CountingFlow(ic_flow(Float64, :f4)); z = ic_to(Float64, IC_START[:f4])
        advance_n(m, cf, z, 0.05, 2); k2 = kick_count(cf)
        reset_counts!(cf); advance_n(m, cf, z, 0.05, 52); k52 = kick_count(cf)
        @test (k52 - k2) == 50 * force_evaluations(m)
        @test force_evaluations(m) == isolated_force_evaluations(m.base) + 2
    end
    # D113: the named cases of the review (S9, leapfrog KDK, SBAB3)
    for (m, per) in ((CalvoSanzSerna4(Float64), 4), (Leapfrog(Float64, :KDK), 1), (Leapfrog(Float64, :DKD), 1), (SBAB(Float64, 3), 3), (ABA864(Float64), 7))
        cf = CountingFlow(ic_flow(Float64, :f4)); z = ic_to(Float64, IC_START[:f4])
        advance_n(m, cf, z, 0.05, 1000)
        @test force_evaluations(m) == per
        @test abs(kick_count(cf) / 1000 - per) <= 0.0101
    end
    # long-run form: adjacent kicks of consecutive steps share one evaluation (S2: one per step; KDK: n + 1 for n steps)
    cf = CountingFlow(ic_flow(Float64, :f4)); z = ic_to(Float64, IC_START[:f4])
    advance_n(Leapfrog(Float64, :DKD), cf, z, 0.05, 100)
    @test kick_count(cf) == 100
    reset_counts!(cf)
    advance_n(Leapfrog(Float64, :KDK), cf, z, 0.05, 100)
    @test kick_count(cf) == 101
    @test force_evaluations(Leapfrog(Float64, :DKD)) == 1 && force_evaluations(Leapfrog(Float64, :KDK)) == 1      # D113: per step in a long run
    @test isolated_force_evaluations(Leapfrog(Float64, :KDK)) == 2
    # merged long run equals the unmerged sequence of steps up to rounding
    zz = z
    for _ in 1:100
        zz = advance(Leapfrog(Float64, :KDK), ic_flow(Float64, :f4), zz, 0.05)
    end
    @test maximum(abs, collect(advance_n(Leapfrog(Float64, :KDK), ic_flow(Float64, :f4), z, 0.05, 100)) .- collect(zz)) <= ic_B(Float64, 1.0, 100, 0.05)
end
