# AC13 to AC16: long runs (bounded energy, angular momentum, Brouwer's law) and allocation / force-evaluation counts.
# Package sizes follow D38 (3 minute ceiling): runs of 1e4 steps instead of 1e5 (AC13, AC14), Brouwer at n up to 1e5 and 100 starts (spec: 3e6 and 200).
# verify/integrator-catalog.jl sets IC_FULL_RUN = true for the full sizes.

# energy samples of a run: H of the stored state in 256 bits at `nchk` equally spaced step counts
function ic_energy_series(m, kind, flow, z0, tau, nsteps, nchk; eps = 1.0)
    ns = [round(Int, j * nsteps / nchk) for j in 1:nchk]
    zs = advance_to(m, flow, ic_to(typeof(tau), z0), tau, ns)
    return setprecision(BigFloat, IC_BITS) do
        H0 = ic_H(kind, ic_to(BigFloat, z0); eps = eps)
        [Float64(abs(ic_H(kind, ic_to(BigFloat, collect(z)); eps = eps) - H0)) for z in zs]
    end
end
ic_lsslope(y) = (n = length(y); x = collect(1.0:n); xm = sum(x) / n; sum((x .- xm) .* (y .- sum(y) / n)) / sum((x .- xm) .^ 2))

@testset "AC13 energy bounded, no secular drift" begin
    _ = (advance_to, force_evaluations)
    nsteps = ic_size(10000, 100000); nchk = 1000
    cases = (("F1", :f1_3, IC_START[:f1_3], 1.0), ("F2e0.1", :f2, ic_kepler_state(0.1, 0.7), 1.0), ("F2e0.5", :f2, IC_START[:f2], 1.0), ("F4", :f4, IC_START[:f4], 1.0))
    for M in IC_S
        M.id in ("S1DK", "S1KD", "S4", "S3") && continue                                # criterion 13: S2 and S5 to S9, S11 to S13
        tau = (startswith(M.id, "S2") || M.id == "S5" || M.id == "S8" || startswith(M.id, "S12")) ? 1 / 20 : 1 / 50   # criterion 13: S2, S5, S8, S12 at 1/20
        for (name, kind, z0, eps) in cases
            E = ic_energy_series(M.mk(Float64), kind, ic_flow(Float64, kind; eps = eps), z0, tau, nsteps, nchk; eps = eps)
            half = nchk ÷ 2; Mn = ic_mmax(z0)
            @test maximum(E[(half + 1):end]) <= 2 * maximum(E[1:half]) + 16 * ic_eps(Float64) * Mn
            @test abs(ic_lsslope(E[(half + 1):end]) * half) <= 0.1 * maximum(E)
        end
    end
    for M in IC_P                                                                       # near-integrable rows on F4 with eps = 1e-2
        E = ic_energy_series(M.mk(Float64), :f4, ic_flow(Float64, :f4; eps = 1e-2), [1.0, 0.3], 1 / 50, nsteps, nchk; eps = 1e-2)
        half = nchk ÷ 2
        @test maximum(E[(half + 1):end]) <= 2 * maximum(E[1:half]) + 16 * ic_eps(Float64) * 1.0
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

@testset "AC16 no allocation after warm-up, force-evaluation counts" begin
    _ = (advance, tangent_advance, CountingFlow, kick_count, reset_counts!, force_evaluations, advance_n)
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
        @test kick_count(cf) == force_evaluations(m) == ic_expected_evals(c, d)         # metadata = measured = derived from the stored table
    end
    # long-run form: adjacent kicks of consecutive steps share one evaluation (S2: one per step; KDK: n + 1 for n steps)
    cf = CountingFlow(ic_flow(Float64, :f4)); z = ic_to(Float64, IC_START[:f4])
    advance_n(Leapfrog(Float64, :DKD), cf, z, 0.05, 100)
    @test kick_count(cf) == 100
    reset_counts!(cf)
    advance_n(Leapfrog(Float64, :KDK), cf, z, 0.05, 100)
    @test kick_count(cf) == 101
    @test force_evaluations(Leapfrog(Float64, :DKD)) == 1 && force_evaluations(Leapfrog(Float64, :KDK)) == 2
    # merged long run equals the unmerged sequence of steps up to rounding
    zz = z
    for _ in 1:100
        zz = advance(Leapfrog(Float64, :KDK), ic_flow(Float64, :f4), zz, 0.05)
    end
    @test maximum(abs, collect(advance_n(Leapfrog(Float64, :KDK), ic_flow(Float64, :f4), z, 0.05, 100)) .- collect(zz)) <= ic_B(Float64, 1.0, 100, 0.05)
end
