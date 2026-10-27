# Recorded lines (no threshold, never a gate; the rli D77 form): the methods of R1 to R8 of the spec, written with @info. Values go to the verification report and the study. Nothing is asserted beyond
# finiteness (that the estimator returns). Default sizes are small (D38); WB_FULL_RUN selects the sizes of the spec (R1: 1000 starts per class and N up to 64000: minutes, stage 04 only).
# Not written here: the cross-tabulation of R1 with L1, FLI, MEGNO, RLI, stretching and the naff-fma D (a study ensemble, stage 04), the device half of R5 (stage 04), R4 (the spec names no panel of
# sticky orbits: plan, Items for the owner).

@testset "R1 class separation on the panel" begin
    _ = (WeightedBirkhoff, CosObservable)
    # K_std = 0.5 regular class and {2 pi, 5, 10}, the five observables, N in {250, 1000, 4000, 16000, 64000} in one pass (H = 5): quantiles (0.1, 0.5, 0.9) of absdig per class and observable at
    # each N; the share of chaotic starts above the regular median or at the floor; the share of regular starts below 5.5 and below 5 (F4); the digits are bimodal or not (the quantile spread).
    # D38: 20 starts per class, N in {250, 1000, 4000}.
    Ns = WB_FULL ? (250, 1000, 4000, 16000, 64000) : (250, 1000, 4000)
    nstart = WB_FULL ? 1000 : 20
    names = ("cos q", "sin q", "cos p", "sin p", "lifted p")
    res = Dict{Symbol, Any}()
    for (cl, sub) in ((:half, 0), (:two_pi, 1), (:five, 6), (:ten, 3))
        c = wb_case(cl, Float64, nstart, sub)
        red = WeightedBirkhoff(; window = Ns, observables = wb_std5())
        dig = zeros(5, length(Ns), nstart)
        fl = falses(5, length(Ns), nstart)
        for j in 1:nstart
            _, _, fin = wb_run((red, wb_Norm()), c.sys, c.x0[j], wb_frame(wb_seed(21, 1), c.ids[j], Float64), 2 * maximum(Ns), 1)
            dig[:, :, j] .= Matrix(fin[1].absdig)
            fl[:, :, j] .= Matrix(fin[1].below_floor)
        end
        res[cl] = (dig = dig, fl = fl)
        @test all(isfinite, dig)
    end
    for (k, nm) in enumerate(names), (h, N) in enumerate(Ns)
        reg = res[:half].dig[k, h, :]
        med = wb_quant(reg, 0.5)
        line = ["regular q10/q50/q90 = $(round.([wb_quant(reg, q) for q in (0.1, 0.5, 0.9)]; digits = 2)), below 5.5: $(round(count(<(5.5), reg) / nstart; digits = 3)), below 5: $(round(count(<(5.0), reg) / nstart; digits = 3))"]
        for cl in (:two_pi, :five, :ten)
            d = res[cl].dig[k, h, :]
            push!(line, "$cl q10/q50/q90 = $(round.([wb_quant(d, q) for q in (0.1, 0.5, 0.9)]; digits = 2)), above regular median or at floor: $(round(count(i -> d[i] > med || res[cl].fl[k, h, i], eachindex(d)) / nstart; digits = 3))")
        end
        @info "R1 $nm N=$N" lines = join(line, " | ")
    end
end

@testset "R2 convergence laws of the five weights" begin
    _ = (weighted_birkhoff, WeightedBirkhoff, CosObservable)
    # Log-log slopes in N of |avg_1| for the weights equal, t (1 - t), sin^2, w^[1], w^[2] on the exact rotation signal cos 2 pi theta (theta_0 = 1/7, alpha = sqrt 2 - 1), Float64 and Double64, against
    # WB11 (the sin^2 spread 2.5 to 4, F7); the N at which |avg_1| first falls to 32 eps (F6 a is the first entry); the median of |Delta| over regular orbits at the same N for w^[1] and w^[2].
    alpha = setprecision(BigFloat, WB_BITS) do
        sqrt(BigFloat(2)) - 1
    end
    th0 = setprecision(BigFloat, WB_BITS) do
        BigFloat(1) / 7
    end
    Ns = (125, 250, 500, 1000, 2000, 4000)
    wts = Dict("equal" => (T, N) -> ones(T, N), "t(1-t)" => (T, N) -> T[(t = T(jj) / T(N); t * (1 - t)) for jj in 0:(N - 1)],
               "sin^2" => (T, N) -> T[sin(T(pi) * T(jj) / T(N))^2 for jj in 0:(N - 1)])
    for T in WB_TYPES
        epsT = wb_eps(T)
        for name in ("equal", "t(1-t)", "sin^2", "w1", "w2")
            errs = Float64[]
            for N in Ns
                cs, _ = wb_rot_big(th0, alpha, 1, 2N)
                f = wb_round_vec(T, cs)
                r = name == "w1" ? weighted_birkhoff(f, N; order = 1) : name == "w2" ? weighted_birkhoff(f, N; order = 2) : weighted_birkhoff(f, N; weights = wts[name](T, N))
                push!(errs, Float64(abs(r.avg1[1])))
            end
            pts = [(log10(N), log10(e)) for (N, e) in zip(Ns, errs) if e > 100 * epsT]
            slope = NaN
            if length(pts) >= 2
                xs, ys = first.(pts), last.(pts)
                slope = sum((xs .- sum(xs) / length(xs)) .* (ys .- sum(ys) / length(ys))) / sum((xs .- sum(xs) / length(xs)) .^ 2)
            end
            nfloor = findfirst(<=(32 * epsT), errs)
            @info "R2 $T weight $name: |avg_1| vs N in $Ns, log-log slope above 100 eps, first N at or below 32 eps" errors = string(round.(errs; sigdigits = 3)) slope = round(slope; digits = 2) nfloor = nfloor === nothing ? nothing : Ns[nfloor]
            @test all(isfinite, errs)
        end
    end
    c = wb_case(:half, Float64, WB_FULL ? 20 : 6, 0)
    for p in (1, 2)
        med = Float64[]
        for N in (250, 1000, 4000)
            red = WeightedBirkhoff(; window = N, observables = (CosObservable(1),), birkhoff_order = p)
            d = [abs(wb_run((red, wb_Norm()), c.sys, c.x0[j], wb_frame(wb_seed(22, 1), c.ids[j], Float64), 2N, 1)[3][1].delta[1]) for j in eachindex(c.x0)]
            push!(med, wb_quant(d, 0.5))
        end
        @info "R2 regular orbits, cos q: median |Delta| at N = 250, 1000, 4000 for w^[$p]" median_delta = string(round.(med; sigdigits = 3))
        @test all(isfinite, med)
    end
end

@testset "R3 observables" begin
    _ = (WeightedBirkhoff, CosObservable, SinObservable, LiftedObservable)
    # Digits under cos, sin of q and p, the lifted p, the wrapped angle (a sawtooth) and the lifted angle (not periodic), per class (K_std = 0.5 regular, 2 pi), N = 1000: median and quartiles of dig and
    # the share of regular starts read as chaotic (dig below 5.5, F4/F5). The spec's dry run: the wrapped angle gives 2.3 to 6.0 digits at N = 250 to 16 000, the lifted angle -2 to -4.
    N = 1000
    nstart = WB_FULL ? 1000 : 20
    obs = (CosObservable(1), SinObservable(1), CosObservable(2), SinObservable(2), LiftedObservable(2, 4), wb_Coord(1), LiftedObservable(1, 3))
    names = ("cos q", "sin q", "cos p", "sin p", "lifted p", "wrapped q", "lifted q")
    for (cl, sub) in ((:half, 0), (:two_pi, 1))
        c = wb_case(cl, Float64, nstart, sub)
        red = WeightedBirkhoff(; window = N, observables = obs)
        dig = zeros(length(obs), nstart)
        for j in 1:nstart
            _, _, fin = wb_run((red, wb_Norm()), c.sys, c.x0[j], wb_frame(wb_seed(23, 1), c.ids[j], Float64), 2N, 1)
            dig[:, j] .= vec(Matrix(fin[1].dig))
        end
        @test all(isfinite, dig)
        for (k, nm) in enumerate(names)
            d = dig[k, :]
            @info "R3 $cl $nm N=$N: dig q25/q50/q75, share below 5.5" quartiles = round.([wb_quant(d, q) for q in (0.25, 0.5, 0.75)]; digits = 2) below_5_5 = round(count(<(5.5), d) / nstart; digits = 3)
        end
    end
end

@testset "R5 cost (CPU)" begin
    _ = (WeightedBirkhoff, CosObservable, reducer_bytes)
    # ns per orbit-step for (K, H) in {(1, 1), (4, 1), (4, 11)} against the stream alone (Benettin, k = 1) and the E2 reader, K_std = 2 pi wrapped, 2000 steps (the largest horizon of H = 11 is 1000),
    # best of 3 runs, one thread. Bytes per orbit: reducer_bytes. The device half is stage 04.
    sys = WrappedStandardMap(2pi)
    x0 = torus_state(SVector(1.0, 0.5))
    W0 = wb_frame(wb_seed(25, 1), 1, Float64)
    n = 2000
    time_chain(chain) = (propagate_orbit(sys, x0, W0, n, 1, chain); minimum(@elapsed(propagate_orbit(sys, x0, W0, n, 1, chain)) for _ in 1:3) / n * 1e9)
    base = time_chain((Benettin(),))
    e2 = time_chain((WeightedPhaseAdvance(; window = 1000, observable = LiftedAngleObservable(1, 3)), Benettin()))
    @info "R5 stream alone (Benettin) and E2 reader, ns per orbit-step" base e2
    for (K, H) in ((1, 1), (4, 1), (4, 11))
        wins = H == 1 ? (1000,) : (16, 32, 64, 128, 256, 400, 500, 600, 700, 800, 1000)
        r = WeightedBirkhoff(; window = wins, observables = ntuple(i -> CosObservable(1, i), K))
        t = time_chain((r, Benettin()))
        @info "R5 reader (K, H) = ($K, $H): ns per orbit-step with Benettin, bytes per orbit" t bytes = reducer_bytes(r, Float64, Val(2), Val(1))
        @test isfinite(t)
    end
end

@testset "R6 precision" begin
    _ = (WeightedBirkhoff,)
    # Float64 against Double64 Delta and digits on the regular class (K_std = 0.5, the five observables, N = 1000; the same starts in each type): median and maximum |Delta|, the share at the floor per
    # type, the Float64 plateau of the digits (F6 c: 14 to 17 because the orbit is rounded, not the sums).
    N = 1000
    nstart = WB_FULL ? 1000 : 20
    for T in WB_TYPES
        c = wb_case(:half, T, nstart, 0)
        red = WeightedBirkhoff(; window = N, observables = wb_std5())
        ds, dg, fl = Float64[], Float64[], Float64[]
        for j in 1:nstart
            _, _, fin = wb_run((red, wb_Norm()), c.sys, c.x0[j], wb_frame(wb_seed(26, 1), c.ids[j], T), 2N, 1)
            f = fin[1]
            append!(ds, Float64.(abs.(vec(Matrix(f.delta)))))
            append!(dg, Float64.(vec(Matrix(f.absdig))))
            append!(fl, Float64.(vec(Matrix(f.below_floor))))
        end
        @test all(isfinite, dg)
        @info "R6 $T regular class: median |Delta|, maximum |Delta|, share at the floor, median / maximum absdig" median_delta = wb_quant(ds, 0.5) max_delta = maximum(ds) at_floor = sum(fl) / length(fl) median_dig = wb_quant(dg, 0.5) max_dig = maximum(dg)
    end
end

@testset "R7 weights: w^[1], w^[2] and the width form" begin
    _ = (weighted_birkhoff, birkhoff_index, CosObservable)
    # w^[1] against w^[2] and the width form g_w(s) = exp(-w/(s (1 - s))) (DuignanMeiss2023 Eq. (22), the constant cancels in the normalised average) for w in {0.5, 1, 2, 5}, through the function with a weight
    # vector, on cos q of the K_std = 0.5 regular class: the median absdig at N in {250, 500, 1000, 2000, 4000}, and the digits as a function of w on a fine grid at N = 1000 (the spikes of
    # DuignanMeiss2023 Fig. 7).
    nstart = WB_FULL ? 100 : 6
    c = wb_case(:half, Float64, nstart, 0)
    Ns = (250, 500, 1000, 2000, 4000)
    sers = [wb_series((CosObservable(1),), wb_states(c.sys, c.x0[j], 2 * maximum(Ns)))[:, 1] for j in 1:nstart]
    width(N, w) = Float64[jj == 0 ? 0.0 : (t = min(jj, N - jj) / N; exp(-w / (t * (1 - t)))) for jj in 0:(N - 1)]
    cases = Pair{String, Any}["w1" => ((f, N) -> weighted_birkhoff(f, N; order = 1)), "w2" => ((f, N) -> weighted_birkhoff(f, N; order = 2))]
    for w in (0.5, 1.0, 2.0, 5.0)
        push!(cases, "g_w, w=$w" => ((f, N) -> weighted_birkhoff(f, N; weights = width(N, w))))
    end
    for (name, fn) in cases
        med = Float64[]
        for N in Ns
            d = Float64[]
            for f in sers
                r = fn(f, N)
                push!(d, birkhoff_index(r.avg1[1], r.avg2[1], r.S[1]).absdig)
            end
            push!(med, wb_quant(d, 0.5))
        end
        @info "R7 weight $name: median absdig at N = $Ns" medians = string(round.(med; digits = 2))
        @test all(isfinite, med)
    end
    spikes = Float64[]
    for w in 0.25:0.25:5.0
        d = [(r = weighted_birkhoff(f, 1000; weights = width(1000, w)); birkhoff_index(r.avg1[1], r.avg2[1], r.S[1]).absdig) for f in sers]
        push!(spikes, wb_quant(d, 0.5))
    end
    @info "R7 median absdig at N = 1000 against the width w on the grid 0.25:0.25:5" digits = string(round.(spikes; digits = 1))
    @test all(isfinite, spikes)
end

@testset "R8 window convention" begin
    _ = (weighted_birkhoff, birkhoff_index, CosObservable)
    # The protocol pair (iterates 1..N and N+1..2N: the samples s_1..s_2N of the handoff series) against the Eq. (9) pair of SanderMeiss2020 (iterates 0..N-1 and N..2N-1, the second window started
    # at f^N(z): the samples s_0..s_{2N-1}, one sample earlier): digits of cos q at N = 1000, median and maximum of |dig_protocol - dig_Eq9| on the regular class and the chaotic panel (K_std = 2 pi).
    N = 1000
    nstart = WB_FULL ? 1000 : 20
    for (cl, sub) in ((:half, 0), (:two_pi, 1))
        c = wb_case(cl, Float64, nstart, sub)
        diffs = Float64[]
        for j in 1:nstart
            sts = wb_states(c.sys, c.x0[j], 2N)
            f = wb_series((CosObservable(1),), vcat([c.x0[j]], sts))[:, 1]          # s_0 .. s_2N
            a = weighted_birkhoff(f[2:(2N + 1)], N)                                    # protocol pair
            b = weighted_birkhoff(f[1:(2N)], N)                                        # Eq. (9) pair
            push!(diffs, abs(birkhoff_index(a.avg1[1], a.avg2[1], a.S[1]).absdig - birkhoff_index(b.avg1[1], b.avg2[1], b.S[1]).absdig))
        end
        @test all(isfinite, diffs)
        @info "R8 $cl: |dig_protocol - dig_Eq9| median / maximum over $nstart starts" median = wb_quant(diffs, 0.5) maximum = maximum(diffs)
    end
end
