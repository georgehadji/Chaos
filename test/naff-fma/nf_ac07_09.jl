# AC7 to AC9: known systems (Y5 rotation block, the K_std = 0 standard map), E2 algebra and the known answer NF15, the maximiser against its independent reference.
# AC7(c) (the series exp(i q) built from the wrapped coordinate) concerns the caller's series builder, not the module: it is gated in nf_selftest.jl (green now by design), see the plan.

@testset "AC7a Y5 rotation block" begin
    _ = (frequency_map, weighted_phase_frequency, phase_increments, frequency_index, circle_distance, WeightedPhaseAdvance, PlanarObservable)
    # z' = exp(-i theta) z (NF14): a single line of frequency -theta/(2 pi) reduced, sign included. Observable q + i p, mean removal off, theta in {0.3, 1.0, 2.5, 3.1}, seeded starts uniform in
    # [-1, 1]^2: E1, E2 (pure function on the phase increments) and E2 through the stream reducer; frequency within 16 eps in both windows, delta at the floor with the flag set.
    # D38: N = 256 with 20 starts (Float64) / 3 (Double64), N = 4096 with 2 starts (FULL: 20 starts, both N, both types).
    for T in NF_TYPES
        epsT = nf_eps(T)
        plan = T === Float64 ? ((256, NF_FULL ? 20 : 20), (4096, NF_FULL ? 20 : 2)) : ((256, NF_FULL ? 20 : 3), (4096, NF_FULL ? 20 : 0))
        bad = String[]
        for theta in (0.3, 1.0, 2.5, 3.1), (N, nstart) in plan
            rng = StableRNG(nf_seed(7, 20 + findfirst(==(theta), (0.3, 1.0, 2.5, 3.1))))
            tru = nf_rot_nu(theta)
            sys = nf_rot_map(T, theta)
            red = WeightedPhaseAdvance(; window = N, observable = PlanarObservable(1, 2))
            for s in 1:nstart
                x0 = nf_rot_start(T, rng, 2)
                st = nf_states(sys, x0, 2N + 1)
                z = nf_obs(nf_obs_planar, st)
                tag = "T=$T theta=$theta N=$N start $s"
                r = frequency_map(z, N)
                (nf_cderr(r.w1[1].nu, tru) <= 16 * epsT && nf_cderr(r.w2[1].nu, tru) <= 16 * epsT) || push!(bad, tag * " E1 frequency")
                (r.delta <= 32 * epsT && r.below_floor) || push!(bad, tag * " E1 floor")
                nu1, nu2 = weighted_phase_frequency(phase_increments(z), N)
                (nf_cderr(nu1, tru) <= 16 * epsT && nf_cderr(nu2, tru) <= 16 * epsT) || push!(bad, tag * " E2 frequency")
                frequency_index(circle_distance(nu1, nu2)).below_floor || push!(bad, tag * " E2 floor")
                _, _, fin = nf_run((red,), sys, x0, nf_frame(nf_seed(7, 2), s, T, 2, 1), 2N + 1, 1)
                f = fin[1]
                (nf_cderr(f.nu1[1], tru) <= 16 * epsT && nf_cderr(f.nu2[1], tru) <= 16 * epsT) || push!(bad, tag * " stream frequency")
                (f.complete && f.below_floor[1] && f.delta[1] <= 32 * epsT) || push!(bad, tag * " stream floor")
            end
        end
        @test nf_report(bad, "AC7a $T")
    end
end

@testset "AC7b K_std = 0 standard map" begin
    _ = (frequency_map, weighted_phase_frequency, WeightedPhaseAdvance, LiftedAngleObservable)
    # NF13: q_j = q_0 + j p_0 modulo 2 pi, exp(i q_j) is a single line of frequency p_0/(2 pi) reduced; the lifted q increases by p_0 per step. p_0 in {0.1, 1.0, 2.0, 3.3, 5.0, 6.1, pi}, random q_0 from
    # StableRNG(seed), wrapped state with the winding counts. E1 on exp(i q) (mean removal off): the reduced p_0/(2 pi) within 16 eps, the two windows at the floor; with mean removal on the window mean
    # of the line is a real DC term, which perturbs the frequency by the leakage of Laskar2003 Eq. (63) (order 12 m/((N/2)^4 (2 pi nu)^3), derived, NOT within 16 eps): gated at 16 eps + 4 times that
    # bound (Items for the owner). E2 on the lifted increments (pure function and stream reducer): p_0/(2 pi) within 16 eps (not reduced), p_0 = pi included. D38: 4 starts, N = 1000
    # (Float64), 2 starts (Double64); FULL: 20 starts, N in {1000, 4096}.
    for T in NF_TYPES
        epsT = nf_eps(T)
        nstart = NF_FULL ? 20 : (T === Float64 ? 4 : 2)
        Ns = NF_FULL ? (1000, 4096) : (1000,)
        rng = StableRNG(nf_seed(7, 10))
        q0s = [T(2pi * rand(rng)) for _ in 1:nstart]
        sys = WrappedStandardMap(nf_K(T, :zero))
        red(N) = WeightedPhaseAdvance(; window = N, observable = LiftedAngleObservable(1, 3))
        bad = String[]
        worst_ratio = 0.0
        for p0 in (0.1, 1.0, 2.0, 3.3, 5.0, 6.1, pi), N in Ns, (s, q0) in enumerate(q0s)
            pT = T(p0)
            x0 = torus_state(SVector{2, T}(q0, pT))
            pt = x0[2]
            tru, truE2 = setprecision(BigFloat, NF_BITS) do
                v = nf_big(pt) / (2 * nf_pi())
                (nf_red(v), v)
            end
            st = nf_states(sys, x0, 2N + 1)
            tag = "T=$T p0=$p0 N=$N start $s"
            z = nf_obs(nf_obs_cis, st[1:(2N)])
            r = frequency_map(z, N)
            (nf_cderr(r.w1[1].nu, tru) <= 16 * epsT && nf_cderr(r.w2[1].nu, tru) <= 16 * epsT) || push!(bad, tag * " E1 frequency")
            (r.delta <= 32 * epsT && r.below_floor) || push!(bad, tag * " E1 floor")
            # mean removal on: derived leakage bound
            rm = frequency_map(z, N; mean_removal = true)
            mabs = setprecision(BigFloat, NF_BITS) do
                th = 2 * nf_pi() * tru
                abs(sin(N * th / 2) / (N * sin(th / 2)))
            end
            B = Float64(192 * mabs / (BigFloat(N)^4 * (2 * nf_pi() * abs(tru))^3))
            for w in (rm.w1[1], rm.w2[1])
                e = nf_cderr(w.nu, tru)
                e <= 16 * epsT + 4 * B || push!(bad, tag * " E1 mean removed: error $e, bound $B")
                B > 0 && (worst_ratio = max(worst_ratio, (e - 16 * epsT) / B))
            end
            # E2 pure and stream
            inc = nf_lifted_incr(st)
            nu1, nu2 = weighted_phase_frequency(inc, N)
            (nf_err(nu1, truE2) <= 16 * epsT && nf_err(nu2, truE2) <= 16 * epsT) || push!(bad, tag * " E2 pure")
            _, _, fin = nf_run((red(N),), sys, x0, nf_frame(nf_seed(7, 11), s, T, 2, 1), 2N + 1, 1)
            f = fin[1]
            (nf_err(f.nu1[1], truE2) <= 16 * epsT && nf_err(f.nu2[1], truE2) <= 16 * epsT) || push!(bad, tag * " E2 stream")
            (f.below_floor[1] && f.delta[1] <= 32 * epsT) || push!(bad, tag * " E2 floor")
        end
        @test nf_report(bad, "AC7b $T")
        @info "AC7b $T: worst (error - 16 eps) / leakage bound with mean removal on (reported; gate is 4)" worst_ratio
    end
end

@testset "AC8 E2 algebra and the known answer" begin
    _ = (birkhoff_weights, weighted_phase_frequency)
    # (a) normalised exponential weights: sum 1 within 4 eps and symmetry under j <-> N - j within 4 eps relative, N in {16, 1000, 10^5}; (b) the weighted average of a constant returns the constant
    # within 2 eps (relative); (c) NF15 (rho = sqrt(2) - 1, alpha = 0.1, beta = 0.2) with increments rounded once: Float64 at N in {500, 1000, 4000}: both windows within 16 eps of rho; Double64 and
    # BigFloat(256) agree with the 512-bit evaluation of the same normalised sum on the same rounded increments within 16 eps_T at N in {500, 1000, 4000, 16000} (D38: up to 4000 for Double64 and BigFloat,
    # FULL adds 16000) and are within 16 eps_T of rho from the smallest N at which the reference sum is within 4 eps_T of rho (that N is reported).
    for T in NF_TYPES
        epsT = nf_eps(T)
        for N in (16, 1000, 100_000)
            w = birkhoff_weights(T, N)
            @test length(w) == N && all(>=(0), w)
            s = setprecision(BigFloat, NF_BITS) do
                sum(nf_big(x) for x in w)
            end
            @test Float64(abs(s - 1)) <= 4 * epsT
            asym = count(j -> abs(w[j + 1] - w[N - j + 1]) > 4 * epsT * max(w[j + 1], w[N - j + 1]), 1:(N - 1))
            @test asym == 0
        end
        # (b)
        N = 1000
        for c in (0.7, -2.3, 3.0)
            cT = T(c)
            nu1, nu2 = weighted_phase_frequency(fill(cT, 2N), N)
            ref = setprecision(BigFloat, NF_BITS) do
                nf_big(cT) / (2 * nf_pi())
            end
            @test Float64(abs(nf_big(nu1) - ref) / abs(ref)) <= 2 * epsT && Float64(abs(nf_big(nu2) - ref) / abs(ref)) <= 2 * epsT
        end
    end
    # (c)
    Nmax = NF_FULL ? 16000 : 4000
    incr, rho = nf_nf15_big(2 * Nmax)
    for T in (Float64, Double64, BigFloat)
        epsT = T === BigFloat ? Float64(setprecision(() -> eps(BigFloat), BigFloat, NF_BITS)) : nf_eps(T)
        Ns = T === Float64 ? (500, 1000, 4000) : (NF_FULL ? (500, 1000, 4000, 16000) : (500, 1000, 4000))
        nstar = nothing
        for N in Ns
            inT = setprecision(BigFloat, NF_BITS) do
                T[nf_round(T, x) for x in incr[1:(2N)]]
            end
            nu = setprecision(BigFloat, NF_BITS) do
                weighted_phase_frequency(inT, N)
            end
            ref = nf_E2_ref(nf_big.(inT), N; bits = 512)
            if T === Float64
                @test Float64(abs(nf_big(nu[1]) - rho)) <= 16 * epsT && Float64(abs(nf_big(nu[2]) - rho)) <= 16 * epsT
            else
                @test Float64(abs(nf_big(nu[1]) - ref[1])) <= 16 * epsT && Float64(abs(nf_big(nu[2]) - ref[2])) <= 16 * epsT
                if nstar === nothing && Float64(abs(ref[1] - rho)) <= 4 * epsT && Float64(abs(ref[2] - rho)) <= 4 * epsT
                    nstar = N
                end
                if nstar !== nothing
                    @test Float64(abs(nf_big(nu[1]) - rho)) <= 16 * epsT && Float64(abs(nf_big(nu[2]) - rho)) <= 16 * epsT
                end
            end
        end
        T === Float64 || @info "AC8c $T: smallest N at which the reference sum is within 4 eps_T of rho (reported for stage 04)" nstar
    end
end

# ---------------------------------------------------------------- AC9

# Series (window means removed in exact arithmetic, rounded once to T) of the criterion 9 set: the K_std = 0.5 regular class and K_std in {2 pi, 5, 10}, observables exp(i q) and q + i p, nstart starts,
# and the signals of criteria 1 to 3 (a line, two lines at 20/N with a_2 = 0.1, a drift alpha N^2 = 1).
function nf_ac9_series(::Type{T}, N::Int, nstart::Int) where {T}
    out = Tuple{String, Vector{Complex{T}}}[]
    for (cl, sub) in ((:half, 0), (:two_pi, 1), (:five, 2), (:ten, 3))
        c = nf_case(cl, T, nstart, sub)
        for j in 1:nstart, (oname, f) in (("cis", nf_obs_cis), ("planar", nf_obs_planar))
            zb = Complex{BigFloat}[nf_big(x) for x in nf_obs(f, nf_states(c.sys, c.x0[j], 2N))]
            push!(out, ("$cl $oname start $j", nf_series(T, nf_demean_big(zb, N))))
        end
    end
    nu = nf_dec("0.1234567891")
    push!(out, ("line", nf_series(T, nf_line_big(nu, BigFloat(1), BigFloat(0.5), 2N))))
    push!(out, ("two lines", nf_series(T, nf_lines_big([(nu, BigFloat(1), BigFloat(0)), (nu + BigFloat(20) / N, nf_dec("0.1"), BigFloat(0))], 2N))))
    push!(out, ("drift", nf_series(T, nf_chirp_big(nu, BigFloat(1) / N^2, 2N))))
    return out
end

# Compare every window of a series set with the BigFloat reference; returns (bad tags, counts by status).
function nf_ac9_check(::Type{T}, sets, N::Int, p::Int) where {T}
    bad = String[]
    counts = Dict(:ok => 0, :tie => 0, :nonunique => 0)
    for (label, z) in sets
        r = frequency_map(z, N; order = p)
        for (k, (w, zw)) in enumerate(((r.w1[1], z[1:N]), (r.w2[1], z[(N + 1):(2N)])))
            st, _ = nf_compare_window(T, zw, w, p, :half, 8)
            haskey(counts, st) ? (counts[st] += 1) : push!(bad, "T=$T N=$N p=$p $label window $k $st")
        end
    end
    return bad, counts
end

@testset "AC9 the maximiser against its independent reference" begin
    _ = (frequency_map,)
    # For every window the frequency equals the BigFloat(256) stationary point of the same |phi|^2 on the same realized series, started from the same coarse maximum, within 16 eps, whenever the reference
    # finds a unique stationary point with negative second derivative in the bracket (the other windows are counted and listed, not failed); the coarse bin equals the reference's bin except ties (a
    # module bin whose |phi| is not below the reference bin's by more than the tie tolerance, 64 eps of Float64), which are counted. D38: 3 starts per class, N in {256, 1000}; FULL: 20 starts, N in {256, 1000, 4096}.
    nstart = NF_FULL ? 20 : 3
    for T in NF_TYPES, N in (NF_FULL ? (256, 1000, 4096) : (256, 1000))
        bad, counts = nf_ac9_check(T, nf_ac9_series(T, N, nstart), N, 1)
        @test nf_report(bad, "AC9 $T N=$N")
        @info "AC9 T=$T N=$N windows by status (listed, not failed)" counts
    end
end
