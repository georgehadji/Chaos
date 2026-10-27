# AC1 to AC3: pure line, two lines, known drift (the exact signals of the windowed Fourier estimator E1 and, for the drift, of the weighted phase advance E2).
# Every testset starts with a line naming the new symbols it uses (UndefVarError when missing). Tolerances are the spec's (16 eps frequency, 8 eps amplitude, 32 eps pair distance).

@testset "AC1 pure line" begin
    _ = (frequency_map,)
    # z_j = a exp(i (2 pi nu j + phi0)), mean removal off, both sample-position conventions, p = 1: every window has dist(nu_est, nu) <= 16 eps, amplitude within 8 eps relative, the "near zero" flag set
    # exactly when |nu| < 2/N, the "mirror" flag clear, and the pair distance of the two windows at most 32 eps. D38: N in {64, 256, 1000, 4096} for Float64 and {64, 256} for Double64 (spec: up to
    # 10001 and 4096).
    nus = ("0.1234567891", "-0.3711", "0.4999", "0.01", "0.2")
    for T in NF_TYPES
        epsT = nf_eps(T)
        Ns = T === Float64 ? (NF_FULL ? (64, 256, 1000, 4096, 10001) : (64, 256, 1000, 4096)) : (NF_FULL ? (64, 256, 1000, 4096) : (64, 256))
        bad = String[]
        for N in Ns, ns in nus, as in ("1", "1e-3"), phi0 in (0.0, 1.3)
            nu, a = nf_dec(ns), nf_dec(as)
            z = nf_series(T, nf_line_big(nu, a, BigFloat(phi0), 2N))
            for pos in (:half, :edge)
                r = frequency_map(z, N; positions = pos)
                for (k, w) in enumerate((r.w1[1], r.w2[1]))
                    tag = "T=$T N=$N nu=$ns a=$as phi0=$phi0 $pos window $k"
                    nf_cderr(w.nu, nu) <= 16 * epsT || push!(bad, tag * " frequency")
                    setprecision(BigFloat, NF_BITS) do
                        abs(nf_big(w.amplitude) / a - 1) <= 8 * epsT || push!(bad, tag * " amplitude")
                    end
                    w.near_zero == (abs(Float64(nu)) < 2 / N) || push!(bad, tag * " near_zero")
                    w.mirror && push!(bad, tag * " mirror")
                    w.status === :ok || push!(bad, tag * " status")
                end
                r.delta <= 32 * epsT || push!(bad, "T=$T N=$N nu=$ns a=$as phi0=$phi0 $pos delta")
            end
        end
        @test nf_report(bad, "AC1 $T")
    end
end

@testset "AC2 two lines" begin
    _ = (frequency_map,)
    # z = exp(i 2 pi nu_1 j) + a_2 exp(i 2 pi nu_2 j), a_2 in {0.5, 0.1, 0.01}, nu_1 in {0.1234567891, -0.3}, separation in {8/N, 20/N, 0.05, 0.3}, and the same with the amplitudes exchanged:
    # (i) the frequency equals the BigFloat(256) stationary point of the same |phi|^2 found from the same coarse maximum within 16 eps; (ii) the returned line is the stronger one (distance to it
    # below 1/(2N)). The error against the true frequency of the stronger line is reported, not gated (R2). D38: N in {256, 1000} (Float64) and {256} (Double64); spec adds 4096.
    for T in NF_TYPES
        epsT = nf_eps(T)
        Ns = T === Float64 ? (NF_FULL ? (256, 1000, 4096) : (256, 1000)) : (NF_FULL ? (256, 1000, 4096) : (256,))
        bad = String[]
        counts = Dict(:ok => 0, :tie => 0, :nonunique => 0)
        worst_true = 0.0
        for N in Ns, a2s in ("0.5", "0.1", "0.01"), n1s in ("0.1234567891", "-0.3"), sep in (BigFloat(8) / N, BigFloat(20) / N, nf_dec("0.05"), nf_dec("0.3")), swap in (false, true)
            n1 = nf_dec(n1s)
            n2 = n1 + sep
            a2 = nf_dec(a2s)
            specs = swap ? [(n1, a2, BigFloat(0)), (n2, BigFloat(1), BigFloat(0))] : [(n1, BigFloat(1), BigFloat(0)), (n2, a2, BigFloat(0))]
            strong = swap ? n2 : n1
            z = nf_series(T, nf_lines_big(specs, 2N))
            r = frequency_map(z, N)
            for (k, (w, zw)) in enumerate(((r.w1[1], z[1:N]), (r.w2[1], z[(N + 1):(2N)])))
                tag = "T=$T N=$N a2=$a2s nu1=$n1s sep=$(Float64(sep)) swap=$swap window $k"
                st, _ = nf_compare_window(T, zw, w, 1, :half, 8)
                st in (:ok, :tie, :nonunique) || push!(bad, tag * " $st")
                haskey(counts, st) && (counts[st] += 1)
                e = nf_cderr(w.nu, strong)
                e < 1 / (2N) || push!(bad, tag * " stronger line")
                worst_true = max(worst_true, e)
            end
        end
        @test nf_report(bad, "AC2 $T")
        @info "AC2 $T: windows by status (listed, not failed: tie, nonunique); worst error against the stronger line's true frequency (R2, reported)" counts worst_true
    end
end

@testset "AC3 known drift" begin
    _ = (frequency_map, weighted_phase_frequency, phase_increments, circle_distance)
    # z_j = exp(2 pi i (nu_0 j + alpha j^2 / 2)), nu_0 = 0.1234567891, alpha N^2 in {+-0.01, +-0.1, +-1, +-5}, windows [0, N) and [N, 2N): for E1 (both position conventions) and for E2 (phase advance,
    # exponential weight) the signed circle difference nu_2 - nu_1 equals alpha N within 16 eps and delta equals |alpha| N within 32 eps. Also (derived, NF11): the frequency of window 1 is the
    # instantaneous frequency of the window centre, nu_0 + alpha (N - 1)/2 for the positions (j + 1/2)/N and nu_0 + alpha N/2 for j/N (E2: nu_0 + alpha (N + 1)/2), within 16 eps.
    # D38: N in {1000, 4000} for Float64 and {1000} for Double64.
    nu0 = nf_dec("0.1234567891")
    for T in NF_TYPES
        epsT = nf_eps(T)
        Ns = T === Float64 ? (1000, 4000) : (1000,)
        bad = String[]
        for N in Ns, aN2 in (-5, -1, -0.1, -0.01, 0.01, 0.1, 1, 5)
            alpha = BigFloat(aN2) / N^2
            zb = nf_chirp_big(nu0, alpha, 2N + 1)
            z = nf_series(T, zb)
            aN = alpha * N
            for pos in (:half, :edge)
                r = frequency_map(z, N; positions = pos)
                tag = "T=$T N=$N aN2=$aN2 $pos"
                n1, n2 = nf_big(r.w1[1].nu), nf_big(r.w2[1].nu)
                d = setprecision(BigFloat, NF_BITS) do
                    x = n2 - n1
                    x - round(x)
                end
                Float64(abs(d - aN)) <= 16 * epsT || push!(bad, tag * " E1 difference")
                Float64(abs(nf_big(r.delta) - abs(aN))) <= 32 * epsT || push!(bad, tag * " E1 delta")
                c = pos === :half ? BigFloat(N - 1) / 2 : BigFloat(N) / 2
                nf_cderr(r.w1[1].nu, nu0 + alpha * c) <= 16 * epsT || push!(bad, tag * " E1 centre")
            end
            nu1, nu2 = weighted_phase_frequency(phase_increments(z), N)
            tag = "T=$T N=$N aN2=$aN2 E2"
            d = setprecision(BigFloat, NF_BITS) do
                x = nf_big(nu2) - nf_big(nu1)
                x - round(x)
            end
            Float64(abs(d - aN)) <= 16 * epsT || push!(bad, tag * " difference")
            Float64(abs(nf_big(circle_distance(nu1, nu2)) - abs(aN))) <= 32 * epsT || push!(bad, tag * " delta")
            nf_cderr(nu1, nu0 + alpha * BigFloat(N + 1) / 2) <= 16 * epsT || push!(bad, tag * " centre")
        end
        @test nf_report(bad, "AC3 $T")
    end
end
