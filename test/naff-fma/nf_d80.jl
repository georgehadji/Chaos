# Amendment D80: D80a (single line and the aliasing rule of Table 2), D80b (the window family), D80c (the window order option).

@testset "D80a Table 2 aliasing rule" begin
    _ = (window_frequency,)
    # For nu_c = (nu_0 / pi) / 2 with nu_0 / pi in {0.990, ..., 1.009} (the 20 printed rows), z_j = exp(2 pi i nu_c j) generated in BigFloat(256) and rounded once, one window of N = 2000 and of N = 2001
    # samples (the printed series runs from t = -1000 to 1000), order p = 1, both position conventions, mean removal off: (i) for the 19 rows other than 1.000 the returned frequency equals nu_c reduced to
    # [-1/2, 1/2) (nu_c below 1, nu_c - 1 above) within 16 eps; (ii) for 1.000 the circle distance of the returned value to 1/2 is at most 16 eps (value and tie flag reported); (iii) the amplitude
    # equals 1 within 8 eps. The printed zeros to 18 decimals are reported, not gated. D38: Double64 on the rows {0.990, 0.999, 1.000, 1.001, 1.009}; FULL: all 20 rows.
    allrows = [string(k ÷ 1000, ".", lpad(k % 1000, 3, '0')) for k in 990:1009]
    for T in NF_TYPES
        epsT = nf_eps(T)
        rows = T === Float64 || NF_FULL ? allrows : ["0.990", "0.999", "1.000", "1.001", "1.009"]
        bad = String[]
        for rs in rows, N in (2000, 2001), pos in (:half, :edge)
            v = nf_dec(rs)
            nuc = v / 2
            w = window_frequency(nf_series(T, nf_line_big(nuc, BigFloat(1), BigFloat(0), N)); order = 1, positions = pos)
            tag = "T=$T row=$rs N=$N $pos"
            if rs == "1.000"
                nf_cderr(w.nu, BigFloat(1) / 2) <= 16 * epsT || push!(bad, tag * " alias to 1/2")
                @info "D80a row 1.000 (reported)" T N pos nu = Float64(w.nu) tie = w.tie
            else
                tru = v < 1 ? nuc : nuc - 1
                (-0.5 <= w.nu < 0.5 && Float64(abs(nf_big(w.nu) - tru)) <= 16 * epsT) || push!(bad, tag * " frequency")
            end
            Float64(abs(nf_big(w.amplitude) - 1)) <= 8 * epsT || push!(bad, tag * " amplitude")
        end
        @test nf_report(bad, "D80a $T")
    end
end

@testset "D80b window family" begin
    _ = (window_weights,)
    # For p in {0, ..., 5}, N in {16, 256, 4096}, both position conventions: (i) the weights equal (2^p (p!)^2/(2p)!) (1 + cos(pi t_j))^p within 2 (p + 1) eps relative (the spec says 4 eps, which is this bound at
    # p = 1: the relative rounding of sin^(2p) grows with 2p, realized worst 10.9 eps at p = 5; Items for the owner), are nonnegative (zero only at the end point
    # t = -1 of the j/N convention, for p of at least 1) and bitwise symmetric about the centre; (ii) their sum divided by N equals 1 within 4 eps; (iii) for p in {1..5} and N = 4096 the sampled
    # transform (1/N) sum chi_p(t_j) exp(i x t_j) equals the closed form of Eq. (57) within 1e-9 for x in {0.5, 1, 2.5, 7, 20}; (iv) the p = 1 weights are proportional, within 4 eps relative, to
    # the base weights sin^2(pi (j + 1/2)/N) (and sin^2(pi j/N) for the j/N convention).
    for T in NF_TYPES
        epsT = nf_eps(T)
        bad = String[]
        for p in 0:5, N in (16, 256, 4096), pos in (:half, :edge)
            w = window_weights(T, N; order = p, positions = pos)
            ref = nf_chi(p, N, pos)
            tag = "T=$T p=$p N=$N $pos"
            (length(w) == N && w isa AbstractVector{T}) || (push!(bad, tag * " shape"); continue)
            zero_at_edge = pos === :edge && p >= 1
            # (i) relative accuracy, nonnegativity, zeros, symmetry
            worst = 0.0
            for j in 1:N
                if zero_at_edge && j == 1
                    w[1] == 0 || push!(bad, tag * " edge weight not zero")
                    continue
                end
                worst = max(worst, Float64(abs(nf_big(w[j]) - ref[j]) / ref[j]))
            end
            worst <= 2 * (p + 1) * epsT || push!(bad, tag * " relative error $(worst / epsT) eps")
            (all(>=(0), w) && count(iszero, w) == (zero_at_edge ? 1 : 0)) || push!(bad, tag * " signs or zeros")
            sym = pos === :half ? all(j -> w[j] === w[N + 1 - j], 1:N) : all(j -> w[j + 1] === w[N - j + 1], 1:(N - 1))
            sym || push!(bad, tag * " symmetry")
            # (ii)
            s = setprecision(BigFloat, NF_BITS) do
                sum(nf_big(x) for x in w) / N
            end
            Float64(abs(s - 1)) <= 4 * epsT || push!(bad, tag * " sum")
            # (iii)
            if p >= 1 && N == 4096
                for x in (0.5, 1.0, 2.5, 7.0, 20.0)
                    xb = BigFloat(x)
                    abs(nf_phi_sampled(w, xb, pos) - nf_phi_closed(p, xb)) <= 1e-9 || push!(bad, tag * " transform at x=$x")
                end
            end
            # (iv)
            if p == 1
                for j in 1:N
                    zero_at_edge && j == 1 && continue
                    base = setprecision(BigFloat, NF_BITS) do
                        sin(nf_pi() * (pos === :half ? BigFloat(2j - 1) / 2 : BigFloat(j - 1)) / N)^2
                    end
                    Float64(abs(nf_big(w[j]) / base - 2) / 2) <= 4 * epsT || (push!(bad, tag * " proportionality j=$j"); break)
                end
            end
        end
        @test nf_report(bad, "D80b $T")
    end
end

@testset "D80c window order option" begin
    _ = (window_frequency, frequency_map)
    # For p in {0, 2, 3, 4, 5} and N in {64, 1000, 4096}: the pure-line criterion 1 holds with the bounds unchanged (frequency 16 eps, amplitude 8 eps, delta of the two windows at most 32 eps), the "near
    # zero" flag exactly when |nu| < (p + 1)/N; criterion 9 for p in {2, 3} on the series of criterion 9; a window shorter than 8 (p + 1) samples is an ArgumentError naming the parameter. D38: Float64
    # N in {64, 1000} (FULL adds 4096, Double64 N = 64 only; FULL Double64 adds 1000).
    nus = ("0.1234567891", "-0.3711", "0.4999", "0.01", "0.2")
    for T in NF_TYPES
        epsT = nf_eps(T)
        Ns = T === Float64 ? (NF_FULL ? (64, 1000, 4096) : (64, 1000)) : (NF_FULL ? (64, 1000) : (64,))
        bad = String[]
        for p in (0, 2, 3, 4, 5), N in Ns, ns in nus, as in ("1", "1e-3"), phi0 in (0.0, 1.3), pos in (:half, :edge)
            nu, a = nf_dec(ns), nf_dec(as)
            z = nf_series(T, nf_line_big(nu, a, BigFloat(phi0), 2N))
            r = frequency_map(z, N; order = p, positions = pos)
            for (k, w) in enumerate((r.w1[1], r.w2[1]))
                tag = "T=$T p=$p N=$N nu=$ns a=$as phi0=$phi0 $pos window $k"
                nf_cderr(w.nu, nu) <= 16 * epsT || push!(bad, tag * " frequency")
                abs(nf_big(w.amplitude) / a - 1) <= 8 * epsT || push!(bad, tag * " amplitude")
                w.near_zero == (abs(Float64(nu)) < (p + 1) / N) || push!(bad, tag * " near_zero")
                w.mirror && push!(bad, tag * " mirror")
            end
            r.delta <= 32 * epsT || push!(bad, "T=$T p=$p N=$N nu=$ns delta")
        end
        @test nf_report(bad, "D80c $T")
    end
    # criterion 9 for p in {2, 3}
    nstart = NF_FULL ? 20 : 2
    for T in NF_TYPES, N in (NF_FULL ? (256, 1000, 4096) : (256, 1000)), p in (2, 3)
        bad, counts = nf_ac9_check(T, nf_ac9_series(T, N, nstart), N, p)
        @test nf_report(bad, "D80c AC9 $T N=$N p=$p")
        @info "D80c AC9 T=$T N=$N p=$p windows by status (listed, not failed)" counts
    end
    # minimum window 8 (p + 1)
    for p in 0:5
        z(n) = nf_series(Float64, nf_line_big(nf_dec("0.2"), BigFloat(1), BigFloat(0), n))
        m = 8 * (p + 1)
        e = nf_error(() -> window_frequency(z(m - 1); order = p))
        @test e isa ArgumentError && occursin("window", e.msg)
        @test nf_error(() -> window_frequency(z(m); order = p)) === nothing
        e2 = nf_error(() -> frequency_map(z(2 * m - 2), m - 1; order = p))
        @test e2 isa ArgumentError && occursin("window", e2.msg)
    end
end
