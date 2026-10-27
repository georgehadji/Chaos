# AC4 to AC6: covariances of the estimator, the index with its floor, the circle distance and the frequency vector.

@testset "AC4 covariances (a) to (d)" begin
    _ = (frequency_map,)
    # Series: a pure line, a two-line signal and K_std = 0.5 regular-class series of exp(i q) (N = 1000; D38: 6 of the 20 starts, FULL 20). (a) multiplication by exp(i 1.3) and by 3.7 changes the
    # frequency by at most 16 eps, by 2^k (k = -20, 7) by 0 (bitwise); (b) conjugation gives -nu within 16 eps; (c) multiplication by exp(2 pi i nu_0 j), nu_0 in {0.1, 0.25, 0.37}, gives nu + nu_0
    # modulo 1 within 16 eps (the products are formed in BigFloat and rounded once, so that the covariance is that of the estimator, not of the test's multiplication); (d) delta is bitwise equal when
    # the two windows are exchanged. Float64 and Double64 (the regular class in Double64 uses 2 starts).
    N = 1000
    nu = nf_dec("0.1234567891")
    for T in NF_TYPES
        epsT = nf_eps(T)
        nreg = T === Float64 ? (NF_FULL ? 20 : 6) : 2
        c = nf_case(:half, T, nreg, 0)
        series = Vector{Complex{T}}[]
        push!(series, nf_series(T, nf_line_big(nu, BigFloat(1), BigFloat(0), 2N)))
        push!(series, nf_series(T, nf_lines_big([(nu, BigFloat(1), BigFloat(0)), (nu + nf_dec("0.05"), nf_dec("0.1"), BigFloat(0))], 2N)))
        for j in 1:nreg
            push!(series, nf_obs(nf_obs_cis, nf_states(c.sys, c.x0[j], 2N)))
        end
        bad = String[]
        for (si, z) in enumerate(series)
            zb = [nf_big(x) for x in z]
            r0 = frequency_map(z, N)
            nus = (r0.w1[1].nu, r0.w2[1].nu)
            tag = "T=$T series $si"
            # (a) unit constant, positive constant, powers of two
            for (name, f, exact) in (("exp(i 1.3)", x -> x * nf_round(T, nf_cis(BigFloat(1.3))), false), ("3.7", x -> x * T(3.7), false), ("2^-20", x -> x * T(2.0)^-20, true), ("2^7", x -> x * T(2.0)^7, true))
                r = frequency_map(map(f, z), N)
                for (k, w) in enumerate((r.w1[1], r.w2[1]))
                    ok = exact ? isequal(w.nu, nus[k]) : nf_cderr(w.nu, nf_big(nus[k])) <= 16 * epsT
                    ok || push!(bad, "$tag $name window $k")
                end
            end
            # (b) conjugation
            rc = frequency_map(conj.(z), N)
            for (k, w) in enumerate((rc.w1[1], rc.w2[1]))
                nf_cderr(w.nu, -nf_big(nus[k])) <= 16 * epsT || push!(bad, "$tag conj window $k")
            end
            # (c) modulation
            for nu0s in ("0.1", "0.25", "0.37")
                nu0 = nf_dec(nu0s)
                mod_b = setprecision(BigFloat, NF_BITS) do
                    w = nf_cis(2 * nf_pi() * nu0)
                    e = one(Complex{BigFloat})
                    out = similar(zb)
                    for j in eachindex(zb)
                        out[j] = zb[j] * e
                        e *= w
                    end
                    out
                end
                rm = frequency_map(nf_series(T, mod_b), N)
                for (k, w) in enumerate((rm.w1[1], rm.w2[1]))
                    nf_cderr(w.nu, nf_big(nus[k]) + nu0) <= 16 * epsT || push!(bad, "$tag shift $nu0s window $k")
                end
            end
            # (d) exchange of the windows
            rx = frequency_map(vcat(z[(N + 1):(2N)], z[1:N]), N)
            isequal(rx.delta, r0.delta) || push!(bad, "$tag exchange")
        end
        @test nf_report(bad, "AC4 $T")
    end
end

@testset "AC4e mirror rule" begin
    _ = (frequency_map,)
    # z_j = cos(2 pi nu j) + e s exp(2 pi i nu j), nu = 0.2, s = +1 in the first window and -1 in the second, N = 1000. e = 0.01: the "mirror" flag is set in both windows (the mirror line has 0.98 of the
    # peak, above 0.9) and the plane distance is the mirror-rule distance, nominally 0; e = 0.2: the flag is clear and the distance is nominally 2 nu (the sign flip is real).
    # The spec states delta <= 32 eps and |delta - 2 nu| <= 32 eps; neither is attainable: the line at -nu leaks into the estimate of the line at +nu (Laskar2003 Eq. (63), order 12 a/((N/2)^4 Omega^3) =
    # 1.2e-11 per unit amplitude ratio at Omega = 2 pi 0.4), and the two windows have different amplitude ratios, so the nominal values are only reached to about 1e-12 (2e-15 for e = 0.01 where the ratios
    # nearly cancel, 6e-12 for e = 0.2). Gated instead (Items for the owner): (i) the flags as in the spec; (ii) delta equals the BigFloat reference delta (stationary points of the two windows, same bins, the
    # mirror rule applied as the flags dictate) within 32 eps; (iii) the nominal value within twice the first-order leakage bound plus 32 eps.
    N = 1000
    nu = nf_dec("0.2")
    for T in NF_TYPES
        epsT = nf_eps(T)
        for (es, flag) in (("0.01", true), ("0.2", false))
            e = nf_dec(es)
            w1 = nf_lines_big([(nu, BigFloat(1) / 2 + e, BigFloat(0)), (-nu, BigFloat(1) / 2, BigFloat(0))], N)
            w2 = nf_lines_big([(nu, BigFloat(1) / 2 - e, BigFloat(0)), (-nu, BigFloat(1) / 2, BigFloat(0))], N)
            z = nf_series(T, vcat(w1, w2))
            r = frequency_map(z, N)
            @test r.w1[1].mirror == flag && r.w2[1].mirror == flag
            # (ii) reference delta
            ref = [nf_stat_ref(z[rng], 1, :half, 8; bin = Int(w.bin)) for (rng, w) in ((1:N, r.w1[1]), ((N + 1):(2N), r.w2[1]))]
            @test ref[1].unique && ref[2].unique
            d_ref = setprecision(BigFloat, NF_BITS) do
                d = nf_cd(ref[1].nu, ref[2].nu)
                flag ? min(d, nf_cd(ref[1].nu, -ref[2].nu)) : d
            end
            @test Float64(abs(nf_big(r.delta) - d_ref)) <= 32 * epsT
            # (iii) the nominal value within the leakage bound
            lk = 12 / ((N / 2)^4 * (2 * pi * 0.4)^3)
            a1, a2, b = 0.5 + parse(Float64, es), 0.5 - parse(Float64, es), 0.5
            nominal = flag ? 0.0 : 0.4
            bound = flag ? 2 * lk * abs(b / a1 - a2 / b) : 2 * lk * (b / a1 + a2 / b)
            @test abs(Float64(nf_big(r.delta)) - nominal) <= bound + 32 * epsT
        end
    end
end

@testset "AC5 index and floor" begin
    _ = (frequency_map, frequency_index)
    # (a) pure line, N in {256, 4096} (D38: {256, 1000}): delta <= 32 eps, the "below floor" flag set, D = log10(32 eps) (bitwise given the floor); (b) drift alpha N in {1e-6, 1e-9}: the flag is clear and
    # D equals log10 of the module's own delta within 2 eps relative; (c) delta is returned unclipped (the criterion 3 value in (b), strictly the rounding value below the floor in (a)); (d) two windows
    # of one line with different phases (phi0 and phi0 + 1) give delta <= 32 eps.
    nu = nf_dec("0.1234567891")
    for T in NF_TYPES
        epsT = nf_eps(T)
        floorD = log10(T(32) * effective_eps(T))
        for N in (T === Float64 && NF_FULL ? (256, 4096) : (256, 1000))
            z = nf_series(T, nf_line_big(nu, BigFloat(1), BigFloat(0.5), 2N))
            r = frequency_map(z, N)
            @test r.delta < 32 * epsT && r.below_floor && r.D == floorD                    # delta is returned as computed, strictly below the floor it is compared with (never replaced by it)
            # (d) phase jump between the windows
            zb = nf_line_big(nu, BigFloat(1), BigFloat(0.5), N)
            zb2 = nf_line_big(nu, BigFloat(1), BigFloat(0.5) + 1, 2N)[(N + 1):(2N)]
            zj = nf_series(T, vcat(zb, zb2))
            rj = frequency_map(zj, N)
            @test rj.delta <= 32 * epsT && rj.below_floor
        end
        N = 1000
        for aN in ("1e-6", "1e-9")
            alpha = nf_dec(aN) / N
            z = nf_series(T, nf_chirp_big(nu, alpha, 2N))
            r = frequency_map(z, N)
            @test !r.below_floor
            @test Float64(abs(nf_big(r.delta) - nf_dec(aN))) <= 32 * epsT                                   # the criterion 3 value
            ref = setprecision(BigFloat, NF_BITS) do
                log10(nf_big(r.delta))
            end
            @test Float64(abs(nf_big(r.D) - ref)) <= 2 * epsT * Float64(abs(ref))
        end
    end
end

@testset "AC6 circle distance and the frequency vector" begin
    _ = (frequency_map, circle_distance)
    # (a) lines nu_1 = 0.4999999 and nu_2 = -0.4999999: delta = 2e-7 within 32 eps (not about 1); (b) two planes from two Y5 rotation blocks (theta_1 = 0.7, theta_2 = 2.3, 20 seeded starts, observables
    # q_k + i p_k): the frequency of plane k is -theta_k / (2 pi) reduced within 16 eps in both windows, delta <= 32 eps; (c) two planes with drifts alpha_k N^2 in {0.1, 1}: delta =
    # sqrt((alpha_1 N)^2 + (alpha_2 N)^2) within 32 eps; (d) a plane with a degenerate series gives no delta for the pair and the status names the plane. D38: N = 256 for (b) with 6 starts in Float64 and 2
    # in Double64 (spec: 20 starts); FULL adds N = 1000 and 20 starts.
    for T in NF_TYPES
        epsT = nf_eps(T)
        # (a)
        for N in (256, 1000)
            n1, n2 = nf_dec("0.4999999"), nf_dec("-0.4999999")
            z = nf_series(T, vcat(nf_line_big(n1, BigFloat(1), BigFloat(0), N), nf_line_big(n2, BigFloat(1), BigFloat(0), N)))
            r = frequency_map(z, N)
            @test Float64(abs(nf_big(r.delta) - nf_dec("2e-7"))) <= 32 * epsT
            @test Float64(nf_cd(nf_big(r.w1[1].nu), n1)) <= 16 * epsT && Float64(nf_cd(nf_big(r.w2[1].nu), n2)) <= 16 * epsT
            @test !r.w1[1].mirror && !r.w2[1].mirror
            @test Float64(abs(nf_big(circle_distance(r.w1[1].nu, r.w2[1].nu)) - nf_dec("2e-7"))) <= 32 * epsT
        end
        # (b)
        th = (0.7, 2.3)
        sys = nf_rot_map(T, th...)
        rng = StableRNG(nf_seed(6, 1))
        nstart = T === Float64 ? (NF_FULL ? 20 : 6) : 2
        bad = String[]
        for N in (NF_FULL ? (256, 1000) : (256,)), s in 1:nstart
            x0 = nf_rot_start(T, rng, 4)
            st = nf_states(sys, x0, 2N)
            planes = [nf_obs(nf_obs_plane(k, 2), st) for k in 1:2]
            r = frequency_map(planes, N)
            for k in 1:2
                tru = nf_rot_nu(th[k])
                (nf_cderr(r.w1[k].nu, tru) <= 16 * epsT && nf_cderr(r.w2[k].nu, tru) <= 16 * epsT) || push!(bad, "T=$T N=$N start $s plane $k")
            end
            r.delta <= 32 * epsT || push!(bad, "T=$T N=$N start $s delta")
        end
        @test nf_report(bad, "AC6b $T")
        # (c)
        N = 1000
        a1, a2 = nf_dec("0.1") / N^2, BigFloat(1) / N^2
        planes = [nf_series(T, nf_chirp_big(nf_dec("0.1234567891"), a1, 2N)), nf_series(T, nf_chirp_big(nf_dec("-0.3"), a2, 2N))]
        r = frequency_map(planes, N)
        ref = setprecision(BigFloat, NF_BITS) do
            sqrt((a1 * N)^2 + (a2 * N)^2)
        end
        @test Float64(abs(nf_big(r.delta) - ref)) <= 32 * epsT
        # (d)
        const_plane = fill(complex(T(1), T(0)), 2N)
        r = frequency_map([const_plane, planes[2]], N; mean_removal = true)
        @test r.status[1] === :degenerate && r.status[2] === :ok && isnan(r.delta) && isnan(r.D) && !r.below_floor
        @test r.w1[1].status === :degenerate && isnan(r.w1[1].nu)
    end
end
