# AC1 to AC3: algebra of the definition, normalisation interval invariance, K = 0 shear.
# Each testset starts by naming the Indicators symbols it uses, so a missing module errors that testset (UndefVarError) and nothing else.

@testset "AC1 algebra of the definition" begin
    _ = (alignment_indices, sali_pair, sali_floor)
    # 1e4 random pairs of non-unit vectors per type (d = 2; 1e3 more in d = 4: A1 to A4 hold in any dimension), half independent, half
    # near-parallel with a perturbation 10^-u (u in [1, 14]), random signs; constructed pairs u2 = rotation of u1 by 10^-m, m = 1..14, parallel
    # and antiparallel. (a) SALI = min(d-, d+) and both indices within B_meas of the BigFloat(256) values; (b) |d-^2 + d+^2 - 4| <= 16 eps (Amendment D67);
    # (c) |GALI_2 - SALI max(d-, d+)/2| <= 8 eps with GALI_2 = d- d+ / 2 (A4); (d) |SALI - SALI_exact| <= B_meas for the stored inputs;
    # (e) SALI bitwise unchanged by a sign change of either vector and by exchanging them; (f) identical and exactly opposite vectors give 0;
    # (g) exactly orthogonal pairs give sqrt 2 within 4 eps. Negative control at the end.
    rng = StableRNG(sl_seed(1))
    for T in SL_TYPES
        epsT, meas = sl_eps(T), sl_meas(T)
        emax = T === Float64 ? 200 : 60
        worst = (a = 0.0, b = 0.0, c = 0.0, d = 0.0)
        bad_min = 0
        bad_bits = 0
        for i in 1:(10_000 + 1_000)
            dim = i <= 10_000 ? 2 : 4
            e1, e2 = rand(rng, -emax:emax), rand(rng, -emax:emax)
            a = sl_randvec(rng, T, dim, e1)
            if iseven(i)
                b = sl_randvec(rng, T, dim, e2)
            else
                delta = T(10.0^(-1 - 13 * rand(rng)))
                b = a + delta * sl_randvec(rng, T, dim, e1)        # near-parallel; the scale of a is kept
            end
            rand(rng, Bool) && (b = -b)
            dm, dp = alignment_indices(a, b)
            s = sali_pair(a, b)
            bad_min += s == min(dm, dp) ? 0 : 1
            dmx, dpx, sx = sl_sali_exact(a, b)
            worst = (a = max(worst.a, sl_err(dm, dmx) / meas, sl_err(dp, dpx) / meas),
                     b = max(worst.b, Float64(abs(dm^2 + dp^2 - 4)) / (16 * epsT)),
                     c = max(worst.c, Float64(abs(dm * dp / 2 - s * max(dm, dp) / 2)) / (8 * epsT)),
                     d = max(worst.d, sl_err(s, sx) / meas))
            bad_bits += (s === sali_pair(-a, b) && s === sali_pair(a, -b) && s === sali_pair(b, a) &&
                         s === sali_pair(-b, -a)) ? 0 : 1
        end
        @test bad_min == 0                                   # (a) the SALI is the smaller of the two indices
        @test worst.a <= 1                                   # (a) d-, d+ against BigFloat, in units of B_meas
        @test worst.b <= 1                                   # (b)
        @test worst.c <= 1                                   # (c)
        @test worst.d <= 1                                   # (d)
        @test bad_bits == 0                                  # (e)
        # constructed pairs: u2 = R(10^-m) u1 and its antiparallel version, m = 1..14
        for m in 1:14, sgn in (1, -1)
            ang = T(10.0)^(-m)
            u1 = SVector{2, T}(cos(T(0.3)), sin(T(0.3)))
            u2 = sgn * SVector{2, T}(cos(T(0.3) + ang), sin(T(0.3) + ang))
            dm, dp = alignment_indices(u1, u2)
            s = sali_pair(u1, u2)
            _, _, sx = sl_sali_exact(u1, u2)
            @test s == min(dm, dp)
            @test sl_err(s, sx) <= meas
            @test s < 3 * 10.0^(-m)                          # the pair really is that close (parallel and antiparallel alike)
            @test abs(dm^2 + dp^2 - 4) <= 16 * epsT
        end
        # (f) identical and exactly opposite vectors
        for _ in 1:50
            a = sl_randvec(rng, T, 2, rand(rng, -emax:emax))
            @test sali_pair(a, a) == 0 && sali_pair(a, -a) == 0 && sali_pair(-a, a) == 0
        end
        # (g) exactly orthogonal pairs (integer entries scaled by a power of two: the inputs are exact)
        for (p, q) in ((SVector{2, T}(1, 0), SVector{2, T}(0, 1)), (SVector{2, T}(3, 4), SVector{2, T}(-4, 3)),
                       (SVector{2, T}(5, 12), SVector{2, T}(-12, 5)), (SVector{2, T}(0, 7), SVector{2, T}(7, 0)))
            for e in (-100, 0, 40)
                @test abs(sali_pair(ldexp.(p, e), ldexp.(q, e)) - sqrt(T(2))) <= 4 * epsT
            end
        end
        # negative control: a test double that keeps only d- reports about 2 for an antiparallel pair, so (a) and (f) discriminate it
        u = SVector{2, T}(cos(T(1)), sin(T(1)))
        @test alignment_indices(u, -u)[1] > 1.9 && sali_pair(u, -u) == 0
    end
end

@testset "AC2 normalisation interval invariance" begin
    _ = (SALI, sali_floor)
    # I1 (A1): the same start and frame with tau in {1, 2, 5, 10}, K in {0.5 (regular class), 2 pi}, n in {10, 100}, 20 starts per type
    # (D38: none reduced, n is short): |SALI_tau(n) - SALI_1(n)| <= B_SALI^tau(n) + B_SALI^1(n) (= 2 B_SALI, both against the same reference)
    # and |L_tau(n) - L_1(n)| <= B_L^tau(n) + B_L^1(n).
    m, N = 20, 100
    for T in SL_TYPES, label in (:half, :two_pi)
        base = sl_rows(label, T, 1, m, N)
        by1 = Dict((r.j, r.n) => r for r in base)
        for tau in (2, 5, 10)
            rows = [r for r in sl_rows(label, T, tau, m, N) if r.n in (10, 100)]
            @test length(rows) == 2m
            wS = maximum(abs(Float64(r.sali - by1[(r.j, r.n)].sali)) / (r.BS + by1[(r.j, r.n)].BS) for r in rows)
            wL = maximum(abs(Float64(r.L - by1[(r.j, r.n)].L)) / (r.BL + by1[(r.j, r.n)].BL) for r in rows)
            @test wS <= 1
            @test wL <= 1
        end
    end
end

# Closed form of A10 for the frame W0 = [[a1, a2], [b1, b2]] (columns (a_i, b_i)): BigFloat SALI and L at step n.
function sl_shear_cf(W0, n::Integer)
    return setprecision(BigFloat, SL_BITS) do
        a1, b1, a2, b2 = sl_big(W0[1, 1]), sl_big(W0[2, 1]), sl_big(W0[1, 2]), sl_big(W0[2, 2])
        det0 = abs(a1 * b2 - a2 * b1)
        w1, w2 = (a1 + n * b1, b1), (a2 + n * b2, b2)
        dotp = abs(w1[1] * w2[1] + w1[2] * w2[2])
        (S = 2 * sin(atan(det0, dotp) / 2), L = log(det0) - log(sqrt(w1[1]^2 + w1[2]^2)) - log(sqrt(w2[1]^2 + w2[2]^2)))
    end
end
# ln SALI_cf(j), j = 1..N, in Float64 (a scale for Lambda_n).
function sl_shear_lnsali(W0, N::Integer)
    a1, b1, a2, b2 = Float64(W0[1, 1]), Float64(W0[2, 1]), Float64(W0[1, 2]), Float64(W0[2, 2])
    det0 = abs(a1 * b2 - a2 * b1)
    return [log(2 * sin(atan(det0, abs((a1 + j * b1) * (a2 + j * b2) + b1 * b2)) / 2)) for j in 1:N]
end

@testset "AC3 K = 0 shear" begin
    _ = (SALI, sali_floor)
    # The tangent is [[1, 1], [0, 1]] everywhere (tangent-core E8), so the answer does not depend on the orbit. Seeded frame and W0 = I, tau in {1, 10},
    # n in {10, 100, 1e4, 1e6}: |SALI(n) - SALI_cf(n)| <= B_SALI(n) with the BigFloat closed form of A10; |L - L_cf| <= B_L; for W0 = I additionally
    # |SALI n - 1| <= 3/(8 n^2) + n B_SALI(n). Lambda_n from the closed form, Psi_n from the constant tangent. Reported: SALI n^2 |b1 b2| / |det W0|.
    # D38: n = 1e6 for Float64 only (Double64 sin/cos at |x| ~ 1e6 cost 10 us per step), 3 starts there; 8 starts at n = 1e4 for Double64.
    seed = sl_seed(3)
    rng = StableRNG(seed)
    M = SMatrix{2, 2, Float64}(1, 0, 1, 1)
    for T in SL_TYPES
        epsT = sl_eps(T)
        sys = StandardMap(exact(T, 0))
        for tau in (1, 10), (case, nid) in (("seeded", 20), ("identity", 1))
            for N in (10_000, 1_000_000)
                (N == 1_000_000 && T !== Float64) && continue
                m = case == "identity" ? 1 : (N == 1_000_000 ? 3 : (T === Double64 ? 8 : nid))
                x0 = [sl_cast(T, SVector(2pi * rand(rng), 2pi * rand(rng))) for _ in 1:m]
                W0 = case == "identity" ? [SMatrix{2, 2, T}(1, 0, 0, 1)] : [sl_frame(seed, j, T) for j in 1:m]
                ns = N == 1_000_000 ? [1_000_000] : [10, 100, 10_000]
                runs = sl_run(SALI(), sys, x0, W0, N, tau; steps = ns)
                for j in 1:m
                    lns = sl_shear_lnsali(W0[j], N)
                    Ms = fill(M, N)
                    for s in runs[j].snaps
                        n = s.n
                        cf = sl_shear_cf(W0[j], n)
                        p1, p2 = get!(SL_CACHE, (:shearpsi, T, case, j, n)) do      # Psi does not depend on tau
                            (lb_Psi(view(Ms, 1:n), W0[j][:, 1]), lb_Psi(view(Ms, 1:n), W0[j][:, 2]))
                        end
                        BS = sl_BSALI(epsT, sl_Lambda(lns, n), Float64(cf.S), p1, p2)
                        @test sl_err(s.sali, cf.S) <= BS
                        @test sl_err(s.L, cf.L) <= sl_BL(epsT, p1, p2)
                        if case == "identity"
                            @test abs(Float64(s.sali) * n - 1) <= 3 / (8 * n^2) + n * BS
                        end
                    end
                    if case == "seeded" && N == 10_000 && j == 1
                        s = runs[j].snaps[end]
                        b1b2 = abs(W0[j][2, 1] * W0[j][2, 2])
                        @info "AC3 reported: SALI n^2 |b1 b2| / |det W0| at n = 1e4" T tau value = Float64(s.sali * s.n^2 * b1b2 / abs(det(W0[j])))
                    end
                end
            end
        end
    end
end
