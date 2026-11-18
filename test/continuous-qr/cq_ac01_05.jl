# AC1 to AC5 and the algebra of the right-hand side: constant Hamiltonian matrices (no orbit), closed forms of the spec (C19) and 256-bit references (C16).
# Every testset starts by naming the new symbols it uses, so a missing module errors that testset (UndefVarError). The "controls" and "oracle" testsets use no new name:
# they run the check against test-side doubles with one defect each (and the closed forms against the 256-bit exponential) and pass on the code as it is today.

# ------------------------------------------------------------------ the algebra of the right-hand side (C5, C8, C10 to C13), pure functions at seeded random points
@testset "ALG right-hand side: on the manifold both forms agree, the penalty and the defect equation of C13, one column (C11)" begin
    _ = (cqr_rhs, cqr_growth)
    rng = StableRNG(cq_seed(0, 1)); d = 4
    for trial in 1:6
        A = randn(rng, d, d); As = SMatrix{d, d, Float64}(A)
        # k = d on the manifold: form 84 equals form 92 (C8 reduces to C5) and the penalty vanishes
        Q = cq_frame(rng, d, d, Float64)
        @test maximum(abs, cqr_rhs(Q, As, 0.0; form = :f84) - cqr_rhs(Q, As, 0.0; form = :f92)) <= 64 * eps(Float64) * norm(As)
        @test maximum(abs, cqr_rhs(Q, As, 3.0; form = :f84) - cqr_rhs(Q, As, 0.0; form = :f84)) <= 16 * eps(Float64) * 3 * norm(Q' * Q - I)  + 1.0e-14
        # the growth is the diagonal of Q^T A Q (C6)
        @test cqr_growth(Q, As) ≈ diag(Q' * As * Q) rtol = 1.0e-13
        # off the manifold (E of largest entry about 0.1): the defect equation of C13, E' = G^T Q + Q^T G
        Qo = Q * (I + 0.05 * (x -> (x + x') / 2)(randn(rng, d, d))); Qo = SMatrix{d, d, Float64}(Qo)
        E = Qo' * Qo - I; U = Qo' * As * Qo; L = [i > j ? U[i, j] : 0.0 for i in 1:d, j in 1:d]; S = L - L'
        for nu in (0.0, 2.5)
            G = cqr_rhs(Qo, As, nu; form = :f84)
            Edot = G' * Qo + Qo' * G
            @test maximum(abs, Edot - (E * S - S * E - 2nu * E * (I + E))) <= 1.0e-12 * (1 + norm(U))
        end
        # k = 2 < d, form 92: E' = -(E U + U^T E) - 2 nu E (I + E) with the k x k blocks (C13)
        Q2 = SMatrix{d, 2, Float64}(Qo[:, 1:2]); E2 = Q2' * Q2 - I; U2 = Q2' * As * Q2
        for nu in (0.0, 2.5)
            G = cqr_rhs(Q2, As, nu; form = :f92)
            @test maximum(abs, (G' * Q2 + Q2' * G) - (-(E2 * U2 + U2' * E2) - 2nu * E2 * (I + E2))) <= 1.0e-12 * (1 + norm(U2))
        end
        # H(Q) Q of Skokos2010 Eq. (92): (A - Q Q^T A + Q S Q^T) Q with the same S (C10) at a non-orthonormal Q2
        L2 = [i > j ? U2[i, j] : 0.0 for i in 1:2, j in 1:2]; S2 = L2 - L2'
        @test maximum(abs, cqr_rhs(Q2, As, 0.0; form = :f92) - (As - Q2 * Q2' * As + Q2 * S2 * Q2') * Q2) <= 1.0e-12 * (1 + norm(U2))
        # one column, on the manifold: q' = A q - (q . A q) q (C11)
        q1 = SMatrix{d, 1, Float64}(Q[:, 1])
        @test maximum(abs, cqr_rhs(q1, As, 0.0; form = :f92) - (As * q1 - dot(q1, As * q1) * q1)) <= 64 * eps(Float64) * norm(As)
    end
    # the form of C5 needs all columns: k < d with form 84 is refused (P1)
    @test_throws ArgumentError cqr_rhs(SMatrix{4, 2, Float64}(1, 0, 0, 0, 0, 1, 0, 0), SMatrix{4, 4, Float64}(I), 0.0; form = :f84)
end

# ------------------------------------------------------------------ AC1: fixed point (C19 (a), I4)
@testset "AC1 fixed point Q = I for upper triangular A: bitwise Q, rounding-only growth" begin
    _ = (ContinuousQR, propagate_cqr_orbit, propagate_cqr_series)
    n = cq_size(400, 2000); h = 0.05
    for (name, A) in (("diag2", cq_diag2()), ("diag4", cq_diag4()), ("tri2", cq_tri2()), ("tri4", cq_tri4())), T in (Float64, Double64), q in (2, 4)
        d = size(A, 1); cq = cq_est(A, T, h, q)
        _, Q, fin, ser = propagate_cqr_series(cq, cq_x0(d, T), cq_eye(d, d, T), n; every = 1)
        @test all(s -> s.Q == cq_eye(d, d, T), ser) && Q == cq_eye(d, d, T)                     # (a) tolerance 0
        worst = setprecision(BigFloat, CQ_BITS) do
            w = 0.0
            for s in ser, i in 1:d
                a = cq_big(T(A[i, i])); t = s.n * cq_big(T(h))
                w = max(w, Float64(abs(cq_big(s.l[i]) - a * t)) / (16 * cq_eps(T) * max(1.0, Float64(abs(a * t)))))
            end
            w
        end
        @test worst <= 1                                                                        # (b) 16 eps max(1, |A_ii| t)
        @test fin.Emax == 0
    end
end

@testset "AC1 negative controls (test-side doubles, no new name)" begin
    n = 400; h = 0.05; T = Float64
    run(A, mode, q) = (d = size(A, 1); cq_scheme(_ -> SMatrix{d, d, T}(A), cq_eye(d, d, T), T(h), q, n; mode = mode))
    I2 = cq_eye(2, 2, T); I4 = cq_eye(4, 4, T)
    for q in (2, 4)
        @test run(cq_diag2(), :ok, q).Q == I2 && run(cq_tri2(), :ok, q).Q == I2 && run(cq_tri4(), :ok, q).Q == I4      # the correct algebra leaves Q = I
        @test run(cq_diag2(), :upper, q).Q == I2 && run(cq_diag4(), :upper, q).Q == I4                                    # S from the upper triangle: fine for diagonal A
        @test run(cq_tri2(), :upper, q).Q != I2 && run(cq_tri4(), :upper, q).Q != I4                                      # but not for the triangular fixtures: fails (a)
        @test run(cq_diag2(), :nuQ, q).Q != I2                                                                            # penalty acting on Q instead of Q E: fails (a)
        r = run(cq_diag4(), :signu, q); @test abs(r.l[1] + 0.7 * n * h) <= 1.0e-12 && abs(r.l[1] - 0.7 * n * h) > 1                # wrong sign of u: fails (b)
    end
end

# ------------------------------------------------------------------ AC2: elliptic degenerate and skew (C19 (d), F7)
@testset "AC2 skew A (oscillator, all omega = 1): l_i = 0 to rounding, Q(t) = exp(A t) Q0" begin
    _ = (ContinuousQR, propagate_cqr_orbit, propagate_cqr_series)
    h = 0.05; tmax = cq_size(40.0, 100.0); n = round(Int, tmax / h)
    for N in (1, 2), q in (2, 4), (T, nf) in ((Float64, 3), (Double64, 1))
        d = 2N; A = cq_skew(N); cq = cq_est(A, T, h, q); rng = StableRNG(cq_seed(2, N))
        frames = vcat([cq_eye(d, d, T)], [cq_frame(rng, d, d, T) for _ in 1:nf])
        for Q0 in frames
            every = max(1, n ÷ (q ÷ 2) ÷ 100)
            _, Q, fin, ser = propagate_cqr_series(cq, cq_x0(d, T), Q0, n; every = every)
            @test all(s -> all(i -> abs(s.l[i]) <= 8 * cq_eps(T) * s.n * h, 1:d), ser)
            AF = cq_AF(A); worst = 0.0
            for s in ser
                t = s.n * T(h)
                Y = SMatrix{d, d, T}(exact_flow_matrix(HarmonicOscillator(SVector{N, T}(ntuple(_ -> one(T), N))), t))
                B = min(cq_Blin(Float64(t), h, q, AF, T), 0.25)                                 # B_lin exceeds the range of the entries at t = 100: capped (plan, calibration)
                worst = max(worst, Float64(maximum(abs, s.Q - Y * Q0)) / B)
            end
            @test worst <= 1
            @test fin.Emax <= CQ_EMAX
        end
    end
end

# ------------------------------------------------------------------ AC3, AC4: hyperbolic closed forms against the 256-bit variational reference (C16), the order of the scheme
cq_hyp_fixtures() = (("a", cq_inv_osc(0.8)), ("b", cq_coupled([0.8, 0.35], CQ_B)), ("c", cq_coupled([0.5, 0.5], CQ_B)))
function cq_hyp_error(A, ::Type{T}, Q0, h, q, t) where {T}
    d = size(A, 1); n = round(Int, t / h); cq = cq_est(A, T, h, q; orthogonality_tol = 0.1)      # h up to 0.1: the stage defect (Delta ||A||)^2 / 2 is above the default tolerance (plan, U6)
    _, _, fin = propagate_cqr_orbit(cq, cq_x0(d, T), Q0, n)
    return setprecision(BigFloat, CQ_BITS) do
        ref = cq_lref(cq.flow.A, Q0, n * cq_big(T(h)))
        (err = [Float64(abs(cq_big(fin.l[i]) - ref[i])) for i in 1:d], ref = Float64.(ref), fin = fin)
    end
end

@testset "AC3 hyperbolic closed forms: l_i against the 256-bit Gram-determinant reference, individually for the degenerate pair" begin
    _ = (ContinuousQR, propagate_cqr_orbit)
    rng = StableRNG(cq_seed(3, 1)); worst = 0.0
    for (name, A) in cq_hyp_fixtures()
        d = size(A, 1); AF = cq_AF(A)
        for (T, hs, ts) in ((Float64, (0.05, 0.025), (2.0, 10.0)), (Double64, (0.05,), (2.0,))), q in (2, 4)
            frames = vcat([cq_eye(d, d, T)], [cq_frame(rng, d, d, T) for _ in 1:(T === Float64 ? 3 : 1)])
            for Q0 in frames, h in hs, t in ts
                r = cq_hyp_error(A, T, Q0, h, q, t); B = cq_Blin(t, h, q, AF, T)
                @test all(r.err .<= B)
                @test r.fin.Emax <= CQ_EMAX
                worst = max(worst, maximum(r.err) / B)
            end
        end
    end
    @info "AC3 reported (R3): worst error / B_lin" worst
end

@testset "AC3 control: the individual values of the degenerate case depend on the frame; a double that sorts the exponents misses them (no new name)" begin
    # In the degenerate case (c) the individual l_1, l_2 depend on the frame; a double that sorts them (or returns them by size) fails on a frame with l_1 < l_2.
    A = cq_coupled([0.5, 0.5], CQ_B); rng = StableRNG(cq_seed(3, 2)); T = Float64; h = 0.025; t = 10.0; n = round(Int, t / h); B = cq_Blin(t, h, 4, cq_AF(A), T); nmiss = 0; best = 0.0
    for _ in 1:20
        Q0 = cq_frame(rng, 4, 4, T)
        ref = setprecision(BigFloat, CQ_BITS) do
            Float64.(cq_lref(SMatrix{4, 4, T}(A), Q0, n * h))
        end
        ok = cq_scheme(_ -> SMatrix{4, 4, T}(A), Q0, T(h), 4, n); res = cq_scheme(_ -> SMatrix{4, 4, T}(A), Q0, T(h), 4, n; mode = :sorted)
        @test abs(ok.l[1] - ref[1]) <= B && abs(ok.l[2] - ref[2]) <= B                          # the correct values pass individually
        best = max(best, ref[2] - ref[1])
        ref[2] - ref[1] > 2B && (nmiss += abs(res.l[1] - ref[1]) > B)
    end
    @test nmiss >= 1                                                                            # at least one seeded frame has l_1 < l_2 by more than 2 B and the sorted double fails there
    @info "AC3 control: largest l_2 - l_1 over 20 seeded frames, frames where the sorted double fails" best bound = B nmiss
end

@testset "AC3 oracle: the 256-bit exponential and Gram determinants reproduce the closed forms of C19 (no new name)" begin
    setprecision(BigFloat, CQ_BITS) do
        lam = BigFloat(8) / 10; t = BigFloat(3)
        Y = cq_expm(Matrix{BigFloat}([0 1; lam^2 0]) * t)
        @test maximum(abs, Y - [cosh(lam * t) sinh(lam * t) / lam; lam * sinh(lam * t) cosh(lam * t)]) <= BigFloat(2)^-230
        # exact for the exact matrix: the stored Float64 matrix differs from the decimal one by rounding only; the Gram-determinant l_1 equals ln ||Y q|| for one column
        Q0 = Matrix{BigFloat}([0.6 -0.8; 0.8 0.6])
        l = cq_lref(Matrix{BigFloat}([0 1; lam^2 0]), Q0, t)
        @test abs(l[1] - log(norm(Y * Q0[:, 1]) / norm(Q0[:, 1]))) <= BigFloat(2)^-230        # Q0 is the Float64 rotation, orthonormal only to rounding: l is relative to the frame itself
        @test abs(l[1] + l[2] - 0) <= BigFloat(2)^-230                                          # det Y = 1 (trace 0)
    end
end

@testset "AC4 order of the Q-scheme: error ratios of consecutive h lie in [2^q / 2, 2 x 2^q] above the rounding floor" begin
    _ = (ContinuousQR, propagate_cqr_orbit)
    t = 10.0; hs = (0.1, 0.05, 0.025, 0.0125); T = Float64; rng = StableRNG(cq_seed(4, 1)); ratios = Dict()
    for (name, A) in cq_hyp_fixtures()[1:2], q in (2, 4)
        d = size(A, 1); Q0 = cq_frame(rng, d, d, T)
        errs = [maximum(cq_hyp_error(A, T, Q0, h, q, t).err) for h in hs]
        floor_ = 1.0e3 * cq_eps(T) * t
        for j in 1:3
            (errs[j] > floor_ && errs[j + 1] > floor_) || continue
            @test 2^q / 2 <= errs[j] / errs[j + 1] <= 2 * 2^q
        end
        ratios[(name, q)] = (errs, errs[1:3] ./ errs[2:4])
    end
    @info "AC4 reported: absolute errors and ratios" ratios
end

# ------------------------------------------------------------------ AC5: rotated frame (C19 (b))
@testset "AC5 rotated frame, A = diag(0.8, -0.8): l_1, the sum and Q(t) against the closed form" begin
    _ = (ContinuousQR, propagate_cqr_orbit, propagate_cqr_series)
    lam = 0.8; A = [lam 0.0; 0.0 -lam]; h = 0.05; tmax = cq_size(10.0, 20.0); n = round(Int, tmax / h); T = Float64; AF = cq_AF(A); worst = 0.0
    for theta in (0.3, 1.0, 1.4), q in (2, 4)
        cq = cq_est(A, T, h, q); Q0 = cq_rot2(theta, T)
        _, _, fin, ser = propagate_cqr_series(cq, cq_x0(2, T), Q0, n; every = max(1, n ÷ (q ÷ 2) ÷ 40))
        for s in ser
            t = s.n * h; B = min(cq_Blin(t, h, q, AF, T), 0.25)
            e1 = abs(s.l[1] - cq_l1_rot(lam, theta, t)); e2 = abs(s.l[1] + s.l[2])
            eQ = maximum(abs, s.Q - cq_rot2(cq_phi_rot(lam, theta, t), T))
            @test e1 <= B && e2 <= B && eQ <= B / max(AF * t, 1)
            worst = max(worst, e1 / B, e2 / B, eQ / (B / max(AF * t, 1)))
        end
    end
    @info "AC5 reported (R3): worst error / bound" worst
end

@testset "AC5 negative control: S of the wrong sign fails l_1 and Q (the sum is blind to it, U3; test-side double, no new name)" begin
    lam = 0.8; A = SMatrix{2, 2, Float64}(lam, 0, 0, -lam); h = 0.05; t = 10.0; n = round(Int, t / h); theta = 1.0; AF = cq_AF(A)
    res = cq_scheme(_ -> A, cq_rot2(theta, Float64), h, 4, n; mode = :Sflip); B = cq_Blin(t, h, 4, AF, Float64)
    @test abs(res.l[1] - cq_l1_rot(lam, theta, t)) > B
    @test abs(res.l[1] + res.l[2]) <= B        # spec point U3: tr(Q^T A Q) = tr A for any orthogonal Q, so the sum cannot detect a wrong S: only l_1 and Q fail
    @test maximum(abs, res.Q - cq_rot2(cq_phi_rot(lam, theta, t), Float64)) > B / (AF * t)
    ok = cq_scheme(_ -> A, cq_rot2(theta, Float64), h, 4, n)
    @test abs(ok.l[1] - cq_l1_rot(lam, theta, t)) <= B && maximum(abs, ok.Q - cq_rot2(cq_phi_rot(lam, theta, t), Float64)) <= B / (AF * t)   # the test-side algebra itself is right
end
