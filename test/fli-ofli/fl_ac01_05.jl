# AC1 to AC5: algebra of the definition, K = 0 shear, linear hyperbolic map, hyperbolic fixed point, rotation and elliptic cases.
# Each testset starts by naming the Indicators symbols it uses, so a missing module errors that testset (UndefVarError) and nothing else.
# Tolerances: B_S(n) = 8 eps (n + G_first) of L1 S3 exactly, except for the two supplied frames of AC3 with a tiny expanding coefficient (|c_1| = 1e-3, 1e-6),
# which use max(B_S(n), 8 eps Psi_n^(i)), the per-vector form of L1 Amendment D35 (the plain B_S is exceeded there, L1 P6; test plan, Deviations item 2).

# ------------------------------------------------------------------ AC1

# Exact partial sums P(m) of c_m (the stored Float64 values), and the exact running maxima of the definition for the diagonal test maps: S_1 = P, S_2 = -P.
function fl_alg_exact(kind::Symbol, N::Int)
    return setprecision(BigFloat, FL_BITS) do
        P = Vector{BigFloat}(undef, N)
        F1, FB = similar(P), similar(P)
        acc = BigFloat(0)
        r1 = rB = BigFloat(0)                       # the s = 0 term
        for m in 1:N
            acc += BigFloat(fl_c(kind, m))
            r1 = max(r1, acc)
            rB = max(rB, acc, -acc)
            P[m], F1[m], FB[m] = acc, r1, rB
        end
        (P = P, F1 = F1, FB = FB)
    end
end
# The criterion as a predicate, so that the negative controls can be run through it: both series within B of the exact running maxima.
fl_alg_ok(f1, fB, ex, B) = all(j -> fl_err(f1[j], ex.F1[j]) <= B[j] && fl_err(fB[j], ex.FB[j]) <= B[j], eachindex(f1))

@testset "AC1 algebra of the definition" begin
    _ = (FLI,)
    # Prescribed diagonal tangent diag(a_m, 1/a_m), a_m = exp(c_m), supplied identity frame, tau = 1, 40 handoffs: (i) c_m = -0.3 for m <= 10 then +0.5; (ii) c_m = -0.2;
    # (iii) c_m alternating +0.2, -0.2. Pass: FLI_1(j) = running maximum of the reducer's own S_1 including the 0 at j = 0, bitwise; FLI_B(j) = running maximum of
    # max(S_1, S_2), bitwise; both within B_S(j) of the exact running maxima of the partial sums of c_m; nondecreasing; exchange of the columns leaves FLI_B bitwise
    # unchanged and gives FLI_1 = the running maximum of the other column's sum. Negative controls: no memory, no j = 0 term, no log.
    for T in FL_TYPES, kind in (:a, :b, :c)
        epsT = fl_eps(T)
        N = 40
        ex = fl_alg_exact(kind, N)
        B = [8 * epsT * (j + 2 + exp(2 * abs(fl_c(kind, 1)))) for j in 1:N]       # 8 eps (j + G_first), G_first = a^2 + a^-2 <= 2 + e^(2|c|)
        sys = fl_Diag(kind)
        x0 = SVector{2, T}(0, 0)
        I2, P2 = SMatrix{2, 2, T}(1, 0, 0, 1), SMatrix{2, 2, T}(0, 1, 1, 0)
        snaps(W) = fl_run((FLI(),), sys, [x0], [W], N, 1; all = true)[1][1].snaps
        sn = snaps(I2)
        @test length(sn) == N
        f1, fB = [s.fli1 for s in sn], [s.fli for s in sn]
        S1, S2 = [s.sum[1] + s.comp[1] for s in sn], [s.sum[2] + s.comp[2] for s in sn]
        @test isequal(f1, accumulate(max, S1; init = zero(T)))                           # own S_1, with the 0 at j = 0
        @test isequal(fB, accumulate(max, max.(S1, S2); init = zero(T)))
        @test issorted(f1) && issorted(fB)
        @test fl_alg_ok(f1, fB, ex, B)
        # the helper against the values the spec prints
        if kind === :a
            @test abs(ex.FB[10] - 3) < 1e-14 && abs(ex.FB[22] - 3) < 1e-14 && ex.F1[16] < 1e-14 && abs(ex.F1[17] - 0.5) < 1e-14
            @test all(j -> f1[j] <= B[j], 1:15) && all(>=(0), f1)                        # FLI_1 = 0 until the sum turns positive (j = 17; at j = 16 the stored c_m sum to +1.1e-16)
            @test fl_err(f1[17], ex.F1[17]) <= B[17] && f1[17] > 0.4
        elseif kind === :b
            @test all(j -> 0 <= f1[j] <= B[j], 1:N)                                      # FLI_1 = 0 for all j
            @test all(j -> fl_err(fB[j], BigFloat(0.2) * j) <= B[j], 1:N)                # FLI_B = 0.2 j
        end
        # exchange of the columns (frame [e2 e1])
        sw = snaps(P2)
        @test isequal([s.fli for s in sw], fB)
        f1sw = [s.fli1 for s in sw]
        @test isequal(f1sw, accumulate(max, S2; init = zero(T)))
        kind === :b && @test f1sw[end] > 0.1 && f1[end] <= B[end]                        # FLI_1 is not symmetric
        # sign change of a column (frame [-e1 e2]): both forms bitwise unchanged
        sg = snaps(SMatrix{2, 2, T}(-1, 0, 0, 1))
        @test isequal([s.fli for s in sg], fB) && isequal([s.fli1 for s in sg], f1)
        # negative controls: each double fails the predicate of the criterion
        if kind !== :b
            @test !fl_alg_ok(S1, max.(S1, S2), ex, B)                                    # the current S: no memory (fails (i) and (iii))
        else
            @test !fl_alg_ok(accumulate(max, S1), fB, ex, B)                             # no j = 0 term: FLI_1 = -0.2 at j = 1 instead of 0
        end
        @test !fl_alg_ok(f1, exp.(fB), ex, B)                                            # the unlogged maximum of the norms
    end
end

# ------------------------------------------------------------------ AC2

# Squared norms of the shear columns after n steps: column (a, b) -> (a + n b, b) (L1 E11), BigFloat from the exact widened frame.
function fl_shear_q(W0, n::Integer)
    return setprecision(BigFloat, FL_BITS) do
        a1, b1, a2, b2 = fl_big(W0[1, 1]), fl_big(W0[2, 1]), fl_big(W0[1, 2]), fl_big(W0[2, 2])
        ((a1 + n * b1)^2 + b1^2, (a2 + n * b2)^2 + b2^2)
    end
end
# FLI closed forms of Q9: the squared norm is convex in n, so the maximum over the handoffs 0, tau, ..., n is at an end point.
function fl_cf(q1, q2)
    return setprecision(BigFloat, FL_BITS) do
        one_ = one(BigFloat)
        (F1 = log(max(one_, q1)) / 2, FB = log(max(one_, q1, q2)) / 2)
    end
end

@testset "AC2 K = 0 shear" begin
    _ = (FLI,)
    # The tangent is [[1, 1], [0, 1]] everywhere (tangent-core E8): the answer does not depend on the orbit. Seeded frame and W0 = I, tau in {1, 10},
    # n in {10, 100, 1e4, 1e6}: |FLI - FLI_cf| <= B_S(n) = 8 eps (n + G_first) for both forms, FLI_cf = ln max(1, norm w(n)) (Q9a); W0 = I: FLI_1 = 0 exactly
    # (tolerance 0) and |FLI_B - (1/2) ln(n^2 + 1)| <= B_S(n). D38: n = 1e6 for Float64 only (3 seeded starts, 1 identity; Double64 sin/cos at |x| ~ 1e6 are slow).
    seed = fl_seed(2)
    rng = StableRNG(seed)
    M = SMatrix{2, 2, Float64}(1, 0, 1, 1)
    # the iterated BigFloat reference of the helper against the closed form (the helper is the oracle of the other criteria)
    let W = fl_frame(seed, 1, Float64), R = fl_ref_steps(fl_mats_const(1, 1, 0, 1), W, 100), cf = fl_shear_q(W, 100)
        @test abs(R.q1[100] - cf[1]) < 1e-60 * cf[1] && abs(R.q2[100] - cf[2]) < 1e-60 * cf[2]
    end
    for T in FL_TYPES
        epsT = fl_eps(T)
        sys = StandardMap(exact(T, 0))
        for tau in (1, 10), (case, nid) in (("seeded", 20), ("identity", 1))
            for N in (10_000, 1_000_000)
                (N == 1_000_000 && T !== Float64) && continue
                m = case == "identity" ? 1 : (N == 1_000_000 ? 3 : nid)
                x0 = [fl_cast(T, SVector(2pi * rand(rng), 2pi * rand(rng))) for _ in 1:m]
                W0 = case == "identity" ? [SMatrix{2, 2, T}(1, 0, 0, 1)] : [fl_frame(seed, j, T) for j in 1:m]
                ns = N == 1_000_000 ? [1_000_000] : [10, 100, 10_000]
                runs = fl_run((FLI(),), sys, x0, W0, N, tau; steps = ns)
                B(n) = 8 * epsT * (n + fl_Gfirst(M, tau))
                for j in 1:m
                    snaps = runs[j][1].snaps
                    @test [s.n for s in snaps] == ns
                    for s in snaps
                        n = s.n
                        cf = fl_cf(fl_shear_q(W0[j], n)...)
                        @test fl_err(s.fli1, cf.F1) <= B(n)
                        @test fl_err(s.fli, cf.FB) <= B(n)
                        if case == "identity"
                            @test s.fli1 == 0                                                       # tolerance 0: S_1 = ln 1
                            @test fl_err(s.fli, setprecision(() -> log(BigFloat(n)^2 + 1) / 2, BigFloat, FL_BITS)) <= B(n)
                        end
                    end
                end
            end
        end
    end
end

# ------------------------------------------------------------------ AC3 and AC4: constant tangents (Y5 and the fixed point of the standard map)

# B_S(n) = 8 eps (n + G_first) of L1 S3 for the constant tangent M (Float64 copy) at step n; with `psi`, max(B_S(n), 8 eps Psi_n^(i)) (Psi memo per frame).
function fl_lin_B(epsT, M::SMatrix, Ms, W0, n::Int, tau::Int, memo, psi::Bool = false)
    psi || return 8 * epsT * (n + fl_Gfirst(M, tau))
    p1, p2 = get!(memo, n) do
        (lb_Psi(view(Ms, 1:n), W0[:, 1]), lb_Psi(view(Ms, 1:n), W0[:, 2]))
    end
    return 8 * epsT * max(n + fl_Gfirst(M, tau), p1, p2)
end

# Eigen data of A = [[2, 1], [1, 1]] in BigFloat (call inside setprecision): mu and the orthonormal v_+, v_- = R90 v_+ (L1 E10).
function fl_y5_eigen()
    mu = (3 + sqrt(BigFloat(5))) / 2
    phi = (1 + sqrt(BigFloat(5))) / 2
    nv = sqrt(phi^2 + 1)
    return (mu = mu, vp = (phi / nv, 1 / nv), vm = (-1 / nv, phi / nv))
end
# Coefficients (c_1, c_2) of the two columns of W0 on (v_+, v_-): c1 = expanding, c2 = contracting (BigFloat, from the exact widened frame).
function fl_y5_coeffs(W0)
    return setprecision(BigFloat, FL_BITS) do
        e = fl_y5_eigen()
        w = ((fl_big(W0[1, 1]), fl_big(W0[2, 1])), (fl_big(W0[1, 2]), fl_big(W0[2, 2])))
        c1 = ntuple(i -> w[i][1] * e.vp[1] + w[i][2] * e.vp[2], 2)
        c2 = ntuple(i -> w[i][1] * e.vm[1] + w[i][2] * e.vm[2], 2)
        (c1 = c1, c2 = c2)
    end
end
# Supplied unit frame whose first column has expanding coefficient c1 (and contracting sqrt(1 - c1^2)); the second column is the first rotated by 90 degrees,
# so its expanding coefficient is -sqrt(1 - c1^2). Built in BigFloat, rounded into T.
function fl_y5_frame(c1::Real, ::Type{T}) where {T}
    return setprecision(BigFloat, FL_BITS) do
        e = fl_y5_eigen()
        c = BigFloat(c1)
        s = sqrt(1 - c^2)
        w1 = (c * e.vp[1] + s * e.vm[1], c * e.vp[2] + s * e.vm[2])
        w2 = (-w1[2], w1[1])
        SMatrix{2, 2, T}(T(w1[1]), T(w1[2]), T(w2[1]), T(w2[2]))
    end
end

@testset "AC3 linear hyperbolic map" begin
    _ = (FLI,)
    # A = [[2, 1], [1, 1]] (Y5), start at the origin (a fixed point), seeded ids 1..100 (ids with |c_1| < 0.05 for column 1 excluded, at most 10) and two supplied unit frames with |c_1| = 1e-3 and 1e-6 for column 1, tau in {1, 5}. (a) n in {1, 2, 4, 8, 20, 100, 200} (tau = 5: the handoffs
    # {5, 10, 20, 100, 200}): both forms within B of the BigFloat closed form; (b) |c_1| = 1e-3, n = 4, tau = 1: FLI_1 = 0 exactly while FLI_B > 0.5 (the initial vector
    # issue, F3); (c) rate at n in {20, 100}: |FLI_B - (n ln mu + ln |c_max|)| <= B + mu^(-4n) (1 + 1 / c_max^2), c_max the larger expanding coefficient of the
    # two columns (the explicit constant of the neglected term).
    seed = fl_seed(3)
    keep = [1, 2, 4, 5, 8, 10, 20, 100, 200]
    for T in FL_TYPES
        epsT = fl_eps(T)
        nid = 100
        sys = LinearSymplecticMap(SMatrix{2, 2, T}(2, 1, 1, 1))
        W0 = vcat([fl_frame(seed, j, T) for j in 1:nid], [fl_y5_frame(1e-3, T), fl_y5_frame(1e-6, T)])
        nfr = length(W0)
        x0 = [SVector{2, T}(0, 0) for _ in 1:nfr]
        M64 = SMatrix{2, 2, Float64}(2, 1, 1, 1)
        Ms = lb_tangents(sys, x0[1], 200)                       # the tangent is A everywhere
        refs = [fl_ref_steps(fl_mats_const(2, 1, 1, 1), W0[j], 200) for j in 1:nfr]
        memo = [Dict{Int, Tuple{Float64, Float64}}() for _ in 1:nfr]
        cs = [fl_y5_coeffs(W0[j]) for j in 1:nfr]
        excluded = [j for j in 1:nid if abs(cs[j].c1[1]) < 0.05]
        @test length(excluded) <= 10
        # supplied frames as intended: |c_1| of column 1, and the second column carries the expanding part
        @test abs(abs(cs[nid + 1].c1[1]) - 1e-3) < 1e-12 && abs(abs(cs[nid + 2].c1[1]) - 1e-6) < 1e-12 && abs(cs[nid + 1].c1[2]) > 0.999
        # the iterated reference against the matrix power A^8 (independent path), first 3 frames
        setprecision(BigFloat, FL_BITS) do
            Y8 = [BigFloat(2) 1; 1 1]^8
            for j in 1:3
                q = fl_qcols(Y8, W0[j])
                @test abs(q[1] - refs[j].q1[8]) < 1e-60 * q[1] && abs(q[2] - refs[j].q2[8]) < 1e-60 * q[2]
            end
        end
        for tau in (1, 5)
            ns = tau == 1 ? [1, 2, 4, 8, 20, 100, 200] : [5, 10, 20, 100, 200]
            run = fl_run((FLI(),), sys, x0, W0, 200, tau; steps = ns)
            for j in 1:nfr
                j in excluded && continue
                at = fl_ref_at(refs[j], tau, ns)
                snaps = run[j][1].snaps
                @test [s.n for s in snaps] == ns
                for (s, r) in zip(snaps, at)
                    n = s.n
                    B = fl_lin_B(epsT, M64, Ms, W0[j], n, tau, memo[j], j > nid)                # the supplied frames: the Psi form
                    @test fl_err(s.fli1, r.F1) <= B                                                 # (a)
                    @test fl_err(s.fli, r.FB) <= B
                    if n in (20, 100)                                                              # (c)
                        asym, tail = setprecision(BigFloat, FL_BITS) do
                            mu = fl_y5_eigen().mu
                            cmax = max(abs(cs[j].c1[1]), abs(cs[j].c1[2]))
                            (n * log(mu) + log(cmax), Float64(mu^(-4n) * (1 + 1 / cmax^2)))
                        end
                        @test fl_err(s.fli, asym) <= B + tail
                    end
                    if j == nid + 1 && tau == 1 && n == 4                                           # (b)
                        @test r.F1 == 0                                                            # the intended case: the exact FLI_1 is the s = 0 term
                        @test s.fli1 == 0
                        @test s.fli > 0.5
                    end
                end
            end
        end
    end
end

@testset "AC4 hyperbolic fixed point" begin
    _ = (FLI,)
    # Fixed point (0, 0) of the standard map, K in {1, 2, 6}, seeded ids 1..100, tau in {1, 5}, n as in AC3 (L1 E12, Q9c): both forms within B
    # of the BigFloat product norm, FLI_cf = ln max(1, norm w_1(n), norm w_2(n)).
    seed = fl_seed(4)
    for T in FL_TYPES, K in (1, 2, 6)
        epsT = fl_eps(T)
        nid = 100
        sys = StandardMap(exact(T, K))
        W0 = [fl_frame(seed, j, T) for j in 1:nid]
        x0 = [SVector{2, T}(0, 0) for _ in 1:nid]
        M64 = SMatrix{2, 2, Float64}(1 + K, K, 1, 1)
        Ms = lb_tangents(sys, x0[1], 200)
        refs = [fl_ref_steps(fl_mats_const(1 + K, 1, K, 1), W0[j], 200) for j in 1:nid]
        memo = [Dict{Int, Tuple{Float64, Float64}}() for _ in 1:nid]
        setprecision(BigFloat, FL_BITS) do
            Y8 = [BigFloat(1 + K) 1; BigFloat(K) 1]^8
            for j in 1:3
                q = fl_qcols(Y8, W0[j])
                @test abs(q[1] - refs[j].q1[8]) < 1e-60 * q[1] && abs(q[2] - refs[j].q2[8]) < 1e-60 * q[2]
            end
        end
        for tau in (1, 5)
            ns = tau == 1 ? [1, 2, 4, 8, 20, 100, 200] : [5, 10, 20, 100, 200]
            run = fl_run((FLI(),), sys, x0, W0, 200, tau; steps = ns)
            for j in 1:nid
                at = fl_ref_at(refs[j], tau, ns)
                for (s, r) in zip(run[j][1].snaps, at)
                    B = fl_lin_B(epsT, M64, Ms, W0[j], s.n, tau, memo[j])
                    @test fl_err(s.fli1, r.F1) <= B
                    @test fl_err(s.fli, r.FB) <= B
                end
            end
        end
    end
end

# ------------------------------------------------------------------ AC5

@testset "AC5 rotation and elliptic cases" begin
    _ = (FLI,)
    # (a) the order-4 rotation [[0, 1], [-1, 0]], supplied identity frame (as L1 criterion 2(a)), n up to 1e4, tau in {1, 3}: FLI = 0 exactly, both forms, at every
    # handoff (Q9d). (b) K = -2 at (0, 0) (Q9e, L1 E13), seeded ids 1..100, tau = 1: for every n >= 1, |FLI_1(n) - max(0, ln norm(M w_1(0)))| <= B
    # and FLI_B <= ln((3 + sqrt 5)/2) + B (and, derived the same way, FLI_B within B of max(0, ln norm(M w_1), ln norm(M w_2))); n in 1..100 and every 10th to 1000
    # (D38; FL_FULL: every n). (c) the same at tau = 2: FLI_1 and FLI_B are 0 within B at every handoff (the odd steps are never sampled).
    seed = fl_seed(5)
    for T in FL_TYPES
        epsT = fl_eps(T)
        rot = LinearSymplecticMap(SMatrix{2, 2, T}(0, -1, 1, 0))
        for tau in (1, 3)
            N = tau == 1 ? 10_000 : 9_999
            snaps = fl_run((FLI(),), rot, [SVector{2, T}(0, 0)], [SMatrix{2, 2, T}(1, 0, 0, 1)], N, tau; all = true)[1][1].snaps
            @test length(snaps) == N ÷ tau
            @test all(s -> s.fli == 0 && s.fli1 == 0, snaps)
        end
        sys = StandardMap(exact(T, -2))
        nid = 100
        x0 = [SVector{2, T}(0, 0) for _ in 1:nid]
        W0 = [fl_frame(seed, j, T) for j in 1:nid]
        nsamp = FL_FULL ? collect(1:1000) : sort(unique(vcat(1:100, 110:10:1000)))
        M64 = SMatrix{2, 2, Float64}(-1, -2, 1, 1)
        Ms = lb_tangents(sys, x0[1], 1000)
        refs = [fl_ref_steps(fl_mats_const(-1, 1, -2, 1), W0[j], 1000) for j in 1:nid]
        memo = [Dict{Int, Tuple{Float64, Float64}}() for _ in 1:nid]
        ln0 = setprecision(() -> log((3 + sqrt(BigFloat(5))) / 2), BigFloat, FL_BITS)
        for tau in (1, 2)
            ns = filter(n -> n % tau == 0, nsamp)
            run = fl_run((FLI(),), sys, x0, W0, 1000, tau; steps = ns)
            for j in 1:nid
                at = fl_ref_at(refs[j], tau, ns)
                # S_i(n) alternates 0 (n even) and ln norm(M w_i(0)) (n odd), so for tau = 1 and n >= 1 the maxima are these constants (Q9e)
                cf = setprecision(BigFloat, FL_BITS) do
                    l1, l2 = log(refs[j].q1[1]) / 2, log(refs[j].q2[1]) / 2
                    e1, e2 = log(sum(abs2, (fl_big(W0[j][1, 1]), fl_big(W0[j][2, 1])))) / 2, log(sum(abs2, (fl_big(W0[j][1, 2]), fl_big(W0[j][2, 2])))) / 2   # even n: ln norm(w_i(0)), O(eps) for a frame unit to roundoff
                    (F1 = max(zero(l1), l1, e1), FB = max(zero(l1), l1, l2, e1, e2))
                end
                tau == 1 && @test abs(at[end].F1 - cf.F1) < 1e-70 && abs(at[end].FB - cf.FB) < 1e-70          # the helper against the closed form
                for (s, r) in zip(run[j][1].snaps, at)
                    B = fl_lin_B(epsT, M64, Ms, W0[j], s.n, tau, memo[j])
                    if tau == 1
                        @test fl_err(s.fli1, cf.F1) <= B
                        @test fl_err(s.fli, cf.FB) <= B
                        @test fl_diff(s.fli, ln0) <= B                                                            # FLI_B <= ln((3 + sqrt 5)/2) + B
                    else
                        @test 0 <= s.fli1 <= B && 0 <= s.fli <= B                                                 # (c)
                    end
                end
            end
        end
    end
end
