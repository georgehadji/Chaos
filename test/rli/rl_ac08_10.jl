# AC8 to AC10: cross-check with the single-orbit stream and L1 (RL5, I5, F8), symmetry and column independence (I3, I4), precision ladder and representation.
# Every testset starts with a line naming the new symbols it uses.

@testset "AC8 cross-check with the single-orbit stream and L1" begin
    _ = (RLI, PairMap)
    # RL5, I5, F8. K in {0.5 (regular class), 2 pi, 10}, Float64 wrapped, tau = 1, 20 starts, n up to 1e4 (the sizes of the spec). (a) the pair orbit points at every handoff are bitwise
    # equal to the points of single-orbit runs started at x*_A and at x*_B, for the chains {empty}, {L1 on k = 1} (single) and {RLI} (pair); (b) the sums s_A, s_B of the RLI
    # reducer are bitwise equal to the L1 k = 1 sums (sum and compensation) of the two single-orbit runs with xi_0 as the supplied frame; (c) post-processing in BigFloat of the
    # two stored L1 series s_i(k) = k X_i(k): DL_post,k = |s_B(k) - s_A(k)| / k and RLI_post its mean: |DL_k - DL_post,k| <= 16 eps max(1, |s_A(k)|, |s_B(k)|) / k + 4 eps DL_post,k
    # at every handoff multiple of 10 and |RLI - RLI_post| <= 16 eps max(1, max_{j<=n,i} |s_i(j)|) H_K / n + 4 eps RLI_post (16 eps because the stored X is rounded, megno
    # Amendment D74); (d) the chain {RLI (non-replacing), MEGNO}: both outputs bitwise equal to their stand-alone values, MEGNO column 1 equal to column 1 of the single-orbit
    # MEGNO of orbit A and column 2 equal to column 1 of the single-orbit MEGNO at x*_B (xi_0 supplied as column 1; the non-replacing position changes nothing).
    T = Float64
    epsT = rl_eps(T)
    N = 10_000
    cps = collect(10:10:N)
    for label in (:half, :two_pi, :ten)
        m = 20
        c = rl_case(label, T, m)
        sys = c.sys
        P = PairMap(sys)
        prs = rl_pairs(c, 1e-10)
        flds = (:Y, :Ybar, :Ybar20, :sum, :comp, :R1, :R2, :Q, :Q2)
        worstD = worstR = 0.0
        for j in 1:m
            pr, xi = prs[j], c.xi0[j]
            W2, W1 = rl_frame(xi), rl_col(xi)
            # (a) the chains without a normalising reducer overflow after about 600 steps on a chaotic orbit (the deviation matrix is never rescaled), so they run N0 = 300 steps and are
            # compared with the first N0 handoffs of the normalising chains
            N0 = 300
            xp0, _, sp0 = propagate_orbit(P, pr.z, W2, N0, 1, (rl_Trace(),); orbit = j)
            xp1, _, sp1 = propagate_orbit(P, pr.z, W2, N, 1, (RLI(), rl_Trace()); orbit = j)
            xa0, _, sa0 = propagate_orbit(sys, pr.zA, W1, N0, 1, (rl_Trace(),); orbit = j)
            xa1, _, sa1 = propagate_orbit(sys, pr.zA, W1, N, 1, (Benettin(), rl_Trace()); orbit = j)
            xb0, _, sb0 = propagate_orbit(sys, pr.zB, W1, N0, 1, (rl_Trace(),); orbit = j)
            xb1, _, sb1 = propagate_orbit(sys, pr.zB, W1, N, 1, (Benettin(), rl_Trace()); orbit = j)
            tp0, tp1 = reducer_final(rl_Trace(), sp0[1]).xs, reducer_final(rl_Trace(), sp1[2]).xs
            ta = (reducer_final(rl_Trace(), sa0[1]).xs, reducer_final(rl_Trace(), sa1[2]).xs)
            tb = (reducer_final(rl_Trace(), sb0[1]).xs, reducer_final(rl_Trace(), sb1[2]).xs)
            @test length(tp0) == N0 && length(tp1) == N && length(ta[2]) == N
            @test isequal(tp0, tp1[1:N0]) && isequal(ta[1], ta[2][1:N0]) && isequal(tb[1], tb[2][1:N0])
            @test isequal(rl_half.(tp1, 1), ta[2]) && isequal(rl_half.(tp1, 2), tb[2]) && isequal(rl_half.(tp0, 1), ta[1]) && isequal(rl_half.(tp0, 2), tb[1])
            @test rl_half(xp0, 1) === xa0 && rl_half(xp0, 2) === xb0 && rl_half(xp1, 1) === xa1 && rl_half(xp1, 2) === xb1
            # (b)
            f = reducer_final(RLI(), sp1[1])
            lA, lB = reducer_final(Benettin(), sa1[1]), reducer_final(Benettin(), sb1[1])
            @test f.n == N && f.handoffs == N
            @test f.sum[1] == lA.sum[1] && f.comp[1] == lA.comp[1] && f.sum[2] == lB.sum[1] && f.comp[2] == lB.comp[1]
            # (c)
            rl = rl_run((RLI(),), P, [pr.z], [W2], N, 1; steps = cps)[1][1].snaps
            XA = rl_run((Benettin(),), sys, [pr.zA], [W1], N, 1; all = true)[1][1].snaps
            XB = rl_run((Benettin(),), sys, [pr.zB], [W1], N, 1; all = true)[1][1].snaps
            sLA, sLB = [k * XA[k].lambda[1] for k in 1:N], [k * XB[k].lambda[1] for k in 1:N]
            smax = accumulate(max, max.(abs.(sLA), abs.(sLB)))
            post = rl_ref_rli(BigFloat.(sLA), BigFloat.(sLB), 1)
            for f in rl
                n = f.n
                tolD = 16 * epsT * max(1.0, abs(sLA[n]), abs(sLB[n])) / n + 4 * epsT * Float64(post.DL[n])
                tolR = 16 * epsT * max(1.0, smax[n]) * rl_H(n) / n + 4 * epsT * Float64(post.RLI[n])
                eD, eR = rl_err(f.DL, post.DL[n]), rl_err(f.RLI, post.RLI[n])
                @test eD <= tolD && eR <= tolR
                worstD, worstR = max(worstD, eD / tolD), max(worstR, eR / tolR)
            end
            # (d)
            for tau in (1, 4)
                _, _, s_r = propagate_orbit(P, pr.z, W2, 400, tau, (RLI(),); orbit = j)
                _, _, s_m = propagate_orbit(P, pr.z, W2, 400, tau, (MEGNO(),); orbit = j)
                _, _, s_c = propagate_orbit(P, pr.z, W2, 400, tau, (RLI(replace = false), MEGNO()); orbit = j)
                @test isequal(reducer_final(RLI(), s_c[1]), reducer_final(RLI(), s_r[1])) && isequal(reducer_final(MEGNO(), s_c[2]), reducer_final(MEGNO(), s_m[1]))
                xi2 = rl_xi2(c.fseed, c.ids[j], T)
                Wm = SMatrix{2, 2, T}(xi[1], xi[2], xi2[1], xi2[2])
                _, _, s_a = propagate_orbit(sys, pr.zA, Wm, 400, tau, (MEGNO(),); orbit = j)
                _, _, s_b = propagate_orbit(sys, pr.zB, Wm, 400, tau, (MEGNO(),); orbit = j)
                fp, fa, fb = reducer_final(MEGNO(), s_c[2]), reducer_final(MEGNO(), s_a[1]), reducer_final(MEGNO(), s_b[1])
                @test all(g -> fp[g][1] == fa[g][1] && fp[g][2] == fb[g][1], flds)
            end
        end
        @info "AC8(c) realized ratio error / tolerance (R2)" label worstD worstR
    end
end

@testset "AC9 symmetry and column independence" begin
    _ = (RLI, PairMap, check_frame)
    # I3, I4. (a) Exchange of the two orbits (start points and block order) with the same xi_0: DL_k and RLI bitwise equal at every handoff, s_A and s_B exchanged; (b) a sign change
    # of xi_0 (both blocks, or one block): all outputs bitwise equal (the norms do not see the sign); (c) the seeded xi_0 does not depend on k (it is the first column of the seeded
    # frame for k = 1 and k = 2) and the initial columns of the block frame are bitwise the same vector in their blocks, with zeros in the other block; check_frame accepts them.
    # K = 2 pi and the K = 0.5 regular class, 5 starts, N = 400, tau in {1, 4}, Float64 and Double64 (ids 1..1000 in (c)).
    for T in RL_TYPES, label in (:two_pi, :half), tau in (1, 4)
        c = rl_case(label, T, 5)
        P = PairMap(c.sys)
        prs = rl_pairs(c, 1e-9)
        Ws = rl_frames(c)
        swap(z) = vcat(rl_half(z, 2), rl_half(z, 1))
        neg(W) = -W
        neg2(W) = SMatrix{4, 2, T}(W[1, 1], W[2, 1], W[3, 1], W[4, 1], -W[1, 2], -W[2, 2], -W[3, 2], -W[4, 2])
        base = rl_run((RLI(),), P, [p.z for p in prs], Ws, 400, tau; all = true)
        ex = rl_run((RLI(),), P, [swap(p.z) for p in prs], Ws, 400, tau; all = true)
        n1 = rl_run((RLI(),), P, [p.z for p in prs], neg.(Ws), 400, tau; all = true)
        n2 = rl_run((RLI(),), P, [p.z for p in prs], neg2.(Ws), 400, tau; all = true)
        for j in 1:5
            a, b = base[j][1].snaps, ex[j][1].snaps
            @test length(a) == 400 ÷ tau
            @test isequal([s.DL for s in a], [s.DL for s in b]) && isequal([s.RLI for s in a], [s.RLI for s in b]) && isequal([s.Q for s in a], [s.Q for s in b])
            @test isequal([reverse(s.sum) for s in a], [s.sum for s in b]) && isequal([reverse(s.comp) for s in a], [s.comp for s in b])
            @test isequal(a, n1[j][1].snaps) && isequal(a, n2[j][1].snaps)
        end
    end
    for T in RL_TYPES
        seed = rl_seed(9, 1)
        nid = T === Float64 ? 1000 : 200
        bad = 0
        for j in 1:nid
            xi = rl_xi0(seed, j, T)
            W = rl_frame(xi)
            same_k = initial_deviations(seed, j, Val(2), Val(1), T)[:, 1] == initial_deviations(seed, j, Val(2), Val(2), T)[:, 1]
            blocks = W[1:2, 1] == xi && W[3:4, 2] == xi && all(iszero, W[3:4, 1]) && all(iszero, W[1:2, 2])
            ok = check_frame(RLI(), W; orbit = j)
            bad += !(same_k && blocks && 0 <= ok <= 8 * rl_eps(T))
        end
        @test bad == 0
    end
end

# One BigFloat(256) pair run in wrapped mode on the start x0 and direction xi (RLI at every handoff 1..N, tau = 1), and the 512-bit reference at the realized points of both
# orbits: returns the errors and bounds at the checkpoints. Everything inside setprecision(BigFloat, 256) (the reference sets its own 512).
function rl_bf_check(label::Symbol, m::Int, N::Int, cps; Dx = 1e-12)
    return setprecision(BigFloat, RL_BITS) do
        T = BigFloat
        c = rl_case(label, T, m)
        epsT = rl_eps(T)
        out = NamedTuple[]
        for j in 1:m
            pr = rl_pair(c.sys, c.x0[j], c.xi0[j], Dx)
            snaps = rl_run((RLI(),), PairMap(c.sys), [pr.z], [rl_frame(c.xi0[j])], N, 1; all = true)[1][1].snaps
            o = rl_pair_ref(c.sys, c.K, pr.zA, pr.zB, c.xi0[j], N; bits = 512)
            ref = rl_ref_rli(o.sA, o.sB, 1; bits = 512)
            for n in cps
                PA, PB = maximum(rl_Psi(o, 1, g) for g in cps if g <= n), maximum(rl_Psi(o, 2, g) for g in cps if g <= n)
                B = rl_BR(epsT, PA, PB, n, n, Float64(ref.RLI[n]))
                push!(out, (j = j, n = n, err = rl_err(snaps[n].RLI, ref.RLI[n]; bits = 512), B = B, rli = snaps[n].RLI))
            end
        end
        out
    end
end

@testset "AC10 precision ladder and representation" begin
    _ = (RLI, PairMap)
    # (a) K in {2 pi, 5}, n = 1e3, and K in {5, 10}, n = 1e4, Dx = 1e-12, wrapped, Float64 against Double64 on the same initial points (exact widening): the ensemble means of RLI
    # agree within 4 combined standard errors (L1 S5). Sizes D38: 300 starts at n = 1e3 and 60 at n = 1e4 (the spec: 1000 and 200). (b) K = 2 pi, 50 starts, n in
    # {10, ..., 60}, BigFloat(256) on its own orbits against the 512-bit reference at its realized points: |RLI - RLI_ref| <= B_R(n). (c) BigFloat(256) against Double64, n = 300,
    # chaotic ensemble level (K in {2 pi, 5}, 100 starts), same rule as (a). (d) K = 10, n = 1e4, Float64: RLI finite, within B_R of the reference, at most 2 L_sigma,
    # and equal to the reference (not clipped or replaced).
    for (label, n, m) in ((:two_pi, 1000, RL_FULL ? 1000 : 300), (:five, 1000, RL_FULL ? 1000 : 300), (:five, 10_000, RL_FULL ? 200 : 60), (:ten, 10_000, RL_FULL ? 200 : 60))
        ens = Dict{Any, Vector{Float64}}()
        for T in RL_TYPES
            c = rl_case(label, T, m)
            prs = rl_pairs(c, 1e-12)
            r = rl_run((RLI(),), PairMap(c.sys), [p.z for p in prs], rl_frames(c), n, 1; steps = [n])
            ens[T] = Float64[f[1].snaps[1].RLI for f in r]
        end
        @test rl_agree(ens[Float64], ens[Double64])
        @info "AC10(a) ensemble means of RLI, Float64 and Double64" label n m rl_mean(ens[Float64]) rl_mean(ens[Double64])
    end
    cps = collect(10:10:60)
    rows = rl_bf_check(:two_pi, 50, 60, cps)
    @test rl_nbad(r -> r.err <= r.B, rows) == 0
    @info "AC10(b) realized ratio max error / B_R, BigFloat against 512 bits (R2)" worst = maximum(r.err / r.B for r in rows)
    for label in (:two_pi, :five)
        m = 100
        bfm = setprecision(BigFloat, RL_BITS) do
            c = rl_case(label, BigFloat, m)
            prs = rl_pairs(c, 1e-12)
            r = rl_run((RLI(),), PairMap(c.sys), [p.z for p in prs], rl_frames(c), 300, 1; steps = [300])
            Float64[Float64(f[1].snaps[1].RLI) for f in r]
        end
        c = rl_case(label, Double64, m)
        r = rl_run((RLI(),), PairMap(c.sys), [p.z for p in rl_pairs(c, 1e-12)], rl_frames(c), 300, 1; steps = [300])
        dm = Float64[Float64(f[1].snaps[1].RLI) for f in r]
        @test rl_agree(bfm, dm)
    end
    let (label, N, m) = last(rl_panel())
        rows = rl_rows(:ten, Float64, 1e-12, 1, m, N)
        top = filter(r -> r.n == N, rows)
        Lsig = log(rl_sigma(10.0))
        @test length(top) == m
        @test rl_nbad(r -> isfinite(r.RLI) && 0 <= r.RLI <= 2 * Lsig && r.eR <= r.BR, top) == 0
        @info "AC10(d) K = 10, n = 1e4: RLI values" rli = [r.RLI for r in top]
    end
end
