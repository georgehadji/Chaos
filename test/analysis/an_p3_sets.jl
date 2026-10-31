# D93d relabelling rule (E49 c), D93e Cal and the sets (E50, I35). Spec: analysis-spec.md, Amendment D93, criteria D93d (S part: the rule, wrapped against lifted) and D93e (S).
# The store-level parts of D93d (items of the P1 and p2 stores, missing data, the change table) are in an_p3_pipeline.jl. Seeds: p3_seed('d', 1) (random labels).

@testset "D93d the reduced rule at every boundary (E49 c): the rule of D89 E23 b plus the diagnostics condition" begin
    p2 = p2_params()                        # label keys of the p2 table: lambda_min 1e-2, growth ratio 5, rotation tolerance 1e-7, estimator bound 1e-8
    thr, rt, tol, bnd = p2.label_lambda_min, p2.label_growth_ratio, p2.label_rotation_tol, p2.label_estimator_bound
    @test (thr, rt, tol, bnd) == (1e-2, 5.0, 1e-7, 1e-8)
    lab(; lam5 = 0.5, tenth = 0.05, lam1 = lam5, rho1 = 0.1, rho2 = 0.1, d5 = false, d1 = false) =
        AN.p3_relabel((lam5 = lam5, lam5_tenth = tenth, lam1 = lam1, rho1 = rho1, rho2 = rho2, diag5 = d5, diag1 = d1), p2)
    # the two classes of a clean orbit
    @test lab().rho === :C && lab().instability == AN.POSITIVE && lab().structure == "irregular"
    @test lab(lam5 = 1e-3, tenth = 1e-2, lam1 = 1e-3).rho === :R && lab(lam5 = 1e-3, tenth = 1e-2, lam1 = 1e-3).instability == AN.NONE_DETECTED && lab(lam5 = 1e-3, tenth = 1e-2, lam1 = 1e-3).structure == "quasiperiodic"
    # lambda_L5(N_1) at label_lambda_min, one unit in the last place below and above (the growth ratio is 10 in all three)
    @test lab(lam5 = thr, tenth = thr / 10, lam1 = thr).rho === :C && lab(lam5 = prevfloat(thr), tenth = thr / 10, lam1 = prevfloat(thr)).rho === :U &&
          lab(lam5 = nextfloat(thr), tenth = thr / 10, lam1 = nextfloat(thr)).rho === :C
    # the growth ratio r = 10 lambda(N) / lambda(N / 10) at label_growth_ratio: r = 5 exactly counts as positive growth, one ulp of lambda(N / 10) above makes it neither
    @test lab(lam5 = 0.5, tenth = 1.0).rho === :C && lab(lam5 = 0.5, tenth = nextfloat(1.0)).rho === :U && lab(lam5 = 0.5, tenth = prevfloat(1.0)).rho === :C
    @test lab(lam5 = 1e-3, tenth = 2e-3, lam1 = 1e-3).rho === :U && lab(lam5 = 1e-3, tenth = nextfloat(2e-3), lam1 = 1e-3).rho === :R     # r = 5 on the regular side: not below 5 (U), just below 5 (R)
    # lambda_L5(N_1 / 10) not positive: the ratio is undefined, both axes unresolved
    @test lab(tenth = 0.0).rho === :U && lab(tenth = -0.1).rho === :U && lab(lam5 = 0.0, tenth = 0.0, lam1 = 0.0).rho === :U && lab(tenth = 0.0).instability == AN.UNRESOLVED
    # |rho_1 - rho_2| at label_rotation_tol (counts as quasiperiodic) and 1.01 times it (does not)
    reg(; kw...) = lab(; lam5 = 1e-3, tenth = 1e-2, lam1 = 1e-3, kw...)
    @test reg(rho1 = 0.0, rho2 = tol).rho === :R && reg(rho1 = 0.0, rho2 = 1.01 * tol).rho === :U && reg(rho1 = 0.0, rho2 = -tol).rho === :R && reg(rho1 = 0.0, rho2 = -1.01 * tol).rho === :U
    @test lab(rho1 = 0.0, rho2 = 1.0).rho === :C                      # a chaotic orbit is irregular whatever its rotation numbers (the structure axis follows the instability axis)
    # |lambda_L1 - lambda_L5| at label_estimator_bound (not above: kept) and one ulp above (both axes unresolved)
    @test lab(lam5 = 0.0, tenth = 0.1, lam1 = bnd).rho === :R && lab(lam5 = 0.0, tenth = 0.1, lam1 = nextfloat(bnd)).rho === :U && lab(lam5 = 0.0, tenth = 0.1, lam1 = nextfloat(bnd)).structure == AN.UNRESOLVED
    @test lab(lam1 = 0.5 + 1e-6).rho === :U && lab(lam1 = 0.5 - 1e-6).rho === :U
    # the diagnostics condition: the flag on the L5 item or on the L1 item read for the orbit leaves both axes unresolved
    @test lab(d5 = true).rho === :U && lab(d1 = true).rho === :U && lab(d5 = true).instability == AN.UNRESOLVED && lab(d5 = true).structure == AN.UNRESOLVED &&
          reg(d5 = true).rho === :U && reg(d1 = true).rho === :U && lab(d5 = false, d1 = false).rho === :C
    # without the flags the rule is the rule of the grid orbits of D89 (the same function of the same numbers), on 400 random inputs; with the flags never a class
    rng = StableRNG(p3_seed('d', 1))
    bad = 0; classes = Dict(:C => 0, :R => 0, :U => 0)
    for _ in 1:400
        lam5 = rand(rng) < 0.5 ? 10.0^(-4 + 3 * rand(rng)) : 10.0^(-1 * rand(rng) - 1.0)
        tenth = lam5 * (rand(rng) < 0.5 ? 1.0 : 10.0^(rand(rng) * 1.5)) * (rand(rng) < 0.05 ? -1 : 1)
        lam1 = lam5 + (rand(rng) < 0.2 ? 1e-7 : rand(rng) < 0.5 ? 1e-9 : 0.0)
        rho1 = rand(rng); rho2 = rho1 + (rand(rng) < 0.3 ? 1e-6 : rand(rng) < 0.5 ? 5e-8 : 0.0)
        d = (lam5 = lam5, lam5_tenth = tenth, lam1 = lam1, rho1 = rho1, rho2 = rho2)
        a = AN.p3_relabel(merge(d, (diag5 = false, diag1 = false)), p2)
        b = AN.p2_grid_label(d, p2)
        c = p3ref_label(merge(d, (diag5 = false, diag1 = false)), (label_lambda_min = thr, label_growth_ratio = rt, label_rotation_tol = tol, label_estimator_bound = bnd))
        classes[a.rho] += 1
        (a.rho === b.rho && a.instability == b.instability && a.structure == b.structure && a.rho === c.rho) || (bad += 1)
        AN.p3_relabel(merge(d, (diag5 = true, diag1 = false)), p2).rho === :U || (bad += 1)
    end
    @test bad == 0 && classes[:C] > 20 && classes[:R] > 20 && classes[:U] > 20
    # a wrapped fixture with winding counts gives the same rotation numbers as the same orbit lifted (P19), to 1e-12; the lifted coordinate of a wrapped value is x + 2 pi winding
    x0, p0, h, N = 0.7, 0.9, 40_000, 80_000
    lifted(n) = x0 + n * p0
    wrapped(n) = (w = floor(lifted(n) / (2pi)); (lifted(n) - 2pi * w, Int(w)))
    xh, wh = wrapped(h); xN, wN = wrapped(N)
    a = AN.p2_rotation_numbers(x0, AN.p2_lift(xh, wh), AN.p2_lift(xN, wN), h, N)
    b = AN.p2_rotation_numbers(x0, lifted(h), lifted(N), h, N)
    @test abs(a[1] - b[1]) < 1e-12 && abs(a[2] - b[2]) < 1e-12 && abs(a[1] - p0 / (2pi)) < 1e-12 && wh > 100 && wN > wh
end

# the orbit rows of the sets: (id, source, kindex, class, control, flags, rho, lambda_ref, lambda_ns); rho is :C, :R, :U or :missing (relabelling data missing)
p3_orb(id; source = :panel, kindex = 2, class = "C", control = false, flags = String[], rho = :C, lref = 0.5, lns = lref) =
    (id = id, source = source, kindex = kindex, class = class, control = control, flags = flags, rho = rho, lambda_ref = lref, lambda_ns = lns)

@testset "D93e Cal and the sets: every condition switched on and off, the boundaries, the conservation of counts (I35)" begin
    q = merge(p3_params(), (; gate_size = 3))
    p2 = p2_params()
    @test p2.gate_k_index == 1 && p2.gate_classes == ["T"]
    orbits = [
        p3_orb(1, lref = 1.0),                                              # S2 at the boundary (>=); not W2, not T_s (lambda_ns = lambda_ref)
        p3_orb(2, lref = prevfloat(1.0)),                                   # one ulp below: not S2
        p3_orb(3, lref = 0.1),                                              # W2 at the boundary (<=)
        p3_orb(4, lref = nextfloat(0.1)),                                   # one ulp above: not W2
        p3_orb(5, lref = 0.5, lns = 0.1),                                   # T_s at the boundary: 0.1 <= 0.2 x 0.5 (<=); T_s is W2 although lambda_ref is above 0.1
        p3_orb(6, lref = 0.5, lns = nextfloat(0.1)),                        # one ulp above the boundary: not T_s
        p3_orb(7, lref = 1.5, lns = 0.1),                                   # T_s and S2 and W2 at once (the overlap of S2 and W2 is reported)
        p3_orb(8, control = true, class = "control"),                       # dropped: control
        p3_orb(9, flags = ["drift"]), p3_orb(10, flags = ["excluded"]),     # dropped: drift, excluded
        p3_orb(11, rho = :U), p3_orb(12, rho = :missing),                   # dropped: unresolved, relabelling data missing
        p3_orb(13, control = true, class = "control", flags = ["drift"]),   # the reason is the first: control
        p3_orb(14, flags = ["drift"], rho = :U),                            # drift before unresolved
        p3_orb(15, rho = :R, class = "T", kindex = 1, lref = 1.4e-4),       # a gate orbit (K index 1, class T) in Cal: in the regular set
        p3_orb(16, rho = :R, class = "T", kindex = 1, lref = 1.4e-4),       # a gate orbit in Cal
        p3_orb(17, rho = :U, class = "T", kindex = 1, lref = 1.4e-4),       # a gate orbit that is dropped (unresolved): still in the gate set (not intersected with Cal)
        p3_orb(18, rho = :R, source = :grid, class = "h4", kindex = 2, lref = 1.4e-4), p3_orb(19, source = :grid, class = "h4", kindex = 1, lref = 2.0),
        p3_orb(20, rho = :R, class = "S", kindex = 1, lref = 1.4e-4)]       # panel, K index 1, class S: regular but not in the gate classes
    S = AN.p3_sets(orbits, q, p2)
    @test S.cal == [1, 2, 3, 4, 5, 6, 7, 15, 16, 18, 19, 20]
    @test S.S2 == [1, 7, 19] && S.T_s == [5, 7] && S.W2 == [3, 5, 7] && S.R_W == [3] && S.regular == [15, 16, 18, 20] && S.gate == [15, 16, 17] && S.overlap == [7]
    @test [(d.id, d.reason) for d in S.dropped] == [(8, "control"), (9, "drift"), (10, "excluded"), (11, "unresolved"), (12, "missing"), (13, "control"), (14, "drift"), (17, "unresolved")]
    # I35: Cal plus the dropped orbits are all the orbits, each once; T_s is a subset of W2, R_W = W2 minus T_s; the sets are subsets of the C orbits of Cal; sorted
    @test sort(vcat(S.cal, [d.id for d in S.dropped])) == 1:20 && issubset(S.T_s, S.W2) && S.R_W == setdiff(S.W2, S.T_s) && issubset(vcat(S.S2, S.W2, S.T_s), S.cal) &&
          all(id -> orbits[id].rho === :C, vcat(S.S2, S.W2, S.T_s)) && all(id -> orbits[id].rho === :R, S.regular) && issorted(S.cal) && issorted(S.W2) && issorted(S.gate)
    # a gate set of another size than gate_size is refused, naming the gate set
    @test an_rejects(() -> AN.p3_sets(orbits, merge(q, (; gate_size = 2)), p2), ["gate"]) && an_rejects(() -> AN.p3_sets(orbits[1:14], q, p2), ["gate"])
    # the keys act: f_s, lam_s2 and lam_w2 of the p3 table (not the p2 names)
    S2b = AN.p3_sets(orbits, merge(q, (; f_s = 0.1)), p2)
    @test S2b.T_s == [7] && S2b.W2 == [3, 7]
    S3 = AN.p3_sets(orbits, merge(q, (; lam_s2 = 1.5, lam_w2 = 0.05)), p2)
    @test S3.S2 == [7, 19] && S3.W2 == [5, 7] && S3.T_s == [5, 7]
    # random tables against the reference rules (the three lines of E50 with the p3 keys)
    rng = StableRNG(p3_seed('e', 1))
    for _ in 1:20
        os = [p3_orb(i, rho = rand(rng) < 0.7 ? :C : :R, lref = rand(rng) < 0.3 ? 0.1 : round(2 * rand(rng); digits = 2), lns = rand(rng) < 0.3 ? 0.05 : round(rand(rng); digits = 2), kindex = rand(rng, 1:2),
                       class = rand(rng) < 0.5 ? "T" : "C", flags = rand(rng) < 0.1 ? ["drift"] : String[]) for i in 1:40]
        g = [o.id for o in os if o.kindex == 1 && o.class == "T"]
        qq = merge(q, (; gate_size = length(g)))
        T = AN.p3_sets(os, qq, p2)
        keep = [o for o in os if !("drift" in o.flags)]
        C = [o for o in keep if o.rho === :C]
        Ts = sort([o.id for o in C if o.lambda_ns <= 0.2 * o.lambda_ref])
        @test T.T_s == Ts && T.S2 == sort([o.id for o in C if o.lambda_ref >= 1.0]) && T.W2 == sort(unique(vcat([o.id for o in C if o.lambda_ref <= 0.1], Ts))) &&
              T.gate == g && T.regular == sort([o.id for o in keep if o.rho === :R]) && T.cal == sort([o.id for o in keep])
    end
end

@testset "D93e sets of 19 and 20 orbits against min_set give 'too small' and a selection (E54 C8, P49)" begin
    q = p3_params(; min_set = 20, e_cal = 1.0, a_cal = 1.0, s_cal = 1.0, c_min = 0.0, delta_tol = 0.0)
    q = merge(q, (; voter_bits = ["SALI", "FLI", "MEGNO", "E1", "E2"], indicators_i7 = ["LE", "SALI", "FLI", "MEGNO", "E1", "E2"], min_voters = 4, w_range = [0], v_range = [0]))
    votes(c) = Dict{Symbol, String}(X => (c ? "c"^10 : "r"^10) for X in (:SALI, :FLI, :MEGNO, :E2))
    function table(nT, nS2, nW)
        rows = NamedTuple[]
        id = 0
        add(; kw...) = (id += 1; push!(rows, p3_row(id; fill = "u"^10, vs = votes(true), kw...)))
        for _ in 1:nS2; add(S2 = true); end
        for _ in 1:nT; add(Ts = true); end
        for _ in 1:nW; add(W2 = true); end
        for _ in 1:3; id += 1; push!(rows, p3_row(id; fill = "u"^10, vs = votes(false), rho = :R)); end
        return rows
    end
    s19 = AN.p3_select(table(19, 20, 5), q)
    s20 = AN.p3_select(table(20, 20, 5), q)
    s2_19 = AN.p3_select(table(20, 19, 5), q)
    @test s19.outcome === :sets_too_small && s19.star === nothing && isempty(s19.tol) && s2_19.outcome === :sets_too_small
    @test s20.outcome === :selected && s20.star !== nothing && s20.outcome === p3ref_select(table(20, 20, 5), q).outcome
end
