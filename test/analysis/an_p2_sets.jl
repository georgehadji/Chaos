# D89k sets (E29, I16) and D89l reference class (E23). Spec: analysis-spec.md, Amendment D89, E23, E24, E29, criteria D89k and D89l.

# one orbit record of the set rules; defaults are a chaotic panel orbit of K index 1 that is in no set
orb(id; source = :panel, kindex = 1, K = 0.5, class = "C", control = false, flags = String[], rho = :C, lref = 0.5, lns = 0.5, sd = false, sl = false) =
    (id = id, source = source, kindex = kindex, K = K, class = class, control = control, flags = flags, rho = rho, lambda_ref = lref, lambda_ns = lns, in_SD = sd, in_SL = sl)

@testset "D89k sets: each condition of E29 switched on and off, boundaries" begin
    p = p2_params()                                   # s2_lambda_min 1.0, w2_lambda_max 0.1, f_s 0.2, gate: K index 1, class T, size 2
    gate = [orb(900, class = "T", rho = :R, lref = 1e-4, lns = 1e-4), orb(901, class = "T", rho = :R, lref = 1e-4, lns = 1e-4)]
    os = vcat(gate, [
        orb(1, lref = 1.0, lns = 1.0),                              # S2 at the boundary (>=)
        orb(2, lref = prevfloat(1.0), lns = prevfloat(1.0)),        # one ulp below: not S2, not W2 (above 0.1), not sticky
        orb(3, lref = 0.1, lns = 0.1),                              # W2 at the boundary (<=)
        orb(4, lref = nextfloat(0.1), lns = nextfloat(0.1)),        # one ulp above: not W2
        orb(5, lref = 0.5, lns = 0.1),                              # T_s at the boundary: 0.1 <= 0.2 x 0.5 (<=); T_s is W2 although lambda_ref is above 0.1
        orb(6, lref = 0.5, lns = nextfloat(0.1)),                   # one ulp above the boundary: not sticky, not W2
        orb(7, lref = 1.5, lns = 0.1),                              # sticky and strong: in T_s, W2 and S2 (the overlap is counted)
        orb(8, rho = :R, lref = 0.05, lns = 0.01),                  # regular: in no chaotic set whatever its numbers
        orb(9, rho = :U, lref = 0.05, lns = 0.01),                  # unresolved: dropped
        orb(10, control = true, rho = :C, lref = 0.05),             # control: dropped
        orb(11, flags = ["drift"], lref = 0.05),                    # flagged: dropped, reason drift
        orb(12, flags = ["excluded"], lref = 0.05),
        orb(13, flags = ["drift"], rho = :U, lref = 0.05),          # flagged and unresolved: listed once, with the reason drift
        orb(14, control = true, flags = ["drift"], rho = :U),       # a control: listed once, reason control
        orb(15, class = "S", lref = 0.4, lns = 0.4),                # design class S, rho C, panel: in T_s_design, not in T_s
        orb(16, class = "S", rho = :R, lref = 0.4, lns = 0.4),      # design class S but regular: not in T_s_design
        orb(100, source = :grid, kindex = 1, class = "h4", lref = 0.05, lns = 0.05),         # grid orbit: in E8
        orb(101, source = :grid, kindex = 2, class = "h4", rho = :R, lref = 0.05, lns = 0.05),
        orb(20, sd = true, sl = true, lref = 0.05, lns = 0.05),     # S_D and S_L
        orb(21, sd = true, rho = :R), orb(22, sl = true, rho = :U, sd = true)])
    s = AN.p2_sets(os, p)
    @test s.E_all == [1, 2, 3, 4, 5, 6, 7, 8, 15, 16, 20, 21, 100, 101, 900, 901]
    @test s.S2 == [1, 7] && s.T_s == [5, 7] && s.W2 == [3, 5, 7, 20, 100]
    @test issubset(s.T_s, s.W2) && s.overlap == [7] && sort(intersect(s.S2, s.W2)) == [7]
    @test s.T_s_design == [15] && s.regular == [8, 16, 21, 101, 900, 901] && s.E8 == [20, 100, 101]
    @test s.S_D == [20, 21] && s.S_L == [20] && s.gate == [900, 901]
    @test [(d.id, d.reason) for d in s.dropped] == [(9, "unresolved"), (10, "control"), (11, "drift"), (12, "excluded"), (13, "drift"), (14, "control"), (22, "unresolved")]
    # I16: E_all plus the dropped orbits are P union G; W2, S2, T_s are subsets of the chaotic orbits of E_all; E8 is a subset of E_all
    ids = sort([o.id for o in os])
    @test sort(vcat(s.E_all, [d.id for d in s.dropped])) == ids && allunique(vcat(s.E_all, [d.id for d in s.dropped]))
    chaotic = [o.id for o in os if o.id in s.E_all && o.rho === :C]
    @test issubset(s.W2, chaotic) && issubset(s.S2, chaotic) && issubset(s.T_s, chaotic) && issubset(s.E8, s.E_all)
    # the gate set is not intersected with E_all or the labels, and its size must be gate_size
    g2 = vcat([orb(900, class = "T", rho = :U, flags = ["drift"])], gate[2:2])
    @test AN.p2_sets(g2, p).gate == [900, 901]
    @test an_rejects(() -> AN.p2_sets(gate[1:1], p), ["gate"]) && an_rejects(() -> AN.p2_sets(vcat(gate, [orb(902, class = "T")]), p), ["gate"])
    # the gate set takes only panel orbits of the gate K index, outside the controls, of the gate classes
    mix = vcat(gate, [orb(903, class = "T", kindex = 2), orb(904, class = "T", control = true), orb(905, class = "X"), orb(906, class = "T", source = :grid)])
    @test AN.p2_sets(mix, p).gate == [900, 901]
    # the thresholds come from the table: another s2_lambda_min and f_s give other sets
    s2 = AN.p2_sets(os, merge(p, (; s2_lambda_min = 1.2, f_s = 0.1)))                         # orbit 5: 0.1 <= 0.1 x 0.5 fails; orbit 7: 0.1 <= 0.15
    @test s2.S2 == [7] && s2.T_s == [7] && s2.W2 == [3, 7, 20, 100]
end

@testset "D89k not decidable below min_set: 19 and 20 orbits (the rule E29 uses the table value)" begin
    p = p2_params(min_set = 20)
    mk(n) = [orb(i, lref = 0.05, lns = 0.05) for i in 1:n]
    for n in (19, 20)
        s = AN.p2_sets(vcat(mk(n), [orb(900, class = "T", rho = :R), orb(901, class = "T", rho = :R)]), p)
        @test length(s.W2) == n && AN.p2_decidable(length(s.W2), p) == (n >= 20)
    end
    @test AN.p2_decidable(19, p) === false && AN.p2_decidable(20, p) === true && AN.p2_decidable(0, p) === false && AN.p2_decidable(5, p2_params()) === true
end

@testset "D89l reference class: the panel rule and the grid rule of E23" begin
    p = p2_params()                                   # label_lambda_min 0.01, growth ratio 5, rotation tolerance 1e-7, estimator bound 1e-8
    POS, NONE, UNR, IRR, QP = "positive growth supported", "none detected on the tested horizon", "unresolved", "irregular", "quasiperiodic"
    @test AN.p2_rho(POS, IRR) === :C && AN.p2_rho(NONE, QP) === :R
    for (i, s) in ((POS, QP), (POS, UNR), (POS, "periodic"), (NONE, IRR), (NONE, UNR), (NONE, "periodic"), (UNR, IRR), (UNR, QP), (UNR, UNR), (UNR, "periodic"))
        @test AN.p2_rho(i, s) === :U
    end
    lab(; lam5 = 0.5, tenth = 0.05, lam1 = lam5, rho1 = 0.0, rho2 = 0.0, q = p) = AN.p2_grid_label((lam5 = lam5, lam5_tenth = tenth, lam1 = lam1, rho1 = rho1, rho2 = rho2), q)
    # instability: lambda_L5(N_1) >= label_lambda_min and r >= 5 (r = 10 lambda(N_1) / lambda(N_1 / 10)); the four cases of E14 on the grid
    @test (l = lab(lam5 = 0.0625, tenth = 0.125); (l.instability, l.structure, l.rho) == (POS, IRR, :C))                     # r = 5 exactly
    @test (l = lab(lam5 = 0.0625, tenth = nextfloat(0.125)); (l.instability, l.structure, l.rho) == (UNR, UNR, :U))        # r just below 5, lambda above the floor
    @test (l = lab(lam5 = 0.01, tenth = 0.002); (l.instability, l.rho) == (POS, :C))                                         # lambda at the floor, r = 50
    @test (l = lab(lam5 = prevfloat(0.01), tenth = 0.002); (l.instability, l.structure, l.rho) == (UNR, UNR, :U))          # lambda below the floor, r above 5
    @test (l = lab(lam5 = 1e-4, tenth = 1.1e-3); (l.instability, l.structure, l.rho) == (NONE, QP, :R))                        # none detected, rotation numbers equal
    @test (l = lab(lam5 = 1e-4, tenth = 0.0); (l.instability, l.structure, l.rho) == (UNR, UNR, :U))                          # lambda_L5(N_1 / 10) <= 0: r undefined
    @test (l = lab(lam5 = 1e-4, tenth = -1e-4); (l.instability, l.rho) == (UNR, :U))
    # structure: quasiperiodic iff none detected and |rho_1 - rho_2| <= tolerance (at the tolerance and 1.01 times it)
    @test (l = lab(lam5 = 1e-4, tenth = 1.1e-3, rho1 = 0.0, rho2 = 1e-7); (l.structure, l.rho) == (QP, :R))
    @test (l = lab(lam5 = 1e-4, tenth = 1.1e-3, rho1 = 0.0, rho2 = 1.01e-7); (l.instability, l.structure, l.rho) == (NONE, UNR, :U))
    @test (l = lab(rho1 = 0.0, rho2 = 0.3); (l.structure, l.rho) == (IRR, :C))                                                  # irregular does not look at the rotation numbers
    # the disagreement test: |lambda_L1 - lambda_L5| above the bound gives unresolved on both axes (2^-26 is exact at 0.5)
    q = merge(p, (; label_estimator_bound = 2.0^-26))
    @test (l = lab(lam5 = 0.5, tenth = 0.05, lam1 = 0.5 + 2.0^-26, q = q); (l.instability, l.structure, l.rho) == (POS, IRR, :C))
    @test (l = lab(lam5 = 0.5, tenth = 0.05, lam1 = 0.5 + 2.0^-26 + 2.0^-52, q = q); (l.instability, l.structure, l.rho) == (UNR, UNR, :U))
    @test (l = lab(lam5 = 0.5, tenth = 0.05, lam1 = 0.5 - 2.0^-26 - 2.0^-53, q = q); (l.instability, l.structure, l.rho) == (UNR, UNR, :U))
    # every grid orbit carries ref_weak in the tables (pipeline test); the rotation numbers: rho_1 = (x_h - x_0) / (2 pi h), rho_2 = (x_N - x_h) / (2 pi (N - h))
    r1, r2 = AN.p2_rotation_numbers(0.5, 0.5 + 0.7 * 40_000, 0.5 + 0.7 * 80_000, 40_000, 80_000)
    @test isapprox(r1, 0.7 / (2 * pi); rtol = 1e-14) && isapprox(r2, 0.7 / (2 * pi); rtol = 1e-14)
    r1, r2 = AN.p2_rotation_numbers(0.0, 10.0, 40.0, 5, 15)                       # (10 - 0) / (2 pi 5) and (40 - 10) / (2 pi 10)
    @test isapprox(r1, 10 / (2 * pi * 5); rtol = 1e-15) && isapprox(r2, 30 / (2 * pi * 10); rtol = 1e-15)
    # a wrapped orbit with winding counts gives the same rotation numbers (to 1e-12) as the same orbit stored lifted
    t = [1, 100, 40_000, 80_000]
    x = 0.3 .+ 0.7 .* t                                                            # lifted
    w = Int64.(floor.(x ./ (2 * pi))); xw = x .- (2 * pi) .* w                     # wrapped value and winding count
    @test all(0 .<= xw .< 2 * pi) && maximum(w) > 1000
    xl = AN.p2_lift(xw, w)
    @test isapprox(xl, x; atol = 1e-8)
    a = AN.p2_rotation_numbers(0.3, x[3], x[4], 40_000, 80_000); b = AN.p2_rotation_numbers(0.3, xl[3], xl[4], 40_000, 80_000)
    @test abs(a[1] - b[1]) < 1e-12 && abs(a[2] - b[2]) < 1e-12
end
