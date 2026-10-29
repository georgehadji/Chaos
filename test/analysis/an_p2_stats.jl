# D89m statistics (F_X, M_X, Q, R_X, phi_X, Delta_phi, delta_dir, delta_prec, F_med, Phi_X) and D89n decisions H1..H6, ROBUST. Pure functions on orbit records.
# Spec: analysis-spec.md, Amendment D89, E32, E33, criteria D89m and D89n. A record is (id, V::Dict{Symbol, 10 verdicts}, T::Dict{Symbol, censored time}).

# a record whose verdicts are regular except at index 7 (budget b_P, k = 6) and index 10 (k = 9): vv[X] = (verdict at b_P, verdict at k = 9)
function p2_rec(id, vv::AbstractDict; T = Dict{Symbol, Int}())
    V = Dict{Symbol, Vector{Symbol}}()
    for (X, (a, b)) in vv
        v = fill(:regular, 10); v[P2_KP] = a; v[10] = b
        V[X] = v
    end
    return (id = id, V = V, T = Dict{Symbol, Int}(T))
end

@testset "D89m F_X, F_med (odd and even number of indicators), and the undefined cases" begin
    p = p2_params()
    # four orbits, F_X = fraction of the orbits decided at b_P and at k = 9 whose verdicts differ
    #   LE: (c,c) (r,c) (c,r) (c,c): 2 of 4 = 0.5       SALI: all equal: 0           FLI: orbit 1 undecided at b_P -> 3 decided, 1 differs: 1/3
    #   MEGNO: 0                                         E1: orbit 4 undecided at k = 9 -> 3 decided, all differ: 1      E2: 0      WBA: 2 of 4 = 0.5
    rows = Dict(:LE => [(c_, c_), (r_, c_), (c_, r_), (c_, c_)], :SALI => fill((r_, r_), 4), :FLI => [(u_, c_), (r_, r_), (c_, r_), (c_, c_)], :MEGNO => fill((c_, c_), 4),
                :E1 => [(r_, c_), (c_, r_), (r_, c_), (c_, u_)], :E2 => fill((r_, r_), 4), :WBA => [(c_, r_), (r_, c_), (r_, r_), (c_, c_)])
    recs = [p2_rec(o, Dict(X => rows[X][o] for X in keys(rows))) for o in 1:4]
    s = AN.p2_stat(:H5, (recs,), p)
    @test s.n == 4 && s.decided == Dict(:LE => 4, :SALI => 4, :FLI => 3, :MEGNO => 4, :E1 => 3, :E2 => 4, :WBA => 4)
    @test s.F[:LE] == 0.5 && s.F[:SALI] == 0.0 && isapprox(s.F[:FLI], 1 / 3; atol = 1e-15) && s.F[:MEGNO] == 0.0 && s.F[:E1] == 1.0 && s.F[:E2] == 0.0 && s.F[:WBA] == 0.5
    # the seven F_X sorted: 0, 0, 0, 1/3, 0.5, 0.5, 1 -> the median of an odd number of indicators is the fourth: 1/3
    @test isapprox(s.F_med, 1 / 3; atol = 1e-15)
    # eight indicators (I_8): RLI with F = 1/4 (one of four differs): sorted 0, 0, 0, 1/4, 1/3, 1/2, 1/2, 1 -> the mean of the two middle values: (1/4 + 1/3) / 2
    recs8 = [p2_rec(o, merge(Dict(X => rows[X][o] for X in keys(rows)), Dict(:RLI => [(c_, r_), (r_, r_), (r_, r_), (c_, c_)][o]))) for o in 1:4]
    s8 = AN.p2_stat(:H5, (recs8,), p; I = P2_X8)
    @test s8.F[:RLI] == 0.25 && isapprox(s8.F_med, (1 / 4 + 1 / 3) / 2; atol = 1e-15)
    # an indicator with no orbit decided at both budgets has no F_X, hence no F_med
    recs0 = [p2_rec(o, Dict(X => (X === :E1 ? (u_, c_) : rows[X][o]) for X in keys(rows))) for o in 1:4]
    s0 = AN.p2_stat(:H5, (recs0,), p)
    @test s0.F[:E1] === nothing && s0.F_med === nothing && s0.decided[:E1] == 0
    # an orbit without the item of an indicator (no V entry) is not counted for it
    rn = [p2_rec(1, Dict(:LE => (c_, r_))), p2_rec(2, Dict(:LE => (c_, c_), :SALI => (c_, c_)))]
    sn = AN.p2_stat(:H5, (rn,), p; I = [:LE, :SALI])
    @test sn.decided == Dict(:LE => 2, :SALI => 1) && sn.F[:LE] == 0.5 && sn.F[:SALI] == 0.0
end

@testset "D89m phi_X and Delta_phi on the sticky set (LE left out, also in the extension)" begin
    p = p2_params()
    # five orbits, verdict at k = 9: SALI r r r c u -> 3 of 4 decided regular; FLI all c -> 0; MEGNO all r -> 1; E1 r r c c c -> 2/5; E2 c c r r c -> 2/5; WBA r r r r c -> 4/5
    col = Dict(:SALI => [r_, r_, r_, c_, u_], :FLI => fill(c_, 5), :MEGNO => fill(r_, 5), :E1 => [r_, r_, c_, c_, c_], :E2 => [c_, c_, r_, r_, c_], :WBA => [r_, r_, r_, r_, c_],
               :LE => fill(c_, 5))
    recs = [p2_rec(o, Dict(X => (r_, col[X][o]) for X in keys(col))) for o in 1:5]
    s = AN.p2_stat(:H4, (recs,), p)
    @test s.n == 5 && !haskey(s.phi, :LE) && Set(keys(s.phi)) == Set(P2_X7[2:end]) && s.decided[:SALI] == 4 && s.decided[:E1] == 5
    @test s.phi[:SALI] == 0.75 && s.phi[:FLI] == 0.0 && s.phi[:MEGNO] == 1.0 && s.phi[:E1] == 0.4 && s.phi[:E2] == 0.4 && s.phi[:WBA] == 0.8
    @test s.dphi == 1.0                                                              # max 1 (MEGNO) minus min 0 (FLI)
    # with RLI (I_8) LE stays out and RLI is in: seven indicators
    col[:RLI] = [c_, c_, c_, c_, r_]
    recs8 = [p2_rec(o, Dict(X => (r_, col[X][o]) for X in keys(col))) for o in 1:5]
    s8 = AN.p2_stat(:H4, (recs8,), p; I = P2_X8)
    @test !haskey(s8.phi, :LE) && s8.phi[:RLI] == 0.2 && length(s8.phi) == 7
    # an indicator undecided for every orbit has no phi, hence no Delta_phi
    recsu = [p2_rec(o, Dict(X => (r_, X === :E2 ? u_ : col[X][o]) for X in keys(col))) for o in 1:5]
    su = AN.p2_stat(:H4, (recsu,), p)
    @test su.phi[:E2] === nothing && su.dphi === nothing && su.decided[:E2] == 0
end

@testset "D89m M_X, Q (H2) and R_X (H3): medians of odd and even sizes" begin
    p = p2_params()
    T = Dict(:S => Dict(:LE => [40, 40, 80, 160, 40960], :SALI => fill(40, 5), :FLI => [80, 160, 160, 320, 640], :MEGNO => [40, 80, 160, 320, 640],
                        :E1 => [640, 1280, 2560, 5120, 40960], :E2 => fill(160, 5), :WBA => [320, 640, 1280, 2560, 5120]),
             :W => Dict(:LE => [160, 320, 640, 40960], :SALI => [40, 80, 80, 160], :FLI => fill(320, 4), :MEGNO => [160, 160, 320, 640],
                        :E1 => [5120, 10240, 40960, 40960], :E2 => fill(1280, 4), :WBA => fill(40960, 4)))
    mk(g, n) = [p2_rec(o, Dict(X => (r_, c_) for X in P2_X7); T = Dict(X => T[g][X][o] for X in P2_X7)) for o in 1:n]
    S, W = mk(:S, 5), mk(:W, 4)
    # M_X(S) (odd, 5): 80, 40, 160, 160, 2560, 160, 1280; Q = 2560 / 40 = 64
    h2 = AN.p2_stat(:H2, (S,), p)
    @test h2.n == 5 && h2.M == Dict(:LE => 80.0, :SALI => 40.0, :FLI => 160.0, :MEGNO => 160.0, :E1 => 2560.0, :E2 => 160.0, :WBA => 1280.0) && h2.Q == 64.0
    # M_X(W) (even, 4): (320 + 640) / 2 = 480; 80; 320; (160 + 320) / 2 = 240; (10240 + 40960) / 2 = 25600; 1280; 40960. R_X = M_W / M_S
    h3 = AN.p2_stat(:H3, (W, S), p)
    @test h3.n_w == 4 && h3.n_s == 5 && h3.M_w == Dict(:LE => 480.0, :SALI => 80.0, :FLI => 320.0, :MEGNO => 240.0, :E1 => 25600.0, :E2 => 1280.0, :WBA => 40960.0)
    @test h3.R == Dict(:LE => 6.0, :SALI => 2.0, :FLI => 2.0, :MEGNO => 1.5, :E1 => 10.0, :E2 => 8.0, :WBA => 32.0) && h3.M_s == h2.M
    # an orbit without the item of an indicator is left out of that M_X
    Sn = vcat(S[1:4], [(id = 9, V = Dict(:SALI => fill(r_, 10)), T = Dict(:SALI => 40960))])
    @test AN.p2_stat(:H2, (Sn,), p; I = [:LE, :SALI]).M[:LE] == 60.0 && AN.p2_stat(:H2, (Sn,), p; I = [:LE, :SALI]).M[:SALI] == 40.0
end

@testset "D89m delta_dir, delta_prec and the flip count Phi_X; dropped and counted" begin
    p = p2_params()
    # direction: four orbits of S_D*, the verdicts at k = 9 in the five directions
    #  LE: orbit 2 differs in one direction -> 1/4; SALI: orbit 3 undecided in one direction is dropped and counted (3 decided, no flip -> 0);
    #  FLI: orbits 1 and 4 not identical -> 2/4; MEGNO: 0
    d5(s) = p2_sym(s)
    dir = [(id = 1, v = Dict(:LE => d5("ccccc"), :SALI => d5("ccccc"), :FLI => d5("crccc"), :MEGNO => d5("rrrrr"))),
           (id = 2, v = Dict(:LE => d5("cccc" * "r"), :SALI => d5("ccccc"), :FLI => d5("ccccc"), :MEGNO => d5("rrrrr"))),
           (id = 3, v = Dict(:LE => d5("rrrrr"), :SALI => d5("rrurr"), :FLI => d5("rrrrr"), :MEGNO => d5("ccccc"))),
           (id = 4, v = Dict(:LE => d5("ccccc"), :SALI => d5("ccccc"), :FLI => d5("rcrrr"), :MEGNO => d5("rrrrr")))]
    # precision: five orbits of S_L*, Float64 against Double64 at k = 9: LE (c,c) (c,r) (r,r) (r,c) (c,c): 2/5; FLI (c,c) (r,r) (c,r) (r,r) (c,u): 4 decided, 1 differs: 1/4
    f(v...) = Dict(zip(P2_X7, v))
    base(o) = Dict(X => r_ for X in P2_X7)
    prec = [(id = o, f64 = merge(base(o), Dict(:LE => a, :FLI => b)), d64 = merge(base(o), Dict(:LE => a2, :FLI => b2))) for (o, a, a2, b, b2) in
            ((1, c_, c_, c_, c_), (2, c_, r_, r_, r_), (3, r_, r_, c_, r_), (4, r_, c_, r_, r_), (5, c_, c_, c_, u_))]
    h = AN.p2_h6_stats(dir, prec, p)
    @test h.dir.n == 4 && h.dir.delta == Dict(:LE => 0.25, :SALI => 0.0, :FLI => 0.5, :MEGNO => 0.0) && h.dir.dropped == Dict(:LE => 0, :SALI => 1, :FLI => 0, :MEGNO => 0)
    @test h.prec.n == 5 && h.prec.delta[:LE] == 0.4 && h.prec.delta[:FLI] == 0.25 && h.prec.dropped[:FLI] == 1 && all(h.prec.delta[X] == 0.0 for X in P2_X7 if X ∉ (:LE, :FLI))
    @test length(h.prec.delta) == 7 && length(h.dir.delta) == 4
    # Phi_X: changes between consecutive indices with both verdicts decided
    @test AN.p2_flips(p2_sym("cccrrcrrcc")) == 4 && AN.p2_flips(p2_sym("cucrrrrrrc")) == 2 && AN.p2_flips(p2_sym("rrrrrrrrrr")) == 0 && AN.p2_flips(p2_sym("crcrcrcrcr")) == 9
    @test AN.p2_flips(p2_sym("uuuuuuuuuu")) == 0 && AN.p2_flips(p2_sym("cuuuuuuuuc")) == 0
end

@testset "D89m statistics against the independent implementation on random record sets" begin
    rng = StableRNG(p2_seed('m', 1))
    p = p2_params()
    syms = (c_, r_, u_)
    bad = String[]
    for t in 1:200
        n = rand(rng, 3:25)
        I = t % 2 == 0 ? P2_X7 : P2_X8
        recs = [(id = o, V = Dict(X => Symbol[syms[rand(rng, 1:3)] for _ in 1:10] for X in I if rand(rng) < 0.95),
                 T = Dict(X => 40 * 2^rand(rng, 0:7) for X in I if rand(rng) < 0.95)) for o in 1:n]
        recs = [(id = r.id, V = r.V, T = Dict(X => v for (X, v) in r.T if haskey(r.V, X))) for r in recs]
        for h in (:H1, :H2, :H4, :H5)
            h === :H2 && !all(X -> any(r -> haskey(r.T, X), recs), I) && continue
            a, b = AN.p2_stat(h, (recs,), p; I = I), p2ref_stat(h, (recs,), p; I = I)
            for k in propertynames(b)
                x, y = getproperty(a, k), getproperty(b, k)
                same(u, v) = u isa AbstractDict ? (keys(u) == keys(v) && all(kk -> same(u[kk], v[kk]), keys(u))) : (u isa Real && v isa Real) ? isapprox(u, v; rtol = 1e-12, atol = 1e-15) : isequal(u, v)
                same(x, y) || push!(bad, "$t $h $k")
            end
        end
    end
    @test isempty(bad)
end

@testset "D89n decisions: thresholds exactly, one ulp on the wrong side, not decidable at min_set - 1" begin
    p = p2_params()                                                                                   # min_set 5
    dec(h, s) = AN.p2_decide(h, s, p).verdict
    # H1: U >= 0.20, lower end >= 0.10, A_min <= 0.80 (all non-strict)
    h1(; n = 5, U = 0.2, lo = 0.1, A = 0.8) = (n = n, U = U, U_lo = lo, A_min = A)
    @test dec(:H1, h1()) === :supported && dec(:H1, h1(U = prevfloat(0.2))) === :not_supported && dec(:H1, h1(lo = prevfloat(0.1))) === :not_supported && dec(:H1, h1(A = nextfloat(0.8))) === :not_supported
    @test dec(:H1, h1(U = nextfloat(0.2), lo = nextfloat(0.1), A = prevfloat(0.8))) === :supported && dec(:H1, h1(n = 4)) === :not_decidable && dec(:H1, h1(lo = nothing)) === :not_decidable
    @test dec(:H1, h1(n = 4, U = 1.0)) === :not_decidable && dec(:H1, h1(A = nothing)) === :not_decidable && dec(:H1, h1(U = nothing)) === :not_decidable
    # H2: Q >= 10 (non-strict), lower end > 3 (strict)
    h2(; n = 5, Q = 10.0, lo = nextfloat(3.0)) = (n = n, Q = Q, Q_lo = lo)
    @test dec(:H2, h2()) === :supported && dec(:H2, h2(lo = 3.0)) === :not_supported && dec(:H2, h2(Q = prevfloat(10.0))) === :not_supported && dec(:H2, h2(n = 4)) === :not_decidable
    @test dec(:H2, h2(lo = nothing)) === :not_decidable
    # H3: at least 5 of the 7 indicators with R_X >= 10 (non-strict) and lower end > 3 (strict); not decidable if W2 or S2 is below min_set
    good = fill(10.0, 7); glo = fill(nextfloat(3.0), 7)
    h3(R, lo; nw = 5, ns = 5) = (n_w = nw, n_s = ns, R = R, R_lo = lo)
    mut(v, i, x) = (w = copy(v); w[i] = x; w)
    R5 = vcat(fill(10.0, 5), 1.0, 1.0); L5 = vcat(fill(nextfloat(3.0), 5), 1.0, 1.0)
    @test dec(:H3, h3(R5, L5)) === :supported && dec(:H3, h3(good, glo)) === :supported
    @test dec(:H3, h3(mut(R5, 1, prevfloat(10.0)), L5)) === :not_supported && dec(:H3, h3(R5, mut(L5, 2, 3.0))) === :not_supported
    @test dec(:H3, h3(vcat(R5, 50.0), vcat(L5, 9.0))) === :supported && dec(:H3, h3(R5, L5; nw = 4)) === :not_decidable && dec(:H3, h3(R5, L5; ns = 4)) === :not_decidable
    @test dec(:H3, h3(Dict(zip(P2_X7, R5)), Dict(zip(P2_X7, L5)))) === :supported          # the statistics may come as dictionaries by indicator
    # the count rule keeps h3_count_min for eight indicators (the extension) and for six (a second version)
    @test dec(:H3, h3(vcat(R5, 1.0), vcat(L5, 1.0))) === :supported && dec(:H3, h3(vcat(fill(10.0, 4), 1.0, 1.0), vcat(fill(nextfloat(3.0), 4), 1.0, 1.0))) === :not_supported
    # H4: Delta_phi >= 0.50 and lower end >= 0.25 (both non-strict)
    h4(; n = 5, d = 0.5, lo = 0.25) = (n = n, dphi = d, dphi_lo = lo)
    @test dec(:H4, h4()) === :supported && dec(:H4, h4(d = prevfloat(0.5))) === :not_supported && dec(:H4, h4(lo = prevfloat(0.25))) === :not_supported && dec(:H4, h4(n = 4)) === :not_decidable
    @test dec(:H4, h4(d = nothing)) === :not_decidable && dec(:H4, h4(lo = nothing)) === :not_decidable
    # H5: F_med >= 0.10 and lower end >= 0.05 (both non-strict); an undefined F_X makes it not decidable
    h5(; n = 5, F = 0.1, lo = 0.05) = (n = n, F_med = F, F_med_lo = lo)
    @test dec(:H5, h5()) === :supported && dec(:H5, h5(F = prevfloat(0.1))) === :not_supported && dec(:H5, h5(lo = prevfloat(0.05))) === :not_supported && dec(:H5, h5(n = 4)) === :not_decidable
    @test dec(:H5, h5(F = nothing)) === :not_decidable
    # H6: SUPPORTED iff a decidable part has its maximum over X >= 0.05; NOT SUPPORTED iff both parts are decidable and neither does; otherwise NOT DECIDABLE
    h6(; na = 5, da = 0.05, nb = 5, db = 0.0) = (dir = (n = na, delta = Dict(:LE => da, :SALI => 0.0, :FLI => 0.0, :MEGNO => 0.0)), prec = (n = nb, delta = Dict(X => (X === :E1 ? db : 0.0) for X in P2_X7)))
    @test dec(:H6, h6()) === :supported && dec(:H6, h6(da = prevfloat(0.05))) === :not_supported && dec(:H6, h6(da = 0.0, db = 0.05)) === :supported
    @test dec(:H6, h6(da = 0.0, db = prevfloat(0.05))) === :not_supported && dec(:H6, h6(da = 0.0, db = 0.3, nb = 4)) === :not_decidable           # a part below min_set does not count
    @test dec(:H6, h6(na = 4, da = 0.0, db = 0.0)) === :not_decidable && dec(:H6, h6(na = 4, da = 0.9, db = 0.0)) === :not_decidable
    @test dec(:H6, h6(na = 4, da = 0.9, db = 0.05)) === :supported && dec(:H6, h6(da = 0.0, na = 4, nb = 4)) === :not_decidable
    # the decision of every hypothesis reports its counts and statistics (the report prints them): the result echoes the inputs
    r = AN.p2_decide(:H2, h2(), p)
    @test r.verdict === :supported && r.n == 5 && r.Q == 10.0
    # ROBUST: both verdicts decidable and equal; NOT ROBUST: both decidable and different; NOT COMPARABLE otherwise
    V = (:supported, :not_supported, :not_decidable)
    exp = Dict((a, b) => (a === :not_decidable || b === :not_decidable) ? :NOT_COMPARABLE : a === b ? :ROBUST : :NOT_ROBUST for a in V, b in V)
    @test all(ab -> AN.p2_robust(ab[1], ab[2]) === exp[ab], keys(exp))
end
