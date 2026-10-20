# AC4 H1 metrics (E4 to E6), AC5 H1 decision (E8), AC6 controls (E9). Spec: analysis-spec.md.

an_lam_all(v = 0.31) = Dict{String, Float64}(m => v for m in vcat(AN_V, "L2"))
const AN_B = ["L1", "L3", "L4", "L4b", "L5", "L5b", "L6"]

@testset "AC4 H1 metrics E4 to E6" begin
    p = an_params()
    r = AN.h1_orbit(an_lam_all(), p)
    @test r.delta_qb == 0 && r.diff_l9 == 0 && r.diff_l2 == 0 && r.pass_a && r.pass_b && r.pass_c
    # prescribed values: B-set spread, L9 and L2 are outside it
    lam = an_lam_all()
    lam["L4"], lam["L3"], lam["L6"] = 0.31 + 3e-9, 0.31 - 2e-9, 0.31 + 1e-9
    lam["L9"], lam["L2"] = 0.31 + 5e-7, 0.31 - 4e-3
    r = AN.h1_orbit(lam, p)
    @test r.delta_qb == maximum(lam[m] for m in AN_B) - minimum(lam[m] for m in AN_B)
    @test r.diff_l9 == abs(lam["L9"] - lam["L5"]) && r.diff_l2 == abs(lam["L2"] - lam["L5"])
    @test r.pass_a && r.pass_b && r.pass_c
    # L9 and L2 are not part of Delta_QB; L9 far out fails only H1b, L2 far out fails only H1c
    lam = an_lam_all(); lam["L9"] = 0.32
    r = AN.h1_orbit(lam, p)
    @test r.delta_qb == 0 && r.pass_a && !r.pass_b && r.pass_c
    lam = an_lam_all(); lam["L2"] = 0.32
    r = AN.h1_orbit(lam, p)
    @test r.pass_a && r.pass_b && !r.pass_c
    # a spread inside B fails H1a only
    lam = an_lam_all(); lam["L5b"] = 0.31 + 2e-8
    r = AN.h1_orbit(lam, p)
    @test !r.pass_a && r.pass_b && r.pass_c
    # an orbit exactly at a bound passes, one unit in the last place above fails (bound set to the computed value)
    lam = an_lam_all(); lam["L4"] = 0.31 + 7e-9
    d = maximum(lam[m] for m in AN_B) - minimum(lam[m] for m in AN_B)
    @test AN.h1_orbit(lam, an_params(h1a_bound = d)).pass_a
    @test !AN.h1_orbit(lam, an_params(h1a_bound = prevfloat(d))).pass_a
    lam = an_lam_all(); lam["L9"] = 0.31 + 6e-7
    d = abs(lam["L9"] - lam["L5"])
    @test AN.h1_orbit(lam, an_params(h1b_bound = d)).pass_b && !AN.h1_orbit(lam, an_params(h1b_bound = prevfloat(d))).pass_b
    lam = an_lam_all(); lam["L2"] = 0.31 + 3e-3
    d = abs(lam["L2"] - lam["L5"])
    @test AN.h1_orbit(lam, an_params(h1c_bound = d)).pass_c && !AN.h1_orbit(lam, an_params(h1c_bound = prevfloat(d))).pass_c
end

an_h1_orbits(n, na, nb, nc) = [(id = i, pass_a = i <= na, pass_b = i <= nb, pass_c = i <= nc) for i in 1:n]

@testset "AC5 H1 decision E8" begin
    p = an_params()
    d = AN.decide_h1(an_h1_orbits(100, 100, 100, 100), p)
    @test d.verdict === :supported && d.n == 100 && d.fractions.a == 1.0 && d.fractions.b == 1.0 && d.fractions.c == 1.0
    # fractions exactly at the limits pass (0.98, 0.98, 0.90), one orbit fewer fails
    @test AN.decide_h1(an_h1_orbits(100, 98, 98, 90), p).verdict === :supported
    @test AN.decide_h1(an_h1_orbits(100, 97, 100, 100), p).verdict === :falsified
    @test AN.decide_h1(an_h1_orbits(100, 100, 97, 100), p).verdict === :falsified
    @test AN.decide_h1(an_h1_orbits(100, 100, 100, 89), p).verdict === :falsified
    # the three sub-rules are independent: each one alone below its limit falsifies, each alone at its limit does not
    for (a, b, c) in ((98, 100, 100), (100, 98, 100), (100, 100, 90))
        @test AN.decide_h1(an_h1_orbits(100, a, b, c), p).verdict === :supported
    end
    r = AN.decide_h1(an_h1_orbits(100, 97, 98, 89), p)
    @test r.verdict === :falsified && r.passes.a == 97 && r.passes.b == 98 && r.passes.c == 89
    # the eligible set below the minimum: not decidable, with its count; at the minimum a verdict
    r = AN.decide_h1(an_h1_orbits(19, 19, 19, 19), p)
    @test r.verdict === :not_decidable && r.n == 19
    @test AN.decide_h1(an_h1_orbits(20, 20, 20, 20), p).verdict === :supported
    @test AN.decide_h1(an_h1_orbits(19, 0, 0, 0), p).verdict === :not_decidable      # neither a pass nor a fail
    # the limits come from the parameter table
    @test AN.decide_h1(an_h1_orbits(100, 90, 100, 100), an_params(h1a_fraction = 0.90)).verdict === :supported
    @test AN.decide_h1(an_h1_orbits(100, 90, 100, 100), p).verdict === :falsified
end

@testset "AC6 exact-exponent controls E9" begin
    p = an_params()
    exact = 0.6931471805599453
    devs(; kw...) = (lam = Dict{String, Float64}(m => exact + 5e-9 for m in AN_V); for (k, v) in kw; lam[string(k)] = exact + v; end; lam)
    # stored: the `control_deviation` attribute of the control's item of each method, as the runner writes it
    row(id, lam; exponent = exact, scope = "bound", stored = :auto) =
        (; id, exponent, scope, lam, stored = stored === :auto ? Dict{String, Float64}(m => abs(lam[m] - exponent) for m in AN_V) : stored)
    # within the bounds (L9 has its own): ok
    r = AN.check_controls([row(871, devs(; L9 = 4e-8))], p)
    @test r.ok && isempty(r.defects) && isempty(r.mismatches)
    # L9 over 5e-8 is a defect naming control and method; L9 exactly at 5e-8 is not
    r = AN.check_controls([row(877, devs(; L9 = 6e-8))], p)
    @test !r.ok && length(r.defects) == 1 && r.defects[1].id == 877 && r.defects[1].method == "L9"
    lam = Dict{String, Float64}(m => exact + 1e-9 for m in AN_V); lam["L9"] = exact + 5e-8
    dev = abs(lam["L9"] - exact)
    @test AN.check_controls([row(877, lam)], an_params(control_bound_L9 = dev)).ok
    @test !AN.check_controls([row(877, lam)], an_params(control_bound_L9 = prevfloat(dev))).ok
    # another method over 1e-8 is a defect (the L9 bound does not apply to it)
    r = AN.check_controls([row(879, devs(; L5 = 1.1e-8))], p)
    @test !r.ok && r.defects[1].method == "L5" && r.defects[1].id == 879
    # 7e-8 on L4b is over its own bound even though L9 could carry it
    @test !AN.check_controls([row(873, devs(; L4b = 7e-8))], p).ok
    # info controls never block, and are listed with their deviations
    r = AN.check_controls([row(878, devs(; L5 = 0.4, L9 = 0.5); scope = "info")], p)
    @test r.ok && isempty(r.defects) && length(r.info) == 1 && r.info[1].id == 878 && r.info[1].deviations["L9"] > 0.4
    # a control with exponent zero is not checked
    r = AN.check_controls([row(872, devs(; L5 = 0.4); exponent = 0.0, scope = "bound")], p)
    @test r.ok && isempty(r.defects)
    # the deviation recomputed from the series against the stored one: more than 1e-12 apart is a mismatch
    lam = devs(; L9 = 4e-8)
    st = Dict{String, Float64}(m => abs(lam[m] - exact) for m in AN_V)
    st1 = merge(st, Dict("L5" => st["L5"] + 5e-13))
    st2 = merge(st, Dict("L5" => st["L5"] + 2e-12))
    @test isempty(AN.check_controls([row(871, lam; stored = st1)], p).mismatches)
    r = AN.check_controls([row(871, lam; stored = st2)], p)
    @test r.mismatches == [871]
    # deviations are absolute
    d = AN.control_deviations(Dict("L5" => exact - 3e-9, "L1" => exact + 2e-9), exact)
    @test isapprox(d["L5"], 3e-9; atol = 1e-15) && isapprox(d["L1"], 2e-9; atol = 1e-15)
end
