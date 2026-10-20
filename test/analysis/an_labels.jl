# AC12 labels (E14): instability and structure axes, the four disagreement tests, the exclusion rule. Spec: analysis-spec.md, E14.

const AN_POS, AN_NONE, AN_UNR = "positive growth supported", "none detected on the tested horizon", "unresolved"
an_ld(; lam_N = 0.3, lam_N10 = 0.3, rho1 = 0.1, rho2 = 0.1, periodic = false, residual = nothing, estimator_fail = false, diag_exceeded = false,
      excluded = false, ladder = nothing, direction = nothing) =
    (; lam_N, lam_N10, rho1, rho2, periodic, residual, estimator_fail, diag_exceeded, excluded, ladder, direction)
an_lad(; f64 = 0.3, d64 = 0.3, bf = 0.3, lf = (AN_POS, "irregular"), ld = (AN_POS, "irregular"), lb = (AN_POS, "irregular")) =
    (; lam_f64 = f64, lam_d64 = d64, lam_bf = bf, labels = (f64 = lf, d64 = ld, bf = lb))

@testset "AC12 growth ratio r(N)" begin
    @test AN.growth_ratio(0.3, 0.3) == 10.0
    @test AN.growth_ratio(0.5, 1.0) == 5.0
    @test AN.growth_ratio(1.0, 0.0) === nothing && AN.growth_ratio(1.0, -0.1) === nothing
end

@testset "AC12 instability axis" begin
    p = an_params()
    base(; kw...) = AN.base_label(an_ld(; kw...), p)
    @test base().instability == AN_POS && base().structure == "irregular"
    # lambda at the threshold: positive (>=); one unit below with r >= 5 is neither rule: unresolved
    @test base(lam_N = 1e-4, lam_N10 = 1e-4).instability == AN_POS
    b = base(lam_N = prevfloat(1e-4), lam_N10 = prevfloat(1e-4))
    @test b.instability == AN_UNR && b.structure == AN_UNR
    # growth ratio: exactly 5 positive, just below 5 with lambda >= threshold unresolved
    @test base(lam_N = 0.5, lam_N10 = 1.0).instability == AN_POS
    @test base(lam_N = 0.3, lam_N10 = 10 * 0.3 / 4.999).instability == AN_UNR
    # lambda below the threshold and r below 5: none detected
    @test base(lam_N = 5e-5, lam_N10 = 2e-4).instability == AN_NONE
    # r undefined: unresolved
    @test base(lam_N10 = 0.0).instability == AN_UNR && base(lam_N10 = -1e-3).instability == AN_UNR
    # thresholds from the parameter table
    @test AN.base_label(an_ld(lam_N = 5e-5, lam_N10 = 5e-5), an_params(lambda_threshold = 1e-5)).instability == AN_POS
    @test AN.base_label(an_ld(lam_N = 0.5, lam_N10 = 1.0), an_params(growth_ratio = 5.0000001)).instability == AN_UNR
end

@testset "AC12 structure axis" begin
    p = an_params()
    base(; kw...) = AN.base_label(an_ld(; kw...), p)
    none = (lam_N = 5e-5, lam_N10 = 2e-4)
    # periodic: only an exact periodic point of the panel, with the return residual at most 1e-10
    @test base(; periodic = true, residual = 1e-10).structure == "periodic"
    @test base(; periodic = true, residual = 1.01e-10).structure == "irregular"           # falls through to the other rules
    @test base(; periodic = true, residual = nothing).structure == "irregular"
    @test base(; periodic = false, residual = 0.0).structure == "irregular"               # not a panel periodic point
    @test base(; none..., periodic = true, residual = 0.0).structure == "periodic"        # periodic comes before quasiperiodic
    @test base(; periodic = true, residual = 0.0).instability == AN_POS                   # the instability axis is separate
    # quasiperiodic: none detected and the two half rotation numbers agree within 1e-7
    @test base(; none..., rho1 = 0.1, rho2 = 0.1).structure == "quasiperiodic"
    @test base(; none..., rho1 = 0.1, rho2 = 0.1 + 1.01e-7).structure == AN_UNR
    d = abs((0.1 + 9e-8) - 0.1)
    @test AN.base_label(an_ld(; none..., rho1 = 0.1, rho2 = 0.1 + 9e-8), an_params(rotation_tol = d)).structure == "quasiperiodic"
    @test AN.base_label(an_ld(; none..., rho1 = 0.1, rho2 = 0.1 + 9e-8), an_params(rotation_tol = prevfloat(d))).structure == AN_UNR
    @test base(; none..., rho1 = nothing, rho2 = nothing).structure == AN_UNR
    # irregular only with positive growth; with instability unresolved the structure is unresolved
    @test base().structure == "irregular"
    @test base(lam_N = 0.3, lam_N10 = 10 * 0.3 / 4.999).structure == AN_UNR
    @test base(; none..., rho1 = 0.1, rho2 = 0.3).structure == AN_UNR                    # none detected but the rotation numbers differ
end

@testset "AC12 disagreement tests force both labels" begin
    p = an_params()
    lab(; kw...) = AN.orbit_label(an_ld(; kw...), p)
    r = lab()
    @test r.instability == AN_POS && r.structure == "irregular" && isempty(r.forced) && r.raw == (instability = AN_POS, structure = "irregular")
    r = lab(estimator_fail = true)
    @test (r.instability, r.structure) == (AN_UNR, AN_UNR) && r.forced == [:estimator]
    @test r.raw.instability == AN_POS                                    # the raw label stays available
    r = lab(diag_exceeded = true)
    @test (r.instability, r.structure) == (AN_UNR, AN_UNR) && r.forced == [:diagnostics]
    r = lab(excluded = true)
    @test (r.instability, r.structure) == (AN_UNR, AN_UNR) && r.forced == [:excluded]
    # a test that does not apply (no ladder data, no direction data) never forces
    @test lab(ladder = nothing, direction = nothing).forced == Symbol[]
    # ladder: lambda difference above delta_lad at N_B (either pair), or labels that differ between precisions
    r = lab(ladder = an_lad(bf = 0.3 + 0.01))
    @test r.forced == [:ladder] && (r.instability, r.structure) == (AN_UNR, AN_UNR) && r.ladder_same === true
    @test lab(ladder = an_lad(d64 = 0.3 + 0.01, bf = 0.3 + 0.01)).forced == [:ladder]                 # Float64 against Double64
    @test lab(ladder = an_lad(f64 = 0.3 + 0.004, d64 = 0.3 + 0.004, bf = 0.3)).forced == Symbol[]      # 4e-3, inside
    r = lab(ladder = an_lad(lb = ("none detected on the tested horizon", "quasiperiodic")))
    @test r.forced == [:ladder] && r.ladder_same === false
    r = lab(ladder = an_lad(lb = (AN_POS, "unresolved")))                                              # only the structure axis differs
    @test r.forced == [:ladder] && r.ladder_same === false
    @test lab(ladder = an_lad()).ladder_same === true && lab().ladder_same === nothing
    # the ladder tolerance exactly: not forced; one unit below: forced
    tol = abs((0.3 + 5e-3) - 0.3)
    ld = an_lad(bf = 0.3 + 5e-3)
    @test AN.orbit_label(an_ld(ladder = ld), an_params(ladder_tol = tol)).forced == Symbol[]
    @test AN.orbit_label(an_ld(ladder = ld), an_params(ladder_tol = prevfloat(tol))).forced == [:ladder]
    # direction: the five labels all equal, or not
    same = [(AN_POS, "irregular") for _ in 1:5]
    @test lab(direction = same).forced == Symbol[]
    r = lab(direction = vcat(same[1:4], [(AN_UNR, AN_UNR)]))
    @test r.forced == [:direction] && (r.instability, r.structure) == (AN_UNR, AN_UNR)
    @test lab(direction = vcat(same[1:4], [(AN_POS, "unresolved")])).forced == [:direction]            # a change on either axis
    # two and all tests together
    r = lab(estimator_fail = true, diag_exceeded = true)
    @test Set(r.forced) == Set([:estimator, :diagnostics]) && (r.instability, r.structure) == (AN_UNR, AN_UNR)
    r = lab(estimator_fail = true, diag_exceeded = true, excluded = true, ladder = an_lad(bf = 0.4), direction = vcat(same[1:4], [(AN_UNR, AN_UNR)]))
    @test Set(r.forced) == Set([:ladder, :estimator, :direction, :diagnostics, :excluded])
    # a periodic orbit is forced too (either axis)
    r = lab(periodic = true, residual = 0.0, diag_exceeded = true)
    @test (r.instability, r.structure) == (AN_UNR, AN_UNR) && r.raw.structure == "periodic"
end
