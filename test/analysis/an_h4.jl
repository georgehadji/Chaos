# AC11 H4 (E13). Spec: analysis-spec.md, Equations E13. Reference values: benchmarks.md rows with source key SalaArtuso2015
# (0.298 at K = pi/2, 1.172 at K = 2 pi) and Shevchenko2014 (Eq. 25 and 26); the test reads them from the parameter table or the rows.

an_h4(; hL1 = 0.298, hL5 = 0.298, tL1 = 1.172, tL5 = 1.172) = (halfpi = (L1 = hL1, L5 = hL5), twopi = (L1 = tL1, L5 = tL5))

@testset "AC11 H4 decision" begin
    p = an_params()
    r = AN.decide_h4(an_h4(), p)
    @test r.verdict === :supported && r.halfpi.pass && r.twopi.pass
    # inside the tolerance: pass
    @test AN.decide_h4(an_h4(hL1 = 0.298 + 0.009, hL5 = 0.298 - 0.009, tL1 = 1.172 + 0.019), p).verdict === :supported
    # exactly at the tolerance passes (tolerance set to the computed distance), one unit below fails
    x = 0.298 + 0.01
    d = abs(x - 0.298)
    @test AN.decide_h4(an_h4(hL1 = x), an_params(h4_tol_halfpi = d)).verdict === :supported
    @test AN.decide_h4(an_h4(hL1 = x), an_params(h4_tol_halfpi = prevfloat(d))).verdict === :did_not_reproduce
    # 1e-9 outside: did not reproduce, with the numbers
    r = AN.decide_h4(an_h4(hL1 = 0.298 + 0.01 + 1e-9), p)
    @test r.verdict === :did_not_reproduce && !r.halfpi.pass && r.twopi.pass
    @test r.halfpi.L1 ≈ 0.308 + 1e-9 atol = 1e-15          # the numbers are reported
    # both methods are required at each K, and both K
    @test AN.decide_h4(an_h4(hL5 = 0.4), p).verdict === :did_not_reproduce
    @test AN.decide_h4(an_h4(tL1 = 1.2), p).verdict === :did_not_reproduce
    @test AN.decide_h4(an_h4(tL5 = 1.2), p).verdict === :did_not_reproduce
    # references and tolerances come from the parameter table
    @test AN.decide_h4(an_h4(tL1 = 1.3, tL5 = 1.3), an_params(h4_ref_2pi = 1.3)).verdict === :supported
end

@testset "AC11 comparators without a verdict" begin
    # Shevchenko2014 fits by branch: Eq. (25) for K < 1, Eq. (26) for 1.1 <= K < 4.4 and for K >= 4.4; no fit in between
    @test AN.shevchenko_fit(0.5) ≈ 0.1333 * 0.5
    @test AN.shevchenko_fit(0.971635406) ≈ 0.1333 * 0.971635406
    @test AN.shevchenko_fit(pi / 2) ≈ 0.469 * sqrt(pi / 2 - 1.037)
    @test AN.shevchenko_fit(4.399) ≈ 0.469 * sqrt(4.399 - 1.037)
    @test AN.shevchenko_fit(4.4) ≈ log(4.4 / 2) + 1 / 4.4^2
    @test AN.shevchenko_fit(2 * pi) ≈ log(pi) + 1 / (2 * pi)^2
    @test AN.shevchenko_fit(10.0) ≈ log(5.0) + 0.01
    @test AN.shevchenko_fit(1.05) === nothing && AN.shevchenko_fit(1.0999) === nothing        # between the two branches
    # Chirikov line
    @test AN.chirikov_line(10.0) == log(5.0) && AN.chirikov_line(2 * pi) == log(pi)
end
