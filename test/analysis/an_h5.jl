# AC13 H5 (E15). Spec: analysis-spec.md, Equations E15 and open point P6.

# one orbit of the label table: in_SL, ladder_same (raw labels at the three precisions, before test (a) acts), the final labels
an_h5o(id; sl = false, same = nothing, inst = AN_POS, str = "irregular", excluded = false) =
    (; id, in_SL = sl, ladder_same = same, instability = inst, structure = str, excluded)

# n orbits: nl in S_L of which `same` have equal labels; nu unresolved on either axis; of the `nn` with "none detected", `ns` have unresolved structure
function an_h5set(; n = 100, nl = 20, same = 20, nu = 0, nn = 40, ns = 0)
    os = []
    for i in 1:n
        sl = i <= nl
        sm = sl ? (i <= same) : nothing
        if i <= nn
            push!(os, an_h5o(i; sl, same = sm, inst = AN_NONE, str = i <= ns ? AN_UNR : "quasiperiodic"))
        else
            push!(os, an_h5o(i; sl, same = sm))
        end
    end
    # the unresolved orbits: take from the end (instability "unresolved" and structure "unresolved")
    for k in 0:(nu - 1)
        i = n - k
        os[i] = an_h5o(i; sl = os[i].in_SL, same = os[i].ladder_same, inst = AN_UNR, str = AN_UNR)
    end
    return os
end

@testset "AC13 H5a ladder stability" begin
    p = an_params()
    r = AN.decide_h5(an_h5set(; nl = 20, same = 20), p)
    @test r.h5a.fraction == 1.0 && r.h5a.n == 20 && r.h5a.pass === true
    @test AN.decide_h5(an_h5set(; nl = 20, same = 18), p).h5a.pass === true         # 0.90 exactly passes
    @test AN.decide_h5(an_h5set(; nl = 20, same = 17), p).h5a.pass === false        # 0.85
    r = AN.decide_h5(an_h5set(; nl = 19, same = 19), p)
    @test r.h5a.verdict === :not_decidable && r.h5a.n == 19 && r.h5a.pass === nothing
    # H5a is read from the raw per-precision labels: an orbit whose ladder labels differ counts as different even though its final labels are unresolved
    os = an_h5set(; nl = 20, same = 20)
    os[1] = an_h5o(1; sl = true, same = false, inst = AN_UNR, str = AN_UNR)
    r = AN.decide_h5(os, p)
    @test r.h5a.fraction == 0.95 && r.h5a.pass === true
    @test AN.decide_h5(os, an_params(h5a_min = 0.96)).h5a.pass === false            # not 1.0 by construction
end

@testset "AC13 H5b and H5c" begin
    p = an_params()
    @test AN.decide_h5(an_h5set(; nu = 30), p).h5b.pass === true && AN.decide_h5(an_h5set(; nu = 30), p).h5b.fraction == 0.3      # 0.30 exactly
    @test AN.decide_h5(an_h5set(; nu = 31), p).h5b.pass === false
    @test AN.decide_h5(an_h5set(; nn = 40, ns = 12), p).h5c.pass === true && AN.decide_h5(an_h5set(; nn = 40, ns = 12), p).h5c.fraction == 0.3
    @test AN.decide_h5(an_h5set(; nn = 40, ns = 13), p).h5c.pass === false
    @test AN.decide_h5(an_h5set(; nu = 30), an_params(h5b_max = 0.2999)).h5b.pass === false
    # H5b over the panel without excluded orbits
    os = an_h5set(; nu = 0)
    os[1] = an_h5o(1; sl = true, same = true, excluded = true, inst = AN_UNR, str = AN_UNR)
    r = AN.decide_h5(os, p)
    @test r.h5b.n == 99 && r.h5b.fraction == 0.0
    # H5c denominator: only the orbits whose instability is "none detected"
    r = AN.decide_h5(an_h5set(; nn = 40, ns = 12), p)
    @test r.h5c.n == 40
    # fewer than 20 in the denominator: not decidable with the count
    r = AN.decide_h5(an_h5set(; nn = 19, ns = 0), p)
    @test r.h5c.verdict === :not_decidable && r.h5c.n == 19
end

@testset "AC13 H5 verdict" begin
    p = an_params()
    @test AN.decide_h5(an_h5set(), p).verdict === :supported
    @test AN.decide_h5(an_h5set(; same = 17), p).verdict === :not_supported                   # H5a fails
    @test AN.decide_h5(an_h5set(; nu = 31), p).verdict === :not_supported                     # H5b fails
    @test AN.decide_h5(an_h5set(; ns = 13), p).verdict === :not_supported                     # H5c fails
    # a sub-rule that is not decidable and the others pass: not decidable overall; with a failing one: not supported
    @test AN.decide_h5(an_h5set(; nl = 19, same = 19), p).verdict === :not_decidable
    @test AN.decide_h5(an_h5set(; nl = 19, same = 19, nu = 31), p).verdict === :not_supported
    # the three results are always reported separately
    r = AN.decide_h5(an_h5set(; same = 17, nu = 31), p)
    @test r.h5a.pass === false && r.h5b.pass === false && r.h5c.pass === true
end
