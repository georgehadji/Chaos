# AC1 to AC3: coefficients (sum conditions, linear order from the stored tables, table identities). Every testset starts by naming the new symbols it uses,
# so a missing module errors that testset (UndefVarError). Tables are read through `coefficients(m)` (drift weights c, kick weights d, drift first in each stage).
# AC3 (c) (package oracle O6) and (d) (REBOUND, O8) and the print comparisons (b2) are stage 04 (verify environment), see the test plan.

const IC_BIB = joinpath(@__DIR__, "..", "..", "_config", "bibliography")
ic_bibkeys() = Set(m.captures[1] for m in eachmatch(r"^@\w+\{([^,\s]+),"m, read(joinpath(IC_BIB, "references.bib"), String)))
# kicks of a table as (position of the kick = cumulative drift before it, weight), zero kicks dropped: invariant under merging of drifts
function ic_kickpos(c, d)
    cs = cumsum(c)
    idx = [i for i in eachindex(d) if !iszero(d[i])]
    return [cs[i] for i in idx], [d[i] for i in idx]
end

@testset "AC1 sum conditions, palindromic symmetry, digits, sources" begin
    _ = (coefficients, stage_count, is_palindromic, table_digits, method_source)
    keys = ic_bibkeys()
    unverified = read(joinpath(IC_BIB, "unverified.md"), String)
    setprecision(BigFloat, IC_BITS) do
        for M in IC_EXPL
            m = M.mk(BigFloat)
            c, d = coefficients(m)
            k = length(c)
            @test length(d) == k && k >= 1 && stage_count(m) == k
            et = max(ic_eps(BigFloat), 10.0^(-M.dtab))
            @test abs(sum(c) - 1) <= 4k * et
            @test abs(sum(d) - 1) <= 4k * et
            @test table_digits(m) == M.dtab
            if M.sym === true
                @test is_palindromic(m)
                @test all(c[i] == c[k + 1 - i] for i in 1:k)           # bitwise after construction
                @test all(d[i] == d[k - i] for i in 1:(k - 1)) && d[k] == 0
            elseif M.sym === false
                @test !is_palindromic(m)
            elseif is_palindromic(m)                                    # S9: not asserted either way, but a declared table must be symmetric
                @test all(c[i] == c[k + 1 - i] for i in 1:k) && all(d[i] == d[k - i] for i in 1:(k - 1))
            end
            # sources are keys of references.bib and none is in unverified.md (never cite from memory)
            src = strip.(split(method_source(m), ','))
            @test !isempty(src) && all(s -> s in keys, src) && !any(s -> occursin(s, unverified), src)
        end
    end
    for T in (Float64, Double64), M in IC_EXPL                          # the same conditions in the lower precisions
        c, d = coefficients(M.mk(T))
        k = length(c)
        et = max(ic_eps(T), 10.0^(-M.dtab))
        @test abs(Float64(sum(ic_big.(c)) - 1)) <= 4k * et && abs(Float64(sum(ic_big.(d)) - 1)) <= 4k * et
    end
end

@testset "AC2 linear order of every explicit table" begin
    _ = (coefficients, nominal_order)
    for M in IC_EXPL
        setprecision(BigFloat, IC_BITS) do
            m = M.mk(BigFloat)
            c, d = coefficients(m)
            p = M.p
            et = max(ic_eps(BigFloat), 10.0^(-M.dtab))
            S = Float64(1 + sum(abs, c) + sum(abs, d))
            Mp = ic_steppoly(c, d, p + 1); R = ic_rotpoly(p + 1)
            dif = [Float64(maximum(abs, Mp[j + 1] - R[j + 1])) for j in 0:(p + 1)]
            @test all(dif[j + 1] <= 64 * et * S^p for j in 0:p)                       # Taylor coefficients through the nominal order match the rotation
            @test nominal_order(m) == p
            dif[p + 2] > 1e3 * et || @info "AC2 reported: linear order higher than nominal (known property of some tables, not a failure)" M.id p dif[p + 2]
        end
    end
end

@testset "AC3a triple jump and Suzuki builders reproduce the closed forms" begin
    _ = (triple_jump, suzuki_jump, Leapfrog, ForestRuth, coefficients)
    setprecision(BigFloat, IC_BITS) do
        th = 1 / (2 - BigFloat(2)^(BigFloat(1) / 3))                                 # theta of Omelyan2002 Eq. (8) (spec item 8)
        ep = ic_eps(BigFloat)
        # 4th order from leapfrog and S5: kick weights (theta, 1 - 2 theta, theta) at drift positions (theta/2, 1/2, 1 - theta/2) (Yoshida1990 Eq. (2.11), derived)
        for m in (triple_jump(Leapfrog(BigFloat), 4), ForestRuth(BigFloat))
            pos, w = ic_kickpos(coefficients(m)...)
            @test length(w) == 3
            @test maximum(abs, w .- [th, 1 - 2th, th]) <= 8 * ep
            @test maximum(abs, pos .- [th / 2, BigFloat(1) / 2, 1 - th / 2]) <= 8 * ep
        end
        # all orders: the F1 step polynomial equals the composition of leapfrog pieces with the triple jump / Suzuki fractions (independent construction)
        deg = 12
        for (builder, wfun, counts) in ((triple_jump, ic_tj_weights, (3, 9, 27)), (suzuki_jump, ic_su_weights, (5, 25, 125)))
            for (order, cnt) in zip((4, 6, 8), counts)
                m = builder(Leapfrog(BigFloat), order)
                c, d = coefficients(m)
                @test count(!iszero, d) == cnt                                          # leapfrog pieces = kicks (Algorithm item 2)
                P = ic_steppoly(c, d, deg); Q = ic_composepoly(wfun(order), deg)
                @test maximum(Float64(maximum(abs, P[j] - Q[j])) for j in 1:(deg + 1)) <= 1e-60
                @test nominal_order(m) == order
            end
        end
        # the builder runs in the element type: Float64 table within 8 eps of the 256-bit one (kick positions and weights)
        pb, wb = ic_kickpos(coefficients(triple_jump(Leapfrog(BigFloat), 8))...)
        pf, wf = ic_kickpos(coefficients(triple_jump(Leapfrog(Float64), 8))...)
        @test length(wf) == length(wb) && maximum(abs, wf .- wb) <= 64 * ic_eps(Float64) && maximum(abs, pf .- pb) <= 64 * ic_eps(Float64)
    end
end

@testset "AC3b Yoshida solution A against Kahan-Li s7odr6, stage counts" begin
    _ = (Yoshida, KahanLi, stage_count, coefficients)
    setprecision(BigFloat, IC_BITS) do
        pa, wa = ic_kickpos(coefficients(Yoshida(BigFloat, 6, :A))...)
        pb, wb = ic_kickpos(coefficients(KahanLi(BigFloat, :s7odr6))...)
        @test length(wa) == 7 && length(wb) == 7
        @test maximum(abs, wa .- wb) <= 1e-14 && maximum(abs, pa .- pb) <= 1e-14             # independent prints of the same method (15 digits)
        @test stage_count(Yoshida(BigFloat, 6, :A)) == 8 && stage_count(Yoshida(BigFloat, 8, :D)) == 16 && stage_count(Yoshida(BigFloat, 8, :A)) == 16   # Yoshida1990 text after Table 2
        @test count(!iszero, coefficients(KahanLi(BigFloat, 6))[2]) == 9 && count(!iszero, coefficients(KahanLi(BigFloat, 8))[2]) == 17   # s9odr6a, s17odr8a
    end
end

@testset "AC3b3 S3: z of the printed formula, weights, error ratio against leapfrog" begin
    _ = (McLachlan2, Leapfrog, coefficients, advance_n)
    setprecision(BigFloat, IC_BITS) do
        y = cbrt(2 * sqrt(BigFloat(326)) - 36)
        z = (y^2 + 6y - 2) / (12y)
        @test abs(Float64(z) - 0.1932) <= 5e-5                                              # McLachlan1995 Table 2 print (spec item 8c)
        c, d = coefficients(McLachlan2(BigFloat))
        pos, w = ic_kickpos(c, d)
        @test abs(sum(c) - 1) <= 8 * ic_eps(BigFloat) && length(w) == 2
        @test maximum(abs, w .- [BigFloat(1) / 2, BigFloat(1) / 2]) <= 8 * ic_eps(BigFloat)
        @test maximum(abs, pos .- [z, 1 - z]) <= 64 * ic_eps(BigFloat)                      # Eq. (4.2): kicks after drifts z and 1 - z
    end
    # effective error ratio at equal work on the pendulum: S3 at tau = 0.1 against leapfrog at tau/2 (T = 8)
    z0 = IC_START[:f4]; ref = ic_ref(:f4, z0, 8); flow = ic_flow(Float64, :f4)
    e3 = ic_err(advance_n(McLachlan2(Float64), flow, ic_to(Float64, z0), 0.1, 80), ref)
    e2 = ic_err(advance_n(Leapfrog(Float64, :DKD), flow, ic_to(Float64, z0), 0.05, 160), ref)
    @info "AC3 reported: S3 / leapfrog error ratio at equal work" e3 e2 e3 / e2
    @test e3 / e2 < 0.6
end

@testset "AC3b4 S13: Eq. (17) of SofroniouSpaletta2005" begin
    _ = (SofroniouSpaletta10, coefficients, table_digits, advance)
    setprecision(BigFloat, IC_BITS) do
        m = SofroniouSpaletta10(BigFloat)
        c, d = coefficients(m)
        @test table_digits(m) == 32
        @test count(!iszero, d) == 35                                                         # 35 stages
        @test abs(Float64(sum(abs, filter(!iszero, d))) - 6.18133) <= 1e-5                    # Table 1 print of the norm (spec item 9a)
        et = 1e-32
        @test abs(sum(c) - 1) <= 4 * length(c) * et && abs(sum(d) - 1) <= 4 * length(c) * et
        # single-step error on the harmonic oscillator (omega = 1.3 is not needed: omega = 1), h = 1/2 .. 1/16: local order 11
        flow = ic_flow(BigFloat, :f1_1); z0 = ic_to(BigFloat, IC_START[:f1_1])
        errs = [ic_err(advance(m, flow, z0, BigFloat(2)^(-j)), ic_f1_exact(:f1_1, IC_START[:f1_1], BigFloat(2)^(-j))) for j in 1:4]
        sl = log2(errs[1] / errs[4]) / 3
        @info "AC3b4 reported: single-step slope" sl errs
        @test abs(sl - 11) <= 0.3
    end
end

@testset "AC3b5 S9: Calvo-Sanz-Serna table, RKN form against the splitting form" begin
    _ = (CalvoSanzSerna4, coefficients, table_digits, advance)
    m = CalvoSanzSerna4(Float64)
    c, d = coefficients(m)
    @test table_digits(m) == 18 && length(c) == 5 && c[1] == 0                                # kick first: gamma_1 = 0 (CalvoSanzSerna1993a Eq. (2.4a))
    et = max(ic_eps(Float64), 1e-18)
    gam = cumsum(c)
    @test abs(sum(d) - 1) <= 4 * et * 5 && abs(gam[5] - 1) <= 4 * et * 5                     # sum b_i = 1, gamma_s = 1 (Eq. (2.4b))
    # RKN form, Eqs. (2.2), (2.3), (2.5): alpha_ij = b_j (gamma_i - gamma_j), beta_i = b_i (1 - gamma_i), Q_i = q + h gamma_i p + h^2 sum_{j<i} alpha_ij f(Q_j)
    function rkn(kind, z, h; eps = 1.0)
        n = length(z) ÷ 2; q = z[1:n]; p = z[(n + 1):end]
        f(Q) = -ic_gradV(kind, Q; eps = eps)
        F = Vector{Vector{Float64}}(undef, 5); s = 5
        for i in 1:s
            Q = q .+ h * gam[i] .* p
            for j in 1:(i - 1)
                Q = Q .+ h^2 * d[j] * (gam[i] - gam[j]) .* F[j]
            end
            F[i] = f(Q)
        end
        pn = p .+ h .* sum(d[i] .* F[i] for i in 1:s)
        qn = q .+ h .* p .+ h^2 .* sum(d[i] * (1 - gam[i]) .* F[i] for i in 1:s)
        return vcat(qn, pn)
    end
    rng = StableRNG(ic_seed(3, 1))
    for kind in (:f1_3, :f4), _ in 1:20
        z0 = ic_draw(rng, kind); flow = ic_flow(Float64, kind)
        a = collect(advance(m, flow, ic_to(Float64, z0), 0.1))
        b = rkn(kind, z0, 0.1)
        @test maximum(abs, a .- b) <= 16 * ic_eps(Float64) * ic_mmax(z0)
    end
end
