# AC5 fixed points (E5), AC6 parabolic band (E4), AC7 period-2 families (E6, oracle O1), AC8 no spurious periodic point (I7),
# AC9 deduplication (E9). Every expected value is an analytic closed form evaluated at 1024 bits (spec, Oracle O1).

# roots x in [0, 2 pi) of g(x) = 2 x + pi c + (K/2) sin x - 2 pi k, by sign change on a grid and bisection at the current precision
function ex_family_c_points(K, c)
    pib = BigFloat(pi)
    pts = Tuple{BigFloat, BigFloat}[]
    for k in -1:3
        g(x) = 2x + pib * c + K / 2 * sin(x) - 2pib * k
        xs = [2pib * i / 4000 for i in 0:4000]
        for i in 1:4000
            a, b = xs[i], xs[i + 1]
            g(a) == 0 && push!(pts, (a, pib * c - K / 2 * sin(a)))
            if g(a) * g(b) < 0
                for _ in 1:1100
                    m = (a + b) / 2
                    g(a) * g(m) <= 0 ? (b = m) : (a = m)
                end
                x0 = (a + b) / 2
                push!(pts, (x0, pib * c - K / 2 * sin(x0)))
            end
        end
    end
    # keep genuine period-2 points: the point is not fixed under T
    return filter(z -> ex_tdist(z, ex_T(K, z)) > 1e-6, unique(pts))
end

@testset "AC5 fixed points" begin
    # K_c = 0.971635406: benchmarks.md row "Golden-mean KAM torus breakup" (also the study's third K value)
    for K in (0.5, 0.971635406, pi / 2, 3.9, 4.0, 4.1, 2pi, 10.0)
        ps = ex_pset(K; q_max = 1, N_g = 16)
        setprecision(BigFloat, EX_BITS) do
            Kb = BigFloat(K)
            for (z, tr_exp) in (((zero(BigFloat), zero(BigFloat)), 2 + Kb), ((BigFloat(pi), zero(BigFloat)), 2 - Kb))
                os = [o for o in ps.orbits if o.q == 1 && any(pt -> ex_tdist(ex_w.(pt), z) < 1e-12, o.points)]
                @test length(os) == 1                                  # listed once after dedup
                isempty(os) && continue
                o = os[1]
                @test abs(ex_w(o.trace) - tr_exp) <= 64 * EX_EPS_BF * max(1, abs(tr_exp))
                cls = z[1] == 0 ? :hyperbolic : (K < 4 ? :elliptic : (K == 4 ? :parabolic : :hyperbolic))
                @test o.class === cls
                want = abs(tr_exp) > 2 && cls === :hyperbolic ? log((abs(tr_exp) + sqrt(tr_exp^2 - 4)) / 2) : zero(BigFloat)
                @test abs(ex_w(o.exponent) - want) <= 1e-60 * max(1, abs(want))
            end
        end
    end
end

@testset "AC6 parabolic band and class from any point of an orbit" begin
    rho = 1e-3
    t = EX.classify_trace
    for s in (1, -1)
        @test t(s * (2 - 2rho), rho) === :elliptic
        @test t(s * (2 - rho / 2), rho) === :parabolic
        @test t(s * 2.0, rho) === :parabolic
        @test t(s * (2 + rho / 2), rho) === :parabolic
        @test t(s * (2 + 2rho), rho) === :hyperbolic
    end
    @test t(0.0, rho) === :elliptic
    # E4: exponent per map step ln|mu_+| / q for |t| > 2 (|mu_+| = (|t| + sqrt(t^2 - 4)) / 2), 0 otherwise
    @test abs(EX.orbit_exponent(3.0, 2) - log((3 + sqrt(5.0)) / 2) / 2) <= 4 * EX_EPS
    @test abs(EX.orbit_exponent(-58.52, 2) - log((58.52 + sqrt(58.52^2 - 4)) / 2) / 2) <= 4 * EX_EPS
    @test EX.orbit_exponent(1.0, 3) == 0 && EX.orbit_exponent(-2.0, 1) == 0
    # the trace, hence the class, is the same at every point of an orbit (tr(AB) = tr(BA))
    ps = ex_pset(5.0)
    setprecision(BigFloat, 256) do
        for o in ps.orbits, i in 1:length(o.points)
            z = (o.points[i][1], o.points[i][2])
            @test abs(ex_orbit_trace(BigFloat(5), z, o.q) - o.trace) <= 1e-15 * max(1, abs(o.trace))
            @test EX.classify_trace(ex_orbit_trace(BigFloat(5), z, o.q), EX_RHO) === o.class
        end
    end
end

@testset "AC7 period-2 families" begin
    for K in (0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 5.0, 7.0, 10.0)
        ps = ex_pset(K)
        setprecision(BigFloat, EX_BITS) do
            Kb, pib = BigFloat(K), BigFloat(pi)
            tol = 1e4 * EX_RHO^2
            # family A: {(0, pi), (pi, pi)}, trace 2 - K^2
            oa = ex_orbit_with(ps, 2, (zero(BigFloat), pib))
            @test oa !== nothing
            if oa !== nothing
                @test any(pt -> ex_tdist(ex_w.(pt), (pib, pib)) <= tol, oa.points)
                tr = 2 - Kb^2
                @test abs(Float64(oa.trace) - Float64(tr)) <= 1e-16 * max(1, abs(Float64(tr)))
                @test oa.class === (K < 2 ? :elliptic : (K == 2 ? :parabolic : :hyperbolic))
            end
            # family B for K > 2 pi: sin x_a = 2 pi / K, trace 2 + 4 pi^2 - K^2
            if K > 2pi
                xa = asin(2pib / Kb)
                for (z0, z1) in (((xa, pib), (xa + pib, pib)), ((pib - xa, pib), (2pib - xa, pib)))
                    ob = ex_orbit_with(ps, 2, (ex_wrap(z0[1]), pib))
                    @test ob !== nothing
                    ob === nothing && continue
                    @test any(pt -> ex_tdist(ex_w.(pt), z1) <= tol, ob.points)
                    tr = 2 + 4pib^2 - Kb^2
                    @test abs(Float64(ob.trace) - Float64(tr)) <= 1e-16 * max(1, abs(Float64(tr)))
                end
            end
            # family C: every root of g (both c) is an orbit found, trace c0^2 + 4 c0 + 2
            if K in (5.0, 7.0, 10.0)
                for c in (0, 1)
                    roots = ex_family_c_points(Kb, c)
                    @test !isempty(roots)
                    for z in roots
                        oc = ex_orbit_with(ps, 2, z)
                        @test oc !== nothing
                        oc === nothing && continue
                        c0 = Kb * cos(z[1])
                        tr = c0^2 + 4c0 + 2
                        @test abs(Float64(oc.trace) - Float64(tr)) <= 1e-16 * max(1, abs(Float64(tr)))
                    end
                end
            end
            # exponent of every period-2 orbit: ln|mu_+| / 2 when hyperbolic, exactly 0 otherwise (E4)
            for o in ps.orbits
                o.q == 2 || continue
                t = Float64(o.trace)
                want = o.class === :hyperbolic ? log((abs(t) + sqrt(t^2 - 4)) / 2) / 2 : 0.0
                o.class === :hyperbolic ? (@test abs(Float64(o.exponent) - want) <= 1e-15 * want) : (@test o.exponent == 0)
            end
        end
    end
    # K = 5: the family C orbit with c = 0 is elliptic and contains the island centre O1 = (0.68002, 0.36002) in unit-torus
    # coordinates (benchmarks.md row "Island centre and unstable 4/9 orbit point", ContopoulosHarsoula2008 Sec. 2 p.4), 1e-4 per unit coordinate
    ps = ex_pset(5.0)
    O1 = (EX_TWO_PI * 0.68002, EX_TWO_PI * 0.36002)
    near(pt) = ex_delta(Float64(pt[1]) - O1[1]) <= EX_TWO_PI * 1e-4 && ex_delta(Float64(pt[2]) - O1[2]) <= EX_TWO_PI * 1e-4
    @test any(o -> o.q == 2 && o.class === :elliptic && any(near, o.points), ps.orbits)
end

@testset "AC8 no spurious periodic point (I7)" begin
    for K in (1.5, 5.0)
        ps = ex_pset(K)
        @test !isempty(ps.orbits)
        setprecision(BigFloat, EX_BITS) do
            Kb = BigFloat(K)
            for o in ps.orbits
                q = o.q
                for z in o.points
                    zb = (ex_w(z[1]), ex_w(z[2]))
                    @test ex_residual(Kb, zb, q) <= EX_RHO * (1 + ldexp(1.0, -200))
                    for q2 in 1:(q - 1)                             # primitive: no proper divisor of q returns the point
                        q % q2 == 0 && @test ex_residual(Kb, zb, q2) > EX_RHO
                    end
                end
                z = (ex_w(o.points[1][1]), ex_w(o.points[1][2]))   # the stored representative carries the windings
                zq = ex_Tq(Kb, z, q)
                tp = 2 * BigFloat(pi)
                @test abs(zq[1] - z[1] - tp * o.m) <= EX_RHO && abs(zq[2] - z[2] - tp * o.j) <= EX_RHO
                @test o.residual <= EX_RHO
            end
        end
        setprecision(BigFloat, 256) do
            for o in ps.orbits
                z = (BigFloat(o.points[1][1]), BigFloat(o.points[1][2]))
                J = ex_J(BigFloat(K), z, o.q)
                @test abs(J[1, 1] * J[2, 2] - J[1, 2] * J[2, 1] - 1) <= 16 * o.q * EX_EPS_BF * max(1, sum(abs2, J))
            end
        end
    end
    # acceptance rule on synthetic residuals; a Newton step does not exist at tr J = 2 and is not attempted
    rho = 1e-6
    @test EX.accept_residual(0.5rho, rho)
    @test !EX.accept_residual(1.5rho, rho)
    @test EX.accept_residual(rho, rho)
    @test !EX.accept_residual(NaN, rho) && !EX.accept_residual(Inf, rho)
    @test EX.newton_correction([1.0 1.0; 0.0 1.0], [1.0, 1.0]) === nothing           # tr = 2, det(J - I) = 0
    J = [2.0 1.0; 1.0 1.0]; F = [0.1, -0.2]                                           # det J = 1, tr J = 3: the E3 closed form applies
    d = EX.newton_correction(J, F)                                                    # Delta = -(J - I)^-1 F
    @test d !== nothing && maximum(abs, (J - I) * collect(d) + F) <= 1e-14
end

@testset "AC9 deduplication" begin
    rho, dup = 1e-6, 1e-8
    sol(q, z, id) = (q = q, z = z, id = id)
    setprecision(BigFloat, 256) do
        pib = BigFloat(pi)
        K = 1.5
        a, b = (zero(BigFloat), pib), (pib, pib)                                       # family A at K = 1.5
        rng = StableRNG(ex_seed(9))
        reps = map(1:5) do _
            sols = [sol(2, a, 4), sol(2, b, 9)]
            sols = sols[sortperm(rand(rng, 2))]
            os = EX.collect_orbits(K, sols; rho = rho, delta_dup = dup)
            @test length(os) == 1 && os[1].q == 2
            os[1].points[1]
        end
        @test all(==(reps[1]), reps)                                                   # same stored representative in any order
        @test Float64(reps[1][1]) == 0.0 && Float64(reps[1][2]) == Float64(pib)        # lexicographically smallest wrapped point
        # (b) a q = 4 solution at the period-2 orbit is stored with period 2 and adds no orbit
        os = EX.collect_orbits(K, [sol(2, a, 1), sol(4, a, 2)]; rho = rho, delta_dup = dup)
        @test length(os) == 1 && os[1].q == 2
        os = EX.collect_orbits(K, [sol(4, b, 2)]; rho = rho, delta_dup = dup)
        @test length(os) == 1 && os[1].q == 2
        # (c) two distinct period-1 orbits near the fixed point (0, 0) at 10.5 delta_dup are kept, at delta_dup / 10 merged
        # (exactly 10 delta_dup sits on the refusal edge of E9; see the plan, ambiguity A3)
        z0 = (zero(BigFloat), zero(BigFloat))
        os = EX.collect_orbits(0.5, [sol(1, z0, -1), sol(1, (BigFloat(10.5dup), zero(BigFloat)), 3)]; rho = rho, delta_dup = dup)
        @test length(os) == 2
        os = EX.collect_orbits(0.5, [sol(1, z0, -1), sol(1, (BigFloat(dup / 10), zero(BigFloat)), 3)]; rho = rho, delta_dup = dup)
        @test length(os) == 1
        # (d) smallest grid id wins; an exact fixed point (id -1) wins over any Newton solution
        os = EX.collect_orbits(0.5, [sol(1, (BigFloat(1e-9), zero(BigFloat)), 7), sol(1, (BigFloat(2e-9), zero(BigFloat)), 3)];
                               rho = rho, delta_dup = dup)
        @test length(os) == 1 && Float64(os[1].points[1][1]) == 2e-9
        os = EX.collect_orbits(0.5, [sol(1, (BigFloat(1e-9), zero(BigFloat)), 1), sol(1, z0, -1)]; rho = rho, delta_dup = dup)
        @test length(os) == 1 && Float64(os[1].points[1][1]) == 0.0
        # (e) the build refuses when distinct orbit points are closer than 10 delta_dup
        @test ex_rejects(() -> EX.collect_orbits(0.5, [sol(1, z0, -1), sol(1, (BigFloat(5dup), zero(BigFloat)), 3)];
                                                 rho = rho, delta_dup = dup), ["delta_dup"])
    end
end
