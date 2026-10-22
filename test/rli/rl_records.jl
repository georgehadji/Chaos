# Recorded lines of the spec that need no threshold and no long run: R1 (the 2D map of SandorMaffione2016 Eq. (2), with the exact identity F7) and R5 (sensitivity to the vector
# and the direction). R2 is the @info of AC1, AC3, AC6, AC7c, AC8 and AC10(b); R3 the @info of AC5; R4, R6 and R7 are stage 04 (long runs, wall time); R8 is the @info of AC7d.
# Every testset starts with a line naming the new symbols it uses. These testsets assert only the identity, the run length and finiteness; the values are reported.

@testset "R1 the map of SandorMaffione2016 Eq. (2) and the identity F7" begin
    _ = (RLI, PairMap)
    # F7: q' = q + y, y' = y - nu sin q' is the repository standard map p' = p + K sin x, x' = x + p' with K = -nu, x = q, p_0 = y_0 + nu sin q_0 and p_{n+1} = y_n. For y_0 = 0 the
    # identity holds bitwise over 50 iterates (the same floating operations in the same order). Then RLI(n) for the ordered start (q, y) = (2, 0) and the chaotic start (3, 0) at nu = 0.5,
    # n in {1e2, 1e3, 1e4, 2e4}, Dx in {1e-14, 1e-12, 1e-10} (wrapped, tau = 1, xi_0 seeded): reported, with the ratio chaotic / ordered in orders of magnitude against the printed 9 to 10
    # (SandorMaffione2016 p.5; separation and vector are not stated there).
    nu = 0.5
    sys = StandardMap(-nu)
    for (q0, y0) in ((2.0, 0.0), (3.0, 0.0))
        z = SVector(q0, y0 + nu * sin(q0))
        q, y = q0, y0
        nbad = 0
        for _ in 1:50
            z = forward(sys, z)
            q = q + y
            nbad += !(z[1] == q && z[2] == y)
            y = y - nu * sin(q)
        end
        @test nbad == 0
    end
    wsys = WrappedStandardMap(-nu)
    ns = (100, 1000, 10_000, 20_000)
    vals = Dict{Tuple{Float64, Float64, Int}, Float64}()
    for (id, (q0, y0)) in enumerate(((2.0, 0.0), (3.0, 0.0))), Dx in (1e-14, 1e-12, 1e-10)
        xi = rl_xi0(rl_seed(16, 1), id, Float64)
        pr = rl_pair(wsys, SVector(q0, y0 + nu * sin(q0)), xi, Dx)
        r = rl_run((RLI(),), PairMap(wsys), [pr.z], [rl_frame(xi)], 20_000, 1; steps = collect(ns))
        snaps = r[1][1].snaps
        @test length(snaps) == length(ns) && all(f -> isfinite(f.RLI) && f.RLI >= 0, snaps)
        for (f, n) in zip(snaps, ns)
            vals[(q0, Dx, n)] = f.RLI
        end
    end
    for Dx in (1e-14, 1e-12, 1e-10), n in ns
        @info "R1 RLI(n) at nu = 0.5: ordered (2, 0), chaotic (3, 0), orders of magnitude between them" Dx n ordered = vals[(2.0, Dx, n)] chaotic = vals[(3.0, Dx, n)] orders = log10(vals[(3.0, Dx, n)] / vals[(2.0, Dx, n)])
    end
end

@testset "R5 sensitivity to the vector and the direction" begin
    _ = (RLI, PairMap, check_frame)
    # P2, P3. Regular-class starts (K = 0.5), each with supplied unit vectors xi_0 of uniformly random orientation (checked by check_frame, as the caller must) and shadow directions u of
    # random orientation: the distribution of RLI / Dx_r at n = 1e3, wrapped, Dx = 1e-9, tau = 1 (reported; D38: 5 starts x 5 vectors x 5 directions; the spec: 20 x 20 and n = 1e3, 1e4).
    # Asserted: every frame is accepted, every run is finite.
    rng = StableRNG(rl_seed(16, 2))
    c = rl_case(:half, Float64, 5)
    rhos = Float64[]
    for j in 1:5, _ in 1:5, _ in 1:5
        a, b = 2pi * rand(rng), 2pi * rand(rng)
        xi, u = SVector(cos(a), sin(a)), SVector(cos(b), sin(b))
        W = rl_frame(xi)
        @test rl_health_error(() -> check_frame(RLI(), W; orbit = j)) === nothing
        pr = rl_pair(c.sys, c.x0[j], u, 1e-9)
        f = rl_run((RLI(),), PairMap(c.sys), [pr.z], [W], 1000, 1; steps = [1000])[1][1].snaps[1]
        @test isfinite(f.RLI) && f.RLI >= 0
        push!(rhos, f.RLI / pr.Dxr)
    end
    @info "R5 distribution of RLI / Dx_r at n = 1e3 (regular class, random vectors and directions)" extrema(rhos) median = sort(rhos)[(length(rhos) + 1) ÷ 2]
end
