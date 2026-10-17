# AC3 screening grid (E10), AC4 torus distance (E8).

@testset "AC3 grid" begin
    for N in (1, 2, 3, 5, 32, 256)
        g = EX.grid_centres(N)
        @test length(g) == N^2                                    # index id + 1 holds the centre of grid cell id = j N + i
        # unique cells, ordered by id = j N + i (x runs fastest)
        @test length(unique(g)) == N^2
        @test all(id -> g[id + 1][1] == g[(id % N) + 1][1] && g[id + 1][2] == g[(id ÷ N) * N + 1][2], 0:(N^2 - 1))
        @test all(c -> 0 < c[1] < EX_TWO_PI && 0 < c[2] < EX_TWO_PI, g)
        # 0 ulp against the binary64 nearest 2 pi (2 i + 1) / (2 N) evaluated at 1024 bits
        ref = setprecision(BigFloat, EX_BITS) do
            [Float64(2 * BigFloat(pi) * (2i + 1) / (2N)) for i in 0:(N - 1)]
        end
        @test all(id -> g[id + 1][1] == ref[(id % N) + 1] && g[id + 1][2] == ref[(id ÷ N) + 1], 0:(N^2 - 1))
        if N >= 2                                                  # N = 1 has a single column: the neighbour is the cell itself
            @test all(i -> abs(EX.torus_distance(g[i + 1], g[i + 2]) - EX_TWO_PI / N) <= 2 * eps(EX_TWO_PI), 0:(N - 2))
        end
    end
end

@testset "AC4 torus distance" begin
    d = EX.torus_distance
    u, v = (0.1, 0.1), (EX_TWO_PI - 0.1, EX_TWO_PI - 0.1)
    # the spec's 4 eps_53 relative is below the rounding of the Float64 input 2 pi - 0.1 itself (ulp 8.9e-16); 4 eps_53 of 2 pi is used
    @test abs(d(u, v) - 0.2 * sqrt(2)) <= 4 * EX_EPS * EX_TWO_PI
    rng = StableRNG(ex_seed(4))
    pts = [(EX_TWO_PI * rand(rng), EX_TWO_PI * rand(rng)) for _ in 1:2000]
    @test all(i -> d(pts[i], pts[i + 1]) == d(pts[i + 1], pts[i]), 1:1999)         # exactly symmetric
    @test all(p -> d(p, p) == 0, pts)
    @test all(i -> pts[i] == pts[i + 1] || d(pts[i], pts[i + 1]) > 0, 1:1999)
    @test all(i -> d(pts[i], pts[i + 1]) <= pi * sqrt(2), 1:1999)
    # triangle inequality on 10^4 seeded triples, slack 8 eps ||.|| with ||.|| <= pi sqrt 2
    slack = 8 * EX_EPS * pi * sqrt(2)
    trip = [(EX_TWO_PI * rand(rng, 6)) for _ in 1:10_000]
    @test all(t -> d((t[1], t[2]), (t[5], t[6])) <= d((t[1], t[2]), (t[3], t[4])) + d((t[3], t[4]), (t[5], t[6])) + slack, trip)
    # translation by 2 pi in either coordinate is the same point of the torus (to the rounding of 2 pi)
    @test all(p -> d(p, (p[1] + EX_TWO_PI, p[2])) <= 4 * EX_EPS * 8, pts[1:200])
    # the distance to the line p = c0 at K = 0 is delta(p - c0): the periodic set of K = 0 is the union of the lines of E7
    ps = EX.periodic_set(0.0; q_max = 3, N_g = 8, rho = EX_RHO, newton_max = 50, delta_dup = EX_DUP)
    for _ in 1:200
        z = (EX_TWO_PI * rand(rng), EX_TWO_PI * rand(rng))
        want = minimum(EX.torus_delta(z[2] - EX_TWO_PI * a / b) for (a, b) in ps.lines)
        @test abs(EX.periodic_distance(z, ps) - want) <= 32 * EX_EPS
    end
    @test EX.torus_delta(EX_TWO_PI + 0.25) ≈ 0.25 atol = 4 * EX_EPS * 8
    @test EX.torus_delta(-0.25) == 0.25
    @test EX.torus_delta(Float64(pi)) <= pi
end
