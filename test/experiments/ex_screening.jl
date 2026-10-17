# AC11 windowed estimate, AC12 K = 0 screening (O2), AC25 screening against L1 (O5), AC26 zero allocation.
# Screening is L5 (Householder), Float64, k = 1, first deviation vector (cos alpha, sin alpha) (E11, P8).

const EX_N = 2000
const EX_W = 100
const EX_ALPHA = pi / 4

ex_scr(K, z; N = EX_N, w = EX_W, alpha = EX_ALPHA, tau = 1) = EX.screen_point(K, z[1], z[2]; N_scr = N, w = w, alpha = alpha, tau = tau)

# the allocation probes call the module through a static reference, never through EX (a property lookup would allocate)
ex_screen!(buf, K, x, p, N, w, a, tau) = Chaoslab.Experiments.screen_point!(buf, K, x, p, N, w, a, tau)
ex_alloc_screen(buf, K, x, p, N, w, a, tau) = @allocated ex_screen!(buf, K, x, p, N, w, a, tau)

@testset "AC11 windowed estimate" begin
    rng = StableRNG(ex_seed(11))
    g = EX.grid_centres(16)
    for K in (0.0, 1.5, 5.0), _ in 1:17                                   # 51 seeded grid cells over the fixture K values
        z = g[1 + floor(Int, 256 * rand(rng))]
        r = ex_scr(K, z)
        @test length(r.blocks) == EX_N ÷ EX_W
        @test abs(sum(r.blocks) / length(r.blocks) - r.lambda) <= 1e-12
        @test abs(EX_W * sum(r.blocks[1:(EX_N ÷ 10 ÷ EX_W)]) - r.S_tenth) <= 1e-9     # S(N_scr / 10) is the sum of the first blocks
        @test abs(EX_N * r.lambda - r.S_N) <= 1e-9
    end
end

@testset "AC12 K = 0 screening: closed form (O2)" begin
    # derived: at K = 0, DT = [[1,1],[0,1]], the deviation vector (cos a, sin a) becomes (cos a + n sin a, sin a), so
    # S(n) = ln sqrt((cos a + n sin a)^2 + sin^2 a) and lambda_scr = S(N)/N. Abs tolerance 8 eps: the S(n) = ... to 8 n eps of the
    # Householder suite (AC4 there); 8 eps relative as written in the spec is below the sum's rounding, see the plan, ambiguity A2.
    g = EX.grid_centres(8)
    rs = [ex_scr(0.0, z) for z in g]
    a, N = EX_ALPHA, EX_N
    closed = log(sqrt((cos(a) + N * sin(a))^2 + sin(a)^2)) / N
    @test allequal(r.lambda for r in rs)                                  # the tangent does not depend on the state at K = 0
    @test all(r -> abs(r.lambda - closed) <= 8 * EX_EPS, rs)
    # the growth ratio of regular decay: S(N)/S(N/10) from the same closed form
    S(n) = log(sqrt((cos(a) + n * sin(a))^2 + sin(a)^2))
    @test abs(rs[1].S_N / rs[1].S_tenth - S(N) / S(N ÷ 10)) <= 1e-10
end

@testset "AC25 screening against L1 (O5)" begin
    seed = ex_seed(25)
    g = EX.grid_centres(16)
    rng = StableRNG(seed)
    a = EX_ALPHA
    dir = SMatrix{2, 1, Float64}(cos(a), sin(a))
    total = 0
    for (K, m) in ((0.0, 67), (1.5, 67), (5.0, 66))                         # 200 seeded grid cells over the fixture K values
        cells = [g[1 + floor(Int, 256 * rand(rng))] for _ in 1:m]
        x0 = [SVector{2, Float64}(c[1], c[2]) for c in cells]
        l1 = propagate(StandardMap(K), x0, EX_N; k = 1, tau = 1, seed = seed, reducers = (Benettin(),), W0 = fill(dir, m))
        for j in 1:m
            r = ex_scr(K, cells[j])
            bs = ex_BS(EX_N, K, ex_Gfirst(K, cells[j][1]))
            @test abs(r.lambda - l1.finals[j][1].lambda[1]) <= 2 * bs / EX_N
            total += 1
        end
    end
    @test total == 200
    # unscored: lambda_scr against L5 at 10 N_scr for 50 chaotic cells at K = 10 (reported, not asserted)
    cells = [g[1 + floor(Int, 256 * rand(rng))] for _ in 1:50]
    x0 = [SVector{2, Float64}(c[1], c[2]) for c in cells]
    long = propagate(StandardMap(10.0), x0, 10 * EX_N; k = 1, tau = 1, seed = seed, reducers = (QRSpectrum(HouseholderQR()),), W0 = fill(dir, 50))
    d = [abs(ex_scr(10.0, cells[j]).lambda - long.finals[j][1].lambda[1]) for j in 1:50]
    @info "AC25 unscored: lambda_scr(N) against L5(10 N) at K = 10, 50 cells" max_abs_difference = maximum(d)
    @test all(isfinite, d)
end

@testset "AC26 zero allocation" begin
    buf = EX.screen_buffer(EX_N, EX_W)
    ex_screen!(buf, 1.5, 1.0, 2.0, EX_N, EX_W, EX_ALPHA, 1)                # warm-up
    @test ex_alloc_screen(buf, 1.5, 1.0, 2.0, EX_N, EX_W, EX_ALPHA, 1) == 0
    @test ex_alloc_screen(buf, 5.0, 3.0, 4.0, EX_N, EX_W, EX_ALPHA, 1) == 0
    # per grid cell: 2 state + 4 tangent numbers (48 bytes), a compensated accumulator (16), N_scr / w block exponents (8 each)
    for (N, w) in ((EX_N, EX_W), (1000, 100), (100_000, 1000))
        @test EX.screen_bytes(N, w) == 48 + 16 + 8 * (N ÷ w)
    end
    @test sizeof(buf.blocks) == 8 * (EX_N ÷ EX_W)
end
