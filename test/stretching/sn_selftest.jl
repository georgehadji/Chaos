# Self-test of the reference code of sn_helpers.jl against hand-computed values. It uses no name of the module under test, so it must pass before and after stage 03 (GREEN in the RED run).

@testset "SELF reference code against hand values" begin
    # cell rule SN4: range [-1, 1], B = 4 (h = 1/2): -1.5 -> 0, -1 -> 1 (z = s_lo), -0.75 -> 1, 0 -> 3, 0.49 -> 3, 0.5 -> 4, 1 -> 4 (last cell closed at s_hi), just above 1 -> 5
    zs = [-1.5, -1.0, -0.75, 0.0, 0.49, 0.5, 1.0, 1.0000001]
    @test sn_cells_ref(zs, -1.0, 1.0, 4) == [0, 1, 1, 3, 3, 4, 4, 5]
    @test sn_hist_ref(zs, -1.0, 1.0, 4) == [1, 2, 0, 2, 2, 1]
    @test sum(sn_hist_ref(zs, -1.0, 1.0, 4)) == length(zs)
    # D of SN5: two single-cell vectors in different cells give sqrt(2 / h); equal vectors give 0; a hand case with k = 4
    @test Float64(sn_D_h([0, 3, 0], [0, 0, 3], 0.5)) == 2.0
    @test Float64(sn_D_h([0, 3, 0], [0, 3, 0], 0.5)) == 0.0
    @test sn_S([0, 1, 3, 0], [1, 1, 1, 1]) == 1 + 0 + 4 + 1
    @test isapprox(Float64(sn_D_h([0, 1, 3, 0], [1, 1, 1, 1], 0.25)), sqrt(6.0) / (4 * sqrt(0.25)); rtol = 1e-15)
    @test isapprox(Float64(sn_D_ref([0, 1, 3, 0], [1, 1, 1, 1], 0.0, 1.0, 2)), sqrt(6.0) / (4 * sqrt(0.5)); rtol = 1e-15)   # B = 2 on [0, 1]: h = 1/2
    # prefix counts and the envelope on one hand sample: z = 0.1 lies on the edge between cells 1 and 2 of [0, 1] with B = 10; the interval 0.1 +- 0.01 meets both
    pc = sn_prefix_counts([0.05, 0.1, 0.95], 0.0, 1.0, 10, [2, 3])
    @test pc[2] == [0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0] && pc[3] == [0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0]
    env = sn_envelope([0.35], [0.001], 0.0, 1.0, 10, [1])[1]
    @test env.U == 0 && env.N == env.Nm == env.Np && env.N[5] == 1
    env = sn_envelope([0.1], [0.01], 0.0, 1.0, 10, [1])[1]
    @test env.U == 1 && env.Nm == zeros(Int, 12) && env.Np[2] == 1 && env.Np[3] == 1 && sum(env.Np) == 2 && sum(env.N) == 1
    # edge condition and the choice of B: z = 0.5 is an edge for B = 4 (h = 1/2) and not for B = 5 (3.75 cells)
    @test !sn_edge_ok([0.5], -1.0, 1.0, 4) && sn_edge_ok([0.5], -1.0, 1.0, 5)
    @test sn_pick_B(4, -1.0, 1.0, [0.5]) == (5, true) && sn_pick_B(5, -1.0, 1.0, [0.5]) == (5, false)
    @test !sn_edge_ok([0.99999999999], 0.0, 1.0, 10; endslack = 1e-6) && sn_edge_ok([0.99999999999], 0.0, 1.0, 10; endslack = 0.0)
    # z from the step-indexed sums: s = (1, 3, 6, 10), tau = 2 gives (3 - 0) / 2 and (10 - 3) / 2
    @test sn_z_from_s(BigFloat.([1, 3, 6, 10]), 2) == [BigFloat(3) / 2, BigFloat(7) / 2]
    # moments: the series (1, 2, 3); the histogram moments of the counts (0, 2, 2, 0) on [0, 2] with B = 2
    mom = sn_prefix_moments(BigFloat.([1, 2, 3]), [3])[3]
    @test mom[1] == 2 && abs(mom[2] - sqrt(BigFloat(2) / 3)) < 1e-70
    m1, sd = sn_hist_moments([0, 2, 2, 0], 0.0, 1.0)
    @test m1 == 1.0 && sd == 0.5
    # the faulty doubles: drop, clip, abs, window
    z4 = [-2.0, -0.5, 0.5, 2.0]
    @test sn_hist_fault(z4, -1.0, 1.0, 2, :drop) == [0, 1, 1, 0] && sn_hist_fault(z4, -1.0, 1.0, 2, :clip) == [0, 2, 2, 0] && sn_hist_fault(z4, -1.0, 1.0, 2, :abs) == [0, 0, 2, 2]
    @test sn_hist_fault(z4, -1.0, 1.0, 2, :window; window = 2) == [0, 0, 1, 1]
    # grids and intervals
    @test sn_grid(25) == vcat(1:9, [10, 20]) && sn_grid(1000)[end] == 1000 && issorted(sn_grid(10_000))
    @test sn_interval(1) == 100 && sn_interval(3) == 102 && sn_interval(4) == 100 && sn_interval(16) == 112
    @test sn_eval_steps(450, 100) == [100, 200, 300, 400, 450] && sn_eval_steps(400, 100) == [100, 200, 300, 400]
    # L_sigma of L1 E14: K = 0 is ln of the golden ratio; the printed values of the spec (F6: 0.725, 2.27, 2.07, 2.70) to the digits printed
    @test abs(sn_Lsigma(0.0) - log((1 + sqrt(BigFloat(5))) / 2)) < 1e-70
    for (K, v) in ((0.5, 0.725), (2pi, 2.27), (5.0, 2.07), (10.0, 2.70))
        @test abs(Float64(sn_Lsigma(K)) - v) < 0.0051
    end
    # T1 closed form (SN8d) against the direct window means of c x_j, x_j = x0 + j p0, in BigFloat (tau = 3, K = 4)
    c, x0, p0, tau = 0.3, 0.1, 0.02, 3
    direct = setprecision(BigFloat, SN_BITS) do
        [sum(BigFloat(c) * (BigFloat(x0) + j * BigFloat(p0)) for j in ((k - 1) * tau):(k * tau - 1)) / tau for k in 1:4]
    end
    cf = sn_T1_z1(c, x0, p0, tau, 4)
    @test maximum(abs.(cf .- direct)) < 1e-70
    @test sn_T1_G(0.0, x0, p0, 3) == 2.0 && sn_T1_G(c, 0.0, 0.0, 3) == 2.0
    # Y5 closed form (SN8b) against the direct power of A = [[2, 1], [1, 1]] in BigFloat (w = (0.6, 0.8), tau = 2, K = 4)
    w = (0.6, 0.8)
    direct = setprecision(BigFloat, SN_BITS) do
        v = (BigFloat(w[1]), BigFloat(w[2]))
        nrm(v) = sqrt(v[1]^2 + v[2]^2)
        prev = nrm(v)
        out = BigFloat[]
        for k in 1:4
            for _ in 1:2
                v = (2 * v[1] + v[2], v[1] + v[2])
            end
            push!(out, log(nrm(v) / prev) / 2)
            prev = nrm(v)
        end
        out
    end
    @test maximum(abs.(sn_y5_z(w, 2, 4) .- direct)) < 1e-70
    @test abs(sn_mu() - (3 + sqrt(BigFloat(5))) / 2) < 1e-70
    @test all(z -> z < log(sn_mu()), sn_y5_z(w, 2, 4)) && abs(sn_y5_z(w, 1, 40)[40] - log(sn_mu())) < 1e-30
    # shear closed form (SN8c) against the direct window norms (a, b) = (0.6, 0.8), tau = 3, K = 5; b = 0 gives exactly 0; tau k z_k -> 1
    direct = setprecision(BigFloat, 128) do
        a, b = BigFloat(0.6), BigFloat(0.8)
        nr(m) = sqrt((a + m * b)^2 + b^2)
        [log(nr(3k) / nr(3 * (k - 1))) / 3 for k in 1:5]
    end
    @test maximum(abs.(sn_shear_z(0.6, 0.8, 3, 5) .- direct)) < 1e-35
    @test all(iszero, sn_shear_z(1.0, 0.0, 3, 5))
    zk = sn_shear_z(0.6, 0.8, 1, 1000)
    @test abs(1000 * zk[1000] - 1) < 0.02
    # the multinomial sampler: mean counts within 5 standard deviations of Z p (seeded)
    rng = StableRNG(sn_seed(0, 40))
    p = sn_dirichlet(rng, 10)
    @test isapprox(sum(p), 1; atol = 1e-14)
    tot = zeros(Int, 10)
    for _ in 1:200
        tot .+= sn_multinomial(rng, 500, cumsum(p))
    end
    @test all(abs.(tot ./ 200 .- 500 .* p) .<= 5 .* sqrt.(500 .* p .* (1 .- p) ./ 200) .+ 1e-9) && sum(tot) == 200 * 500
    # a reference orbit on a tiny case: the K = 0 shear, column (1, 0) stays and column (0, 1) gives (m, 1): s2(m) = ln sqrt(m^2 + 1)
    sys = StandardMap(0.0)
    o = sn_ref_orbit(sys, 0.0, SVector(0.3, 0.7), SMatrix{2, 2, Float64}(1, 0, 0, 1), 6)
    @test all(iszero, o.s1) && abs(o.s2[3] - log(sqrt(BigFloat(10)))) < 1e-70 && length(o.Ms) == 6
    @test issorted(sn_Pmax(o, 2, 6)) && all(>(0), sn_Pmax(o, 2, 6))
end
