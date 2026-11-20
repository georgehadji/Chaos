# Amendment D126, stage 02: the finite-time forms of the E58 pairing threshold (d > 2) and of the exact-exponent control tolerance.
# Spec: experiments-spec.md, "Amendment D126". Both testsets run on lower-layer functions only (propagate, QRSpectrum), so they pass before stage 03.
# Constants: the D126 defaults, pinned from the K24 calibration (plan, "D126"); not changed after the data. Seeds: d125_seed('x', .) and d125_seed('y', .) (64200 to 64299).
# D126_PAIR_CONST and D126_CTL_CONST (pair_const, control_const) are defined in ex_d125_helpers.jl (ex_d125_items.jl uses the first)

# E58 as restated by D126, rounding term replaced by its lower bound 64 d eps_53 (the E58 term 16 d eps (G_max / tau + M_max) is not below it for G_max >= 0, M_max >= 1)
d126_pair_bound(d, n, tau, h) = D126_PAIR_CONST / (n * tau * h) + 64 * d * EX_EPS
d126_ctl_bound(d, n, tau, h, ex) = D126_CTL_CONST / (n * tau * h) + 64 * d * EX_EPS * max(1.0, abs(ex))

@testset "D126a pairing |lambda_i + lambda_(d+1-i)| <= C_pair / (N tau h) + rounding, d > 2 (E58)" begin
    rng = StableRNG(d125_seed('x'))
    kinds = Any[("y3", 4, () -> d125_system("froeschle", ["0.5", "0.7", "0.3"], Float64), 1.0),
                ("ch3r", 6, () -> d125_system("chain", ["3", "0.8", "0.35", "ring"], Float64), 1.0),
                ("ch5r", 10, () -> d125_system("chain", ["5", "1.1", "0.2", "ring"], Float64), 1.0),
                ("hh-lf", 4, () -> d125_flow_map(d125_integrators()[1], Float64), 0.0625)]
    worst = 0.0
    for (nm, d, sf, h) in kinds, N in (200, 1000), j in 1:25
        sys = sf()
        x0 = if nm == "hh-lf"
            z = nothing
            while z === nothing
                q2 = -0.3 + 0.7 * rand(rng); p2 = -0.25 + 0.5 * rand(rng)
                z, st = Chaoslab.energy_start(Chaoslab.HenonHeilesFlow(), 1 / 12, q2, p2); st === :ok || (z = nothing)
            end
            Float64.(collect(z))
        else
            [2pi * rand(rng) for _ in 1:d]
        end
        lam = d125_lambda(sys, x0, N, "L5", d, EX_SEED + 65000 + j, j) ./ h                  # per unit time for the flow, per step for the maps
        pair = maximum(abs(lam[i] + lam[d + 1 - i]) for i in 1:(d ÷ 2))
        worst = max(worst, pair / d126_pair_bound(d, N, 1, h))
        @test pair <= d126_pair_bound(d, N, 1, h)
    end
    @test worst < 1.0
    # the rounding-level threshold of D125 is exceeded by orders of magnitude (the reason for D126): the same sample against 64 d eps alone
    sys = d125_system("froeschle", ["0.5", "0.7", "0.3"], Float64)
    lam = d125_lambda(sys, [1.0, 2.0, 3.0, 4.0], 1000, "L5", 4, EX_SEED + 65001, 1)
    @test maximum(abs(lam[i] + lam[5 - i]) for i in 1:2) > 1e6 * 64 * 4 * EX_EPS
end

@testset "D126b controls |lambda_i - exact_i| <= C_ctl / (N tau h) + 64 d eps max(1, |exact_i|) (E54, R42 item 4)" begin
    function ymap(mus)
        d = 2 * length(mus); A = zeros(d, d)
        setprecision(BigFloat, 256) do
            for (i, m) in enumerate(mus)
                A[2i-1:2i, 2i-1:2i] .= Float64.(EX.hyperbolic_block(BigFloat, log(d125_big(m))))
            end
        end
        return Chaoslab.LinearSymplecticMap(Chaoslab.SMatrix{d, d, Float64}(A))
    end
    cases = [("y5-d4", () -> ymap(["3", "5/2"]), d125_exact_y5(["3", "5/2"])), ("y5-d6", () -> ymap(["2", "3/2", "4"]), d125_exact_y5(["2", "3/2", "4"])),
             ("chain2", () -> d125_system("chain", ["2", "1.5", "0", "open"], Float64), d125_exact_chain_origin(2, "1.5")),
             ("chainK0", () -> d125_system("chain", ["2", "0", "0", "open"], Float64), zeros(4))]
    for N in (100, 160, 1000), (nm, sf, ex) in cases
        sys = sf(); d = length(ex)
        for id in 1:12, (m, k) in (("L5", d), ("L1", 1), ("L3", d))
            lam = d125_lambda(sys, zeros(d), N, m, k, EX_RUN_SEED, id)
            @test all(i -> abs(lam[i] - ex[i]) <= d126_ctl_bound(d, N, 1, 1, ex[i]), 1:k)
        end
    end
    # the rounding-level tolerance of D125 is not attainable by a finite-time exponent: N = 100, y5 d = 4, frame 1
    sys = cases[1][2](); lam = d125_lambda(sys, zeros(4), 100, "L5", 4, EX_RUN_SEED, 1)
    @test maximum(abs.(lam .- cases[1][3])) > 1e8 * 64 * 4 * EX_EPS * 3
end

# D126e (owner: "Ναι πρόσθεσε", stage 02 checkpoint): negative controls. The finite-time bounds of D126a and D126b must still reject a real defect.
# Detection limit by construction: a per-step defect delta is seen once delta > C / (N tau h), so the assertions sit at N = 1e5 (bound 3.7e-4 and 5.6e-4)
# and the same defects are below the bound at N = 1000 (0.037, 0.056), which is recorded, not hidden. No random draw: linear maps, start 0, frames from EX_RUN_SEED.
function d126_blockmap(mus, scales)
    d = 2 * length(mus); A = zeros(d, d)
    setprecision(BigFloat, 256) do
        for (i, m) in enumerate(mus)
            A[2i-1:2i, 2i-1:2i] .= scales[i] .* Float64.(EX.hyperbolic_block(BigFloat, log(d125_big(m))))
        end
    end
    return Chaoslab.LinearSymplecticMap(Chaoslab.SMatrix{d, d, Float64}(A))   # the constructor does not check symplecticity: the defect is in the tangent map
end

@testset "D126e negative controls: a non-symplectic tangent map and a wrong exact exponent are rejected at N = 1e5 (E58, E54)" begin
    s = 1.01; delta = 2 * log(s)                                                    # pairing defect per step, 0.0199
    # blocks scaled by s and 1/s: determinant 1, so only the pairing sees the defect (the volume and sum checks do not)
    bad = d126_blockmap(["3", "5/2"], (s, 1 / s))
    good = d126_blockmap(["3", "5/2"], (1.0, 1.0))
    for N in (1000, 100_000), id in 1:4
        pb = d125_lambda(bad, zeros(4), N, "L5", 4, EX_RUN_SEED, id)
        pg = d125_lambda(good, zeros(4), N, "L5", 4, EX_RUN_SEED, id)
        pairb = maximum(abs(pb[i] + pb[5 - i]) for i in 1:2)
        pairg = maximum(abs(pg[i] + pg[5 - i]) for i in 1:2)
        @test abs(abs(sum(pb)) - 0.0) <= 1e-12 * N                                 # the sum is blind to this defect (det 1)
        @test abs(pairb - delta) <= D126_PAIR_CONST / N                             # the measured defect is the planted one
        @test pairg <= d126_pair_bound(4, N, 1, 1)                                  # the correct map passes at both N
        if N == 100_000
            @test pairb > d126_pair_bound(4, N, 1, 1)                               # rejected at the study scale
        else
            @test pairb <= d126_pair_bound(4, N, 1, 1)                              # below the detection limit C_pair / N at N = 1000 (recorded)
        end
    end
    # control with a wrong exact exponent: map with mu = 3, 5/2; the control claims mu = 301/100 (offset ln(301/300) = 3.3e-3 on two exponents)
    sys = d126_blockmap(["3", "5/2"], (1.0, 1.0))
    exr = d125_exact_y5(["3", "5/2"]); exw = d125_exact_y5(["301/100", "5/2"])
    for N in (1000, 100_000), id in 1:4
        lam = d125_lambda(sys, zeros(4), N, "L5", 4, EX_RUN_SEED, id)
        @test all(i -> abs(lam[i] - exr[i]) <= d126_ctl_bound(4, N, 1, 1, exr[i]), 1:4)
        caught = any(i -> abs(lam[i] - exw[i]) > d126_ctl_bound(4, N, 1, 1, exw[i]), 1:4)
        @test caught == (N == 100_000)                                              # rejected at 1e5, below the limit C_ctl / N at 1000
    end
end
