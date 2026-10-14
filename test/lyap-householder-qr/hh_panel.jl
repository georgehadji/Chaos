# AC7 and AC13: standard map against a 256-bit run on the same orbit points, tau sweep.

# 256-bit reference exponents (windows of one step; S_i(n) does not depend on tau, I2) along the stored orbit points zs[1..n]:
# per step M = [[1 + c, 1], [c, 1]], c = K cos x (E8), one MGS step in BigFloat, product of the R_ii, one log at the end.
function hh_ref_std(K, zs, Q0, n)
    return setprecision(BigFloat, 256) do
        Kb = hh_big(K)
        a1, b1 = hh_big(Q0[1, 1]), hh_big(Q0[2, 1])
        a2, b2 = hh_big(Q0[1, 2]), hh_big(Q0[2, 2])
        P1 = P2 = one(Kb)
        for i in 1:n
            c = Kb * cos(hh_big(zs[i][1]))
            w1x, w1y = (1 + c) * a1 + b1, c * a1 + b1
            w2x, w2y = (1 + c) * a2 + b2, c * a2 + b2
            r11 = sqrt(w1x^2 + w1y^2)
            a1, b1 = w1x / r11, w1y / r11
            r12 = a1 * w2x + b1 * w2y
            vx, vy = w2x - r12 * a1, w2y - r12 * b1
            r22 = sqrt(vx^2 + vy^2)
            a2, b2 = vx / r22, vy / r22
            P1 *= r11
            P2 *= r22
        end
        [log(P1) / n, log(P2) / n]
    end
end

# stored orbit points z_0 .. z_{n-1} in T and the largest tangent entry over them
function hh_orbit(sys, x0, n)
    zs = [x0]
    m = 0.0
    for _ in 1:(n - 1)
        m = max(m, Float64(maximum(abs, tangent(sys, zs[end]))))
        push!(zs, forward(sys, zs[end]))
    end
    m = max(m, Float64(maximum(abs, tangent(sys, zs[end]))))
    return zs, m
end

const HH_AC7_KS = (0.5, 0.9716354, 2.0, 5.0, 10.0)   # K_test = the Float64 nearest the decimal, widened exactly to every type

function hh_ac7(::Type{T}, K::Float64, starts, seed) where {T}
    eps_T = Float64(effective_eps(T))
    m, n = length(starts), 1000
    sys = StandardMap(T(K))
    x0s = [hh_cast(T, z) for z in starts]
    W0s = [initial_deviations(seed, j, Val(2), Val(2), T) for j in 1:m]
    orb = [hh_orbit(sys, x0s[j], n) for j in 1:m]
    ref = [hh_ref_std(T(K), orb[j][1], W0s[j], n) for j in 1:m]
    Psi = [hh_Psi(sys, orb[j][1], W0s[j]) for j in 1:m]                                         # Psi_n of the orbit and the first column (L1 S3, Amendment D35)
    cprime = 0.0
    for tau in (1, 2, 4)
        lam, thetas = Dict(), Dict()
        for alg in hh_algs(T)
            res = propagate(sys, x0s, n; k = 2, tau = tau, seed = seed, W0 = W0s, reducers = (QRSpectrum(alg),))
            f = [res.finals[j][1] for j in 1:m]
            th = [hh_theta(T, f[j].G_max, tau, orb[j][2]) for j in 1:m]
            thp = [hh_theta_psi(T, f[j].G_max, tau, orb[j][2], Psi[j], n) for j in 1:m]         # theta_psi(n), Amendment D35
            errs = setprecision(BigFloat, 256) do
                [maximum(Float64(abs(hh_big(f[j].lambda[i]) - ref[j][i])) for i in 1:2) for j in 1:m]
            end
            @test all(errs .<= thp)                                                             # (a) vs the 256-bit run
            @test all(j -> Float64(abs(f[j].lambda[1] + f[j].lambda[2])) <= th[j], 1:m)         # (c)
            for j in 1:m
                cprime = max(cprime, Float64(abs(f[j].lambda[1] + f[j].lambda[2])) / (eps_T * (Float64(f[j].G_max) / tau + orb[j][2])))
            end
            lam[hh_name(alg)], thetas[hh_name(alg)] = [f[j].lambda for j in 1:m], thp
        end
        names = collect(keys(lam))
        for a in 1:length(names), b in (a + 1):length(names)                                    # (b) every pair of variants
            @test all(j -> maximum(Float64(abs(lam[names[a]][j][i] - lam[names[b]][j][i])) for i in 1:2) <= thetas[names[a]][j] + thetas[names[b]][j], 1:m)
        end
    end
    return cprime
end

@testset "AC7 standard map against a 256-bit run" begin
    # K in {0.5, 0.9716354, 2, 5, 10}, seeded starts (K = 0.5 screened: regular class), n = 1e3, k = 2, tau in {1, 2, 4}, every variant.
    # Float64 50 starts, Double64 10 starts (shrunk for the 3 minute budget; same starts as the first 10 of Float64).
    # (a) max_i |lambda_i - ref| <= theta_psi(n) = theta_sum + 8 eps Psi_n/n, (b) pairs of variants within theta_psi^A + theta_psi^B,
    # (c) |lambda_1 + lambda_2| <= theta_sum (all Amendment D35; no tail rule: every orbit is held to the bound)
    seed = hh_seed(7)
    cprime = 0.0
    for (iK, K) in enumerate(HH_AC7_KS)
        rng = StableRNG(seed + iK)
        starts = K == 0.5 ? hh_regular_starts(K, 50, rng)[1] : hh_starts(rng, 50)
        cprime = max(cprime, hh_ac7(Float64, K, starts, seed + 10 + iK))
        cprime = max(cprime, hh_ac7(Double64, K, starts[1:10], seed + 20 + iK))
    end
    @info "lyap-householder-qr AC7: largest measured c' of E11 (|lambda_1 + lambda_2| / (eps (G_max/tau + M_max)))" cprime
end

@testset "AC13 tau sweep" begin
    # K = 2 pi, k = 2, one seeded orbit, n = 1e4, tau in {1, 2, 4}: |lambda_1 + lambda_2| <= theta_sum and
    # |lambda_i(tau) - lambda_i(1)| <= theta_sum(tau) + theta_sum(1) (I2); both T, every variant
    seed = hh_seed(13)
    n = 10_000
    for T in (Float64, Double64)
        sys = StandardMap(T(HH_TWO_PI))
        x0 = hh_cast(T, hh_starts(StableRNG(seed), 1)[1])
        Mmax = hh_Mmax(sys, x0, n)
        for alg in hh_algs(T)
            f = Dict(tau => propagate(sys, [x0], n; k = 2, tau = tau, seed = seed, reducers = (QRSpectrum(alg),)).finals[1][1] for tau in (1, 2, 4))
            th = Dict(tau => hh_theta(T, f[tau].G_max, tau, Mmax) for tau in (1, 2, 4))
            for tau in (1, 2, 4)
                @test Float64(abs(f[tau].lambda[1] + f[tau].lambda[2])) <= th[tau]
                for i in 1:2
                    @test Float64(abs(f[tau].lambda[i] - f[1].lambda[i])) <= th[tau] + th[1]
                end
            end
        end
    end
end
