# AC6 to AC9: pathwise agreement with L1 (Benettin) on the identical orbit and direction, d0 sweep (U-shape), Double64, tau sweep.
# Wrapped mode, M_state = 2 pi, eps = effective eps of the element type. Orbit ids are the positions in the start vector.

function l2_ac6(K, starts, ns, seed)
    sys = WrappedStandardMap(K)
    x0s = [torus_state(z) for z in starts]
    Gf = [l2_Gfirst(sys, x, 1) for x in x0s]
    eps64 = Float64(effective_eps(Float64))
    outside = 0
    worst, worst_old = 0.0, 0.0
    for n in ns
        X = l2_X(sys, x0s, n, seed)
        Psi = l2_Psis(sys, x0s, seed, Float64, n)
        for m in 3:12
            d0 = 10.0^(-m)
            if !l2_valid(K, d0, 1, L2_TWO_PI, eps64)
                outside += 1
                continue
            end
            X2 = l2_X2(sys, x0s, n, seed, d0)
            D = l2_D2T(K, d0, L2_TWO_PI, eps64)
            bound = [l2_pathbound(D, Psi[j], n, l2_BS(Float64, n, Gf[j], K)) for j in eachindex(x0s)]
            @test all(j -> abs(X2[j] - X[j]) <= bound[j], eachindex(x0s))
            worst = max(worst, maximum(abs.(X2 .- X) ./ bound))
            worst_old = max(worst_old, maximum(abs.(X2 .- X)) / D)
        end
    end
    return outside, worst, worst_old
end

@testset "AC6 pathwise agreement with L1" begin
    # K in {0.5 (regular class, S7 screening), 2, 5, 2 pi}, 50 starts, Float64, wrapped, tau = 1, d0 = 10^-m (m = 3..12),
    # n in {100, 1e3} and 1e4 for K = 0.5: |X2 - X| <= Delta_2T(d0) Psi_n/n + B_S^L1(n)/n (T5a, Amendment D35) whenever the validity limits of T5 hold
    seed = l2_seed(6)
    for (iK, K) in enumerate((0.5, 2.0, 5.0, L2_TWO_PI))
        starts = K == 0.5 ? l2_regular_starts(K, 50, StableRNG(seed + iK)) : l2_starts(StableRNG(seed + iK), 50)
        outside, worst, worst_old = l2_ac6(K, starts, K == 0.5 ? (100, 1000, 10_000) : (100, 1000), seed + 10 + iK)
        @info "lyap-two-trajectory AC6: K = $K, grid points outside the validity limits (reported) = $outside, largest realized |X2 - X| / bound (T5a) = $worst, largest |X2 - X| / Delta_2T = $worst_old"
    end
end

# median over starts of |X2 - X| for d0 = 10^-m (Float64) or exact(Double64, "1e-m"), wrapped K, n = 1e3, tau = 1; cached
const L2_GRID_CACHE = Dict{Any, Any}()
function l2_grid(::Type{T}, K, nstart, ms, seed) where {T}
    return get!(L2_GRID_CACHE, (T, K, nstart, ms, seed)) do
        sys = WrappedStandardMap(T(K))
        x0s = [l2_state(:wrapped, T, z) for z in l2_starts(StableRNG(seed), nstart)]
        X = l2_X(sys, x0s, 1000, seed + 1)
        Dict(m => l2_median([Float64(abs(a - b)) for (a, b) in zip(l2_X2(sys, x0s, 1000, seed + 1, T === Float64 ? 10.0^(-m) : exact(T, "1e-$m")), X)]) for m in ms)
    end
end

@testset "AC7 d0 sweep U-shape" begin
    # Float64, K = 5, wrapped, 200 starts, n = 1e3, m = 2..13, e(m) = median over starts of |X2 - X|: (a) slope of log10 e against log10 d0 over m in {2,3,4}
    # lies in [0.8, 2.2]; (b) over m in {10,...,13} in [-1.3, -0.7]; (c) the grid minimum lies at d0 in [d0*/100, 100 d0*] and is at most Delta_2T(d0*)
    seed = l2_seed(7)
    eps64 = Float64(effective_eps(Float64))
    e = l2_grid(Float64, 5.0, 200, 2:13, seed)
    slope(ms) = l2_slope([-Float64(m) for m in ms], [log10(e[m]) for m in ms])
    @test 0.8 <= slope(2:4) <= 2.2
    @test -1.3 <= slope(10:13) <= -0.7
    mstar = argmin(Dict(m => e[m] for m in 2:13))
    d0star = l2_d0star(5.0, L2_TWO_PI, eps64)
    @test d0star / 100 <= 10.0^(-mstar) <= 100 * d0star
    @test e[mstar] <= l2_D2T(5.0, d0star, L2_TWO_PI, eps64)
end

@testset "AC8 Double64" begin
    # Double64, K in {5, 2 pi}, 50 starts, n = 1e3, d0 = 10^-m (m = 3..28): criterion 6 (T5a) with eps = 2^-100; at K = 5 the grid minimum of the
    # median error is <= 1e-9 and >= 1e3 times below the Float64 grid minimum of AC7
    seed = l2_seed(8)
    eps = Float64(effective_eps(Double64))
    meds = Dict{Int, Float64}()       # K = 5: median error per m, for the grid minimum
    for (iK, K) in enumerate((5.0, L2_TWO_PI))
        sys = WrappedStandardMap(Double64(K))
        starts = l2_starts(StableRNG(seed + iK), 50)
        x0s = [l2_state(:wrapped, Double64, z) for z in starts]
        n = 1000
        X = l2_X(sys, x0s, n, seed + 10 + iK)
        Gf = [l2_Gfirst(sys, x, 1) for x in x0s]
        Psi = l2_Psis(sys, x0s, seed + 10 + iK, Double64, n)
        for m in 3:28
            d0 = exact(Double64, "1e-$m")
            l2_valid(K, Float64(d0), 1, L2_TWO_PI, eps) || continue
            X2 = l2_X2(sys, x0s, n, seed + 10 + iK, d0)
            D = l2_D2T(K, Float64(d0), L2_TWO_PI, eps)
            @test all(j -> Float64(abs(X2[j] - X[j])) <= l2_pathbound(D, Psi[j], n, l2_BS(Double64, n, Gf[j], K)), eachindex(x0s))
            K == 5.0 && (meds[m] = l2_median([Float64(abs(X2[j] - X[j])) for j in eachindex(x0s)]))
        end
    end
    ef = l2_grid(Float64, 5.0, 200, 2:13, l2_seed(7))
    emin, efmin = minimum(values(meds)), minimum(values(ef))
    @info "lyap-two-trajectory AC8: grid minimum of the median error, Double64 $emin, Float64 $efmin, ratio $(efmin / emin)"
    @test emin <= 1e-9
    @test efmin / emin >= 1e3
end

@testset "AC9 tau sweep" begin
    # Float64, d0 = 1e-8, 50 starts, n = 1e3: |X2 - X| <= sigma^(2(tau-1)) Delta_2T(d0) Psi_n/n + B_S^L1(n)/n (T6a, Amendment D35); K = 0.5 (regular) with tau in {1,2,4}, K = 5 with tau in {1,2};
    # K = 5 with tau in {4, 8} is reported only; negative control T7: K = 10, d0 = 1e-3, tau = 16, N = 1008 (a multiple of 16): X2 <= ln(pi sqrt 2/d0)/16
    seed = l2_seed(9)
    eps64 = Float64(effective_eps(Float64))
    n, d0 = 1000, 1e-8
    for (K, taus, reported) in ((0.5, (1, 2, 4), ()), (5.0, (1, 2), (4, 8)))
        starts = K == 0.5 ? l2_regular_starts(K, 50, StableRNG(seed)) : l2_starts(StableRNG(seed + 1), 50)
        sys = WrappedStandardMap(K)
        x0s = [torus_state(z) for z in starts]
        X = l2_X(sys, x0s, n, seed)
        Gf = [l2_Gfirst(sys, x, 1) for x in x0s]
        Psi = l2_Psis(sys, x0s, seed, Float64, n)
        for tau in (taus..., reported...)
            X2 = l2_X2(sys, x0s, n, seed, d0; tau = tau)
            D = l2_D2T(K, d0, L2_TWO_PI, eps64)
            bound = [l2_pathbound(D, Psi[j], n, l2_BS(Float64, n, Gf[j], K); factor = l2_sigma(K)^(2 * (tau - 1))) for j in eachindex(x0s)]
            ratio = maximum(abs.(X2 .- X) ./ bound)
            if tau in taus
                @test ratio <= 1
                @info "lyap-two-trajectory AC9 (gated): K = $K, tau = $tau, largest |X2 - X| / bound = $ratio"
            else
                @info "lyap-two-trajectory AC9 (reported, not gated): K = $K, tau = $tau, largest |X2 - X| / bound = $ratio"
            end
        end
    end
    sys = WrappedStandardMap(10.0)
    x0s = [torus_state(z) for z in l2_starts(StableRNG(seed + 2), 50)]
    N, d0c = 1008, 1e-3
    X2 = l2_X2(sys, x0s, N, seed, d0c; tau = 16)
    X = l2_X(sys, x0s, N, seed)
    cap = log(pi * sqrt(2) / d0c) / 16
    @test all(X2 .<= cap + 1e-12)
    @test all(j -> X[j] < 1.1 || X[j] - X2[j] >= 0.5, eachindex(X))
    @info "lyap-two-trajectory AC9 control: starts with X >= 1.1 = $(count(>=(1.1), X)) of 50, cap = $cap"
end
