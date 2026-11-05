# AC1 to AC4: evaluation accuracy and route agreement, closed forms by prescribed Gram matrix, algebra of the definition, k = 2 against SALI.
# Each testset starts by naming the new symbols it uses, so a missing module errors that testset (UndefVarError) and nothing else.
# The negative controls (test-side doubles, no new symbol) are separate testsets and pass now.

const GA_PAIRS = ((2, 2), (3, 2), (4, 2), (8, 2), (4, 3), (4, 4), (6, 3), (6, 6), (8, 5))   # (d, K): K columns in d dimensions

@testset "AC1 evaluation accuracy and route agreement" begin
    _ = (gali_values, gali_floor)
    # Random sets of (d, K) in GA_PAIRS, raw column norms spread over 1e-3 .. 1e3. (a) in BigFloat(256) the four routes of B2 (minors, Gram determinant,
    # singular values, QR) agree to relative 1e-60; (b) the module evaluation of all orders 2..K from one K-column set satisfies |G_T - exact| <= 2 B_meas(T, k)
    # for every k, exact = BigFloat(256) evaluation of the same rounded input. D38: 20 sets per pair for the routes (the minors route of (8, 5) is 56
    # determinants) and 100 sets per pair and type for (b) (spec: 1000).
    rng = StableRNG(ga_seed(1))
    worst_route = 0.0
    for (d, K) in GA_PAIRS, _ in 1:ga_size(20, 1000)
        W = ga_randset(rng, Float64, d, K)
        setprecision(BigFloat, GA_BITS) do
            U = ga_unitcols(BigFloat[ga_big(W[i, j]) for i in 1:d, j in 1:K])
            v = ga_route_qr(U)
            for r in (ga_route_minors(U), ga_route_gram(U), ga_route_svd(U))
                worst_route = max(worst_route, Float64(abs(r - v) / v))
            end
        end
    end
    @test worst_route <= 1e-60
    for T in GA_TYPES
        worst = 0.0
        for (d, K) in GA_PAIRS, _ in 1:ga_size(100, 1000)
            W = ga_randset(rng, T, d, K)
            g = gali_values(W)
            @test g isa SVector{K - 1, T}
            ex = ga_exact(W)
            for k in 2:K
                worst = max(worst, ga_err(g[k - 1], ex[k - 1]) / (2 * ga_meas(T, k)))
            end
        end
        @test worst <= 1                                       # in units of 2 B_meas(T, k)
        @info "AC1 reported: worst error / (2 B_meas)" T worst
    end
    @test gali_floor(Float64, 2) == ga_meas(Float64, 2) && gali_floor(Double64, 4) == ga_meas(Double64, 4)
end

@testset "AC1c negative control: skip-normalisation double" begin
    # A double that evaluates the volume of the raw columns (no normalisation) fails the bound of AC1(b) at the spread sets, for nearly every set.
    rng = StableRNG(ga_seed(1, 1))
    fail, total = 0, 0
    for (d, K) in GA_PAIRS, _ in 1:20, k in 2:K
        W = ga_randset(rng, Float64, d, K)
        ex = ga_exact(W)
        fail += ga_err(ga_ctl_noscale(W, k), ex[k - 1]) > 2 * ga_meas(Float64, k) ? 1 : 0
        total += 1
    end
    @test fail >= 0.9 * total
end

@testset "AC2 closed forms by prescribed Gram matrix" begin
    _ = (gali_values,)
    # Unit-diagonal Gram matrices, vectors from a BigFloat(256) Cholesky factor embedded in d in {k, 4, 8} by a random orthogonal map, rounded to T;
    # pass |G_T - sqrt(det Gram)| <= 3 B_meas(T, k). (a) k = 2, theta = 1e-m and pi/2 - 1e-m (m = 1..14); (b) k = 3, equal cosine c in {0, 0.5, -0.4, 0.9}
    # and 1 - 1e-m (m = 1..7), and c = -1/2 (coplanar, closed form 0); (c) k = 4, c in {0, 0.3, -0.3}; (d) orthonormal sets give 1, k = d random sets give
    # the absolute determinant of the unit matrix. The cosines are the Float64 numbers, widened exactly.
    rng = StableRNG(ga_seed(2))
    for T in GA_TYPES
        worst = Ref(0.0)
        function chk(Gram, k, d, closed)
            local W, ex
            setprecision(BigFloat, GA_BITS) do
                V = ga_gram_vectors(Gram, d, ga_randorth(rng, d))
                W = ga_sm([ga_to(T, V[i, j]) for i in 1:d, j in 1:k])
                ex = closed(Gram)
            end
            g = gali_values(W)[k - 1]
            worst[] = max(worst[], ga_err(g, ex) / (3 * ga_meas(T, k)))
        end
        eqgram(k, c) = (G = fill(c, k, k); for i in 1:k; G[i, i] = one(c); end; G)
        for m in 1:14, d in (2, 4, 8)
            th = setprecision(() -> BigFloat(10.0)^(-m), BigFloat, GA_BITS)
            for ang in (th, setprecision(() -> BigFloat(pi) / 2 - BigFloat(10.0)^(-m), BigFloat, GA_BITS))
                setprecision(BigFloat, GA_BITS) do
                    c = cos(ang)
                    chk(BigFloat[1 c; c 1], 2, d, G -> sin(ang))
                end
            end
        end
        cs3 = setprecision(BigFloat, GA_BITS) do
            vcat(BigFloat.([0.0, 0.5, -0.4, 0.9, -0.5]), [1 - BigFloat(10.0)^(-m) for m in 1:7])
        end
        for c in cs3, d in (3, 4, 8)
            setprecision(BigFloat, GA_BITS) do
                chk(eqgram(3, c), 3, d, G -> sqrt(max((1 + 2c) * (1 - c)^2, zero(c))))
            end
        end
        for c in setprecision(() -> BigFloat.([0.0, 0.3, -0.3]), BigFloat, GA_BITS), d in (4, 8)
            setprecision(BigFloat, GA_BITS) do
                chk(eqgram(4, c), 4, d, G -> sqrt((1 + 3c) * (1 - c)^3))
            end
        end
        for (d, k) in ((4, 4), (6, 6), (8, 8))                     # (d) k = d random sets: the absolute determinant of the unit matrix
            for _ in 1:5
                W = ga_randset(rng, T, d, k)
                ex = ga_exact(W)
                worst[] = max(worst[], ga_err(gali_values(W)[k - 1], ex[k - 1]) / (3 * ga_meas(T, k)))
            end
        end
        @test worst[] <= 1                                         # in units of 3 B_meas(T, k)
    end
end

@testset "AC3 algebra and invariants of the definition" begin
    _ = (gali_values,)
    # The random sets of AC1 (D38: 30 per pair and type): (a) values in [0, 1 + B_meas]; (b) nesting G_{k+1} <= G_k + 2 B_meas(T, k+1); (c) the value of
    # order k changes by at most 2 B_meas under a column sign flip, a permutation of the first k columns and a random orthogonal map of the ambient space;
    # (d) a last column equal, opposite or a random combination of the others gives at most 2 B_meas(T, K); (e) columns multiplied by 2^+-498 (about 1e+-150,
    # a power of two so the scaling itself is exact) give the unscaled value within 2 B_meas and no error (Float64).
    rng = StableRNG(ga_seed(3))
    for T in GA_TYPES
        bad = (range = 0, nest = 0, sign = 0, perm = 0, orth = 0, dep = 0, scale = 0)
        for (d, K) in GA_PAIRS, _ in 1:ga_size(30, 300)
            W = ga_randset(rng, T, d, K)
            g = gali_values(W)
            bad = merge(bad, (range = bad.range + count(x -> !(0 <= x <= 1 + ga_meas(T, K)), g),))
            bad = merge(bad, (nest = bad.nest + count(k -> g[k] > g[k - 1] + 2 * ga_meas(T, k + 1), 2:(K - 1)),))
            i = rand(rng, 1:K)
            Wf = ga_sm(hcat((j == i ? -W[:, j] : W[:, j] for j in 1:K)...))
            gf = gali_values(Wf)
            bad = merge(bad, (sign = bad.sign + count(k -> abs(gf[k - 1] - g[k - 1]) > 2 * ga_meas(T, k), 2:K),))
            nperm = 0
            for k in 2:K
                p = vcat(sortperm(rand(rng, k)), (k + 1):K)
                nperm += abs(gali_values(ga_sm(W[:, p]))[k - 1] - g[k - 1]) > 2 * ga_meas(T, k) ? 1 : 0
            end
            bad = merge(bad, (perm = bad.perm + nperm,))
            Wo = setprecision(BigFloat, GA_BITS) do
                Q = ga_randorth(rng, d)
                B = Q * BigFloat[ga_big(W[a, b]) for a in 1:d, b in 1:K]
                ga_sm([ga_to(T, B[a, b]) for a in 1:d, b in 1:K])
            end
            go = gali_values(Wo)
            bad = merge(bad, (orth = bad.orth + count(k -> abs(go[k - 1] - g[k - 1]) > 2 * ga_meas(T, k), 2:K),))
            kind = rand(rng, 1:3)
            coef = kind == 1 ? [1.0] : kind == 2 ? [-2.5] : [2 * rand(rng) - 1 for _ in 1:(K - 1)]
            last = sum(T(coef[min(j, length(coef))]) * W[:, j] for j in 1:(kind == 3 ? K - 1 : 1))
            Wd = ga_sm(hcat(W[:, 1:(K - 1)], last))
            bad = merge(bad, (dep = bad.dep + (gali_values(Wd)[end] > 2 * ga_meas(T, K) ? 1 : 0),))
        end
        for (d, K) in GA_PAIRS, _ in 1:ga_size(10, 100)
            T === Float64 || break
            W = ga_randset(rng, Float64, d, K; spread = false)
            sc = [ldexp(1.0, rand(rng, Bool) ? 498 : -498) for _ in 1:K]
            Ws = ga_sm(hcat((W[:, j] * sc[j] for j in 1:K)...))
            g, gs = gali_values(W), gali_values(Ws)
            bad = merge(bad, (scale = bad.scale + count(k -> !isfinite(gs[k - 1]) || abs(gs[k - 1] - g[k - 1]) > 2 * ga_meas(T, k), 2:K),))
        end
        @test bad == (range = 0, nest = 0, sign = 0, perm = 0, orth = 0, dep = 0, scale = 0)
    end
end

# ---- AC4: k = 2 against SALI ----
function ga_pair(rng, ::Type{T}, d::Int) where {T}
    a = ga_randvec(rng, T, d) * T(10.0^(6 * rand(rng) - 3))
    if iseven(rand(rng, 1:2))
        b = ga_randvec(rng, T, d) * T(10.0^(6 * rand(rng) - 3))
    else
        s = T(10.0^(4 * rand(rng) - 2))
        b = a * s + (norm(a) * s * T(10.0^(-1 - 13 * rand(rng)))) * ga_randvec(rng, T, d)    # near-parallel: perturbation 10^-u, u in [1, 14]
    end
    rand(rng, Bool) && (b = -b)
    return a, b
end
function ga_constructed_pairs(::Type{T}, d::Int) where {T}
    out = Tuple{Vector{T}, Vector{T}}[]
    for m in 1:14, sgn in (1, -1), base in (T(10)^(-m), T(pi) / 2 - T(10)^(-m))
        u1 = vcat(T[cos(T(0.3)), sin(T(0.3))], zeros(T, d - 2))
        u2 = sgn * vcat(T[cos(T(0.3) + base), sin(T(0.3) + base)], zeros(T, d - 2))
        push!(out, (u1, u2))
    end
    return out
end

@testset "AC4 k = 2 against SALI and the conversions" begin
    _ = (gali_values, sali_pair, alignment_indices, sali_floor)
    # d in {2, 3, 4, 8}, random pairs (half independent, half near-parallel, random signs) and constructed pairs at angles 1e-m, pi/2 - 1e-m (m = 1..14) and
    # antiparallel. With g = G_2 of the module and SALI, d-, d+ of the sali module, the input SALI carrying its own floor 4 eps (sali_floor), the tolerance of
    # every comparison is 2 B_meas(T, 2) + sali_floor(T) (the spec states 2 B_meas; the second term is the error of the SALI input): (a) |g - SALI max(d-, d+)/2|;
    # (b) |g - SALI sqrt(4 - SALI^2)/2|; (c) the inverse conversion sqrt(2) g / sqrt(1 + sqrt(1 - g^2)) against SALI within the tolerance / cos(phi) for g <= 0.9;
    # (d) SALI/sqrt 2 <= g <= SALI within the tolerance, and SALI <= 1.0013 g + 2 tolerance for 10 B_meas <= g <= 0.1 (the relative form is untestable at
    # g = 10 B_meas where g itself is uncertain by 10 percent); (e) the threshold map: signs of SALI - tau_S and g - tau_g agree except within the tolerance of the
    # threshold, tau_S in {0.1, 1e-3, 1e-6}, with pairs constructed at SALI = tau_S (1 + delta), delta in +-{1e-2, 1e-6, 1e-10}; the count of disagreements is
    # reported. D38: 300 random pairs per d and type (spec: 1e4).
    rng = StableRNG(ga_seed(4))
    for T in GA_TYPES
        tol = 2 * ga_meas(T, 2) + Float64(sali_floor(T))
        worst = zeros(4)
        bad = 0
        count_sign, disagree = 0, 0
        for d in (2, 3, 4, 8)
            pairs = vcat([ga_pair(rng, T, d) for _ in 1:ga_size(300, 10_000)], ga_constructed_pairs(T, d))
            for tS in (0.1, 1e-3, 1e-6), delta in (1e-2, -1e-2, 1e-6, -1e-6, 1e-10, -1e-10)
                phi = 2 * asin(T(tS * (1 + delta)) / 2)                   # SALI of the pair = tS (1 + delta)
                push!(pairs, (vcat(T[1, 0], zeros(T, d - 2)), vcat(T[cos(phi), sin(phi)], zeros(T, d - 2))))
            end
            for (a, b) in pairs
                g = gali_values(ga_sm(hcat(a, b)))[1]
                dm, dp = alignment_indices(a, b)
                S = sali_pair(a, b)
                worst[1] = max(worst[1], Float64(abs(g - S * max(dm, dp) / 2)) / tol)
                worst[2] = max(worst[2], Float64(abs(g - S * sqrt(4 - S^2) / 2)) / tol)
                if g <= 0.9
                    cphi = Float64(sqrt(1 - g^2))
                    worst[3] = max(worst[3], Float64(abs(sqrt(T(2)) * g / sqrt(1 + sqrt(1 - g^2)) - S)) / (tol / cphi))
                end
                bad += (g <= S + tol && g >= S / sqrt(T(2)) - tol) ? 0 : 1
                if 10 * ga_meas(T, 2) <= g <= 0.1
                    bad += S <= 1.0013 * g + 2 * tol ? 0 : 1
                end
                for tS in (0.1, 1e-3, 1e-6)
                    tg = tS * sqrt(4 - tS^2) / 2
                    count_sign += 1
                    if sign(S - tS) != sign(g - tg) && abs(Float64(g) - tg) > tol
                        disagree += 1
                    end
                end
            end
        end
        @test worst[1] <= 1 && worst[2] <= 1 && worst[3] <= 1
        @test bad == 0
        @test disagree == 0
        @info "AC4 reported" T worst count_sign
    end
end

@testset "AC4f negative control: SALI returned as GALI_2" begin
    # A double that returns SALI as GALI_2 fails the conversion check (b) by far (the gap is SALI^3/8) for pairs with SALI >= 1e-3.
    rng = StableRNG(ga_seed(4, 1))
    fail, total = 0, 0
    for _ in 1:200
        a, b = ga_randvec(rng, Float64, 4), ga_randvec(rng, Float64, 4)
        S = ga_sali(a, b)
        S >= 1e-3 || continue
        total += 1
        fail += abs(S - S * sqrt(4 - S^2) / 2) > 2 * ga_meas(Float64, 2) + 4 * ga_eps(Float64) ? 1 : 0
    end
    @test total > 100 && fail == total
end
