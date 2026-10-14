# AC1 to AC3: factor identities on seeded synthetic matrices, near-parallel columns, variant agreement.
# Spec: lyap-householder-qr-spec.md, Acceptance criteria 1 to 3. Orthonormalize is the pure form of the API contract.

const HH_NS = 200   # seeded matrices per (d, k, T)

hh_panel(::Type{T}, vd, vk, seed) where {T} = (rng = StableRNG(seed); [hh_randW(rng, T, vd, vk) for _ in 1:HH_NS])

# Exact (256-bit) evaluation of the three AC1 ratios, each already divided by its tolerance (pass: <= 1).
function hh_eval(W::SMatrix{d, k, T}, Q, R, Rref, eps_T) where {d, k, T}
    return setprecision(BigFloat, 256) do
        Wb, Qb, Rb = hh_bigmat(W), hh_bigmat(Q), hh_bigmat(R)
        res = Qb * Rb - Wb
        ra = maximum(Float64(sqrt(sum(abs2, res[:, j]))) / (8 * d * eps_T * hh_colnorm(W, j)) for j in 1:k)
        rb = Float64(maximum(abs, Qb' * Qb - I)) / (8 * k * eps_T)
        rd = maximum(Float64(abs(Rb[i, i] - Rref[i, i])) / (8 * d * eps_T * hh_kappa(W, i - 1) * hh_colnorm(W, i)) for i in 1:k)
        (ra, rb, rd)
    end
end

function hh_ac1(::Type{T}, vd::Val{d}, vk::Val{k}, seed) where {T, d, k}
    eps_T = Float64(effective_eps(T))
    algs = hh_algs(T)
    worst = [(0.0, 0.0, 0.0) for _ in algs]
    structural = zeros(Int, length(algs))   # status, shapes, zeros below the diagonal, positive diagonal, counts
    for W in hh_panel(T, vd, vk, seed)
        Rref = setprecision(BigFloat, 256) do
            hh_ref_qr(W)[2]
        end
        for (ia, alg) in enumerate(algs)
            Q, R, cnt, status, col = orthonormalize(alg, W)
            ok = status == QR_OK && Q isa SMatrix{d, k, T} && R isa SMatrix{k, k, T} && cnt == hh_count(alg, d, k) &&
                 all(R[i, j] == 0 for j in 1:k for i in (j + 1):k) && all(R[i, i] > 0 for i in 1:k)
            structural[ia] += ok ? 0 : 1
            worst[ia] = map(max, worst[ia], hh_eval(W, Q, R, Rref, eps_T))
        end
    end
    for (ia, alg) in enumerate(algs)
        @test structural[ia] == 0                 # (c) exact zeros below the diagonal, R_ii > 0; (e) counts r(d,k), g(d,k); status
        @test worst[ia][1] <= 1                   # (a) column residual <= 8 d eps ||W e_j||
        @test worst[ia][2] <= 1                   # (b) max |Q'Q - I| <= 8 k eps
        @test worst[ia][3] <= 1                   # (d) |R_ii - R_ii^ref| <= 8 d eps kappa_{i-1} ||W e_i||
    end
end

@testset "AC1 factor identities" begin
    # (d,k) in the seven pairs, 200 seeded W each, column scales 10^u; V-H library (Float64 only), V-H, V-A, V-G kernels; both T
    for (it, T) in enumerate((Float64, Double64)), (idk, (vd, vk)) in enumerate(HH_DK)
        hh_ac1(T, vd, vk, hh_seed(1, 10 * it + idk))
    end
end

@testset "AC2 near-parallel columns" begin
    # d = 4, k = 2, w2 = w1 + theta w1_perp rounded to T: loss <= 8 k eps for every theta; |R_22 - ref| <= 8 d eps ||w2|| (E10, no 1/theta growth)
    for T in (Float64, Double64), theta in (2.0^-10, 2.0^-20, 2.0^-40)
        eps_T = Float64(effective_eps(T))
        s, c = sincos(T(0.3))
        w1 = SVector{4, T}(c, s, 0, 0)
        w2 = w1 + T(theta) * SVector{4, T}(-s, c, 0, 0)
        W = hcat(w1, w2)
        R22ref = setprecision(BigFloat, 256) do
            hh_ref_qr(W)[2][2, 2]
        end
        for alg in hh_algs(T)
            Q, R, cnt, status, col = orthonormalize(alg, W)
            @test status == QR_OK
            loss, err = setprecision(BigFloat, 256) do
                Qb = hh_bigmat(Q)
                (Float64(maximum(abs, Qb' * Qb - I)), Float64(abs(hh_big(R[2, 2]) - R22ref)))
            end
            @test loss <= 8 * 2 * eps_T
            @test err <= 8 * 4 * eps_T * hh_colnorm(W, 2)
        end
    end
end

function hh_ac3(::Type{T}, ::Val{d}, ::Val{k}, seed) where {T, d, k}
    eps_T = Float64(effective_eps(T))
    algs = hh_algs(T)
    npair = length(algs) * (length(algs) - 1) ÷ 2
    worstR, worstQ = zeros(npair), zeros(npair)
    for W in hh_panel(T, Val(d), Val(k), seed)
        res = [orthonormalize(alg, W) for alg in algs]
        p = 0
        for a in 1:length(algs), b in (a + 1):length(algs)
            p += 1
            Qa, Ra = res[a][1], res[a][2]
            Qb, Rb = res[b][1], res[b][2]
            rR = maximum(Float64(abs(Ra[i, i] - Rb[i, i])) / (16 * d * eps_T * hh_kappa(W, i - 1) * hh_colnorm(W, i)) for i in 1:k)
            rQ = maximum(Float64(maximum(abs, Qa[:, i] - Qb[:, i])) / (16 * d * eps_T * hh_kappa(W, i)) for i in 1:k)
            worstR[p] = max(worstR[p], rR)
            worstQ[p] = max(worstQ[p], rQ)
        end
    end
    @test all(<=(1), worstR)
    @test all(<=(1), worstQ)
end

@testset "AC3 variant agreement" begin
    # data of AC1; pairs of variants: max_i |R_ii^A - R_ii^B| <= 16 d eps kappa_{i-1} ||W e_i||, max |Q^A e_i - Q^B e_i| <= 16 d eps kappa_i
    for (it, T) in enumerate((Float64, Double64)), (idk, (vd, vk)) in enumerate(HH_DK)
        hh_ac3(T, vd, vk, hh_seed(1, 10 * it + idk))
    end
end
