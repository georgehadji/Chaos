# Stage 04 verification of module lyap-householder-qr (HouseholderQR, AugmentedQR, GivensQR, LapackQR).
# Run from the repo root:  julia +release -t 2 --project=verify verify/lyap-householder-qr.jl
# Spec: workspaces/module-build/stages/01-spec/output/lyap-householder-qr-spec.md (criteria 1 to 17, Amendment D35); package suite: test/lyap-householder-qr/.
# Section 1 ("cross"): the full cross-product of the package-suite testsets (decision D38, HH_FULL_RUN = true, 2035 assertions, about 21 minutes).
# The other sections print one line per criterion:  AC<n> | quantity | value | bound | PASS/FAIL (value = the measured number; ratios are error / bound, pass at <= 1).
# The measurement code reuses the helpers of the package suite (hh_helpers.jl and the definitions in the hh_*.jl files, loaded without their testsets).
# Seeds: HH_SEED + 44000 + 50 ac + sub (the package suite uses HH_SEED + 24000 .. 24999; verify/ scripts of other slugs use other blocks).
# HV_ONLY="AC7,D35" runs only the named sections (cross, AC1-3, AC4-6, AC7, AC8, AC9-12, AC13, AC15, alg2, dev, D35, AC17, R, figure);
# HV_LOG sets the log file; HV_N6 (default 1000000) and HV_N6_REF (default 2) set the AC7 long-horizon run; HV_AC7_STARTS (default 50); HV_D35_PER_K (default 3000).
using Printf, Random, LinearAlgebra, Statistics, Test
using Chaoslab, StaticArrays, DoubleFloats, StableRNGs
import ChaosTools as CT
import CairoMakie as CM
import BenchmarkTools: @belapsed

const HH_FULL_RUN = true
const ROOT = dirname(@__DIR__)
const OUT = joinpath(ROOT, "workspaces", "module-build", "stages", "04-verify", "output")
const TESTDIR = joinpath(ROOT, "test", "lyap-householder-qr")
mkpath(OUT)
const ONLY = filter(!isempty, split(get(ENV, "HV_ONLY", ""), ","))
const LOG = open(get(ENV, "HV_LOG", joinpath(OUT, "lyap-householder-qr-verify-log.txt")), "w")
const N6 = parse(Int, get(ENV, "HV_N6", "1000000"))
const N6_REF = parse(Int, get(ENV, "HV_N6_REF", "2"))
const AC7_STARTS = parse(Int, get(ENV, "HV_AC7_STARTS", "50"))
const D35_PER_K = parse(Int, get(ENV, "HV_D35_PER_K", "3000"))

fmt(x::AbstractFloat) = @sprintf("%.4g", x)
fmt(x) = string(x)
function say(s)
    println(s); println(LOG, s); flush(stdout); flush(LOG)
end
emit(ac, quantity, value, bound, verdict) = say(string(ac, " | ", quantity, " | ", fmt(value), " | ", bound, " | ", verdict))
pf(ok::Bool) = ok ? "PASS" : "FAIL"
ratio_line(ac, quantity, r) = emit(ac, quantity, r, "<= 1", pf(r <= 1))
count_line(ac, quantity, c) = emit(ac, quantity, c, "0", pf(c == 0))
selected(name) = isempty(ONLY) || name in ONLY
vseed(ac, sub = 0) = 20261012 + 44000 + 50 * ac + sub

function section(f, name)
    selected(name) || return
    t0 = time()
    try
        f()
    catch e
        say(string(name, " | ERROR | ", first(sprint(showerror, e), 400), " | - | FAIL"))
    end
    say(@sprintf("# section %s took %.1f s", name, time() - t0))
end

# Definitions of the package suite without its testsets (one home for helpers, references and bounds).
function include_defs(path)
    ex = Meta.parseall(read(path, String); filename = path)
    for e in ex.args
        e isa LineNumberNode && continue
        e isa Expr && e.head === :macrocall && e.args[1] === Symbol("@testset") && continue
        Core.eval(Main, e)
    end
end
load_defs() = foreach(f -> include_defs(joinpath(TESTDIR, f)), ("hh_helpers.jl", "hh_factor.jl", "hh_exact.jl", "hh_panel.jl", "hh_cross.jl", "hh_runtime.jl", "hh_backend.jl"))

const T0 = time()
say("# lyap-householder-qr verification, julia $(VERSION), threads $(Threads.nthreads()), base seed 20261012 + 44000 (+ 50 ac + sub), time $(round(time())) (unix)")
say("# timings in this log were taken on a loaded machine (other agents were running concurrently); runtime is informative only")

# ---------------------------------------------------------------------------------------------------
# Section 1: full cross-product of the package-suite testsets (D38)
if selected("cross")
    tc = time()
    res = try
        @testset verbose = true "lyap-householder-qr full cross-product" begin
            include(joinpath(TESTDIR, "runtests.jl"))
        end
    catch e
        e
    end
    if res isa Test.TestSetException
        emit("CROSS", "full cross-product (D38) assertions pass / fail / error / broken", "$(res.pass) / $(res.fail) / $(res.error) / $(res.broken)", "fail = error = 0", pf(res.fail == 0 && res.error == 0))
    elseif res isa Test.AbstractTestSet
        c = Test.get_test_counts(res)
        emit("CROSS", "full cross-product (D38) assertions pass / fail / error / broken", "$(c.cumulative_passes) / $(c.cumulative_fails) / $(c.cumulative_errors) / $(c.cumulative_broken)", "fail = error = 0", pf(c.cumulative_fails == 0 && c.cumulative_errors == 0))
    else
        say(string("CROSS | ERROR | ", first(sprint(showerror, res), 400), " | - | FAIL"))
    end
    say(@sprintf("# section cross took %.1f s", time() - tc))
end
@isdefined(hh_eval) || load_defs()

const Lyr = Chaoslab.Lyapunov
const RES = Dict{Symbol, Any}()        # numbers shared with the figures and the oracle summary
const eps64 = Float64(effective_eps(Float64))

# ---------------------------------------------------------------------------------------------------
# AC1 to AC3: factor identities (200 seeded W per shape, as the spec asks), near-parallel columns, variant agreement
function ac13_shape(::Type{T}, vd::Val{d}, vk::Val{k}, seed) where {T, d, k}
    eps_T = Float64(effective_eps(T))
    algs = hh_algs(T)
    na = length(algs)
    rng = StableRNG(seed)
    worst = zeros(3, na); structural = zeros(Int, na)
    prs = [(a, b) for a in 1:na for b in (a + 1):na]
    wR, wQ = zeros(length(prs)), zeros(length(prs))
    for _ in 1:200
        W = hh_randW(rng, T, vd, vk)
        Rref = setprecision(BigFloat, 256) do
            hh_ref_qr(W)[2]
        end
        res = [orthonormalize(alg, W) for alg in algs]
        for ia in 1:na
            Q, R, cnt, status, col = res[ia]
            ok = status == QR_OK && Q isa SMatrix{d, k, T} && R isa SMatrix{k, k, T} && cnt == hh_count(algs[ia], d, k) &&
                 all(R[i, j] == 0 for j in 1:k for i in (j + 1):k) && all(R[i, i] > 0 for i in 1:k)
            structural[ia] += ok ? 0 : 1
            worst[:, ia] .= max.(worst[:, ia], collect(hh_eval(W, Q, R, Rref, eps_T)))
        end
        for (p, (a, b)) in enumerate(prs)
            Qa, Ra, Qb, Rb = res[a][1], res[a][2], res[b][1], res[b][2]
            wR[p] = max(wR[p], maximum(Float64(abs(Ra[i, i] - Rb[i, i])) / (16 * d * eps_T * hh_kappa(W, i - 1) * hh_colnorm(W, i)) for i in 1:k))
            wQ[p] = max(wQ[p], maximum(Float64(maximum(abs, Qa[:, i] - Qb[:, i])) / (16 * d * eps_T * hh_kappa(W, i)) for i in 1:k))
        end
    end
    return (algs = algs, worst = worst, structural = structural, prs = prs, wR = wR, wQ = wQ)
end

function ac2_rows(::Type{T}) where {T}
    eps_T = Float64(effective_eps(T))
    rows = []
    for theta in (2.0^-10, 2.0^-20, 2.0^-40)
        s, c = sincos(T(0.3))
        w1 = SVector{4, T}(c, s, 0, 0)
        w2 = w1 + T(theta) * SVector{4, T}(-s, c, 0, 0)
        W = hcat(w1, w2)
        R22ref = setprecision(BigFloat, 256) do
            hh_ref_qr(W)[2][2, 2]
        end
        loss, err = 0.0, 0.0
        for alg in hh_algs(T)
            Q, R, cnt, status, col = orthonormalize(alg, W)
            status == QR_OK || (loss = Inf)
            l, e = setprecision(BigFloat, 256) do
                Qb = hh_bigmat(Q)
                (Float64(maximum(abs, Qb' * Qb - I)), Float64(abs(hh_big(R[2, 2]) - R22ref)))
            end
            loss, err = max(loss, l), max(err, e)
        end
        push!(rows, (theta, loss / eps_T, err / (eps_T * hh_colnorm(W, 2)), 8 * 2, 8 * 4))
    end
    return rows
end

section("AC1-3") do
    fig = Dict()
    sa = Dict(Float64 => 0.0, Double64 => 0.0); sd = copy(sa); s3R = copy(sa); s3Q = copy(sa); sc = Dict(Float64 => 0, Double64 => 0)
    for (it, T) in enumerate((Float64, Double64)), (idk, (vd, vk)) in enumerate(HH_DK)
        r = ac13_shape(T, vd, vk, vseed(1, 10 * it + idk))
        d, k = typeof(vd).parameters[1], typeof(vk).parameters[1]
        rb = vec(r.worst[2, :])
        ratio_line("AC1", "(b) max|Q'Q - I| / (8 k eps), worst over 200 W, $T, (d,k) = ($d,$k); per variant " *
                          join(["$(hh_name(a)) $(fmt(x))" for (a, x) in zip(r.algs, rb)], ", ") * "; realized loss of the worst variant $(fmt(8 * maximum(rb))) k eps", maximum(rb))
        for (ia, alg) in enumerate(r.algs)
            fig[(T, idk, hh_name(alg))] = Tuple(r.worst[:, ia])
        end
        sa[T] = max(sa[T], maximum(r.worst[1, :])); sd[T] = max(sd[T], maximum(r.worst[3, :])); sc[T] += sum(r.structural)
        s3R[T] = max(s3R[T], maximum(r.wR)); s3Q[T] = max(s3Q[T], maximum(r.wQ))
        say("#   AC1 $T ($d,$k): (a) per variant " * join([fmt(x) for x in r.worst[1, :]], ", ") * "; (d) per variant " * join([fmt(x) for x in r.worst[3, :]], ", ") *
            "; AC3 pairs (R, Q) worst " * fmt(maximum(r.wR)) * ", " * fmt(maximum(r.wQ)))
    end
    RES[:ac1] = fig
    for T in (Float64, Double64)
        ratio_line("AC1", "(a) column residual ||W e_j - Q R e_j|| / (8 d eps ||W e_j||), worst over 7 shapes x 200 W x variants, $T", sa[T])
        ratio_line("AC1", "(d) |R_ii - R_ii^ref(256 bit)| / (8 d eps kappa_(i-1) ||W e_i||), worst, $T", sd[T])
        count_line("AC1", "(c)+(e) failures of: status QR_OK, shapes, exact zeros below diagonal, R_ii > 0, counts r(d,k) and g(d,k) (7 shapes x 200 W x variants), $T", sc[T])
        ratio_line("AC3", "variant agreement, R_ii: max_i |R_ii^A - R_ii^B| / (16 d eps kappa_(i-1) ||W e_i||), all variant pairs, $T", s3R[T])
        ratio_line("AC3", "variant agreement, Q: max|Q^A e_i - Q^B e_i| / (16 d eps kappa_i), all variant pairs, $T", s3Q[T])
    end
    for T in (Float64, Double64), (theta, loss, err, bl, be) in ac2_rows(T)
        ratio_line("AC2", "(a) max|Q'Q - I| / (8 k eps) at theta = 2^$(Int(round(log2(theta)))), all variants, $T (loss $(fmt(loss)) eps)", loss / bl)
        ratio_line("AC2", "(b) |R_22 - ref| / (8 d eps ||w_2||) at theta = 2^$(Int(round(log2(theta)))), all variants, $T (error $(fmt(err)) eps ||w_2||)", err / be)
    end
    RES[:ac2] = ac2_rows(Float64)
end

function y5_ratios(::Type{T}, mus, nonnormal::Bool, ::Val{d}, seed) where {T, d}
    N = d ÷ 2
    A0 = block_diagonal(map(m -> hyperbolic_block(T, m), mus)...)
    A = nonnormal ? conjugated(A0, N == 2 ? SMatrix{2, 2, T}(1, 2, 2, 1) : SMatrix{3, 3, T}(1, 2, 0, 2, 1, 1, 0, 1, 1)) : A0
    sys = LinearSymplecticMap(A)
    Q0 = initial_deviations(seed, 1, Val(d), Val(d), T)
    n = 500
    Af = Float64[Float64(A[i, j]) for i in 1:d, j in 1:d]
    sv = svdvals(Af)
    theta = hh_theta(T, sv[1] / sv[end], 1, maximum(abs, Af))
    lnvol, Lam, Bterm = setprecision(BigFloat, 2048) do
        Y = hh_matpow(hh_bigmat(A), n) * hh_bigmat(Q0)
        lv = [log(hh_det(Y[:, 1:i]' * Y[:, 1:i])) / 2 for i in 1:d]
        lst = Tuple{BigFloat, Int, Int}[]
        for (j, m) in enumerate(mus)
            push!(lst, (BigFloat(m), j, 1)); push!(lst, (-BigFloat(m), j, -1))
        end
        sort!(lst; by = t -> -t[1])
        V = zeros(BigFloat, d, d)
        for (c, (_, j, sg)) in enumerate(lst)
            V[j, c] = 1 / sqrt(BigFloat(2)); V[N + j, c] = sg / sqrt(BigFloat(2))
        end
        Q0b = hh_bigmat(Q0)
        (lv, cumsum([t[1] for t in lst]), [log(abs(hh_det(V[:, 1:i]' * Q0b[:, 1:i]))) for i in 1:d])
    end
    ra, rb, rc = 0.0, 0.0, 0.0
    for alg in hh_algs(T)
        f = propagate(sys, [zero(SVector{d, T})], n; k = d, tau = 1, seed = seed, W0 = [Q0], reducers = (QRSpectrum(alg),)).finals[1][1]
        setprecision(BigFloat, 2048) do
            C = cumsum([hh_big(f.sum[j]) + hh_big(f.comp[j]) for j in 1:d])
            for i in 1:d
                ra = max(ra, Float64(abs(C[i] - lnvol[i])) / (i * n * theta))
                nonnormal || (rb = max(rb, Float64(abs(C[i] / n - Lam[i] - Bterm[i] / n)) / (i * theta)))
            end
        end
        rc = max(rc, Float64(abs(sum(f.lambda))) / theta)
    end
    return ra, rb, rc
end

# ---------------------------------------------------------------------------------------------------
# AC4 to AC6: K = 0, cat map, Y5
section("AC4-6") do
    # AC4(a)
    for T in (Float64, Double64)
        bad = 0
        for alg in hh_algs(T)
            W0 = SMatrix{2, 2, T}(1, 0, 0, 1)
            shear = LinearSymplecticMap(SMatrix{2, 2, T}(1, 0, 1, 1))
            res = propagate(shear, [zero(SVector{2, T})], 10_000; k = 2, tau = 1, seed = vseed(4), W0 = [W0], reducers = (QRSpectrum(alg), hh_Dev()))
            f, dv = res.finals[1]
            bad += !(all(iszero, f.sum) && all(iszero, f.comp) && all(iszero, f.lambda) && dv.maxdev == 0 && dv.count == 10_000 && res.W[1] == W0)
        end
        count_line("AC4", "(a) shear, W_0 = I, tau = 1, n = 1e4: variants (of $(length(hh_algs(T)))) with any nonzero S_i, lambda or |W - I|, $T", bad)
    end
    # AC4(b)
    for T in (Float64, Double64)
        eps_T = Float64(effective_eps(T))
        for n in (10, 1000, 10_000)
            w, lnb = 0.0, 0.0
            for alg in hh_algs(T)
                f = propagate(StandardMap(T(0)), [SVector{2, T}(T(0.25), T(0.5))], n; k = 1, tau = 1, seed = vseed(4), W0 = [SMatrix{2, 1, T}(0, 1)], reducers = (QRSpectrum(alg),)).finals[1][1]
                S = f.sum[1] + f.comp[1]
                eb = setprecision(BigFloat, 256) do
                    Float64(abs(hh_big(S) - log(BigFloat(n)^2 + 1) / 2))
                end
                w = max(w, eb / (8 * n * eps_T))
                lnb = max(lnb, Float64(S) / n - log(n) / n)
            end
            ratio_line("AC4", "(b) |S_1(n) - (1/2) ln(n^2+1)| / (8 n eps), K = 0, k = 1, n = $n, all variants, $T", w)
            n == 10_000 && emit("AC4", "(b) max S_1/n - ln(n)/n at n = 1e4, $T (E12)", lnb, "<= 1e-12", pf(lnb <= 1e-12))
        end
    end
    # AC5
    lnmu, lnc, r0, q = setprecision(BigFloat, 256) do
        mub = (3 + sqrt(BigFloat(5))) / 2
        cc = 1 / sqrt(1 + ((sqrt(BigFloat(5)) - 1) / 2)^2)
        (log(mub), log(cc), (1 - cc^2) / cc^2, mub^-4)
    end
    mu = Float64((3 + sqrt(5)) / 2)
    ref1 = Dict(n => setprecision(BigFloat, 256) do
        log(sqrt(BigFloat(hh_fibbig(2n + 1))^2 + BigFloat(hh_fibbig(2n))^2))
    end for n in (10, 100, 1000, 10_000))
    for T in (Float64, Double64)
        eps_T = Float64(effective_eps(T))
        ra, rbS, rc, rd, obs = 0.0, 0.0, 0.0, 0.0, 0.0
        for tau in (1, 2, 5), n in (10, 100, 1000, 10_000), alg in hh_algs(T)
            sys = LinearSymplecticMap(SMatrix{2, 2, T}(2, 1, 1, 1))
            f = propagate(sys, [zero(SVector{2, T})], n; k = 2, tau = tau, seed = vseed(5), W0 = [SMatrix{2, 2, T}(1, 0, 0, 1)], reducers = (QRSpectrum(alg),)).finals[1][1]
            S1, S2 = f.sum[1] + f.comp[1], f.sum[2] + f.comp[2]
            Gfirst = Float64(hh_fibbig(2tau + 1)^2 + 2 * hh_fibbig(2tau)^2 + hh_fibbig(2tau - 1)^2)
            BS = hh_BS(T, n, Gfirst)
            theta = hh_theta(T, mu^(2tau), tau, 2.0)
            e1, e2 = setprecision(BigFloat, 256) do
                (Float64(abs(hh_big(S1) - ref1[n])), Float64(abs(hh_big(f.lambda[1]) - (lnmu + (lnc + log(1 + r0 * q^n) / 2) / n))))
            end
            ra = max(ra, e1 / BS); rbS = max(rbS, Float64(abs(S1 + S2)) / (n * theta)); rc = max(rc, e2 / (BS / n)); rd = max(rd, Float64(abs(f.lambda[1] + f.lambda[2])) / theta)
            tau == 1 && (obs = max(obs, e1 / (n * eps_T)))
        end
        ratio_line("AC5", "(a) |S_1 - ln||A^n e_1|| (exact Fibonacci)| / B_S, tau in {1,2,5}, n in {10..1e4}, all variants, $T (tau = 1: max $(fmt(obs)) n eps; spec observed 2.03)", ra)
        ratio_line("AC5", "(b) |S_1 + S_2| / (n theta_sum), same grid, $T", rbS)
        ratio_line("AC5", "(c) |lambda_1 - exact E10 form| / (B_S / n), same grid, $T", rc)
        ratio_line("AC5", "(d) |lambda_1 + lambda_2| / theta_sum, same grid, $T", rd)
    end
    # AC6
    for T in (Float64, Double64)
        ra, rb, rc = 0.0, 0.0, 0.0
        for (im, mus) in enumerate(((0.45, 0.2), (0.45, 0.3, 0.1))), (inn, nn) in enumerate((false, true))
            a, b, c = y5_ratios(T, mus, nn, Val(2 * length(mus)), vseed(6, 10 * (T === Float64 ? 1 : 2) + 2 * im + inn))
            ra, rb, rc = max(ra, a), max(rb, b), max(rc, c)
            say("#   AC6 $T N=$(length(mus)) " * (nn ? "non-normal A_S" : "normal A") * ": (a) $(fmt(a)), (b) $(nn ? "n/a" : fmt(b)), (c) $(fmt(c))")
        end
        ratio_line("AC6", "(a) |C_i - ln vol_i (2048-bit Gram det)| / (i n theta_sum), N in {2,3}, A and A_S, n = 500, all variants, $T", ra)
        ratio_line("AC6", "(b) |C_i/n - Lambda_i - B_i/n| / (i theta_sum), A normal, $T", rb)
        ratio_line("AC6", "(c) |sum lambda_i| / theta_sum, $T", rc)
    end
end

# ---------------------------------------------------------------------------------------------------
# AC7: standard map against a 256-bit run on the same orbit points; n = 1e3 (Float64 and Double64, 50 starts, tau in {1,2,4}, all variants)
# and the stage 04 horizon n = N6 (Float64; (b),(c) on the full panel, (a) on the first N6_REF starts per K)
function ac7_stats(::Type{T}, K::Float64, starts, seed, n, taus, nref) where {T}
    eps_T = Float64(effective_eps(T))
    m = length(starts)
    sys = StandardMap(T(K))
    x0s = [hh_cast(T, z) for z in starts]
    W0s = [initial_deviations(seed, j, Val(2), Val(2), T) for j in 1:m]
    if n <= 10_000
        orb = [hh_orbit(sys, x0s[j], n) for j in 1:m]
        Mmax = [o[2] for o in orb]
        Psi = [hh_Psi(sys, orb[j][1], W0s[j]) for j in 1:m]
        zsref = j -> orb[j][1]
    else
        Mmax = [hh_Mmax(sys, x0s[j], n) for j in 1:m]
        Psi = [lb_Psi(sys, x0s[j], Float64.(W0s[j][:, 1]), n) for j in 1:m]
        zsref = j -> hh_orbit(sys, x0s[j], n)[1]
    end
    ref = Vector{Any}(undef, nref)
    tref = @elapsed Threads.@threads for j in 1:nref
        ref[j] = hh_ref_std(T(K), zsref(j), W0s[j], n)
    end
    ra, rb, rsum, cprime, e1 = 0.0, 0.0, 0.0, 0.0, 0.0
    for tau in taus
        lam, thetas = Dict(), Dict()
        for alg in hh_algs(T)
            res = propagate(sys, x0s, n; k = 2, tau = tau, seed = seed, W0 = W0s, reducers = (QRSpectrum(alg),))
            f = [res.finals[j][1] for j in 1:m]
            th = [hh_theta(T, f[j].G_max, tau, Mmax[j]) for j in 1:m]
            thp = [hh_theta_psi(T, f[j].G_max, tau, Mmax[j], Psi[j], n) for j in 1:m]
            setprecision(BigFloat, 256) do
                for j in 1:nref
                    errs = [Float64(abs(hh_big(f[j].lambda[i]) - ref[j][i])) for i in 1:2]
                    ra = max(ra, maximum(errs) / thp[j]); e1 = max(e1, errs[1] / (eps_T * Psi[j] / n))
                end
            end
            for j in 1:m
                s = Float64(abs(f[j].lambda[1] + f[j].lambda[2]))
                rsum = max(rsum, s / th[j]); cprime = max(cprime, s / (eps_T * (Float64(f[j].G_max) / tau + Mmax[j])))
            end
            lam[hh_name(alg)], thetas[hh_name(alg)] = [f[j].lambda for j in 1:m], thp
        end
        names = collect(keys(lam))
        for a in 1:length(names), b in (a + 1):length(names), j in 1:m
            rb = max(rb, maximum(Float64(abs(lam[names[a]][j][i] - lam[names[b]][j][i])) for i in 1:2) / (thetas[names[a]][j] + thetas[names[b]][j]))
        end
    end
    return (ra = ra, rb = rb, rsum = rsum, cprime = cprime, e1 = e1, tref = tref, m = m, nref = nref)
end

section("AC7") do
    cprime = 0.0
    nst = AC7_STARTS
    for (iK, K) in enumerate((0.5, 0.9716354, 2.0, 5.0, 10.0))
        rng = StableRNG(vseed(7, iK))
        starts = K == 0.5 ? hh_regular_starts(K, nst, rng)[1] : hh_starts(rng, nst)
        for (T, ns, tag) in ((Float64, nst, "Float64"), (Double64, nst, "Double64"))
            s = ac7_stats(T, K, starts[1:ns], vseed(7, 10 + iK), 1000, (1, 2, 4), ns)
            cprime = max(cprime, s.cprime)
            ratio_line("AC7", "(a) max_i |lambda_i - ref256| / theta_psi, K = $K, n = 1e3, tau in {1,2,4}, 4 variants, $ns starts, $tag (lambda_1 error / (eps Psi/n): $(fmt(s.e1)); spec observed 0.81)", s.ra)
            ratio_line("AC7", "(b) pairs of variants / (theta_psi^A + theta_psi^B), K = $K, n = 1e3, $tag", s.rb)
            ratio_line("AC7", "(c) |lambda_1 + lambda_2| / theta_sum, K = $K, n = 1e3, $tag (c' = $(fmt(s.cprime)))", s.rsum)
            T === Float64 && (RES[Symbol("ac7_", K)] = s)
        end
    end
    emit("AC7", "(c) largest measured c' of E11 = |lambda_1 + lambda_2| / (eps (G_max/tau + M_max)) over all K, T, tau, variants (n = 1e3)", cprime, "< 4 (diagnostics reference, P8)", pf(cprime < 4))
    # stage 04 horizon
    cp6 = 0.0
    for (iK, K) in enumerate((0.5, 0.9716354, 2.0, 5.0, 10.0))
        rng = StableRNG(vseed(7, 20 + iK))
        starts = K == 0.5 ? hh_regular_starts(K, nst, rng)[1] : hh_starts(rng, nst)
        s = ac7_stats(Float64, K, starts, vseed(7, 30 + iK), N6, (1, 2, 4), min(N6_REF, nst))
        cp6 = max(cp6, s.cprime)
        ratio_line("AC7", "(a) max_i |lambda_i - ref256| / theta_psi, K = $K, n = $N6, tau in {1,2,4}, 4 variants, first $(s.nref) of $nst starts (reference time $(fmt(s.tref)) s), Float64 (lambda_1 error / (eps Psi/n): $(fmt(s.e1)))", s.ra)
        ratio_line("AC7", "(b) pairs of variants / (theta_psi^A + theta_psi^B), K = $K, n = $N6, full $nst starts, Float64", s.rb)
        ratio_line("AC7", "(c) |lambda_1 + lambda_2| / theta_sum, K = $K, n = $N6, full $nst starts, Float64 (c' = $(fmt(s.cprime)))", s.rsum)
    end
    emit("AC7", "(c) largest c' at n = $N6, Float64", cp6, "< 4 (P8)", pf(cp6 < 4))
end

# ---------------------------------------------------------------------------------------------------
# AC8: cross-estimator agreement on the identical orbit
section("AC8") do
    seed = vseed(8)
    n = 10_000
    wa = Dict(hh_name(a) => 0.0 for a in hh_algs(Float64))
    for (iK, K) in enumerate((0.5, 2.0, 5.0, HH_TWO_PI))
        sys = StandardMap(K)
        starts = hh_starts(StableRNG(seed + iK), 50)
        Psi = [lb_Psi(sys, starts[j], initial_deviations(seed, j, Val(2), Val(1), Float64)[:, 1], n) for j in 1:50]
        for tau in (1, 4)
            X = propagate(sys, starts, n; k = 1, tau = tau, seed = seed, reducers = (Benettin(),))
            lam = [X.finals[j][1].lambda[1] for j in 1:50]
            for alg in hh_algs(Float64)
                Q = propagate(sys, starts, n; k = 1, tau = tau, seed = seed, reducers = (QRSpectrum(alg),))
                wa[hh_name(alg)] = max(wa[hh_name(alg)], maximum(abs(Q.finals[j][1].lambda[1] - lam[j]) / (16 * eps64 * Psi[j] / n) for j in 1:50))
            end
        end
    end
    ratio_line("AC8", "(a) |lambda_1^QR - X_Benettin| / (16 eps Psi_n/n), k = 1, K in {0.5,2,5,2pi}, 50 orbits, n = 1e4, tau in {1,4}; per variant " * join(["$k $(fmt(v))" for (k, v) in sort(collect(wa))], ", "), maximum(values(wa)))
    seed = vseed(8, 10)
    n = 1000
    wb = Dict("CGS" => 0.0, "MGS" => 0.0, "RGS" => 0.0); orthQR, orthGS = 0.0, 0.0
    for (iK, K) in enumerate((0.5, 2.0, 5.0, 10.0, HH_TWO_PI))
        sys = StandardMap(K)
        starts = hh_starts(StableRNG(seed + iK), 50)
        Mmax = [hh_Mmax(sys, z, n) for z in starts]
        Psi = [lb_Psi(sys, starts[j], initial_deviations(seed, j, Val(2), Val(2), Float64)[:, 1], n) for j in 1:50]
        for tau in (1, 2, 4)
            runq(alg) = propagate(sys, starts, n; k = 2, tau = tau, seed = seed, reducers = (QRSpectrum(alg),)).finals
            theta(f, j) = hh_theta_psi(Float64, f[j][1].G_max, tau, Mmax[j], Psi[j], n)
            gs = [(alg, runq(alg)) for alg in (CGS(), MGS(), RGS())]
            for (_, fg) in gs
                orthGS = max(orthGS, maximum(fg[j][1].orth_max for j in 1:50) / (2 * eps64))
            end
            for alg in hh_algs(Float64)
                fq = runq(alg)
                orthQR = max(orthQR, maximum(fq[j][1].orth_max for j in 1:50) / (2 * eps64))
                for (ag, fg) in gs
                    wb[hh_name(ag)] = max(wb[hh_name(ag)], maximum(maximum(abs(fq[j][1].lambda[i] - fg[j][1].lambda[i]) for i in 1:2) / (theta(fq, j) + theta(fg, j)) for j in 1:50))
                end
            end
        end
    end
    ratio_line("AC8", "(b) max_i |lambda_i^QR - lambda_i^GS| / (theta_psi^QR + theta_psi^GS), k = 2, K in {0.5,2,5,10,2pi}, 50 orbits, n = 1e3, tau in {1,2,4}, 4 QR variants; per GS variant " * join(["$k $(fmt(v))" for (k, v) in sort(collect(wb))], ", "), maximum(values(wb)))
    say("#   AC8(b) observed orthogonality loss max|Q'Q - I| / (k eps): QR variants $(fmt(orthQR)) (bound 8), Gram-Schmidt variants (CGS, MGS, RGS) $(fmt(orthGS))")
    seed = vseed(8, 20)
    ws, wsum = 0.0, 0.0
    for (iK, K) in enumerate((0.5, 2.0, 5.0, HH_TWO_PI))
        sys = StandardMap(K)
        starts = hh_starts(StableRNG(seed + iK), 50)
        Mmax = [hh_Mmax(sys, z, n) for z in starts]
        Psi = [lb_Psi(sys, starts[j], initial_deviations(seed, j, Val(2), Val(2), Float64)[:, 1], n) for j in 1:50]
        for tau in (1, 2, 4)
            fs = propagate(sys, starts, n; k = 2, tau = tau, seed = seed, reducers = (SVDSpectrum(),)).finals
            ths = [hh_theta_psi(Float64, fs[j][1].G_max, tau, Mmax[j], Psi[j], n) for j in 1:50]
            wsum = max(wsum, maximum(abs(fs[j][1].lambda[1] + fs[j][1].lambda[2]) / hh_theta(Float64, fs[j][1].G_max, tau, Mmax[j]) for j in 1:50))
            for alg in hh_algs(Float64)
                fq = propagate(sys, starts, n; k = 2, tau = tau, seed = seed, reducers = (QRSpectrum(alg),)).finals
                for j in 1:50
                    f = fs[j][1]
                    thq = hh_theta_psi(Float64, fq[j][1].G_max, tau, Mmax[j], Psi[j], n)
                    mid = f.lambda[1] - fq[j][1].lambda[1] - log(1 + f.c^2) / (2n)
                    lo, hi = -(thq + ths[j]), f.u / (2n) + thq + ths[j]
                    ws = max(ws, abs(mid - (lo + hi) / 2) / ((hi - lo) / 2))
                end
            end
        end
    end
    ratio_line("AC8", "(c) SVD V7 sandwich: |mid - centre| / half-width (<= 1 inside), K in {0.5,2,5,2pi}, 50 orbits, tau in {1,2,4}, n = 1e3, 4 variants", ws)
    ratio_line("AC8", "(c) |lambda_1^SVD + lambda_2^SVD| / theta_sum", wsum)
    # (d) ensemble and pathwise against the two-trajectory estimator
    seed = vseed(8, 30)
    sysW = WrappedStandardMap(HH_TWO_PI)
    x0w = [torus_state(z) for z in hh_starts(StableRNG(seed), 1000)]
    d0 = hh_d0star(HH_TWO_PI, HH_TWO_PI, eps64)
    q = propagate(sysW, x0w, 10_000; k = 2, tau = 1, seed = seed, reducers = (QRSpectrum(HouseholderQR()),))
    s = propagate_shadow(sysW, x0w, 10_000; d0 = d0, tau = 1, seed = seed)
    lq = [q.finals[j][1].lambda[1] for j in 1:1000]; x2 = [s.finals[j][1].X2 for j in 1:1000]
    se = hh_se(lq, x2) / 4
    ratio_line("AC8", "(d) ensemble K = 2pi wrapped, 1000 starts, n = 1e4, d0 = d0*: |mean lambda_1^QR - mean X2| / (4 combined SE); means $(fmt(hh_mean(lq))) and $(fmt(hh_mean(x2))), SE $(fmt(se))", abs(hh_mean(lq) - hh_mean(x2)) / hh_se(lq, x2))
    K = 0.5
    sysK = WrappedStandardMap(K)
    st, _ = hh_regular_starts(K, 50, StableRNG(seed + 1))
    x0k = [torus_state(z) for z in st]
    q = propagate(sysK, x0k, 1000; k = 2, tau = 1, seed = seed, reducers = (QRSpectrum(HouseholderQR()),))
    s = propagate_shadow(sysK, x0k, 1000; d0 = 1e-6, tau = 1, seed = seed)
    D = hh_D2T(K, 1e-6, HH_TWO_PI, eps64)
    Psik = [lb_Psi(sysK, x0k[j], initial_deviations(seed, j, Val(2), Val(2), Float64)[:, 1], 1000) for j in 1:50]
    ratio_line("AC8", "(d) pathwise K = 0.5, 50 starts, n = 1e3, d0 = 1e-6: |lambda_1^QR - X2| / (Delta_2T Psi_n/n + 3 B_S/n)",
               maximum(abs(q.finals[j][1].lambda[1] - s.finals[j][1].X2) / (D * Psik[j] / 1000 + 3 * hh_BS(Float64, 1000, hh_Gfirst(sysK, x0k[j], 1), K) / 1000) for j in 1:50))
end

# ---------------------------------------------------------------------------------------------------
# AC9 to AC12: invariance, zero allocation, health and scaling, state independence
function ac9_mismatch(alg, sys, x0, W0, seed, n, chunks)
    reducers = (QRSpectrum(alg),)
    serial = [propagate_orbit(sys, x0[j], W0[j], n, 1, reducers; orbit = j) for j in eachindex(x0)]
    sx, sW, sf = [s[1] for s in serial], [s[2] for s in serial], [map(reducer_final, reducers, s[3]) for s in serial]
    bad = 0
    for chunk in chunks
        res = propagate(sys, x0, n; k = 2, tau = 1, seed = seed, reducers = reducers, chunk = chunk, W0 = W0)
        bad += !isequal(res.x, sx) + !isequal(res.W, sW) + !isequal(res.finals, sf)
    end
    return bad
end

function health_stream_ok(alg, M, col)
    I2 = SMatrix{2, 2, Float64}(1, 0, 0, 1)
    sys = hh_Degen(12, M)
    x0 = [SVector(0.0, i == 3 ? 1.0 : 0.0) for i in 1:5]
    r = (QRSpectrum(alg),)
    ok = true
    for err in (hh_health_error(() -> propagate(sys, x0, 20; k = 2, tau = 5, seed = 1, W0 = fill(I2, 5), reducers = r, chunk = 2)),
                hh_health_error(() -> propagate_orbit(sys, x0[3], I2, 20, 5, r; orbit = 3)))
        ok &= err isa NumericalHealthError && err.what == "orthonormalisation breakdown" && err.orbit == 3 && err.step == 15 && err.column == col && err.backend === :cpu && err.value == 0.0
    end
    good = propagate(sys, x0[[1, 2, 4, 5]], 20; k = 2, tau = 5, seed = 1, W0 = fill(I2, 4), reducers = r)
    return ok & all(f -> all(iszero, f[1].sum) && f[1].n == 20, good.finals)
end

function health_direct_ok(alg, Wbad::SMatrix{d, k, Float64}, Wgood, col) where {d, k}
    r = QRSpectrum(alg); x = zero(SVector{d, Float64})
    st1, _ = reducer_step(r, reducer_init(r, Float64, Val(d), Val(k)), 1, x, Wgood)
    st2, W2 = reducer_step(r, st1, 2, x, Wbad)
    f1, f2 = reducer_final(r, st1), reducer_final(r, st2)
    e = reducer_error(r, st2, 7, :cpu)
    thrown = try; reducer_update(r, st1, 2, x, Wbad); false; catch ex; ex isa NumericalHealthError; end
    return !reducer_failed(r, st1) && reducer_failed(r, st2) && W2 === Wbad && f2.sum == f1.sum && f2.comp == f1.comp &&
           e.what == "orthonormalisation breakdown" && e.orbit == 7 && e.column == col && e.step == 2 && e.value == 0.0 && thrown
end

function scale_ratios(::Type{T}, vd::Val{d}, vk::Val{k}, s, alg, seed) where {T, d, k}
    eps_T = Float64(effective_eps(T))
    Qt = initial_deviations(seed, 1, vd, vk, T)
    W = T(s) * Qt
    Q, R, cnt, status, col = orthonormalize(alg, W)
    r = QRSpectrum(alg)
    st, _ = reducer_step(r, reducer_init(r, T, vd, vk), 1, zero(SVector{d, T}), W)
    fin = reducer_final(r, st)
    ok = status == QR_OK && all(isfinite, Q) && all(isfinite, R) && !reducer_failed(r, st) && isfinite(fin.G_max)
    return (ok, Float64(maximum(abs(R[i, i] / T(s) - 1) for i in 1:k)) / (8 * d * eps_T), Float64(maximum(abs, Q - Qt)) / (16 * d * eps_T),
            maximum(abs(Float64(fin.lambda[i]) - log(s)) / (8 * d * eps_T + 4 * eps_T * abs(log(s))) for i in 1:k))
end

section("AC9-12") do
    seed = vseed(9)
    m, n = 4096, 1000
    sys = StandardMap(HH_TWO_PI)
    x0 = hh_starts(StableRNG(seed), m)
    W0 = [initial_deviations(seed, j, Val(2), Val(2), Float64) for j in 1:m]
    for alg in hh_kernel_algs()
        count_line("AC9", "non-identical results (x, final Q, finals) of propagate vs serial propagate_orbit, 4096 orbits, n = 1e3, chunks 64 / 1000 / default, $(hh_name(alg)), $(Threads.nthreads()) threads", ac9_mismatch(alg, sys, x0, W0, seed, n, (64, 1000, nothing)))
    end
    c9 = ac9_mismatch(LapackQR(), sys, x0, W0, seed, n, (64, 1000, nothing))
    emit("AC9", "non-identical results, library path LapackQR (record; owner is asked if nonzero), chunks 64 / 1000 / default", c9, "0", pf(c9 == 0))
    # AC10
    bytes = Dict()
    for T in (Float64, Double64), alg in hh_kernel_algs(), (i, (vd, vk)) in enumerate(HH_AC10_DK)
        bytes[(T, hh_name(alg), i)] = hh_ac10(alg, T, vd, vk, seed + i)
    end
    for (i, (vd, vk)) in enumerate(HH_AC10_DK)
        bytes[(Float64, "LapackQR", i)] = hh_ac10(LapackQR(), Float64, vd, vk, seed + 10 + i)
    end
    for T in (Float64, Double64)
        sysT = StandardMap(T(HH_TWO_PI))
        for alg in hh_kernel_algs()
            bytes[(T, hh_name(alg), 0)] = hh_orbit_bytes(sysT, SVector{2, T}(T(1.0), T(2.0)), SMatrix{2, 2, T}(1, 0, 0, 1), 100, 1, (QRSpectrum(alg),))
        end
    end
    kern = [v for (k, v) in bytes if k[2] != "LapackQR"]; lib = [v for (k, v) in bytes if k[2] == "LapackQR"]
    emit("AC10", "max bytes allocated by one reducer_step (kernel variants, Float64 and Double64, (d,k) in {(2,1),(2,2),(4,4),(8,8)}) and by propagate_orbit n = 100, $(length(kern)) measurements", maximum(kern), "0", pf(maximum(kern) == 0))
    emit("AC10", "max bytes allocated by one reducer_step, library path (LapackQR, Float64, 4 shapes)", maximum(lib), "0", pf(maximum(lib) == 0))
    # AC11
    cases = (("zero column 2", SMatrix{2, 2, Float64}(1, 2, 0, 0), 2), ("identical columns", SMatrix{2, 2, Float64}(1, 0, 1, 0), 2), ("zero column 1", SMatrix{2, 2, Float64}(0, 0, 1, 2), 1))
    Wgood = SMatrix{2, 2, Float64}(2, 1, 0.5, 3)
    bad = 0; total = 0
    for alg in hh_algs(Float64), (name, M, col) in cases
        bad += !health_stream_ok(alg, M, col) + !health_direct_ok(alg, M, Wgood, col); total += 2
    end
    for alg in hh_algs(Float64)
        bad += !health_direct_ok(alg, SMatrix{4, 2, Float64}(1, 0, 0, 0, 1, 0, 0, 0), SMatrix{4, 2, Float64}(1, 2, 0, 1, 0, 1, 3, 0.5), 2); total += 1
    end
    count_line("AC11", "(a)(b) health cases failing (of $total: 4 variants x {zero column 2, identical columns, zero column 1} x {stream, direct} + d = 4 identical columns); kind observed \"orthonormalisation breakdown\" (spec text says \"deviation matrix\", superseded by D34), orbit 3, step 15, column, value 0, accumulators unchanged", bad)
    for T in (Float64, Double64)
        ok, ra, rq, rl = true, 0.0, 0.0, 0.0
        for (i, (vd, vk)) in enumerate(((Val(2), Val(2)), (Val(4), Val(4)), (Val(4), Val(2)))), s in (1e-200, 1e200), alg in hh_algs(T)
            o, a, b, c = scale_ratios(T, vd, vk, s, alg, seed + 20 + i)
            ok &= o; ra, rq, rl = max(ra, a), max(rq, b), max(rl, c)
        end
        ratio_line("AC11", "(c) W = s Q_test, s in {1e-200, 1e200}: max_i |R_ii/s - 1| / (8 d eps), (d,k) in {(2,2),(4,4),(4,2)}, all variants, $T; status OK, finite Q, R, G_max: $ok", ra)
        ratio_line("AC11", "(c) max|Q - Q_test| / (16 d eps), $T", rq)
        ratio_line("AC11", "(c) |lambda_i - ln s| / (8 d eps + 4 eps |ln s|) (no overflow or underflow of the norm), $T", rl)
        ok || emit("AC11", "(c) status or finiteness failure at s = 1e+-200, $T", false, "true", "FAIL")
    end
    # AC12
    bad12 = 0
    for T in (Float64, Double64), k in (1, 2), alg in hh_algs(T)
        K = T(HH_TWO_PI)
        z2 = [hh_cast(T, z) for z in hh_starts(StableRNG(vseed(12)), 20)]
        z4 = [SVector(z[1], z[2], zero(T), zero(T)) for z in z2]
        long = propagate(hh_Long(K), z4, 500; k = k, tau = 2, seed = vseed(12), reducers = (QRSpectrum(alg),))
        short = propagate(hh_Short(K), z2, 500; k = k, tau = 2, seed = vseed(12), reducers = (QRSpectrum(alg),))
        bad12 += !(isequal(long.W, short.W) && isequal(long.finals, short.finals))
    end
    count_line("AC12", "state length 4 vs 2: runs whose lambda, G or final Q differ bitwise (of 2 T x 2 k x variants), N = 500, tau = 2, 20 starts", bad12)
end

# ---------------------------------------------------------------------------------------------------
# AC13: tau sweep; R5 above tau_max
section("AC13") do
    seed = vseed(13)
    n = 10_000
    for T in (Float64, Double64)
        sys = StandardMap(T(HH_TWO_PI))
        x0 = hh_cast(T, hh_starts(StableRNG(seed), 1)[1])
        Mmax = hh_Mmax(sys, x0, n)
        rs, rd = 0.0, 0.0
        for alg in hh_algs(T)
            f = Dict(tau => propagate(sys, [x0], n; k = 2, tau = tau, seed = seed, reducers = (QRSpectrum(alg),)).finals[1][1] for tau in (1, 2, 4))
            th = Dict(tau => hh_theta(T, f[tau].G_max, tau, Mmax) for tau in (1, 2, 4))
            for tau in (1, 2, 4)
                rs = max(rs, Float64(abs(f[tau].lambda[1] + f[tau].lambda[2])) / th[tau])
                for i in 1:2
                    rd = max(rd, Float64(abs(f[tau].lambda[i] - f[1].lambda[i])) / (th[tau] + th[1]))
                end
            end
        end
        ratio_line("AC13", "|lambda_1 + lambda_2| / theta_sum(tau), K = 2pi, n = 1e4, tau in {1,2,4}, all variants, $T", rs)
        ratio_line("AC13", "|lambda_i(tau) - lambda_i(1)| / (theta_sum(tau) + theta_sum(1)), $T", rd)
    end
end

# ---------------------------------------------------------------------------------------------------
# AC15: backend resolution and provenance (exact strings), counted as failed checks
section("AC15") do
    seed = vseed(15)
    sys = StandardMap(2.0)
    x0 = hh_starts(StableRNG(seed), 8)
    run15(alg, req; extra = ()) = propagate(sys, x0, 20; k = 2, tau = 2, seed = seed, reducers = (QRSpectrum(alg), extra...), backend = req)
    none, present = hh_F.NoDevice(), hh_FakeDevice(true, true)
    checks = Bool[]
    for dev in (none, present), req in (:auto, :gpu)
        p = hh_with_device(() -> run15(LapackQR(), req), dev).provenance
        push!(checks, p.backend === :cpu && p.device == "cpu" && p.backend_fallback_reason == "reducer is CPU-only")
    end
    p = hh_with_device(() -> run15(HouseholderQR(), :gpu), hh_FakeDevice(false, false)).provenance
    push!(checks, p.backend === :cpu && p.backend_fallback_reason == "CUDA.jl not loadable")
    p = hh_with_device(() -> run15(HouseholderQR(), :gpu), hh_FakeDevice(true, false)).provenance
    push!(checks, p.backend === :cpu && p.backend_fallback_reason == "no functional device")
    for alg in (HouseholderQR(), AugmentedQR(), GivensQR())
        p = hh_with_device(() -> run15(alg, :auto), none).provenance
        push!(checks, p.backend === :cpu && p.backend_fallback_reason in ("CUDA.jl not loadable", "no functional device"))
    end
    setprecision(BigFloat, 128) do
        sysB = StandardMap(BigFloat(2))
        xB = [SVector{2, BigFloat}(BigFloat(0.5) + j, BigFloat(0.25)) for j in 1:3]
        for dev in (none, present), alg in (HouseholderQR(), GivensQR(), AugmentedQR())
            p = hh_with_device(() -> propagate(sysB, xB, 4; k = 2, tau = 2, seed = seed, reducers = (QRSpectrum(alg),), backend = :gpu), dev).provenance
            push!(checks, p.backend === :cpu && p.backend_fallback_reason == "kernel does not compile for the element type")
        end
        push!(checks, resolve_backend(:gpu, BigFloat; cpu_only_reducer = true).reason == "kernel does not compile for the element type")
    end
    push!(checks, try; reducer_init(QRSpectrum(LapackQR()), Double64, Val(2), Val(2)); false; catch e; e isa ArgumentError; end)
    for alg in hh_algs(Float64)
        push!(checks, run15(alg, :cpu).provenance.backend_fallback_reason == "")
        name, impl = HH_PROV[hh_name(alg)]
        expected = (estimator = "qr", orthonormalizer = name, iteration_path = "diagonal", implementation = impl, orthogonality_tol = 8.0)
        push!(checks, hh_nt_eq(reducer_provenance(QRSpectrum(alg)), expected))
        for res in (run15(alg, :cpu), run15(alg, :auto; extra = (ResidualMonitor(sys, 4.0, Inf),)))
            p = res.provenance
            push!(checks, all(k -> getproperty(p, k) == getproperty(expected, k), keys(expected)) && p.backend in (:cpu, :gpu) && p.k == 2 && p.tau == 2 && p.seed == seed &&
                          p.precision == "Float64" && all(k -> haskey(p, k), (:device, :backend_fallback_reason, :threads, :cpu, :chunk_size, :headroom, :gpu_chunks, :gpu_chunk_size)))
        end
    end
    push!(checks, reducer_provenance(QRSpectrum(HouseholderQR(); orthogonality_tol = 4)).orthogonality_tol === 4.0)
    count_line("AC15", "failed exact-string / provenance checks (of $(length(checks))): reason strings for library path, no runtime, no functional device, BigFloat; provenance keys and values of all 4 variants; ArgumentError for LapackQR Double64", count(!, checks))
end

# ---------------------------------------------------------------------------------------------------
# Algorithm item 2: the Gram-Schmidt AC7 case (A^20 in Float64): GS raises, Householder and Givens continue
section("alg2") do
    F(m) = Float64(hh_fibbig(m))
    A20 = SMatrix{2, 2, Float64}(F(41), F(40), F(40), F(39))          # A^20, A = [[2,1],[1,1]], exact in Float64
    sys = LinearSymplecticMap(A20)
    W0 = initial_deviations(vseed(2), 1, Val(2), Val(2), Float64)
    Mm = F(41)
    exact = setprecision(BigFloat, 256) do
        Wb = hh_bigmat(A20) * hh_bigmat(W0)
        log(sqrt(Wb[1, 1]^2 + Wb[2, 1]^2))
    end
    for alg in (CGS(), MGS(), RGS(), HouseholderQR(), AugmentedQR(), GivensQR(), LapackQR())
        r = QRSpectrum(alg)
        err = hh_health_error(() -> propagate_orbit(sys, zero(SVector{2, Float64}), W0, 1, 1, (r,); orbit = 1))
        if err isa NumericalHealthError
            emit("ALG2", "$(hh_name(alg)) at W = A^20 Q_0 (Float64, one handoff): raises, kind \"$(err.what)\", step $(err.step)", 1, "GS raises, HH/Givens continue", pf(alg isa Union{CGS, MGS, RGS}))
        else
            st = propagate_orbit(sys, zero(SVector{2, Float64}), W0, 1, 1, (r,))
            f = reducer_final(r, st[3][1])
            orth = Float64(maximum(abs, st[2]' * st[2] - I))
            e1 = setprecision(BigFloat, 256) do
                Float64(abs(hh_big(f.sum[1] + f.comp[1]) - exact))
            end
            emit("ALG2", "$(hh_name(alg)) at W = A^20 Q_0 (Float64): continues; max|Q'Q - I| / (8 k eps) = $(fmt(orth / (16 * eps64))), |S_1 - ln||A^20 q_1||| = $(fmt(e1 / eps64)) eps, |lambda_1 + lambda_2| = $(fmt(Float64(abs(f.lambda[1] + f.lambda[2])))) against theta_sum = $(fmt(hh_theta(Float64, f.G_max, 1, Mm))) (R_22 carries an O(1) relative error, the inherent conditioning of E10), G_max = $(fmt(Float64(f.G_max)))",
                 orth / (16 * eps64), "<= 1", pf(!(alg isa Union{CGS, MGS, RGS}) && orth <= 16 * eps64 && Float64(abs(f.lambda[1] + f.lambda[2])) <= hh_theta(Float64, f.G_max, 1, Mm)))
        end
    end
end

# ---------------------------------------------------------------------------------------------------
# Deviations from the impl notes, checked with numbers
function old_Gw(R::SMatrix{k, k, T}) where {k, T}       # window_Gw before the 1e100 rescale (copy of the committed body without its first two lines)
    Ri = zero(SMatrix{k, k, T})
    for j in 1:k
        Ri = setindex(Ri, inv(R[j, j]), j + (j - 1) * k)
        for i in (j - 1):-1:1
            s = zero(T)
            for l in (i + 1):j
                s += R[i, l] * Ri[l, j]
            end
            Ri = setindex(Ri, -s / R[i, i], i + (j - 1) * k)
        end
    end
    return sqrt(sum(abs2, R)) * sqrt(sum(abs2, Ri))
end

section("dev") do
    rng = StableRNG(vseed(20))
    diffs = 0; nsame = 0
    for _ in 1:20_000, k in (2, 4)
        sc = 10.0^(rand(rng) * 160 - 80)
        R = sc * (k == 2 ? SMatrix{2, 2, Float64}(1 + rand(rng), 0.0, randn(rng), 1 + rand(rng)) :
                  SMatrix{4, 4, Float64}(ntuple(l -> (i = mod(l - 1, 4) + 1; j = div(l - 1, 4) + 1; i > j ? 0.0 : (i == j ? 1 + rand(rng) : randn(rng))), 16)))
        a, b = Lyr.window_Gw(R), old_Gw(R)
        nsame += 1; diffs += !(a === b || isequal(a, b))
    end
    count_line("DEV", "G_w rescale (orthonormalize.jl): triangles with entries inside [1e-100, 1e100] where window_Gw differs bitwise from the pre-rescale code ($nsame draws, k in {2,4})", diffs)
    bad = 0; rel = 0.0
    for alg in (CGS(), MGS(), RGS(), HouseholderQR(), AugmentedQR(), GivensQR(), LapackQR()), s in (1e-200, 1e200)
        Qt = initial_deviations(vseed(20), 1, Val(2), Val(2), Float64)
        r = QRSpectrum(alg)
        W = s * Qt * SMatrix{2, 2, Float64}(1, 0, 0.5, 3)
        Ws = Qt * SMatrix{2, 2, Float64}(1, 0, 0.5, 3)
        f2 = reducer_final(r, reducer_step(r, reducer_init(r, Float64, Val(2), Val(2)), 1, zero(SVector{2, Float64}), W)[1])
        f3 = reducer_final(r, reducer_step(r, reducer_init(r, Float64, Val(2), Val(2)), 1, zero(SVector{2, Float64}), Ws)[1])
        bad += !isfinite(f2.G_max); rel = max(rel, abs(f2.G_max - f3.G_max) / f3.G_max)
    end
    count_line("DEV", "G_w rescale: variants x s in {1e-200, 1e200} with non-finite G_max at W = s Q B (7 cases x 2); max rel. difference of G_max to the unscaled W = Q B: $(fmt(rel))", bad)
    nt = Threads.nthreads()
    before = LinearAlgebra.BLAS.get_num_threads()
    reducer_init(QRSpectrum(LapackQR()), Float64, Val(2), Val(2))
    after = LinearAlgebra.BLAS.get_num_threads()
    emit("DEV", "LapackQR reducer_init sets BLAS threads globally: Julia threads $nt, BLAS threads before $before, after $after", after, "1 when Julia threads > 1", pf(nt == 1 || after == 1))
    W = hh_randW(StableRNG(vseed(21)), Float64, Val(4), Val(4))
    orthonormalize(LapackQR(), W)
    by = @allocated orthonormalize(LapackQR(), W)
    emit("DEV", "pure LapackQR orthonormalize (fresh workspace per call) allocates bytes (documented deviation; reducer_step is 0, see AC10)", by, "> 0 expected", "NOTE")
    Wb = setprecision(BigFloat, 128) do
        Wq = SMatrix{4, 2, BigFloat}(ntuple(i -> BigFloat(sin(i)), 8))
        Q, R, cnt, status, col = orthonormalize(HouseholderQR(), Wq)
        (status == QR_OK, Float64(maximum(abs, Q' * Q - I)))
    end
    emit("DEV", "generic BigFloat(128) HouseholderQR 4x2: status OK, max|Q'Q - I|", Wb[2], "status OK ($(Wb[1]))", pf(Wb[1] && Wb[2] < 1e-30))
    emit("DEV", "H1 (treppen path) not implemented; AC16 H1 versus H2 excluded by D30", "-", "-", "EXCLUDED (D30)")
end

# ---------------------------------------------------------------------------------------------------
# D35 reproduction: 6000 chaotic orbits, K in {2, 5}, n = 1e3, tau = 1, kernel Householder, 256-bit reference
section("D35") do
    n = 1000
    all_err, all_thp, all_sum, all_ths, all_K = Float64[], Float64[], Float64[], Float64[], Float64[]
    e1max, over = 0.0, 0
    for (iK, K) in enumerate((2.0, 5.0))
        starts = hh_starts(StableRNG(vseed(35, iK)), D35_PER_K)
        sys = StandardMap(K)
        W0s = [initial_deviations(vseed(35, 10 + iK), j, Val(2), Val(2), Float64) for j in 1:D35_PER_K]
        orb = [hh_orbit(sys, starts[j], n) for j in 1:D35_PER_K]
        Psi = [hh_Psi(sys, orb[j][1], W0s[j]) for j in 1:D35_PER_K]
        ref = Vector{Any}(undef, D35_PER_K)
        Threads.@threads for j in 1:D35_PER_K
            ref[j] = hh_ref_std(K, orb[j][1], W0s[j], n)
        end
        f = [r[1] for r in propagate(sys, starts, n; k = 2, tau = 1, seed = vseed(35, 10 + iK), W0 = W0s, reducers = (QRSpectrum(HouseholderQR()),)).finals]
        for j in 1:D35_PER_K
            thp = hh_theta_psi(Float64, f[j].G_max, 1, orb[j][2], Psi[j], n); th = hh_theta(Float64, f[j].G_max, 1, orb[j][2])
            errs = setprecision(BigFloat, 256) do
                [Float64(abs(hh_big(f[j].lambda[i]) - ref[j][i])) for i in 1:2]
            end
            push!(all_err, maximum(errs)); push!(all_thp, thp); push!(all_sum, Float64(abs(f[j].lambda[1] + f[j].lambda[2]))); push!(all_ths, th); push!(all_K, K)
            e1max = max(e1max, errs[1] / (eps64 * Psi[j] / n)); over += maximum(errs) > th
        end
    end
    RES[:d35] = (err = all_err, thp = all_thp, sum = all_sum, ths = all_ths, K = all_K)
    N = length(all_err)
    ratio_line("D35", "max_i |lambda_i - ref256| / theta_psi over $N chaotic orbits (spec: 0.076)", maximum(all_err ./ all_thp))
    ratio_line("D35", "|lambda_1 + lambda_2| / theta_sum over $N orbits (spec: 0.0035)", maximum(all_sum ./ all_ths))
    emit("D35", "max |lambda_1 - ref| / (eps Psi_n/n) over $N orbits (spec observed 0.81 on its 6000 orbits); the bound is B_S^Psi = 8 eps Psi_n/n, ratio to the bound $(fmt(e1max / 8))", e1max, "<= 8", pf(e1max <= 8))
    emit("D35", "orbits whose error exceeds theta_sum (these need the Psi term; spec found 5 of 6000)", over, "informational", "INFO")
end

# ---------------------------------------------------------------------------------------------------
# AC17: ChaosTools full-spectrum oracle (lifted rule, hand-coded Jacobian, identical start and Q_0 = I)
ct_rule(u, p, t) = (pp = u[2] + p[1] * sin(u[1]); SVector(u[1] + pp, pp))
ct_jac(u, p, t) = (c = p[1] * cos(u[1]); SMatrix{2, 2, Float64}(1 + c, c, 1.0, 1.0))
function ct_tands(K)
    ds = CT.DeterministicIteratedMap(ct_rule, SVector(0.1, 0.2), [K])
    return CT.TangentDynamicalSystem(ds; J = ct_jac, Q0 = SMatrix{2, 2, Float64}(1, 0, 0, 1))
end

section("AC17") do
    I2 = SMatrix{2, 2, Float64}(1, 0, 0, 1)
    K = 0.5
    z0, _ = hh_regular_starts(K, 1, StableRNG(vseed(17)))
    z0 = z0[1]
    n = 100_000
    lam_ct = CT.lyapunovspectrum(ct_tands(K), n; u0 = z0)
    for alg in hh_algs(Float64)
        f = propagate(StandardMap(K), [z0], n; k = 2, tau = 1, seed = vseed(17), W0 = [I2], reducers = (QRSpectrum(alg),)).finals[1][1]
        d = [abs(f.lambda[i] - lam_ct[i]) for i in 1:2]
        emit("AC17", "(a) regular K = 0.5, n = 1e5, $(hh_name(alg)): |lambda_i - ChaosTools| for i = 1, 2 (lambda = $(fmt(f.lambda[1])), $(fmt(f.lambda[2])); ChaosTools $(fmt(lam_ct[1])), $(fmt(lam_ct[2])); ln(n)/n = $(fmt(log(n) / n)))", maximum(d), "<= 1e-4", pf(maximum(d) <= 1e-4))
    end
    nst, nn = 1000, 10_000
    sysT = StandardMap(HH_TWO_PI)
    starts = hh_starts(StableRNG(vseed(17, 1)), nst)
    tands = ct_tands(HH_TWO_PI)
    lct = [CT.lyapunovspectrum(tands, nn; u0 = z)[1] for z in starts]
    for alg in (HouseholderQR(), AugmentedQR(), GivensQR(), LapackQR())
        res = propagate(sysT, starts, nn; k = 2, tau = 1, seed = vseed(17, 1), W0 = fill(I2, nst), reducers = (QRSpectrum(alg),))
        lq = [res.finals[j][1].lambda[1] for j in 1:nst]
        ratio_line("AC17", "(b) chaotic K = 2pi, 1000 starts, n = 1e4, $(hh_name(alg)): |mean lambda_1 - mean ChaosTools| / (4 combined SE); means $(fmt(hh_mean(lq))), $(fmt(hh_mean(lct))); max pathwise |difference| $(fmt(maximum(abs.(lq .- lct))))", abs(hh_mean(lq) - hh_mean(lct)) / hh_se(lq, lct))
    end
end

# ---------------------------------------------------------------------------------------------------
# Recorded lines R4, R5, R6 and the criteria outside this script
section("R") do
    I2 = SMatrix{2, 2, Float64}(1, 0, 0, 1)
    for K in (HH_TWO_PI, 10.0)
        starts = hh_starts(StableRNG(vseed(40, K > 7 ? 1 : 0)), 1000)
        res = propagate(StandardMap(K), starts, 10_000; k = 2, tau = 1, seed = vseed(40), W0 = fill(I2, 1000), reducers = (QRSpectrum(HouseholderQR()),))
        lq = [res.finals[j][1].lambda[1] for j in 1:1000]
        say("R4 | ensemble mean lambda_1, K = $(fmt(K)), 1000 starts, n = 1e4 (V-H) | $(fmt(hh_mean(lq))) +- $(fmt(sqrt(hh_var(lq) / 1000))) (SE) | benchmark rows: K = 2pi 1.172 (single orbit, T = 1e9), large K ln(K/2) = $(fmt(log(K / 2))) | RECORDED")
    end
    sys = StandardMap(HH_TWO_PI); x0 = hh_starts(StableRNG(vseed(41)), 1)[1]; Mm = hh_Mmax(sys, x0, 10_000)
    lam1 = propagate(sys, [x0], 10_000; k = 2, tau = 1, seed = vseed(41), reducers = (QRSpectrum(HouseholderQR()),)).finals[1][1]
    for tau in (8, 10, 16)
        f = nothing
        t = @elapsed (try
            f = propagate(sys, [x0], 10_000; k = 2, tau = tau, seed = vseed(41), reducers = (QRSpectrum(HouseholderQR()),)).finals[1][1]
        catch e
            say("R5 | tau = $tau, K = 2pi, n = 1e4 | raises: $(sprint(showerror, e)) (window product numerically rank deficient above tau_max) | no threshold | RECORDED")
        end)
        f === nothing && continue
        say("R5 | tau = $tau, K = 2pi, n = 1e4 | |lambda_1 + lambda_2| = $(fmt(Float64(abs(f.lambda[1] + f.lambda[2])))), theta_sum = $(fmt(hh_theta(Float64, f.G_max, tau, Mm))), G_max = $(fmt(Float64(f.G_max))), |lambda_1 - lambda_1(tau=1)| = $(fmt(Float64(abs(f.lambda[1] - lam1.lambda[1])))), time $(fmt(t)) s | no threshold | RECORDED")
    end
    for (d, vd) in ((2, Val(2)), (8, Val(8)))
        W = hh_randW(StableRNG(vseed(42)), Float64, vd, vd)
        ts = Dict()
        for alg in (HouseholderQR(), AugmentedQR(), GivensQR(), LapackQR(), CGS(), MGS(), RGS())
            orthonormalize(alg, W)
            ts[hh_name(alg)] = @belapsed orthonormalize($alg, $W)
        end
        base = ts["HouseholderQR"]
        say("R6 | orthonormalize time d = k = $d (loaded machine) | " * join(["$k $(fmt(1e6 * v)) us ($(fmt(v / base))x)" for (k, v) in sort(collect(ts); by = last)], ", ") * " | no threshold; H1 not built (D30) | RECORDED")
    end
    emit("AC14", "kernel path CPU vs GPU, backend resolution on the device, performance record for K15", "-", "-", "STAGE04-GPU (verify/gpu.jl, local GPU run)")
    emit("AC16", "H1 versus H2", "-", "-", "EXCLUDED (D30)")
end

# ---------------------------------------------------------------------------------------------------
# Figures
section("figure") do
    if haskey(RES, :ac1)
        fig = CM.Figure(size = (1500, 480))
        shapes = ["($(typeof(vd).parameters[1]),$(typeof(vk).parameters[1]))" for (vd, vk) in HH_DK]
        names = ["HouseholderQR", "AugmentedQR", "GivensQR", "LapackQR"]
        cols = [:crimson, :steelblue, :darkgreen, :darkorange]
        for (ip, (lab, idx)) in enumerate((("(a) column residual / (8 d eps ||w||)", 1), ("(b) max|Q'Q - I| / (8 k eps)", 2), ("(d) |R_ii - ref256| / (8 d eps kappa ||w||)", 3)))
            ax = CM.Axis(fig[1, ip]; xlabel = "(d,k)", ylabel = "worst ratio to bound over 200 W", title = lab, xticks = (1:7, shapes), yscale = log10)
            for (ia, nm) in enumerate(names), (T, mk) in ((Float64, :circle), (Double64, :utriangle))
                ys = [get(RES[:ac1], (T, idk, nm), (NaN, NaN, NaN))[idx] for idk in 1:7]
                any(!isnan, ys) && CM.scatter!(ax, (1:7) .+ 0.08 * (ia - 2.5), max.(ys, 1e-3); color = cols[ia], marker = mk, markersize = 11, label = "$nm $T")
            end
            CM.hlines!(ax, [1.0]; color = :black, linestyle = :dash)
            CM.ylims!(ax, 1e-2, 3)
            ip == 1 && CM.axislegend(ax; position = :lb, labelsize = 8, nbanks = 2)
        end
        CM.save(joinpath(OUT, "lyap-householder-qr-factor-error.png"), fig)
        say("# figure saved: lyap-householder-qr-factor-error.png")
    end
    if haskey(RES, :d35)
        r = RES[:d35]
        fig = CM.Figure(size = (1100, 480))
        ax = CM.Axis(fig[1, 1]; xscale = log10, yscale = log10, xlabel = "theta_psi(n) of the orbit", ylabel = "max_i |lambda_i - ref256|", title = "$(length(r.err)) chaotic orbits, n = 1e3, tau = 1, V-H kernel")
        for (K, c) in ((2.0, :crimson), (5.0, :steelblue))
            m = r.K .== K
            CM.scatter!(ax, r.thp[m], max.(r.err[m], 1e-20); color = (c, 0.35), markersize = 4, label = "K = $K")
        end
        lo, hi = extrema(r.thp)
        CM.lines!(ax, [lo, hi], [lo, hi]; color = :black, linestyle = :dash, label = "error = theta_psi")
        CM.axislegend(ax; position = :lt)
        ax2 = CM.Axis(fig[1, 2]; xscale = log10, yscale = log10, xlabel = "theta_sum of the orbit", ylabel = "|lambda_1 + lambda_2|", title = "sum identity against theta_sum")
        CM.scatter!(ax2, r.ths, max.(r.sum, 1e-20); color = (:darkgreen, 0.35), markersize = 4)
        lo, hi = extrema(r.ths)
        CM.lines!(ax2, [lo, hi], [lo, hi]; color = :black, linestyle = :dash)
        CM.save(joinpath(OUT, "lyap-householder-qr-theta-psi.png"), fig)
        say("# figure saved: lyap-householder-qr-theta-psi.png")
    end
end

say(@sprintf("# total wall time %.1f s", time() - T0))
close(LOG)
