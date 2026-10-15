# Stage 04 verification of module lyap-svd-stabilized (SVDSpectrum, d = k = 2; catalog L9).
# Run from the repo root:  julia +release -t 2 --project=verify verify/lyap-svd-stabilized.jl
# Spec: workspaces/module-build/stages/01-spec/output/lyap-svd-stabilized-spec.md (criteria 1 to 12, Amendment D35).
# One printed line per criterion:  AC<n> | quantity | value | bound | PASS/FAIL/INFO/STAGE04-GPU
# Values are "worst ratio of error to bound" (pass when <= 1) unless the quantity says otherwise.
# Seeds: StableRNG / run seed 20261012 + 50000 + 100 x criterion (package tests use 30100..31200, GS 20100..21400, standard-map verify 40000+).
# The measurement functions are re-written here (the package criteria are @testsets); the helpers (gs_*, sv_*) are included from the
# test folders (one home). Every 256-bit reference is built on the run's own element-type orbit (risk K18).
# SV_ONLY="AC3,AC11" runs only the named sections and writes the log to SV_LOG (debugging aid); the default runs everything.
using Printf, Random, LinearAlgebra, Statistics
using Chaoslab, StaticArrays, DoubleFloats, StableRNGs
import ChaosTools as CT
import DynamicalSystemsBase as DSB
import CairoMakie as CM

const ROOT = dirname(@__DIR__)
const OUT = joinpath(ROOT, "workspaces", "module-build", "stages", "04-verify", "output")
mkpath(OUT)
const ONLY = filter(!isempty, split(get(ENV, "SV_ONLY", ""), ","))
const LOG = open(get(ENV, "SV_LOG", joinpath(OUT, "lyap-svd-stabilized-verify-log.txt")), "w")

fmt(x::AbstractFloat) = @sprintf("%.4g", x)
fmt(x) = string(x)
function say(s)
    println(s); println(LOG, s); flush(stdout); flush(LOG)
end
emit(ac, quantity, value, bound, verdict) = say(string(ac, " | ", quantity, " | ", fmt(value), " | ", bound, " | ", verdict))
pf(ok::Bool) = ok ? "PASS" : "FAIL"
ratio_line(ac, quantity, r) = emit(ac, quantity, r, "<= 1", pf(r <= 1))

function section(f, name)
    isempty(ONLY) || name in ONLY || return
    t0 = time()
    try
        f()
    catch e
        say(string(name, " | ERROR | ", sprint(showerror, e)[1:min(end, 400)], " | - | FAIL"))
    end
    say(@sprintf("# section %s took %.1f s", name, time() - t0))
end

const TESTDIR = joinpath(ROOT, "test")
include(joinpath(TESTDIR, "lyap-gram-schmidt", "gs_helpers.jl"))
include(joinpath(TESTDIR, "lyap-svd-stabilized", "svd_helpers.jl"))

VS(ac) = GS_SEED + 50000 + 100 * ac       # this script's offsets: 50000 .. 51200
const Lyap = Chaoslab.Lyapunov
const T0 = time()
say("# lyap-svd-stabilized verification, julia $(VERSION), threads $(Threads.nthreads()), base seed $GS_SEED, offsets VS(ac) = base + 50000 + 100 ac, time $(round(time())) (unix)")
say("# ChaosTools $(pkgversion(CT)), DynamicalSystemsBase $(pkgversion(DSB)); effective eps: Float64 $(gs_eps(Float64)), Double64 $(gs_eps(Double64)) (= 2^-100)")
say("# timings in this log were taken on a loaded machine (other agents were running concurrently)")
say("# format: AC<n> | quantity | value | bound | verdict")

lname(T) = string(T)

# ---------------------------------------------------------------------------------------------------
# API V6: closed form against the 2 x 2 singular values at 256 bits (full size: 4 values of c x 20 for both types)
section("V6") do
    rng = StableRNG(VS(0))
    for T in gs_types()
        worst = 0.0
        for c in (0.0, 1.27, -0.19, 3.0e3), _ in 1:20
            S = SVector{2, T}(T(15 * rand(rng)), T(-15 * rand(rng)))
            out = log_singular_values(S, T(c))
            ref = setprecision(BigFloat, 256) do
                a, dl = exp(gs_big(S[1])), exp(gs_big(S[2]))
                b = BigFloat(c) * a
                F, Dt = a^2 + b^2 + dl^2, a * dl
                s1 = sqrt((F + sqrt(F^2 - 4 * Dt^2)) / 2)
                (log(s1), log(Dt / s1))
            end
            tol = 16 * gs_eps(T) * (1 + abs(Float64(S[1])) + abs(Float64(S[2])) + log1p(c^2))
            worst = max(worst, gs_absdiff(out[1], ref[1]) / tol, gs_absdiff(out[2], ref[2]) / tol)
        end
        ratio_line("V6", "max ratio |(ln sigma_1, ln sigma_2) - 256-bit| / 16 eps (1 + |S_1| + |S_2| + ln(1 + c^2)), c in {0, 1.27, -0.19, 3e3} x 20, $(lname(T))", worst)
    end
    out = log_singular_values(SVector(1.0e5, -1.0e5 + 3.0), 1.5)
    ok = all(isfinite, out) && out[1] >= out[2] && out[1] >= 1.0e5
    emit("V6", "S_1 = 1e5: finite, ordered; |sum - 3| / (16 eps 1e5)", abs(out[1] + out[2] - 3.0) / (16 * eps(Float64) * 1.0e5), "<= 1 and finite", pf(ok && abs(out[1] + out[2] - 3.0) <= 16 * eps(Float64) * 1.0e5))
    ok0 = all(isfinite, log_singular_values(SVector(0.0, 0.0), 0.0))
    emit("V6", "sigma_1 = sigma_2 (S = 0, c = 0) finite (V3)", ok0, "true", pf(ok0))
end

# ---------------------------------------------------------------------------------------------------
# AC1: symmetric map (cat map), tau in {1, 2, 4}, n in {40, 400}, 20 seeded frames
function ac1(T)
    seed = VS(1); lnmu = gs_lnmu(); sys = gs_cat(T)
    x0s = [gs_zero(T, 2) for _ in 1:20]
    w1 = w2 = wq = wp = 0.0; nord = 0
    for n in (40, 400)
        l1 = Dict{Tuple{Int, Int}, Tuple{T, Float64}}()
        for tau in (1, 2, 4)
            res = sv_run(sys, x0s, n; tau, seed)
            for j in 1:20
                f = sv_fin(res, j); th = gs_theta(T, f.G_max, tau, 2.0)
                q1 = initial_deviations(seed, j, Val(2), Val(2), T)[:, 1]
                w1 = max(w1, gs_absdiff(f.lambda[1], lnmu) / th)
                w2 = max(w2, gs_absdiff(f.lambda[2], -lnmu) / th)
                wq = max(wq, gs_absdiff(f.lambda_qr[1], gs_cat_lambda1(q1, n)) / th)
                nord += f.lambda[1] >= f.lambda[2] ? 0 : 1
                l1[(tau, j)] = (f.lambda[1], th)
            end
        end
        for a in keys(l1), b in keys(l1)
            a < b && (wp = max(wp, Float64(abs(l1[a][1] - l1[b][1])) / (l1[a][2] + l1[b][2])))
        end
    end
    ratio_line("AC1", "max ratio |lambda_1^SVD - ln mu| / theta_sum, cat map, tau 1,2,4, n 40,400, 20 frames, $(lname(T))", w1)
    ratio_line("AC1", "max ratio |lambda_2^SVD + ln mu| / theta_sum, $(lname(T))", w2)
    ratio_line("AC1", "max ratio |lambda_1^QR - G12a offset formula| / theta_sum, $(lname(T))", wq)
    ratio_line("AC1", "max ratio |lambda_1^SVD(frame/tau a) - lambda_1^SVD(b)| / (theta_a + theta_b), all pairs, $(lname(T))", wp)
    emit("AC1", "order violations lambda_1 < lambda_2, $(lname(T))", nord, "0", pf(nord == 0))
end
section("AC1") do
    for T in gs_types(); ac1(T); end
end

# ---------------------------------------------------------------------------------------------------
# AC2: K = 0
function ac2(T)
    seed = VS(2); sys = StandardMap(zero(T)); x0 = SVector{2, T}(T(0.3), T(0.7))
    frames = (SMatrix{2, 2, T}(0, 1, 1, 0), SMatrix{2, 2, T}(1, 0, 0, 1), initial_deviations(seed, 1, Val(2), Val(2), T))
    r1 = r2 = r3 = 0.0; mar = -Inf
    for (fi, W0) in enumerate(frames), n in (1, 10, 100, 1000, 10^4)
        f = sv_fin(sv_run(sys, [x0], n; tau = 1, seed, W0 = [W0]), 1)
        th = gs_theta(T, f.G_max, 1, 1.0)
        l1 = sv_V9b_lnsigma1(n)
        r1 = max(r1, gs_absdiff(f.lambda[1], l1 / n) / th)
        r2 = max(r2, Float64(abs(f.lambda[2] + f.lambda[1])) / th)
        if fi == 1
            lq = setprecision(BigFloat, 256) do
                log(BigFloat(n)^2 + 1) / 2
            end
            r3 = max(r3, gs_absdiff(f.lambda[1] - f.lambda_qr[1], (l1 - lq) / n) / th)
        end
        n == 10^4 && (mar = max(mar, Float64(f.lambda[1]) - log(1e4) / 1e4))
    end
    ratio_line("AC2", "max ratio |lambda_1^SVD - ln((n + sqrt(n^2 + 4))/2)/n| / theta_sum, 3 frames, n 1..1e4, $(lname(T))", r1)
    ratio_line("AC2", "max ratio |lambda_1 + lambda_2| / theta_sum, $(lname(T))", r2)
    ratio_line("AC2", "frame (e_2, e_1): max ratio |(SVD - QR) - [ln sigma_1 - (1/2) ln(n^2 + 1)]/n| / theta_sum, $(lname(T))", r3)
    emit("AC2", "lambda_1^SVD(1e4) - ln(1e4)/1e4 (derived margin 1.0e-12), max over frames, $(lname(T))", mar, "<= 2e-12", pf(mar <= 2e-12))
end
section("AC2") do
    for T in gs_types(); ac2(T); end
end

# ---------------------------------------------------------------------------------------------------
# AC3: identity with Gram-Schmidt, sandwich V7, c convergence by class (Amendment D35); Float64 and Double64 (Double64 is the stage 04 part)
const BINS = ("S1(500)<0", "0<=S1<2", "2<=S1<5", "5<=S1<25", "S1>=25")
binof(s) = s < 0 ? 1 : s < 2 ? 2 : s < 5 ? 3 : s < 25 ? 4 : 5
function ac3(T)
    seed = VS(3); n = 1000; epsT = gs_eps(T)
    nbit = nbad_a = 0; lo = hi = 0.0; nord = 0; nsand = 0
    cnt = zeros(Int, 5); relmax = zeros(5)
    hiwr = 0.0; nfail_hi = 0; fails = NTuple{4, Float64}[]; lowr = 0.0; nfail_lo = 0
    const_rows = String[]
    for K in gs_KS(), tau in (1, 2, 4)
        sys = gs_stdmap(T, K)
        x0s = gs_starts(T, seed + 1, 50)
        pts = [gs_points(sys, x0s[j], n) for j in 1:50]
        res = sv_run(sys, x0s, n; tau, seed)
        ref = gs_run(sys, x0s, n; k = 2, tau, seed, alg = RGS())
        res500 = sv_run(sys, x0s, 500; tau, seed)
        kv = Float64[]; cv = Float64[]
        for j in 1:50
            f, g = sv_fin(res, j), gs_fin(ref, j)
            nbit += 1
            nbad_a += (isequal(f.sum, g.sum) && isequal(f.comp, g.comp) && isequal(f.lambda_qr, g.lambda)) ? 0 : 1
            th = gs_theta(T, f.G_max, tau, gs_Mmax(sys, pts[j]))
            # the sandwich terms are evaluated at 256 bits (a Float64 log1p would put a 1e-19 error under a Double64 theta_sum)
            gapb, up = setprecision(BigFloat, 256) do
                cb = gs_big(f.c)
                (gs_big(f.lambda[1]) - gs_big(f.lambda_qr[1]) - log(1 + cb^2) / (2n), gs_big(f.u) / (2n))
            end
            gap = Float64(gapb)
            lo = max(lo, -gap / th)
            hi = max(hi, (gap - Float64(up)) / th)
            nord += (f.lambda[1] >= f.lambda[2] && f.lambda[1] + th >= f.lambda_qr[1]) ? 0 : 1
            push!(kv, log1p(Float64(f.c)^2) / (2n)); push!(cv, Float64(f.c))
            h = sv_fin(res500, j)
            S1_500 = Float64(h.sum[1] + h.comp[1])
            dc = Float64(abs(f.c - h.c))
            b = binof(S1_500); cnt[b] += 1
            relmax[b] = max(relmax[b], dc / max(1.0, Float64(abs(h.c))))
            if b == 5
                r = dc / (8 * epsT * max(1.0, Float64(abs(f.c))))
                hiwr = max(hiwr, r)
                if r > 1
                    nfail_hi += 1; push!(fails, (S1_500, dc, r, Float64(K isa Symbol ? 2pi : K)))
                end
            elseif b == 4
                r = dc / (0.06 * max(1.0, Float64(abs(h.c))))
                lowr = max(lowr, r); nfail_lo += r > 1 ? 1 : 0
            end
        end
        push!(const_rows, "K=$K tau=$tau: (1/(2n)) ln(1 + c_n^2) min/median/max = $(fmt(minimum(kv)))/$(fmt(median(kv)))/$(fmt(maximum(kv))); c median $(fmt(median(cv))), min $(fmt(minimum(cv))), max $(fmt(maximum(cv)))")
    end
    emit("AC3a", "QR-family output (sum, comp, lambda_QR) not bitwise equal to QRSpectrum(RGS()): count of $nbit orbits, $(lname(T))", nbad_a, "0", pf(nbad_a == 0))
    ratio_line("AC3b", "max ratio of the lower sandwich side -(lambda_1^SVD - lambda_1^QR - ln(1 + c^2)/(2n)) / theta_sum, $nbit orbits, $(lname(T))", lo)
    ratio_line("AC3b", "max ratio of the upper side (... - u/(2n)) / theta_sum, $(lname(T))", hi)
    emit("AC3b", "order/sandwich violations (lambda_1 >= lambda_2, lambda_1 + theta >= lambda_1^QR), $(lname(T))", nord, "0", pf(nord == 0))
    say("# AC3(c) $(lname(T)): classes by S_1(500) over $nbit orbits (K x tau x 50): " * join(["$(BINS[i]): $(cnt[i]) orbits, worst relative change $(fmt(relmax[i]))" for i in 1:5], "; "))
    emit("AC3c", "class S_1(500) >= 25: worst |c_1000 - c_500| / (8 eps_eff max(1, |c|)) over $(cnt[5]) orbits; failing orbits: $nfail_hi; $(lname(T))", hiwr, "<= 1", pf(hiwr <= 1))
    if !isempty(fails)
        s1s = [x[1] for x in fails]
        say("# AC3(c) $(lname(T)) failing S_1(500) >= 25 orbits: $(length(fails)); S_1(500) range [$(fmt(minimum(s1s))), $(fmt(maximum(s1s)))], worst |dc| $(fmt(maximum(x[2] for x in fails))), by K: " *
            join(["K=$(fmt(k)): $(count(x -> x[4] == k, fails))" for k in sort(unique(x[4] for x in fails))], ", "))
    end
    emit("AC3c", "class 5 <= S_1(500) < 25: worst |c_1000 - c_500| / (0.06 max(1, |c_500|)) over $(cnt[4]) orbits; failing: $nfail_lo; $(lname(T))", lowr, "<= 1", pf(lowr <= 1))
    emit("AC3c", "class S_1(500) < 5: $(cnt[1] + cnt[2] + cnt[3]) orbits reported, worst relative change $(fmt(maximum(relmax[1:3]))), $(lname(T))", maximum(relmax[1:3]), "reported, not scored", "INFO")
    for r in const_rows
        say("# AC3 measured constant $(lname(T)) $r")
    end
end
section("AC3") do
    for T in gs_types(); ac3(T); end
end

# ---------------------------------------------------------------------------------------------------
# AC4: 256-bit recomputation on the identical orbit, FULL spec size for both types (the package suite cut Double64 to 4 orbits, 2 K, 2 tau)
function ac4(T)
    seed = VS(4); ns = (12, 20, 100, 1000); m = 20
    worst = 0.0; werr = 0.0; per = Dict{Any, Float64}()
    for K in gs_KS()
        sys = gs_stdmap(T, K)
        x0s = gs_starts(T, seed + 1, m)
        pts = [gs_points(sys, x0s[j], maximum(ns)) for j in 1:m]
        Ys = [gs_ref_std_Y(sys.K, pts[j], ns) for j in 1:m]
        wk = 0.0
        for n in ns, tau in (1, 2, 4)
            res = sv_run(sys, x0s, n; tau, seed)
            for j in 1:m
                f = sv_fin(res, j)
                th = gs_theta(T, f.G_max, tau, gs_Mmax(sys, pts[j][1:n]))
                l1, l2 = sv_ref_logsigma(Ys[j][n])
                e = max(gs_absdiff(f.lambda[1], l1 / n), gs_absdiff(f.lambda[2], l2 / n))
                wk = max(wk, e / th); werr = max(werr, e / gs_eps(T))
            end
        end
        per[K] = wk; worst = max(worst, wk)
    end
    ratio_line("AC4", "max ratio |lambda_i^SVD - lambda_i^ref| / theta_sum, 256-bit product at the run's own orbit, K 4 values, 20 orbits, tau 1,2,4, n 12,20,100,1000, $(lname(T)) (per K: " * join(["$(k): $(fmt(per[k]))" for k in gs_KS()], ", ") * ")", worst)
    emit("AC4", "max absolute error in lambda / eps_eff (informational; scratch value 8 eps at n = 2000), $(lname(T))", werr, "-", "INFO")
end
section("AC4") do
    for T in gs_types(); ac4(T); end
end

# ---------------------------------------------------------------------------------------------------
# AC5: negative control (L10b) and overflow, K = 2 pi, tau = 1; plus the D35 60-start table
function naive_logsig2_dd(Y, l1ref)    # Double64 plain product: ln sigma_2 from the naive determinant (informational)
    d = abs(Y[1, 1] * Y[2, 2] - Y[1, 2] * Y[2, 1])
    return iszero(d) ? -Inf : Float64(log(d)) - l1ref
end
section("AC5") do
    seed = VS(5)
    for T in gs_types()
        sys = gs_stdmap(T, :twopi)
        x0s = gs_starts(T, seed + 1, 60)
        lnsig1(x0) = Float64(sv_ref_logsigma(gs_ref_std_Y(sys.K, gs_points(sys, x0, 20), (20,))[20])[1])
        ls = [lnsig1(x) for x in x0s]
        j = findfirst(>=(25), ls)
        j === nothing && (emit("AC5", "no start with ln sigma_1(20) >= 25 among 60, $(lname(T))", "-", "found", "FAIL"); continue)
        x0 = x0s[j]
        pts = gs_points(sys, x0, 2000)
        Ys = gs_ref_std_Y(sys.K, pts, (20, 2000))
        l2ref = sv_ref_logsigma(Ys[20])[2]
        if T === Float64
            gaps = Float64[]
            for x in x0s
                p20 = gs_points(sys, x, 20)
                Y20 = gs_ref_std_Y(sys.K, p20, (20,))[20]
                naive = svdvals(Matrix(sv_naive_product(sys, p20, 20)))
                push!(gaps, naive[2] == 0 ? Inf : abs(log(naive[2]) - Float64(sv_ref_logsigma(Y20)[2])))
            end
            big_ = gaps .> 10
            say("# AC5 D35 table, Float64, 60 starts at K = 2 pi, n = 20: ln sigma_1(20) from $(fmt(minimum(ls))) to $(fmt(maximum(ls))), median $(fmt(median(ls))); gap > 10 (or exact zero singular value) for $(count(big_)) of 60; " *
                "all orbits with ln sigma_1 >= 24.9 have gap > 10: $(all(big_[ls .>= 24.9])); smallest ln sigma_1 with gap > 10: $(fmt(minimum(ls[big_]))); largest ln sigma_1 with finite gap <= 10: $(fmt(maximum(ls[.!big_]))); starts with ln sigma_1 >= 25: $(count(>=(25), ls)); chosen start index $j (ln sigma_1(20) = $(fmt(ls[j])))")
            naive = svdvals(Matrix(sv_naive_product(sys, pts, 20)))
            gap = abs(log(naive[2]) - Float64(l2ref))
            emit("AC5", "plain Float64 product + svdvals at n = 20: |ln sigma_2^naive - ln sigma_2^ref| (naive $(fmt(log(naive[2]))), ref $(fmt(Float64(l2ref))); scratch 13.3)", gap, "> 10", pf(gap > 10))
        else
            Yn = sv_naive_product(sys, pts, 20)
            emit("AC5", "plain Double64 product at n = 20: |ln sigma_2^naive - ln sigma_2^ref| from the naive determinant (informational: Double64 keeps ~32 digits, so no gap is required; Float64 control is the L10b statement)", abs(naive_logsig2_dd(Yn, Float64(-l2ref)) - Float64(l2ref)), "-", "INFO")
        end
        nfirst = findfirst(n -> any(!isfinite, sv_naive_product(sys, pts, n)), 1:700)
        emit("AC5", "plain $(lname(T)) product has a non-finite entry first at n = $(nfirst === nothing ? "none up to 700" : nfirst) (derived 607, measured 618 in one orbit)", nfirst === nothing ? 700 : nfirst, "< 700", pf(nfirst !== nothing))
        wr = 0.0; allfin = true
        for n in (20, 2000)
            f = sv_fin(sv_run(sys, [x0], n; tau = 1, seed), 1)
            th = gs_theta(T, f.G_max, 1, gs_Mmax(sys, pts[1:n]))
            l1, l2 = sv_ref_logsigma(Ys[n])
            allfin &= all(isfinite, f.lambda) && isfinite(f.c)
            wr = max(wr, gs_absdiff(f.lambda[1], l1 / n) / th, gs_absdiff(f.lambda[2], l2 / n) / th)
        end
        emit("AC5", "estimator at n = 20 and 2000: finite = $allfin; max ratio |lambda_i - 256-bit| / theta_sum, $(lname(T))", wr, "<= 1", pf(allfin && wr <= 1))
    end
end

# ---------------------------------------------------------------------------------------------------
# AC6: rotation 0.7311 (exponents 0, 0)
section("AC6") do
    seed = VS(6)
    for T in gs_types(), tau in (1, 10)
        A = rotation_block(T, 0.7311); sys = LinearSymplecticMap(A)
        res = sv_run(sys, [gs_zero(T, 2) for _ in 1:4], 10^4; tau, seed)
        th = gs_theta_lin(T, A, tau)
        w = 0.0; fin = true
        for j in 1:4
            f = sv_fin(res, j)
            fin &= all(isfinite, f.lambda) && isfinite(f.c) && all(isfinite, res.W[j])
            w = max(w, maximum(i -> Float64(abs(f.lambda[i])) / th, 1:2))
        end
        emit("AC6", "rotation 0.7311, n = 1e4, tau = $tau, 4 frames, $(lname(T)): finite = $fin; max ratio |lambda_i| / theta_sum", w, "<= 1", pf(fin && w <= 1))
    end
end

# ---------------------------------------------------------------------------------------------------
# AC7: pairing and sum, both types (the package suite runs the standard map in Float64 only)
section("AC7") do
    seed = VS(7)
    function stats(f, T, th)
        s = f.lambda[1] + f.lambda[2]; S = f.sum .+ f.comp
        (Float64(abs(s)) / th, gs_absdiff(s, gs_big(S[1] + S[2]) / f.n) / (8 * gs_eps(T)),
         (eps_pair(f.lambda) == eps_sum(f.lambda) && f.lambda[1] >= f.lambda[2]) ? 0 : 1)
    end
    for T in gs_types()
        sys = gs_cat(T)
        res = sv_run(sys, [gs_zero(T, 2) for _ in 1:4], 10^5; tau = 1, seed)
        st = [stats(sv_fin(res, j), T, gs_theta(T, sv_fin(res, j).G_max, 1, 2.0)) for j in 1:4]
        ratio_line("AC7", "cat map n = 1e5, 4 frames, max |lambda_1 + lambda_2| / theta_sum, $(lname(T))", maximum(first, st))
        ratio_line("AC7", "cat map, max |lambda_1 + lambda_2 - (S_1 + S_2)/n| / 8 eps, $(lname(T))", maximum(x -> x[2], st))
        emit("AC7", "cat map, eps_pair != eps_sum or order violations, $(lname(T))", sum(x -> x[3], st), "0", pf(sum(x -> x[3], st) == 0))
        for K in gs_KS()
            sysK = gs_stdmap(T, K)
            x0s = gs_starts(T, seed + 1, 10)
            pts = [gs_points(sysK, x, 10^4) for x in x0s]
            resK = sv_run(sysK, x0s, 10^4; tau = 1, seed)
            stK = [stats(sv_fin(resK, j), T, gs_theta(T, sv_fin(resK, j).G_max, 1, gs_Mmax(sysK, pts[j]))) for j in 1:10]
            emit("AC7", "standard map K = $K n = 1e4, 10 orbits, $(lname(T)): max |sum|/theta_sum; max |sum - (S_1+S_2)/n|/8 eps; eps_pair/order violations $(sum(x -> x[3], stK))",
                 string(fmt(maximum(first, stK)), " ; ", fmt(maximum(x -> x[2], stK))), "<= 1 ; <= 1 ; 0",
                 pf(maximum(first, stK) <= 1 && maximum(x -> x[2], stK) <= 1 && sum(x -> x[3], stK) == 0))
        end
    end
end

# ---------------------------------------------------------------------------------------------------
# AC8: wrapped state (state length 4), both types
section("AC8") do
    seed = VS(8)
    for T in gs_types()
        nbad = 0; ntot = 0
        for K in (2.0, 5.0), tau in (1, 3)
            sys = WrappedStandardMap(T(K))
            rng = StableRNG(seed + 1)
            x0s = [SVector{4, T}(T(2 * Float64(pi) * rand(rng)), T(2 * Float64(pi) * rand(rng)), zero(T), zero(T)) for _ in 1:20]
            res = propagate(sys, x0s, 1200; k = 2, tau, seed, reducers = (SVDSpectrum(),), backend = :cpu)
            r = SVDSpectrum()
            for j in 1:20
                W0 = initial_deviations(seed, j, Val(2), Val(2), T)
                xm, Wm, st = gs_manual_loop(r, sys, x0s[j], W0, 1200, tau)
                ntot += 1
                nbad += (isequal(sv_fin(res, j), reducer_final(r, st)) && isequal(res.x[j], xm) && isequal(res.W[j], Wm) && length(res.x[j]) == 4) ? 0 : 1
            end
        end
        emit("AC8", "wrapped standard map K {2,5} x tau {1,3}, n = 1200, 20 orbits: final record, x, W not bitwise equal to the plain per-orbit loop: $nbad of $ntot, $(lname(T))", nbad, "0", pf(nbad == 0))
        sys = WrappedStandardMap(T(5))
        rng = StableRNG(seed + 2)
        x0s = [SVector{4, T}(T(2 * Float64(pi) * rand(rng)), T(2 * Float64(pi) * rand(rng)), zero(T), zero(T)) for _ in 1:8]
        t = @elapsed (res = propagate(sys, x0s, 10^6; k = 2, tau = 1, seed, reducers = (SVDSpectrum(),), backend = :cpu))
        w = 0.0; fin = true
        for j in 1:8
            f = sv_fin(res, j)
            fin &= all(isfinite, f.lambda) && isfinite(f.c)
            w = max(w, Float64(abs(f.lambda[1] + f.lambda[2])) / gs_theta(T, f.G_max, 1, 6.0))
        end
        emit("AC8", "10^6 steps, K = 5, 8 orbits ($(fmt(t)) s, loaded machine), $(lname(T)): finite = $fin, no health error; max |lambda_1 + lambda_2| / theta_sum (M_max = 6)", w, "<= 1", pf(fin && w <= 1))
        W0 = initial_deviations(seed, 1, Val(2), Val(2), T)
        al = gs_alloc_bytes(sys, x0s[1], W0, 10^4, 1, (SVDSpectrum(),))
        emit("AC8", "bytes allocated by propagate_orbit (10^4 steps, SVDSpectrum, wrapped map) after warm-up, $(lname(T))", al, "0", pf(al == 0))
    end
end

# ---------------------------------------------------------------------------------------------------
# AC9: thread and chunk invariance, 4096 orbits, K = 2 pi, n = 1000, tau = 1, both types
section("AC9") do
    seed = VS(9); n_orb, N = 4096, 1000
    for T in gs_types()
        sys = gs_stdmap(T, :twopi)
        x0s = gs_starts(T, seed + 1, n_orb)
        W0 = [initial_deviations(seed, j, Val(2), Val(2), T) for j in 1:n_orb]
        reducers = (SVDSpectrum(),)
        serial = [propagate_orbit(sys, x0s[j], W0[j], N, 1, reducers; orbit = j) for j in 1:n_orb]
        xs, Ws = [s[1] for s in serial], [s[2] for s in serial]
        fs = [map(reducer_final, reducers, s[3]) for s in serial]
        for chunk in (64, 1000, nothing)
            res = propagate(sys, x0s, N; k = 2, tau = 1, seed, reducers, chunk, W0, backend = :cpu)
            ok = isequal(res.x, xs) && isequal(res.W, Ws) && isequal(res.finals, fs)
            emit("AC9", "propagate ($(Threads.nthreads()) threads, chunk $(chunk === nothing ? "default" : chunk)) vs serial propagate_orbit: x, W, finals (S_i, compensation, c) bitwise, $(lname(T))", ok, "true (tolerance 0)", pf(ok))
        end
        W = initial_deviations(seed, 1, Val(2), Val(2), T)
        al = gs_alloc_bytes(sys, x0s[1], W, 10^4, 1, (SVDSpectrum(),))
        emit("AC9", "bytes allocated per orbit body after warm-up, $(lname(T))", al, "0", pf(al == 0))
    end
    say("# AC9 note: thread counts compared are serial (1 thread) against the $(Threads.nthreads())-thread kernel; a 4-thread run was not made")
end

# ---------------------------------------------------------------------------------------------------
# AC12: scope guard and provenance, both types
section("AC12") do
    seed = VS(12)
    for T in gs_types()
        nbad = 0
        for (sys, d, k) in ((sv_Boom4(), 4, 4), (sv_Boom4(), 4, 2), (sv_Boom4(), 4, 1), (gs_cat(T), 2, 1))
            x0s = [gs_zero(T, d)]
            err = try
                propagate(sys, x0s, 4; k, tau = 2, seed, reducers = (SVDSpectrum(),), backend = :cpu); nothing
            catch e
                e isa ArgumentError ? e : rethrow()
            end
            nbad += (err isa ArgumentError && err.msg == SV_GUARD_MESSAGE) ? 0 : 1
            err = try
                reducer_init(SVDSpectrum(), T, Val(d), Val(k)); nothing
            catch e
                e isa ArgumentError ? e : rethrow()
            end
            nbad += (err isa ArgumentError && err.msg == SV_GUARD_MESSAGE) ? 0 : 1
        end
        emit("AC12", "guard: exact message from propagate and reducer_init for (d,k) in {(4,4),(4,2),(4,1),(2,1)}; deviations from the exact string or not rejected before work (4-d map throws if stepped): $nbad of 8, $(lname(T))", nbad, "0", pf(nbad == 0))
        nprov = 0
        for (tol, tau) in ((8.0, 1), (4.0, 3), (Inf, 3))
            nprov += reducer_provenance(SVDSpectrum(; orthogonality_tol = tol)) == (estimator = "svd", orthonormalizer = "RGS", d = 2, orthogonality_tol = tol) ? 0 : 1
            r = SVDSpectrum(; orthogonality_tol = tol)
            res = propagate(gs_stdmap(T, 2.0), gs_starts(T, seed + 1, 3), 12; k = 2, tau, seed, reducers = (r,), backend = :cpu)
            p = res.provenance
            ok = p.estimator == "svd" && p.orthonormalizer == "RGS" && p.d == 2 && p.orthogonality_tol === tol &&
                 p.k == 2 && p.tau == tau && p.seed == seed && p.precision == precision_label(T) &&
                 p.backend === :cpu && p.device == "cpu" && p.backend_fallback_reason == "" && !haskey(p, :iteration_path) && !haskey(p, :implementation)
            nprov += ok ? 0 : 1
        end
        resm = propagate(gs_stdmap(T, 2.0), gs_starts(T, seed + 1, 3), 12; k = 2, tau = 1, seed, reducers = (SVDSpectrum(), ResidualMonitor(gs_stdmap(T, 2.0), 4.0, Inf)), backend = :cpu)
        nprov += (resm.provenance.estimator == "svd" && !haskey(resm.provenance, :iteration_path)) ? 0 : 1
        emit("AC12", "provenance record (estimator, orthonormalizer, d, orthogonality_tol, k, tau, seed, precision = $(precision_label(T)), backend, device, fallback reason; chained ResidualMonitor): mismatches $nprov of 4 cases, $(lname(T))", nprov, "0", pf(nprov == 0))
    end
end

# ---------------------------------------------------------------------------------------------------
# AC10: CPU vs GPU is run by verify/gpu.jl (D29); this script does not touch the GPU (shared with concurrent agents)
emit("AC10", "CPU vs GPU lambda_i^SVD (tangent-core AC12 tolerances) and |c_cpu - c_gpu| <= 1e-12 at n = 1e4, K = 10, 1024 orbits with S_1(500) >= 25", "not run here; verify/gpu.jl has no SVDSpectrum section", "AC12 tolerances; 1e-12", "STAGE04-GPU")

# ---------------------------------------------------------------------------------------------------
# AC11: ChaosTools full-spectrum oracle (lyapunovspectrum), lifted standard-map rule, hand-coded Jacobian, identical start and frame
# Facts read from the installed source (ChaosTools 3.5.4, DynamicalSystemsBase 3.19.3):
#   lyapunovspectrum(tands, N; u0) calls reinit!(tands, u0), whose default Q0 is default_deviations = the first k identity columns,
#   so a Q0 given at construction is overwritten: the oracle always starts from the identity frame. The Chaoslab side therefore
#   runs with W0 = identity (lyap-benettin S2 identical frame), not with the seeded Haar frame.
ct_rule(u, p, n) = (pn = u[2] + p[1] * sin(u[1]); SVector(u[1] + pn, pn))
ct_jac(u, p, n) = (c = p[1] * cos(u[1]); SMatrix{2, 2}(1 + c, c, 1, 1))
const I2 = SMatrix{2, 2, Float64}(1.0, 0.0, 0.0, 1.0)

# lambda_1, lambda_2 of ChaosTools for each start, threaded with one tangent system per task
function ct_spectra(K, x0s, N)
    out = Vector{SVector{2, Float64}}(undef, length(x0s))
    nt = Threads.nthreads()
    Threads.@threads :static for t in 1:nt
        ds = CT.DeterministicIteratedMap(ct_rule, x0s[1], [K])
        tands = DSB.TangentDynamicalSystem(ds; J = ct_jac, k = 2)
        for j in t:nt:length(x0s)
            out[j] = SVector{2, Float64}(CT.lyapunovspectrum(tands, N; u0 = x0s[j]))
        end
    end
    return out
end
chaoslab_svd(K, x0s, N) = propagate(StandardMap(K), x0s, N; k = 2, tau = 1, seed = VS(11), reducers = (SVDSpectrum(),), W0 = fill(I2, length(x0s)), backend = :cpu).finals

section("AC11") do
    seed = VS(11)
    # frame-reading check: the oracle result does not depend on Q0 given at construction (reinit! resets it)
    ds = CT.DeterministicIteratedMap(ct_rule, SVector(1.0, 2.0), [2.0])
    Qrot = SMatrix{2, 2}(0.0, 1.0, 1.0, 0.0)
    la = CT.lyapunovspectrum(DSB.TangentDynamicalSystem(ds; J = ct_jac, k = 2), 50; u0 = SVector(1.0, 2.0))
    lb = CT.lyapunovspectrum(DSB.TangentDynamicalSystem(ds; J = ct_jac, Q0 = Qrot), 50; u0 = SVector(1.0, 2.0))
    emit("AC11", "FINDING: lyapunovspectrum overwrites the constructor Q0 with the identity (reinit! default): result with Q0 = identity vs Q0 = swap are equal", la == lb, "true (so W0 = identity on the Chaoslab side)", "INFO")
    # (a) regular class K = 0.5, S7 screening: X(1e5) <= 5e-4 (lambda_1^QR of the estimator, identity frame, tau = 1), replaced starts counted
    sysR = StandardMap(0.5); N = 100_000
    rng = StableRNG(seed + 1)
    cands = SVector{2, Float64}[]; reg = Int[]
    while length(reg) < 20
        lo = length(cands) + 1
        append!(cands, [SVector(2 * Float64(pi) * rand(rng), 2 * Float64(pi) * rand(rng)) for _ in 1:40])
        fin = chaoslab_svd(0.5, cands[lo:end], N)
        append!(reg, (lo - 1) .+ findall(f -> f[1].lambda_qr[1] <= 5.0e-4, fin))
    end
    keep = reg[1:20]; x0r = cands[keep]
    say("# AC11 regular class: candidates screened $(length(cands)) (first 20 regular accepted at candidate index $(keep[end]); replaced = $(keep[end] - 20))")
    fin = chaoslab_svd(0.5, x0r, N)
    ct = ct_spectra(0.5, x0r, N)
    r_svd = [Float64(abs(fin[j][1].lambda[1] - ct[j][1])) / (1.0e-4 + log1p(Float64(fin[j][1].c)^2) / (2N)) for j in 1:20]
    pq = [Float64(abs(fin[j][1].lambda_qr[1] - ct[j][1])) for j in 1:20]
    emit("AC11", "K = 0.5 regular, first orbit (n = 1e5): lambda_1^SVD $(fmt(Float64(fin[1][1].lambda[1]))), oracle $(fmt(ct[1][1])), c_n $(fmt(Float64(fin[1][1].c))); |diff| / (1e-4 + ln(1 + c^2)/(2n))", r_svd[1], "<= 1", pf(r_svd[1] <= 1))
    ratio_line("AC11", "K = 0.5 regular, 20 orbits, n = 1e5: max |lambda_1^SVD - oracle| / (1e-4 + ln(1 + c_n^2)/(2n))", maximum(r_svd))
    emit("AC11", "K = 0.5 regular, 20 orbits: max |lambda_1^SVD - oracle| and (informational) max |lambda_1^QR - oracle| (pathwise, same orbit and frame)", string(fmt(maximum(Float64(abs(fin[j][1].lambda[1] - ct[j][1])) for j in 1:20)), " ; ", fmt(maximum(pq))), "-", "INFO")
    emit("AC11", "K = 0.5 regular, lambda_2: max |lambda_2^SVD - lambda_2^oracle| and |lambda_2^QR - lambda_2^oracle| (informational)", string(fmt(maximum(Float64(abs(fin[j][1].lambda[2] - ct[j][2])) for j in 1:20)), " ; ", fmt(maximum(Float64(abs(fin[j][1].lambda_qr[2] - ct[j][2])) for j in 1:20))), "-", "INFO")
    # (b) chaotic K = 2 pi, 1000 seeded orbits, n = 1e4: ensemble means (S5)
    K = 2 * Float64(pi); NC = 10^4
    x0c = gs_starts(Float64, seed + 2, 1000)
    t = @elapsed (finc = chaoslab_svd(K, x0c, NC))
    t2 = @elapsed (ctc = ct_spectra(K, x0c, NC))
    a = [Float64(f[1].lambda[1]) for f in finc]; b = [c[1] for c in ctc]; q = [Float64(f[1].lambda_qr[1]) for f in finc]
    se = sqrt(var(a) / 1000 + var(b) / 1000)
    emit("AC11", "K = 2 pi chaotic, 1000 orbits, n = 1e4 (Chaoslab $(fmt(t)) s, ChaosTools $(fmt(t2)) s): mean lambda_1^SVD $(fmt(mean(a))), oracle $(fmt(mean(b))), combined s.e. $(fmt(se)); |difference of means| / (4 s.e.)", abs(mean(a) - mean(b)) / (4 * se), "<= 1", pf(abs(mean(a) - mean(b)) <= 4 * se))
    emit("AC11", "K = 2 pi: pathwise max |lambda_1^SVD - oracle| and max |lambda_1^QR - oracle| over 1000 orbits (informational)", string(fmt(maximum(abs.(a .- b))), " ; ", fmt(maximum(abs.(q .- b)))), "-", "INFO")
    emit("AC11", "K = 2 pi: mean of (lambda_1^SVD - lambda_1^QR) over the ensemble, and max (the V7 offset; informational)", string(fmt(mean(a .- q)), " ; ", fmt(maximum(a .- q))), "-", "INFO")
end

# ---------------------------------------------------------------------------------------------------
# Deviations from the impl notes, checked
section("DEV") do
    r = SVDSpectrum()
    st0 = reducer_init(r, Float64, Val(2), Val(2))
    f0 = reducer_final(r, st0)
    names = (:n, :sum, :comp, :lambda_qr, :lambda, :c, :u, :G_max, :orth_max, :count_max)
    emit("DEV", "reducer_final field order equals the contract (n, sum, comp, lambda_qr, lambda, c, u, G_max, orth_max, count_max)", keys(f0) == names, "true", pf(keys(f0) == names))
    emit("DEV", "lambda is zero at n = 0 on a fresh state", f0.lambda == zero(SVector{2, Float64}), "true", pf(f0.lambda == zero(SVector{2, Float64})))
    # c unchanged when the QR step fails: a rank-deficient window (equal columns) and a zero window
    stg, Wg = reducer_step(r, st0, 1, SVector(0.0, 0.0), SMatrix{2, 2, Float64}(2.0, 1.0, 1.0, 1.0))   # one good step: c becomes nonzero
    stg, Wg = reducer_step(r, stg, 2, SVector(0.0, 0.0), SMatrix{2, 2, Float64}(2.0, 1.0, 1.0, 1.0) * Wg)
    Wbad = SMatrix{2, 2, Float64}(1.0, 2.0, 2.0, 4.0)
    st1, W1 = reducer_step(r, stg, 3, SVector(0.0, 0.0), Wbad)
    sw = reducer_step(r, st1, 4, SVector(0.0, 0.0), Wbad)
    emit("DEV", "failed QR step (rank-deficient window) after 2 good steps: code $(st1.qr.code), c unchanged ($(st1.c) == $(stg.c), nonzero: $(stg.c != 0)), W returned unchanged", st1.qr.code != 0 && stg.c != 0 && st1.c == stg.c && W1 == Wbad, "true", pf(st1.qr.code != 0 && stg.c != 0 && st1.c == stg.c && W1 == Wbad))
    emit("DEV", "step after a failure leaves state untouched (st.code != 0 shortcut)", isequal(sw[1], st1), "true", pf(isequal(sw[1], st1)))
    emit("DEV", "_qr_step still defined in qr_spectrum.jl (the shared step the SVD reducer calls)", isdefined(Lyap, :_qr_step), "true", pf(isdefined(Lyap, :_qr_step)))
    # reducer_failed / reducer_error route through the QR state
    emit("DEV", "reducer_failed(SVDSpectrum, failed state) is true and a fresh state is not failed", reducer_failed(r, st1) && !reducer_failed(r, st0), "true", pf(reducer_failed(r, st1) && !reducer_failed(r, st0)))
end

# ---------------------------------------------------------------------------------------------------
# Figure: (a) error of lambda_i^SVD against the 256-bit value in units of theta_sum, (b) the V7 sandwich with c settling
section("figure") do
    seed = VS(13)
    nmax = 10_000; nref = 2000
    refns = sort(unique(round.(Int, exp10.(range(log10(2), log10(nref); length = 40)))))
    fig = CM.Figure(size = (1200, 520))
    axa = CM.Axis(fig[1, 1]; xscale = log10, yscale = log10, xlabel = "n (map steps, tau = 1)", ylabel = "|lambda_i^SVD(n) - lambda_i^ref(n)| / theta_sum", title = "(a) singular-value exponents against 256 bits")
    axb = CM.Axis(fig[1, 2]; xscale = log10, yscale = log10, xlabel = "n", ylabel = "lambda_1^SVD - lambda_1^QR", title = "(b) V7: (1/(2n)) ln(1 + c_n^2) <= gap <= that + u/(2n)")
    say("# figure data (a): n up to $nref, worst ratio per curve below")
    for (T, K, col, lab) in ((Float64, :twopi, :crimson, "K = 2 pi Float64"), (Double64, :twopi, :darkorange, "K = 2 pi Double64"), (Float64, 0.5, :steelblue, "K = 0.5 Float64"))
        sys = gs_stdmap(T, K); x0 = gs_starts(T, seed + 1, 1)[1]
        pts = gs_points(sys, x0, nmax)
        W0 = initial_deviations(seed, 1, Val(2), Val(2), T)
        r = SVDSpectrum(); st = reducer_init(r, T, Val(2), Val(2)); W = W0; x = x0
        Mrun = 0.0; rows = Dict{Int, Any}()
        for n in 1:nmax
            W = tangent(sys, x) * W; Mrun = max(Mrun, Float64(maximum(abs, tangent(sys, x)))); x = forward(sys, x)
            st, W = reducer_step(r, st, n, x, W)
            (n in refns || n % 100 == 0 || n <= 20) && (rows[n] = (reducer_final(r, st), Mrun))
        end
        Ys = gs_ref_std_Y(sys.K, pts, Tuple(refns))
        e1 = Float64[]; e2 = Float64[]; nn = Int[]
        for n in refns
            f, M = rows[n]; th = gs_theta(T, f.G_max, 1, M)
            l1, l2 = sv_ref_logsigma(Ys[n])
            push!(nn, n); push!(e1, max(gs_absdiff(f.lambda[1], l1 / n) / th, 1e-4)); push!(e2, max(gs_absdiff(f.lambda[2], l2 / n) / th, 1e-4))
        end
        CM.lines!(axa, nn, e1; color = col, linewidth = 2, label = "$lab, i = 1")
        CM.lines!(axa, nn, e2; color = col, linestyle = :dash, label = "$lab, i = 2")
        say("# figure (a) $lab: max over n <= $nref of error / theta_sum: i=1 $(fmt(maximum(e1))), i=2 $(fmt(maximum(e2))) (values below 1e-4 are drawn at 1e-4)")
        T === Float64 && begin
            ns = sort(collect(keys(rows)))
            gap = [max(Float64(rows[n][1].lambda[1] - rows[n][1].lambda_qr[1]), 1e-18) for n in ns]
            lo = [max(log1p(Float64(rows[n][1].c)^2) / (2n), 1e-18) for n in ns]
            hi = [max(log1p(Float64(rows[n][1].c)^2) / (2n) + Float64(rows[n][1].u) / (2n), 1e-18) for n in ns]
            CM.lines!(axb, ns, gap; color = col, linewidth = 2, label = "$lab: measured gap")
            CM.lines!(axb, ns, lo; color = col, linestyle = :dash, label = "$lab: lower, ln(1 + c^2)/(2n)")
            CM.lines!(axb, ns, hi; color = col, linestyle = :dot, label = "$lab: upper, + u/(2n)")
            cs = [Float64(rows[n][1].c) for n in ns]
            say("# figure (b) $lab: c_n at n = 1000, 10000: $(fmt(Float64(rows[1000][1].c))), $(fmt(Float64(rows[10000][1].c))); min over n of gap - lower (negative = rounding below the bound) $(fmt(minimum(gap .- lo)))")
        end
    end
    CM.axislegend(axa; position = :rt, labelsize = 10)
    CM.axislegend(axb; position = :lb, labelsize = 9)
    CM.save(joinpath(OUT, "lyap-svd-stabilized-exponent-error.png"), fig)
    say("# figure saved: lyap-svd-stabilized-exponent-error.png")
end

# ---------------------------------------------------------------------------------------------------
# Allocation recheck at the end of the run (W7 of the standard-map report was context dependent)
section("alloc-recheck") do
    for T in gs_types()
        sys = WrappedStandardMap(T(5)); x0 = SVector{4, T}(T(1), T(2), zero(T), zero(T))
        W0 = initial_deviations(VS(8), 1, Val(2), Val(2), T)
        al = gs_alloc_bytes(sys, x0, W0, 10^4, 1, (SVDSpectrum(),))
        emit("AC8", "recheck at end of run: bytes allocated by propagate_orbit (wrapped, 10^4 steps) after warm-up, $(lname(T))", al, "0", pf(al == 0))
    end
end

# ---------------------------------------------------------------------------------------------------
say("# H1 to H3 (general-d algorithm, held, D30/D32): not scored in M2")
say("# AC3(c) class table above reproduces Amendment D35 with fresh seeds; the per-orbit bound ratio of D35 (largest 0.104) is derivation only, not scored, and is not recomputed here")
say(@sprintf("# total wall time %.1f s", time() - T0))
close(LOG)
