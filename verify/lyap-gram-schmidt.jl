# Stage 04 verification of module lyap-gram-schmidt (QRSpectrum with CGS, MGS, RGS).
# Run from the repo root:  julia +release -t 2 --project=verify verify/lyap-gram-schmidt.jl
# Spec: workspaces/module-build/stages/01-spec/output/lyap-gram-schmidt-spec.md (criteria 1 to 14, Amendment D35).
# One printed line per measurement:  AC<n> | quantity | value | bound | PASS/FAIL (INFO, DEFERRED, NOTE are not scored).
# Most "value" entries are the worst ratio error / bound (pass when <= 1), as in the package tests; the quantity text names
# the case that attains the worst value ("at ...").
# Seeds: StableRNG(20261012 + 50000 + 100 ac + sub); the package tests use 20261012 + 20000 + 100 ac, standard-map verify 40000+,
# lyap-two-trajectory verify 45000+, so the samples here are fresh.
# Environment: GS_ONLY="ac1,ac4" runs only the named sections (log goes to GS_LOG or the temp dir), GS_QUICK=1 shrinks the samples
# (smoke test only, never for the report).
using Printf, Random, LinearAlgebra, Statistics
using Chaoslab, StaticArrays, DoubleFloats, StableRNGs
import ChaosTools as CT
import CairoMakie as CM

const ROOT = dirname(@__DIR__)
const OUT = joinpath(ROOT, "workspaces", "module-build", "stages", "04-verify", "output")
mkpath(OUT)
const ONLY = filter(!isempty, split(get(ENV, "GS_ONLY", ""), ","))
const QUICK = !isempty(get(ENV, "GS_QUICK", ""))
const LOGPATH = isempty(ONLY) && !QUICK ? joinpath(OUT, "lyap-gram-schmidt-verify-log.txt") :
                get(ENV, "GS_LOG", joinpath(tempdir(), "gs-verify-partial-log.txt"))
const LOG = open(LOGPATH, "w")

fmt(x::AbstractFloat) = @sprintf("%.4g", x)
fmt(x) = string(x)
function say(s)
    println(s); println(LOG, s); flush(stdout); flush(LOG)
end
emit(ac, quantity, value, bound, verdict) = say(string(ac, " | ", quantity, " | ", fmt(value), " | ", bound, " | ", verdict))
pf(ok::Bool) = ok ? "PASS" : "FAIL"

# Reuse the helpers of the package tests (one home): gs_helpers.jl whole, gs_s_loop sliced from gs_ac06_09.jl.
const TESTDIR = joinpath(ROOT, "test", "lyap-gram-schmidt")
include(joinpath(TESTDIR, "gs_helpers.jl"))
let txt = replace(read(joinpath(TESTDIR, "gs_ac06_09.jl"), String), "\r\n" => "\n")
    i = first(findfirst("function gs_s_loop", txt))
    j = last(findnext("\nend\n", txt, i))
    include_string(Main, txt[i:j], "gs_s_loop")
end

const VS = GS_SEED + 50000
vs(ac, sub = 0) = VS + 100 * ac + sub
const Lyap = Chaoslab.Lyapunov
const TWOPI = 2 * Float64(pi)

# Worst-value accumulator: keeps the largest ratio and the label of the case that produced it.
mutable struct Worst
    v::Float64
    at::String
end
Worst() = Worst(-Inf, "")
function upd!(w::Worst, v, at)
    x = Float64(v)
    isnan(x) && (x = Inf)
    x > w.v && (w.v = x; w.at = string(at))
    return w
end
emitw(ac, quantity, w::Worst; bound = 1.0, verdict = nothing) =
    emit(ac, string(quantity, " [worst at ", w.at, "]"), w.v, string("<= ", fmt(bound)), verdict === nothing ? pf(w.v <= bound) : verdict)

# Accumulators shared by sections: orthogonality ratios (AC6 a, b), the O3 constant, RGS pass counts, CGS trips.
const ORTH = Dict{Tuple{String,String},Worst}()      # (type, class) => worst orth_max / bound
const O3 = Dict{String,Worst}()                      # type => worst |sum lambda| / (eps (G_max/tau + M_max))
const PASSES = Dict{Tuple{String,Int},Int}()         # (type, count_max) => number of RGS final records
const TRIPS = Tuple{String,String,Int,String}[]      # (where, type, tau, alg) of CGS runs that raised the orthogonality error
tname(::Type{T}) where {T} = string(T)

function track!(T, alg, f, k, at)
    cls = alg isa RGS ? "RGS (8 k eps)" : "CGS/MGS (8 k G_w eps)"
    if alg isa RGS || k == 2     # the spec scopes (b) to k = 2 for CGS and MGS
        upd!(get!(ORTH, (tname(T), cls), Worst()), Float64(f.orth_max) / gs_orth_bound(alg, T, k, f.G_max), at)
    end
    alg isa RGS && (PASSES[(tname(T), f.count_max)] = get(PASSES, (tname(T), f.count_max), 0) + 1)
    return nothing
end
# O3 constant c = |sum lambda_i| / (eps (G_max/tau + M_max)) = 16 |sum| / theta_sum (pairing_residual_tol = 16)
track_o3!(T, sum_abs, theta, at) = upd!(get!(O3, tname(T), Worst()), 16 * sum_abs / theta, at)

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

sz(full, quick) = QUICK ? quick : full
const T0 = time()
say("# lyap-gram-schmidt verification, julia $(VERSION), threads $(Threads.nthreads()), base seed $GS_SEED, seed block +50000 (packages tests +20000), time $(round(time())) (unix)$(QUICK ? " QUICK (smoke test)" : "")")
say("# timings in this log were taken on a loaded machine (other agents were running concurrently)")
say("# effective eps: Float64 $(fmt(gs_eps(Float64))), Double64 $(fmt(gs_eps(Double64))) (= 2^$(round(Int, log2(gs_eps(Double64)))))")

# ---------------------------------------------------------------------------------------------------
# AC1 cat map (G12a)
section("ac1") do
    seed = vs(1)
    for T in gs_types()
        w1, w2, w3, w4 = Worst(), Worst(), Worst(), Worst()
        for alg in gs_algs(), tau in (1, 2, 4), n in (40, 400)
            at = "$(gs_name(alg)) tau=$tau n=$n"
            sys = gs_cat(T)
            x0s = [gs_zero(T, 2) for _ in 1:20]
            res = gs_run(sys, x0s, n; k = 2, tau, seed, alg)
            for j in 1:20
                f = gs_fin(res, j)
                th = gs_theta(T, f.G_max, tau, 2.0)
                q1 = initial_deviations(seed, j, Val(2), Val(2), T)[:, 1]
                upd!(w1, gs_absdiff(f.lambda[1], gs_cat_lambda1(q1, n)) / th, at)
                upd!(w2, Float64(abs(f.lambda[1] + f.lambda[2])) / th, at)
                track!(T, alg, f, 2, "AC1 " * at)
            end
            f = gs_fin(gs_run(sys, x0s[1:1], n; k = 2, tau, seed, alg, W0 = [gs_cat_eigframe(T)]), 1)
            upd!(w3, gs_absdiff(f.lambda[1], gs_lnmu()) / gs_theta(T, f.G_max, tau, 2.0), at)
            upd!(w4, Float64(f.orth_max) / gs_orth_bound(alg, T, 2, f.G_max), at)
        end
        emitw("AC1", "|lambda_1 - G12a| / theta_sum, 20 frames, $T, all variants, tau in {1,2,4}, n in {40,400}", w1)
        emitw("AC1", "|lambda_1 + lambda_2| / theta_sum, $T", w2)
        emitw("AC1", "eigenbasis frame |lambda_1 - ln mu| / theta_sum, $T", w3)
        emitw("AC1", "orth_max / orthogonality bound (eigenbasis frame runs), $T", w4)
    end
end

# ---------------------------------------------------------------------------------------------------
# AC2 K = 0 (G12b)
section("ac2") do
    seed = vs(2)
    rng = StableRNG(seed + 1)
    for T in gs_types()
        wS1, wS2, wlam = Worst(), Worst(), Worst()
        for alg in gs_algs()
            sys = StandardMap(zero(T))
            x0 = SVector{2,T}(T(TWOPI * rand(rng)), T(TWOPI * rand(rng)))
            W0 = SMatrix{2,2,T}(0, 1, 1, 0)
            for n in (1, 10, 1000, 10^4)
                f = gs_fin(gs_run(sys, [x0], n; k = 2, tau = 1, seed, alg, W0 = [W0]), 1)
                S = gs_S(f)
                ref = setprecision(BigFloat, 256) do
                    log(BigFloat(n)^2 + 1) / 2
                end
                at = "$(gs_name(alg)) n=$n"
                upd!(wS1, gs_absdiff(S[1], ref) / (8 * n * gs_eps(T)), at)
                upd!(wS2, Float64(abs(S[2] + S[1])) / (8 * n * gs_eps(T)), at)
                n == 10^4 && upd!(wlam, Float64(f.lambda[1]) / (log(1e4) / 1e4 + 1e-12), at)
            end
        end
        emitw("AC2", "|S_1(n) - (1/2) ln(n^2+1)| / (8 n eps), n in {1,10,1000,1e4}, $T", wS1)
        emitw("AC2", "|S_2 + S_1| / (8 n eps), $T", wS2)
        emitw("AC2", "lambda_1(1e4) / (ln(1e4)/1e4 + 1e-12), $T", wlam)
    end
end

# ---------------------------------------------------------------------------------------------------
# AC3 bitwise structure
section("ac3") do
    seed = vs(3)
    for T in gs_types()
        mism = 0; cases = 0; total = 0; dS = 0.0
        for K in (0.5, 2.0, 5.0, 10.0), tau in (1, 3)
            sys = gs_stdmap(T, K)
            m = T === Float64 ? sz(100, 10) : sz(20, 5)
            x0s = gs_starts(T, seed + 1, m)
            rc = gs_run(sys, x0s, 1200; k = 2, tau, seed, alg = CGS())
            rm = gs_run(sys, x0s, 1200; k = 2, tau, seed, alg = MGS())
            mism += count(j -> !(isequal(gs_fin(rc, j).sum, gs_fin(rm, j).sum) && isequal(gs_fin(rc, j).comp, gs_fin(rm, j).comp)), 1:m)
            mism += isequal(rc.W, rm.W) ? 0 : 1
            total += m + 1; cases += 1
            dS = max(dS, maximum(j -> Float64(maximum(abs, gs_S(gs_fin(rc, j)) - gs_S(gs_fin(rm, j)))), 1:m))
        end
        emit("AC3", "(a) standard map K in {0.5,2,5,10}, tau in {1,3}, n=1200, k=2: orbit comparisons (of $total in $cases cases) where CGS and MGS differ in S_i, comp or final Q, $T", mism, "0", pf(mism == 0))
        emit("AC3", "(a) magnitude: max over orbits of |S_i(CGS) - S_i(MGS)| / eps, $T", dS / gs_eps(T), "(information)", "INFO")
    end
    rng = StableRNG(seed + 2)
    for (T, d) in ((Float64, 2), (Float64, 4), (Double64, 2))
        bad = 0; N = sz(600, 50)
        for _ in 1:N
            W = SMatrix{d,2,T}(ntuple(_ -> T(randn(rng)), Val(2d)))
            fc, fm = orthonormalize(CGS(), W), orthonormalize(MGS(), W)
            bad += (isequal(fc[1], fm[1]) && isequal(fc[2], fm[2])) ? 0 : 1
        end
        emit("AC3", "(a) orthonormalize Q and R, CGS vs MGS, $N random $d x 2 matrices, $T: mismatches", bad, "0", pf(bad == 0))
    end
    # Diagnosis of a (a) mismatch: R12 = <q1, w2> is dot(q1, w2) in both variants; StaticArrays' dot is an @inbounds @simd loop, whose
    # reassociation/contraction is compiler's choice. Count which rounding each variant's R12 equals (plain = two roundings, fma = fused).
    rng2 = StableRNG(seed + 3)
    cnt = Dict(nm => [0, 0, 0, 0] for nm in ("CGS", "MGS"))   # plain, fma (q1 first), fma (q2 first), none of the three
    NW = sz(2000, 200)
    for _ in 1:NW
        W = SMatrix{2,2,Float64}(ntuple(_ -> randn(rng2), Val(4)))
        for alg in (CGS(), MGS())
            Q, R = orthonormalize(alg, W)
            q, w = Q[:, 1], W[:, 2]
            plain = q[1] * w[1] + q[2] * w[2]
            f1 = fma(q[1], w[1], q[2] * w[2]); f2 = fma(q[2], w[2], q[1] * w[1])
            c = cnt[gs_name(alg)]
            c[1] += R[1, 2] == plain; c[2] += R[1, 2] == f1; c[3] += R[1, 2] == f2
            c[4] += !(R[1, 2] == plain || R[1, 2] == f1 || R[1, 2] == f2)
        end
    end
    for nm in ("CGS", "MGS")
        say("AC3 | diagnosis: of $NW random 2x2 Float64 matrices, R12 of $nm equals plain q1*w1+q2*w2 / fma(q1,w1,q2*w2) / fma(q2,w2,q1*w1) / none: $(join(cnt[nm], " / ")) (equal counts mean the two roundings coincide on that matrix) | recorded | - | INFO")
    end
    for T in gs_types()
        A3 = gs_AS(T)^3
        Qc = Qm = initial_deviations(seed + 4, 1, Val(4), Val(4), T)
        differ = 0
        for _ in 1:200
            Qc = orthonormalize(CGS(), A3 * Qc)[1]
            Qm = orthonormalize(MGS(), A3 * Qm)[1]
            differ += Qc != Qm
        end
        emit("AC3", "(b) A_S k=4 tau=3: windows (of 200) where the CGS and MGS chains differ in Q, $T", differ, ">= 1", pf(differ >= 1))
    end
    sys = LinearSymplecticMap(gs_AS(Float64))
    rc = gs_run(sys, [gs_zero(Float64, 4)], 600; k = 4, tau = 3, seed, alg = CGS(), tol = Inf)
    rm = gs_run(sys, [gs_zero(Float64, 4)], 600; k = 4, tau = 3, seed, alg = MGS(), tol = Inf)
    emit("AC3", "(b) final Q of CGS and MGS streams differ (A_S, tau=3, Float64)", rc.W != rm.W, "true", pf(rc.W != rm.W))
end

# ---------------------------------------------------------------------------------------------------
# AC4 256-bit recomputation on the identical orbit
section("ac4") do
    seed = vs(4)
    n = 1000
    for (T, Ks, m) in ((Float64, gs_KS(), sz(20, 3)), (Double64, (2.0, 5.0, :twopi), sz(8, 2)))
        w, wraw, wo = Worst(), Worst(), Worst()
        for K in Ks
            sys = gs_stdmap(T, K)
            x0s = gs_starts(T, seed + 1, m)
            pts = [gs_points(sys, x0s[j], n) for j in 1:m]
            Mmax = [gs_Mmax(sys, pts[j]) for j in 1:m]
            ref = [gs_ref_std_S(sys.K, pts[j], initial_deviations(seed, j, Val(2), Val(2), T), n) for j in 1:m]
            for alg in gs_algs(), tau in (1, 2, 4)
                res = gs_run(sys, x0s, n; k = 2, tau, seed, alg)
                for j in 1:m
                    f = gs_fin(res, j)
                    th = gs_theta(T, f.G_max, tau, Mmax[j])
                    at = "K=$K $(gs_name(alg)) tau=$tau"
                    for i in 1:2
                        e = gs_absdiff(f.lambda[i], ref[j][i] / n)
                        upd!(w, e / th, at)
                        upd!(wraw, e / (th / 16), at)
                    end
                    track!(T, alg, f, 2, "AC4 " * at)
                end
            end
        end
        emitw("AC4", "|lambda_i - lambda_i^256|/theta_sum, i=1,2, K panel, n=1000, all variants, tau in {1,2,4}, $(length(Ks)) K x $m orbits, $T (eps_eff = $(fmt(gs_eps(T))))", w)
        emit("AC4", "scratch constant: max |lambda_i - ref| / (eps (G_max/tau + M_max)) [worst at $(wraw.at)], $T", wraw.v, "<= 16 (spec scratch 0.9)", pf(wraw.v <= 16))
    end
end

# ---------------------------------------------------------------------------------------------------
# AC5 cross variant and cross tau
function ac5_std(T, seed, m)
    n = 1000
    w, wsum = Worst(), Worst()
    for K in gs_KS()
        sys = gs_stdmap(T, K)
        x0s = gs_starts(T, seed + 1, m)
        Mmax = [gs_Mmax(sys, gs_points(sys, x0s[j], n)) for j in 1:m]
        runs = Dict((alg, tau) => gs_run(sys, x0s, n; k = 2, tau, seed, alg) for alg in gs_algs(), tau in (1, 2, 4))
        for ((alg, tau), r) in runs, j in 1:m
            f = gs_fin(r, j)
            track!(T, alg, f, 2, "AC5 K=$K $(gs_name(alg)) tau=$tau")
            track_o3!(T, Float64(abs(f.lambda[1] + f.lambda[2])), gs_theta(T, f.G_max, tau, Mmax[j]), "K=$K $(gs_name(alg)) tau=$tau")
        end
        ks = collect(keys(runs))
        for a in 1:length(ks), b in (a + 1):length(ks), j in 1:m
            fa, fb = gs_fin(runs[ks[a]], j), gs_fin(runs[ks[b]], j)
            th = gs_theta(T, fa.G_max, ks[a][2], Mmax[j]) + gs_theta(T, fb.G_max, ks[b][2], Mmax[j])
            for i in 1:2
                upd!(w, Float64(abs(fa.lambda[i] - fb.lambda[i])) / th,
                     "K=$K $(gs_name(ks[a][1])) tau=$(ks[a][2]) vs $(gs_name(ks[b][1])) tau=$(ks[b][2]) i=$i orbit $j")
            end
        end
    end
    return w
end

section("ac5") do
    seed = vs(5)
    for T in gs_types()
        m = sz(50, 5)
        w = ac5_std(T, seed, m)
        emitw("AC5", "pair bound |lambda_i^A - lambda_i^B| / (theta^A + theta^B), K panel x $m orbits, n=1000, 36 (variant,tau) pairs, i=1,2, $T", w)
    end
    for T in gs_types()
        A = gs_AS(T)
        runs = Dict{Tuple{String,Int},Any}()
        nm = Dict(gs_name(a) => a for a in gs_algs())
        for alg in gs_algs(), tau in (1, 2)
            res, err = gs_try_run(LinearSymplecticMap(A), [gs_zero(T, 4) for _ in 1:10], 600; k = 4, tau, seed = seed + 2, alg)
            if err === nothing
                runs[(gs_name(alg), tau)] = res
                for j in 1:10
                    f = gs_fin(res, j)
                    track!(T, alg, f, 4, "AC5 Y5 $(gs_name(alg)) tau=$tau")
                    track_o3!(T, Float64(abs(sum(f.lambda))), gs_theta_lin(T, A, tau), "Y5 $(gs_name(alg)) tau=$tau")
                end
            else
                push!(TRIPS, ("AC5 Y5", tname(T), tau, gs_name(alg)))
                say("AC5 | A_S k=4 n=600 $T $(gs_name(alg)) tau=$tau raised health error '$(err.what)' at step $(err.step) (reported, not scored when CGS)| - | - | INFO")
            end
        end
        w = Worst()
        ks = collect(keys(runs))
        for a in 1:length(ks), b in (a + 1):length(ks), j in 1:10
            fa, fb = gs_fin(runs[ks[a]], j), gs_fin(runs[ks[b]], j)
            th = gs_theta_lin(T, A, ks[a][2]) + gs_theta_lin(T, A, ks[b][2])
            for i in 1:4
                upd!(w, Float64(abs(fa.lambda[i] - fb.lambda[i])) / th, "$(ks[a]) vs $(ks[b]) i=$i frame $j")
            end
        end
        emitw("AC5", "A_S k=4 n=600 pair bound, $(length(ks)) completed (variant,tau) runs of 6, 10 frames, $T", w)
    end
end

# ---------------------------------------------------------------------------------------------------
# AC6 orthogonality
function gv_windows(alg, sys, x0::SVector{s,T}, W0::SMatrix{d,k,T}, N, tau) where {s,d,k,T}
    x, W = x0, W0
    Gw, loss, passes = Float64[], Float64[], Int[]
    for _ in 1:(N ÷ tau)
        for _ in 1:tau
            W = tangent(sys, x) * W
            x = forward(sys, x)
        end
        Q, R, c, st, _ = orthonormalize(alg, W)
        st == QR_OK || break
        push!(Gw, Float64(window_Gw(R))); push!(loss, Float64(Lyap.orth_loss(Q))); push!(passes, c)
        W = Q
    end
    return (Gw = Gw, loss = loss, passes = passes)
end

section("ac6") do
    seed = vs(6)
    for T in gs_types()
        w = Dict("RGS" => Worst(), "CGS" => Worst(), "MGS" => Worst())
        ratio = Dict("CGS" => Worst(), "MGS" => Worst())
        for K in gs_KS(), tau in (1, 2, 4), alg in gs_algs()
            sys = gs_stdmap(T, K)
            m = sz(20, 4)
            x0s = gs_starts(T, seed + 1, m)
            res = gs_run(sys, x0s, 1000; k = 2, tau, seed, alg)
            for j in 1:m
                f = gs_fin(res, j)
                upd!(w[gs_name(alg)], Float64(f.orth_max) / gs_orth_bound(alg, T, 2, f.G_max), "K=$K tau=$tau")
                track!(T, alg, f, 2, "AC6 K=$K tau=$tau $(gs_name(alg))")
            end
            if !(alg isa RGS)   # measured loss / (eps G_w) per window
                for j in 1:min(m, 4)
                    wd = gv_windows(alg, sys, x0s[j], initial_deviations(seed, j, Val(2), Val(2), T), 1000, tau)
                    isempty(wd.loss) || upd!(ratio[gs_name(alg)], maximum(wd.loss ./ (gs_eps(T) .* wd.Gw)), "K=$K tau=$tau")
                end
            end
        end
        emitw("AC6", "(a) RGS orth_max / (8 k eps), K panel x tau {1,2,4} x 20 orbits, n=1000, $T", w["RGS"])
        emitw("AC6", "(b) CGS orth_max / (8 k (G_max+1) eps), same runs, $T", w["CGS"])
        emitw("AC6", "(b) MGS orth_max / (8 k (G_max+1) eps), same runs, $T", w["MGS"])
        for a in ("CGS", "MGS")
            emit("AC6", "recorded: max loss / (eps G_w) per window, $a, k=2, [worst at $(ratio[a].at)], $T", ratio[a].v, "<= 8 k = 16 (spec scratch 1.05)", pf(ratio[a].v <= 16))
        end
    end
    # (c) discrimination on A_S, k = 4, Float64, 20 seeded frames
    sys = LinearSymplecticMap(gs_AS(Float64))
    x0 = [gs_zero(Float64, 4)]
    m = sz(20, 4)
    raised = Dict("CGS" => 0, "MGS" => 0, "RGS" => 0)
    minc, minm = Inf, Inf
    order1 = 0; minc1, minm1 = Inf, Inf
    for j in 1:m
        W0 = [initial_deviations(seed + 2, j, Val(4), Val(4), Float64)]
        loss3 = Dict(gs_name(a) => Float64(gs_fin(gs_run(sys, x0, 600; k = 4, tau = 3, seed, alg = a, tol = Inf, W0), 1).orth_max) for a in gs_algs())
        minc = min(minc, loss3["CGS"] / loss3["MGS"]); minm = min(minm, loss3["MGS"] / loss3["RGS"])
        for a in gs_algs()
            e = gs_health_error(() -> gs_run(sys, x0, 600; k = 4, tau = 3, seed, alg = a, W0))
            raised[gs_name(a)] += (e isa NumericalHealthError && e.what == "orthogonality") ? 1 : 0
        end
        loss1 = Dict(gs_name(a) => Float64(gs_fin(gs_run(sys, x0, 200; k = 4, tau = 1, seed, alg = a, tol = Inf, W0), 1).orth_max) for a in gs_algs())
        order1 += (loss1["CGS"] > loss1["MGS"] > loss1["RGS"]) ? 1 : 0
        minc1 = min(minc1, loss1["CGS"] / loss1["MGS"]); minm1 = min(minm1, loss1["MGS"] / loss1["RGS"])
    end
    emit("AC6", "(c) A_S k=4 tau=3: CGS orthogonality errors raised (of $m frames), key 8", raised["CGS"], "= $m", pf(raised["CGS"] == m))
    emit("AC6", "(c) A_S k=4 tau=3: MGS / RGS orthogonality errors raised (of $m frames each)", string(raised["MGS"], " / ", raised["RGS"]), "0 / 0", pf(raised["MGS"] == 0 && raised["RGS"] == 0))
    emit("AC6", "(c) tau=3: min over frames of loss(CGS)/loss(MGS)", minc, ">= 100 (spec scratch 407)", pf(minc >= 100))
    emit("AC6", "(c) tau=3: min over frames of loss(MGS)/loss(RGS)", minm, ">= 100 (spec scratch 2000)", pf(minm >= 100))
    emit("AC6", "(c) tau=1: frames (of $m) with loss(CGS) > loss(MGS) > loss(RGS)", order1, "= $m", pf(order1 == m))
    emit("AC6", "(c) tau=1: min loss(CGS)/loss(MGS) and min loss(MGS)/loss(RGS) (recorded; spec scratch 3.2 and 7.3)", string(fmt(minc1), " / ", fmt(minm1)), "> 1", pf(minc1 > 1 && minm1 > 1))
end

# ---------------------------------------------------------------------------------------------------
# AC7 breakdown
function ac7_frames(T, pbad, alg, seed, n)
    sysb = LinearSymplecticMap(gs_Apow(T, pbad)); A = gs_Apow(T, pbad)
    r = QRSpectrum(alg)
    raised = fields = npass = ncond = 0
    minegw = Inf
    for j in 1:n
        W0 = initial_deviations(seed, j, Val(2), Val(2), T)
        err = gs_health_error(() -> propagate_orbit(sysb, gs_zero(T, 2), W0, 1, 1, (r,); orbit = j))
        if err isa NumericalHealthError && err.what == "orthonormalisation breakdown"
            raised += 1
            fields += (err.orbit == j && err.step == 1 && err.backend === :cpu) ? 1 : 0
            err.threshold == Float64(Lyap.MAX_PASSES) ? (npass += 1) : (ncond += 1)
        end
        _, R, _, st, _ = orthonormalize(alg, A * W0)
        st == QR_OK && (minegw = min(minegw, Float64(window_Gw(R)) * gs_eps(T)))
    end
    return (raised = raised, fields = fields, npass = npass, ncond = ncond, minegw = minegw)
end

section("ac7") do
    seed = vs(7)
    n = sz(200, 20)
    for (T, pbad, pok) in ((Float64, 20, 16), (Double64, 40, 32)), alg in gs_algs()
        r = ac7_frames(T, pbad, alg, seed, n)
        pb = pbad
        if T === Double64 && r.raised < n    # spec fallback: the power 44
            say("AC7 | NOTE: Double64 A^40 raised in only $(r.raised) of $n frames for $(gs_name(alg)); spec fallback power 44 also run | - | - | NOTE")
            r = ac7_frames(T, 44, alg, seed, n); pb = 44
        end
        emit("AC7", "breakdown raised at handoff 1 with orbit, step=1, backend=:cpu: A^$pb $T $(gs_name(alg)); frames raised / correct fields (of $n)", string(r.raised, " / ", r.fields), "$n / $n", pf(r.raised == n && r.fields == n))
        emit("AC7", "  cause split: fifth-pass frames / conditioning-rule frames; min eps_eff G_w over the frames whose factorisation completed", string(r.npass, " / ", r.ncond, " ; ", fmt(r.minegw)), "min eps G_w >= 1/2 expected on completed frames (spec scratch 0.99 CGS/MGS, 1.74 RGS)", r.minegw >= 0.5 ? "PASS" : "FAIL")
        syso = LinearSymplecticMap(gs_Apow(T, pok))
        res = gs_serial(syso, [gs_zero(T, 2) for _ in 1:20], sz(200, 20); k = 2, tau = 1, seed, reducers = (QRSpectrum(alg),))
        G = maximum(j -> Float64(gs_fin(res, j).G_max), 1:20)
        ratio = maximum(j -> Float64(gs_fin(res, j).orth_max) / gs_orth_bound(alg, T, 2, G), 1:20)
        emit("AC7", "A^$pok no breakdown, n=200, 20 frames: orth_max / bound, G_max=$(fmt(G)), eps G_w ~ $(fmt(G * gs_eps(T))), $T $(gs_name(alg))", ratio, "<= 1", pf(ratio <= 1))
    end
    ok = true
    syso = LinearSymplecticMap(gs_Apow(Float64, 16)); x0s = [gs_zero(Float64, 2) for _ in 1:10]
    rc, rm = (gs_run(syso, x0s, 200; k = 2, tau = 1, seed, alg) for alg in (CGS(), MGS()))
    emit("AC7", "I6: CGS and MGS finals and Q bitwise equal at k=2 on A^16, Float64", isequal(rc.finals, rm.finals) && isequal(rc.W, rm.W), "true", pf(isequal(rc.finals, rm.finals) && isequal(rc.W, rm.W)))
    for T in gs_types(), alg in gs_algs()
        r = QRSpectrum(alg)
        sys4 = LinearSymplecticMap(one(SMatrix{4,4,T})); x4 = gs_zero(T, 4)
        okc = true; vals = String[]
        for (zc, k) in ((3, 4), (1, 4))
            cols = ntuple(c -> c == zc ? zero(SVector{4,T}) : SVector{4,T}(ntuple(i -> i == c ? one(T) : zero(T), Val(4))), Val(k))
            W0 = SMatrix{4,k,T}(reduce(vcat, cols)...)
            e1 = gs_health_error(() -> gs_serial(sys4, [x4], 2; k, tau = 2, seed, reducers = (r,), W0 = [W0]))
            st = reducer_init(r, T, Val(4), Val(k))
            e2 = gs_health_error(() -> reducer_update(r, st, 5, x4, W0))
            okc &= e1 isa NumericalHealthError && e1.what == "orthonormalisation breakdown" && e1.column == zc && e1.step == 2 && e1.orbit == 1 && e1.backend === :cpu
            okc &= e2 isa NumericalHealthError && e2.what == "orthonormalisation breakdown" && e2.column == zc && e2.step == 5 && e2.orbit == 0
            e2 isa NumericalHealthError && push!(vals, "col$zc value=$(fmt(e2.value)) threshold=$(fmt(e2.threshold))")
        end
        Wn = SMatrix{4,4,T}(1, 0, 0, 0, 0, 1, 0, 0, NaN, 0, 0, 0, 0, 0, 0, 1)
        e3 = gs_health_error(() -> gs_serial(sys4, [x4], 2; k = 4, tau = 2, seed, reducers = (r,), W0 = [Wn]))
        st = reducer_init(r, T, Val(4), Val(4))
        e4 = gs_health_error(() -> reducer_update(r, st, 2, x4, Wn))
        okn = e3 isa NumericalHealthError && e3.what == "deviation matrix" && e3.step == 2 &&
              e4 isa NumericalHealthError && e4.what == "orthonormalisation breakdown" && e4.column == 3 && isinf(e4.threshold)
        emit("AC7", "injected zero column 3 and 1 (propagate and reducer_update): kind, column, step, orbit correct; $T $(gs_name(alg)); error fields: $(join(vals, "; "))", okc, "true", pf(okc))
        emit("AC7", "injected NaN: propagate fires '$(e3 isa NumericalHealthError ? e3.what : "none")' first, reducer_update alone '$(e4 isa NumericalHealthError ? e4.what : "none")' (column $(e4 isa NumericalHealthError ? e4.column : -1), threshold $(e4 isa NumericalHealthError ? fmt(e4.threshold) : "-")), $T $(gs_name(alg))", okn, "true", pf(okn))
    end
end

# ---------------------------------------------------------------------------------------------------
# AC8 wrapped state
function gs_try_run_wrapped(sys, x0s, N, seed)
    try
        return propagate(sys, x0s, N; k = 2, tau = 1, seed, reducers = (QRSpectrum(RGS()),)), nothing
    catch e
        e isa NumericalHealthError ? (nothing, e) : rethrow()
    end
end

section("ac8") do
    seed = vs(8)
    for T in gs_types()
        mism = Dict(nm => 0 for nm in ("stream vs manual loop", "x, W vs manual", "S_i vs orthonormalize+add loop"))
        total = 0; w = Worst()
        for K in (2.0, 5.0), tau in (1, 3), alg in gs_algs()
            sys = WrappedStandardMap(T(K))
            rng = StableRNG(seed + 1)
            m = sz(20, 4)
            x0s = [SVector{4,T}(T(TWOPI * rand(rng)), T(TWOPI * rand(rng)), zero(T), zero(T)) for _ in 1:m]
            res = propagate(sys, x0s, 1200; k = 2, tau, seed, reducers = (QRSpectrum(alg),))
            r = QRSpectrum(alg)
            for j in 1:m
                W0 = initial_deviations(seed, j, Val(2), Val(2), T)
                xm, Wm, st = gs_manual_loop(r, sys, x0s[j], W0, 1200, tau)
                acc = gs_s_loop(alg, sys, x0s[j], W0, 1200, tau)
                total += 1
                isequal(gs_fin(res, j), reducer_final(r, st)) || (mism["stream vs manual loop"] += 1)
                (isequal(res.x[j], xm) && isequal(res.W[j], Wm) && length(res.x[j]) == 4) || (mism["x, W vs manual"] += 1)
                all(i -> isequal(gs_fin(res, j).sum[i], acc[i].s) && isequal(gs_fin(res, j).comp[i], acc[i].c), 1:2) || (mism["S_i vs orthonormalize+add loop"] += 1)
                track!(T, alg, gs_fin(res, j), 2, "AC8 K=$K tau=$tau $(gs_name(alg))")
            end
        end
        for (nm, c) in mism
            emit("AC8", "wrapped state (length 4), K in {2,5}, tau in {1,3}, 3 variants: orbits (of $total) with mismatch, $nm, $T", c, "0", pf(c == 0))
        end
        sys = WrappedStandardMap(T(5.0))
        rng = StableRNG(seed + 2)
        nlong = sz(10^6, 10^4)
        x0s = [SVector{4,T}(T(TWOPI * rand(rng)), T(TWOPI * rand(rng)), zero(T), zero(T)) for _ in 1:8]
        t = @elapsed ((res, err) = gs_try_run_wrapped(sys, x0s, nlong, seed))
        if err !== nothing
            emit("AC8", "$nlong-step run raised health error '$(err.what)' at step $(err.step), $T", "error", "no error", "FAIL")
        else
            for j in 1:8
                f = gs_fin(res, j)
                upd!(w, Float64(abs(f.lambda[1] + f.lambda[2])) / gs_theta(T, f.G_max, 1, 6.0), "orbit $j")
                track_o3!(T, Float64(abs(f.lambda[1] + f.lambda[2])), gs_theta(T, f.G_max, 1, 6.0), "AC8 orbit $j")
                track!(T, RGS(), f, 2, "AC8 long orbit $j")
            end
            fin = all(j -> all(isfinite, gs_fin(res, j).lambda) && all(isfinite, res.W[j]), 1:8)
            emitw("AC8", "$nlong steps, K=5, 8 orbits, tau=1, RGS: |lambda_1 + lambda_2| / theta_sum (M_max = 1+K), all finite=$fin, no health error, $(fmt(t)) s, $T", w; verdict = pf(w.v <= 1 && fin))
        end
        W0 = initial_deviations(seed, 1, Val(2), Val(2), T)
        for alg in gs_algs()
            b = gs_alloc_bytes(sys, x0s[1], W0, 10^4, 1, (QRSpectrum(alg),))
            emit("AC8", "bytes allocated by propagate_orbit after warm-up (wrapped K=5, 1e4 steps), $T $(gs_name(alg))", b, "0", pf(b == 0))
        end
    end
end

# ---------------------------------------------------------------------------------------------------
# AC9 thread and chunk invariance, zero allocation
section("ac9") do
    seed = vs(9)
    N = 1000
    for T in gs_types()
        sys = gs_stdmap(T, :twopi)
        for alg in gs_algs()
            n_orb = sz(4096, 256)
            x0s = gs_starts(T, seed + 1, n_orb)
            W0 = [initial_deviations(seed, j, Val(2), Val(2), T) for j in 1:n_orb]
            reducers = (QRSpectrum(alg),)
            serial = [propagate_orbit(sys, x0s[j], W0[j], N, 1, reducers; orbit = j) for j in 1:n_orb]
            xs, Ws = [s[1] for s in serial], [s[2] for s in serial]
            fs = [map(reducer_final, reducers, s[3]) for s in serial]
            for chunk in (64, 1000, nothing)
                res = propagate(sys, x0s, N; k = 2, tau = 1, seed, reducers, chunk, W0)
                bad = count(j -> !(isequal(res.x[j], xs[j]) && isequal(res.W[j], Ws[j]) && isequal(res.finals[j], fs[j])), 1:n_orb)
                emit("AC9", "propagate ($(Threads.nthreads()) threads, chunk=$(chunk === nothing ? "default" : chunk)) vs serial propagate_orbit: orbits (of $n_orb) with x, W, S_i, comp, record not bitwise equal, K=2pi n=1000, $T $(gs_name(alg))", bad, "0", pf(bad == 0))
            end
            for j in 1:min(n_orb, 64)
                track!(T, alg, fs[j][1], 2, "AC9 $(gs_name(alg))")
            end
        end
        sys4 = LinearSymplecticMap(gs_AS(T))
        W4 = initial_deviations(seed, 1, Val(4), Val(4), T)
        for alg in gs_algs()
            b = gs_alloc_bytes(sys4, gs_zero(T, 4), W4, 10^4, 1, (QRSpectrum(alg),))
            emit("AC9", "bytes allocated by propagate_orbit after warm-up, A_S d=k=4, 1e4 steps, $T $(gs_name(alg))", b, "0", pf(b == 0))
        end
    end
    say("# AC9 note: this run has $(Threads.nthreads()) threads; the spec compares 1 and 4 threads. The comparison here is serial propagate_orbit versus the threaded propagate at three chunk sizes; the tangent-core AC6 child-process test covers 1 versus 4 threads.")
end

# ---------------------------------------------------------------------------------------------------
# AC6 summary over all runs of AC1, AC4, AC5, AC6, AC8, AC9 (needs those sections in the same process)
section("ac6sum") do
    for ((t, cls), w) in sort(collect(ORTH); by = first)
        emit("AC6", "(a)/(b) $cls: max orth_max / bound over all runs of AC1, AC4, AC5, AC6, AC8, AC9 [worst at $(w.at)], $t", w.v, "<= 1", pf(w.v <= 1))
    end
    for (t, h) in sort(collect(PASSES); by = first)
        say("# RGS pass histogram (count_max of final records, all RGS runs above) $(t[1]): maximum passes $(t[2]): $h records")
    end
    for (t, w) in O3
        emit("O3", "max |sum_i lambda_i| / (eps (G_max/tau + M_max)) over AC5 and AC8 runs [worst at $(w.at)], $t", w.v, "<= 4 (diagnostics c; spec scratch 0.015 / 0.056)", pf(w.v <= 4))
    end
    for tr in TRIPS
        say("# CGS orthogonality trip recorded: $tr")
    end
end

# ---------------------------------------------------------------------------------------------------
# AC10 agreement with Benettin (k = 1)
function gv_gfirst(sys, x0, tau)
    P = tangent(sys, x0)
    x = forward(sys, x0)
    for _ in 2:tau
        P = tangent(sys, x) * P
        x = forward(sys, x)
    end
    return Float64(sum(abs2, P))
end

section("ac10") do
    seed = vs(10)
    n = sz(10^4, 10^3)
    for T in gs_types()
        w1, w2 = Worst(), Worst()
        nbit = 0; ntot = 0; dmax = 0.0
        for K in gs_KS(), tau in (1, 4)
            sys = gs_stdmap(T, K)
            m = sz(50, 5)
            x0s = gs_starts(T, seed + 1, m)
            rb = propagate(sys, x0s, n; k = 1, tau, seed, reducers = (Benettin(),))
            for alg in gs_algs()
                r1 = gs_run(sys, x0s, n; k = 1, tau, seed, alg)
                r2 = gs_run(sys, x0s[1:min(m, 10)], n; k = 2, tau, seed, alg)
                for j in 1:m
                    X = rb.finals[j][1]
                    bound = 2 * gs_BS(T, n, gv_gfirst(sys, x0s[j], tau), K) / n
                    f1 = gs_fin(r1, j)
                    d1 = Float64(abs(f1.lambda[1] - X.lambda[1]))
                    at = "K=$K tau=$tau $(gs_name(alg))"
                    upd!(w1, d1 / bound, at)
                    dmax = max(dmax, d1)
                    ntot += 1
                    nbit += (isequal(f1.sum[1], X.sum[1]) && isequal(f1.comp[1], X.comp[1])) ? 0 : 1
                    j <= 10 && upd!(w2, Float64(abs(gs_fin(r2, j).lambda[1] - X.lambda[1])) / bound, at)
                end
            end
        end
        emitw("AC10", "k=1: |lambda_1^GS - X| / ((B_S^GS + B_S^L1)/n), K panel x tau {1,4} x $(sz(50, 5)) orbits, n=$n, all variants, $T", w1)
        emit("AC10", "k=1: runs (of $ntot) whose sum or compensation is not bitwise equal to Benettin; realised max |lambda_1^GS - X| = $(fmt(dmax)), $T", nbit, "0", pf(nbit == 0))
        emitw("AC10", "k=2 same first column, first 10 orbits: |lambda_1 - X| / bound, $T", w2)
    end
end

# ---------------------------------------------------------------------------------------------------
# AC11 volume and sum identities, and the G8 sandwich
function gv_det(M::Matrix{BigFloat})
    n = size(M, 1); A = copy(M); d = one(BigFloat)
    for c in 1:n
        pv = argmax([abs(A[r, c]) for r in c:n]) + c - 1
        if pv != c
            A[[c, pv], :] = A[[pv, c], :]; d = -d
        end
        d *= A[c, c]
        for r in (c + 1):n
            f = A[r, c] / A[c, c]
            for k in c:n
                A[r, k] -= f * A[c, k]
            end
        end
    end
    return d
end

# Cyclic Jacobi eigensolver for a symmetric BigFloat matrix; returns eigenvalues in descending order and the eigenvectors.
function gv_jacobi(A0::Matrix{BigFloat}; sweeps = 40)
    n = size(A0, 1); A = copy(A0); V = Matrix{BigFloat}(I, n, n)
    for _ in 1:sweeps
        for p in 1:(n - 1), q in (p + 1):n
            iszero(A[p, q]) && continue
            theta = (A[q, q] - A[p, p]) / (2 * A[p, q])
            t = sign(theta == 0 ? one(theta) : theta) / (abs(theta) + sqrt(theta^2 + 1))
            c = 1 / sqrt(t^2 + 1); s = t * c
            for k in 1:n
                akp, akq = A[k, p], A[k, q]
                A[k, p] = c * akp - s * akq; A[k, q] = s * akp + c * akq
            end
            for k in 1:n
                apk, aqk = A[p, k], A[q, k]
                A[p, k] = c * apk - s * aqk; A[q, k] = s * apk + c * aqk
            end
            for k in 1:n
                vkp, vkq = V[k, p], V[k, q]
                V[k, p] = c * vkp - s * vkq; V[k, q] = s * vkp + c * vkq
            end
        end
    end
    ord = sortperm([A[i, i] for i in 1:n]; rev = true)
    return [A[i, i] for i in ord], V[:, ord]
end

# Bounds (1/n)[sum_{i<=p} ln sigma_i + ln|det(V_p' X)|] and (1/n) sum ln sigma_i of G8 for Y_n = A^n, X = first p columns of Q0; bits-bit arithmetic.
function gv_sandwich(Aint, n, Q0, bits)
    setprecision(BigFloat, bits) do
        Y = BigFloat.(gs_bigint_power(Aint, n))
        G = [sum(Y[l, i] * Y[l, j] for l in 1:4) for i in 1:4, j in 1:4]
        ev, V = gv_jacobi(G)
        resid = maximum(abs, G * V - V * Diagonal(ev)) / maximum(abs, G)
        orth = maximum(abs, V' * V - I)
        X = gs_bigmat(Q0)
        lnsig = log.(ev) ./ 2
        out = map(1:4) do p
            Mp = [sum(V[l, i] * X[l, j] for l in 1:4) for i in 1:p, j in 1:p]
            hi = sum(lnsig[1:p]); lo = hi + log(abs(gv_det(Mp)))
            (Float64(lo / n), Float64(hi / n))
        end
        return out, Float64(resid), Float64(orth)
    end
end

section("ac11") do
    seed = vs(11)
    n = 60
    m = sz(10, 3)
    for T in gs_types()
        A = gs_AS(T); Aint = gs_intmat(A)
        sys = LinearSymplecticMap(A)
        Q0s = [initial_deviations(seed, j, Val(4), Val(4), T) for j in 1:m]
        refs = [[gs_exact_logvol(Aint, n, Q0s[j], p) / n for p in 1:4] for j in 1:m]
        wv, wt = Worst(), Worst()
        for alg in gs_algs(), tau in (1, 2, 3)
            res, err = gs_try_run(sys, [gs_zero(T, 4) for _ in 1:m], n; k = 4, tau, seed, alg)
            if err !== nothing
                push!(TRIPS, ("AC11", tname(T), tau, gs_name(alg)))
                say("AC11 | A_S k=4 n=60 $T $(gs_name(alg)) tau=$tau raised health error '$(err.what)' at step $(err.step) (reported, not scored when CGS) | - | - | INFO")
                continue
            end
            th = gs_theta_lin(T, A, tau)
            for j in 1:m
                f = gs_fin(res, j)
                at = "$(gs_name(alg)) tau=$tau"
                for p in 1:4
                    upd!(wv, gs_absdiff(sum(f.lambda[1:p]), refs[j][p]) / (p * th), at * " p=$p")
                end
                upd!(wt, Float64(abs(sum(f.lambda))) / th, at)
                track_o3!(T, Float64(abs(sum(f.lambda))), th, "Y5 $at")
                track!(T, alg, f, 4, "AC11 $at")
            end
        end
        emitw("AC11", "|sum_{i<=p} lambda_i - (1/n) ln vol_p| / (p theta_sum), p=1..4, A_S n=60, $m frames, tau in {1,2,3}, completed runs, $T", wv)
        emitw("AC11", "|sum_i lambda_i| / theta_sum, $T", wt)
    end
    # G8 sandwich, n = 12, 512 bits
    n8 = 12
    for T in gs_types()
        A = gs_AS(T); Aint = gs_intmat(A)
        sys = LinearSymplecticMap(A)
        viol, resmax, orthmax, minwidth = Worst(), 0.0, 0.0, Inf
        exactviol = Worst()
        for j in 1:sz(10, 3)
            Q0 = initial_deviations(seed + 1, j, Val(4), Val(4), T)
            bounds, resid, orth = gv_sandwich(Aint, n8, Q0, 512)
            resmax = max(resmax, resid); orthmax = max(orthmax, orth)
            for p in 1:4
                lo, hi = bounds[p]
                ex = Float64(gs_exact_logvol(Aint, n8, Q0, p) / n8)
                upd!(exactviol, max(lo - ex, ex - hi), "p=$p frame $j")
            end
            for alg in gs_algs(), tau in (1, 2, 3)
                res, err = gs_try_run(sys, [gs_zero(T, 4)], n8; k = 4, tau, seed = seed + 1, alg, W0 = [Q0])
                err === nothing || continue
                f = gs_fin(res, 1); th = gs_theta_lin(T, A, tau)
                for p in 1:4
                    s = Float64(sum(f.lambda[1:p])); lo, hi = bounds[p]
                    upd!(viol, max(lo - s, s - hi) / (p * th), "$(gs_name(alg)) tau=$tau p=$p frame $j")
                    p < 4 && (minwidth = min(minwidth, (hi - lo) / (p * th)))
                end
            end
        end
        say("# AC11 sandwich: Jacobi eigen-residual max |Y'Y V - V L| / max|Y'Y| = $(fmt(resmax)), max |V'V - I| = $(fmt(orthmax)) (512 bits, n=12), $T")
        emitw("AC11", "G8 sandwich: (lower bound - sum_p lambda) or (sum_p lambda - upper bound), divided by p theta_sum (negative = inside), n=12, 512 bits, A_S, $T", viol)
        emitw("AC11", "G8 sandwich: exact (1/n) ln vol_p, violation of the lower/upper bound (absolute, negative = inside), $T", exactviol; bound = 1e-12)
        emit("AC11", "G8 sandwich width (upper - lower) / (p theta_sum), smallest over p=1..3 (information: the sandwich is wider than the tolerance, so the check is not vacuous only if this is > 1), $T", minwidth, "> 1 (information)", minwidth > 1 ? "INFO" : "NOTE")
    end
end

# ---------------------------------------------------------------------------------------------------
# AC14 provenance
section("ac14") do
    seed = vs(14)
    for T in gs_types()
        bad = 0; total = 0
        for alg in gs_algs(), (tol, tau) in ((8.0, 1), (4.0, 3), (Inf, 3))
            r = QRSpectrum(alg; orthogonality_tol = tol)
            res = propagate(gs_stdmap(T, 2.0), gs_starts(T, seed + 1, 3), 12; k = 2, tau, seed, reducers = (r,))
            p = res.provenance
            checks = (p.estimator == "qr", p.orthonormalizer == gs_name(alg), p.orthogonality_tol === tol, p.k == 2, p.tau == tau, p.seed == seed,
                      p.precision == precision_label(T), p.backend === :cpu, p.device == "cpu", p.backend_fallback_reason == "",
                      !haskey(p, :iteration_path), !haskey(p, :implementation), !haskey(p, :d),
                      reducer_provenance(r) == (estimator = "qr", orthonormalizer = gs_name(alg), orthogonality_tol = tol))
            total += length(checks); bad += count(!, checks)
        end
        emit("AC14", "provenance record fields (estimator, orthonormalizer, orthogonality_tol, k, tau, seed, precision, backend, device, fallback reason, absent Householder keys): failed checks of $total, $T", bad, "0", pf(bad == 0))
    end
    emit("AC14", "QRSpectrum(RGS()).orthogonality_tol === 8.0 (key value)", QRSpectrum(RGS()).orthogonality_tol === 8.0, "true", pf(QRSpectrum(RGS()).orthogonality_tol === 8.0))
end

# ---------------------------------------------------------------------------------------------------
# AC12 (CPU vs GPU): not run by this script
section("ac12") do
    emit("AC12", "CPU vs GPU, RGS and MGS, S_i/n for i=1..k (tangent-core AC12 tolerances), Float64 and Double64", "not run by this script", "tangent-core AC12", "DEFERRED")
end

# ---------------------------------------------------------------------------------------------------
# Oracle: ChaosTools lyapunovspectrum (AC13, d = 4 case, R3)
# Facts (ChaosTools 3.5.4, src/chaosdetection/lyapunovs/lyapunovspectrum.jl, DynamicalSystemsBase 3.19.3 tangent_space.jl):
#   lyapunovspectrum(tands, N; u0) calls reinit!(tands, u0), whose Q0 defaults to default_deviations = identity columns, so a seeded Q0 given
#   to the TangentDynamicalSystem constructor is overwritten before the first step. Both sides therefore start from W_0 = I
#   (probed in the first lines of the oracle section). Each step: step!(tands), Q, R = qr(deviations), lambda_i += log(abs(R_ii)).
ct_rule(x, par, n) = (p = x[2] + par[1] * sin(x[1]); SVector(x[1] + p, p))
ct_jac(x, par, n) = (c = par[1] * cos(x[1]); SMatrix{2,2}(1 + c, c, 1.0, 1.0))
ct_tands(K, x0) = CT.TangentDynamicalSystem(CT.DeterministicIteratedMap(ct_rule, x0, [K]); J = ct_jac)
ct_spec(tands, x0, N) = SVector{2,Float64}(CT.lyapunovspectrum(tands, N; u0 = x0))
const I2 = SMatrix{2,2,Float64}(1.0, 0.0, 0.0, 1.0)
const I4 = SMatrix{4,4,Float64}(I)

function gv_regular_start(K, seed)     # S7-style screening: first seeded start whose lambda_1(1e5) <= 5e-4 (estimator under test, identity frame)
    rng = StableRNG(seed)
    sys = StandardMap(K)
    for c in 1:200
        x0 = SVector(TWOPI * rand(rng), TWOPI * rand(rng))
        f = gs_fin(gs_run(sys, [x0], 10^5; k = 2, tau = 1, seed, alg = RGS(), W0 = [I2]), 1)
        f.lambda[1] <= 5e-4 && return x0, c - 1
    end
    error("no regular start found")
end

section("oracle") do
    seed = vs(13)
    # probe: seeded Q0 is ignored by lyapunovspectrum
    let x0 = SVector(1.25, 2.5), K = 2.0
        Q0 = SMatrix{2,2}(0.6, 0.8, -0.8, 0.6)
        ds = CT.DeterministicIteratedMap(ct_rule, x0, [K])
        a = CT.lyapunovspectrum(CT.TangentDynamicalSystem(ds; J = ct_jac, Q0 = Q0), 1; u0 = x0)
        b = CT.lyapunovspectrum(CT.TangentDynamicalSystem(ds; J = ct_jac), 1; u0 = x0)
        M = ct_jac(x0, [K], 0) * Q0; R = qr(M).R
        emit("ORACLE", "lyapunovspectrum with seeded Q0 vs identity frame, N=1 (equal means Q0 is ignored); seeded-frame manual QR gives $(fmt(log(abs(R[1, 1]))))", maximum(abs, a .- b), "0 (finding)", "INFO")
        sys = StandardMap(K); xo = x0
        ds2 = CT.DeterministicIteratedMap(ct_rule, x0, [K])
        for _ in 1:1000
            CT.step!(ds2); xo = forward(sys, xo)
        end
        emit("ORACLE", "ChaosTools lifted rule vs Chaoslab forward, K=2, 1000 steps: max |difference| (same arithmetic)", maximum(abs, SVector{2,Float64}(CT.current_state(ds2)) - xo), "0 (information)", "INFO")
    end
    # AC13 regular
    K = 0.5
    x0, rejected = gv_regular_start(K, seed)
    tands = ct_tands(K, x0)
    lam_or = ct_spec(tands, x0, sz(10^5, 10^4))
    N = sz(10^5, 10^4)
    for alg in gs_algs()
        f = gs_fin(gs_run(StandardMap(K), [x0], N; k = 2, tau = 1, seed, alg, W0 = [I2]), 1)
        th = gs_theta(Float64, f.G_max, 1, 1 + K)
        emit("AC13", "regular K=0.5 (S7 screen, $rejected candidates rejected), 1 orbit, n=$N, identity frame, $(gs_name(alg)): |lambda_1 - lambda_1^oracle| (lambda_1 = $(fmt(f.lambda[1])), oracle $(fmt(lam_or[1])))", Float64(abs(f.lambda[1] - lam_or[1])), "<= 1e-4", pf(abs(f.lambda[1] - lam_or[1]) <= 1e-4))
        emit("AC13", "regular K=0.5, $(gs_name(alg)): |lambda_1 + lambda_2| / theta_sum (oracle sum $(fmt(lam_or[1] + lam_or[2])))", Float64(abs(f.lambda[1] + f.lambda[2])) / th, "<= 1", pf(abs(f.lambda[1] + f.lambda[2]) / th <= 1))
        emit("AC13", "regular K=0.5, $(gs_name(alg)): |lambda_2 - lambda_2^oracle| / (2 theta_sum) (pathwise, information beyond the spec)", Float64(abs(f.lambda[2] - lam_or[2])) / (2 * th), "<= 1", pf(abs(f.lambda[2] - lam_or[2]) / (2 * th) <= 1))
    end
    # AC13 chaotic ensemble
    K = TWOPI
    m = sz(1000, 50); N = sz(10^4, 10^3)
    x0s = gs_starts(Float64, seed + 1, m)
    tands = ct_tands(K, x0s[1])
    or = [ct_spec(tands, x0s[j], N) for j in 1:m]
    mean_se(v) = (mean(v), std(v) / sqrt(length(v)))
    for alg in gs_algs()
        res = gs_run(StandardMap(K), x0s, N; k = 2, tau = 1, seed, alg, W0 = [I2 for _ in 1:m])
        for i in 1:2
            a = [Float64(gs_fin(res, j).lambda[i]) for j in 1:m]; b = [or[j][i] for j in 1:m]
            m1, s1 = mean_se(a); m2, s2 = mean_se(b)
            z = abs(m1 - m2) / sqrt(s1^2 + s2^2 + 1e-300)
            emit("AC13", "chaotic K=2pi, $m orbits, n=$N, $(gs_name(alg)), lambda_$i: mean $(fmt(m1)) (SE $(fmt(s1))) vs oracle $(fmt(m2)) (SE $(fmt(s2))); |diff|/combined SE; max pathwise |diff| = $(fmt(maximum(abs.(a .- b))))", z, "<= 4", pf(z <= 4))
        end
        wsum = Worst()
        for j in 1:m
            f = gs_fin(res, j)
            upd!(wsum, Float64(abs(f.lambda[1] + f.lambda[2])) / gs_theta(Float64, f.G_max, 1, 1 + K), "orbit $j")
        end
        emitw("AC13", "chaotic K=2pi $(gs_name(alg)): |lambda_1 + lambda_2| / theta_sum (M_max = 1+K)", wsum)
    end
end

section("oracle4") do
    # d = 4 linear symplectic map A_S with known exponents, fixed point 0, identity frame, against the exact partial sums (G7) and ChaosTools
    A = gs_AS(Float64); Aint = gs_intmat(A)
    sys = LinearSymplecticMap(A)
    rule4(x, par, n) = A * x
    jac4(x, par, n) = A
    ds = CT.DeterministicIteratedMap(rule4, zero(SVector{4,Float64}), [0.0])
    tands = CT.TangentDynamicalSystem(ds; J = jac4)
    for n in (60, sz(600, 120))
        lam_or = SVector{4,Float64}(CT.lyapunovspectrum(tands, n; u0 = zero(SVector{4,Float64})))
        ex = [gs_exact_logvol(Aint, n, I4, p) / n for p in 1:4]
        th = gs_theta_lin(Float64, A, 1)
        wo = Worst()
        for p in 1:4
            upd!(wo, gs_absdiff(sum(lam_or[1:p]), ex[p]) / (p * th), "p=$p")
        end
        emitw("ORACLE-D4", "ChaosTools lyapunovspectrum, A_S d=4, n=$n, identity frame: |sum_{i<=p} lambda_i^CT - (1/n) ln vol_p| / (p theta_sum)", wo)
        for alg in gs_algs()
            f = gs_fin(gs_run(sys, [gs_zero(Float64, 4)], n; k = 4, tau = 1, seed = vs(13), alg, W0 = [I4]), 1)
            w, wp = Worst(), Worst()
            for p in 1:4
                upd!(w, gs_absdiff(sum(f.lambda[1:p]), ex[p]) / (p * th), "p=$p")
                upd!(wp, Float64(abs(sum(f.lambda[1:p]) - sum(lam_or[1:p]))) / (2 * p * th), "p=$p")
            end
            emitw("ORACLE-D4", "$(gs_name(alg)) tau=1, A_S d=4, n=$n: |sum_{i<=p} lambda_i - exact| / (p theta_sum)", w)
            emitw("ORACLE-D4", "$(gs_name(alg)) vs ChaosTools, n=$n: |sum_{i<=p} (lambda_i - lambda_i^CT)| / (2 p theta_sum)", wp)
            n == 60 || say("# ORACLE-D4 n=$n $(gs_name(alg)) lambda = $(join(map(fmt, f.lambda), ", ")); ChaosTools $(join(map(fmt, lam_or), ", ")); asymptotic +-0.9624236501, +-1.3169578969 (G12c, finite-time ln kappa(S)/n offset expected)")
        end
    end
end

# R3: ChaosTools k = 2 sum at 1e6 steps, lifted rule, next to theta_sum and the Chaoslab |lambda_1 + lambda_2| of the same orbit
section("r3") do
    seed = vs(15)
    N = sz(10^6, 10^4)
    x0s = gs_starts(Float64, seed + 1, 1)
    for K in (0.5, 0.9716354, 2.0, 5.0, 10.0)
        tands = ct_tands(K, x0s[1])
        lam_or = ct_spec(tands, x0s[1], N)
        sys = StandardMap(K)
        pts = gs_points(sys, x0s[1], N)
        Mmax = gs_Mmax(sys, pts)
        for alg in gs_algs()
            f = gs_fin(gs_run(sys, x0s, N; k = 2, tau = 1, seed, alg, W0 = [I2]), 1)
            th = gs_theta(Float64, f.G_max, 1, Mmax)
            track_o3!(Float64, Float64(abs(f.lambda[1] + f.lambda[2])), th, "R3 K=$K")
            emit("R3", "K=$K, n=$N, $(gs_name(alg)): ChaosTools lambda_1+lambda_2 = $(fmt(lam_or[1] + lam_or[2])); Chaoslab |lambda_1+lambda_2| = $(fmt(Float64(abs(f.lambda[1] + f.lambda[2])))); theta_sum = $(fmt(th)) (M_max = $(fmt(Mmax)), G_max = $(fmt(Float64(f.G_max)))); lambda_1 Chaoslab $(fmt(f.lambda[1])) vs CT $(fmt(lam_or[1]))", abs(lam_or[1] + lam_or[2]), "(recorded, no threshold)", "INFO")
        end
    end
    # loose sanity: lambda_1 at K = 10 against ln(K/2) (benchmarks.md), non-scored
    K = 10.0
    f = gs_fin(gs_run(StandardMap(K), x0s, N; k = 2, tau = 1, seed, alg = RGS(), W0 = [I2]), 1)
    emit("SANITY", "K=10 lambda_1 = $(fmt(f.lambda[1])) vs ln(K/2) = $(fmt(log(K / 2))) (single orbit; non-scored, risk K18)", Float64(f.lambda[1]) - log(K / 2), "none", "INFO")
end

# ---------------------------------------------------------------------------------------------------
# Impl-note deviations and code facts
section("deviations") do
    r = QRSpectrum(RGS())
    x4 = gs_zero(Float64, 4)
    # (1) the pivot is read from R[column, column] on a failed status
    Wz = SMatrix{4,4,Float64}(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1)
    Wn = SMatrix{4,4,Float64}(1, 0, 0, 0, 0, 1, 0, 0, NaN, 0, 0, 0, 0, 0, 0, 1)
    for alg in gs_algs()
        _, Rz, _, sz_, cz = orthonormalize(alg, Wz)
        _, Rn, _, sn, cn = orthonormalize(alg, Wn)
        ez = gs_health_error(() -> reducer_update(QRSpectrum(alg), reducer_init(QRSpectrum(alg), Float64, Val(4), Val(4)), 1, x4, Wz))
        en = gs_health_error(() -> reducer_update(QRSpectrum(alg), reducer_init(QRSpectrum(alg), Float64, Val(4), Val(4)), 1, x4, Wn))
        ok = sz_ == QR_PIVOT && cz == 3 && Rz[3, 3] == 0 && sn == QR_PIVOT && cn == 3 && isnan(Rn[3, 3]) &&
             ez.value == 0 && isnan(en.value) && isinf(en.threshold)
        emit("DEV1", "$(gs_name(alg)): failed status leaves the pivot in R[col,col] (zero column: R33 = $(Rz[3, 3]), error value $(ez.value); NaN column: error value $(en.value), threshold $(en.threshold))", ok, "true", pf(ok))
    end
    # (2) Double64: health error fields are Float64
    T = Double64
    rD = QRSpectrum(RGS())
    WzD = SMatrix{4,4,Double64}(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1)
    eD = gs_health_error(() -> reducer_update(rD, reducer_init(rD, T, Val(4), Val(4)), 1, gs_zero(T, 4), WzD))
    emit("DEV2", "Double64 health error field types (value, threshold)", string(typeof(eD.value), ", ", typeof(eD.threshold)), "Float64, Float64", pf(eD.value isa Float64 && eD.threshold isa Float64))
    # (3) G_max, orth_max, count_max, sums unchanged on failure
    for alg in gs_algs()
        rr = QRSpectrum(alg)
        st0 = reducer_init(rr, Float64, Val(4), Val(4))
        Wg = SMatrix{4,4,Float64}(1.0, 0.1, 0.0, 0.2, 0.3, 1.0, 0.1, 0.0, 0.0, 0.2, 1.0, 0.3, 0.1, 0.0, 0.2, 1.0)
        st1, _ = reducer_step(rr, st0, 1, x4, Wg)
        st2, _ = reducer_step(rr, st1, 2, x4, Wz)
        same = st2.G_max == st1.G_max && st2.orth_max == st1.orth_max && st2.count_max == st1.count_max && st2.sums == st1.sums && st2.n == st1.n && st2.code != 0 && st2.step == 2
        emit("DEV3", "$(gs_name(alg)): after a failing handoff G_max, orth_max, count_max, sums, n are those of the last good handoff; code=$(st2.code), step=$(st2.step)", same, "true", pf(same))
    end
    # (4) window_G equals _G_from_Gw(window_Gw) and the exact condition number at k = 2; upper bound at k = 3
    rng = StableRNG(vs(16))
    wg = Worst(); wb = Worst(); same = true
    for _ in 1:2000
        R2 = SMatrix{2,2,Float64}(10.0^(6 * (rand(rng) - 0.5)), 0.0, 10.0^(6 * (rand(rng) - 0.5)) * (rand(rng) - 0.5), 10.0^(6 * (rand(rng) - 0.5)))
        g = window_G(R2); c = cond(Matrix(R2))
        same &= g == Lyap._G_from_Gw(window_Gw(R2), Val(2))
        upd!(wg, abs(g - c) / c, "R=$(R2)")
        R3 = SMatrix{3,3,Float64}(10.0^(4 * rand(rng)), 0.0, 0.0, randn(rng), 10.0^(4 * rand(rng)), 0.0, randn(rng), randn(rng), 10.0^(4 * rand(rng)))
        upd!(wb, cond(Matrix(R3)) / window_G(R3), "k=3")
    end
    emit("DEV4", "window_G == _G_from_Gw(window_Gw) (2000 random 2x2)", same, "true", pf(same))
    emit("DEV4", "k=2: max relative difference of window_G and cond_2(R) (SVD) over 2000 random triangles with entries over 6 decades [worst at $(wg.at)]", wg.v, "<= 1e-6", pf(wg.v <= 1e-6))
    emit("DEV4", "k=3: max cond_2(R) / window_G (upper bound claim, must be <= 1+1e-12)", wb.v, "<= 1", pf(wb.v <= 1 + 1e-12))
    # (5) window_Gw rescale guard (the impl notes say there is none)
    R = SMatrix{2,2,Float64}(3.0, 0.0, 2.0, 0.5)
    g0 = window_Gw(R)
    emit("DEV5", "window_Gw scale invariance with the rescale guard: max |G_w(1e200 R)/G_w(R) - 1|, |G_w(1e-200 R)/G_w(R) - 1| (notes: 'no overflow guard'; code has one at lines 111-112)", max(abs(window_Gw(1e200 * R) / g0 - 1), abs(window_Gw(1e-200 * R) / g0 - 1)), "<= 1e-12", pf(max(abs(window_Gw(1e200 * R) / g0 - 1), abs(window_Gw(1e-200 * R) / g0 - 1)) <= 1e-12))
    # (6) stale header comment in orthonormalize.jl
    hdr = readlines(joinpath(ROOT, "src", "Lyapunov", "orthonormalize.jl"))[3]
    emit("DEV6", "orthonormalize.jl line 3 reads: $hdr", occursin("MMatrix", hdr), "(information: stale if it mentions MMatrix scratch; the notes say it was replaced)", "NOTE")
    # (7) allocations of orthonormalize
    for T in gs_types(), alg in gs_algs()
        W = SMatrix{4,4,T}(ntuple(i -> T(1) / T(mod(i - 1, 4) + div(i - 1, 4) + 1), Val(16)))
        b = gs_alloc_orthonormalize(alg, W)
        emit("DEV7", "bytes allocated by orthonormalize (4x4) after warm-up, $T $(gs_name(alg))", b, "0", pf(b == 0))
    end
    # (8) k = 1: orth_max stays 0
    f = gs_fin(gs_run(gs_stdmap(Float64, 2.0), gs_starts(Float64, vs(16), 3), 100; k = 1, tau = 1, seed = vs(16), alg = CGS()), 1)
    emit("DEV8", "k=1: orth_max", f.orth_max, "0", pf(f.orth_max == 0))
end

# ---------------------------------------------------------------------------------------------------
# P1 vacuous-check band and loss / (eps G_w) on the powers of the cat map
section("p1band") do
    seed = vs(17)
    for (T, p) in ((Float64, 16), (Double64, 32), (Float64, 10))
        sys = LinearSymplecticMap(gs_Apow(T, p))
        for alg in (CGS(), MGS(), RGS())
            wd = gv_windows(alg, sys, gs_zero(T, 2), initial_deviations(seed, 1, Val(2), Val(2), T), 200, 1)
            isempty(wd.loss) && continue
            Gw = maximum(wd.Gw); loss = maximum(wd.loss)
            thr = alg isa RGS ? 8 * 2 * gs_eps(T) : 8 * 2 * gs_eps(T) * Gw
            say("P1 | A^$p $T $(gs_name(alg)): G_w = $(fmt(Gw)) (eps G_w = $(fmt(Gw * gs_eps(T))); band eps^-1/2 < G_w < 1/(2 eps) is $(fmt(gs_eps(T)^-0.5)) .. $(fmt(0.5 / gs_eps(T)))), max loss = $(fmt(loss)), threshold_orth = $(fmt(thr)), loss/(eps G_w) = $(fmt(loss / (gs_eps(T) * Gw))) | recorded | INFO")
        end
    end
end

# ---------------------------------------------------------------------------------------------------
# Figures (CairoMakie last: its first call is slow)
section("figures") do
    seed = vs(18)
    cols = Dict("CGS" => :crimson, "MGS" => :darkorange, "RGS" => :steelblue)
    mk = Dict("CGS" => :circle, "MGS" => :utriangle, "RGS" => :rect)
    fig = CM.Figure(size = (1250, 520))
    panels = (
        ("standard map K = 2 pi, k = 2, tau = 1..8, 200 windows", 2, [(alg, tau) for alg in gs_algs(), tau in 1:8],
         (alg, tau) -> gv_windows(alg, gs_stdmap(Float64, :twopi), gs_starts(Float64, seed + 1, 1)[1], initial_deviations(seed, 1, Val(2), Val(2), Float64), 200 * tau, tau)),
        ("A_S (Y5 pair), k = 4, tau = 1..4, 200 windows", 4, [(alg, tau) for alg in gs_algs(), tau in 1:4],
         (alg, tau) -> gv_windows(alg, LinearSymplecticMap(gs_AS(Float64)), gs_zero(Float64, 4), initial_deviations(seed, 2, Val(4), Val(4), Float64), 200 * tau, tau)))
    for (i, (title, k, cases, f)) in enumerate(panels)
        ax = CM.Axis(fig[1, i]; xscale = log10, yscale = log10, xlabel = "window G_w", ylabel = "orthogonality loss max |Q'Q - I|", title = title)
        done = Set{String}()
        for (alg, tau) in vec(cases)
            wd = f(alg, tau); isempty(wd.loss) && continue
            nm = gs_name(alg)
            CM.scatter!(ax, wd.Gw, max.(wd.loss, 1e-18); color = (cols[nm], 0.5), marker = mk[nm], markersize = 6, label = nm in done ? nothing : nm)
            push!(done, nm)
        end
        g = 10.0 .^ range(0, 16; length = 50)
        CM.lines!(ax, g, fill(8 * k * gs_eps(Float64), length(g)); color = :steelblue, linestyle = :dash, label = "8 k eps (RGS bound)")
        CM.lines!(ax, g, min.(8 * k * gs_eps(Float64) .* g, 1e2); color = :black, linestyle = :dash, label = "8 k G_w eps (CGS, MGS bound)")
        CM.ylims!(ax, 1e-17, 10); CM.xlims!(ax, 1, 1e16)
        i == 1 && CM.axislegend(ax; position = :lt, labelsize = 11)
    end
    CM.save(joinpath(OUT, "lyap-gram-schmidt-orthogonality.png"), fig)
    say("# figure saved: lyap-gram-schmidt-orthogonality.png (windows with a failed factorisation are omitted)")

    # oracle figure: (a) |lambda_1 - lambda_1^CT| at K = 2 pi against n for the three variants, (b) A_S d = 4 partial sum p = 2 error against the exact value
    fig2 = CM.Figure(size = (1250, 520))
    K = TWOPI; x0 = gs_starts(Float64, seed + 2, 1)[1]
    ns = sz([10, 30, 100, 300, 1000, 3000, 10^4, 3 * 10^4, 10^5], [10, 30, 100, 300, 1000])
    ax = CM.Axis(fig2[1, 1]; xscale = log10, yscale = log10, xlabel = "n (map steps)", ylabel = "|lambda_1 - lambda_1^ChaosTools|", title = "standard map K = 2 pi, one orbit, identity frame")
    tands = ct_tands(K, x0)
    orc = [ct_spec(tands, x0, n)[1] for n in ns]
    thetas = Float64[]
    for alg in gs_algs()
        d = Float64[]
        for (i, n) in enumerate(ns)
            f = gs_fin(gs_run(StandardMap(K), [x0], n; k = 2, tau = 1, seed, alg, W0 = [I2]), 1)
            push!(d, max(Float64(abs(f.lambda[1] - orc[i])), 1e-18)); alg isa RGS && push!(thetas, gs_theta(Float64, f.G_max, 1, 1 + K))
        end
        CM.scatterlines!(ax, ns, d; color = cols[gs_name(alg)], marker = mk[gs_name(alg)], label = gs_name(alg))
    end
    CM.lines!(ax, ns, thetas; color = :black, linestyle = :dash, label = "theta_sum (RGS run, M_max = 1 + K)")
    CM.axislegend(ax; position = :lb, labelsize = 11)
    A = gs_AS(Float64); Aint = gs_intmat(A)
    ns2 = sz([10, 20, 40, 80, 160, 320, 600], [10, 20, 40])
    ax2 = CM.Axis(fig2[1, 2]; xscale = log10, yscale = log10, xlabel = "n (map steps)", ylabel = "|lambda_1 + lambda_2 - (1/n) ln vol_2|", title = "A_S, d = 4, identity frame, tau = 1")
    ds = CT.DeterministicIteratedMap((x, par, n) -> A * x, zero(SVector{4,Float64}), [0.0])
    t4 = CT.TangentDynamicalSystem(ds; J = (x, par, n) -> A)
    ex = [gs_exact_logvol(Aint, n, I4, 2) / n for n in ns2]
    dct = [max(gs_absdiff(sum(CT.lyapunovspectrum(t4, n; u0 = zero(SVector{4,Float64}))[1:2]), ex[i]), 1e-18) for (i, n) in enumerate(ns2)]
    CM.scatterlines!(ax2, ns2, dct; color = :gray30, marker = :diamond, label = "ChaosTools")
    for alg in gs_algs()
        d = [max(gs_absdiff(sum(gs_fin(gs_run(LinearSymplecticMap(A), [gs_zero(Float64, 4)], n; k = 4, tau = 1, seed, alg, W0 = [I4]), 1).lambda[1:2]), ex[i]), 1e-18) for (i, n) in enumerate(ns2)]
        CM.scatterlines!(ax2, ns2, d; color = cols[gs_name(alg)], marker = mk[gs_name(alg)], label = gs_name(alg))
    end
    CM.lines!(ax2, ns2, fill(2 * gs_theta_lin(Float64, A, 1), length(ns2)); color = :black, linestyle = :dash, label = "2 theta_sum")
    CM.axislegend(ax2; position = :lb, labelsize = 11)
    CM.save(joinpath(OUT, "lyap-gram-schmidt-oracle.png"), fig2)
    say("# figure saved: lyap-gram-schmidt-oracle.png")
end

say(@sprintf("# total wall time %.1f s", time() - T0))
close(LOG)
