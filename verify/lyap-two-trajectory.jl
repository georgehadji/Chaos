# Stage 04 verification of module lyap-two-trajectory (propagate_shadow, separation helpers).
# Run from the repo root:  julia +release -t 2 --project=verify verify/lyap-two-trajectory.jl
# Spec: workspaces/module-build/stages/01-spec/output/lyap-two-trajectory-spec.md (criteria 1 to 17, Amendment D35).
# Section "suite": the package suite at full size (D38, L2_FULL_RUN = true), 568 assertions, about 4 minutes.
# Other sections: one printed line per criterion  AC<n> | quantity | value | bound | PASS/FAIL  with the MEASURED number, reusing the helpers of
# test/lyap-two-trajectory/l2_helpers.jl (one home for tolerance formulas). Values are "worst ratio error / bound" (pass when <= 1) unless stated.
# Seeds: L2_SEED + 45000 + 50 ac + sub (package tests use L2_SEED + 25000 .. 25999), so the samples differ from the package tests.
# L2_ONLY="AC6,oracle" runs only the named sections (debugging aid) and L2_LOG names the log file; the default runs everything.
using Printf, Test, Random, LinearAlgebra, Statistics
using Chaoslab, StaticArrays, DoubleFloats, StableRNGs
import ChaosTools as CT
import CairoMakie as CM

const L2_FULL_RUN = true
const ROOT = dirname(@__DIR__)
const TESTDIR = joinpath(ROOT, "test", "lyap-two-trajectory")
include(joinpath(TESTDIR, "l2_helpers.jl"))

vseed(ac, sub = 0) = L2_SEED + 45000 + 50 * ac + sub

# Thread-count invariance (AC11): the same X2 series hashed in this process and in child processes with other thread counts.
function l2_thread_hash()
    m, n, d0 = 4096, 1000, 1e-8
    sys = WrappedStandardMap(L2_TWO_PI)
    x0 = [l2_state(:wrapped, Float64, z) for z in l2_starts(StableRNG(vseed(11)), m)]
    res = propagate_shadow(sys, x0, n; d0 = d0, tau = 1, seed = vseed(11))
    return hash((reinterpret(UInt64, [f[1].X2 for f in res.finals]), reinterpret(UInt64, [f[1].sum for f in res.finals]), res.x, res.W))
end
if get(ENV, "L2_CHILD", "") == "1"
    println("HASH ", l2_thread_hash(), " THREADS ", Threads.nthreads())
    exit(0)
end

const OUT = joinpath(ROOT, "workspaces", "module-build", "stages", "04-verify", "output")
mkpath(OUT)
const ONLY = filter(!isempty, split(get(ENV, "L2_ONLY", ""), ","))
const LOG = open(get(ENV, "L2_LOG", joinpath(OUT, "lyap-two-trajectory-verify-log.txt")), "w")
const VERDICTS = String[]

fmt(x::AbstractFloat) = @sprintf("%.4g", x)
fmt(x) = string(x)
function say(s)
    println(s); println(LOG, s); flush(stdout); flush(LOG)
end
function emit(ac, quantity, value, bound, verdict)
    push!(VERDICTS, verdict)
    say(string(ac, " | ", quantity, " | ", fmt(value), " | ", bound, " | ", verdict))
end
pf(ok::Bool) = ok ? "PASS" : "FAIL"
ratio_line(ac, quantity, r) = emit(ac, quantity, r, "<= 1", pf(r <= 1))

function section(f, name)
    isempty(ONLY) || name in ONLY || return
    t0 = time()
    try
        f()
    catch e
        say(string(name, " | ERROR | ", sprint(showerror, e)[1:min(end, 300)], " | - | FAIL"))
        push!(VERDICTS, "FAIL")
    end
    say(@sprintf("# section %s took %.1f s", name, time() - t0))
end

function l2_alloc(sys, x0, u0, N, tau, d0)
    propagate_shadow_orbit(sys, x0, u0, N, tau, d0)
    return @allocated propagate_shadow_orbit(sys, x0, u0, N, tau, d0)
end

const T0 = time()
const DATA = Dict{Any, Any}()
const EPS64 = Float64(effective_eps(Float64))
const EPSDD = Float64(effective_eps(Double64))
say("# lyap-two-trajectory verification, julia $(VERSION), threads $(Threads.nthreads()), base seed $L2_SEED (verify offset 45000), time $(round(time())) (unix)")
say("# timings in this log were taken on a loaded machine (other agents were running concurrently)")

# ---------------------------------------------------------------------------------------------------
# The package suite at full size (D38)
section("suite") do
    local ts
    t = @elapsed try
        ts = @testset verbose = true "lyap-two-trajectory verify" begin
            include(joinpath(TESTDIR, "runtests.jl"))
        end
    catch e
        e isa Test.TestSetException || rethrow()
        ts = nothing
        say("# suite: TestSetException: $(e.pass) passed, $(e.fail) failed, $(e.error) errors")
    end
    ok = ts !== nothing
    c = ok ? Test.get_test_counts(ts) : nothing
    npass = ok ? c.passes + c.cumulative_passes : -1
    nbad = ok ? c.fails + c.errors + c.cumulative_fails + c.cumulative_errors : -1
    emit("SUITE", "package suite at full size (L2_FULL_RUN), assertions passed (fail + error = $nbad); wall time $(round(t; digits = 1)) s", npass, "568 pass, 0 fail, 0 error", pf(ok && nbad == 0 && npass == 568))
end

# ---------------------------------------------------------------------------------------------------
# AC1 to AC5: analytic cases
section("AC1-5") do
    # AC1: linear map at the origin
    for T in (Float64, Double64)
        w = 0.0
        for tau in (1, 5, 10), n in (20, 100, 1000)
            sys = LinearSymplecticMap(SMatrix{2, 2, T}(2, 1, 1, 1))
            res = propagate_shadow(sys, [zero(SVector{2, T})], n; d0 = T(2.0^-20), tau = tau, seed = vseed(1), U0 = [SVector{2, T}(1, 0)])
            f = res.finals[1][1]
            Gfirst = Float64(l2_fibbig(2tau + 1)^2 + 2 * l2_fibbig(2tau)^2 + l2_fibbig(2tau - 1)^2)
            err = setprecision(BigFloat, 256) do
                Float64(abs(l2_big(f.sum + f.comp) - l2_ln_cat(n)))
            end
            w = max(w, err / l2_BS(T, n, Gfirst))
        end
        ratio_line("AC1", "max |n X2 - ln||A^n u0|| | / B_S(n), tau in {1,5,10}, n in {20,100,1000}, d0 = 2^-20, $T", w)
    end
    # AC2: K = 0 shear
    rng = StableRNG(vseed(2))
    angles = [2pi * rand(rng) for _ in 1:20]
    for T in (Float64, Double64)
        w = 0.0
        for d0 in (T(2.0^-20), T === Float64 ? 1e-9 : exact(Double64, "1e-9")), tau in (1, 10), n in (10, 100, 10_000)
            U0 = [SVector{2, T}(cos(th), sin(th)) for th in angles]
            res = propagate_shadow(StandardMap(T(0)), [zero(SVector{2, T}) for _ in angles], n; d0 = d0, tau = tau, seed = vseed(2), U0 = U0)
            errs = setprecision(BigFloat, 256) do
                [Float64(abs(l2_big(res.finals[j][1].sum + res.finals[j][1].comp) -
                             log(sqrt((l2_big(U0[j][1]) + n * l2_big(U0[j][2]))^2 + l2_big(U0[j][2])^2) / sqrt(l2_big(U0[j][1])^2 + l2_big(U0[j][2])^2)))) for j in eachindex(angles)]
            end
            w = max(w, maximum(errs) / l2_BS(T, n, 2 + tau^2))
        end
        ratio_line("AC2", "max |n X2 - ln sqrt((a+nb)^2+b^2)| / B_S(n), 20 seeded u0, d0 in {2^-20, 1e-9}, tau in {1,10}, n up to 1e4, $T", w)
    end
    # AC3: hyperbolic fixed point
    rng = StableRNG(vseed(3))
    angles = [2pi * rand(rng) for _ in 1:100]
    U0 = [SVector(cos(th), sin(th)) for th in angles]
    w = 0.0
    for K in (1.0, 2.0, 6.0), d0 in (1e-6, 1e-9), n in (20, 100, 1000)
        res = propagate_shadow(StandardMap(K), [zeros(SVector{2, Float64}) for _ in angles], n; d0 = d0, tau = 1, seed = vseed(3), U0 = U0)
        bound = l2_BS(Float64, n, (1 + K)^2 + 1 + K^2 + 1) + 2 * n * l2_Enl(K, d0)
        errs = setprecision(BigFloat, 256) do
            Kb = BigFloat(K); s = sqrt(Kb * (Kb + 4)); mu = (2 + Kb + s) / 2; vu = sqrt(1 + ((s - Kb) / 2)^2)
            [Float64(abs(l2_big(res.finals[j][1].sum + res.finals[j][1].comp) -
                         (n * log(mu) + log(abs(BigFloat(U0[j][1]) * (Kb + s) / 2 + BigFloat(U0[j][2])) / s * vu / sqrt(BigFloat(U0[j][1])^2 + BigFloat(U0[j][2])^2))))) for j in eachindex(angles)]
        end
        w = max(w, maximum(errs) / bound)
    end
    ratio_line("AC3", "max |n X2 - [n ln mu+ + ln(|c1| ||vu||)]| / (B_S + 2 n E_nl), K in {1,2,6}, d0 in {1e-6,1e-9}, 100 u0, n in {20,100,1000}", w)
    # AC4: minimal image
    d0 = 2.0^-20
    sys = WrappedStandardMap(0.0)
    Ero = l2_Ero(0.0, d0, L2_TWO_PI, EPS64)
    u0x = SVector(1.0, 0.0)
    strad = torus_state(SVector(-d0 / 2, 0.0))
    plain = torus_state(SVector(1.0, 0.0))
    wx, wd = 0.0, 0.0
    for n in (1, 10, 100)
        X2 = [f[1].X2 for f in propagate_shadow(sys, [strad, plain], n; d0 = d0, tau = 1, seed = vseed(4), U0 = [u0x, u0x]).finals]
        wx = max(wx, maximum(abs, X2) / (2 * Ero)); wd = max(wd, abs(X2[1] - X2[2]) / (2 * Ero))
    end
    ratio_line("AC4", "straddling pair and pair at x = 1: max |X2| / (2 E_ro(d0)), n in {1,10,100} (E_ro = $(fmt(Ero)))", wx)
    ratio_line("AC4", "same two pairs: |X2(straddle) - X2(x=1)| / (2 E_ro(d0))", wd)
    res = propagate_shadow(sys, [strad], 1; d0 = d0, tau = 1, seed = vseed(4), U0 = [u0x])
    xref, xsh = res.x[1], res.W[1]
    pl = SVector(xsh[1] - xref[1], xsh[2] - xref[2])
    sep = separation(sys, xref, xsh)
    lnart = log(sqrt(pl[1]^2 + pl[2]^2) / d0)
    emit("AC4", "plain wrapped difference of the straddling pair: ln(norm / d0), the artifact of T3", lnart, ">= 10", pf(lnart >= 10))
    r = abs(sqrt(sep[1]^2 + sep[2]^2) - d0) / (4 * Ero * d0)
    emit("AC4", "separation(): max |entry| <= pi ($(fmt(maximum(abs, sep)))) and ||sep|| = d0 within 4 E_ro d0, ratio", r, "<= 1", pf(r <= 1 && maximum(abs, sep) <= pi))
    ri = maximum(abs, separation(sys, xref, displaced(sys, xref, sep)) - sep) / (4 * EPS64 * L2_TWO_PI)
    ratio_line("AC4", "separation(displaced(sep)) - sep over 4 eps 2 pi", ri)
    double = l2_Plain(sys)
    err = l2_health_error(() -> propagate_shadow(double, [strad], 1; d0 = d0, tau = 1, seed = vseed(4), U0 = [u0x]))
    emit("AC4", "plain-difference double on the straddling pair: rejected before any work, kind/step/orbit/backend", err === nothing ? "no error" : string(err.what, "/", err.step, "/", err.orbit, "/", err.backend), "separation/0/1/cpu",
         pf(err isa NumericalHealthError && err.what == "separation" && err.step == 0 && err.orbit == 1 && err.backend == :cpu))
    ok = propagate_shadow(double, [plain], 10; d0 = d0, tau = 1, seed = vseed(4), U0 = [u0x]).finals[1][1].X2
    good = propagate_shadow(sys, [plain], 10; d0 = d0, tau = 1, seed = vseed(4), U0 = [u0x]).finals[1][1].X2
    ratio_line("AC4", "double vs correct hook at x = 1: |difference| / (2 E_ro)", abs(ok - good) / (2 * Ero))
    mirror = torus_state(SVector(1.0, -d0 / 2))
    wm = 0.0
    for n in (1, 10, 100)
        f = propagate_shadow(sys, [mirror], n; d0 = d0, tau = 1, seed = vseed(4), U0 = [SVector(0.0, 1.0)]).finals[1][1]
        err = setprecision(BigFloat, 256) do
            Float64(abs(l2_big(f.sum + f.comp) - log(sqrt(BigFloat(n)^2 + 1))))
        end
        wm = max(wm, err / (l2_BS(Float64, n, 3) + 2 * n * Ero))
    end
    ratio_line("AC4", "mirror case across p = 0, u0 = (0,1): max |n X2 - ln sqrt(n^2+1)| / (B_S + 2 n E_ro)", wm)
    # AC5: identical reference orbit, full grid
    for T in (Float64, Double64)
        bad, cells = 0, 0
        for mode in (:lifted, :wrapped), K in (0.5, L2_TWO_PI), tau in (1, 4)
            sys = l2_sys(mode, T(K))
            m = T === Float64 ? 100 : 20
            x0s = [l2_state(mode, T, z) for z in l2_starts(StableRNG(vseed(5)), m)]
            for N in (4, 8, 40, 400, 1000)
                plainx = map(x -> (z = x; for _ in 1:N; z = forward(sys, z); end; z), x0s)
                shadow = propagate_shadow(sys, x0s, N; d0 = T(1e-8), tau = tau, seed = vseed(5))
                l1 = propagate(sys, x0s, N; k = 1, tau = tau, seed = vseed(5), reducers = (Benettin(),))
                bad += !isequal(shadow.x, plainx) + !isequal(l1.x, plainx)
                cells += 1
            end
        end
        emit("AC5", "cells where the reference orbit differs bitwise from the plain forward loop or the L1 stream orbit ($cells cells: mode x K x tau x N, $T)", bad, "0", pf(bad == 0))
    end
end

# ---------------------------------------------------------------------------------------------------
# AC6: pathwise agreement with L1
section("AC6") do
    for (iK, K) in enumerate((0.5, 2.0, 5.0, L2_TWO_PI))
        starts = K == 0.5 ? l2_regular_starts(K, 50, StableRNG(vseed(6) + iK)) : l2_starts(StableRNG(vseed(6) + iK), 50)
        sys = WrappedStandardMap(K)
        x0s = [torus_state(z) for z in starts]
        Gf = [l2_Gfirst(sys, x, 1) for x in x0s]
        outside, worst, worst_old = 0, 0.0, 0.0
        perm = Dict{Int, Float64}()
        for n in (K == 0.5 ? (100, 1000, 10_000) : (100, 1000))
            X = l2_X(sys, x0s, n, vseed(6) + 10 + iK)
            Psi = l2_Psis(sys, x0s, vseed(6) + 10 + iK, Float64, n)
            for m in 3:12
                d0 = 10.0^(-m)
                if !l2_valid(K, d0, 1, L2_TWO_PI, EPS64)
                    outside += 1
                    continue
                end
                X2 = l2_X2(sys, x0s, n, vseed(6) + 10 + iK, d0)
                D = l2_D2T(K, d0, L2_TWO_PI, EPS64)
                bound = [l2_pathbound(D, Psi[j], n, l2_BS(Float64, n, Gf[j], K)) for j in eachindex(x0s)]
                r = maximum(abs.(X2 .- X) ./ bound)
                worst = max(worst, r); worst_old = max(worst_old, maximum(abs.(X2 .- X)) / D)
                perm[m] = max(get(perm, m, 0.0), r)
            end
        end
        ratio_line("AC6", "K = $K: max |X2 - X| / (Delta_2T Psi_n/n + B_S/n) (T5a), 50 starts, d0 = 10^-m m = 3..12, n up to $(K == 0.5 ? "1e4" : "1e3"); points outside validity (reported) = $outside", worst)
        say("#   K = $K: old ratio max |X2 - X| / Delta_2T = $(fmt(worst_old)); T5a ratio per decade m: " * join(["$m:$(fmt(perm[m]))" for m in sort(collect(keys(perm)))], " "))
    end
end

# ---------------------------------------------------------------------------------------------------
# AC7 and AC8: U-shape and Double64
function ushape_data()
    return get!(DATA, :ushape) do
        seed = vseed(7)
        sys = WrappedStandardMap(5.0)
        x0s = [l2_state(:wrapped, Float64, z) for z in l2_starts(StableRNG(seed), 200)]
        X = l2_X(sys, x0s, 1000, seed + 1)
        e, emax = Dict{Int, Float64}(), Dict{Int, Float64}()
        for m in 2:13
            err = [abs(a - b) for (a, b) in zip(l2_X2(sys, x0s, 1000, seed + 1, 10.0^(-m)), X)]
            e[m] = l2_median(err); emax[m] = maximum(err)
        end
        (e = e, emax = emax)
    end
end

function dd_data()
    return get!(DATA, :dd) do
        seed = vseed(8)
        meds, worst = Dict{Int, Float64}(), Dict{Float64, Float64}()
        for (iK, K) in enumerate((5.0, L2_TWO_PI))
            sys = WrappedStandardMap(Double64(K))
            x0s = [l2_state(:wrapped, Double64, z) for z in l2_starts(StableRNG(seed + iK), 50)]
            n = 1000
            X = l2_X(sys, x0s, n, seed + 10 + iK)
            Gf = [l2_Gfirst(sys, x, 1) for x in x0s]
            Psi = l2_Psis(sys, x0s, seed + 10 + iK, Double64, n)
            w = 0.0
            for m in 3:28
                d0 = exact(Double64, "1e-$m")
                l2_valid(K, Float64(d0), 1, L2_TWO_PI, EPSDD) || continue
                X2 = l2_X2(sys, x0s, n, seed + 10 + iK, d0)
                D = l2_D2T(K, Float64(d0), L2_TWO_PI, EPSDD)
                w = max(w, maximum(j -> Float64(abs(X2[j] - X[j])) / l2_pathbound(D, Psi[j], n, l2_BS(Double64, n, Gf[j], K)), eachindex(x0s)))
                K == 5.0 && (meds[m] = l2_median([Float64(abs(X2[j] - X[j])) for j in eachindex(x0s)]))
            end
            worst[K] = w
        end
        (meds = meds, worst = worst)
    end
end

section("AC7-8") do
    u = ushape_data()
    e = u.e
    slope(ms) = l2_slope([-Float64(m) for m in ms], [log10(e[m]) for m in ms])
    sa, sb = slope(2:4), slope(10:13)
    emit("AC7", "slope of log10 median|X2-X| against log10 d0, m in {2,3,4} (E_nl linear in d0), K = 5, 200 starts, n = 1e3", sa, "[0.8, 2.2]", pf(0.8 <= sa <= 2.2))
    emit("AC7", "same, m in {10,...,13} (E_ro proportional to 1/d0)", sb, "[-1.3, -0.7]", pf(-1.3 <= sb <= -0.7))
    mstar = argmin(e)
    d0star = l2_d0star(5.0, L2_TWO_PI, EPS64)
    emit("AC7", "grid minimum at m = $mstar: d0/d0* with d0* = $(fmt(d0star)) (T8)", 10.0^(-mstar) / d0star, "[0.01, 100]", pf(d0star / 100 <= 10.0^(-mstar) <= 100 * d0star))
    D = l2_D2T(5.0, d0star, L2_TWO_PI, EPS64)
    ratio_line("AC7", "grid minimum of the median error $(fmt(e[mstar])) over Delta_2T(d0*) = $(fmt(D))", e[mstar] / D)
    say("#   AC7 median |X2 - X| per m: " * join(["$m:$(fmt(e[m]))" for m in 2:13], " "))
    d = dd_data()
    for K in (5.0, L2_TWO_PI)
        ratio_line("AC8", "Double64 K = $(fmt(K)): max |X2 - X| / (Delta_2T Psi_n/n + B_S/n) (T5a, eps 2^-100), 50 starts, n = 1e3, m = 3..28", d.worst[K])
    end
    emin, efmin = minimum(values(d.meds)), minimum(values(e))
    emit("AC8", "Double64 K = 5 grid minimum of the median error", emin, "<= 1e-9", pf(emin <= 1e-9))
    emit("AC8", "Float64 grid minimum / Double64 grid minimum (T8 predicts 1.7e7)", efmin / emin, ">= 1e3", pf(efmin / emin >= 1e3))
    say("#   AC8 Double64 median |X2 - X| per m (K = 5): " * join(["$m:$(fmt(d.meds[m]))" for m in sort(collect(keys(d.meds)))], " "))
end

# ---------------------------------------------------------------------------------------------------
# AC9: tau sweep
function tau_data()
    return get!(DATA, :tau) do
        seed = vseed(9)
        n, d0 = 1000, 1e-8
        out = Dict{Tuple{Float64, Int}, NamedTuple}()
        for (K, taus) in ((0.5, (1, 2, 4, 8)), (5.0, (1, 2, 4, 8)))
            starts = K == 0.5 ? l2_regular_starts(K, 50, StableRNG(seed)) : l2_starts(StableRNG(seed + 1), 50)
            sys = WrappedStandardMap(K)
            x0s = [torus_state(z) for z in starts]
            X = l2_X(sys, x0s, n, seed)
            Gf = [l2_Gfirst(sys, x, 1) for x in x0s]
            Psi = l2_Psis(sys, x0s, seed, Float64, n)
            for tau in taus
                n % tau == 0 || continue
                X2 = l2_X2(sys, x0s, n, seed, d0; tau = tau)
                D = l2_D2T(K, d0, L2_TWO_PI, EPS64)
                bound = [l2_pathbound(D, Psi[j], n, l2_BS(Float64, n, Gf[j], K); factor = l2_sigma(K)^(2 * (tau - 1))) for j in eachindex(x0s)]
                out[(K, tau)] = (ratio = maximum(abs.(X2 .- X) ./ bound), maxerr = maximum(abs.(X2 .- X)), mederr = l2_median(abs.(X2 .- X)))
            end
        end
        out
    end
end

section("AC9") do
    t = tau_data()
    for (K, tau) in ((0.5, 1), (0.5, 2), (0.5, 4), (5.0, 1), (5.0, 2))
        r = t[(K, tau)]
        ratio_line("AC9", "K = $K, tau = $tau: max |X2 - X| / (sigma^(2(tau-1)) Delta_2T Psi_n/n + B_S/n), d0 = 1e-8, 50 starts, n = 1e3 (gated)", r.ratio)
    end
    for (K, tau) in ((5.0, 4), (5.0, 8))
        haskey(t, (K, tau)) || continue
        emit("AC9", "K = $K, tau = $tau: same ratio (reported, not gated)", t[(K, tau)].ratio, "reported", "INFO")
    end
    sys = WrappedStandardMap(10.0)
    x0s = [torus_state(z) for z in l2_starts(StableRNG(vseed(9) + 2), 50)]
    N, d0c = 1008, 1e-3
    X2 = l2_X2(sys, x0s, N, vseed(9), d0c; tau = 16)
    X = l2_X(sys, x0s, N, vseed(9))
    cap = log(pi * sqrt(2) / d0c) / 16
    emit("AC9", "negative control (T7), K = 10, d0 = 1e-3, tau = 16: max X2 against ln(pi sqrt 2/d0)/16", maximum(X2), "<= $(fmt(cap))", pf(all(X2 .<= cap + 1e-12)))
    idx = findall(>=(1.1), X)
    gap = isempty(idx) ? NaN : minimum(X[idx] .- X2[idx])
    emit("AC9", "negative control: min (X - X2) over the $(length(idx)) of 50 starts with X >= 1.1", gap, ">= 0.5", pf(!isempty(idx) && gap >= 0.5))
end

# ---------------------------------------------------------------------------------------------------
# AC10 and AC11
section("AC10-11") do
    seed = vseed(10)
    d0 = 1e-8
    errs = []
    push!(errs, l2_health_error(() -> propagate_shadow(l2_Boom(), [SVector(1.0, 2.0), SVector(1.0, 2.0)], 4; d0 = 1e-20, tau = 1, seed = seed)))
    push!(errs, l2_health_error(() -> propagate_shadow_orbit(l2_Boom(), SVector(1.0, 2.0), SVector(1.0, 0.0), 4, 1, 1e-20; orbit = 5)))
    oka = all(e -> e isa NumericalHealthError && e.what == "separation" && e.step == 0 && e.backend === :cpu && e.value == 0.0, errs)
    emit("AC10", "(a) d0 = 1e-20 at magnitude ~1: rejected before any work (the map throws if called), kind separation, step 0, backend cpu, both entry points", oka, "true", pf(oka))
    x0 = [SVector(0.5, i == 3 ? 1.0 : 0.0) for i in 1:5]
    U0 = fill(SVector(0.6, 0.8), 5)
    sysc = l2_Fault(:collapse, 2.25)
    e1 = l2_health_error(() -> propagate_shadow(sysc, x0, 4; d0 = d0, tau = 1, seed = seed, U0 = U0, chunk = 2))
    e2 = l2_health_error(() -> propagate_shadow_orbit(sysc, x0[3], U0[3], 4, 1, d0; orbit = 3))
    okb = all(e -> e isa NumericalHealthError && e.what == "separation" && e.orbit == 3 && e.step == 2, (e1, e2)) && e1.backend === :cpu
    emit("AC10", "(b) map collapsing the pair at the second handoff, orbit 3: error with kind separation, orbit 3, step 2, nothing stored", okb, "true", pf(okb))
    okc = true
    for kind in (:inf_ref, :inf_shadow, :inf_both)
        sysi = l2_Fault(kind, 2.25)
        for call in (() -> propagate_shadow(sysi, x0, 4; d0 = d0, tau = 1, seed = seed, U0 = U0, chunk = 2), () -> propagate_shadow_orbit(sysi, x0[3], U0[3], 4, 1, d0; orbit = 3))
            e = l2_health_error(call)
            okc &= e isa NumericalHealthError && e.what == "state" && e.orbit == 3 && e.step == 2 && e.threshold == Inf && !isfinite(e.value)
        end
    end
    emit("AC10", "(c) non-finite in reference, shadow or both: kind state, orbit 3, step 2, threshold Inf (6 calls)", okc, "true", pf(okc))
    okd = true
    for (T, mode) in ((Float64, :lifted), (Float64, :wrapped), (Double64, :wrapped))
        sys = l2_sys(mode, T(5))
        x0s = [l2_state(mode, T, z) for z in l2_starts(StableRNG(seed), 8)]
        dT = T === Float64 ? d0 : exact(Double64, "1e-8")
        p = propagate_shadow(sys, x0s, 20; d0 = dT, tau = 2, seed = seed).provenance
        okd &= p.estimator == "two-trajectory" && p.state_mode == string(mode) && p.k == 0 && p.d0 == Float64(dT) && p.tau == 2 && p.seed == seed && p.precision == string(T)
    end
    emit("AC10", "(d) run record carries d0, tau, estimator, state_mode, seed, precision (3 type/mode cases)", okd, "true", pf(okd))
    sys = StandardMap(2.0)
    x1 = [SVector(0.3, 0.4)]
    thr(f) = try; f(); false; catch e; e isa ArgumentError; end
    okx = thr(() -> propagate_shadow(sys, x1, 7; d0 = d0, tau = 5, seed = seed)) && thr(() -> propagate_shadow(sys, x1, 0; d0 = d0, tau = 1, seed = seed)) &&
          all(bad -> thr(() -> propagate_shadow(sys, x1, 4; d0 = bad, tau = 1, seed = seed)), (0.0, -1e-6, NaN, Inf)) &&
          setprecision(BigFloat, 256) do
              thr(() -> propagate_shadow(sys, x1, 4; d0 = BigFloat("1e-8"), tau = 1, seed = seed))
          end
    emit("AC10", "contract checks: N not a multiple of tau, N = 0, d0 in {0, <0, NaN, Inf}, BigFloat d0 narrowing: ArgumentError", okx, "true", pf(okx))

    # AC11
    m, n, d0 = 4096, 1000, 1e-8
    sysw = WrappedStandardMap(L2_TWO_PI)
    x0w = [l2_state(:wrapped, Float64, z) for z in l2_starts(StableRNG(vseed(11)), m)]
    u0s = [SVector{2, Float64}(initial_deviations(vseed(11), j, Val(2), Val(1), Float64)[:, 1]) for j in 1:m]
    serial = [propagate_shadow_orbit(sysw, x0w[j], u0s[j], n, 1, d0; orbit = j) for j in 1:m]
    for chunk in (64, 1000, nothing)
        res = propagate_shadow(sysw, x0w, n; d0 = d0, tau = 1, seed = vseed(11), chunk = chunk)
        ok = isequal(res.x, [s[1] for s in serial]) && isequal(res.W, [s[2] for s in serial]) && isequal(res.finals, [(s[3],) for s in serial])
        emit("AC11", "propagate_shadow equals the serial map of propagate_shadow_orbit bitwise, 4096 orbits, K = 2 pi, n = 1e3, chunk = $(chunk === nothing ? "default" : chunk)", ok, "true", pf(ok))
    end
    hashes = Dict{Int, String}(Threads.nthreads() => string(l2_thread_hash()))
    for nt in (1, 4)
        nt == Threads.nthreads() && continue
        cmd = addenv(`$(Base.julia_cmd()) -t $nt --project=$(dirname(Base.active_project())) $(@__FILE__)`, "L2_CHILD" => "1")
        txt = read(pipeline(cmd; stderr = devnull), String)
        mt = match(r"HASH (\d+) THREADS (\d+)", txt)
        hashes[nt] = mt === nothing ? "child failed: $(first(txt, 100))" : mt.captures[1]
    end
    okt = length(unique(values(hashes))) == 1
    emit("AC11", "X2 series, sums, x and W bitwise identical for thread counts $(sort(collect(keys(hashes)))) (hashes of all outputs, separate processes)", okt, "true", pf(okt))
    for T in (Float64, Double64), mode in (:wrapped, :lifted)
        dT = T === Float64 ? d0 : exact(Double64, "1e-8")
        s = l2_sys(mode, T(L2_TWO_PI))
        a = l2_alloc(s, l2_state(mode, T, SVector(1.0, 2.0)), SVector{2, T}(T(0.6), T(0.8)), 10_000, 1, dT)
        emit("AC11", "bytes allocated by propagate_shadow_orbit after warm-up, n = 1e4, $T, $mode", a, "0", pf(a == 0))
    end
end

# ---------------------------------------------------------------------------------------------------
# AC12: wrapped versus lifted
section("AC12") do
    seed = vseed(12)
    K, d0 = 5.0, 1e-8
    maxerr = Dict{Symbol, Float64}(:wrapped => 0.0, :lifted => 0.0)
    for (n, nstart) in ((1000, 100), (100_000, 100))
        starts = l2_starts(StableRNG(seed), nstart)
        for mode in (:wrapped, :lifted)
            sys = l2_sys(mode, K)
            x0s = [l2_state(mode, Float64, z) for z in starts]
            X = l2_X(sys, x0s, n, seed)
            w, tested, reported, errored = 0.0, 0, 0, 0
            for j in eachindex(x0s)
                M = mode === :wrapped ? L2_TWO_PI : l2_Mstate(sys, x0s[j], n)
                Ero = l2_Ero(K, d0, M, EPS64)
                try
                    u = SVector{2, Float64}(initial_deviations(seed, j, Val(2), Val(1), Float64)[:, 1])
                    X2 = propagate_shadow(sys, x0s[j:j], n; d0 = d0, tau = 1, seed = seed, U0 = [u]).finals[1][1].X2
                    if Ero <= 0.5
                        tested += 1
                        Psi = lb_Psi(sys, x0s[j], u, n)
                        w = max(w, abs(X2 - X[j]) / l2_pathbound(l2_D2T(K, d0, M, EPS64), Psi, n, l2_BS(Float64, n, l2_Gfirst(sys, x0s[j], 1), K)))
                        n == 100_000 && (maxerr[mode] = max(maxerr[mode], abs(X2 - X[j])))
                    else
                        reported += 1
                    end
                catch e
                    e isa NumericalHealthError && e.what == "separation" || rethrow()
                    Ero > 0.5 ? (reported += 1) : (errored += 1)
                end
            end
            ratio_line("AC12", "$mode, n = $(fmt(Float64(n))): max |X2 - X| / (Delta_2T Psi_n/n + B_S/n) over $tested tested starts of $nstart (E_ro > 0.5, reported: $reported; unexpected separation errors: $errored), K = 5, d0 = 1e-8", w)
            (n == 1000 && mode === :lifted) && emit("AC12", "lifted n = 1e3: all 100 starts valid and tested, no unexpected errors", (tested, errored), "(100, 0)", pf(tested == nstart && errored == 0))
        end
    end
    emit("AC12", "recorded (no threshold): max |X2 - X| at n = 1e5 over valid starts, lifted / wrapped = $(fmt(maxerr[:lifted])) / $(fmt(maxerr[:wrapped]))", maxerr[:lifted] / maxerr[:wrapped], "recorded", "INFO")
end

# ---------------------------------------------------------------------------------------------------
# Oracle (D30): ChaosTools lyapunov (two-trajectory, Benettin 1976 per its docstring) on the standard map
# Convention (read from ChaosTools 3.5.4 src/chaosdetection/lyapunovs/lyapunov.jl, see the report):
#   state (theta, p), parameter k: p += k sin(theta); theta += p. Lifted rule: no wrap. Wrapped rule: mod2pi after every step (docs rule).
#   lyapunov(ds, T; u0, d0, inittest, d0_lower = 1e-3 d0, d0_upper = 1e3 d0, Dt = 1, Ttr = 0): test state u1 + inittest; test orbit
#   rescaled to distance d0 along the separation only when the plain Euclidean distance leaves [d0_lower, d0_upper]; lambda = sum ln(d/d0) / T.
ct_lifted(x, par, n) = (p = x[2] + par[1] * sin(x[1]); SVector(x[1] + p, p))
ct_wrapped(x, par, n) = mod2pi.(ct_lifted(x, par, n))

function ct_lyap(rule, K, x0, udir, n, d0; kw...)
    ds = CT.DeterministicIteratedMap(rule, x0, [K])
    return CT.lyapunov(ds, n; u0 = x0, d0 = d0, inittest = (u1, d) -> u1 .+ d .* udir, kw...)
end

seeded_dirs(seed, m) = [SVector{2, Float64}(initial_deviations(seed, j, Val(2), Val(1), Float64)[:, 1]) for j in 1:m]
mean_se(v) = (mean(v), std(v) / sqrt(length(v)))

# Oracle and Chaoslab estimates on the same starts, lifted mode, same d0, oracle test orbit placed along the seeded direction of Chaoslab
function oracle_pair(K, starts, n, d0, seed; kw...)
    sys = StandardMap(K)
    x0s = [SVector{2, Float64}(z) for z in starts]
    X2 = l2_X2(sys, x0s, n, seed, d0)
    dirs = seeded_dirs(seed, length(x0s))
    CTv = [ct_lyap(ct_lifted, K, x0s[j], dirs[j], n, d0; kw...) for j in eachindex(x0s)]
    return X2, CTv
end

function ensemble_line(ac, label, X2, CTv)
    m1, s1 = mean_se(X2); m2, s2 = mean_se(CTv)
    se = sqrt(s1^2 + s2^2)
    z = abs(m1 - m2) / se
    emit(ac, "$label: mean X2 $(fmt(m1)) (SE $(fmt(s1))), mean oracle $(fmt(m2)) (SE $(fmt(s2))), |diff| = $(fmt(abs(m1 - m2))); |diff| / combined SE", z, "<= 4", pf(z <= 4))
    say("#   pathwise |X2 - oracle|: median $(fmt(l2_median(abs.(X2 .- CTv)))), max $(fmt(maximum(abs.(X2 .- CTv))))")
    return (m1, s1, m2, s2)
end

section("oracle") do
    d0 = 1e-6
    # convention check: both sides step the same lifted orbit bitwise
    z0 = SVector(1.25, 2.5)
    ds = CT.DeterministicIteratedMap(ct_lifted, z0, [5.0])
    zc = z0; sysl = StandardMap(5.0); zo = z0
    for _ in 1:1000
        CT.step!(ds); zo = forward(sysl, zo)
    end
    emit("ORACLE", "orbit of ChaosTools lifted rule vs Chaoslab forward, K = 5, 1000 steps: max |difference| (same arithmetic)", maximum(abs, SVector{2, Float64}(CT.current_state(ds)) - zo), "0 (information)", "INFO")

    # AC13: regular class, pathwise
    starts = l2_regular_starts(0.5, 20, StableRNG(vseed(13)))
    X2, CTv = oracle_pair(0.5, starts, 100_000, d0, vseed(13))
    dif = maximum(abs.(X2 .- CTv))
    emit("AC13", "K = 0.5 regular class, 20 starts, n = 1e5, d0 = 1e-6, lifted: max |X2 - lambda_oracle| (mean X2 $(fmt(mean(X2))), mean oracle $(fmt(mean(CTv))), Delta_2T(1e-6) = $(fmt(l2_D2T(0.5, 1e-6, L2_TWO_PI, EPS64))))", dif, "<= 1e-4", pf(dif <= 1e-4))
    # AC13 diagnostics (not gated): which side departs from the L1 tangent value X, and the oracle with narrower rescaling thresholds
    X13 = l2_X(StandardMap(0.5), [SVector{2, Float64}(z) for z in starts], 100_000, vseed(13))
    emit("AC13", "info: max |X2 - X_L1| and max |oracle - X_L1| (same orbit and direction, lifted), default thresholds", string(fmt(maximum(abs.(X2 .- X13))), " / ", fmt(maximum(abs.(CTv .- X13)))), "info", "INFO")
    sys13 = StandardMap(0.5)
    x13 = [SVector{2, Float64}(z) for z in starts]
    Ms = [l2_Mstate(sys13, x13[j], 100_000) for j in eachindex(x13)]
    Psi13 = l2_Psis(sys13, x13, vseed(13), Float64, 100_000)
    r13 = maximum(j -> abs(X2[j] - X13[j]) / l2_pathbound(l2_D2T(0.5, d0, Ms[j], EPS64), Psi13[j], 100_000, l2_BS(Float64, 100_000, l2_Gfirst(sys13, x13[j], 1), 0.5)), eachindex(x13))
    emit("AC13", "info: lifted-mode T5a ratio |X2 - X_L1| / (Delta_2T(M_state) Psi_n/n + B_S/n) with the measured M_state (max over starts $(fmt(maximum(Ms))), max Psi_n/n $(fmt(maximum(Psi13) / 100_000)), E_ro(max M) = $(fmt(l2_Ero(0.5, d0, maximum(Ms), EPS64))))", r13, "info (vacuous bound)", "INFO")
    for d0x in (1e-5, 1e-4)
        X2x, CTx = oracle_pair(0.5, starts, 100_000, d0x, vseed(13))
        emit("AC13", "info: same comparison at d0 = $d0x (default thresholds): max |X2 - oracle| / max |X2 - X_L1|", string(fmt(maximum(abs.(X2x .- CTx))), " / ", fmt(maximum(abs.(X2x .- X13)))), "info", "INFO")
    end
    for (lo, up) in ((d0 / 2, 2d0), (d0 / 10, 10d0), (d0 / 100, 100d0))
        _, CTn = oracle_pair(0.5, starts, 100_000, d0, vseed(13); d0_lower = lo, d0_upper = up)
        emit("AC13", "info: oracle thresholds [d0*$(fmt(lo / d0)), d0*$(fmt(up / d0))]: max |X2 - oracle| / max |oracle - X_L1|", string(fmt(maximum(abs.(X2 .- CTn))), " / ", fmt(maximum(abs.(CTn .- X13)))), "info (gate 1e-4 applies to the default)", "INFO")
    end
    # AC14
    rng = StableRNG(vseed(14))
    starts = l2_starts(rng, 1000)
    X2, CTv = oracle_pair(L2_TWO_PI, starts, 10_000, d0, vseed(14))
    ensemble_line("AC14", "K = 2 pi, 1000 starts, n = 1e4, d0 = 1e-6", X2, CTv)
    # informational: oracle with a tight schedule (rescaled when leaving [d0/2, 2 d0], close to the periodic schedule)
    _, CTt = oracle_pair(L2_TWO_PI, starts[1:200], 10_000, d0, vseed(14); d0_lower = d0 / 2, d0_upper = 2d0)
    m1, s1 = mean_se(X2[1:200]); m2, s2 = mean_se(CTt)
    emit("AC14", "info: first 200 starts, oracle with tight thresholds [d0/2, 2 d0]: |diff of means| / combined SE (X2 mean $(fmt(m1)), oracle $(fmt(m2)))", abs(m1 - m2) / sqrt(s1^2 + s2^2), "info", "INFO")
    # AC15
    for K in (5.0, 10.0)
        starts = l2_starts(StableRNG(vseed(15) + Int(K)), 200)
        X2, CTv = oracle_pair(K, starts, 10_000, d0, vseed(15) + Int(K))
        ensemble_line("AC15", "K = $(Int(K)), 200 starts, n = 1e4, d0 = 1e-6", X2, CTv)
    end
    # AC16: cross-estimator, wrapped, d0 = d0*
    K = L2_TWO_PI
    d0s = l2_d0star(K, L2_TWO_PI, EPS64)
    sysw = WrappedStandardMap(K)
    x0s = [torus_state(z) for z in l2_starts(StableRNG(vseed(16)), 1000)]
    X = l2_X(sysw, x0s, 10_000, vseed(16))
    X2w = l2_X2(sysw, x0s, 10_000, vseed(16), d0s)
    ensemble_line("AC16", "K = 2 pi wrapped, 1000 starts, n = 1e4, d0 = d0* = $(fmt(d0s)): X2 (first) against L1 X (second, in the 'oracle' slot)", X2w, X)
    say("#   AC16 pathwise |X2 - X| (same orbit and direction): median $(fmt(l2_median(abs.(X2w .- X)))), max $(fmt(maximum(abs.(X2w .- X))))")
    # AC17 stays with the GPU environment
    emit("AC17", "CPU vs GPU on X2: needs a CUDA device and the shadow kernel; verify/gpu.jl has no L2 case", "not run here", "tangent-core AC12 tolerances", "STAGE04-GPU")
    # R2: wrapped (mod 2 pi) oracle rule, lifted oracle rule, K = 0.5 and 5, d0 = 1e-6, n = 1e5
    for K in (0.5, 5.0)
        starts = K == 0.5 ? l2_regular_starts(0.5, 40, StableRNG(vseed(17))) : l2_starts(StableRNG(vseed(17) + 1), 40)
        x0s = [SVector{2, Float64}(z) for z in starts]
        dirs = seeded_dirs(vseed(17), length(x0s))
        lif = [ct_lyap(ct_lifted, K, x0s[j], dirs[j], 100_000, 1e-6) for j in eachindex(x0s)]
        wrp = [ct_lyap(ct_wrapped, K, x0s[j], dirs[j], 100_000, 1e-6) for j in eachindex(x0s)]
        ml, sl = mean_se(lif); mw, sw = mean_se(wrp)
        emit("R2", "oracle with wrapped rule vs lifted rule, K = $K, 40 starts, n = 1e5, d0 = 1e-6: mean wrapped $(fmt(mw)) (SE $(fmt(sw))), mean lifted $(fmt(ml)) (SE $(fmt(sl))), artifact = wrapped - lifted", mw - ml, "recorded", "INFO")
    end
end

# ---------------------------------------------------------------------------------------------------
# Figures
section("figure") do
    u = ushape_data(); d = dd_data(); t = tau_data()
    fig = CM.Figure(size = (1200, 480))
    ax = CM.Axis(fig[1, 1]; xscale = log10, yscale = log10, xlabel = "d0", ylabel = "median |X2 - X| over starts", title = "AC7, AC8: error against d0, K = 5, n = 1e3, wrapped")
    ms = 2:13
    CM.scatterlines!(ax, [10.0^(-m) for m in ms], [u.e[m] for m in ms]; color = :crimson, label = "Float64, 200 starts, median")
    CM.scatterlines!(ax, [10.0^(-m) for m in ms], [u.emax[m] for m in ms]; color = :crimson, linestyle = :dash, label = "Float64, max")
    md = sort(collect(keys(d.meds)))
    CM.scatterlines!(ax, [10.0^(-m) for m in md], [max(d.meds[m], 1e-20) for m in md]; color = :navy, label = "Double64, 50 starts, median")
    dgrid = 10.0 .^ range(-28, -2; length = 200)
    CM.lines!(ax, dgrid, [l2_D2T(5.0, x, L2_TWO_PI, EPS64) for x in dgrid]; color = :black, label = "Delta_2T, Float64 (T5, per-step)")
    CM.lines!(ax, dgrid, [l2_D2T(5.0, x, L2_TWO_PI, EPSDD) for x in dgrid]; color = :gray40, linestyle = :dot, label = "Delta_2T, Double64")
    CM.vlines!(ax, [l2_d0star(5.0, L2_TWO_PI, EPS64)]; color = :steelblue, linestyle = :dashdot, label = "d0* Float64 (T8)")
    CM.ylims!(ax, 1e-20, 1e2)
    CM.axislegend(ax; position = :ct, labelsize = 10)
    ax2 = CM.Axis(fig[1, 2]; yscale = log10, xlabel = "tau", ylabel = "max |X2 - X| over 50 starts", title = "AC9: tau sweep, d0 = 1e-8, n = 1e3", xticks = [1, 2, 4, 8])
    for (K, c) in ((0.5, :darkorange), (5.0, :seagreen))
        taus = sort([tau for (kk, tau) in keys(t) if kk == K])
        CM.scatterlines!(ax2, taus, [t[(K, tau)].maxerr for tau in taus]; color = c, label = "K = $K")
        CM.lines!(ax2, taus, [l2_sigma(K)^(2 * (tau - 1)) * l2_D2T(K, 1e-8, L2_TWO_PI, EPS64) for tau in taus]; color = c, linestyle = :dash, label = "K = $K, sigma^(2(tau-1)) Delta_2T (no Psi)")
    end
    CM.axislegend(ax2; position = :lt, labelsize = 10)
    CM.save(joinpath(OUT, "lyap-two-trajectory-ushape.png"), fig)
    say("# figure saved: lyap-two-trajectory-ushape.png")
end

section("figure-oracle") do
    d0 = 1e-6
    fig = CM.Figure(size = (1200, 480))
    ax = CM.Axis(fig[1, 1]; xlabel = "oracle lambda (ChaosTools, lifted)", ylabel = "Chaoslab X2", title = "AC13: K = 0.5 regular, 20 starts, n = 1e5, d0 = 1e-6")
    starts = l2_regular_starts(0.5, 20, StableRNG(vseed(13)))
    X2, CTv = oracle_pair(0.5, starts, 100_000, d0, vseed(13))
    CM.scatter!(ax, CTv, X2; color = :crimson)
    lim = extrema(vcat(CTv, X2))
    CM.lines!(ax, [lim...], [lim...]; color = :black)
    ax2 = CM.Axis(fig[1, 2]; xticks = (1:4, ["2 pi", "5", "10", "2 pi (L1)"]), ylabel = "ensemble mean of exponent estimate", title = "AC14 to AC16: ensemble means, n = 1e4, +- 1 SE")
    specs = ((L2_TWO_PI, 1000, vseed(14), 1), (5.0, 200, vseed(15) + 5, 2), (10.0, 200, vseed(15) + 10, 3))
    for (K, m, sd, pos) in specs
        starts = K == L2_TWO_PI ? l2_starts(StableRNG(vseed(14)), 1000) : l2_starts(StableRNG(sd), m)
        a, b = oracle_pair(K, starts, 10_000, d0, sd)
        m1, s1 = mean_se(a); m2, s2 = mean_se(b)
        CM.errorbars!(ax2, [pos - 0.12], [m1], [s1]; color = :crimson); CM.scatter!(ax2, [pos - 0.12], [m1]; color = :crimson, label = pos == 1 ? "Chaoslab X2" : nothing)
        CM.errorbars!(ax2, [pos + 0.12], [m2], [s2]; color = :navy); CM.scatter!(ax2, [pos + 0.12], [m2]; color = :navy, label = pos == 1 ? "ChaosTools lyapunov" : nothing)
    end
    sysw = WrappedStandardMap(L2_TWO_PI)
    x0s = [torus_state(z) for z in l2_starts(StableRNG(vseed(16)), 1000)]
    X = l2_X(sysw, x0s, 10_000, vseed(16))
    X2w = l2_X2(sysw, x0s, 10_000, vseed(16), l2_d0star(L2_TWO_PI, L2_TWO_PI, EPS64))
    m1, s1 = mean_se(X2w); m2, s2 = mean_se(X)
    CM.errorbars!(ax2, [4 - 0.12], [m1], [s1]; color = :crimson); CM.scatter!(ax2, [4 - 0.12], [m1]; color = :crimson)
    CM.errorbars!(ax2, [4 + 0.12], [m2], [s2]; color = :darkorange); CM.scatter!(ax2, [4 + 0.12], [m2]; color = :darkorange, label = "Chaoslab L1 X")
    CM.axislegend(ax2; position = :lt, labelsize = 10)
    CM.save(joinpath(OUT, "lyap-two-trajectory-oracle.png"), fig)
    say("# figure saved: lyap-two-trajectory-oracle.png")
end

say(@sprintf("# verdicts: %d PASS, %d FAIL, %d other", count(==("PASS"), VERDICTS), count(==("FAIL"), VERDICTS), count(v -> v != "PASS" && v != "FAIL", VERDICTS)))
say(@sprintf("# total wall time %.1f s", time() - T0))
close(LOG)
