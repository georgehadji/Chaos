# Stage 04 verification of module gali (GALI_k reducer, src/Indicators/gali.jl; catalog I3).
# Run from the repo root, one process per group (the groups are independent; GA_LOG names the log, GA_APPEND=1 appends):
#   GA_ONLY=AC1-4   julia --project=verify -t 8 verify/gali.jl      (also AC5-8, AC9-12, AC13-16, AC17-20 = the package test files at the FULL sizes)
#   GA_ONLY=AC19    julia --project=verify -t 8 verify/gali.jl      (ChaosTools gali, criterion 19)
#   GA_ONLY=AC20t   julia --project=verify -t 8 verify/gali.jl      (scale timing, criterion 20 V)
#   GA_ONLY=R       julia --project=verify -t 8 verify/gali.jl      (recorded lines R1, R2, R4, R5 and the AC2 / AC9 realised ratios)
#   GA_ONLY=figure  julia --project=verify -t 8 verify/gali.jl      (gali-laws.png)
# Spec: workspaces/module-build/stages/01-spec/output/gali-spec.md (criteria 1 to 21, Amendment D107). Seeds: 20261012 + 40000..40999 (ga_seed), the verify-only
# sections use 20261012 + 40900 .. 40999.
# One printed line per measurement:  AC<n> | quantity | value | bound | PASS/FAIL/INFO
# The package test files (test/gali/) are included with GA_FULL_RUN = true: the sizes marked D38 there take their full (spec) values, the assertions and tolerances
# are the package ones. Each testset is reported as one line with its counts; the @info lines of the tests (realised ratios, slopes, counts) are copied into the log.
const ONLY0 = filter(!isempty, split(get(ENV, "GA_ONLY", ""), ","))
using Printf, Random, LinearAlgebra, Statistics, Serialization, Logging
using Test
const GA_FULL_RUN = true
using Chaoslab, StaticArrays, DoubleFloats, StableRNGs, BenchmarkTools
import ChaosTools as CT
if isempty(ONLY0) || "figure" in ONLY0
    import CairoMakie as CM
end

const ROOT = dirname(@__DIR__)
const OUT = joinpath(ROOT, "workspaces", "module-build", "stages", "04-verify", "output")
mkpath(OUT)
const ONLY = filter(!isempty, split(get(ENV, "GA_ONLY", ""), ","))
const LOG = open(get(ENV, "GA_LOG", joinpath(OUT, "gali-verify-log.txt")), get(ENV, "GA_APPEND", "") == "1" ? "a" : "w")

fmt(x::AbstractFloat) = @sprintf("%.4g", x)
fmt(x) = string(x)
const LOGLOCK = ReentrantLock()
function say(s)
    lock(LOGLOCK) do
        println(s); println(LOG, s); flush(stdout); flush(LOG)
    end
end
emit(ac, quantity, value, bound, verdict) = say(string(ac, " | ", quantity, " | ", fmt(value), " | ", bound, " | ", verdict))
pf(ok::Bool) = ok ? "PASS" : "FAIL"
ratio_line(ac, quantity, r) = emit(ac, quantity, r, "<= 1", pf(r <= 1))

# @info of the tests -> the log
struct TeeLog <: Logging.AbstractLogger end
Logging.min_enabled_level(::TeeLog) = Logging.Info
Logging.shouldlog(::TeeLog, args...) = true
Logging.catch_exceptions(::TeeLog) = true
function Logging.handle_message(::TeeLog, level, message, _mod, group, id, file, line; kwargs...)
    kv = join(("$k = $(first(repr(v), 700))" for (k, v) in kwargs), "; ")
    say(string("# ", level, " ", message, " | ", kv))
end
global_logger(TeeLog())

function section(f, name)
    isempty(ONLY) || name in ONLY || return
    t0 = time()
    try
        f()
    catch e
        say(string(name, " | ERROR | ", sprint(showerror, e, catch_backtrace())[1:min(end, 900)], " | - | FAIL"))
    end
    say(@sprintf("# section %s took %.1f s", name, time() - t0))
end

const T0 = time()
say("# gali verification, julia $(VERSION), threads $(Threads.nthreads()), base seed 20261012 + 40000, GA_FULL_RUN = true, time $(round(time())) (unix)")
say("# timings in this log were taken on a machine that may be loaded by other processes")

# ---------------------------------------------------------------------------------------------------
# The package suite at full size, one section per test file. Counts per testset (cumulative of its children).
include(joinpath(ROOT, "test", "gali", "ga_helpers.jl"))
say("# GA_FULL = $GA_FULL (the full sizes of ga_size are in force)")

# Pass counts are only known for the whole file (DefaultTestSet rolls the counts of finished children into the parent); children are reported with their time and failures.
function nfail(ts)
    n = 0
    for r in ts.results
        r isa Test.AbstractTestSet ? (n += nfail(r)) : (r isa Union{Test.Fail, Test.Error} && (n += 1))
    end
    return n
end
function report_ts(ac, ts; top = true)
    el = ts.time_end - ts.time_start
    if top
        c = Test.get_test_counts(ts)
        ok = c.cumulative_fails == 0 && c.cumulative_errors == 0
        emit(ac, string("file total (", c.cumulative_passes, " pass, ", c.cumulative_fails, " fail, ", c.cumulative_errors, " error, ", c.cumulative_broken, " broken, ", @sprintf("%.1f", el), " s)"),
             c.cumulative_passes, "0 fail, 0 error", pf(ok))
    end
    for ch in ts.results
        ch isa Test.AbstractTestSet || continue
        nf = nfail(ch)
        emit(ac, string("testset ", ch.description, " (", nf, " fail or error, ", @sprintf("%.1f", ch.time_end - ch.time_start), " s)"), nf, "0", pf(nf == 0))
        report_ts(ac, ch; top = false)
    end
end

# Definitions (function, const, assignment) of an earlier test file, without running its testsets: later files use helpers defined there.
function include_defs(file)
    for e in Meta.parseall(read(file, String)).args
        e isa Expr && e.head in (:function, :const, :(=)) && Core.eval(Main, e)
    end
end

function suite_section(name, ac, file; defs = String[])
    section(name) do
        for f in defs
            include_defs(joinpath(ROOT, "test", "gali", f))
        end
        outer = Test.@testset "outer $name" begin
            ts = Test.@testset "$name" begin
                include(joinpath(ROOT, "test", "gali", file))
            end
            report_ts(ac, ts)
        end
    end
end

# `@testset` at the top level throws on failure after reporting; the section catches it and the failure shows in the counts above.
suite_section("AC1-4", "AC1-4", "ga_ac01_04.jl")
suite_section("AC5-8", "AC5-8", "ga_ac05_08.jl")
suite_section("AC9-12", "AC9-12", "ga_ac09_12.jl")
suite_section("AC13-16", "AC13-16", "ga_ac13_16.jl"; defs = ["ga_ac05_08.jl", "ga_ac09_12.jl"])
suite_section("AC17-20", "AC17-20", "ga_ac17_20.jl"; defs = ["ga_ac05_08.jl", "ga_ac09_12.jl", "ga_ac13_16.jl"])
suite_section("D108", "D108", "ga_d108.jl")      # D108j: subnormal column norm is a health failure


# ---------------------------------------------------------------------------------------------------
# From here: sections that are not the package testsets. They need the helper functions of all five test files (definitions only).
const DEFS_DONE = Ref(false)
function need_defs()
    DEFS_DONE[] && return
    for f in ("ga_ac01_04.jl", "ga_ac05_08.jl", "ga_ac09_12.jl", "ga_ac13_16.jl", "ga_ac17_20.jl")
        include_defs(joinpath(ROOT, "test", "gali", f))
    end
    DEFS_DONE[] = true
end
# Top level (so that the definitions are in the world age of the sections below); the suite groups include the test files themselves.
const SUITE = ("AC1-4", "AC5-8", "AC9-12", "AC13-16", "AC17-20")
any(x -> x in SUITE, ONLY) || need_defs()
const VSEED = 20261012 + 40900            # verify-only seeds: 40900 .. 40999
const PERM = [1, 3, 2, 4]                 # module labelling (q1, q2, p1, p2) <-> item-2 labelling (x1, x2 | x3, x4 = pair 1 | pair 2), an involution

# ---------------------------------------------------------------------------------------------------
# DIAG: re-evaluation of the criteria whose package testsets fail at full size, with the cause isolated (reference saturation, missing (1 + a_max) factor, window).
# AC9 (same-orbit bounds). ga_check_orbit of the suite, with three changes: (i) the direct form is gated only where the 256-bit reference is resolved
# (G_ref >= 1e-70, i.e. above the saturation 2^-256 of unit columns rounded to 256 bits; below it the reference volume can be exactly 0 and the bound is NaN);
# (ii) D108a: the log-form bound carries the factor (1 + a_max), a_max the largest coordinate of the lifted orbit up to n (the tangent is evaluated at the rounded
# sums q + p of a lifted state: coupled-maps criterion 1 bound 16 eps (1 + a_max)); the error is also reported in units of eps sum Psi (the old multiplier);
# (iii) D108b: the floor control (d) allows 2 B_meas + B_G^k, as in the suite.
function chk9(b, snaps, K::Int, ::Type{T}, amax::Vector{Float64}) where {T}
    d, epsT = b.d, ga_eps(T)
    bits = T === BigFloat ? 2 * GA_BITS : GA_BITS
    thr = T === BigFloat ? 1.0e-140 : 1.0e-70
    r = (G = 0.0, L = 0.0, Lmult = 0.0, Lamax = 0.0, agree = 0.0, floorbad = 0, nonfinite = 0, skipped = 0, nres = 0)
    for s in snaps
        ic = findfirst(==(s.n), b.cps)
        psi, lam, Gref = b.psi[ic], b.lam[ic], b.ref.G[s.n]
        for k in 2:K
            gr = Float64(Gref[k - 1])
            if gr >= thr
                B = ga_BG(epsT, k, lam[k], gr, psi)
                r = merge(r, (G = max(r.G, ga_err(s.gali[k - 1], Gref[k - 1]; bits = bits) / B), nres = r.nres + 1))
            else
                r = merge(r, (skipped = r.skipped + 1,))
            end
            Bk = ga_BG(epsT, k, lam[k], gr, psi)
            gr < ga_meas(T, k) && gr > 0 && isfinite(Bk) && (r = merge(r, (floorbad = r.floorbad + (s.gali[k - 1] <= 2 * ga_meas(T, k) + Bk ? 0 : 1),)))
        end
        if K == d
            isfinite(s.L) || (r = merge(r, (nonfinite = r.nonfinite + 1,)))
            e = ga_err(s.L, ga_Lref(b.ref, s.n; bits = bits); bits = bits)
            r = merge(r, (L = max(r.L, e / ga_BL(epsT, psi, amax[s.n])), Lmult = max(r.Lmult, e / (epsT * sum(psi))), Lamax = max(r.Lamax, e / (epsT * sum(psi) * (1 + amax[s.n])))))
            gd = Float64(Gref[d - 1])
            if gd >= max(10 * ga_meas(T, d), thr)
                Bd = ga_BG(epsT, d, lam[d], gd, psi) + 8 * epsT * gd * sum(psi)
                r = merge(r, (agree = max(r.agree, Float64(abs(s.gali[d - 1] - exp(s.L))) / Bd),))
            end
        end
    end
    return r
end

section("DIAG-AC9") do
    need_defs()
    say("# DIAG-AC9: cases of the suite at full size (ga_ac9_cases): orbit C n = 1e4, orbit D n = 1e3, 20 chaotic Y3 starts (2, 0.5, 0.1) n = 1e3, 20 chaotic Y4 N = 2 starts n = 1e3, Y4 N = 3 ring, Y4 N = 2 regular n = 1e4, Y4 N = 5 n = 1e3; both types; every multiple of 10")
    worst = Dict{Any, Any}()
    for T in GA_TYPES, (j, (label, c, z, N)) in enumerate(ga_ac9_cases())
        d = ga_d(c)
        b = ga_bundle(c, T, z, (ga_seed(9, 1), j), N)
        pts, _ = ga_pts(ga_sys(T, c), b.x0, N)
        amax = accumulate(max, [Float64(maximum(abs, Float64.(p))) for p in pts])
        for K in unique([4, d])
            r = ga_run1(GALI(), ga_sys(T, c), b.x0, b.W0[:, 1:K] |> ga_sm, N, 1; steps = b.cps, orbit = j)
            q = chk9(b, r.snaps, K, T, amax)
            key = (T, K, endswith(label, "chaotic") ? "chaotic" : "regular/orbit")
            w = get(worst, key, (G = 0.0, L = 0.0, Lmult = 0.0, Lamax = 0.0, agree = 0.0, floorbad = 0, nonfinite = 0, skipped = 0, nres = 0, n = 0))
            worst[key] = (G = max(w.G, q.G), L = max(w.L, q.L), Lmult = max(w.Lmult, q.Lmult), Lamax = max(w.Lamax, q.Lamax), agree = max(w.agree, q.agree),
                          floorbad = w.floorbad + q.floorbad, nonfinite = w.nonfinite + q.nonfinite, skipped = w.skipped + q.skipped, nres = w.nres + q.nres, n = w.n + 1)
        end
    end
    for key in sort(collect(keys(worst)); by = string)
        w = worst[key]; T, K, cl = key
        ratio_line("AC9", "(a) direct form, $T, K = $K, $cl ($(w.n) runs): worst |G - G_ref| / B_G^k over the $(w.nres) resolved (G_ref >= 1e-70) order-handoffs; $(w.skipped) order-handoffs with an unresolved 256-bit reference skipped", w.G)
        emit("AC9", "(d) floor control $T K = $K $cl: handoffs where G_ref < B_meas and the value exceeds 2 B_meas + B_G^k (D108b)", w.floorbad, "0", pf(w.floorbad == 0))
        if w.Lmult > 0
            ratio_line("AC9", "(b) log form, $T, K = d = $K, $cl: worst |L - L_ref| / B_L (B_L = 32 eps (1 + a_max) sum Psi, D108a); nonfinite $(w.nonfinite)", w.L)
            emit("AC9", "(b) same: worst |L - L_ref| / (eps sum Psi) (before D108a: the multiplier B_L would have needed instead of 32)", w.Lmult, "info", "INFO")
            emit("AC9", "(b) same: worst |L - L_ref| / (eps sum Psi (1 + a_max)), a_max the largest lifted coordinate up to n", w.Lamax, "info", "INFO")
            ratio_line("AC9", "(c) $T K = $K $cl: worst |G_d,direct - exp(L)| / (B_G^d + 8 eps G_ref sum Psi) at resolved handoffs", w.agree)
        end
    end
end

# AC11 by window and size of the chain, AC12a by chain length: the failures of the suite at full size are N = 10 (window) and N = 5 (not elliptic).
section("DIAG-AC11") do
    need_defs()
    cases = vcat([(ga_y3(6.0, 0.1, 0.3), "Y3 nu=6"), (ga_y3(4.04, 0.1, 0.3), "Y3 nu=4.04")],
                 [(ga_y4(N, 0.05, 0.04; ring = ring), "Y4 N=$N ring=$ring") for N in (2, 3, 5, 10) for ring in (false, true) if !(ring && N < 3)])
    ns = [100, 200, 400, 800]
    for (lab, sel) in (("N <= 5 and Y3 (d <= 10)", c -> ga_d(c) <= 10), ("N = 10 (d = 20)", c -> ga_d(c) == 20))
        gated_win = lab == "N = 10 (d = 20)" ? "(400, 800)" : "(100, 200)"     # D108d
        for (win, (i1, i2)) in (("(100, 200)", (1, 2)), ("(400, 800)", (3, 4)))
            tot = (ok = 0, bad = 0, unres = 0)
            for (ci, (c, label)) in enumerate(cases)
                sel(c) || continue
                d = ga_d(c); K = d <= 10 ? d : 6
                l = ga_exponents(c)
                E = Float64[Float64(ga_Ek(l, k)) for k in 1:K]
                for j in 1:5, T in (BigFloat, Double64)
                    T === BigFloat && d > 6 && continue
                    r = ga_fixed_run(c, T, K, (ga_seed(11, 1), 10 * ci + j), ns)
                    ok, bad, unres = ga_law_counts(r.snaps[i1], r.snaps[i2], E, K, T)
                    tot = (ok = tot.ok + ok, bad = tot.bad + bad, unres = tot.unres + unres)
                end
            end
            isg = win == gated_win
            emit("AC11", "$lab, window n = $win ($(isg ? "gated, D108d" : "not gated, reported")), 5 frames per system, BigFloat (d <= 6) and Double64: orders within 0.5 / outside / unresolved = $(tot.ok) / $(tot.bad) / $(tot.unres)", tot.bad, isg ? "0 outside" : "reported", isg ? pf(tot.bad == 0) : "INFO")
        end
    end
end

section("DIAG-AC12") do
    need_defs()
    for N in (2, 3, 5)
        ep = N == 5 ? 0.2 : 0.3                                            # D108e
        g = [-1.0 + ep * (2 - 2cos(pi * m / N)) for m in 0:(N - 1)]       # open-chain mode values g_m = K + eps lambda_m, lambda_m = 2 - 2 cos(pi m / N)
        emit("AC12", "Y4 N=$N K = -1, eps = $ep: largest mode value g = K + eps lambda_max (elliptic iff -4 < g < 0)", maximum(g), "< 0 for a stable origin", pf(maximum(g) < 0))
    end
    for eps in (0.2,)
        c = ga_y4(5, -1.0, eps)
        for j in 1:5
            d = ga_d(c)
            r = ga_fixed_run(c, Float64, d, (ga_seed(12, 1), 50 + j), collect(1:2000))
            inr = all(s -> all(g -> 0.05 <= g <= 1 + ga_meas(Float64, d), s.gali), r.snaps)
            sl = ga_slopes(filter(s -> s.n % 10 == 0, r.snaps), 100, 2000, d, Float64)
            emit("AC12", "Y4 N=5 eps=$eps (D108e) frame $j: G_k in [0.05, 1 + B_meas] for all k, n <= 2000; max |slope| over [100, 2000] = $(fmt(maximum(abs, sl)))", maximum(abs, sl), "G in range, |slope| <= 0.2", pf(inr && all(s -> abs(s) <= 0.2, sl)))
        end
    end
end

# AC13: only the offsets whose singular values from the long 256-bit product are valid (the suite's ga_logsv at n = 800 gives a sum of the log singular values far from 0
# for some strongly chaotic starts: symplectic pairing requires 0; the module log form equals the 512-bit reference there).
section("DIAG-AC13") do
    need_defs()
    rng = StableRNG(ga_seed(13))
    cases = Any[("orbit D", GA_Y3A, GA_ORBIT_D, [100, 400, 800, 1600])]
    for _ in 1:20; push!(cases, ("Y3 chaotic", GA_Y3B, ga_draw_chaotic(rng, GA_Y3B), [100, 400, 800])); end
    for N in (2, 3), _ in 1:20; push!(cases, ("Y4 N=$N chaotic", ga_y4(N, 3.0, 0.6102), ga_draw_chaotic(rng, ga_y4(N, 3.0, 0.6102)), [100, 400, 800])); end
    ninvalid, nchk, worstc, worstd, ngated = 0, 0, 0.0, 0.0, 0
    for (j, (label, c, z, ns)) in enumerate(cases)
        res = setprecision(BigFloat, GA_BITS) do
            sys, x0 = ga_sys(BigFloat, c), ga_state(BigFloat, c, z)
            pts, _ = ga_pts(sys, x0, maximum(ns))
            lsv = ga_logsv([ga_tanB(sys, c, p) for p in pts], ns)
            Dict(n => Float64(sum(lsv[n])) for n in ns)
        end
        o = ga_offsets(c, z, (ga_seed(13, 1), j), ns)
        valid = [n for n in ns if abs(res[n]) < 1.0e-6]
        length(valid) < length(ns) && (ninvalid += 1)
        o.gap[400] <= 0.01 && continue
        ngated += 1
        for n in valid, v in o.c[n]
            isnan(v) || (worstc = max(worstc, abs(v)))
        end
        if 100 in valid
            for n2 in valid
                n2 == 100 && continue
                for (a, b) in zip(o.c[100], o.c[n2])
                    (isnan(a) || isnan(b)) && continue
                    nchk += 1; worstd = max(worstd, abs(b - a))
                end
            end
        end
    end
    emit("AC13", "orbits with an invalid singular-value sum at some n in the suite's ga_logsv (|sum ln sigma_i| >= 1e-6, symplectic pairing violated; those n are dropped), of $(length(cases)) ($ngated gated)", ninvalid, "info", "INFO")
    emit("AC13", "offsets at the valid n: worst |c_k(n)| ($nchk offset pairs checked)", worstc, "<= 15", pf(worstc <= 15))
    emit("AC13", "offsets at the valid n: worst |c_k(n_2) - c_k(100)| (D108f)", worstd, "<= 5.0", pf(worstd <= 5.0))
end

# AC2 realised ratio (the suite asserts it without printing): worst |G_T - sqrt(det Gram)| / (3 B_meas), as the suite's loop.
section("R3-AC2") do
    need_defs()
    rng = StableRNG(VSEED + 1)
    for T in GA_TYPES
        worst = Ref(0.0)
        function chk(Gram, k, d, closed)
            local W, ex
            setprecision(BigFloat, GA_BITS) do
                V = ga_gram_vectors(Gram, d, ga_randorth(rng, d))
                W = ga_sm([ga_to(T, V[i, j]) for i in 1:d, j in 1:k])
                ex = closed(Gram)
            end
            worst[] = max(worst[], ga_err(gali_values(W)[k - 1], ex) / (3 * ga_meas(T, k)))
        end
        eqgram(k, c) = (G = fill(c, k, k); for i in 1:k; G[i, i] = one(c); end; G)
        for m in 1:14, d in (2, 4, 8)
            th = setprecision(() -> BigFloat(10.0)^(-m), BigFloat, GA_BITS)
            for ang in (th, setprecision(() -> BigFloat(pi) / 2 - BigFloat(10.0)^(-m), BigFloat, GA_BITS))
                setprecision(BigFloat, GA_BITS) do
                    c = cos(ang); chk(BigFloat[1 c; c 1], 2, d, G -> sin(ang))
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
        ratio_line("AC2", "R3 worst |G_T - sqrt(det Gram)| / (3 B_meas(T, k)), $T, constructed Gram matrices (a) to (c)", worst[])
    end
end

# ---------------------------------------------------------------------------------------------------
# Recorded lines R1, R2, R4
section("R1") do
    need_defs()
    c = GA_Y3B; N = 300; m = 1000
    rng = StableRNG(VSEED + 2)
    thr = (1.0e-3, 1.0e-6, 1.0e-10)
    first_below(v, t) = (i = findfirst(<(t), v); i === nothing ? 0 : i)
    tauG(t) = t * sqrt(4 - t^2) / 2                 # B4 threshold map, SALI threshold -> GALI_2 threshold
    sali_of(g) = sqrt(2) * g / sqrt(1 + sqrt(max(1 - g^2, 0.0)))
    cnt = Dict(t => 0 for t in thr); maxd = Dict(t => 0 for t in thr); nocross = Dict(t => 0 for t in thr); mapdiff = Dict(t => 0 for t in thr)
    sys = ga_sys(Float64, c)
    for j in 1:m
        z = ga_draw_chaotic(rng, c)
        x0 = ga_state(Float64, c, z); W0 = ga_frame(VSEED, j, Float64, 4, 2)
        rg = ga_run1(GALI(), sys, x0, W0, N, 1; all = true, orbit = j)
        rs = ga_run1(SALI(), sys, x0, W0, N, 1; all = true, orbit = j)
        g2 = [Float64(s.gali[1]) for s in rg.snaps]; sa = [Float64(s.sali) for s in rs.snaps]
        for t in thr
            a = first_below(g2, tauG(t)); b = first_below(sa, t); c2 = first_below(sali_of.(g2), t)
            a == 0 && b == 0 && (nocross[t] += 1)
            a != b && (cnt[t] += 1; maxd[t] = max(maxd[t], abs(a - b)))
            c2 != b && (mapdiff[t] += 1)
        end
    end
    for t in thr
        emit("R1", "$m chaotic Y3 starts (2, 0.5, 0.1), n <= $N, tau_S = $t: crossing step of G_2 (threshold tau_g = tau_S sqrt(4 - tau_S^2)/2) vs crossing step of SALI: starts with different steps / max |difference| / starts with no crossing; G_2 mapped to SALI by B4 against tau_S: starts that differ", string(cnt[t], " / ", maxd[t], " / ", nocross[t], "; ", mapdiff[t]), "recorded", "INFO")
    end
end

section("R2") do
    need_defs()
    rng = StableRNG(VSEED + 3)
    function cross(v, t)
        i = findfirst(<(t), v); i === nothing ? 0 : i
    end
    for (label, c, starts, N, dense) in (("orbit D", GA_Y3A, [GA_ORBIT_D], 20_000, false), ("orbit C", GA_Y3A, [GA_ORBIT_C], 20_000, false),
                                         ("100 chaotic Y3 starts (2, 0.5, 0.1)", GA_Y3B, [ga_draw_chaotic(rng, GA_Y3B) for _ in 1:100], 3000, true),
                                         ("100 regular small-amplitude Y3 starts (0.5, 0.1, 1e-3)", GA_Y3A, [ga_draw_regular(rng, GA_Y3A) for _ in 1:100], 20_000, false))
        res = Dict((k, t) => Int[] for k in (3, 4), t in (1.0e-8, 1.0e-16))
        for (j, z) in enumerate(starts)
            x0 = ga_state(Float64, c, z); W0 = ga_frame(VSEED, j, Float64, 4, 4)
            steps = dense ? Int[] : collect(10:10:N)
            r = ga_run1(GALI(), ga_sys(Float64, c), x0, W0, N, 1; steps = steps, all = dense, orbit = j)
            g3 = [Float64(s.gali[2]) for s in r.snaps]; Lg = [Float64(s.L) for s in r.snaps]          # G_3 direct, ln G_4 by the log form
            ns = [s.n for s in r.snaps]
            a = cross(g3, 1.0e-8); push!(res[(3, 1.0e-8)], a == 0 ? 0 : ns[a])
            a = cross(Lg, log(1.0e-8)); push!(res[(4, 1.0e-8)], a == 0 ? 0 : ns[a])
            a = cross(Lg, log(1.0e-16)); push!(res[(4, 1.0e-16)], a == 0 ? 0 : ns[a])
            rD = ga_run1(GALI(), ga_sys(Double64, c), ga_state(Double64, c, z), ga_frame(VSEED, j, Double64, 4, 4), N, 1; steps = steps, all = dense, orbit = j)
            g3d = [Float64(s.gali[2]) for s in rD.snaps]
            a = cross(g3d, 1.0e-16); push!(res[(3, 1.0e-16)], a == 0 ? 0 : ns[a])
        end
        for k in (3, 4), t in (1.0e-8, 1.0e-16)
            v = res[(k, t)]; cr = filter(>(0), v)
            how = k == 4 ? "ln G_4 from the log form, Float64" : (t == 1.0e-8 ? "direct form, Float64" : "direct form, Double64 (the Float64 floor 1.3e-15 is above the threshold)")
            emit("R2", "$label, n <= $N ($(dense ? "every step" : "every 10th step")): crossing step of G_$k at $t ($how): crossed $(length(cr)) of $(length(v)); min / median / max step = $(isempty(cr) ? "-" : string(minimum(cr), " / ", sort(cr)[(end + 1) ÷ 2], " / ", maximum(cr)))", length(cr), "recorded", "INFO")
        end
    end
end

section("R4") do
    need_defs()
    # orbit D: the depth to which the direct form tracks the exact value (F5): computed G_k against the 256-bit value (k = 4 also exp(L))
    c = GA_Y3A; N = 3000; cps = [400, 800, 1500, 3000]
    z = GA_ORBIT_D
    W0 = ga_frame(VSEED, 1, Float64, 4, 4)
    pts, _ = ga_pts(ga_sys(Float64, c), ga_state(Float64, c, z), N)
    ref = ga_ref(c, pts, W0, cps; all_steps = false)
    for T in GA_TYPES
        r = ga_run1(GALI(), ga_sys(T, c), ga_state(T, c, z), T === Float64 ? W0 : ga_frame(VSEED, 1, T, 4, 4), N, 1; steps = cps, orbit = 1)
        for s in r.snaps, k in 2:4
            gr = ref.G[s.n][k - 1]
            emit("R4", "orbit D, $T direct form G_$k at n = $(s.n): computed $(fmt(Float64(s.gali[k - 1]))) vs exact (256-bit) $(fmt(Float64(gr))); floor B_meas $(fmt(ga_meas(T, k)))" * (k == 4 ? "; exp(L) = $(fmt(exp(Float64(s.L))))" : ""), abs(Float64(s.gali[k - 1]) - Float64(gr)), "recorded (the exact value is below 1e-77 beyond n about 1500 for k = 4: saturation of the 256-bit unit columns)", "INFO")
        end
    end
end

# ---------------------------------------------------------------------------------------------------
# AC19: ChaosTools gali
ct3_rule(u, p, n) = (nu = p[1]; ka = p[2]; mu = p[3]; s1 = u[1] + u[2]; s2 = u[3] + u[4]; c = mu * (1 - cos(u[1] + u[2] + u[3] + u[4]));
                     SVector(s1, u[2] - nu * sin(s1) - c, s2, u[4] - ka * sin(s2) - c))
function ct3_jac(u, p, n)
    nu, ka, mu = p[1], p[2], p[3]
    S1 = u[1] + u[2]; S2 = u[3] + u[4]; sT = sin(u[1] + u[2] + u[3] + u[4])
    a = -nu * cos(S1) - mu * sT; b = -mu * sT; cc = -ka * cos(S2) - mu * sT
    return @SMatrix [1.0 1.0 0.0 0.0; a 1+a b b; 0.0 0.0 1.0 1.0; b b cc 1+cc]
end
ct3_tds(p, K) = CT.TangentDynamicalSystem(CT.DeterministicIteratedMap(ct3_rule, SVector(0.1, 0.1, 0.1, 0.1), collect(p)); J = ct3_jac, k = K)
# oracle series on one start: GALI_K(t), t = 0..n, and the orbit points before each step (module labelling)
function ct_series(p, K, z, W0, n)
    tds = ct3_tds(p, K)
    Q0 = SMatrix{4, K, Float64}(W0[PERM, 1:K])
    CT.reinit!(tds, SVector{4, Float64}(z[PERM]); Q0 = Q0)
    states = Vector{SVector{4, Float64}}(undef, n)
    for i in 1:n
        states[i] = SVector{4, Float64}(CT.current_state(tds)[PERM])
        CT.step!(tds)
    end
    CT.reinit!(tds, SVector{4, Float64}(z[PERM]); Q0 = Q0)
    g = CT.gali(tds, n; threshold = 0)[1]
    return g, states
end

section("AC19") do
    need_defs()
    say("# AC19: ChaosTools $(pkgversion(CT)) gali(tands, T; threshold = 0, Delta t = 1), TangentDynamicalSystem with the hand-coded Jacobian ct3_jac (item-2 labelling of Y3, permuted by PERM), seeded frame by reinit!(Q0 = first K columns of the module frame, rows permuted)")
    # (b) confirmation first: K = 0 shear (nu = kappa = mu = 0) with W0 = I, closed forms G_3 = 1/sqrt(n^2 + 1), G_4 = 1/(n^2 + 1); and orbit C slopes
    ok_b = true
    for K in (3, 4), n in (10, 100, 1000)
        g, _ = ct_series((0.0, 0.0, 0.0), K, zeros(4), Matrix{Float64}(I, 4, 4), n)
        exact = K == 3 ? 1 / sqrt(n^2 + 1.0) : 1 / (n^2 + 1.0)
        e = abs(g[n + 1] - exact); okk = e <= 4 * ga_meas(Float64, K)
        ok_b &= okk
        emit("AC19", "(b) confirmation, K = 0 shear, W0 = I, ChaosTools GALI_$K at n = $n: |value - closed form| (closed form $(fmt(exact)))", e, "<= 4 B_meas = $(fmt(4 * ga_meas(Float64, K)))", pf(okk))
    end
    ns = collect(1000:10:20_000)
    for k in (3, 4)
        gC, _ = ct_series((0.5, 0.1, 1.0e-3), k, [0.5, 0.5, 0.0, 0.0], ga_frame(VSEED, 2, Float64, 4, 4), 20_000)
        sl = ga_slope(log.(Float64.(ns)), log.(gC[ns .+ 1]))
        okk = k == 3 ? -2.3 <= sl <= -1.7 : -4.3 <= sl <= -3.7
        ok_b &= okk
        emit("AC19", "(b) confirmation, orbit C, ChaosTools GALI_$k slope of ln G against ln n over [1e3, 2e4] (B7: $(k == 3 ? -2 : -4))", sl, k == 3 ? "[-2.3, -1.7]" : "[-4.3, -3.7]", pf(okk))
    end
    emit("AC19", "(b) the behaviour holds (ChaosTools gali accepts k = 3, 4 on the hand-written Y3 and reproduces the closed forms and the regular-orbit slopes)", ok_b, "true, else withdrawn", ok_b ? "PASS" : "WITHDRAWN")
    cases = Any[("orbit C", GA_Y3A, GA_ORBIT_C, (0.5, 0.1, 1.0e-3)), ("orbit D", GA_Y3A, GA_ORBIT_D, (0.5, 0.1, 1.0e-3))]
    rng = StableRNG(VSEED + 4)
    for _ in 1:20; push!(cases, ("chaotic Y3 (2, 0.5, 0.1)", GA_Y3B, ga_draw_chaotic(rng, GA_Y3B), (2.0, 0.5, 0.1))); end
    for _ in 1:20; push!(cases, ("regular Y3 (0.5, 0.1, 1e-3)", GA_Y3A, ga_draw_regular(rng, GA_Y3A), (0.5, 0.1, 1.0e-3))); end
    nmax = 1000; cps = [10, 100, 1000]
    stats = Dict{Any, Any}()
    nbit = Ref(0); ntot = Ref(0); maxorb = Ref(0.0)
    for (j, (label, c, z, p)) in enumerate(cases)
        sysF = ga_sys(Float64, c)
        x0 = ga_state(Float64, c, z); W0 = ga_frame(VSEED + 7, j, Float64, 4, 4)
        pts_m, _ = ga_pts(sysF, x0, nmax)
        gO, states = ct_series(p, 4, z, W0, nmax)
        ntot[] += 1; nbit[] += pts_m == states ? 1 : 0
        maxorb[] = max(maxorb[], maximum(maximum(abs, pts_m[i] - states[i]) for i in 1:nmax))
        # module bundle (own orbit) and oracle reference (at the oracle's own orbit points), K = 4 columns of the seeded frame
        bm = ga_bundle_core(sysF, c, x0, W0, nmax, cps)
        ref_o = ga_ref(c, states, W0, cps)
        Ms = [Float64.(Matrix(tangent(sysF, q))) for q in states]
        psi_o = ga_psis(Ms, W0, cps)
        lam_o = [[ga_Lambda(ref_o.lnG[k, :], n) for k in 1:4] for n in cps]
        for K in (2, 3, 4)
            ok_b || K == 2 || continue
            gO_K = K == 4 ? gO : ct_series(p, K, z, W0, nmax)[1]
            r = ga_run1(GALI(), sysF, x0, W0[:, 1:K] |> ga_sm, nmax, 1; steps = cps, orbit = j)
            cls = (label == "orbit D" || startswith(label, "chaotic")) ? :chaotic : :regular
            for (ic, n) in enumerate(cps)
                k = K                       # ChaosTools gali returns the order K = number of deviation vectors only; orders k < K come from runs with K = k
                grm = Float64(bm.ref.G[n][k - 1]); gro = Float64(ref_o.G[n][k - 1])
                stats[(K, k, :tot, cls)] = get(stats, (K, k, :tot, cls), 0) + 1
                if !(grm >= 10 * ga_meas(Float64, k) && gro >= 10 * ga_meas(Float64, k))
                    stats[(K, k, :unres, cls)] = get(stats, (K, k, :unres, cls), 0) + 1
                    continue
                end
                Bm = ga_BG(ga_eps(Float64), k, bm.lam[ic][k], grm, bm.psi[ic]); Bo = ga_BG(ga_eps(Float64), k, lam_o[ic][k], gro, psi_o[ic])
                gm = Float64(r.snaps[ic].gali[k - 1]); go = gO_K[n + 1]
                rm = abs(gm - grm) / Bm; ro = abs(go - gro) / Bo; pair = abs(gm - go) / (Bm + Bo)
                st = get(stats, (K, k), (rm = 0.0, ro = 0.0, pair = 0.0, n = 0, inb = 0, maxd = 0.0))
                stats[(K, k)] = (rm = max(st.rm, rm), ro = max(st.ro, ro), pair = max(st.pair, pair), n = st.n + 1, inb = st.inb + (pair <= 1 ? 1 : 0), maxd = max(st.maxd, abs(gm - go)))
            end
        end
    end
    emit("AC19", "oracle orbit equals the stream orbit bitwise over n <= $nmax: starts with a bitwise equal orbit / max coordinate difference $(fmt(maxorb[]))", string(nbit[], " of ", ntot[]), "bitwise, else own-orbit references and the pair bound", "INFO")
    for K in (2, 3, 4)
        k = K
        haskey(stats, (K, k)) || continue
        st = stats[(K, k)]; un = sum(get(stats, (K, k, :unres, cl), 0) for cl in (:chaotic, :regular))
        for cl in (:chaotic, :regular)
            t = get(stats, (K, k, :tot, cl), 0); u = get(stats, (K, k, :unres, cl), 0)
            reg = cl == :regular        # D108h: the 10 percent allowance gates the regular starts; chaotic starts are reported
            emit("AC19", "order k = $k: excluded as unresolved (reference below 10 B_meas), $cl starts: $u of $t start-handoffs ($(round(100 * u / max(t, 1); digits = 1)) percent)", u / max(t, 1), reg ? "<= 10 percent (regular starts, D108h)" : "reported (chaotic: a property of the Float64 floor, D108h)", reg ? pf(u / max(t, 1) <= 0.1) : "INFO")
        end
        tag = k == 2 ? "(a)" : "(b)"
        ratio_line("AC19", "$tag run K = $K, order k = $k, $(st.n) resolved start-handoffs ($un excluded as unresolved, $(round(100 * un / (st.n + un); digits = 1)) percent): module max |G - G_ref| / B_G^k at its own orbit", st.rm)
        ratio_line("AC19", "$tag same for the ChaosTools value against the exact value at the oracle's own orbit (B_oracle = the same form, eps(Float64))", st.ro)
        emit("AC19", "$tag same: pairs with |G_module - G_oracle| <= B_module + B_oracle: $(st.inb) of $(st.n); max pair ratio $(fmt(st.pair)); max absolute difference $(fmt(st.maxd))", st.pair, "<= 1 where the orbits are bitwise equal; else the orbit divergence term (reported)", st.inb == st.n ? "PASS" : "INFO")
    end
end

# ---------------------------------------------------------------------------------------------------
# AC20 (V): scale and time against the L1 k = 1 pass, Y4 N = 100
section("AC20t") do
    need_defs()
    for ring in (false, true)
        c = ga_y4(100, 3.0, 0.2; ring = ring)
        nm = ring ? "ring" : "open"
        sys = ga_sys(Float64, c)
        rng = StableRNG(VSEED + (ring ? 6 : 5))
        norb, n = 200, 1000
        x0 = [ga_state(Float64, c, [(2 * rand(rng) - 1) * 4 * pi for _ in 1:200]) for _ in 1:norb]
        Wg = [ga_randset(rng, Float64, 200, 4; spread = false) for _ in 1:norb]      # Gaussian frames (D107l): the 200 x 200 Haar QR of the default frame is out of scope
        kw(k, red) = (k = k, tau = 1, seed = VSEED, reducers = red, backend = :cpu, W0 = [w[:, 1:k] |> ga_sm for w in Wg])
        t0 = time(); propagate(sys, x0[1:2], 3; kw(4, (GALI(),))..., W0 = kw(4, (GALI(),)).W0[1:2]); say(@sprintf("# AC20 Y4 N = 100 %s: first call GALI K = 4 (compile) %.1f s", nm, time() - t0))
        t0 = time(); propagate(sys, x0[1:2], 3; kw(1, (Benettin(),))..., W0 = kw(1, (Benettin(),)).W0[1:2]); say(@sprintf("# AC20 Y4 N = 100 %s: first call Benettin k = 1 (compile) %.1f s", nm, time() - t0))
        t0 = time(); propagate(sys, x0[1:2], 3; kw(2, (GALI(),))..., W0 = kw(2, (GALI(),)).W0[1:2]); say(@sprintf("# AC20 Y4 N = 100 %s: first call GALI K = 2 (compile) %.1f s", nm, time() - t0))
        tg = @elapsed rg = propagate(sys, x0, n; kw(4, (GALI(),))...)
        tl = @elapsed rl = propagate(sys, x0, n; kw(1, (Benettin(),))...)
        t2 = @elapsed r2 = propagate(sys, x0, n; kw(2, (GALI(),))...)
        fin = [f[1] for f in rg.finals]
        inrange = all(f -> all(g -> 0 <= g <= 1 + ga_meas(Float64, 4), f.gali), fin)
        nested = all(f -> f.gali[1] + 2 * ga_meas(Float64, 3) >= f.gali[2] && f.gali[2] + 2 * ga_meas(Float64, 4) >= f.gali[3], fin)
        finite = all(f -> all(isfinite, f.gali), fin)
        bitw = count(j -> rg.finals[j][1].gali[1] === r2.finals[j][1].gali[1], 1:norb)
        emit("AC20", "Y4 N = 100 $nm, K = 4, $norb seeded orbits, n = $n, Float64 CPU ($(Threads.nthreads()) threads): G_k in [0, 1 + B_meas] / nested / finite", string(inrange, " / ", nested, " / ", finite), "true / true / true", pf(inrange && nested && finite))
        emit("AC20", "$nm: G_2 of the K = 4 run against the K = 2 run on the identical orbit (nesting): bitwise equal orbits", string(bitw, " of ", norb), "within 2 B_meas (bitwise reported)", pf(all(j -> abs(rg.finals[j][1].gali[1] - r2.finals[j][1].gali[1]) <= 2 * ga_meas(Float64, 2), 1:norb)))
        emit("AC20", "$nm: time per orbit-step, GALI K = 4 / GALI K = 2 / Benettin k = 1 (ns)", string(fmt(tg / (norb * n) * 1e9), " / ", fmt(t2 / (norb * n) * 1e9), " / ", fmt(tl / (norb * n) * 1e9)), "reported", "INFO")
        emit("AC20", "$nm: ratio of the K = 4 GALI pass to the L1 k = 1 pass on the identical orbit (cm criterion 18 style)", tg / tl, "reported", "INFO")
    end
end

# ---------------------------------------------------------------------------------------------------
# Figure: ln G_k against n (orbit D) and ln n (orbit C)
section("figure") do
    need_defs()
    c = GA_Y3A
    W0 = ga_frame(VSEED, 1, Float64, 4, 4)
    # orbit D, n <= 3000: k = 2, 3 direct Float64 where resolved (>= 10 B_meas), k = 4 by the log form (K = d); predicted slope -E_k from the exponents of the 256-bit product
    ND = 3000
    rD = ga_run1(GALI(), ga_sys(Float64, c), ga_state(Float64, c, GA_ORBIT_D), W0, ND, 1; steps = collect(10:10:ND), orbit = 1)
    lsv = setprecision(BigFloat, GA_BITS) do
        sys = ga_sys(BigFloat, c); x0 = ga_state(BigFloat, c, GA_ORBIT_D)
        pts, _ = ga_pts(sys, x0, ND)
        ga_logsv([ga_tanB(sys, c, p) for p in pts], [ND])[ND]
    end
    l = Float64.(lsv) ./ ND
    Ek(k) = sum(l[1] - l[j] for j in 2:k)
    say(@sprintf("# figure: orbit D exponents at n = 3000 from the 256-bit product: l = (%.5g, %.5g, %.5g, %.5g); E_2 = %.5g, E_3 = %.5g, E_4 = %.5g (E_4 = 4 l_1 = %.5g)", l..., Ek(2), Ek(3), Ek(4), 4 * l[1]))
    nsD = [s.n for s in rD.snaps]
    NC = 20_000
    rC = ga_run1(GALI(), ga_sys(Float64, c), ga_state(Float64, c, GA_ORBIT_C), W0, NC, 1; steps = collect(10:10:NC), orbit = 2)
    nsC = [s.n for s in rC.snaps]
    fig = CM.Figure(size = (1300, 520))
    ax1 = CM.Axis(fig[1, 1]; xlabel = "n", ylabel = "ln G_k(n)", title = "orbit D (chaotic): ln G_k against n")
    cols = (:steelblue, :darkorange, :crimson)
    for (i, k) in enumerate(2:4)
        v = k == 4 ? [Float64(s.L) for s in rD.snaps] : [Float64(s.gali[k - 1]) >= 10 * ga_meas(Float64, k) ? log(Float64(s.gali[k - 1])) : NaN for s in rD.snaps]
        CM.lines!(ax1, nsD, v; color = cols[i], linewidth = 2.2, label = k == 4 ? "ln G_4 (log form)" : "ln G_$k (direct, above 10 B_meas)")
        n0 = 500; i0 = findfirst(==(n0), nsD); off = v[i0] + Ek(k) * n0
        CM.lines!(ax1, nsD, off .- Ek(k) .* nsD; color = cols[i], linestyle = :dash, linewidth = 1.2, label = @sprintf("slope -E_%d = %.4f (offset at n = 500)", k, -Ek(k)))
    end
    CM.axislegend(ax1; position = :lb, labelsize = 10)
    ax2 = CM.Axis(fig[1, 2]; xscale = log10, xlabel = "n", ylabel = "ln G_k(n)", title = "orbit C (regular): ln G_k against n, slopes 0, -2, -4")
    for (i, k) in enumerate(2:4)
        v = [log(Float64(s.gali[k - 1])) for s in rC.snaps]
        CM.lines!(ax2, nsC, v; color = cols[i], linewidth = 1.6, label = "ln G_$k")
        i0 = findfirst(==(1000), nsC); sl = (0.0, -2.0, -4.0)[i]
        CM.lines!(ax2, nsC[i0:end], v[i0] .+ sl .* (log.(Float64.(nsC[i0:end])) .- log(1000.0)); color = cols[i], linestyle = :dash, linewidth = 1.2, label = "slope $(Int(sl)) in ln n")
    end
    CM.hlines!(ax2, [log(ga_meas(Float64, 4) * 10)]; color = :gray50, linestyle = :dot, label = "10 B_meas(Float64, 4)")
    CM.axislegend(ax2; position = :lb, labelsize = 10)
    CM.save(joinpath(OUT, "gali-laws.png"), fig)
    say("# figure saved: gali-laws.png")
end

say(@sprintf("# total wall time %.1f s", time() - T0))
close(LOG)
