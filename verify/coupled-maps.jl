# Stage 04 verification of module coupled-maps (Systems: Y3 Froeschle 4D map, Y4 chain of N coupled standard maps, wrapped forms, D105 fixed-size states).
# Run from the repo root (main run, 8 threads):   julia --project=verify -t 8 verify/coupled-maps.jl
# Thread comparison (D104c, criterion 18):         CM_ONLY=threads CM_APPEND=1 julia --project=verify -t 1 verify/coupled-maps.jl
#                                                   CM_ONLY=compare CM_APPEND=1 julia --project=verify -t 1 verify/coupled-maps.jl
# (the main run writes the per-orbit hashes of the 8-thread run, the -t 1 run those of the 1-thread run, "compare" compares them bitwise)
# Spec: workspaces/module-build/stages/01-spec/output/coupled-maps-spec.md (criteria 1 to 21, Amendments D104, D105).
# One printed line per measurement:  AC<n> | quantity | value | bound | PASS/FAIL/INFO
# Most "value" entries are the worst ratio error / bound (pass when <= 1), as in the package tests (test/coupled-maps/).
# The measurement functions of the package suite are reused (the test files are included with @testset / @test turned into no-ops), so the
# assertions and tolerances are identical; only the sample sizes are the full D38 sizes. Seeds: StableRNG(20261012 + offset), offsets from 41000.
# Environment: CM_ONLY="AC1-4,AC20" runs only the named sections; CM_LOG overrides the log file; CM_APPEND=1 appends; CM_SCALE=0.01 scales every sample
# size (smoke test only); CM_DUMP is the directory of the thread-comparison dumps (default: a temp directory, never inside the repo).
using Printf, Random, LinearAlgebra, Statistics, Serialization
using Chaoslab, StaticArrays, DoubleFloats, StableRNGs, ForwardDiff, BenchmarkTools
import ChaosTools as CT
import CairoMakie as CM

const ROOT = dirname(@__DIR__)
const OUT = joinpath(ROOT, "workspaces", "module-build", "stages", "04-verify", "output")
mkpath(OUT)
const ONLY = filter(!isempty, split(get(ENV, "CM_ONLY", ""), ","))
const SC = parse(Float64, get(ENV, "CM_SCALE", "1"))
sz(n) = max(1, round(Int, n * SC))
const LOG = open(get(ENV, "CM_LOG", joinpath(OUT, "coupled-maps-verify-log.txt")), get(ENV, "CM_APPEND", "") == "1" ? "a" : "w")
const DUMP = get(ENV, "CM_DUMP", joinpath(tempdir(), "cm-threads"))

fmt(x::AbstractFloat) = @sprintf("%.4g", x)
fmt(x) = string(x)
function say(s)
    println(s); println(LOG, s); flush(stdout); flush(LOG)
end
emit(ac, quantity, value, bound, verdict) = say(string(ac, " | ", quantity, " | ", fmt(value), " | ", bound, " | ", verdict))
pf(ok::Bool) = ok ? "PASS" : "FAIL"
ratio_line(ac, quantity, r) = emit(ac, quantity, r, "<= 1", pf(r <= 1))

# Reuse the measurement functions of the package tests: the test files are included with the test macros turned into no-ops.
macro testset(args...) nothing end
macro test(args...) nothing end
macro test_throws(args...) nothing end
const TESTDIR = joinpath(ROOT, "test", "coupled-maps")
for f in ("cm_helpers.jl", "cm_ac01_04.jl", "cm_ac05_06.jl", "cm_ac07_12.jl", "cm_ac13_16.jl", "cm_ac17_19.jl", "cm_d105.jl")
    include(joinpath(TESTDIR, f))
end

function section(f, name)
    isempty(ONLY) || name in ONLY || return
    t0 = time()
    try
        f()
    catch e
        say(string(name, " | ERROR | ", sprint(showerror, e, catch_backtrace())[1:min(end, 600)], " | - | FAIL"))
    end
    say(@sprintf("# section %s took %.1f s", name, time() - t0))
end

const T0 = time()
say("# coupled-maps verification, julia $(VERSION), threads $(Threads.nthreads()), base seed $CM_SEED, scale $SC, time $(round(time())) (unix)")
say("# timings in this log were taken on a machine that may be loaded by other processes")
say("# n_H (Y3 published) = $(cm_nH(cm_y3(0.5, 0.1, 1.0e-3))), n_H (Y3 (2, 0.5, 0.1)) = $(cm_nH(cm_y3(2.0, 0.5, 0.1))); sigma published = $(fmt(cm_sigma(cm_y3(0.5, 0.1, 1.0e-3)))) (spec: 2.07, n_H 19)")

# extra N = 100 cases (criteria 1, 2, 4, 10 name N = 100; the package panel stops at N = 10)
cm_big_cases(rng, m) = [(cm_y4(100, K, e; ring), cm_draw(rng, 200, 4 * Float64(pi))) for ring in (false, true), (K, e) in ((3.0, 0.6102), (5.0, 1.0)), _ in 1:m] |> vec

# ---------------------------------------------------------------------------------------------------
# AC1 to AC4 (1000 points per case, 200 orbits per case; N = 100 added at 5 points / orbits per case)
section("AC1-4") do
    nf = sz(1000)
    cF = cm_cases(StableRNG(CM_SEED + 41001), nf)
    cD = cm_cases(StableRNG(CM_SEED + 41002), nf; maxN = 5)
    cB = cm_big_cases(StableRNG(CM_SEED + 41003), sz(5))
    say("# AC1 sample: Float64 $(length(cF)) points (Y3: 5 sets + uniform, $nf each; Y4 panel of 8 configs N <= 10, $nf each), Double64 $(length(cD)), N = 100 Float64 $(length(cB))")
    for (nm, cs) in (("N <= 10", cF), ("N = 100", cB))
        ad, cs_, rf = cm_ac1_ratios(Float64, cs; complexstep = true)
        ratio_line("AC1", "max ratio |M - J_ForwardDiff| / 16 eps (1 + a_max), Float64, $nm", ad)
        ratio_line("AC1", "max ratio |M - J_complexstep(1e-30)| / 16 eps (1 + a_max), Float64, $nm", cs_)
        ratio_line("AC1", "max ratio |M - plain-formula tangent| / 16 eps (1 + a_max), Float64, $nm", rf)
    end
    ad, _, rf = cm_ac1_ratios(Double64, cD; complexstep = false)
    ratio_line("AC1", "max ratio |M - J_ForwardDiff| / 16 eps (1 + a_max), Double64 (eps 2^-100), $(length(cD)) pts", ad)
    ratio_line("AC1", "max ratio |M - plain-formula tangent|, Double64", rf)
    emit("AC1", "same Double64 ForwardDiff ratio against nominal eps 2^-104 (ratio x16)", 16 * ad, "<= 1", pf(16 * ad <= 1))
    for (nm, T, cs) in (("Float64, N <= 10", Float64, cm_cases(StableRNG(CM_SEED + 42001), nf)), ("Double64, N <= 5", Double64, cm_cases(StableRNG(CM_SEED + 42002), nf; maxN = 5)), ("Float64, N = 100", Float64, cB))
        a, b, r = cm_ac2_ratios(T, cs)
        ratio_line("AC2", "max ratio |M M_inv - I| / 16 eps (1 + a_max)^2, $nm", a)
        ratio_line("AC2", "max ratio |M_inv - J_ForwardDiff(inverse)| / 16 eps (1 + a_max), $nm", b)
        ratio_line("AC2", "max ratio |M_inv - plain-formula inverse tangent| / 16 eps (1 + a_max), $nm", r)
    end
    for (T, seed) in ((Float64, 43001), (Double64, 43002))
        a, b, r, sy = cm_ac3_ratios(T, cm_cases(StableRNG(CM_SEED + seed), nf; maxN = T === Float64 ? 10 : 5))
        ratio_line("AC3", "max ratio |det M - 1| (256 bits) / 8 d eps ||M||_F ||M^-1||_F, $T", a)
        ratio_line("AC3", "max ratio |M' Om M - Om| (256 bits) / 8 d eps ||M||_max^2, $T", b)
        ratio_line("AC3", "max ratio of module det_residual / symplectic_residual to 64 d / 32 d times the same scales, $T", r)
        emit("AC3", "lower-left block (A or G) not symmetric to the bit: count flag, $T", sy, "0", pf(sy == 0))
    end
    for (T, cs, seed) in ((Float64, 200, 44001), (Double64, 200, 44002))
        m = sz(cs)
        r = cm_ac4_ratio(T, cm_cases(StableRNG(CM_SEED + seed), m; maxN = T === Float64 ? 10 : 5))
        ratio_line("AC4", "max ratio round-trip difference / tau(n), n in {1, n_H}, both orders, $m orbits per case, $T", r)
    end
    ratio_line("AC4", "same, Float64, N = 100 ($(length(cB)) orbits)", cm_ac4_ratio(Float64, cB))
    # the module orbit and inverse vs the independent plain formulas at 256 bits (as the package test, 20 orbits per case)
    setprecision(BigFloat, CM_BITS) do
        worst = 0.0
        for (c, z) in cm_cases(StableRNG(CM_SEED + 44003), sz(20); uniform = false)
            nH = cm_nH(c)
            sysb = cm_sys(BigFloat, c); z0 = cm_state(BigFloat, c, z)
            a = cm_orbit(sysb, z0, nH)
            b = cm_ref_orbit(c, cm_params(BigFloat, c), z0, nH)
            worst = max(worst, cm_absdiff(a[end], b[end]) / cm_tau(ldexp(1.0, 1 - CM_BITS), c, a, nH))
            worst = max(worst, cm_absdiff(inverse(sysb, a[end]), cm_ref_inverse(c, cm_params(BigFloat, c), a[end])) / cm_tau(ldexp(1.0, 1 - CM_BITS), c, a, nH))
        end
        ratio_line("AC4", "Chaoslab 256-bit orbit and inverse vs independent plain 256-bit formulas, 20 orbits per case", worst)
    end
end

# ---------------------------------------------------------------------------------------------------
# AC5, AC6: wrapped system (200 orbits per case, Float64 and Double64; counts 50 teacher-forced steps)
section("AC5-6") do
    m = sz(200)
    for (T, seed) in ((Float64, 45001), (Double64, 45002))
        cs = cm_cases(StableRNG(CM_SEED + seed), m; maxN = 5, wrapped = true)
        ok, wt, wl = cm_ac5_ratios(T, cs)
        emit("AC5", "range [0, 2 pi) and integer counts hold at every step, $T, $(length(cs)) orbits", ok, "true", pf(ok))
        ratio_line("AC5", "max ratio torus distance wrapped state vs lifted orbit / tau(k), k <= n_H, $T", wt)
        ratio_line("AC5", "max ratio |lift(state) at 256 bits - lifted orbit| / tau(k), $T", wl)
        okc = cm_ac5_counts(T, cs)
        emit("AC5", "counts equal the 256-bit increments for 50 teacher-forced steps, $T, $(length(cs)) orbits", okc, "true", pf(okc))
    end
    rng = StableRNG(CM_SEED + 45011)
    nbad = 0
    for c in (cm_y3(0.5, 0.1, 1.0e-3), cm_y4(2, 3.0, 0.2), cm_y4(5, 1.0, 0.6102; ring = true), cm_y4(100, 3.0, 0.2)), _ in 1:sz(100)
        d = cm_d(c); wsys = cm_wsys(Float64, c); zr = cm_draw(rng, d, 100.0)
        s = torus_state(wsys, cm_state(Float64, c, zr))
        ok = length(s) == 2d && all(j -> 0 <= s[j] < 2 * Float64(pi), 1:d) && all(j -> s[d + j] == round(s[d + j]), 1:d) && all(abs.(lift(wsys, s) .- zr) .<= 16 * eps(100.0))
        nbad += ok ? 0 : 1
    end
    emit("AC5", "torus_state / lift round trip |lift(torus_state(z)) - z| <= 16 eps(100), z in [-100, 100]^d, 4 systems incl. N = 100: failures", nbad, "0", pf(nbad == 0))
    for (T, seed) in ((Float64, 46001), (Double64, 46002))
        cs = cm_cases(StableRNG(CM_SEED + seed), m; maxN = 5, wrapped = true)
        counts_ok, worst = cm_ac6_ratio(T, cs)
        emit("AC6", "counts exactly equal after n forward / n inverse (both orders), n in {1, n_H}, $T, $(length(cs)) orbits", counts_ok, "true", pf(counts_ok))
        ratio_line("AC6", "max ratio torus distance round trip / tau(n), $T", worst)
    end
    res = cm_overflow_checks()
    emit("AC6", "overflow guard: 10 checks (NaN in all components one below 2^53, finite and exact two below)", string(count(res), "/", length(res)), "10/10", pf(length(res) == 10 && all(res)))
end

# ---------------------------------------------------------------------------------------------------
# AC7, AC8: K = 0
section("AC7-8") do
    for c in (cm_y3(0, 0, 0), cm_y4(2, 0, 0), cm_y4(3, 0, 0), cm_y4(100, 0, 0; ring = true))
        d1 = cm_ac7_orbit_defects(Float64, c, StableRNG(CM_SEED + 47001), sz(1000))
        d2 = cm_ac7_orbit_defects(Double64, c, StableRNG(CM_SEED + 47002), sz(200))
        emit("AC7", "orbit defects (q_j != q_j(0) + n p_j(0) or p_j != p_j(0)), 1000 steps, $(c.kind === :y3 ? "Y3" : "Y4 N=$(c.N)$(c.ring ? " ring" : "")"): Float64 $(sz(1000)) orbits / Double64 $(sz(200)) orbits", string(d1, " / ", d2), "0", pf(d1 == 0 && d2 == 0))
    end
    for (T, cs) in ((Float64, (cm_y3(0, 0, 0), cm_y4(3, 0, 0), cm_y4(5, 0, 0; ring = true), cm_y4(10, 0, 0; ring = true))), (Double64, (cm_y3(0, 0, 0), cm_y4(3, 0, 0))))
        for c in cs
            bad = sum(cm_ac7_product_defects(T, c, StableRNG(CM_SEED + 47010 + i)) for i in 1:3)
            emit("AC7", "product of n tangents != [[I, nI], [0, I]] for n in {1, 2, 10, 100, 1000, 1e4}, $(c.kind === :y3 ? "Y3" : "Y4 N=$(c.N)$(c.ring ? " ring" : "")"), $T, 3 starts", bad, "0", pf(bad == 0))
        end
    end
    ns = (1, 10, 100, 1000, 10_000, 1_000_000)
    for (T, cs) in ((Float64, (cm_y3(0, 0, 0), cm_y4(3, 0, 0), cm_y4(5, 0, 0; ring = true), cm_y4(10, 0, 0; ring = true))), (Double64, (cm_y3(0, 0, 0), cm_y4(3, 0, 0))))
        for c in cs
            sys = cm_sys(T, c); d, h = cm_d(c), cm_h(c)
            z = cm_state(T, c, cm_dyadic(StableRNG(CM_SEED + 48001), T, d))
            M0 = Matrix(tangent(sys, z)); same = true
            for _ in 1:100
                z = forward(sys, z); same &= Matrix(tangent(sys, z)) == M0
            end
            wp = zeros(T, d); wp[h + 1] = one(T); wq = zeros(T, d); wq[1] = one(T)
            exact_ok, worst, posx = true, 0.0, 0.0
            for n in ns
                P = M0^n
                exact_ok &= P == cm_shear(T, d, n)
                X = cm_expo(P, wp, n)
                setprecision(BigFloat, CM_BITS) do
                    worst = max(worst, Float64(abs(cm_big(X) - log(sqrt(BigFloat(n)^2 + 1)) / n)) / (8 * cm_eps(T)))
                end
                posx = max(posx, abs(Float64(cm_expo(P, wq, n))))
            end
            nm = "$(c.kind === :y3 ? "Y3" : "Y4 N=$(c.N)$(c.ring ? " ring" : "")"), $T"
            emit("AC8", "K=0 tangent state independent on 100 orbit points and M0^n == [[I, nI], [0, I]] for n up to 1e6, $nm", same && exact_ok, "true", pf(same && exact_ok))
            ratio_line("AC8", "max ratio |X(n) - ln(sqrt(n^2+1))/n| / 8 eps, n up to 1e6, w0 = one momentum component, $nm", worst)
            emit("AC8", "max |X(n)| for a position component (exactly 0), $nm", posx, "0", pf(posx == 0))
        end
    end
end

# ---------------------------------------------------------------------------------------------------
# AC9, AC10
section("AC9-10") do
    m = sz(200)
    for (T, seed) in ((Float64, 49001), (Double64, 49002))
        a = cm_ac9a_defects(T, m, CM_SEED + seed)
        emit("AC9", "(a) Y3 mu = 0 mapped by C vs Y1 with K = -nu, -kappa: steps with a difference, $m orbits per (nu, kappa) x 5 sets x 1000 steps, $T", a, "0", pf(a == 0))
    end
    for (T, seed) in ((Float64, 49003), (Double64, 49004))
        b = cm_ac9b_defects(T, m, CM_SEED + seed)
        emit("AC9", "(b) Y4 eps = 0, N in {2, 5, 10}, K in {1, 3}, open and ring: steps with a site differing from Y1, $m orbits per config x 1000 steps, $T", b, "0", pf(b == 0))
    end
    for (T, seed) in ((Float64, 49005), (Double64, 49006))
        c = cm_ac9c_defects(T, sz(100), CM_SEED + seed)
        emit("AC9", "(c) independence of the unperturbed pair / sites, $(sz(100)) perturbations per case x 1000 steps, $T", c, "0", pf(c == 0))
    end
    for (T, seed) in ((Float64, 50001), (Double64, 50002))
        r = cm_ac10_ratio(T, cm_cases(StableRNG(CM_SEED + seed), m; maxN = T === Float64 ? 10 : 5))
        ratio_line("AC10", "max ratio torus distance of the shifted orbit at negated on-site parameters / tau(k), k <= n_H, $m orbits per case, $T", r)
    end
    ratio_line("AC10", "same, Float64, N = 100 (20 orbits)", cm_ac10_ratio(Float64, cm_big_cases(StableRNG(CM_SEED + 50003), sz(5))))
end

# ---------------------------------------------------------------------------------------------------
# AC11, AC12
section("AC11-12") do
    for T in (Float64, Double64), p in (CM_Y3_SETS..., (-6.0, 3.0, 2.5), (4.0, 4.0, 5.0))
        sys = cm_sys(T, cm_y3(p...)); z = zero(SVector{4, T}); ok = true
        for _ in 1:10_000
            z = forward(sys, z); ok &= all(iszero, z)
        end
        emit("AC11", "origin stays for 1e4 steps, Y3 $p, $T", ok, "true", pf(ok))
    end
    vals = (-6.0, -2.0, 1.0, 3.0, 4.5, 6.0)
    for T in (Float64, Double64)
        worst_tr, worst_ev, wrong_class, off = 0.0, 0.0, 0, 0.0
        for nu in vals, ka in vals                                   # all 36 pairs (the package test uses 6)
            c = cm_y3(nu, ka, 1.0e-3)
            M = Matrix(tangent(cm_sys(T, c), zero(SVector{4, T})))
            off = max(off, Float64(maximum(abs, M[[1, 3], [2, 4]])), Float64(maximum(abs, M[[2, 4], [1, 3]])))
            for (pair, par) in (([1, 3], nu), ([2, 4], ka))
                tr = M[pair[1], pair[1]] + M[pair[2], pair[2]]
                setprecision(BigFloat, CM_BITS) do
                    trx = 2 - BigFloat(par)
                    kx, mx = cm_class(trx)
                    tol = 8 * cm_eps(T) * (1 + abs(nu) + abs(ka))
                    worst_tr = max(worst_tr, Float64(abs(cm_big(tr) - trx)) / tol)
                    kind, mag = cm_class(cm_big(tr))
                    wrong_class += kind == kx ? 0 : 1
                    worst_ev = max(worst_ev, Float64(abs(mag - mx)) / (32 * cm_eps(T) * (1 + abs(par)) * max(1, 1 / Float64(sqrt(abs(trx^2 - 4))))))
                end
            end
        end
        emit("AC11", "off-pair entries of M(0) (exactly 0), 36 (nu, kappa) pairs, $T", off, "0", pf(off == 0))
        ratio_line("AC11", "max ratio |tr block - (2 - par)| / 8 eps (1 + |nu| + |kappa|), $T", worst_tr)
        ratio_line("AC11", "max ratio |angle or ln lambda_+ - closed form| / 32 eps (1 + |par|) cond, $T", worst_ev)
        emit("AC11", "wrong classes (elliptic / hyperbolic), $T", wrong_class, "0", pf(wrong_class == 0))
    end
    for T in (Float64, Double64)
        M = Matrix(tangent(cm_sys(T, cm_y3(4, 0.1, 1.0e-3)), zero(SVector{4, T})))
        B = M[[1, 3], [1, 3]]; P = B; worst = 0.0
        for n in 2:1000
            P = B * P
            n in (10, 100, 1000) && (worst = max(worst, abs(Float64(maximum(abs, P)) / n - 4)))
        end
        emit("AC11", "nu = 4 Jordan block: |max entry of M^n / n - 4|, n in {10, 100, 1000}, $T", worst, "<= 1e-12", pf(worst <= 1.0e-12))
    end
    for (T, c) in ((Float64, cm_y4(3, 1, 0.6102; ring = true)), (Float64, cm_y4(100, 3, 0.2)), (Float64, cm_y4(100, -2, -0.5; ring = true)), (Double64, cm_y4(5, 3, -0.5)))
        sys = cm_sys(T, c); z = zeros(T, cm_d(c)); ok = true
        for _ in 1:10_000
            z = forward(sys, z); ok &= all(iszero, z)
        end
        emit("AC12", "origin stays for 1e4 steps, Y4 N=$(c.N)$(c.ring ? " ring" : ""), $T", ok, "true", pf(ok))
    end
    worst_ev, worst_off, worst_exp, wrong_class, tr_err, nmodes, nband = 0.0, 0.0, 0.0, 0, 0.0, 0, 0
    worst_off_base = 0.0
    for N in (2, 3, 5, 10, 50, 100), ring in (false, true), K in (-2.0, 0.5, 1.0, 3.0), e in (0.2, 0.6102, -0.5)
        (ring && N < 3) && continue
        c, M, G, ev, gx = cm_ac12_data(N, ring, K, e)
        epsT = cm_eps(Float64); a = 1 + cm_amax(c)
        worst_ev = max(worst_ev, maximum(abs, ev.values .- gx) / (16 * epsT * (abs(K) + 4 * abs(e))))
        U = ev.vectors; Z = zeros(N, N); UU = [U Z; Z U]; Mt = UU' * M * UU
        offmask = trues(2N, 2N)
        for k in 1:N, (i, j) in ((k, k), (k, N + k), (N + k, k), (N + k, N + k))
            offmask[i, j] = false
        end
        worst_off = max(worst_off, maximum(abs.(Mt[offmask])) / (16 * N * epsT * a))
        worst_off_base = max(worst_off_base, maximum(abs.(Mt[offmask])) / (16 * epsT * a))
        for k in 1:N
            g = ev.values[k]; trx = 2 + gx[k]; band = abs(abs(trx) - 2)
            nmodes += 1
            tr_err = max(tr_err, abs((2 + g) - trx) / (32 * epsT * a))
            if band > 64 * epsT * a
                kx, mx = cm_class(trx); kind, mag = cm_class(2 + g)
                wrong_class += kind == kx ? 0 : 1
                band > 1.0e-6 && (worst_exp = max(worst_exp, abs(mag - mx) / (32 * epsT * a * max(1.0, 1 / sqrt(abs(trx^2 - 4))))))
            else
                nband += 1
            end
        end
    end
    ratio_line("AC12", "max ratio |eig(G(0)) - (K + eps l_k)| / 16 eps (|K| + 4|eps|), 6 N x 2 boundaries x 4 K x 3 eps grid, Float64", worst_ev)
    ratio_line("AC12", "max ratio off-block entries of (U+U)' M(0) (U+U) / 16 N eps (1 + a_max) (D104a)", worst_off)
    emit("AC12", "same off-block ratio with the base bound 16 eps (1 + a_max) of the spec before D104a (informational; spec D104a: 0.02, 0.6, 20, 44 at N = 2, 10, 50, 100)", worst_off_base, "informational", "INFO")
    ratio_line("AC12", "max ratio |trace per mode - (2 + K + eps l_k)| / 32 eps (1 + a_max)", tr_err)
    ratio_line("AC12", "max ratio |ln lambda_+ or angle - closed form| / 32 eps (1 + a_max) cond, outside the parabolic band (D104b)", worst_exp)
    emit("AC12", "wrong classes outside the parabolic band (D104b); modes inside the band (trace within 64 eps (1 + a_max) of +-2, class not asserted): $nband of $nmodes", wrong_class, "0", pf(wrong_class == 0))
end

# ---------------------------------------------------------------------------------------------------
# AC13 to AC16
section("AC13-15") do
    cases = cm_cases(StableRNG(CM_SEED + 51001), sz(10); maxN = 5, uniform = false)
    wp, ws = cm_ac13_worst(cases, (5,))
    wp2, ws2 = 0.0, 0.0
    for (c, z) in cases
        a, b = cm_ac13_worst([(c, z)], (cm_nH(c),)); wp2, ws2 = max(wp2, a), max(ws2, b)
    end
    ratio_line("AC13", "max ratio |lambda_i(n) + lambda_{d+1-i}(n)| / (B(n)/n), n = 5, $(length(cases)) orbits (Y3 sets, Y4 N in {2, 3, 5})", wp)
    ratio_line("AC13", "max ratio |sum lambda_i(n)| / (B(n)/n), n = 5", ws)
    ratio_line("AC13", "same pairing ratio at n = n_H", wp2)
    ratio_line("AC13", "same sum ratio at n = n_H", ws2)
    m = sz(200)
    r = cm_ac14_ratio(Float64, m, CM_SEED + 52001)
    emit("AC14", "total momentum, Float64: max |sum p(k) - sum p(0)| / 8 N eps M_k k over $m orbits per config (N in {2,3,5,10}, eps in 4, open / ring), k <= n_H", r, "<= 1", pf(r <= 1))
    r = cm_ac14_ratio(Double64, m, CM_SEED + 52002)
    emit("AC14", "same, Double64", r, "<= 1", pf(r <= 1))
    rng = StableRNG(CM_SEED + 53001)
    worst = 0.0
    for beta in (0.1, 0.3051), K in (1.0, 3.0), _ in 1:m
        c = cm_y4(2, K, 2 * beta)
        z = cm_draw(rng, 4, 4 * Float64(pi)); nH = cm_nH(c)
        zs = cm_orbit(cm_sys(Float64, c), cm_state(Float64, c, z), nH)
        pr = cm_printed_orbit(K, beta, (z[1], z[3], z[2], z[4]) ./ (2 * pi), nH)
        for k in 1:nH
            conv = [2 * pi * pr[k + 1][1], 2 * pi * pr[k + 1][3], 2 * pi * pr[k + 1][2], 2 * pi * pr[k + 1][4]]
            worst = max(worst, cm_tdist(conv, zs[k + 1]) / cm_tau(cm_eps(Float64), c, zs, k))
        end
    end
    ratio_line("AC15", "max ratio torus distance printed unit-torus map (symplectic reading, converted) vs chain / tau(k), N = 2, beta in {0.1, 0.3051}, K in {1, 3}, $m orbits each", worst)
    # the literal reading of the printed lines has determinant 1 - dF/dq, not 1: numbers for the finding (not scored)
    let K = 3.0, beta = 0.1
        lit(u) = [u[1] + u[2], u[2] + (K / (2pi)) * sin(2pi * u[1]) - (beta / pi) * sin(2pi * (u[3] - u[1])), u[3] + u[4], u[4] + (K / (2pi)) * sin(2pi * u[3]) - (beta / pi) * sin(2pi * (u[1] - u[3]))]
        dets = [det(ForwardDiff.jacobian(lit, [0.3 + 0.7i, 0.2, 1.1 + 0.5i, 0.4])) for i in 0:2]
        say("# AC15 finding (not scored): determinant of the literal reading of the printed lines (force at old positions, position update with the old momentum) at 3 points: " * join(fmt.(dets), ", ") * " (symplectic reading: 1)")
    end
end

section("AC16") do
    m = sz(200)
    rng = StableRNG(CM_SEED + 54001)
    big100 = Tuple{CMCase, Vector{Float64}}[]
    for K in (1.0, 5.0), e in (0.2, 1.0), _ in 1:sz(5)
        push!(big100, (cm_y4(100, K, e), cm_draw(rng, 200, 4 * Float64(pi))))
    end
    cases = vcat(cm_cases(StableRNG(CM_SEED + 54002), m), big100)
    setprecision(BigFloat, CM_BITS) do
        r = cm_ac16_ratio(Float64, cases, cm_eps(Float64))
        ratio_line("AC16", "Float64 vs 256-bit at n_H: max ratio (lifted and torus difference) / tau(n_H), $(length(cases)) orbits ($m per case; N = 100: $(length(big100)) orbits)", r)
        r100 = cm_ac16_ratio(Float64, big100, cm_eps(Float64))
        ratio_line("AC16", "same restricted to the N = 100 orbits", r100)
        dd = vcat(cm_cases(StableRNG(CM_SEED + 54003), m; maxN = 5), big100)
        xs = Double64[]
        for (c, z) in dd
            for s in cm_orbit(cm_sys(Double64, c), cm_state(Double64, c, z), cm_nH(c))
                append!(xs, cm_trig_args(c, s))
            end
        end
        rng2 = StableRNG(CM_SEED + 54004)
        rand_dd(lim) = (u = (2 * rand(rng2) - 1) * lim; Double64(u) + Double64((rand(rng2) - 0.5) * eps(u)))
        xu = [rand_dd(4 * Float64(pi)) for _ in 1:sz(50_000)]
        xw = [rand_dd(1.0e3) for _ in 1:sz(20_000)]
        eO, eU, eW = cm_dd_trig_error(xs), cm_dd_trig_error(xu), cm_dd_trig_error(xw)
        say(@sprintf("# Double64 sin/cos max relative error vs BigFloat-256: orbit arguments (%d values) %.3g, uniform [-4pi,4pi] (%d) %.3g, uniform [-1e3,1e3] (%d) %.3g",
                     length(xs), eO, length(xu), eU, length(xw), eW))
        eff = max(CM_EPS_DD_FLOOR, eO, eU, eW)
        emit("AC16", "effective Double64 eps (max relative sin/cos error over the three samples; floor 2^-104)", eff, "<= 2^-100 = $(fmt(ldexp(1.0, -100)))", pf(eff <= CM_EPS_DD))
        emit("AC16", "same, log2", log2(eff), "<= -100", pf(eff <= CM_EPS_DD))
        rdd = cm_ac16_ratio(Double64, dd, eff)
        ratio_line("AC16", "Double64 vs 256-bit at n_H with eps = effective value above, $(length(dd)) orbits", rdd)
        rdd_nom = cm_ac16_ratio(Double64, dd, CM_EPS_DD_FLOOR)
        emit("AC16", "same ratio with the nominal eps 2^-104 (informational)", rdd_nom, "<= 1", pf(rdd_nom <= 1))
        rdd100 = cm_ac16_ratio(Double64, filter(t -> t[1].N == 100, dd), eff)
        ratio_line("AC16", "Double64 vs 256-bit, restricted to the N = 100 orbits", rdd100)
    end
end

# ---------------------------------------------------------------------------------------------------
# AC17: banded frame action (package part, more draws) and timing (V)
section("AC17") do
    rng = StableRNG(CM_SEED + 55001)
    worst = 0.0
    for N in (10, 50, 100), ring in (false, true), (K, e) in ((3.0, 0.6102), (5.0, 1.0)), _ in 1:sz(20)
        worst = max(worst, cm_ac17_ratios(cm_y4(N, K, e; ring), rng))
    end
    ratio_line("AC17", "max ratio |banded action of M and M^-1 on 4 vectors - dense product| / 16 eps (1 + a_max), N in {10, 50, 100}, open and ring, 20 draws per case", worst)
end

function ac17_times(::Val{N}, ring::Bool) where {N}
    c = cm_y4(N, 3.0, 0.6102; ring)
    sys = cm_sys(Float64, c)
    rng = StableRNG(CM_SEED + 55100 + N)
    z = cm_state(Float64, c, cm_draw(rng, 2N, 4 * Float64(pi)))
    W = reshape(cm_draw(rng, 2N * 4, 1.0), 2N, 4)
    out = zeros(2N, 4)
    zs = SVector{2N, Float64}(z); Ws = SMatrix{2N, 4, Float64}(W)
    ta = @belapsed tangent_action!($out, $sys, $z, $W)
    ia = @belapsed inverse_tangent_action!($out, $sys, $z, $W)
    ts = @belapsed tangent($sys, $zs) * $Ws
    td = @belapsed mul!($out, $(Matrix(tangent(sys, z))), $W)
    return (ta, ia, ts, td)
end

section("AC17-time") do
    say("# AC17 (V) timing: BenchmarkTools minimum over samples, frame of 4 vectors, 5 rounds interleaved N = 50 / 100; timer resolution on this host is 100 ns (see the standard-map notes)")
    for ring in (false, true)
        for N in (10, 20, 50, 100)                              # warm-up and compile
            ac17_times(Val(N), ring)
        end
        res = Dict(N => [Inf, Inf, Inf, Inf] for N in (10, 20, 50, 100))
        for _ in 1:5, N in (10, 20, 50, 100)
            res[N] = min.(res[N], collect(ac17_times(Val(N), ring)))
        end
        names = ("tangent_action! (vector state)", "inverse_tangent_action!", "tangent(sys, SVector) * SMatrix (fixed-size, band built per call)", "dense M * W (mul!, tangent excluded; quadratic, informational)")
        for (i, nm) in enumerate(names)
            r = res[100][i] / res[50][i]
            line = "N=10 $(fmt(res[10][i] * 1e9)) ns, N=20 $(fmt(res[20][i] * 1e9)) ns, N=50 $(fmt(res[50][i] * 1e9)) ns, N=100 $(fmt(res[100][i] * 1e9)) ns"
            i <= 3 ? emit("AC17", "time ratio N=100 / N=50, $nm, $(ring ? "ring" : "open"): $line", r, "<= 2.6 (expected 2)", pf(r <= 2.6)) :
                     emit("AC17", "time ratio N=100 / N=50, $nm, $(ring ? "ring" : "open"): $line", r, "informational (expected 4)", "INFO")
        end
    end
end

# ---------------------------------------------------------------------------------------------------
# D105: fixed-size Y4 states, full parameter grid of criterion 1 (25 (K, eps) pairs), plus the stream at 2000 orbits
section("D105") do
    Ks, Es = (0.0, 0.5, 1.0, 3.0, 5.0), (0.0, 0.2, 0.6102, 1.0, -0.5)
    rng = StableRNG(CM_SEED + 56001)
    for T in (Float64, Double64)
        bad, ncase = 0, 0
        for N in (2, 3, 5, 10, 100), ring in (false, true)
            ring && N < 3 && continue
            for K in Ks, e in Es, _ in 1:(N == 100 ? 1 : 2)
                if N == 100 && !((K, e) in ((3.0, 0.6102), (5.0, 1.0), (0.0, 0.0), (1.0, -0.5)))
                    continue
                end
                bad += cm105_mismatches(cm_y4(N, K, e; ring), T, rng); ncase += 1
            end
        end
        emit("D105a", "fixed-size vs vector-state results not bitwise equal (forward, inverse, tangent and inverse tangent entries, frame products vs banded actions, torus_state, lift; lifted and wrapped), $T, N in {2,3,5,10,100} open / ring, 25 (K, eps) pairs (N=100: 4), $ncase cases", bad, "0", pf(bad == 0))
    end
    for N in (2, 5), wrapped in (false, true)
        c = cm_y4(N, 3.0, 0.2)
        sys = wrapped ? cm_wsys(Float64, c) : cm_sys(Float64, c)
        n, d = cm_nH(c), 2N
        X0 = [cm105_sv(wrapped ? cm105_wstate(Float64, rng, d) : cm_draw(rng, d, 4 * Float64(pi))) for _ in 1:2000]
        W0 = [cm105_sm(cm105_frame(Float64, rng, d)) for _ in 1:2000]
        res = propagate(sys, X0, n; k = 4, tau = n, seed = CM_SEED, W0, backend = :cpu)
        ref = [cm105_serial(sys, X0[j], W0[j], n) for j in eachindex(X0)]
        bx = count(j -> !isequal(res.x[j], ref[j][1]), eachindex(X0)); bw = count(j -> !isequal(res.W[j], ref[j][2]), eachindex(X0))
        typed = tangent(sys, X0[1]) isa ChainTangent && all(r -> r[1] isa SVector && r[2] isa SMatrix, ref)
        emit("D105b", "stream (CPU backend, $(Threads.nthreads()) threads) vs plain serial loop, N=$N $(wrapped ? "wrapped" : "lifted"), 2000 orbits, k = 4, n = $n: orbits with different final state / frame; fixed-size types kept", string(bx, " / ", bw, "; ", typed), "0 / 0; true", pf(bx == 0 && bw == 0 && typed))
    end
    for T in (Float64, Double64), N in (2, 10, 100)
        c = cm_y4(N, 3.0, 0.2); d = 2N
        Ws = cm105_sm(cm105_frame(T, rng, d))
        a = cm105_alloc(cm_sys(T, c), cm105_sv(cm_state(T, c, cm_draw(rng, d, 4.0))), Ws)
        b = cm105_alloc(cm_wsys(T, c), cm105_sv(cm105_wstate(T, rng, d)), Ws)
        emit("D105c", "bytes allocated by forward, inverse, tangent * W, inverse_tangent * W after warm-up, N=$N, $T: lifted / wrapped; isbits system values", string(a, " / ", b, "; ", isbitstype(typeof(cm_sys(T, c))), "/", isbitstype(typeof(cm_wsys(T, c)))), "(0,0,0,0) / (0,0,0,0); true/true", pf(all(==(0), a) && all(==(0), b) && isbitstype(typeof(cm_sys(T, c))) && isbitstype(typeof(cm_wsys(T, c)))))
    end
end

# ---------------------------------------------------------------------------------------------------
# AC18 (CPU parts) and AC19 (re-run) ; the GPU parts are in verify/gpu.jl
function plain_frames(sys, Z, W0, n)
    out = similar(W0)
    for j in eachindex(Z)
        z = Z[j]; W = W0[j]
        for _ in 1:n
            W = tangent(sys, z) * W; z = forward(sys, z)
        end
        out[j] = W
    end
    return out
end

section("AC18") do
    rng = StableRNG(CM_SEED + 57001)
    lim = 4 * Float64(pi)
    for p in ((0.5, 0.1, 1.0e-3), (2.0, 0.5, 0.1))
        c = cm_y3(p...); sys = cm_sys(Float64, c); n = cm_nH(c)
        Z = [SVector{4, Float64}(ntuple(_ -> (2 * rand(rng) - 1) * lim, 4)) for _ in 1:sz(100_000)]
        W0 = fill(one(SMatrix{4, 4, Float64}), length(Z))
        plain = cm_run_serial(sys, Z, n); thr = cm_run_threaded(sys, Z, n)
        t = @elapsed (res = propagate(sys, Z, n; k = 4, tau = 1, seed = CM_SEED, W0 = W0, backend = :cpu))
        res2 = propagate(sys, Z, n; k = 4, tau = 1, seed = CM_SEED, W0 = W0, backend = :cpu, chunk = 777)
        Wp = plain_frames(sys, Z, W0, n)
        emit("AC18", "Y3 $p n=$n: $(length(Z)) orbits, threaded forward loop ($(Threads.nthreads()) threads) vs serial loop, bitwise", isequal(plain, thr), "true", pf(isequal(plain, thr)))
        emit("AC18", "Y3 $p: propagate (CPU kernel) x and W vs plain loops, bitwise (wall $(fmt(t)) s)", isequal(res.x, plain) && isequal(res.W, Wp), "true", pf(isequal(res.x, plain) && isequal(res.W, Wp)))
        emit("AC18", "Y3 $p: propagate with chunk 777 vs default chunk, bitwise", isequal(res.x, res2.x) && isequal(res.W, res2.W), "true", pf(isequal(res.x, res2.x) && isequal(res.W, res2.W)))
        Zw = [SVector{8, Float64}(ntuple(i -> i <= 4 ? 2 * Float64(pi) * rand(rng) : 0.0, 8)) for _ in 1:sz(10_000)]
        wsys = cm_wsys(Float64, c)
        pw = cm_run_serial(wsys, Zw, n); rw = propagate(wsys, Zw, n; k = 4, tau = 1, seed = CM_SEED, W0 = fill(one(SMatrix{4, 4, Float64}), length(Zw)), backend = :cpu)
        emit("AC18", "wrapped Y3 $p: $(length(Zw)) orbits threaded vs serial; propagate x vs plain loop, bitwise", isequal(pw, cm_run_threaded(wsys, Zw, n)) && isequal(rw.x, pw), "true", pf(isequal(pw, cm_run_threaded(wsys, Zw, n)) && isequal(rw.x, pw)))
    end
    for N in (2, 5)
        c = cm_y4(N, 3.0, 0.2); sys = cm_sys(Float64, c); n = cm_nH(c); d = 2N
        Zv = [cm_draw(rng, d, lim) for _ in 1:sz(10_000)]
        Zs = [SVector{d, Float64}(z) for z in Zv]
        W0 = fill(SMatrix{d, 4, Float64}(ntuple(i -> Float64(((i - 1) % d) == ((i - 1) ÷ d)), d * 4)), length(Zs))
        a = isequal(cm_run_serial(sys, Zv, n), cm_run_threaded(sys, Zv, n))
        plain = cm_run_serial(sys, Zs, n)
        res = propagate(sys, Zs, n; k = 4, tau = 1, seed = CM_SEED, W0 = W0, backend = :cpu)
        Wp = plain_frames(sys, Zs, W0, n)
        emit("AC18", "Y4 N=$N n=$n: $(length(Zv)) orbits, Vector states threaded vs serial; SVector stream x and W vs plain loops, bitwise; SVector forward vs Vector forward", string(a, "; ", isequal(res.x, plain) && isequal(res.W, Wp), "; ", isequal(plain, [SVector{d}(z) for z in cm_run_serial(sys, Zv, n)])), "true; true; true", pf(a && isequal(res.x, plain) && isequal(res.W, Wp) && isequal(plain, [SVector{d}(z) for z in cm_run_serial(sys, Zv, n)])))
    end
    emit("AC18", "GPU vs CPU (Y3 1e5, Y4 N in {2, 5} 1e4, Double64 gate): see verify/gpu.jl, log gpu-coupled-maps-verify-log.txt", "-", "tau(n)", "STAGE04-GPU")
end

# 1 vs 8 threads (D104c): per-orbit hashes of the results, written by each run, compared by the "compare" section
section("threads") do
    mkpath(DUMP)
    res = Dict{String, Vector{UInt64}}()
    rng = StableRNG(CM_SEED + 58001)
    lim = 4 * Float64(pi)
    hs(x) = hash(x.data)
    for p in ((0.5, 0.1, 1.0e-3), (2.0, 0.5, 0.1))
        c = cm_y3(p...); sys = cm_sys(Float64, c); n = cm_nH(c)
        Z = [SVector{4, Float64}(ntuple(_ -> (2 * rand(rng) - 1) * lim, 4)) for _ in 1:sz(100_000)]
        res["Y3 $p threaded forward loop"] = UInt64[hs(z) for z in cm_run_threaded(sys, Z, n)]
        r = propagate(sys, Z, n; k = 4, tau = 1, seed = CM_SEED, W0 = fill(one(SMatrix{4, 4, Float64}), length(Z)), backend = :cpu)
        res["Y3 $p propagate x and W"] = UInt64[hash((r.x[j].data, r.W[j].data)) for j in eachindex(Z)]
        r2 = propagate(sys, Z, n; k = 4, tau = 1, seed = CM_SEED, backend = :cpu)         # random initial frames from the stream keyed by the orbit index
        res["Y3 $p propagate x and W, random W0 (seed stream)"] = UInt64[hash((r2.x[j].data, r2.W[j].data)) for j in eachindex(Z)]
        Zw = [SVector{8, Float64}(ntuple(i -> i <= 4 ? 2 * Float64(pi) * rand(rng) : 0.0, 8)) for _ in 1:sz(10_000)]
        rw = propagate(cm_wsys(Float64, c), Zw, n; k = 4, tau = 1, seed = CM_SEED, backend = :cpu)
        res["wrapped Y3 $p propagate x and W"] = UInt64[hash((rw.x[j].data, rw.W[j].data)) for j in eachindex(Zw)]
    end
    for N in (2, 5), wrapped in (false, true)
        c = cm_y4(N, 3.0, 0.2); d = 2N; sys = wrapped ? cm_wsys(Float64, c) : cm_sys(Float64, c); n = cm_nH(c)
        Zs = [SVector{wrapped ? 2d : d, Float64}(wrapped ? cm105_wstate(Float64, rng, d) : cm_draw(rng, d, lim)) for _ in 1:sz(10_000)]
        r = propagate(sys, Zs, n; k = 4, tau = 1, seed = CM_SEED, backend = :cpu)
        res["Y4 N=$N $(wrapped ? "wrapped" : "lifted") stream x and W (random W0)"] = UInt64[hash((r.x[j].data, r.W[j].data)) for j in eachindex(Zs)]
    end
    file = joinpath(DUMP, "cm-threads-t$(Threads.nthreads()).bin")
    serialize(file, res)
    for k in sort(collect(keys(res)))
        say("# threads=$(Threads.nthreads()) | $k | $(length(res[k])) orbits | digest $(hash(res[k]))")
    end
end

"compare" in ONLY && section("compare") do
    f1, f8 = joinpath(DUMP, "cm-threads-t1.bin"), joinpath(DUMP, "cm-threads-t8.bin")
    (isfile(f1) && isfile(f8)) || error("missing dump: $f1 or $f8")
    a, b = deserialize(f1), deserialize(f8)
    for k in sort(collect(keys(a)))
        nd = count(j -> a[k][j] != b[k][j], eachindex(a[k]))
        emit("AC18", "1 vs 8 threads, per-orbit result hash (x and W bits), $k, $(length(a[k])) orbits: orbits that differ", nd, "0", pf(nd == 0 && length(a[k]) == length(b[k])))
    end
end

# ---------------------------------------------------------------------------------------------------
# O9: ForwardDiff of the n-step map vs the product of analytic tangents, n = n_H, 10^4 seeded Y3 points and Y4 N in {2, 5}
function o9_point(c::CMCase, z)
    n = cm_nH(c); sys = cm_sys(Float64, c); d = cm_d(c)
    z0 = cm_state(Float64, c, z)
    P = Matrix{Float64}(I, d, d); zz = z0; fro = 1.0
    pts = Vector{typeof(z0)}(undef, n)
    for k in 1:n
        pts[k] = zz
        M = Matrix(tangent(sys, zz)); P = M * P; fro *= norm(M); zz = forward(sys, zz)
    end
    J = Matrix(ForwardDiff.jacobian(v -> (w = v; for _ in 1:n; w = forward(sys, w); end; w), z0))
    bound = 8 * n * ldexp(1.0, -52) * fro * sqrt(Float64(d))
    rB, rJB = setprecision(BigFloat, CM_BITS) do
        sysB = cm_sys(BigFloat, c)
        Pb = Matrix{BigFloat}(I, d, d)
        for k in 1:n
            Pb = Matrix(tangent(sysB, cm_state(BigFloat, c, pts[k]))) * Pb
        end
        (Float64(maximum(abs(BigFloat(P[i, j]) - Pb[i, j]) for i in 1:d, j in 1:d)) / bound, Float64(maximum(abs(BigFloat(J[i, j]) - Pb[i, j]) for i in 1:d, j in 1:d)) / bound)
    end
    return (Float64(maximum(abs, J - P)) / bound, rB, rJB)
end

section("O9") do
    rng = StableRNG(CM_SEED + 59001)
    cs = cm_cases(rng, sz(1667); y4 = false)                          # 6 groups x 1667 = 10002 Y3 points: 5 sets + uniform
    c4 = [(c, z) for (c, z) in cm_cases(rng, sz(1000); y3 = false) if c.N in (2, 5)]
    for (nm, set) in (("Y3", cs), ("Y4 N in {2, 5} (4 configs, $(sz(1000)) points each)", c4))
        res = Vector{NTuple{3, Float64}}(undef, length(set))
        Threads.@threads for i in eachindex(set)
            res[i] = o9_point(set[i]...)
        end
        say("# O9 sample $nm: $(length(set)) points, n = n_H, z0 uniform in [-4 pi, 4 pi]^d; bound = 8 n eps prod_k ||M_k||_F ||I||_F")
        ratio_line("O9", "$nm: max ratio |J_ForwardDiff(n-step) - prod M_k| / bound", maximum(first, res))
        ratio_line("O9", "$nm: max ratio |prod M_k (Float64) - product at 256 bits, same realised points| / bound", maximum(r -> r[2], res))
        ratio_line("O9", "$nm: max ratio |J_ForwardDiff - product at 256 bits| / bound", maximum(last, res))
    end
end

# ---------------------------------------------------------------------------------------------------
# AC20 (orbits C and D, published behaviour), O7 (independent QR spectrum), O6b (ChaosTools spectra), figure
const Y3PUB = (0.5, 0.1, 1.0e-3)
const ORB_C = SVector(0.5, 0.5, 0.0, 0.0)          # item 1 labelling, positions first (item 2 labelling: (0.5, 0, 0.5, 0))
const ORB_D = SVector(3.0, 0.5, 0.0, 0.0)          # (item 2: (3, 0, 0.5, 0))
const W_SALI = SMatrix{4, 2, Float64}(1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0)   # deviation vectors (1,1,1,1) and (1,0,0,0)

# One pass along the orbit: SALI frame (two vectors, renormalised every step), finite-time exponent L_N from the first of them, and a 4-vector QR spectrum
# (modified Gram-Schmidt, identity start). Returns the series at the checkpoints ns and the SALI of every step in [lo, hi] summarised.
function orbit_series(sys, x0, ns, nmax; salirange = (10_000, 100_000))
    x = x0; V = W_SALI; Q = one(SMatrix{4, 4, Float64})
    s1 = 0.0; s2 = 0.0; S = zeros(4)
    sali = zeros(length(ns)); L = zeros(length(ns)); lam = zeros(4, length(ns))
    ci = 1
    salisum = 0.0; saliN = 0; salimin = Inf; salimax = 0.0
    first12 = 0; sali_at = Dict{Int, Float64}()
    for n in 1:nmax
        M = tangent(sys, x)
        x = forward(sys, x)
        V = M * V
        a1 = norm(V[:, 1]); a2 = norm(V[:, 2])
        u1, u2 = V[:, 1] / a1, V[:, 2] / a2
        s1 += log(a1); s2 += log(a2)
        V = hcat(u1, u2)
        sa = min(norm(u1 + u2), norm(u1 - u2))
        # modified Gram-Schmidt on M * Q
        B = M * Q
        c1 = B[:, 1]; r1 = norm(c1); q1 = c1 / r1
        c2 = B[:, 2] - dot(q1, B[:, 2]) * q1; r2 = norm(c2); q2 = c2 / r2
        c3 = B[:, 3] - dot(q1, B[:, 3]) * q1; c3 = c3 - dot(q2, c3) * q2; r3 = norm(c3); q3 = c3 / r3
        c4 = B[:, 4] - dot(q1, B[:, 4]) * q1; c4 = c4 - dot(q2, c4) * q2; c4 = c4 - dot(q3, c4) * q3; r4 = norm(c4); q4 = c4 / r4
        Q = hcat(q1, q2, q3, q4)
        S[1] += log(r1); S[2] += log(r2); S[3] += log(r3); S[4] += log(r4)
        if salirange[1] <= n <= salirange[2]
            salisum += sa; saliN += 1; salimin = min(salimin, sa); salimax = max(salimax, sa)
        end
        first12 == 0 && sa < 1.0e-12 && (first12 = n)
        n in (1000, 5000, 10_000, 20_000, 100_000, 1_000_000) && (sali_at[n] = sa)
        if ci <= length(ns) && n == ns[ci]
            sali[ci] = sa; L[ci] = s1 / n; lam[:, ci] = S ./ n; ci += 1
        end
    end
    return (; sali, L, lam, salimean = saliN == 0 ? NaN : salisum / saliN, salimin, salimax, first12, sali_at, s1, s2)
end

# Independent Y3 in the labelling of Skokos2003a Eq. (9) (pair ordering), written without Chaoslab; plain Matrix and LinearAlgebra.qr
function indep_spectrum(p, u0, N)
    nu, ka, mu = p
    u = collect(u0)
    Q = Matrix{Float64}(I, 4, 4)[:, PERM]; S = zeros(4)   # e_j of the module labelling is e_PERM[j] here: same Gram-Schmidt order as the module run
    for _ in 1:N
        S1 = u[1] + u[2]; S2 = u[3] + u[4]; T = u[1] + u[2] + u[3] + u[4]
        sT, cT = sin(T), cos(T)
        a = -nu * cos(S1) - mu * sT; b = -mu * sT; cc = -ka * cos(S2) - mu * sT
        J = [1.0 1.0 0.0 0.0; a 1 + a b b; 0.0 0.0 1.0 1.0; b b cc 1 + cc]
        u = [S1, u[2] - nu * sin(S1) - mu * (1 - cT), S2, u[4] - ka * sin(S2) - mu * (1 - cT)]
        F = qr(J * Q)
        Q = Matrix(F.Q); R = F.R
        for j in 1:4
            S[j] += log(abs(R[j, j]))
        end
    end
    return S ./ N
end
const PERM = [1, 3, 2, 4]
const WPERM = SMatrix{4, 4, Float64}(Matrix{Float64}(I, 4, 4)[:, PERM])    # item 2 labelling (x1, x2 | x3, x4 = pair 1 | pair 2) <-> item 1 labelling (positions first): an involution

ct3_rule(u, p, n) = (nu = p[1]; ka = p[2]; mu = p[3]; s1 = u[1] + u[2]; s2 = u[3] + u[4]; c = mu * (1 - cos(u[1] + u[2] + u[3] + u[4]));
                     SVector(s1, u[2] - nu * sin(s1) - c, s2, u[4] - ka * sin(s2) - c))
function ct3_jac(u, p, n)
    nu, ka, mu = p[1], p[2], p[3]
    S1 = u[1] + u[2]; S2 = u[3] + u[4]; sT = sin(u[1] + u[2] + u[3] + u[4])
    a = -nu * cos(S1) - mu * sT; b = -mu * sT; cc = -ka * cos(S2) - mu * sT
    return @SMatrix [1.0 1.0 0.0 0.0; a 1+a b b; 0.0 0.0 1.0 1.0; b b cc 1+cc]
end
function ct4_rule(::Val{N}) where {N}
    return function (u, p, n)
        K, e = p[1], p[2]
        q = ntuple(j -> u[j], Val(N)); pm = ntuple(j -> u[N + j], Val(N))
        F = ntuple(Val(N)) do j
            f = K * sin(q[j])
            j > 1 && (f += e * sin(q[j] - q[j - 1]))
            j < N && (f += e * sin(q[j] - q[j + 1]))
            f
        end
        pn = ntuple(j -> pm[j] + F[j], Val(N)); qn = ntuple(j -> q[j] + pn[j], Val(N))
        return SVector{2N, Float64}(qn..., pn...)
    end
end
function ct4_jac(::Val{N}) where {N}
    return function (u, p, n)
        K, e = p[1], p[2]
        G = zeros(MMatrix{N, N, Float64})
        for j in 1:N
            G[j, j] = K * cos(u[j])
            if j > 1
                G[j, j] += e * cos(u[j] - u[j - 1]); G[j, j - 1] = -e * cos(u[j] - u[j - 1])
            end
            if j < N
                G[j, j] += e * cos(u[j] - u[j + 1]); G[j, j + 1] = -e * cos(u[j] - u[j + 1])
            end
        end
        Gs = SMatrix{N, N, Float64}(G); Iv = one(SMatrix{N, N, Float64}); Z = zero(SMatrix{N, N, Float64})
        return SMatrix{2N, 2N, Float64}(vcat(hcat(Iv + Gs, Iv), hcat(Gs, Iv)))
    end
end

section("AC20") do
    sys = FroeschleMap(Y3PUB...)
    ns = sort(unique(vcat(round.(Int, 10 .^ range(0, 6, length = 241)), [1000, 5000, 10_000, 20_000, 100_000, 1_000_000])))
    serC = orbit_series(sys, ORB_C, ns, 1_000_000)
    serD = orbit_series(sys, ORB_D, ns, 1_000_000)
    DATA = (; ns, serC, serD)
    Serialization.serialize(joinpath(tempdir(), "cm-ac20-series.bin"), DATA)
    say("# AC20: parameters (nu, kappa, mu) = $Y3PUB, item 1 form; C = $(Vector(ORB_C)) (item 2 labelling (0.5, 0, 0.5, 0)), D = $(Vector(ORB_D)) (item 2 labelling (3, 0, 0.5, 0)); deviation vectors (1,1,1,1) and (1,0,0,0); binary64; one pass of 1e6 steps each")
    say(@sprintf("# AC20 SALI of C over n in [1e4, 1e5]: mean %.4f, min %.4f, max %.4f (reported 0.28; drafting check 0.14 to 0.39 for single times)", serC.salimean, serC.salimin, serC.salimax))
    emit("AC20", "(a) time average of SALI of C over n in [1e4, 1e5]", serC.salimean, "in [0.1, 0.5] (reported 0.28)", pf(0.1 <= serC.salimean <= 0.5))
    emit("AC20", "(b) SALI of D: first n with SALI < 1e-12 (reported about 4.7e3); SALI at n = 1e3 / 5e3 / 1e4 / 2e4: $(join(fmt.([serD.sali_at[k] for k in (1000, 5000, 10_000, 20_000)]), " / "))", serD.first12, "some n in [1e3, 2e4]", pf(1000 <= serD.first12 <= 20_000))
    LD, LC = serD.L[end], serC.L[end]
    emit("AC20", "(c) finite-time maximal exponent L_N of D at n = 1e6 (first deviation vector (1,1,1,1)); Skokos2010 Fig. 7 about 8e-3, Skokos2003a text 5e-2, Fig. 2b about 1e-2", LD, "in [5e-3, 1.5e-2]", pf(5.0e-3 <= LD <= 1.5e-2))
    emit("AC20", "(d) finite-time maximal exponent L_N of C at n = 1e6 (reported 1.6e-6 at 1e7); decay law ~ ln(n)/n", LC, "<= 1e-4", pf(LC <= 1.0e-4))
    i5, i6 = findfirst(==(100_000), ns), findfirst(==(1_000_000), ns)
    emit("AC20", "info: L_N of C at n = 1e5 / 1e6 (1/n law: ratio should be about 10 to 14): $(fmt(serC.L[i5])) / $(fmt(serC.L[i6]))", serC.L[i5] / serC.L[i6], "info", "INFO")
    emit("AC20", "info: L_N of D at n = 1e3 / 5e3 / 1e4 / 2e4 / 1e5 / 1e6 (drafting check 0.0150, 0.0090, 0.0085, 0.0083, 0.0093): $(join(fmt.([serD.L[findfirst(==(k), ns)] for k in (1000, 5000, 10_000, 20_000, 100_000, 1_000_000)]), " / "))", LD, "info", "INFO")
    for (nm, ser) in (("C", serC), ("D", serD)), (i, n) in ((i5, 100_000), (i6, 1_000_000))
        lam = ser.lam[:, i]
        say(@sprintf("# AC20 QR spectrum of %s at n = %d: (%.5g, %.5g, %.5g, %.5g); pairing |l1+l4| = %.3g, |l2+l3| = %.3g; sum = %.3g (QR exponents; the exact pairing is for singular values, criterion 13)",
                     nm, n, lam[1], lam[2], lam[3], lam[4], abs(lam[1] + lam[4]), abs(lam[2] + lam[3]), sum(lam)))
    end
    # cross-check of the hand loop against the package estimators on the same orbit (SALI reducer, Benettin, QR spectrum through propagate)
    for (nm, x0, ser) in (("C", ORB_C, serC), ("D", ORB_D, serD)), n in (10_000, 100_000, 1_000_000)
        i = findfirst(==(n), ns)
        rb = propagate(sys, [x0], n; k = 1, tau = 1, seed = CM_SEED, reducers = (Benettin(),), W0 = [SMatrix{4, 1, Float64}(1.0, 1.0, 1.0, 1.0)], backend = :cpu)
        rs = propagate(sys, [x0], n; k = 2, tau = 1, seed = CM_SEED, reducers = (SALI(),), W0 = [W_SALI], backend = :cpu)
        rq = propagate(sys, [x0], n; k = 4, tau = 1, seed = CM_SEED, reducers = (QRSpectrum(HouseholderQR()),), W0 = [one(SMatrix{4, 4, Float64})], backend = :cpu)
        lb = rb.finals[1][1].lambda[1]; sa = rs.finals[1][1].sali; lq = rq.finals[1][1].lambda
        dB = abs(lb - ser.L[i]); dS = abs(sa - ser.sali[i]); dQ = maximum(abs.(collect(lq) .- ser.lam[:, i]))
        emit("AC20", "cross-check at n = $n, orbit $nm: |Benettin L_N - hand loop| / |SALI reducer - hand loop| / max |QRSpectrum(HouseholderQR) - hand-loop MGS spectrum|: $(fmt(dB)) / $(fmt(dS)) / $(fmt(dQ))", max(dB, dQ), "<= 1e-9 (SALI: informational)", pf(dB <= 1e-9 && dQ <= 1e-9))
    end
end

section("O7-O6b") do
    sys = FroeschleMap(Y3PUB...)
    # O7: independent re-implementation (item 2 labelling, plain Matrix, LinearAlgebra.qr) vs Chaoslab QRSpectrum, same frame (identity of the module labelling,
    # i.e. the permuted identity in the independent labelling). Short horizons and the regular orbit C are scored; the chaotic orbit D at long n is informational
    # because two orbits that differ by rounding separate exponentially, so its finite-time exponents differ by sampling, not by formula.
    for (nm, x0) in (("C", ORB_C), ("D", ORB_D)), n in (100, 1000, 100_000, 1_000_000)
        li = indep_spectrum(Y3PUB, x0[PERM], n)
        rq = propagate(sys, [x0], n; k = 4, tau = 1, seed = CM_SEED, reducers = (QRSpectrum(HouseholderQR()),), W0 = [one(SMatrix{4, 4, Float64})], backend = :cpu)
        d = maximum(abs.(li .- collect(rq.finals[1][1].lambda)))
        scored = nm == "C" || n <= 1000
        emit("O7", "orbit $nm, n = $n: independent re-implementation (labelling of Skokos2003a Eq. 9, qr) spectrum $(join(fmt.(li), ", ")) vs Chaoslab QRSpectrum; max difference", d, scored ? "<= 1e-9" : "informational (chaotic, long n)", scored ? pf(d <= 1.0e-9) : "INFO")
    end
    rng = StableRNG(CM_SEED + 60001)
    worst = 0.0; worst_long = 0.0
    for _ in 1:sz(10)
        x0 = SVector{4, Float64}(ntuple(_ -> (2 * rand(rng) - 1) * 4 * Float64(pi), 4))
        for (n, isshort) in ((200, true), (20_000, false))
            li = indep_spectrum(Y3PUB, x0[PERM], n)
            rq = propagate(sys, [x0], n; k = 4, tau = 1, seed = CM_SEED, reducers = (QRSpectrum(HouseholderQR()),), W0 = [one(SMatrix{4, 4, Float64})], backend = :cpu)
            d = maximum(abs.(li .- collect(rq.finals[1][1].lambda)))
            isshort ? (worst = max(worst, d)) : (worst_long = max(worst_long, d))
        end
    end
    emit("O7", "$(sz(10)) seeded random starts, n = 200: max difference of the spectra, independent vs Chaoslab", worst, "<= 1e-9", pf(worst <= 1.0e-9))
    emit("O7", "same starts, n = 2e4 (chaotic starts separate by rounding growth: informational)", worst_long, "informational", "INFO")
    # O6b: ChaosTools lyapunovspectrum with the hand-written discrete system (item 2 labelling), loose
    tds = CT.TangentDynamicalSystem(CT.DeterministicIteratedMap(ct3_rule, SVector(0.5, 0.0, 0.5, 0.0), collect(Y3PUB)); J = ct3_jac, k = 4)
    for (nm, x0) in (("C", ORB_C), ("D", ORB_D)), n in (100_000, 1_000_000)
        lc = collect(CT.lyapunovspectrum(tds, n; u0 = SVector{4, Float64}(x0[PERM])))
        rq = propagate(sys, [x0], n; k = 4, tau = 1, seed = CM_SEED, reducers = (QRSpectrum(HouseholderQR()),), W0 = [WPERM], backend = :cpu)
        lq = collect(rq.finals[1][1].lambda)
        say(@sprintf("# O6b ChaosTools lyapunovspectrum, orbit %s, n = %d: (%.5g, %.5g, %.5g, %.5g); pairing |l1+l4| = %.3g, |l2+l3| = %.3g; sum = %.3g", nm, n, lc[1], lc[2], lc[3], lc[4], abs(lc[1] + lc[4]), abs(lc[2] + lc[3]), sum(lc)))
        d = maximum(abs.(lc .- lq))
        emit("O6b", "orbit $nm, n = $n: max |ChaosTools - Chaoslab QRSpectrum| over the four exponents (loose, K18; same frame and exponent order; the regular orbit C agrees to rounding, the chaotic orbit D differs by the separation of the two rounding-perturbed orbits)", d, "<= 1e-3", pf(d <= 1.0e-3))
    end
    # Y4 N = 2 chain check of the spec drafting check (K = 3, eps = 0.2, one random orbit, 2e5 steps): (0.749, 0.461, -0.461, -0.749)
    tdc = CT.TangentDynamicalSystem(CT.DeterministicIteratedMap(ct4_rule(Val(2)), SVector(1.0, 2.0, 0.5, -0.5), [3.0, 0.2]); J = ct4_jac(Val(2)), k = 4)
    lc = collect(CT.lyapunovspectrum(tdc, 200_000; u0 = SVector{4, Float64}(1.0, 2.0, 0.5, -0.5)))
    csys = StandardMapChain(2, 3.0, 0.2)
    rq = propagate(csys, [SVector{4, Float64}(1.0, 2.0, 0.5, -0.5)], 200_000; k = 4, tau = 1, seed = CM_SEED, reducers = (QRSpectrum(HouseholderQR()),), W0 = [one(SMatrix{4, 4, Float64})], backend = :cpu)
    lq = collect(rq.finals[1][1].lambda)
    say(@sprintf("# O6b chain N = 2, K = 3, eps = 0.2, start (1, 2, 0.5, -0.5), n = 2e5: ChaosTools (%.4g, %.4g, %.4g, %.4g), Chaoslab QR via the stream with fixed-size states (%.4g, %.4g, %.4g, %.4g); sum %.2g / %.2g; spec drafting check (0.749, 0.461, -0.461, -0.749) for one other orbit",
                 lc..., lq..., sum(lc), sum(lq)))
    d = maximum(abs.(sort(lc; rev = true) .- sort(lq; rev = true)))
    emit("O6b", "chain N = 2 spectrum: max |ChaosTools - Chaoslab stream QRSpectrum|", d, "<= 1e-3", pf(d <= 1.0e-3))
end

# ---------------------------------------------------------------------------------------------------
# AC21: ChaosTools hand-written discrete systems vs this module, 100 seeded orbits per system
function ct_compare(label, c::CMCase, mk_tds, tostate, fromQ, jac, par, m, seed)
    rng = StableRNG(seed)
    sys = cm_sys(Float64, c); d = cm_d(c); nH = cm_nH(c)
    worst_orb, worst_prod, worst_tf = 0.0, 0.0, 0.0
    Id = one(SMatrix{d, d, Float64})
    tds = mk_tds()
    for _ in 1:m
        z = cm_draw(rng, d, 4 * Float64(pi))
        z0 = c.kind === :y3 ? SVector{4, Float64}(ntuple(i -> z[i], 4)) : SVector{d, Float64}(ntuple(i -> z[i], d))
        CT.reinit!(tds, tostate(z0); Q0 = Id)
        zs = cm_orbit(sys, z0, nH)
        P = Matrix{Float64}(I, d, d); fro = 1.0
        Ptf = Matrix{Float64}(I, d, d)
        for k in 1:nH
            M = Matrix(tangent(sys, zs[k])); P = M * P; fro *= norm(M)
            CT.step!(tds)
            zo = CT.current_state(tds); Q = fromQ(Matrix(CT.current_deviations(tds)))
            worst_orb = max(worst_orb, cm_tdist(fromQ === identity ? zo : zo[PERM], zs[k + 1]) / cm_tau(cm_eps(Float64), c, zs, k))
            worst_prod = max(worst_prod, maximum(abs, Q - P) / (8 * k * cm_eps(Float64) * fro))
        end
        # teacher-forced: the oracle Jacobian evaluated at the module's own orbit points (isolates the Jacobian conversion from the orbit divergence)
        Pt = Matrix{Float64}(I, d, d); frot = 1.0
        for k in 1:nH
            u = tostate(zs[k])
            Mo = fromQ(Matrix(jac(u, par, 0))); Pt = Mo * Pt
            frot *= norm(Matrix(tangent(sys, zs[k])))
            Pm = Matrix{Float64}(I, d, d)
            for i in 1:k
                Pm = Matrix(tangent(sys, zs[i])) * Pm
            end
            worst_tf = max(worst_tf, maximum(abs, Pt - Pm) / (8 * k * cm_eps(Float64) * frot))
        end
    end
    return worst_orb, worst_prod, worst_tf
end

section("AC21") do
    m = sz(100)
    say("# AC21: ChaosTools $(pkgversion(CT)) DeterministicIteratedMap + TangentDynamicalSystem(J = hand-coded), $m seeded orbits per system, z0 uniform in [-4 pi, 4 pi]^d, n <= n_H; Q0 = identity")
    say("# AC21 convention: Y3 hand system in the labelling of Skokos2003a Eq. (9) (x1, x2 = position, momentum of map 1; x3, x4 of map 2; right sides at the pre-step values); conversion to the module labelling (positions first) by the permutation [1, 3, 2, 4] of states and of the Jacobian product (rows and columns); Y4 hand system in radians, open chain, same ordering (q_1..q_N, p_1..p_N), eps = 2 beta")
    for p in (Y3PUB, (2.0, 0.5, 0.1))
        c = cm_y3(p...)
        mk() = CT.TangentDynamicalSystem(CT.DeterministicIteratedMap(ct3_rule, SVector(0.1, 0.1, 0.1, 0.1), collect(p)); J = ct3_jac, k = 4)
        o, pr, tf = ct_compare("Y3 $p", c, mk, z -> SVector{4, Float64}(z[PERM]), Q -> Q[PERM, PERM], ct3_jac, collect(p), m, CM_SEED + 61001)
        ratio_line("AC21", "Y3 $p: max ratio torus distance (ChaosTools orbit vs module) / tau(k), k <= n_H", o)
        ratio_line("AC21", "Y3 $p: max ratio |Q_ChaosTools - prod M_k| / 8 k eps prod ||M_k||_F, k <= n_H (each side along its own orbit; the bound of the spec)", pr)
        ratio_line("AC21", "Y3 $p: same with the oracle Jacobian evaluated at the module's orbit points (teacher forced; isolates the formula)", tf)
    end
    for N in (2, 3), K in (1.0, 3.0), e in (0.2, 0.6102)
        c = cm_y4(N, K, e)
        mk() = CT.TangentDynamicalSystem(CT.DeterministicIteratedMap(ct4_rule(Val(N)), zero(SVector{2N, Float64}), [K, e]); J = ct4_jac(Val(N)), k = 2N)
        o, pr, tf = ct_compare("Y4", c, mk, z -> SVector{2N, Float64}(z), identity, ct4_jac(Val(N)), [K, e], m, CM_SEED + 61002 + N)
        ratio_line("AC21", "Y4 N=$N K=$K eps=$e: max ratio torus distance / tau(k), k <= n_H", o)
        ratio_line("AC21", "Y4 N=$N K=$K eps=$e: max ratio tangent product difference / 8 k eps prod ||M_k||_F (own orbits)", pr)
        ratio_line("AC21", "Y4 N=$N K=$K eps=$e: same, oracle Jacobian at the module's orbit points (teacher forced)", tf)
    end
end

# ---------------------------------------------------------------------------------------------------
# Figure: AC20 series of orbits C and D (finite-time exponent, SALI, pairing)
section("figure") do
    f = joinpath(tempdir(), "cm-ac20-series.bin")
    isfile(f) || error("run section AC20 first (series file $f)")
    DATA = Serialization.deserialize(f)
    ns, C, D = DATA.ns, DATA.serC, DATA.serD
    fig = CM.Figure(size = (1500, 480))
    ax1 = CM.Axis(fig[1, 1]; xscale = log10, yscale = log10, xlabel = "n (iterations)", ylabel = "finite-time maximal exponent L(n)", title = "Y3 orbits C (ordered) and D (chaotic)")
    CM.lines!(ax1, ns, max.(C.L, 1e-12); color = :steelblue, linewidth = 2, label = "C (0.5, 0.5, 0, 0)")
    CM.lines!(ax1, ns, D.L; color = :crimson, linewidth = 2, label = "D (3, 0.5, 0, 0)")
    ref = [log(n) / n for n in ns]
    CM.lines!(ax1, ns, ref .* (C.L[end] / ref[end]); color = :black, linestyle = :dash, label = "ln(n)/n law, scaled to C at 1e6")
    CM.hlines!(ax1, [8e-3]; color = :gray50, linestyle = :dot, label = "8e-3 (Skokos2010 Fig. 7)")
    CM.hspan!(ax1, 5e-3, 1.5e-2; color = (:gray, 0.15))
    CM.axislegend(ax1; position = :lb, labelsize = 11)
    ax2 = CM.Axis(fig[1, 2]; xscale = log10, yscale = log10, xlabel = "n (iterations)", ylabel = "SALI(n)", title = "SALI, deviation vectors (1,1,1,1), (1,0,0,0)")
    CM.lines!(ax2, ns, max.(C.sali, 1e-17); color = :steelblue, linewidth = 2, label = "C")
    CM.lines!(ax2, ns, max.(D.sali, 1e-17); color = :crimson, linewidth = 2, label = "D")
    CM.hlines!(ax2, [1e-12]; color = :black, linestyle = :dash, label = "1e-12 threshold of (b)")
    CM.hspan!(ax2, 0.1, 0.5; color = (:steelblue, 0.12))
    CM.vspan!(ax2, 1e3, 2e4; color = (:crimson, 0.08))
    CM.axislegend(ax2; position = :lb, labelsize = 11)
    ax3 = CM.Axis(fig[1, 3]; xscale = log10, yscale = log10, xlabel = "n (iterations)", ylabel = "pairing defect and sum (QR exponents)", title = "orbit D: lambda_i(n) + lambda_5-i(n), sum")
    CM.lines!(ax3, ns, max.(abs.(D.lam[1, :] .+ D.lam[4, :]), 1e-17); color = :crimson, linewidth = 2, label = "|l1 + l4|")
    CM.lines!(ax3, ns, max.(abs.(D.lam[2, :] .+ D.lam[3, :]), 1e-17); color = :darkorange, linewidth = 2, label = "|l2 + l3|")
    CM.lines!(ax3, ns, max.(abs.(vec(sum(D.lam; dims = 1))), 1e-17); color = :black, linewidth = 1.5, label = "|sum l_i|")
    CM.lines!(ax3, ns, 1 ./ ns; color = :gray50, linestyle = :dot, label = "1/n")
    CM.axislegend(ax3; position = :lb, labelsize = 11)
    CM.save(joinpath(OUT, "coupled-maps-orbits-cd.png"), fig)
    say("# figure saved: coupled-maps-orbits-cd.png")
end

say(@sprintf("# total wall time %.1f s", time() - T0))
close(LOG)
