# Stage 04 verification of module standard-map (Systems: Y1 standard map, Y5 linear symplectic maps, D28 wrapped map).
# Run from the repo root:  julia --project=verify -t 4 verify/standard-map.jl
# Spec: workspaces/module-build/stages/01-spec/output/standard-map-spec.md (criteria 1 to 19, amendment D28 criteria W1 to W8).
# One printed line per criterion:  AC<n> | quantity | value | bound | PASS/FAIL/STAGE04-GPU
# Most "value" entries are the worst ratio error / bound (pass when <= 1), as in the package tests.
# Seeds: StableRNG(20261012 + offset); the offsets used here start at 40000 so the sample differs from the package tests.
using Printf, Random, LinearAlgebra, Statistics
using Chaoslab, StaticArrays, DoubleFloats, StableRNGs, ForwardDiff, BenchmarkTools
import ChaosTools as CT
import CairoMakie as CM

const ROOT = dirname(@__DIR__)
const OUT = joinpath(ROOT, "workspaces", "module-build", "stages", "04-verify", "output")
mkpath(OUT)
# SM_ONLY="W,AC18" runs only the named sections and writes the log to SM_LOG (debugging aid); the default runs everything
const ONLY = filter(!isempty, split(get(ENV, "SM_ONLY", ""), ","))
const LOG = open(get(ENV, "SM_LOG", joinpath(OUT, "standard-map-verify-log.txt")), "w")

fmt(x::AbstractFloat) = @sprintf("%.4g", x)
fmt(x) = string(x)
function say(s)
    println(s); println(LOG, s); flush(stdout); flush(LOG)
end
function emit(ac, quantity, value, bound, verdict)
    say(string(ac, " | ", quantity, " | ", fmt(value), " | ", bound, " | ", verdict))
end
pf(ok::Bool) = ok ? "PASS" : "FAIL"
ratio_line(ac, quantity, r) = emit(ac, quantity, r, "<= 1", pf(r <= 1))

# Reuse the measurement functions of the package tests: include each test file up to its first @testset.
function include_defs(path)
    txt = read(path, String)
    i = findfirst("\n@testset", txt)
    include_string(Main, i === nothing ? txt : txt[1:first(i)], path)
end
const TESTDIR = joinpath(ROOT, "test", "standard-map")
include(joinpath(TESTDIR, "sm_helpers.jl"))
for f in ("sm_ac01_06.jl", "sm_ac07_13.jl", "sm_ac14.jl", "sm_ac18.jl", "sm_y5.jl", "sm_wrapped.jl")
    include_defs(joinpath(TESTDIR, f))
end

# Run a section; an exception is reported as an ERROR line and does not stop the other sections.
function section(f, name)
    isempty(ONLY) || name in ONLY || return
    t0 = time()
    try
        f()
    catch e
        say(string(name, " | ERROR | ", sprint(showerror, e)[1:min(end, 300)], " | - | FAIL"))
    end
    say(@sprintf("# section %s took %.1f s", name, time() - t0))
end

const T0 = time()
say("# standard-map verification, julia $(VERSION), threads $(Threads.nthreads()), base seed $SM_SEED, time $(round(time())) (unix)")
say("# timings in this log were taken on a loaded machine (other agents were running concurrently)")
const NH = Dict(K => sm_nH(K) for K in SM_KS)
say("# n_H(K) = " * join(["$K:$(NH[K])" for K in sort(SM_KS)], " "))

# ---------------------------------------------------------------------------------------------------
# AC1 to AC6: map invariants, re-measured (Float64 500 per K, Double64 100 per K; AC3/5/6: 200 orbits per listed K)
section("AC1-6") do
    nf, nd = 500, 100
    cF = sm_cases(StableRNG(SM_SEED + 40001), nf)
    cD = sm_cases(StableRNG(SM_SEED + 40002), nd)
    ad, cs = sm_ac1_ratios(Float64, cF; complexstep = true)
    ad2, _ = sm_ac1_ratios(Double64, cD; complexstep = false)
    ratio_line("AC1", "max ratio |DT - J_ForwardDiff| / 8 eps (1+|K|), Float64, $(length(cF)) pts", ad)
    ratio_line("AC1", "max ratio |DT - J_complexstep(1e-30)| / 8 eps (1+|K|), Float64, $(length(cF)) pts", cs)
    ratio_line("AC1", "max ratio |DT - J_ForwardDiff| / 8 eps (1+|K|), Double64 (eps 2^-100), $(length(cD)) pts", ad2)
    emit("AC1", "same Double64 ratio against nominal eps 2^-104 (ratio x16)", 16 * ad2, "<= 1", pf(16 * ad2 <= 1))
    r2f = sm_ac2_ratios(Float64, cF)
    r2d = sm_ac2_ratios(Double64, cD)
    ratio_line("AC2", "max ratio of |det DT - 1| and |DT' Om DT - Om| to 8 eps (1+|K|), Float64", r2f)
    ratio_line("AC2", "same, Double64 (eps 2^-100)", r2d)
    emit("AC2", "same Double64 ratio against nominal eps 2^-104 (ratio x16)", 16 * r2d, "<= 1", pf(16 * r2d <= 1))
    c3 = sm_cases(StableRNG(SM_SEED + 40003), 200; uniform = false)
    ratio_line("AC3", "max ratio round-trip error / tau(n), n in {1, n_H}, both orders, Float64, $(length(c3)) orbits", sm_ac3_ratio(Float64, c3))
    ratio_line("AC3", "same, Double64", sm_ac3_ratio(Double64, c3))
    r4f = sm_ac4_ratio(Float64, cF)
    r4d = sm_ac4_ratio(Double64, cD)
    ratio_line("AC4", "max ratio |DT DT_inv(Tz) - I| / 8 eps (1+|K|)^2, Float64", r4f)
    ratio_line("AC4", "same, Double64 (eps 2^-100)", r4d)
    c5 = sm_cases(StableRNG(SM_SEED + 40005), 200; uniform = false)
    ratio_line("AC5", "max ratio torus distance lifted-wrapped vs wrap-each-step / tau(k), k <= n_H, Float64, $(length(c5)) orbits", sm_ac5_ratio(Float64, c5))
    ratio_line("AC5", "same, Double64", sm_ac5_ratio(Double64, c5))
    c6 = sm_cases(StableRNG(SM_SEED + 40006), 200; uniform = false)
    ratio_line("AC6", "max ratio torus distance K / -K symmetry / tau(k), k <= n_H, Float64, $(length(c6)) orbits", sm_ac6_ratio(Float64, c6))
    ratio_line("AC6", "same, Double64", sm_ac6_ratio(Double64, c6))
end

# ---------------------------------------------------------------------------------------------------
# AC7 to AC13: known answers
section("AC7-13") do
    d1 = sm_ac7_orbit_defects(Float64, StableRNG(SM_SEED + 40071), 1000)
    d2 = sm_ac7_orbit_defects(Double64, StableRNG(SM_SEED + 40072), 200)
    p1 = sm_ac7_product_defects(Float64, StableRNG(SM_SEED + 40073), 20)
    p2 = sm_ac7_product_defects(Double64, StableRNG(SM_SEED + 40074), 20)
    emit("AC7", "orbit defects x_n != x_0+n p_0 or p_n != p_0 (1000 steps; Float64 1000 orbits, Double64 200 orbits)", string(d1, " / ", d2), "0", pf(d1 == 0 && d2 == 0))
    emit("AC7", "tangent product defects vs [[1,n],[0,1]] for n in {1,2,10,100,1000,1e4} (20 starts each type)", string(p1, " / ", p2), "0", pf(p1 == 0 && p2 == 0))

    ns = (1, 10, 100, 1000, 10_000, 1_000_000)
    for T in (Float64, Double64)
        prods = sm_ac8_products(T, ns)
        wr, wx, wl = 0.0, 0.0, 0.0
        setprecision(BigFloat, SM_BITS) do
            for n in ns
                P = prods[n]
                X = sm_expo(P, SVector(zero(T), one(T)), n)
                ref = log(sqrt(BigFloat(n)^2 + 1)) / n
                wr = max(wr, Float64(abs(sm_big(X) - ref)) / (8 * sm_eps(T)))
                wx = max(wx, Float64(abs(sm_expo(P, SVector(one(T), zero(T)), n))))
            end
        end
        X4 = sm_expo(prods[10_000], SVector(zero(T), one(T)), 10_000)
        wl = abs(Float64(10_000 * X4 / log(T(10_000))) - 1)
        ratio_line("AC8", "max ratio |X(n) - ln(sqrt(n^2+1))/n| / 8 eps, n up to 1e6, w0=(0,1), $T", wr)
        emit("AC8", "|n X(n)/ln n - 1| at n = 1e4, $T (expected 5.4e-10)", wl, "<= 1e-9", pf(wl <= 1e-9))
        emit("AC8", "max |X(n)| for w0=(1,0), $T", wx, "0", pf(wx == 0))
    end

    for T in (Float64, Double64)
        wt, we, wc = sm_ac9_ratios(T)
        jd = sm_ac9_jordan_defects(T)
        ratio_line("AC9", "max ratio |tr - (2 + K cos x)| / 8 eps (1+|K|), $T", wt)
        ratio_line("AC9", "max ratio |lambda_+ or angle - closed form| / 32 eps (1+|K|), $T", we)
        emit("AC9", "wrong elliptic/hyperbolic classes (12 cases), $T", wc, "0", pf(wc == 0))
        emit("AC9", "Jordan blocks: |max entry of DT^n / n - (4 or 1)| for n in {10,100,1000}, $T", jd, "<= 1e-12", pf(jd <= 1e-12))
    end

    for T in (Float64, Double64)
        bad = 0
        for K in SM_KS
            z = SVector(zero(T), zero(T)); sys = StandardMap(T(K))
            for _ in 1:10_000
                z = forward(sys, z)
                bad += (z[1] == 0 && z[2] == 0) ? 0 : 1
            end
        end
        emit("AC10", "steps that left (0,0) in 1e4 steps, 17 K values, $T", bad, "0", pf(bad == 0))
    end

    for T in (Float64, Double64)
        worst0, worstp, wrongc = 0.0, 0.0, 0
        for Ks in ("1", "3", "3.9", "4.1", "5")
            K = exact(T, Ks)
            sys = StandardMap(K)
            tol = 8 * sm_eps(T) * (1 + abs(Float64(K)))
            xc = 2 * exact(T, pi) * exact(T, 1 // 2)
            DTc = tangent(sys, SVector(xc, zero(T)))
            DT0 = tangent(sys, SVector(zero(T), zero(T)))
            trc = DTc[1, 1] + DTc[2, 2]
            worstp = max(worstp, sm_absdiff((trc,), (2 - K,)) / tol)
            worst0 = max(worst0, sm_absdiff((DT0[1, 1] + DT0[2, 2],), (2 + K,)) / tol)
            wrongc += sm_class(trc)[1] == (Float64(K) < 4 ? :elliptic : :hyperbolic) ? 0 : 1
        end
        ratio_line("AC11", "max ratio |tr(DT)(pi,0) - (2-K)| / 8 eps (1+|K|), K in {1,3,3.9,4.1,5}, $T", worstp)
        ratio_line("AC11", "max ratio |tr(DT)(0,0) - (2+K)| / 8 eps (1+|K|), $T", worst0)
        emit("AC11", "wrong classes (elliptic K<4, hyperbolic K>4), $T", wrongc, "0", pf(wrongc == 0))
    end

    K = 1.0
    s = sqrt(K * (K + 4))
    DT = tangent(StandardMap(K), SVector(0.0, 0.0))
    tr = DT[1, 1] + DT[2, 2]
    lp = (tr + sqrt(tr^2 - 4)) / 2; lm = (tr - sqrt(tr^2 - 4)) / 2
    ru = (lp - DT[1, 1]) / DT[1, 2]; rs = (lm - DT[1, 1]) / DT[1, 2]
    e1 = max(abs(ru - (s - K) / 2), abs(rs + (K + s) / 2))
    e2 = max(maximum(abs, DT * SVector(1.0, ru) - lp * SVector(1.0, ru)), maximum(abs, DT * SVector(1.0, rs) - lm * SVector(1.0, rs)))
    e3 = max(abs(0.62 - abs(ru)), abs(0.62 - 1 / abs(rs)))
    emit("AC12", "v_p/v_x unstable $(fmt(ru)), stable $(fmt(rs)); max error vs derived (s-K)/2 and -(K+s)/2", e1, "<= 1e-12", pf(e1 <= 1e-12))
    emit("AC12", "eigen-equation residual", e2, "<= 1e-12", pf(e2 <= 1e-12))
    emit("AC12", "printed 0.62 vs derived magnitudes 0.618 (loose check)", e3, "<= 0.005", pf(e3 <= 0.005))

    sys = StandardMap(5.0)
    tp = 2 * Float64(pi); tol = tp * 1.0e-4
    O1 = SVector(tp * 0.68002, tp * 0.36002); P1 = SVector(tp * 0.642118, tp * 0.369344)
    dO = sm_tdist(sm_orbit(sys, O1, 2)[3], O1)
    zsP = sm_orbit(sys, P1, 18)
    dP18 = sm_tdist(zsP[19], P1); dP9 = sm_tdist(zsP[10], P1)
    emit("AC13", "torus distance O1 after 2 steps (K=5), tolerance 2 pi 1e-4 = $(fmt(tol))", dO, "<= $(fmt(tol))", pf(dO <= tol))
    emit("AC13", "torus distance P1 after 18 steps", dP18, "<= $(fmt(tol))", pf(dP18 <= tol))
    emit("AC13", "torus distance P1 after 9 steps (must exceed tolerance)", dP9, "> $(fmt(tol))", pf(dP9 > tol))
end

# ---------------------------------------------------------------------------------------------------
# AC14 heavy part: larger 256-bit sample (1000 orbits per listed K) and the effective Double64 accuracy of sin, cos
rand_dd(rng, lim) = (u = (2 * rand(rng) - 1) * lim; Double64(u) + Double64((rand(rng) - 0.5) * eps(u)))

function dd_trig_stats(xs)
    rels = Float64[]; args = Double64[]; which = String[]
    setprecision(BigFloat, SM_BITS) do
        for x in xs
            b = sm_big(x)
            for (nm, f) in (("sin", sin), ("cos", cos))
                r = f(b)
                iszero(r) && continue
                push!(rels, Float64(abs(sm_big(f(x)) - r) / abs(r))); push!(args, x); push!(which, nm)
            end
        end
    end
    i = argmax(rels)
    return (max = rels[i], arg = args[i], fn = which[i], q999 = quantile(rels, 0.999), median = median(rels), n = length(rels))
end

section("AC14") do
    n_per = 1000
    cases = sm_cases(StableRNG(SM_SEED + 40141), n_per; uniform = false)
    setprecision(BigFloat, SM_BITS) do
        r64 = sm_ac14_ratio(Float64, cases, sm_eps(Float64))
        ratio_line("AC14", "max ratio |orbit_Float64 - orbit_256| / tau(n_H), lifted and torus, $(length(cases)) orbits", r64)
        # effective Double64 eps: orbit arguments of the case set, uniform sample on [-4pi, 4pi], uniform on [-1e3, 1e3]
        xs_orb = Double64[]
        for (K, z) in cases[1:3:end]
            append!(xs_orb, first.(sm_orbit(StandardMap(Double64(K)), sm_cast(Double64, z), sm_nH(K))))
        end
        rng = StableRNG(SM_SEED + 40142)
        xs_u = [rand_dd(rng, 4 * Float64(pi)) for _ in 1:50_000]
        xs_w = [rand_dd(rng, 1.0e3) for _ in 1:20_000]
        sO, sU, sW = dd_trig_stats(xs_orb), dd_trig_stats(xs_u), dd_trig_stats(xs_w)
        for (nm, st) in (("orbit arguments", sO), ("uniform [-4pi,4pi]", sU), ("uniform [-1e3,1e3]", sW))
            say(@sprintf("# Double64 sin/cos relative error vs BigFloat-256, %s: %d values, max %.3g (%s at x=%s), 99.9%% %.3g, median %.3g",
                         nm, st.n, st.max, st.fn, string(Float64(st.arg)), st.q999, st.median))
        end
        eff = max(SM_EPS_DD_FLOOR, sO.max, sU.max, sW.max)
        emit("AC14", "effective Double64 eps = max relative sin/cos error over the three samples ($(sO.n + sU.n + sW.n) evaluations)", eff, "<= 2^-100 = $(fmt(ldexp(1.0, -100)))", pf(eff <= ldexp(1.0, -100)))
        emit("AC14", "same, as a fraction of 2^-100 and as a power of two (log2)", string(fmt(eff / ldexp(1.0, -100)), " ; ", fmt(log2(eff))), "<= 1", pf(eff <= ldexp(1.0, -100)))
        rdd = sm_ac14_ratio(Double64, cases, eff)
        ratio_line("AC14", "max ratio |orbit_Double64 - orbit_256| / tau(n_H) with eps = effective value above", rdd)
        rdd_nom = sm_ac14_ratio(Double64, cases, SM_EPS_DD_FLOOR)
        emit("AC14", "same ratio with the nominal eps 2^-104 (informational)", rdd_nom, "<= 1", pf(rdd_nom <= 1))
        # sanity: Chaoslab 256-bit orbit vs independent plain-formula 256-bit orbit
        worst = 0.0
        for (K, z) in cases[1:400]
            nH = sm_nH(K)
            a = sm_orbit(StandardMap(BigFloat(K)), sm_cast(BigFloat, z), nH)
            b = sm_ref_orbit(BigFloat(K), sm_cast(BigFloat, z), nH)
            worst = max(worst, sm_absdiff(a[end], b[end]) / sm_tau(ldexp(1.0, 1 - SM_BITS), K, a, nH))
        end
        ratio_line("AC14", "Chaoslab 256-bit orbit vs independent plain 256-bit orbit (400 orbits)", worst)
    end
end

# ---------------------------------------------------------------------------------------------------
# AC15 to AC17: Y5
section("AC15-17") do
    rng = StableRNG(SM_SEED + 40151)
    for T in (Float64, Double64)
        w, ok, omok = 0.0, true, true
        for N in (1, 2, 3, 5), _ in 1:10
            a, b, c = sm_ac15_ratios(T, rng, Val(N))
            w = max(w, a); ok &= b; omok &= c
        end
        ratio_line("AC15", "max ratio of |A' Om A - Om| and |det A - 1| to 16 eps ||A||^p, A and A_S, N in {1,2,3,5}, 10 draws, $T (eps 2^-100 for Double64)", w)
        emit("AC15", "map contract (dimension, tangent == A, forward = A z) and Omega block form hold, $T", ok && omok, "true", pf(ok && omok))
    end
    rng = StableRNG(SM_SEED + 40161)
    wn, ws = 0.0, 0.0
    for N in (1, 2, 3, 5), _ in 1:10
        a, b = sm_ac16_ratios(rng, Val(N))
        wn = max(wn, a); ws = max(ws, b)
    end
    ratio_line("AC16", "max ratio |(1/n) ln smax(A^n) - lambda_1|, |(1/n) ln smax(A^-n) + lambda_2N| to 16 eps (block-diagonal), n in {10,50}, Float64", wn)
    ratio_line("AC16", "same for A_S with bound 16 eps + ln(kappa(S))/n", ws)
    rng = StableRNG(SM_SEED + 40171)
    for T in (Float64, Double64)
        w = 0.0
        for N in (1, 2, 3, 5), _ in 1:10
            w = max(w, sm_ac17_ratios(T, rng, Val(N)))
        end
        ratio_line("AC17", "max ratio of inverse residuals (A^-1 A, A A^-1, Om^-1 A' Om vs block inverse) to 16 eps ||A||^2, $T", w)
    end
end

# ---------------------------------------------------------------------------------------------------
# O6: ForwardDiff n-step Jacobian vs the product of analytic tangents, 10^4 seeded points, n = n_H(K)
function o6_point(K, z)
    n = sm_nH(K)
    sys = StandardMap(K)
    P = SMatrix{2, 2}(1.0, 0.0, 0.0, 1.0); zz = z; fro = 1.0
    pts = Vector{typeof(z)}(undef, n)
    for k in 1:n
        pts[k] = zz
        DT = tangent(sys, zz)
        P = DT * P; fro *= norm(DT); zz = forward(sys, zz)
    end
    J = ForwardDiff.jacobian(v -> (w = v; for _ in 1:n; w = forward(sys, w); end; w), z)
    bound = 8 * n * ldexp(1.0, -52) * fro * sqrt(2.0)
    rJ = Float64(maximum(abs, SMatrix{2, 2}(J) - P)) / bound
    rB, rJB = setprecision(BigFloat, SM_BITS) do
        sysB = StandardMap(BigFloat(K))
        Pb = BigFloat[1 0; 0 1]
        for k in 1:n
            Mb = tangent(sysB, SVector(BigFloat(pts[k][1]), BigFloat(pts[k][2])))
            Pb = Mb * Pb
        end
        dP = maximum(abs(BigFloat(P[i, j]) - Pb[i, j]) for i in 1:2, j in 1:2)
        dJ = maximum(abs(BigFloat(J[i, j]) - Pb[i, j]) for i in 1:2, j in 1:2)
        (Float64(dP) / bound, Float64(dJ) / bound)
    end
    return (rJ, rB, rJB)
end

section("O6") do
    cases = sm_cases(StableRNG(SM_SEED + 40601), 556)
    res = Vector{NTuple{3, Float64}}(undef, length(cases))
    Threads.@threads for i in eachindex(cases)
        K, z = cases[i]
        res[i] = o6_point(K, z)
    end
    a = maximum(first, res); b = maximum(r -> r[2], res); c = maximum(last, res)
    say("# O6 sample: $(length(cases)) points = 556 for each of the 17 listed K and 556 with K uniform in [-20, 20]; z0 uniform in [-4 pi, 4 pi]^2; n = n_H(K); bound = 8 n eps ||W0||_F prod_k ||DT_k||_F (||W0||_F = sqrt 2), derived (spec O6 states no bound)")
    ratio_line("AC1(O6)", "max ratio |J_ForwardDiff(n-step) - prod DT_k| / bound, $(length(cases)) points", a)
    ratio_line("AC1(O6)", "max ratio |prod DT_k (Float64) - product at 256 bits, same realised points| / bound", b)
    ratio_line("AC1(O6)", "max ratio |J_ForwardDiff - product at 256 bits| / bound", c)
end

# ---------------------------------------------------------------------------------------------------
# AC18: CPU backend of the kernel vs the plain loop, threaded vs serial; GPU parts are STAGE04-GPU
section("AC18") do
    rng = StableRNG(SM_SEED + 40181)
    lim = 4 * Float64(pi)
    Z = [SVector((2 * rand(rng) - 1) * lim, (2 * rand(rng) - 1) * lim) for _ in 1:100_000]
    I2 = SMatrix{2, 2, Float64}(1.0, 0.0, 0.0, 1.0)
    W0 = fill(I2, length(Z))
    for K in (1.0, 5.0)
        sys = StandardMap(K); n = sm_nH(K)
        plain = sm_run_serial(sys, Z, n)
        thr = sm_run_threaded(sys, Z, n)
        t = @elapsed (res = propagate(sys, Z, n; k = 2, tau = 1, seed = SM_SEED, W0 = W0, backend = :cpu))
        res2 = propagate(sys, Z, n; k = 2, tau = 1, seed = SM_SEED, W0 = W0, backend = :cpu, chunk = 777)
        ser = [propagate_orbit(sys, Z[j], I2, n, 1, ())[1] for j in 1:1000]
        Wplain = map(1:length(Z)) do j
            z = Z[j]; W = I2
            for _ in 1:n
                W = tangent(sys, z) * W; z = forward(sys, z)
            end
            W
        end
        bit_thr = isequal(plain, thr)
        bit_kernel = isequal(res.x, plain) && isequal(res.W, Wplain)
        bit_chunk = isequal(res.x, res2.x) && isequal(res.W, res2.W)
        bit_ser = isequal(res.x[1:1000], ser)
        maxd = maximum(sm_tdist(res.x[j], plain[j]) for j in 1:2000)
        emit("AC18", "K=$K n=$n: $(length(Z)) orbits, propagate (CPU kernel, $(Threads.nthreads()) threads) vs plain loop: bitwise equal x and W; torus distance max over first 2000 orbits = $(fmt(maxd)); propagate wall $(fmt(t)) s (loaded machine)", bit_kernel, "true (bound tau(n))", pf(bit_kernel))
        emit("AC18", "K=$K: threaded forward loop vs serial loop, bitwise", bit_thr, "true", pf(bit_thr))
        emit("AC18", "K=$K: propagate with chunk 777 vs default chunk, and vs serial propagate_orbit (first 1000), bitwise", bit_chunk && bit_ser, "true", pf(bit_chunk && bit_ser))
    end
    emit("AC18", "GPU vs CPU max torus distance, 1e5 orbits, K in {1,5}; Double64 kernel gate (compile or reason)", "see verify/gpu.jl", "tau(n)", "STAGE04-GPU")
    emit("AC18", "literal 1 vs 8 thread comparison: only serial vs $(Threads.nthreads())-thread executed here (orbits are independent; each orbit is one work item)", "n/a", "bitwise", "NOTE")
end

# ---------------------------------------------------------------------------------------------------
# AC19 (oracle O5): ChaosTools standard map rule
# Convention (read from the package source of the verify environment, see the report):
#   ChaosTools 3.5.4, test/chaosdetection/01test.jl lines 3-12 and docs/src/chaos_detection.md (rule `standardmap_rule`):
#     state (theta, p), parameter k = par[1]:  p += k sin(theta);  theta += p;  both reduced to [0, 2 pi) after every step.
#   Canonical form of this module (spec Equations item 1): p' = p + K sin x, x' = x + p'. So x = theta, p = p, K = k:
#   the conversion is the identity and the state ordering is the same (x, p). The only difference is wrapping.
#   No predefined standard map exists in DynamicalSystemsBase 3.19.3 / ChaosTools 3.5.4 (PredefinedDynamicalSystems is not in the manifest).
function ct_rule_mod2pi(x, par, n)
    theta = x[1]; p = x[2]
    p += par[1] * sin(theta)
    theta += p
    return mod2pi.(SVector(theta, p))
end
function ct_rule_while(x, par, n)
    theta = x[1]; p = x[2]
    p += par[1] * sin(theta)
    theta += p
    while theta >= 2π; theta -= 2π; end
    while theta < 0; theta += 2π; end
    while p >= 2π; p -= 2π; end
    while p < 0; p += 2π; end
    return SVector(theta, p)
end

function ct_orbit(rule, K, z0, n)
    ds = CT.DeterministicIteratedMap(rule, z0, [K])
    out = Vector{SVector{2, Float64}}(undef, n + 1)
    out[1] = SVector{2, Float64}(CT.current_state(ds))
    for k in 1:n
        CT.step!(ds)
        out[k + 1] = SVector{2, Float64}(CT.current_state(ds))
    end
    return out
end

section("AC19") do
    rng = StableRNG(SM_SEED + 40191)
    lim = 4 * Float64(pi)
    Ks = (0.5, 1.0, 2.0, 5.0)
    sample = [(K, SVector((2 * rand(rng) - 1) * lim, (2 * rand(rng) - 1) * lim)) for K in Ks for _ in 1:25]
    say("# AC19 sample: 100 seeded orbits, 25 per K in $(Ks), z0 uniform in [-4 pi, 4 pi]^2 handed unwrapped to both; ChaosTools 3.5.4 / DynamicalSystemsBase 3.19.3")
    for (nm, rule) in (("mod2pi rule (docs/src/chaos_detection.md)", ct_rule_mod2pi), ("while-loop rule (test/chaosdetection/01test.jl)", ct_rule_while))
        worst = 0.0; worst_perK = Dict(K => 0.0 for K in Ks); sanity = 0.0; wrapworst = 0.0
        for (K, z0) in sample
            nH = sm_nH(K)
            ours = sm_orbit(StandardMap(K), z0, 100)
            ct = ct_orbit(rule, K, z0, 100)
            for k in 1:nH
                r = sm_tdist(ours[k + 1], ct[k + 1]) / sm_tau(sm_eps(Float64), K, ours, k)
                worst = max(worst, r); worst_perK[K] = max(worst_perK[K], r)
                # secondary: our lifted orbit wrapped in Float64 by Chaoslab.wrap first (same quantity in the module's own wrap)
                wrapworst = max(wrapworst, sm_tdist(wrap(ours[k + 1]), ct[k + 1]) / sm_tau(sm_eps(Float64), K, ours, k))
            end
            sanity = max(sanity, sm_tdist(ours[101], ct[101]))
        end
        ratio_line("AC19", "$nm: max ratio torus distance (ChaosTools vs Chaoslab, k <= n_H) / tau(k), 100 orbits", worst)
        say("#   per-K worst ratio: " * join(["K=$K: $(fmt(worst_perK[K]))" for K in Ks], ", "))
        emit("AC19", "$nm: same with Chaoslab.wrap(lifted state) in Float64 instead of the exact reduction", wrapworst, "<= 1", pf(wrapworst <= 1))
        emit("AC19", "$nm: unscored sanity, max torus distance at n = 100 (orbits decorrelated, chaotic)", sanity, "unscored", "INFO")
    end
end

# ---------------------------------------------------------------------------------------------------
# D28 amendment criteria W1 to W8
section("W") do
    ok1, w1 = smw_w1(Float64, 500, SM_SEED + 40501)
    ok1d, w1d = smw_w1(Double64, 200, SM_SEED + 40502)
    ratio_line("W1", "max ratio torus distance wrapped state vs lifted orbit / tau(k), range and integrality ok=$ok1, Float64, 500 orbits per K in {0.5,2,5,10}", w1)
    ratio_line("W1", "same, Double64 (eps 2^-100), range and integrality ok=$ok1d, 200 orbits per K", w1d)
    ratio_line("W2", "max ratio |lift(state) at 256 bits - lifted orbit| / tau(k), Float64, 500 orbits per K", smw_w2(Float64, 500, SM_SEED + 40511))
    ratio_line("W2", "same, Double64, 200 orbits per K", smw_w2(Double64, 200, SM_SEED + 40512))

    t3 = @elapsed ((inrange, ints, bounded, first51) = smw_w3_run(5.0, SVector(1.0, 2.0, 0.0, 0.0), 1_000_000))
    zs = sm_orbit(StandardMap(5.0), SVector(1.0, 2.0), sm_nH(5.0))
    forced = smw_w3_forced(5.0, first51)
    wfree, same = smw_w3_free(5.0, first51, zs, sm_nH(5.0))
    emit("W3", "1e6 steps at K=5 (Float64): x, p in [0,2 pi), m, n integer, |m|,|n| < 2^53 ($(fmt(t3)) s, loaded machine)", string(inrange, "/", ints, "/", bounded), "true/true/true", pf(inrange && ints && bounded))
    emit("W3", "teacher-forced 256-bit increments a, b equal binary64 increments, first 50 steps", forced, "true", pf(forced))
    ratio_line("W3", "free-running 256-bit orbit vs binary64 wrapped state, n <= n_H(5): max ratio torus distance / tau", wfree)
    emit("W3", "free-running 256-bit counts equal binary64 counts, n <= n_H(5)", same, "true", pf(same))
    # extra: counts at 1e6 steps, magnitude
    s = SVector(1.0, 2.0, 0.0, 0.0); sysw = WrappedStandardMap(5.0)
    for _ in 1:1_000_000
        s = forward(sysw, s)
    end
    say("# W3 info: state after 1e6 steps at K=5: m = $(s[3]), n = $(s[4])")

    w4 = smw_w4(SM_SEED + 40541)
    emit("W4", "K=0: p, n unchanged and m_k = floor((x0 + k p0)/2pi) at 256 bits, k <= 1000, 20 dyadic starts", w4, "true", pf(w4))
    a, b = smw_w5(Float64, 200, SM_SEED + 40551)
    ratio_line("W5", "max ratio wrapped tangent product vs 256-bit product / (8 n eps prod||DT||_F ||W0||_F), Float64, 200 orbits per K", a)
    ratio_line("W5", "max ratio wrapped vs lifted tangent product / (2 x bound + orbit term), Float64", b)
    a, b = smw_w5(Double64, 100, SM_SEED + 40552)
    ratio_line("W5", "same, Double64 (eps 2^-100), first ratio", a)
    ratio_line("W5", "same, Double64, second ratio", b)
    ok, w = smw_w6(Float64, 200, SM_SEED + 40561)
    ratio_line("W6", "max ratio round trip torus distance / tau(n), counts exactly equal=$ok, Float64, 200 orbits per K", w)
    ok, w = smw_w6(Double64, 200, SM_SEED + 40562)
    ratio_line("W6", "same, Double64, counts exactly equal=$ok", w)
    for T in (Float64, Double64)
        sys = WrappedStandardMap(T(5))
        s = SVector{4, T}(T(1), T(2), zero(T), zero(T))
        W0 = SMatrix{2, 2, T}(one(T), zero(T), zero(T), one(T))
        al = smw_allocs(sys, s, W0)
        emit("W7", "bytes allocated by forward, inverse, tangent, inverse_tangent, propagate_orbit (100 steps) after warm-up, $T", string(al), "(0, 0, 0, 0, 0)", pf(al == (0, 0, 0, 0, 0)))
    end
    lim = maxintfloat(Float64)
    W0 = SMatrix{2, 2, Float64}(1.0, 0.0, 0.0, 1.0)
    sm_up = SVector(6.0, 3.0, lim - 1, 0.0); sn_up = SVector(Float64(pi) / 2, 5.0, 0.0, lim - 1)
    nan_m = all(isnan, forward(WrappedStandardMap(0.0), sm_up))
    nan_n = all(isnan, forward(WrappedStandardMap(3.0), sn_up))
    nan_i = all(isnan, inverse(WrappedStandardMap(0.0), SVector(1.0, 3.0, -(lim - 1), 0.0)))
    e1 = smw_health(() -> propagate_orbit(WrappedStandardMap(0.0), sm_up, W0, 1, 1, ()))
    e2 = smw_health(() -> propagate(WrappedStandardMap(3.0), [sn_up], 1; k = 2, seed = 1))
    okf = forward(WrappedStandardMap(0.0), SVector(6.0, 3.0, lim - 2, 0.0))
    w8 = nan_m && nan_n && nan_i && e1 isa NumericalHealthError && e1.what == "state" && e2 isa NumericalHealthError && e2.what == "state" && all(isfinite, okf) && okf[3] == lim - 1
    emit("W8", "NaN in all components at the count limit (m up, n up, inverse), health error kind 'state' from propagate_orbit and propagate, two below the limit stays finite and exact", w8, "true", pf(w8))
end

# ---------------------------------------------------------------------------------------------------
# Moved from the package suite (D38). The package suite was cut to stay under 3 minutes; the sizes it used before are run here
# with the package seeds (SM_SEED + 1001 ...), so the earlier package-suite assertions stay checked at their earlier size.
# Not repeated here because a section above already runs the same claim at equal or larger size (different seeds):
#   AC3, AC5, AC6 (200 orbits per listed K in "AC1-6"), AC7 to AC13 (same sizes in "AC7-13"), AC14 (1000 per K in "AC14"),
#   AC15 to AC17 (N in {1, 2, 3, 5}, 10 draws, both types in "AC15-17"), W1 (500 / 200), W2 (500 / 200), W5 (200 / 100),
#   W6 (200 / 200) in "W".
section("moved-D38") do
    cF = sm_cases(StableRNG(SM_SEED + 1001), 1000)
    cD = sm_cases(StableRNG(SM_SEED + 1002), 200)
    ad, cs = sm_ac1_ratios(Float64, cF; complexstep = true)
    ad2, _ = sm_ac1_ratios(Double64, cD; complexstep = false)
    ratio_line("AC1", "moved from package suite: max ratio |DT - J_ForwardDiff| / 8 eps (1+|K|), Float64, $(length(cF)) pts", ad)
    ratio_line("AC1", "moved from package suite: max ratio |DT - J_complexstep(1e-30)| / 8 eps (1+|K|), Float64, $(length(cF)) pts", cs)
    ratio_line("AC1", "moved from package suite: max ratio |DT - J_ForwardDiff| / 8 eps (1+|K|), Double64 (eps 2^-100), $(length(cD)) pts", ad2)
    ratio_line("AC2", "moved from package suite: max ratio of |det DT - 1| and |DT' Om DT - Om| to 8 eps (1+|K|), Float64, $(length(cF)) pts",
               sm_ac2_ratios(Float64, sm_cases(StableRNG(SM_SEED + 2001), 1000)))
    ratio_line("AC2", "moved from package suite: same, Double64 (eps 2^-100), $(length(cD)) pts", sm_ac2_ratios(Double64, sm_cases(StableRNG(SM_SEED + 2002), 200)))
    ratio_line("AC4", "moved from package suite: max ratio |DT DT_inv(Tz) - I| / 8 eps (1+|K|)^2, Float64, $(length(cF)) pts",
               sm_ac4_ratio(Float64, sm_cases(StableRNG(SM_SEED + 4001), 1000)))
    ratio_line("AC4", "moved from package suite: same, Double64 (eps 2^-100), $(length(cD)) pts", sm_ac4_ratio(Double64, sm_cases(StableRNG(SM_SEED + 4002), 200)))
end

# ---------------------------------------------------------------------------------------------------
# Timings (loaded machine) for the record
section("timings") do
    for T in (Float64, Double64)
        sys = StandardMap(T(5)); z = SVector(T(1.25), T(-0.5))
        tf = @belapsed forward($sys, $z)
        tt = @belapsed tangent($sys, $z)
        say(@sprintf("# timing (loaded machine) %s: forward %.1f ns, tangent %.1f ns at (1.25, -0.5), K = 5", string(T), 1e9 * tf, 1e9 * tt))
    end
end

# ---------------------------------------------------------------------------------------------------
# Figure: orbit agreement error (Float64 vs 256-bit) against n, with the bound tau(n), K in {1, 5}
section("figure") do
    nmax = 40; m = 1000
    fig = CM.Figure(size = (1100, 460))
    for (j, K) in enumerate((1.0, 5.0))
        rng = StableRNG(SM_SEED + 40700 + j)
        lim = 4 * Float64(pi)
        errs = zeros(m, nmax); taus = zeros(m, nmax); tor = zeros(m, nmax)
        setprecision(BigFloat, SM_BITS) do
            for i in 1:m
                z0 = SVector((2 * rand(rng) - 1) * lim, (2 * rand(rng) - 1) * lim)
                ref = sm_orbit(StandardMap(BigFloat(K)), sm_cast(BigFloat, z0), nmax)
                low = sm_orbit(StandardMap(K), z0, nmax)
                for n in 1:nmax
                    errs[i, n] = sm_absdiff(low[n + 1], ref[n + 1])
                    tor[i, n] = sm_tdist(low[n + 1], ref[n + 1])
                    taus[i, n] = sm_tau(ldexp(1.0, -52), K, ref, n)
                end
            end
        end
        nH = sm_nH(K)
        ax = CM.Axis(fig[1, j]; yscale = log10, xlabel = "step n", ylabel = "max-norm error vs 256-bit orbit",
                  title = "K = $K, $m orbits, Float64 (horizon n_H = $nH)", yminorticksvisible = true)
        ns = 1:nmax
        pos(v) = [x > 0 ? x : NaN for x in v]
        CM.lines!(ax, ns, pos(vec(maximum(errs; dims = 1))); color = :crimson, linewidth = 2, label = "lifted error, ensemble max")
        CM.lines!(ax, ns, pos(vec(median(errs; dims = 1))); color = :crimson, linestyle = :dash, label = "lifted error, ensemble median")
        CM.lines!(ax, ns, pos(vec(maximum(tor; dims = 1))); color = :darkorange, linewidth = 1.5, label = "torus distance, ensemble max")
        CM.lines!(ax, ns, vec(minimum(taus; dims = 1)); color = :black, linewidth = 2, label = "tau(n), ensemble min (bound)")
        CM.lines!(ax, ns, vec(maximum(taus; dims = 1)); color = :gray40, linestyle = :dot, label = "tau(n), ensemble max")
        CM.vlines!(ax, [nH]; color = :steelblue, linestyle = :dashdot, label = "n_H(K)")
        CM.ylims!(ax, 1e-18, 1e12)
        j == 1 && CM.axislegend(ax; position = :lt, labelsize = 11)
        # numbers behind the figure
        okh = all(errs[:, 1:nH] .<= taus[:, 1:nH]) && all(tor[:, 1:nH] .<= taus[:, 1:nH])
        say("# figure K=$K: n <= n_H max over orbits of error/tau = " * fmt(maximum(max.(errs[:, 1:nH], tor[:, 1:nH]) ./ taus[:, 1:nH])) * " (all below tau: $okh); first n where ensemble-max lifted error exceeds ensemble-min tau: " *
            string(something(findfirst(n -> maximum(errs[:, n]) > minimum(taus[:, n]), 1:nmax), 0)))
    end
    CM.save(joinpath(OUT, "standard-map-agreement.png"), fig)
    say("# figure saved: standard-map-agreement.png")
end

say(@sprintf("# total wall time %.1f s", time() - T0))
close(LOG)
