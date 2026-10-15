# Stage 04 verification of module lyap-benettin (catalog L1, Benettin single-vector estimator).
# Run from the repo root:  julia +release -t 2 --project=verify verify/lyap-benettin.jl
# Spec: workspaces/module-build/stages/01-spec/output/lyap-benettin-spec.md (criteria 1 to 19, Amendment D35).
# One printed line per criterion:  AC<n> | quantity | value | bound | PASS/FAIL/INFO/STAGE04-GPU
# Most "value" entries are the worst ratio error / bound (pass when <= 1), as in the package tests.
# Seeds: LB_SEED + 41000 + 50 * criterion + sub-index (the package tests use LB_SEED + 21000 .. 21999, the sibling scripts 44000 and 45000).
# Package-test drivers that hard-code a package seed (lb_ac5_data, lb_ac6_data, lb_regular) are included with a textual patch of the seed
# and the sample size (include_patched checks that every patch target exists), so the measurement code itself is the package's.
using Printf, Random, LinearAlgebra, Statistics
using Chaoslab, StaticArrays, DoubleFloats, StableRNGs
import ChaosTools as CT
import CairoMakie as CM
const HAVE_CUDA = try
    @eval using CUDA
    true
catch
    false
end

const ROOT = dirname(@__DIR__)
const OUT = joinpath(ROOT, "workspaces", "module-build", "stages", "04-verify", "output")
mkpath(OUT)
# LBV_ONLY="AC1,Y1" runs only the named sections and writes the log to LBV_LOG (debugging aid); the default runs everything
const ONLY = filter(!isempty, split(get(ENV, "LBV_ONLY", ""), ","))
const LOG = open(get(ENV, "LBV_LOG", joinpath(OUT, "lyap-benettin-verify-log.txt")), "w")

fmt(x::AbstractFloat) = @sprintf("%.5g", x)
fmt(x) = string(x)
function say(s)
    println(s); println(LOG, s); flush(stdout); flush(LOG)
end
function emit(ac, quantity, value, bound, verdict)
    say(string(ac, " | ", quantity, " | ", fmt(value), " | ", bound, " | ", verdict))
end
pf(ok::Bool) = ok ? "PASS" : "FAIL"
ratio_line(ac, quantity, r) = emit(ac, quantity, r, "<= 1", pf(r <= 1))

function include_patched(path, pairs...)
    txt = read(path, String)
    i = findfirst("\n@testset", txt)
    txt = i === nothing ? txt : txt[1:first(i)]
    for (a, b) in pairs
        occursin(a, txt) || error("patch target not found: $a in $path")
        txt = replace(txt, a => b)
    end
    include_string(Main, txt, path)
end

const TESTDIR = joinpath(ROOT, "test", "lyap-benettin")
include_patched(joinpath(TESTDIR, "lb_helpers.jl"),
    "seed = LB_SEED + LB_OFF + 70" => "seed = LB_SEED + 41350",
    "while length(ok) < 50" => "while length(ok) < 100",
    "keep = ok[1:50]" => "keep = ok[1:100]",
    "last(keep) - 50" => "last(keep) - 100")
include_patched(joinpath(TESTDIR, "lb_ac01_04.jl"))
include_patched(joinpath(TESTDIR, "lb_ac05_07.jl"),
    "seed = LB_SEED + LB_OFF + 50" => "seed = LB_SEED + 41250",
    "T === Double64 ? (N > 1000 ? 6 : 10) : 20" => "20",
    "seed = LB_SEED + LB_OFF + 60" => "seed = LB_SEED + 41300")

vs(ac, sub = 0) = LB_SEED + 41000 + 50 * ac + sub
const KLAB = Dict(:half => "0.5", :two_pi => "2 pi", :ten => "10")

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
say("# lyap-benettin verification, julia $(VERSION), threads $(Threads.nthreads()), base seed $LB_SEED, time $(round(time())) (unix), CUDA loaded: $HAVE_CUDA")
say("# timings in this log were taken on a loaded machine (other agents were running concurrently)")

# ---------------------------------------------------------------------------------------------------
# External oracle: ChaosTools lyapunovspectrum (tangent method, QR of the deviation vectors at every step), lifted standard map.
ct_rule(u, p, n) = (pp = u[2] + p[1] * sin(u[1]); SVector(u[1] + pp, pp))
ct_jac(u, p, n) = (c = p[1] * cos(u[1]); SMatrix{2, 2}(1 + c, c, 1, 1))     # [[1 + K cos x, 1], [K cos x, 1]], state (x, p)
ct_make(K) = CT.TangentDynamicalSystem(CT.DeterministicIteratedMap(ct_rule, SVector(0.1, 0.1), [K]); J = ct_jac, k = 1)
ct_lambda(tds, z0, N) = CT.lyapunovspectrum(tds, N; u0 = SVector{2, Float64}(z0))[1]
const E1 = SMatrix{2, 1, Float64}(1.0, 0.0)          # ChaosTools 3.5.4 restarts every call from the first unit vector (see report)

function ct_lambdas(K, x0, N)
    out = zeros(length(x0))
    Threads.@threads for j in eachindex(x0)
        out[j] = ct_lambda(ct_make(K), x0[j], N)
    end
    return out
end
# Chaoslab Benettin lambda for an ensemble (CPU backend explicitly: CUDA is loaded here, :auto would pick the GPU).
function lb_lambdas(sys, x0, N; seed, tau = 1, W0 = nothing)
    kw = W0 === nothing ? (;) : (; W0 = W0)
    res = propagate(sys, x0, N; k = 1, tau = tau, seed = seed, reducers = (Benettin(),), backend = :cpu, kw...)
    return Float64[f[1].lambda[1] for f in res.finals]
end
# Rule S5: difference of ensemble means against 4 combined standard errors; returns (difference, 4 SE, |difference| / 4 SE).
function s5(a, b)
    d = mean(a) - mean(b)
    se = sqrt(var(a) / length(a) + var(b) / length(b))
    return (d, 4 * se, abs(d) / (4 * se))
end
uniform_starts(seed, m) = lb_starts(StableRNG(seed), m)

# ---------------------------------------------------------------------------------------------------
# AC1: Y5 hyperbolic map
section("AC1") do
    seed = vs(1)
    ns = (20, 100, 10_000)
    for T in LB_TYPES
        epsT = lb_eps(T); nid = 100
        sys = LinearSymplecticMap(SMatrix{2, 2, T}(2, 1, 1, 1))
        x0 = [SVector{2, T}(0, 0) for _ in 1:nid]
        w0s = [lb_w0(seed, j, T) for j in 1:nid]
        wa, wb, excluded, nbad = 0.0, 0.0, 0, 0
        setprecision(BigFloat, LB_BITS) do
            A = (BigFloat(2), BigFloat(1), BigFloat(1), BigFloat(1))
            mu = (3 + sqrt(BigFloat(5))) / 2; phi = (1 + sqrt(BigFloat(5))) / 2; nv = sqrt(phi^2 + 1)
            vp = (phi / nv, 1 / nv); vm = (-1 / nv, phi / nv)
            powers = Dict(n => lb_matpow(A, n) for n in ns)
            for tau in (1, 5, 10)
                G = sum(abs2, Float64.(lb_matpow(A, tau)))
                for n in ns
                    res = lb_run(sys, x0, n; k = 1, tau = tau, seed = seed, reducers = (Benettin(),))
                    P = powers[n]
                    for j in 1:nid
                        w = (lb_big(w0s[j][1]), lb_big(w0s[j][2]))
                        c1 = w[1] * vp[1] + w[2] * vp[2]; c2 = w[1] * vm[1] + w[2] * vm[2]
                        if abs(c1) < 0.05
                            tau == 1 && n == ns[1] && (excluded += 1)
                            continue
                        end
                        y = (P[1] * w[1] + P[2] * w[2], P[3] * w[1] + P[4] * w[2])
                        n0 = log(w[1]^2 + w[2]^2) / 2
                        refa = log(y[1]^2 + y[2]^2) / 2 - n0
                        refb = n * log(mu) + log(abs(c1)) + log(1 + (c2 / c1)^2 * mu^(-4n)) / 2 - n0
                        S, C = lb_nX(lb_fin(res, j))
                        B = lb_BS(epsT, n, G)
                        wa = max(wa, lb_err(S, C, refa) / B); wb = max(wb, lb_err(S, C, refb) / B)
                        nbad += lb_fin(res, j).n == n ? 0 : 1
                    end
                end
            end
        end
        ratio_line("AC1", "(a) max |nX - ln(|A^n w0|/|w0|)_ref| / B_S, tau {1,5,10}, n {20,100,1e4}, $nid ids, $T", wa)
        ratio_line("AC1", "(b) max |nX - closed form E10| / B_S, $T (excluded |c1| < 0.05: $excluded of $nid ids, limit 10; final n wrong: $nbad)", wb)
    end
end

# AC2: rotation and elliptic point
section("AC2") do
    seed = vs(2)
    for T in LB_TYPES
        epsT = lb_eps(T)
        rot = LinearSymplecticMap(SMatrix{2, 2, T}(0, -1, 1, 0))
        nz = 0
        for tau in (1, 3), n in (12, 120, 1200, 9996)
            res = lb_run(rot, [SVector{2, T}(0, 0)], n; k = 1, tau = tau, seed = seed, reducers = (Benettin(),), W0 = [SMatrix{2, 1, T}(1, 0)])
            f = lb_fin(res, 1)
            nz += (f.sum[1] == 0 && f.comp[1] == 0 && f.lambda[1] == 0) ? 0 : 1
        end
        emit("AC2", "(a) rot4, w0 = (1,0), tau {1,3}, n up to 9996: cases with sum, comp or lambda != 0, $T", nz, "0", pf(nz == 0))
        sys = StandardMap(exact(T, -2))
        nid = 100
        x0 = [SVector{2, T}(0, 0) for _ in 1:nid]
        w0s = [lb_w0(seed, j, T) for j in 1:nid]
        bound_odd = log((3 + sqrt(5.0)) / 2)
        Mpow = setprecision(BigFloat, LB_BITS) do
            Mb = (BigFloat(-1), BigFloat(1), BigFloat(-2), BigFloat(1))
            Dict(n => lb_matpow(Mb, n) for n in (1, 2, 3, 4, 10, 11, 100, 101, 1000, 1001))
        end
        wref, weven, wodd = 0.0, 0.0, 0.0
        for (tau, ns) in ((1, (1, 2, 3, 4, 10, 11, 100, 101, 1000, 1001)), (2, (2, 4, 10, 100, 1000)))
            G = tau == 1 ? 7.0 : 2.0
            for n in ns
                res = lb_run(sys, x0, n; k = 1, tau = tau, seed = seed, reducers = (Benettin(),))
                B = lb_BS(epsT, n, G)
                for j in 1:nid
                    S, C = lb_nX(lb_fin(res, j))
                    ref = setprecision(BigFloat, LB_BITS) do
                        P = Mpow[n]; w1, w2 = lb_big(w0s[j][1]), lb_big(w0s[j][2])
                        y1, y2 = P[1] * w1 + P[2] * w2, P[3] * w1 + P[4] * w2
                        (log(y1^2 + y2^2) - log(w1^2 + w2^2)) / 2
                    end
                    wref = max(wref, lb_err(S, C, ref) / B)
                    v = abs(Float64(S + C))
                    iseven(n) ? (weven = max(weven, v / B)) : (wodd = max(wodd, v / (bound_odd + B)))
                end
            end
        end
        ratio_line("AC2", "(b) K = -2 at (0,0), max |nX - ln(|M^n w0|/|w0|)_ref| / B_S, tau {1,2}, $nid ids, $T", wref)
        ratio_line("AC2", "(b) even n: max |nX| / B_S, $T", weven)
        ratio_line("AC2", "(b) odd n: max |nX| / (ln((3+sqrt5)/2) + B_S), $T", wodd)
    end
end

# AC3: K = 0 shear, including the Double64 n = 1e6 case that the package suite does not carry (D38 / test plan deviation 3)
section("AC3") do
    seed = vs(3)
    rng = StableRNG(seed)
    dy() = rand(rng, -2048:2048) / 1024
    starts = [(dy(), dy()) for _ in 1:20]
    for T in LB_TYPES
        epsT = lb_eps(T)
        sys = StandardMap(exact(T, 0))
        wr, wlt, wr6 = 0.0, 0.0, 0.0
        for tau in (1, 10), n in (10, 100, 10_000, 1_000_000)
            m = n == 1_000_000 ? 3 : 20
            x0 = [SVector{2, T}(T(s[1]), T(s[2])) for s in starts[1:m]]
            res = lb_run(sys, x0, n; k = 1, tau = tau, seed = seed, reducers = (Benettin(),))
            B = lb_BS(epsT, n, 2.0 + tau^2)
            for j in 1:m
                S, C = lb_nX(lb_fin(res, j))
                w = lb_w0(seed, j, T)
                ref = setprecision(BigFloat, LB_BITS) do
                    a, b = lb_big(w[1]), lb_big(w[2])
                    (log((a + n * b)^2 + b^2) - log(a^2 + b^2)) / 2
                end
                r = lb_err(S, C, ref) / B
                n == 1_000_000 ? (wr6 = max(wr6, r)) : (wr = max(wr, r))
                wlt = max(wlt, Float64(S + C) / log(n + 1.0))
            end
        end
        w01 = 0.0
        for tau in (1, 10)
            res = lb_run(sys, [SVector{2, T}(T(0.25), T(0.5))], 10_000; k = 1, tau = tau, seed = seed, reducers = (Benettin(),), W0 = [SMatrix{2, 1, T}(0, 1)])
            S, C = lb_nX(lb_fin(res, 1))
            ref = setprecision(BigFloat, LB_BITS) do
                log(BigFloat(10_000)^2 + 1) / 2
            end
            w01 = max(w01, lb_err(S, C, ref) / (8 * 10_000 * epsT))
        end
        ratio_line("AC3", "max |nX - ln(|w_n|/|w0|)| / B_S, tau {1,10}, n {10,100,1e4} x 20 ids, $T", wr)
        ratio_line("AC3", "same, n = 1e6, 3 ids, $T $(T === Double64 ? "(moved from the package suite)" : "")", wr6)
        ratio_line("AC3", "max nX / ln(n+1) (strict bound < 1), all n, $T", wlt)
        ratio_line("AC3", "w0 = (0,1), n = 1e4: |S_n - 0.5 ln(n^2+1)| / (8 n eps), $T", w01)
    end
end

# AC4: hyperbolic fixed point
section("AC4") do
    seed = vs(4)
    for T in LB_TYPES
        epsT = lb_eps(T)
        x0 = [SVector{2, T}(0, 0) for _ in 1:100]
        w0s = [lb_w0(seed, j, T) for j in 1:100]
        worst = Dict(K => 0.0 for K in (1, 2, 6))
        for K in (1, 2, 6)
            sys = StandardMap(exact(T, K))
            Mf = SMatrix{2, 2, Float64}(1 + K, K, 1, 1)
            for tau in (1, 5), n in (20, 100, 1000)
                B = lb_BS(epsT, n, sum(abs2, Mf^tau))
                res = lb_run(sys, x0, n; k = 1, tau = tau, seed = seed, reducers = (Benettin(),))
                setprecision(BigFloat, LB_BITS) do
                    Kb = BigFloat(K); s = sqrt(Kb * (Kb + 4)); mu = (2 + Kb + s) / 2
                    vu = (one(Kb), (s - Kb) / 2); vsb = (one(Kb), -(Kb + s) / 2)
                    cross(a, b) = a[1] * b[2] - a[2] * b[1]
                    for j in 1:100
                        w = (lb_big(w0s[j][1]), lb_big(w0s[j][2]))
                        c1 = cross(w, vsb) / cross(vu, vsb); c2 = cross(vu, w) / cross(vu, vsb)
                        nrm2 = c1^2 * mu^(2n) * (vu[1]^2 + vu[2]^2) + 2 * c1 * c2 * (vu[1] * vsb[1] + vu[2] * vsb[2]) + c2^2 * mu^(-2n) * (vsb[1]^2 + vsb[2]^2)
                        ref = (log(nrm2) - log(w[1]^2 + w[2]^2)) / 2
                        S, C = lb_nX(lb_fin(res, j))
                        worst[K] = max(worst[K], lb_err(S, C, ref) / B)
                    end
                end
            end
        end
        for K in (1, 2, 6)
            ratio_line("AC4", "max |nX - closed form E12| / B_S, K = $K, tau {1,5}, n {20,100,1e3}, 100 ids, $T", worst[K])
        end
    end
end

# AC5, AC6, AC18: the package drivers with fresh seeds; Double64 with 20 starts at every n (moved)
section("AC5-6-18") do
    rows = lb_ac5_data()
    for T in LB_TYPES, (label, N) in LB_AC5_CASES
        rs = filter(r -> r.T === T && r.label == label, rows)
        w = maximum(r.ratio for r in rs)
        g = count(r -> r.G > 1e4, rs)
        okx = all(r -> r.xeq && r.nsnaps == div(r.N, 100), rs)
        emit("AC5", "max |nX - ln(|Y_n w0|/|w0|)_ref| / B_S at every n multiple of 100, K = $(KLAB[label]), n = $N, $T, $(length(rs)) runs (20 starts x tau {1,4}); G_first > 1e4: $g; orbit == forward and all handoffs read: $okx",
             w, "<= 1", pf(w <= 1 && okx))
    end
    for T in LB_TYPES, (label, N) in LB_AC5_CASES
        rs = filter(r -> r.T === T && r.label == label, rows)
        w = maximum(r.rv for r in rs)
        ratio_line("AC18", "max alpha_j range violation ratio (E14), AC5 runs K = $(KLAB[label]), n = $N, $T", w)
    end
    rows6 = lb_ac6_data()
    for T in LB_TYPES
        rs = filter(r -> r.T === T, rows6)
        w = maximum(r.diff / r.bound for r in rs)
        ratio_line("AC6", "max |n X_tau - n X_1| / (B_S^tau + B_S^1), tau {2,5,10}, K {0.5, 2pi}, n {100,1e3}, $(length(rs)) comparisons, $T", w)
        ratio_line("AC18", "max alpha_j range violation ratio (E14), AC6 runs, $T", maximum(r.rv for r in rs))
    end
end

# AC7: full BigFloat orbit, regular class (screened, fresh seed)
section("AC7") do
    ids, starts, rejected, seed = lb_regular()
    say("# S7 screening (fresh seed $seed): first 100 regular starts, $rejected candidates replaced")
    for T in LB_TYPES, n in (10, 100)
        epsT = lb_eps(T)
        sys = StandardMap(exact(T, 0.5))
        excluded, w, wprint = 0, 0.0, 0.0
        for (id, s) in zip(ids[1:50], starts[1:50])
            W0 = initial_deviations(seed, id, Val(2), Val(1), T)
            f = lb_final(sys, lb_cast(T, s), W0, n, 1)
            err, rho, M = lb_ref_full(0.5, s, W0[:, 1], n, f.sum[1], f.comp[1])
            if rho > 50
                excluded += 1
                continue
            end
            w = max(w, err / (32 * 0.5 * M * n^2 * epsT * rho))
            wprint = max(wprint, err / (32 * 0.5 * n^2 * epsT * rho))      # the printed form without M_state (superseded, informational)
        end
        ratio_line("AC7", "max |nX_T - nX_ref| / (32 K M_state n^2 eps rho_n), n = $n, 50 starts, $T (rho_n > 50 excluded: $excluded, limit 5)", w)
        emit("AC7", "same without the M_state factor (printed form before Amendment D35, informational), n = $n, $T", wprint, "<= 1", "INFO")
    end
end

# AC8: norm evaluation
section("AC8") do
    rng = StableRNG(vs(8))
    for T in LB_TYPES
        epsT = lb_eps(T)
        relerr(v, ref) = setprecision(BigFloat, LB_BITS) do
            Float64(abs(lb_big(safe_norm(v)) - ref) / ref)
        end
        w = 0.0; nonfinite = 0
        for k in (600, -600, 0, 1000, -1000)
            c = ldexp(one(T), k)
            v = SVector{2, T}(3c, 4c)
            nonfinite += isfinite(safe_norm(v)) ? 0 : 1
            w = max(w, relerr(v, 5 * BigFloat(2)^k) / (2 * epsT))
        end
        for _ in 1:2000
            d = rand(rng, (2, 4)); e = rand(rng, -900:900)
            v = d == 2 ? SVector{2, T}(ldexp(T(2 * rand(rng) - 1), e), ldexp(T(2 * rand(rng) - 1), e)) :
                SVector{4, T}(ntuple(_ -> ldexp(T(2 * rand(rng) - 1), e), 4))
            ref = setprecision(BigFloat, LB_BITS) do
                sqrt(sum(x -> lb_big(x)^2, v))
            end
            iszero(ref) && continue
            w = max(w, relerr(v, ref) / (2 * epsT))
        end
        z = safe_norm(SVector{2, T}(0, 0)) == 0 && !isfinite(safe_norm(SVector{2, T}(T(NaN), 1))) && !isfinite(safe_norm(SVector{2, T}(T(Inf), 1)))
        ratio_line("AC8", "max relative error of safe_norm / 2 eps, c (3,4) c in 2^{600,-600,0,1000,-1000} and 2000 random d {2,4} at scale 2^e, |e| <= 900, $T (non-finite results: $nonfinite)", w)
        emit("AC8", "zero gives 0, NaN and Inf give a non-finite norm, $T", z, "true", pf(z))
    end
end

# AC9: health (checks counted)
section("AC9") do
    seed = vs(9)
    x0s(flagged) = [SVector{3, Float64}(0, 0, j in flagged ? 1 : 0) for j in 1:6]
    zero_map, nan_map = lb_Fault(:zero, 0), lb_Fault(:nan, 3)
    checks = Bool[]
    ok = propagate(zero_map, x0s(()), 12; k = 1, tau = 3, seed = seed, reducers = (Benettin(),), chunk = 2, backend = :cpu)
    push!(checks, length(ok.finals) == 6)
    res = nothing
    err = lb_health_error(() -> (res = propagate(zero_map, x0s((4,)), 12; k = 1, tau = 3, seed = seed, reducers = (Benettin(),), chunk = 2, backend = :cpu)))
    append!(checks, [err isa NumericalHealthError, res === nothing, err.what == "deviation matrix", err.orbit == 4, err.step == 3, err.backend === :cpu, err.column == 1, err.value == 0.0])
    res = nothing
    err = lb_health_error(() -> (res = propagate(zero_map, x0s((5, 2)), 12; k = 1, tau = 3, seed = seed, reducers = (Benettin(),), chunk = 2, backend = :cpu)))
    append!(checks, [err isa NumericalHealthError && err.orbit == 2 && err.step == 3, res === nothing])
    res = nothing
    err = lb_health_error(() -> (res = propagate(nan_map, x0s((3,)), 20; k = 1, tau = 5, seed = seed, reducers = (Benettin(),), chunk = 2, backend = :cpu)))
    append!(checks, [err isa NumericalHealthError, res === nothing, err.what == "deviation matrix", err.orbit == 3, err.step == 5, err.column == 0, err.threshold == Inf, !isfinite(err.value)])
    for T in LB_TYPES
        x = SVector{3, T}(0, 0, 1)
        W0 = initial_deviations(seed, 1, Val(2), Val(1), T)
        e = lb_health_error(() -> propagate_orbit(zero_map, x, W0, 12, 3, (Benettin(),); orbit = 4))
        append!(checks, [e isa NumericalHealthError && e.what == "deviation matrix" && e.orbit == 4 && e.step == 3, e.backend === :cpu, e.column == 1, e.value == 0.0])
        e = lb_health_error(() -> propagate_orbit(nan_map, x, W0, 20, 5, (Benettin(),); orbit = 3))
        append!(checks, [e isa NumericalHealthError && e.what == "deviation matrix" && e.orbit == 3 && e.step == 5, e.column == 0, e.threshold == Inf, !isfinite(e.value)])
    end
    ex = try
        propagate(lb_Boom(), [SVector(0.0, 0.0)], 4; k = 2, tau = 2, seed = seed, reducers = (Benettin(),), backend = :cpu)
        nothing
    catch e
        e
    end
    push!(checks, ex isa ArgumentError && occursin("Benettin needs k = 1, got k = 2", ex.msg))
    emit("AC9", "health: zero tangent (kind, orbit, step, backend, column, value, no result stored, lowest id first), NaN tangent (column 0, threshold Inf), Double64 per-orbit body, k != 1 ArgumentError: checks passed / total", string(count(checks), "/", length(checks)), "all", pf(all(checks)))
end

# AC10: orbit and execution invariance; Double64 thread invariance is the moved part (1024 orbits, runtime)
section("AC10") do
    seed = vs(10)
    ok_a, tot_a = 0, 0
    for tau in (1, 4)
        sys = StandardMap(lb_K(Float64, :two_pi))
        mon = ResidualMonitor(sys, LB_DET_TOL, Inf)
        for (j, s) in enumerate(lb_starts(StableRNG(seed + 1), 5))
            x0 = lb_cast(Float64, s)
            W0 = initial_deviations(seed, j, Val(2), Val(1), Float64)
            x_a, _, st_a = propagate_orbit(sys, x0, W0, 400, tau, (lb_Trace(),); orbit = j)
            x_b, _, st_b = propagate_orbit(sys, x0, W0, 400, tau, (Benettin(), lb_Trace()); orbit = j)
            x_c, _, st_c = propagate_orbit(sys, x0, W0, 400, tau, (Benettin(), mon, lb_Trace()); orbit = j)
            tr_a = reducer_final(lb_Trace(), st_a[1]).xs; tr_b = reducer_final(lb_Trace(), st_b[2]).xs; tr_c = reducer_final(lb_Trace(), st_c[3]).xs
            tot_a += 1
            ok_a += (length(tr_a) == 400 ÷ tau && isequal(tr_a, tr_b) && isequal(tr_a, tr_c) && x_a === x_b && x_a === x_c) ? 1 : 0
        end
    end
    emit("AC10", "(a) orbit points at every handoff bitwise equal for chains {}, {L1}, {L1, det}: cases equal / total (5 starts x tau {1,4})", string(ok_a, "/", tot_a), "all", pf(ok_a == tot_a))
    nthr = Threads.nthreads()
    for T in LB_TYPES
        m, N = T === Float64 ? 4096 : 1024, 1000
        okb, totb = 0, 0
        for tau in (1, 4)
            sys = StandardMap(lb_K(T, :two_pi))
            reducers = (Benettin(), ResidualMonitor(sys, LB_DET_TOL, Inf))
            x0 = [lb_cast(T, s) for s in lb_starts(StableRNG(seed + 2), m)]
            W0 = [initial_deviations(seed, j, Val(2), Val(1), T) for j in 1:m]
            serial = [propagate_orbit(sys, x0[j], W0[j], N, tau, reducers; orbit = j) for j in 1:m]
            xs = [s[1] for s in serial]; Ws = [s[2] for s in serial]; fs = [map(reducer_final, reducers, s[3]) for s in serial]
            for chunk in (64, 1000)
                res = propagate(sys, x0, N; k = 1, tau = tau, seed = seed, reducers = reducers, chunk = chunk, W0 = W0, backend = :cpu)
                totb += 1
                okb += (isequal(res.x, xs) && isequal(res.W, Ws) && isequal(res.finals, fs) && res.provenance.estimator == "benettin") ? 1 : 0
            end
        end
        emit("AC10", "(b) propagate ($nthr threads, chunks 64 and 1000) bitwise equal to the serial map of propagate_orbit (x, W, finals, provenance), $m orbits, K = 2 pi, n = 1e3, tau {1,4}, $T $(T === Double64 ? "(moved from the package suite)" : ""): equal / total", string(okb, "/", totb), "all", pf(okb == totb))
    end
    for T in LB_TYPES, tau in (1, 4)
        sys = StandardMap(lb_K(T, :two_pi))
        reducers = (Benettin(), ResidualMonitor(sys, LB_DET_TOL, Inf))
        b = lb_alloc_bytes(sys, SVector{2, T}(1, 2), initial_deviations(seed, 1, Val(2), Val(1), T), 10_000, tau, reducers)
        emit("AC10", "(c) bytes allocated per orbit after warm-up, n = 1e4, chain {L1, det}, tau = $tau, $T", b, "0", pf(b == 0))
    end
end

# AC11: wrapped vs lifted
section("AC11") do
    seed = vs(11)
    rng = StableRNG(seed)
    dy() = rand(rng, 0:2048) / 1024
    for T in LB_TYPES
        okc, totc = 0, 0
        for tau in (1, 3)
            zs = [SVector{2, T}(T(dy()), T(dy())) for _ in 1:50]
            resL = lb_run(StandardMap(exact(T, 0)), zs, 300; k = 1, tau = tau, seed = seed)
            resW = lb_run(WrappedStandardMap(exact(T, 0)), torus_state.(zs), 300; k = 1, tau = tau, seed = seed)
            okc += count(j -> isequal(resL.finals[j][1], resW.finals[j][1]), 1:50); totc += 50
        end
        emit("AC11", "(a) K = 0, dyadic starts: wrapped and lifted finals bitwise identical, $T: identical / total", string(okc, "/", totc), "all", pf(okc == totc))
    end
    ids, starts, rejected, rseed = lb_regular()
    for T in LB_TYPES, n in (10, 100)
        W0 = [initial_deviations(rseed, id, Val(2), Val(1), T) for id in ids[1:50]]
        st = [lb_cast(T, s) for s in starts[1:50]]
        resL = lb_run(StandardMap(exact(T, 0.5)), st, n; k = 1, tau = 1, seed = seed, W0 = W0)
        resW = lb_run(WrappedStandardMap(exact(T, 0.5)), torus_state.(st), n; k = 1, tau = 1, seed = seed, W0 = W0)
        excluded, w = 0, 0.0
        for j in 1:50
            fL, fW = resL.finals[j][1], resW.finals[j][1]
            _, rho, M = lb_ref_full(0.5, starts[j], Float64.(W0[j][:, 1]), n, 0.0, 0.0)
            if rho > 50
                excluded += 1
                continue
            end
            diff = Float64(abs((fW.sum[1] + fW.comp[1]) - (fL.sum[1] + fL.comp[1])))
            w = max(w, diff / (64 * 0.5 * M * n^2 * lb_eps(T) * rho))
        end
        ratio_line("AC11", "(b) max |nX_wrapped - nX_lifted| / (64 K M_state n^2 eps rho_n), K = 0.5 regular class, n = $n, 50 starts, $T (rho_n > 50 excluded: $excluded, limit 5)", w)
    end
end

# AC17: regular-orbit decay
section("AC17") do
    ids, starts, rejected, seed = lb_regular()
    ns = unique(round.(Int, 10 .^ range(4, 6; length = 21)))
    sys = StandardMap(0.5)
    good, slopes = 0, Float64[]
    for (id, s) in zip(ids[1:20], starts[1:20])
        probe = lb_Probe(Benettin(), ns; alphas = false)
        W0 = initial_deviations(seed, id, Val(2), Val(1), Float64)
        _, _, st = propagate_orbit(sys, SVector(s[1], s[2]), W0, 1_000_000, 1, (probe,); orbit = id)
        fin = reducer_final(probe, st[1])
        X = [(S + C) / n for (n, S, C) in fin.snaps]
        all(>(0), X) || (push!(slopes, NaN); continue)
        u = log.(Float64[n for (n, _, _) in fin.snaps]); v = log.(X)
        slope = sum((u .- mean(u)) .* (v .- mean(v))) / sum((u .- mean(u)) .^ 2)
        push!(slopes, slope)
        -1.1 <= slope <= -0.8 && (good += 1)
    end
    emit("AC17", "starts with slope of ln X vs ln n over n in [1e4,1e6] in [-1.1,-0.8] (of 20, K = 0.5 regular class, n up to 1e6); slope min / median / max = $(fmt(minimum(filter(!isnan, slopes)))) / $(fmt(median(filter(!isnan, slopes)))) / $(fmt(maximum(filter(!isnan, slopes))))", good, ">= 16", pf(good >= 16))
end

# AC19: symplectic dimension 2
section("AC19") do
    seed = vs(19)
    for (ci, (label, N)) in enumerate(LB_AC5_CASES), T in LB_TYPES, tau in (1, 4)
        K = lb_K(T, label)
        sys = StandardMap(K)
        x0 = [lb_cast(T, s) for s in lb_starts(StableRNG(seed + ci), 20)]
        res = lb_run(sys, x0, N; k = 1, tau = tau, seed = seed, reducers = (Benettin(), ResidualMonitor(sys, LB_DET_TOL, Inf)), probe = false)
        dmax = maximum(f[2].det_max for f in res.finals)
        bound = LB_DET_TOL * effective_eps(T) * (2 + K)
        emit("AC19", "(a) max det residual / threshold 4 eps (2 + K) over 20 starts, K = $(KLAB[label]), n = $N, tau = $tau, $T (nonzero: $(dmax > 0))", Float64(dmax / bound), "<= 1", pf(dmax <= bound && dmax > 0))
    end
    for T in LB_TYPES
        A0 = SMatrix{2, 2, T}(2, 1, 1, 1)
        theta = Float64(det_threshold(LB_DET_TOL, T, A0))
        above = [exp2(j) for j in -130:20 if 2theta <= exp2(j) < 4theta]
        below = [exp2(j) for j in -130:20 if theta / 4 < exp2(j) <= theta / 2]
        runf(delta, N) = begin
            sys = LinearSymplecticMap(SMatrix{2, 2, T}(2 + T(delta), 1, 1, 1))
            mon = ResidualMonitor(sys, LB_DET_TOL, Inf)
            lb_health_error(() -> lb_run(sys, [SVector{2, T}(0, 0)], N; k = 1, tau = 1, seed = seed, reducers = (Benettin(), mon), W0 = [SMatrix{2, 1, T}(1, 0)], probe = false))
        end
        raised = count(delta -> (e = runf(delta, 1); e isa NumericalHealthError && e.what == "det" && e.value == delta && e.orbit == 1 && e.step == 1), above)
        quiet = count(delta -> runf(delta, 100) === nothing, below)
        emit("AC19", "(b) injected det fault: raised for delta in [2 theta, 4 theta) (theta = $(fmt(theta))) / quiet for delta in (theta/4, theta/2], $T", string(raised, "/", length(above), " ; ", quiet, "/", length(below)), "all", pf(raised == length(above) && quiet == length(below)))
    end
end

# ---------------------------------------------------------------------------------------------------
# AC12: precision agreement of the series (Float64 vs Double64 vs BigFloat(256)). Lifted state, CPU backend.
section("AC12") do
    for K in (2pi, 5.0)
        x0 = uniform_starts(vs(12, 1), 1000)
        sysF = StandardMap(K); sysD = StandardMap(exact(Double64, K))
        a = lb_lambdas(sysF, x0, 10_000; seed = vs(12))
        b = lb_lambdas(sysD, [lb_cast(Double64, z) for z in x0], 10_000; seed = vs(12))
        d, se4, r = s5(a, b)
        emit("AC12", "(a) Float64 vs Double64 (lifted), K = $(fmt(K)), 1000 starts, n = 1e4: mean X $(fmt(mean(a))) vs $(fmt(mean(b))), difference $(fmt(d)), 4 combined SE $(fmt(se4)); value = |difference| / 4 SE", r, "<= 1", pf(r <= 1))
    end
    x0 = uniform_starts(vs(12, 2), 100)
    Kt(T) = lb_K(T, :two_pi)
    res = Dict{Any, Vector{Float64}}()
    for T in (Float64, Double64, BigFloat)
        setprecision(BigFloat, LB_BITS) do
            K = T === BigFloat ? exact(BigFloat, pi) * 2 : Kt(T)
            sys = StandardMap(K)
            W0 = [initial_deviations(vs(12), j, Val(2), Val(1), T) for j in 1:100]
            r = lb_run(sys, [lb_cast(T, z) for z in x0], 1000; k = 1, tau = 1, seed = vs(12), W0 = W0)
            res[T] = Float64[f[1].lambda[1] for f in r.finals]
        end
    end
    for (A, B) in ((Float64, BigFloat), (Double64, BigFloat))
        d, se4, r = s5(res[A], res[B])
        emit("AC12", "(b) $A vs BigFloat(256) (lifted), K = 2 pi, 100 starts, n = 1e3: mean X $(fmt(mean(res[A]))) vs $(fmt(mean(res[B]))), difference $(fmt(d)), 4 SE $(fmt(se4)); value = |difference| / 4 SE", r, "<= 1", pf(r <= 1))
    end
    ids, starts, rejected, seed = lb_regular()
    W0F = [initial_deviations(seed, id, Val(2), Val(1), Float64) for id in ids]
    W0D = [initial_deviations(seed, id, Val(2), Val(1), Double64) for id in ids]
    a = lb_lambdas(StandardMap(0.5), [SVector(s[1], s[2]) for s in starts], 100_000; seed = seed, W0 = W0F)
    b = lb_lambdas(StandardMap(exact(Double64, 0.5)), [lb_cast(Double64, s) for s in starts], 100_000; seed = seed, W0 = W0D)
    m = maximum(abs.(a .- b))
    emit("AC12", "(c) regular class K = 0.5 (S7 screening, $rejected candidates replaced), 100 starts, n = 1e5 (lifted): max |X_64 - X_D64| (X values $(fmt(minimum(a))) .. $(fmt(maximum(a))))", m, "<= 1e-4", pf(m <= 1e-4))
end

# AC13 to AC15: ChaosTools oracle
section("AC13-15") do
    # orbit identity: is the ChaosTools lifted orbit bitwise the Chaoslab orbit?
    for K in (0.5, 2pi, 5.0, 10.0)
        sys = StandardMap(K); same, tot = 0, 0
        for z0 in uniform_starts(vs(13, 9), 10)
            ds = CT.DeterministicIteratedMap(ct_rule, SVector{2, Float64}(z0), [K])
            z = SVector{2, Float64}(z0)
            for _ in 1:1000
                CT.step!(ds); z = forward(sys, z)
                same += CT.current_state(ds) == z ? 1 : 0; tot += 1
            end
        end
        emit("AC13", "ChaosTools lifted orbit bitwise equal to Chaoslab forward, K = $(fmt(K)), 10 starts x 1000 steps: equal / total", string(same, "/", tot), "all", "INFO")
    end
    # AC13 loose: regular class K = 0.5, 20 starts, n = 1e5
    ids, starts, rejected, seed = lb_regular()
    st = [SVector(s[1], s[2]) for s in starts[1:20]]
    W0 = [initial_deviations(seed, id, Val(2), Val(1), Float64) for id in ids[1:20]]
    lc = lb_lambdas(StandardMap(0.5), st, 100_000; seed = seed, W0 = W0)
    lt = ct_lambdas(0.5, st, 100_000)
    m = maximum(abs.(lc .- lt))
    emit("AC13", "regular class K = 0.5 (S7), 20 starts, n = 1e5: max |X Chaoslab (seeded w0) - lambda_1 ChaosTools| (X range $(fmt(minimum(lc))) .. $(fmt(maximum(lc))))", m, "<= 1e-4", pf(m <= 1e-4))
    le = lb_lambdas(StandardMap(0.5), st, 100_000; seed = seed, W0 = fill(E1, 20))
    emit("AC13", "same with the identical initial vector e1 on both sides: max |X - lambda_1 ChaosTools|", maximum(abs.(le .- lt)), "<= 1e-4", pf(maximum(abs.(le .- lt)) <= 1e-4))
    # pathwise line: identical w0 = e1, lifted rule, K = 0.5, n = 100, 50 starts, tau = 1
    w, excluded = 0.0, 0
    for (id, s) in zip(ids[1:50], starts[1:50])
        z0 = SVector(s[1], s[2]); n = 100
        f = lb_final(StandardMap(0.5), z0, SMatrix{2, 1, Float64}(1.0, 0.0), n, 1)
        _, rho, M = lb_ref_full(0.5, s, [1.0, 0.0], n, f.sum[1], f.comp[1])
        if rho > 50
            excluded += 1
            continue
        end
        lam = ct_lambda(ct_make(0.5), z0, n)
        w = max(w, abs(Float64(f.sum[1] + f.comp[1]) - n * lam) / (64 * 0.5 * M * n^2 * eps(Float64) * rho))
    end
    ratio_line("AC13", "pathwise, w0 = e1 on both sides, K = 0.5 (S7), n = 100, 50 starts, tau = 1: max |nX - n lambda_ChaosTools| / (64 K M_state n^2 eps rho_n) (rho_n > 50 excluded: $excluded, limit 5)", w)
    # pathwise, chaotic: orbits are bitwise equal, so compare sums at n = 1e3 in units of 8 eps Psi_n
    wps = 0.0
    sys = StandardMap(2pi)
    for (j, s) in enumerate(uniform_starts(vs(13, 3), 20))
        z0 = SVector(s[1], s[2]); n = 1000
        f = lb_final(sys, z0, SMatrix{2, 1, Float64}(1.0, 0.0), n, 1)
        psi = lb_Psi(sys, z0, [1.0, 0.0], n)
        lam = ct_lambda(ct_make(2pi), z0, n)
        wps = max(wps, abs(Float64(f.sum[1] + f.comp[1]) - n * lam) / (8 * eps(Float64) * psi))
    end
    emit("AC13", "pathwise, w0 = e1, K = 2 pi, n = 1e3, 20 starts: max |nX - n lambda_ChaosTools| / (8 eps Psi_n) (two estimators, so up to 2 expected; informational)", wps, "<= 2", "INFO")
    # AC14 and AC15: chaotic ensembles
    for (ac, K, m, N) in (("AC14", 2pi, 1000, 10_000), ("AC15", 5.0, 200, 10_000), ("AC15", 10.0, 200, 10_000))
        x0 = uniform_starts(vs(14, round(Int, K)), m)
        a = lb_lambdas(StandardMap(K), x0, N; seed = vs(14))
        c = ct_lambdas(K, x0, N)
        d, se4, r = s5(a, c)
        emit(ac, "chaotic ensemble K = $(fmt(K)), $m starts, n = $N: mean X Chaoslab $(fmt(mean(a))) vs ChaosTools $(fmt(mean(c))), difference $(fmt(d)), 4 combined SE $(fmt(se4)); value = |difference| / 4 SE", r, "<= 1", pf(r <= 1))
        e = lb_lambdas(StandardMap(K), x0, N; seed = vs(14), W0 = fill(E1, m))
        emit(ac, "same starts, w0 = e1 on both sides: max per-start |X - lambda_ChaosTools| (chaotic orbits are bitwise equal; informational)", maximum(abs.(e .- c)), "-", "INFO")
    end
end

# AC16: CPU vs GPU, same tolerances as verify/gpu.jl AC12c and AC12a (tangent-core AC12, S6); gpu.jl itself has no Benettin line
section("AC16") do
    if !HAVE_CUDA
        emit("AC16", "CPU vs GPU", "-", "tangent-core AC12", "STAGE04-GPU (CUDA not loadable here)")
        return
    end
    ongpu(args...; kw...) = (r = propagate(args...; backend = :gpu, kw...); (r, r.provenance.backend === :gpu ? nothing : r.provenance.backend_fallback_reason))
    Xs(res) = Float64[f[1].lambda[1] for f in res.finals]
    for T in LB_TYPES
        # chaotic class K = 10, 1024 orbits, n = 1e4, tau = 10
        rng = StableRNG(vs(16))
        x0 = [lb_cast(T, SVector(2pi * rand(rng), 2pi * rand(rng))) for _ in 1:1024]
        sys = StandardMap(exact(T, 10))
        kw = (k = 1, tau = 10, seed = vs(16), reducers = (Benettin(),))
        cpu = propagate(sys, x0, 10_000; kw..., backend = :cpu)
        gpu, why = ongpu(sys, x0, 10_000; kw...)
        if why !== nothing
            emit("AC16", "chaotic class K = 10, 1024 orbits, n = 1e4, tau = 10, $T", "-", "6 sigma; 4 SE", "FALLBACK($why)")
        else
            sc, sg = Xs(cpu), Xs(gpu); sig = std(sc)
            maxd = maximum(abs.(sc .- sg)) / sig
            dm = abs(mean(sc) - mean(sg)) / sqrt((var(sc) + var(sg)) / length(sc))
            emit("AC16", "chaotic class K = 10, 1024 orbits, n = 1e4, tau = 10, $T: max per-orbit |X_cpu - X_gpu| / sigma (sigma $(fmt(sig))); |difference of means| / SE", string(fmt(maxd), "; ", fmt(dm)), "6; 4", pf(maxd <= 6 && dm <= 4))
        end
        T === Float64 || continue
        # regular class K = 0.5, 1024 orbits, n = 4e4, tau = 10, against a Double64 CPU reference (D33 symmetric form)
        rng = StableRNG(vs(16, 1))
        x0 = [SVector(2pi * rand(rng), 2pi * rand(rng)) for _ in 1:1024]
        sys = StandardMap(0.5)
        cpu = propagate(sys, x0, 40_000; kw..., backend = :cpu)
        gpu, why = ongpu(sys, x0, 40_000; kw...)
        if why !== nothing
            emit("AC16", "regular class K = 0.5, 1024 orbits, n = 4e4, tau = 10, $T", "-", "count diff <= 3; gap <= 2e-6", "FALLBACK($why)")
        else
            ref = propagate(StandardMap(Double64(0.5)), [SVector{2, Double64}(z) for z in x0], 40_000; kw..., backend = :cpu)
            sc, sg, sr = Xs(cpu), Xs(gpu), Xs(ref)
            dcr, dgr, dcg = abs.(sc .- sr), abs.(sg .- sr), abs.(sc .- sg)
            ncpu, ngpu = count(>(1e-6), dcr), count(>(1e-6), dgr)
            both = (dcr .<= 1e-6) .& (dgr .<= 1e-6)
            gap = maximum(dcg[both])
            emit("AC16", "regular class K = 0.5, 1024 orbits, n = 4e4, tau = 10, $T: orbits above 1e-6 vs the Double64 reference, CPU / GPU; max CPU-GPU gap over the $(count(both)) orbits where both agree", string(ncpu, " / ", ngpu, "; ", fmt(gap)), "count diff <= 3; gap <= 2e-6", pf(abs(ncpu - ngpu) <= 3 && gap <= 2e-6))
        end
    end
end

# ---------------------------------------------------------------------------------------------------
# Y1: large-K known answer (loose sanity lines, K18: reported, never pass/fail). Float64, lifted, CPU, tau = 1.
section("Y1") do
    N = 1_000_000
    for K in (pi / 2, 2pi, 10.0, 20.0)
        x0 = uniform_starts(vs(1, 700 + round(Int, 10K)), 100)
        a = lb_lambdas(StandardMap(K), x0, N; seed = vs(1, 700))
        se = std(a) / sqrt(length(a))
        chirikov = log(K / 2); shev = log(K / 2) + 1 / K^2
        ref = K >= 4.4 ? "ln(K/2) = $(fmt(chirikov)), ln(K/2) + 1/K^2 = $(fmt(shev))" : "ln(K/2) = $(fmt(chirikov)) (not asymptotic at this K)"
        say(string("Y1 | ensemble K = $(fmt(K)), 100 uniform starts, n = 1e6: mean X $(fmt(mean(a))) +- $(fmt(se)) (SE), median $(fmt(median(a))), min $(fmt(minimum(a))), max $(fmt(maximum(a))); ", ref,
                   "; mean - ln(K/2) = $(fmt(mean(a) - chirikov))", K >= 4.4 ? ", mean - (ln(K/2) + 1/K^2) = $(fmt(mean(a) - shev)) = $(fmt((mean(a) - shev) / se)) SE" : "",
                   " | - | INFO"))
    end
    # single orbit at the E9-converted start (x, p) = (1e-3, -1e-3), n = 1e7, against the SalaArtuso2015 benchmark values
    for (K, ref) in ((pi / 2, 0.298), (2pi, 1.172))
        a = lb_lambdas(StandardMap(K), [SVector(1e-3, -1e-3)], 10_000_000; seed = vs(1, 800))
        say("Y1 | single orbit K = $(fmt(K)), start (x, p) = (1e-3, -1e-3) (E9 of (x,y) = (1,2) 1e-3), n = 1e7: X = $(fmt(a[1])) vs SalaArtuso2015 $ref (their n = 1e9); difference $(fmt(a[1] - ref)) | - | INFO")
    end
end

# ---------------------------------------------------------------------------------------------------
# PR-trig: Double64 sin, cos error for large |x| (precision spec Amendment D35): absolute error against BigFloat(256), by decade of |x|
section("PR-trig") do
    rng = StableRNG(vs(30))
    epsEff = Float64(effective_eps(Double64))
    worst_abs, worst_c = 0.0, 0.0
    for e in 0:7
        lo, hi = 10.0^e, 10.0^(e + 1)
        mabs, mrel, mc = 0.0, 0.0, 0.0
        xmax = 0.0
        setprecision(BigFloat, LB_BITS) do
            for _ in 1:1500
                u = lo * (hi / lo)^rand(rng)
                x = Double64(u) + Double64((rand(rng) - 0.5) * eps(u))
                b = lb_big(x)
                for (f, g) in ((sin, sin), (cos, cos))
                    r = g(b)
                    err = Float64(abs(lb_big(f(x)) - r))
                    mabs = max(mabs, err / epsEff)
                    mc = max(mc, err / (epsEff * u))
                    iszero(r) || (mrel = max(mrel, err / Float64(abs(r)) / epsEff))
                end
                xmax = max(xmax, u)
            end
        end
        worst_abs = max(worst_abs, mabs); worst_c = max(worst_c, mc)
        say("PR-trig | Double64 sin, cos, |x| in [1e$e, 1e$(e + 1)), 3000 evaluations: max abs error $(fmt(mabs)) x eps_eff, max relative error $(fmt(mrel)) x eps_eff, max abs error / (eps_eff |x|) = $(fmt(mc)) | eps_eff = 2^-100 | INFO")
    end
    say("PR-trig | summary: Double64 trig floor over |x| in [1, 1e8): max abs error = $(fmt(worst_abs)) eps_eff = $(fmt(worst_abs * epsEff)) (2^$(fmt(log2(worst_abs * epsEff)))); fit err <= c eps_eff |x| with c = $(fmt(worst_c)) | eps_eff = effective_eps(Double64) = $(fmt(epsEff)) | INFO")
    # The same on realized lifted orbit points (K = 2 pi, Double64, 20 starts, n = 1e4, every 7th point): where do the arguments go, and is the error larger there?
    sys = StandardMap(lb_K(Double64, :two_pi))
    xs = Double64[]
    for s in uniform_starts(vs(30, 1), 20)
        pts, _ = lb_pts(sys, lb_cast(Double64, s), 10_000)
        append!(xs, [p[1] for p in pts[1:7:end]])
    end
    bins = Dict{Int, Float64}(); cnt = Dict{Int, Int}(); over = Dict{Int, Int}(); cmax = 0.0
    setprecision(BigFloat, LB_BITS) do
        for x in xs
            e = x == 0 ? 0 : floor(Int, log10(abs(Float64(x))))
            b = lb_big(x)
            err = max(Float64(abs(lb_big(sin(x)) - sin(b))), Float64(abs(lb_big(cos(x)) - cos(b)))) / epsEff
            bins[e] = max(get(bins, e, 0.0), err); cnt[e] = get(cnt, e, 0) + 1
            over[e] = get(over, e, 0) + (err > 1 ? 1 : 0)
            abs(Float64(x)) >= 1e3 && (cmax = max(cmax, err / abs(Float64(x))))
        end
    end
    say("PR-trig | realized orbit arguments, K = 2 pi lifted Double64, $(length(xs)) points, max |x| = $(fmt(Float64(maximum(abs, xs)))): max abs error of sin, cos in units of eps_eff by decade of |x| (count): " *
        join(["1e$e: $(fmt(bins[e])) ($(cnt[e]) points, $(over[e]) above 1)" for e in sort(collect(keys(bins)))], "; ") * "; for |x| >= 1e3 max abs error / (eps_eff |x|) = $(fmt(cmax)) | eps_eff | INFO")
    # Attribution of the Double64 lifted excess over the old linear bound (Amendment D35): error of nX against (i) the 256-bit product at the
    # realized points (E1, the reference of criterion 5) and (ii) the 256-bit product of the tangents as stored in Double64.
    rows = NamedTuple[]
    for (j, s) in enumerate(uniform_starts(vs(30, 2), 6))
        x0 = lb_cast(Double64, s); n = 10_000
        W0 = initial_deviations(vs(30), j, Val(2), Val(1), Double64)
        pts, _ = lb_pts(sys, x0, n)
        f = lb_final(sys, x0, W0, n, 1)
        ref = lb_ref_std(lb_K(Double64, :two_pi), pts, W0[:, 1], n)[end]
        stored = setprecision(BigFloat, LB_BITS) do
            w1, w2 = lb_big(W0[1]), lb_big(W0[2]); n0 = log(w1^2 + w2^2) / 2
            tmax = 0.0
            Kb = lb_big(sys.K)
            for z in pts
                M = tangent(sys, z)
                m11, m12, m21, m22 = lb_big(M[1, 1]), lb_big(M[1, 2]), lb_big(M[2, 1]), lb_big(M[2, 2])
                c = Kb * cos(lb_big(z[1]))
                tmax = max(tmax, Float64(abs(m11 - 1 - c)))
                w1, w2 = m11 * w1 + m12 * w2, m21 * w1 + m22 * w2
            end
            (log(w1^2 + w2^2) / 2 - n0, tmax)
        end
        Ms = [SMatrix{2, 2, Float64}(map(Float64, tangent(sys, z))) for z in pts]
        G = lb_Gfirst(Ms, 1)
        e1 = lb_err(f.sum[1], f.comp[1], ref); e2 = lb_err(f.sum[1], f.comp[1], stored[1])
        old = lb_BS(lb_eps(Double64), n, G)
        push!(rows, (; e1_old = e1 / old, e2_old = e2 / old, tmax = stored[2] / epsEff, x = Float64(maximum(abs(p[1]) for p in pts))))
    end
    say("PR-trig | attribution, K = 2 pi lifted Double64, n = 1e4, 6 starts: max over starts of error / old linear bound 8 eps (n + G_first): vs 256-bit product at realized points $(fmt(maximum(r.e1_old for r in rows))); vs product of the stored Double64 tangents $(fmt(maximum(r.e2_old for r in rows))); max |stored tangent entry - K cos x at 256 bits| = $(fmt(maximum(r.tmax for r in rows))) eps_eff; max |x| on the orbit $(fmt(maximum(r.x for r in rows))) | - | INFO")
end

# ---------------------------------------------------------------------------------------------------
# Figure: convergence with N of the exponent against the ChaosTools oracle (K = 2 pi, K = 0.5) and the large-K answer (K = 10)
section("figure") do
    ns_c = [10, 30, 100, 300, 1000, 3000, 10_000]
    ns_r = [100, 300, 1000, 3000, 10_000, 30_000, 100_000]
    ns_y = [100, 300, 1000, 3000, 10_000, 30_000, 100_000]
    function series(K, x0, ns; seed, W0 = nothing)
        sys = StandardMap(K)
        cl = zeros(length(x0), length(ns)); ct = zeros(length(x0), length(ns))
        Threads.@threads for j in eachindex(x0)
            W = W0 === nothing ? initial_deviations(seed, j, Val(2), Val(1), Float64) : W0[j]
            probe = lb_Probe(Benettin(), ns; alphas = false)
            _, _, st = propagate_orbit(sys, SVector{2, Float64}(x0[j]), W, last(ns), 1, (probe,); orbit = j)
            fin = reducer_final(probe, st[1])
            for (i, (n, S, C)) in enumerate(fin.snaps)
                cl[j, i] = (S + C) / n
            end
            tds = ct_make(K)
            for (i, n) in enumerate(ns)
                ct[j, i] = ct_lambda(tds, x0[j], n)
            end
        end
        return cl, ct
    end
    fig = CM.Figure(size = (1500, 480))
    # panel 1: chaotic K = 2 pi
    xc = uniform_starts(vs(40, 1), 200)
    cl, ct = series(2pi, xc, ns_c; seed = vs(40))
    ax = CM.Axis(fig[1, 1]; xscale = log10, xlabel = "n", ylabel = "ensemble mean of X(n)", title = "K = 2 pi, 200 starts: Chaoslab vs ChaosTools")
    for (m, col, lab, ls) in ((cl, :crimson, "Chaoslab Benettin (seeded w0)", :solid), (ct, :steelblue, "ChaosTools lyapunovspectrum (w0 = e1)", :dash))
        mu = vec(mean(m; dims = 1)); se = vec(std(m; dims = 1)) ./ sqrt(size(m, 1))
        CM.band!(ax, ns_c, mu .- 4se, mu .+ 4se; color = (col, 0.2))
        CM.lines!(ax, ns_c, mu; color = col, linestyle = ls, linewidth = 2, label = lab)
    end
    CM.axislegend(ax; position = :rt, labelsize = 11)
    d = vec(mean(cl; dims = 1)) .- vec(mean(ct; dims = 1))
    se4 = 4 .* sqrt.(vec(var(cl; dims = 1)) ./ 200 .+ vec(var(ct; dims = 1)) ./ 200)
    say("# figure K=2pi: |difference of means| / 4 SE at n = $(ns_c): " * join(fmt.(abs.(d) ./ se4), ", "))
    # panel 2: regular K = 0.5
    ids, starts, _, rseed = lb_regular()
    st = [SVector(s[1], s[2]) for s in starts[1:20]]
    W0 = [initial_deviations(rseed, id, Val(2), Val(1), Float64) for id in ids[1:20]]
    rl, rt = series(0.5, st, ns_r; seed = rseed, W0 = W0)
    ax2 = CM.Axis(fig[1, 2]; xscale = log10, yscale = log10, xlabel = "n", ylabel = "X(n)", title = "K = 0.5 regular class, 20 starts")
    for j in 1:20
        CM.lines!(ax2, ns_r, [x > 0 ? x : NaN for x in rl[j, :]]; color = (:crimson, 0.5))
        CM.lines!(ax2, ns_r, [x > 0 ? x : NaN for x in rt[j, :]]; color = (:steelblue, 0.5), linestyle = :dash)
    end
    CM.lines!(ax2, ns_r, log.(ns_r) ./ ns_r; color = :black, linestyle = :dot, linewidth = 2, label = "ln n / n")
    CM.lines!(ax2, ns_r, vec(maximum(abs.(rl .- rt); dims = 1)); color = :darkorange, linewidth = 2, label = "max |X - X_ChaosTools|")
    CM.hlines!(ax2, [1e-4]; color = :gray40, linestyle = :dashdot, label = "bound 1e-4")
    CM.axislegend(ax2; position = :lb, labelsize = 11)
    say("# figure K=0.5: max |X Chaoslab - X ChaosTools| over 20 starts at n = $(ns_r): " * join(fmt.(vec(maximum(abs.(rl .- rt); dims = 1))), ", "))
    # panel 3: Y1, K = 10
    xy = uniform_starts(vs(40, 3), 200)
    yl, yt = series(10.0, xy, ns_y; seed = vs(40, 3))
    ax3 = CM.Axis(fig[1, 3]; xscale = log10, xlabel = "n", ylabel = "ensemble mean of X(n)", title = "K = 10 (Y1), 200 starts")
    CM.lines!(ax3, ns_y, vec(mean(yl; dims = 1)); color = :crimson, linewidth = 2, label = "Chaoslab Benettin")
    CM.lines!(ax3, ns_y, vec(mean(yt; dims = 1)); color = :steelblue, linestyle = :dash, linewidth = 2, label = "ChaosTools")
    CM.hlines!(ax3, [log(5.0)]; color = :black, linestyle = :dot, label = "ln(K/2)")
    CM.hlines!(ax3, [log(5.0) + 0.01]; color = :gray40, linestyle = :dashdot, label = "ln(K/2) + 1/K^2")
    CM.axislegend(ax3; position = :rb, labelsize = 11)
    say("# figure K=10: ensemble mean X(n) at n = $(ns_y): " * join(fmt.(vec(mean(yl; dims = 1))), ", ") * "; ln(K/2) = $(fmt(log(5.0)))")
    CM.save(joinpath(OUT, "lyap-benettin-convergence.png"), fig)
    say("# figure saved: lyap-benettin-convergence.png")
end

say(@sprintf("# total wall time %.1f s", time() - T0))
close(LOG)
