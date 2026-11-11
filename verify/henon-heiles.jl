# Stage 04 verification of module henon-heiles (catalog row F3: Chaoslab.Systems henon_heiles.jl, Chaoslab.Integrators section.jl and escape.jl; spec D116).
# Run from the repo root, one process per group (HH_ONLY; HH_LOG names the log, HH_APPEND=1 appends):
#   HH_ONLY=FULL    julia --project=verify -t 8 verify/henon-heiles.jl   (the package test files of test/henon-heiles/ at the FULL sizes, HH_FULL_RUN = true; FULL1 to FULL4 = one file each)
#   HH_ONLY=O1      (BigFloat 256-bit closed forms: criteria 1, 4, 14)
#   HH_ONLY=O2      (ForwardDiff and complex-step derivatives, the tangent step against the AD derivative of the one-step map: criterion 2, tangent of 14/17)
#   HH_ONLY=O3      (TaylorIntegration BigFloat references; O3a crossings 7, O3b section Jacobian 8, O3c escape times 12b, O3d spectrum at T = 100 13d)
#   HH_ONLY=O4      (ChaosTools Lyapunov spectrum at tolerance 1e-14: criterion 13d)
#   HH_ONLY=O5      (portrait at E = 1/8: criterion 10; saves the portrait for FIG)
#   HH_ONLY=O6      (escape continuation to 2 R_esc with Taylor, tau/8 self-convergence: criterion 12a)
#   HH_ONLY=AC16T   (1 against N CPU threads, bitwise; spawns a one-thread child)
#   HH_ONLY=FIG     (henon-heiles-energy.png, henon-heiles-portrait.png; needs O5 first)
# Not run: the GPU part of criterion 16 (verify/gpu.jl writes data/; owner decision). Nothing here writes under data/.
# Spec: workspaces/module-build/stages/01-spec/output/henon-heiles-spec.md; test plan points U1 to U10 are the operative reading. Seeds 20261012 + 42000 .. 42999 (hh_seed); verify-only 42900 ...
# One printed line per measurement:  AC<n> | quantity | value | bound | PASS/FAIL/INFO
const ONLY = filter(!isempty, split(get(ENV, "HH_ONLY", ""), ","))
using Printf, Random, LinearAlgebra, Statistics, Serialization, Logging
using Test
const HH_FULL_RUN = true
using Chaoslab, StaticArrays, DoubleFloats, StableRNGs, ForwardDiff, TOML
const NEEDTAY = any(o -> startswith(o, "O3") || startswith(o, "O6"), ONLY)
NEEDTAY && (import TaylorIntegration as TI)
"O4" in ONLY && (import ChaosTools as CT)
"FIG" in ONLY && (import CairoMakie as CM)

const ROOT = dirname(@__DIR__)
const OUT = joinpath(ROOT, "workspaces", "module-build", "stages", "04-verify", "output")
const SCR = get(ENV, "HH_SCRATCH", joinpath(tempdir(), "hh-verify-scratch"))
mkpath(SCR)
const LOG = open(get(ENV, "HH_LOG", joinpath(OUT, "henon-heiles-verify-log.txt")), get(ENV, "HH_APPEND", "") == "1" ? "a" : "w")

fmt(x::AbstractFloat) = @sprintf("%.4g", x)
fmt(x) = string(x)
const LOGLOCK = ReentrantLock()
function say(s)
    lock(LOGLOCK) do
        println(s); println(LOG, s); flush(stdout); flush(LOG)
    end
end
emit(ac, quantity, value, bound, verdict) = say(string("AC", ac, " | ", quantity, " | ", fmt(value), " | ", bound, " | ", verdict))
pf(ok::Bool) = ok ? "PASS" : "FAIL"

struct TeeLog <: Logging.AbstractLogger end
Logging.min_enabled_level(::TeeLog) = Logging.Info
Logging.shouldlog(::TeeLog, args...) = true
Logging.catch_exceptions(::TeeLog) = true
function Logging.handle_message(::TeeLog, level, message, _mod, group, id, file, line; kwargs...)
    kv = join(("$k = $(first(repr(v), 900))" for (k, v) in kwargs), "; ")
    say(string("# ", level, " ", message, " | ", kv))
end
global_logger(TeeLog())

want(name) = isempty(ONLY) || any(o -> startswith(name, o), ONLY)
function section(f, name)
    want(name) || return
    t0 = time()
    try
        f()
    catch e
        say(string(name, " | ERROR | ", sprint(showerror, e, catch_backtrace())[1:min(end, 1800)], " | - | FAIL"))
    end
    say(@sprintf("# section %s took %.1f s", name, time() - t0))
end

say("# henon-heiles verification, julia $(VERSION), threads $(Threads.nthreads()), base seed 20261012 + 42000, HH_FULL_RUN = true, groups $(isempty(ONLY) ? "all" : join(ONLY, ",")), time $(round(time())) (unix)")
say("# timings in this log were taken on a machine that may be loaded by other processes")

const TD = joinpath(ROOT, "test", "henon-heiles")
include(joinpath(TD, "hh_helpers.jl"))
const VSEED = 20261012 + 42900
const BITS = 256

# ---------------------------------------------------------------------------------------------------
# FULL: the package test files at the full sizes
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
        nf = nfail(ch); c = Test.get_test_counts(ch)
        emit(ac, string("testset ", ch.description, " (", c.cumulative_passes, " pass, ", nf, " fail or error, ", c.cumulative_broken, " broken, ", @sprintf("%.1f", ch.time_end - ch.time_start), " s)"), nf, "0", pf(nf == 0))
    end
end
function suite_section(name, file)
    section(name) do
        Test.@testset "outer $name" begin
            ts = Test.@testset "$name" begin
                include(joinpath(TD, file))
            end
            report_ts(name, ts)
        end
    end
end
suite_section("FULL1", "hh_ac01_05.jl")      # criteria 1 to 5
suite_section("FULL2", "hh_ac06_09.jl")      # criteria 6 to 9
suite_section("FULL3", "hh_ac10_12.jl")      # criteria 10 to 12
suite_section("FULL4", "hh_ac13_17.jl")      # criteria 13 to 17

# ---------------------------------------------------------------------------------------------------
# independent test-side physics in BigFloat (call inside setprecision)
bV(q1, q2) = (q1^2 + q2^2) / 2 + q1^2 * q2 - q2^3 / 3
bg(q1, q2) = (q1 + 2q1 * q2, q2 + q1^2 - q2^2)
bu(q2) = q2^2 / 2 - q2^3 / 3
function bbisect(f, lo, hi)                      # sign change on [lo, hi]
    flo = f(lo) > 0
    for _ in 1:(BITS + 8)
        mid = (lo + hi) / 2
        (f(mid) > 0) == flo ? (lo = mid) : (hi = mid)
    end
    return (lo + hi) / 2
end
function binterval(E)                            # Equation 7 roots a < 0 < b of u = E
    E == BigFloat(1) / 6 && return (-BigFloat(1) / 2, BigFloat(1))
    return (bbisect(q -> bu(q) - E, -BigFloat(1), BigFloat(0)), bbisect(q -> bu(q) - E, BigFloat(0), BigFloat(1)))
end
function barea(E; n = 2000)                      # A(E) = int_a^b 2 sqrt(2 (E - u)) dq2, q2 = (a+b)/2 - (b-a)/2 cos(phi), midpoint rule in phi
    a, b = binterval(E); s = BigFloat(0)
    for k in 1:n
        phi = BigFloat(pi) * (k - BigFloat(1) / 2) / n
        q2 = (a + b) / 2 - (b - a) / 2 * cos(phi)
        s += 2 * sqrt(max(2 * (E - bu(q2)), BigFloat(0))) * (b - a) / 2 * sin(phi)
    end
    return s * BigFloat(pi) / n
end

# ------------------------------------------------------------------------------------------------ O1
section("O1") do
    flow = HenonHeilesFlow(); tol = BigFloat(2)^-200
    setprecision(BigFloat, BITS) do
        s3 = sqrt(BigFloat(3)); h = BigFloat(1) / 2
        pts = ((BigFloat(0), BigFloat(0)), (BigFloat(0), BigFloat(1)), (s3 / 2, -h), (-s3 / 2, -h))
        vals = (BigFloat(0), BigFloat(1) / 6, BigFloat(1) / 6, BigFloat(1) / 6)
        eigs = ((BigFloat(1), BigFloat(1)), (BigFloat(3), -BigFloat(1)), (BigFloat(3), -BigFloat(1)), (BigFloat(3), -BigFloat(1)))
        gmax = BigFloat(0); vmax = BigFloat(0); emax = BigFloat(0); cpmax = BigFloat(0); fgmax = BigFloat(0)
        cpb = critical_points(flow, BigFloat)
        for k in 1:4
            q1, q2 = pts[k]; g = bg(q1, q2)
            gmax = max(gmax, abs(g[1]), abs(g[2])); vmax = max(vmax, abs(bV(q1, q2) - vals[k]))
            tr = (1 + 2q2) + (1 - 2q2); dt = (1 + 2q2) * (1 - 2q2) - (2q1)^2; disc = sqrt(tr^2 - 4dt)
            emax = max(emax, abs((tr + disc) / 2 - eigs[k][1]), abs((tr - disc) / 2 - eigs[k][2]))
            cpmax = max(cpmax, abs(cpb[k][1] - q1), abs(cpb[k][2] - q2))
            fg = force(flow, SVector{2, BigFloat}(q1, q2)); fgmax = max(fgmax, abs(fg[1] - g[1]), abs(fg[2] - g[2]))
        end
        emit("1", "O1: max |grad V| at O, S0, S1, S2 (own formula, 256 bits)", Float64(gmax), "2^-200 = 6.2e-61", pf(gmax <= tol))
        emit("1", "O1: max |V - (0, 1/6, 1/6, 1/6)| (256 bits)", Float64(vmax), "2^-200", pf(vmax <= tol))
        emit("1", "O1: max |Hessian eigenvalue - (1,1), (3,-1) x 3| (256 bits)", Float64(emax), "2^-200", pf(emax <= tol))
        emit("1", "O1: module critical_points(BigFloat) against the closed forms; module force against own gradient at them", string(fmt(Float64(cpmax)), " / ", fmt(Float64(fgmax))), "2^-200", pf(cpmax <= tol && fgmax <= tol))
        gn = [sqrt(sum(abs2, bg(sg * s3 / 2, h))) for sg in (1, -1)]
        emit("1", "O1: gradient norm at the benchmarks.md points (+-sqrt(3)/2, +1/2) (not critical)", Float64(minimum(gn)), ">= 1", pf(minimum(gn) >= 1))
        # exact rationals: q1 enters V and g only through q1^2 = 3/4 at S1,2 (g1 = q1 (1 + 2 q2) = 0 at q2 = -1/2)
        q1sq = 3 // 4; q2 = -1 // 2
        Vex = (q1sq + q2^2) / 2 + q1sq * q2 - q2^3 / 3; g2ex = q2 + q1sq - q2^2
        Vs0 = (0 + 1 // 1) / 2 - (1 // 1) / 3; g2s0 = (1 // 1) + 0 - (1 // 1)
        emit("1", "O1: exact rationals: V(S1,2) = 1/6, g2(S1,2) = 0, V(S0) = 1/6, g2(S0) = 0", string(Vex, ", ", g2ex, ", ", Vs0, ", ", g2s0), "1//6, 0, 1//6, 0", pf(Vex == 1 // 6 && g2ex == 0 && Vs0 == 1 // 6 && g2s0 == 0))
        for T in (Float64, Double64)
            c = critical_points(flow, T); w = 0.0
            for k in 1:4, i in 1:2
                w = max(w, Float64(abs(hh_big(c[k][i]) - pts[k][i])))
            end
            emit("1", "O1: critical_points($(T)) rounding error against 256-bit points", w, "<= eps($(T)) = $(fmt(hh_eps(T)))", pf(w <= hh_eps(T)))
        end
    end
    # criterion 4: harmonic limit of A(E), quadrature convergence, module grid area fraction against the 256-bit area
    setprecision(BigFloat, BITS) do
        defects = BigFloat[]
        Es = [BigFloat(10)^-k for k in 2:6]
        for E in Es
            A2 = barea(E; n = 2000); A4 = barea(E; n = 4000)
            push!(defects, A4 / (2 * BigFloat(pi) * E) - 1)
            emit("4", @sprintf("O1: quadrature convergence of A(E) at E = %.0e (n = 2000 against 4000), relative", Float64(E)), Float64(abs(A2 - A4) / A4), "<= 2^-100 = 7.9e-31", pf(abs(A2 - A4) / A4 <= BigFloat(2)^-100))
        end
        for (E, d) in zip(Es, defects)
            emit("4", @sprintf("O1: A(E) / (2 pi E) - 1 at E = %.0e, and (defect)/E", Float64(E)), string(fmt(Float64(d)), " ; ", fmt(Float64(d / E))), "defect -> 0 (harmonic limit)", "INFO")
        end
        c = Float64.(defects ./ Es)
        emit("4", "O1: harmonic limit: |A/(2 pi E) - 1| decreases by a factor >= 9 per decade of E, and (defect)/E converges (last two within 1e-3 relative)",
             string(fmt(Float64(abs(defects[end]))), " ; ratio ", fmt(Float64(abs(defects[end - 1] / defects[end]))), " ; (defect/E) ", fmt(c[end - 1]), " -> ", fmt(c[end])),
             "ratio >= 9; converged 1e-3", pf(all(i -> abs(defects[i]) / abs(defects[i + 1]) >= 9, 1:(length(defects) - 1)) && abs(c[end] - c[end - 1]) <= 1.0e-3 * abs(c[end])))
        for E in (1 / 12, 1 / 8, 1 / 6)
            zs = energy_grid(flow, E, 1000); a, b = section_interval(flow, E); pm = sqrt(2E)
            Eb = BigFloat(1) / 12 * (E == 1 / 12) + BigFloat(1) / 8 * (E == 1 / 8) + BigFloat(1) / 6 * (E == 1 / 6)
            Ab = barea(Eb; n = 3000)
            frac = BigFloat(length(zs)) * BigFloat(b - a) * 2 * BigFloat(pm) / 1000^2 / Ab
            emit("4", @sprintf("O1: energy_grid 1000 x 1000 at E = %.6g: n points = %d, area fraction against the 256-bit A(E)", E, length(zs)), Float64(abs(frac - 1)), "<= 1e-3", pf(abs(frac - 1) <= 1.0e-3))
        end
    end
    # criterion 14: tangent at S0 and O against closed forms; the matrix exponential by its Taylor series in 256 bits is the independent check of the closed form
    setprecision(BigFloat, BITS) do
        t = BigFloat(5); s3 = sqrt(BigFloat(3))
        A = zeros(BigFloat, 4, 4); A[1, 3] = 1; A[2, 4] = 1; A[3, 1] = -3; A[4, 2] = 1     # tangent flow at S0: Hess = diag(3, -1)
        E = Matrix{BigFloat}(I, 4, 4); term = Matrix{BigFloat}(I, 4, 4)
        for n in 1:200
            term = term * (t * A) / n; E += term
        end
        Mcl = zeros(BigFloat, 4, 4)
        Mcl[1, 1] = cos(s3 * t); Mcl[1, 3] = sin(s3 * t) / s3; Mcl[3, 1] = -s3 * sin(s3 * t); Mcl[3, 3] = cos(s3 * t)
        Mcl[2, 2] = cosh(t); Mcl[2, 4] = sinh(t); Mcl[4, 2] = sinh(t); Mcl[4, 4] = cosh(t)
        d1 = maximum(abs.(E .- Mcl))
        emit("14", "O1: closed-form tangent at S0 (T = 5) against the 256-bit exponential series of the constant matrix", Float64(d1), "<= 1e-60", pf(d1 <= BigFloat(10)^-60))
        flowc = HenonHeilesFlow(); I4 = one(SMatrix{4, 4, Float64})
        cp = critical_points(flowc, Float64); zS = SVector(cp[2][1], cp[2][2], 0.0, 0.0); zO = SVector(cp[1][1], cp[1][2], 0.0, 0.0)
        m = KahanLi(Float64, 8); tau = 0.05; z = zS; W = I4
        for _ in 1:100
            z, W = tangent_advance(m, flowc, z, W, tau)
        end
        e14 = maximum(abs.(Float64.(Mcl) .- Matrix(W)))
        emit("14", "O1: S0, T = 5, order-8 table, Float64, tau = 1/20: max |tangent - closed form|", e14, "<= 1e-9", pf(e14 <= 1.0e-9))
        mb = KahanLi(BigFloat, 8); fb = HenonHeilesFlow(); zb = SVector{4, BigFloat}(0, 1, 0, 0); Wb = one(SMatrix{4, 4, BigFloat})
        for _ in 1:100
            zb, Wb = tangent_advance(mb, fb, zb, Wb, BigFloat(1) / 20)
        end
        eb = maximum(abs.(Mcl .- Matrix(Wb)))
        emit("14", "O1: S0, T = 5, order-8 table in BigFloat, tau = 1/20: max |tangent - closed form| (truncation error times e^T, informational)", Float64(eb), "<= 1e-9", pf(Float64(eb) <= 1.0e-9))
        # O: rotation with unit frequencies
        n = round(Int, 10 / tau); z = zO; W = I4
        for _ in 1:n
            z, W = tangent_advance(m, flowc, z, W, tau)
        end
        c10 = cos(BigFloat(10)); s10 = sin(BigFloat(10)); R = zeros(BigFloat, 4, 4)
        R[1, 1] = c10; R[2, 2] = c10; R[3, 3] = c10; R[4, 4] = c10; R[1, 3] = s10; R[2, 4] = s10; R[3, 1] = -s10; R[4, 2] = -s10
        eO = maximum(abs.(Float64.(R) .- Matrix(W)))
        emit("14", "O1: O, T = 10: max |tangent - 256-bit unit-frequency rotation|", eO, "<= B(n) = $(fmt(hh_B(Float64, 1.0, n, tau)))", pf(eO <= hh_B(Float64, 1.0, n, tau)))
        lam = sort(collect(hh_lyap(m, flowc, zS, tau, 100.0).lambda); rev = true)
        emit("14", "O1: S0 spectrum at T = 100 (sorted)", join(fmt.(lam), ", "), "(1, 0, 0, -1) within 0.05", pf(maximum(abs.(lam .- (1.0, 0.0, 0.0, -1.0))) <= 0.05))
        lo = maximum(abs, hh_lyap(m, flowc, zO, tau, 1000.0).lambda)
        emit("14", "O1: O spectrum at T = 1000, max |lambda|", lo, "<= 1e-12", pf(lo <= 1.0e-12))
    end
end

# ------------------------------------------------------------------------------------------------ O2
section("O2") do
    flow = HenonHeilesFlow()
    for (T, seed) in ((Float64, 1), (Double64, 2))
        rng = StableRNG(VSEED + seed); et = hh_eps(T); wg = 0.0; wh = 0.0; wcs = 0.0
        for _ in 1:1000
            x = 3 .* rand(rng, 2) .- 1.5; q = SVector{2, T}(T(x[1]), T(x[2]))
            gad = ForwardDiff.gradient(w -> hh_V(w[1], w[2]), T[q[1], q[2]]); had = ForwardDiff.hessian(w -> hh_V(w[1], w[2]), T[q[1], q[2]])
            fad = ForwardDiff.jacobian(w -> collect(force(flow, SVector{2}(w[1], w[2]))), T[q[1], q[2]])
            g = force(flow, q); Hs = force_jacobian(flow, q)
            for i in 1:2
                wg = max(wg, Float64(abs(g[i] - gad[i])) / (et * max(1, abs(Float64(gad[i])))))
                for j in 1:2
                    wh = max(wh, Float64(abs(Hs[i, j] - had[i, j])) / (et * max(1, abs(Float64(had[i, j])))), Float64(abs(Hs[i, j] - fad[i, j])) / (et * max(1, abs(Float64(fad[i, j])))))
                end
            end
            if T === Float64
                h = 1.0e-20
                for j in 1:2
                    e = SVector(j == 1 ? 1.0 : 0.0, j == 2 ? 1.0 : 0.0)
                    gcs = imag(hh_V(q[1] + im * h * e[1], q[2] + im * h * e[2])) / h
                    fcs = imag.(force(flow, SVector{2, ComplexF64}(q[1] + im * h * e[1], q[2] + im * h * e[2]))) ./ h
                    wcs = max(wcs, abs(g[j] - gcs) / (et * max(1, abs(gcs))), maximum(abs(Hs[i, j] - fcs[i]) / (et * max(1, abs(fcs[i]))) for i in 1:2))
                end
            end
        end
        emit("2", "O2: $(T), 1000 seeded points in [-1.5, 1.5]^2: worst |force - ForwardDiff grad V| / (eps max(1, |entry|))", wg, "<= 8", pf(wg <= 8))
        emit("2", "O2: $(T): worst |force_jacobian - ForwardDiff (hessian of V and jacobian of force)| / (eps max(1, |entry|))", wh, "<= 8", pf(wh <= 8))
        T === Float64 && emit("2", "O2: Float64: worst complex-step (V gradient, force Hessian columns) / (eps max(1, |entry|)), h = 1e-20", wcs, "<= 8", pf(wcs <= 8))
    end
    # the tangent step against the AD derivative of the one-step map (integrator-catalog O3 form: 16 eps (1 + |M|) k per step)
    rng = StableRNG(VSEED + 3)
    for (nm, m) in (("order-8 KahanLi", KahanLi(Float64, 8)), ("order-4 ForestRuth", ForestRuth(Float64)), ("leapfrog KDK", Leapfrog(Float64, :KDK)))
        k = max(1, stage_count(m)); tau = 0.05
        for nsteps in (1, 20)
            worst = 0.0
            for _ in 1:200
                zs, _ = random_energy_start(flow, 1 / 8, (-0.4, 0.4), (-0.5, 0.5), rng)
                zs = SVector{4, Float64}(zs[1] + 0.1 * (rand(rng) - 0.5), zs[2], zs[3], zs[4])
                W = one(SMatrix{4, 4, Float64}); z = zs
                for _ in 1:nsteps
                    z, W = tangent_advance(m, flow, z, W, tau)
                end
                B = ForwardDiff.jacobian(w -> advance_n(m, flow, w, tau, nsteps), zs)
                worst = max(worst, maximum(abs, Matrix(W) .- B) / (16 * hh_eps(Float64) * (1 + maximum(abs, W)) * k * nsteps))
            end
            emit("2", "O2: tangent step ($(nsteps) steps) against ForwardDiff of the one-step map, $(nm), worst error / (16 eps (1 + |M|) k n)", worst, "<= 1", pf(worst <= 1))
        end
    end
end

# ------------------------------------------------------------------------------------------------ Taylor helpers (O3, O6): hand-written vector field, independent of Chaoslab's
function hh_f!(dx, x, p, s)
    dx[1] = x[3]; dx[2] = x[4]; dx[3] = -(x[1] + 2 * x[1] * x[2]); dx[4] = -(x[2] + x[1]^2 - x[2]^2)
    return nothing
end
function var_f!(dx, x, p, s)                    # state (4) and the 4 x 4 frame, column j at 4j+1 .. 4j+4: W' = [[0, I], [-Hess V, 0]] W
    q1 = x[1]; q2 = x[2]
    dx[1] = x[3]; dx[2] = x[4]; dx[3] = -(q1 + 2 * q1 * q2); dx[4] = -(q2 + q1^2 - q2^2)
    h11 = 1 + 2 * q2; h12 = 2 * q1; h22 = 1 - 2 * q2
    for j in 1:4
        o = 4 * j
        dx[o + 1] = x[o + 3]; dx[o + 2] = x[o + 4]
        dx[o + 3] = -(h11 * x[o + 1] + h12 * x[o + 2]); dx[o + 4] = -(h12 * x[o + 1] + h22 * x[o + 2])
    end
    return nothing
end
tay(f!, x0, T; order = 40, maxsteps = 10^6) = TI.taylorinteg(f!, x0, zero(BigFloat), BigFloat(T), order, ldexp(BigFloat(1), -200); maxsteps = maxsteps)
# root of the dense polynomial g(s) = 0 on [0, dt] with g(0) < 0 <= g(dt), bisection to 2^-250
function poly_root(g, dt)
    lo = zero(dt); hi = dt
    for _ in 1:250
        mid = (lo + hi) / 2
        g(mid) < 0 ? (lo = mid) : (hi = mid)
    end
    return (lo + hi) / 2
end
# upward crossings (q1: - to +, p1 > 0) of the dense Taylor solution: (t, q2, p2) as BigFloat vectors; a start on the section is the first point at t = 0
function tay_section(sol; start = true)
    ts = BigFloat[]; q2s = BigFloat[]; p2s = BigFloat[]; x = sol.x; t = sol.t
    if start && x[1, 1] == 0 && x[1, 3] > 0
        push!(ts, BigFloat(0)); push!(q2s, x[1, 2]); push!(p2s, x[1, 4])
    end
    for i in 1:(length(t) - 1)
        if x[i, 1] < 0 && x[i + 1, 1] >= 0 && x[i + 1, 3] > 0
            P = sol.p[i, 1]; s = poly_root(w -> P(w), t[i + 1] - t[i])
            push!(ts, t[i] + s); push!(q2s, sol.p[i, 2](s)); push!(p2s, sol.p[i, 4](s))
        end
    end
    return (t = ts, q2 = q2s, p2 = p2s)
end
# the next upward crossing from the section point (E, q2, p2) in 256 bits: (q2', p2')
function tay_return(E, q2, p2)
    p1 = sqrt(2 * (E - p2^2 / 2 - bu(q2))); x = [BigFloat(0), q2, p1, p2]
    for _ in 1:40
        sol = tay(hh_f!, x, 8)
        c = tay_section(sol; start = false)
        isempty(c.t) || return (c.q2[1], c.p2[1])
        x = copy(sol.x[end, :])
    end
    return (BigFloat(NaN), BigFloat(NaN))
end
# first time r > R (BigFloat), horizon tmax, from the Float64 start z0: time or nothing
function tay_escape(z0, R, tmax)
    sol = tay(hh_f!, BigFloat.(collect(z0)), tmax; maxsteps = 20000); x = sol.x; t = sol.t; R2 = BigFloat(R)^2
    for i in 1:(length(t) - 1)
        if x[i, 1]^2 + x[i, 2]^2 <= R2 < x[i + 1, 1]^2 + x[i + 1, 2]^2
            P1 = sol.p[i, 1]; P2 = sol.p[i, 2]
            return t[i] + poly_root(w -> P1(w)^2 + P2(w)^2 - R2 + 0 * w, t[i + 1] - t[i])
        end
    end
    return nothing
end

# ------------------------------------------------------------------------------------------------ O3a: crossings (criterion 7)
section("O3a") do
    flow = HenonHeilesFlow(); m = KahanLi(Float64, 8); tau = 0.05; TREF = 200.0; NST = 20; TAG = 60.0
    say("# O3a: reference Taylor order 40, abstol 2^-200, 256 bits; $(NST) seeded starts per energy (the starts of AC7a of the package suite); reference horizon $(TREF) (spec: 1000: a chaotic start diverges from any other integrator by exp(lambda t), spec point U1), agreement window t <= $(TAG)")
    setprecision(BigFloat, BITS) do
        for (k, E) in enumerate((1 / 12, 1 / 8, 1 / 6))
            rng = StableRNG(hh_seed(7, k)); a, b = section_interval(flow, E); pm = sqrt(2E)
            nbad = 0; wt = 0.0; wq = 0.0; wp = 0.0; hor = Float64[]; ncr = Int[]
            t0 = time()
            for _ in 1:NST
                z0, _ = random_energy_start(flow, E, (a, b), (-pm, pm), rng)
                rec = section_map(m, flow, z0, tau, round(Int, TREF / tau))
                ref = tay_section(tay(hh_f!, BigFloat.(collect(z0)), TREF))
                ir = findall(<=(TAG), Float64.(ref.t)); ic = findall(<=(TAG), rec.t)
                push!(ncr, length(ir))
                if length(ir) != length(ic)
                    nbad += 1
                else
                    wt = max(wt, maximum(abs.(rec.t[ic] .- Float64.(ref.t[ir])))); wq = max(wq, maximum(abs.([rec.z[i][2] for i in ic] .- Float64.(ref.q2[ir]))))
                    wp = max(wp, maximum(abs.([rec.z[i][4] for i in ic] .- Float64.(ref.p2[ir]))))
                end
                # agreement horizon: time of the first crossing where count or location differ by more than 1e-8
                h = TREF
                for i in 1:min(length(rec.t), length(ref.t))
                    d = max(abs(rec.t[i] - ref.t[i]), abs(rec.z[i][2] - ref.q2[i]), abs(rec.z[i][4] - ref.p2[i]))
                    d > 1.0e-8 && (h = Float64(ref.t[i]); break)
                end
                length(rec.t) != length(ref.t) && (h = min(h, Float64(ref.t[min(length(rec.t), length(ref.t))])))
                push!(hor, h)
            end
            Es = @sprintf("%.6g", E)
            emit("7", "O3a: E = $(Es), $(NST) starts, crossings t <= $(TAG): starts whose count differs from the Taylor reference (counts per start $(minimum(ncr)) to $(maximum(ncr)))", nbad, "0", pf(nbad == 0))
            emit("7", "O3a: E = $(Es): max |t_cross - reference| over the window", wt, "<= 1e-8", pf(wt <= 1.0e-8))
            emit("7", "O3a: E = $(Es): max |q2 - reference|, max |p2 - reference|", string(fmt(wq), " ; ", fmt(wp)), "<= 1e-8", pf(wq <= 1.0e-8 && wp <= 1.0e-8))
            emit("7", "O3a: E = $(Es): agreement horizon (first crossing differing by more than 1e-8; reference to $(TREF)): min / median over starts", string(fmt(minimum(hor)), " ; ", fmt(median(hor))), "info (window bound 60)", "INFO")
            say(@sprintf("# O3a E = %s took %.1f s", Es, time() - t0))
        end
    end
end

# ------------------------------------------------------------------------------------------------ O3b: section Jacobian (criterion 8)
section("O3b") do
    flow = HenonHeilesFlow(); m = KahanLi(Float64, 8); tau = 0.05; et = hh_eps(Float64); dtol = hh_defaults()["det_residual_tol"]; im = IntegratorMap(m, flow, tau)
    NST = 10
    say("# O3b: reference Jacobian = central difference (step 2^-100) of the Taylor return map (order 40, abstol 2^-200, 256 bits), legs 2 to 4 of $(NST) starts per energy")
    setprecision(BigFloat, BITS) do
        h = ldexp(BigFloat(1), -100)
        for (k, E) in enumerate((1 / 12, 1 / 8, 1 / 6))
            rng = StableRNG(hh_seed(8, k)); a, b = section_interval(flow, E); pm = sqrt(2E)
            wj = 0.0; wd = 0.0; nl = 0; ok_det = true
            for _ in 1:NST
                z0, _ = random_energy_start(flow, E, (a, b), (-pm, pm), rng)
                zb = SVector{4, BigFloat}(BigFloat.(z0)); Eb = (zb[3]^2 + zb[4]^2) / 2 + bV(zb[1], zb[2])
                rec = section_map(m, flow, z0, tau, round(Int, 30 / tau); jacobian = true)
                for i in 2:min(length(rec.t), 4)
                    J = rec.J[i]; zp = rec.z[i - 1]; q2 = BigFloat(zp[2]); p2 = BigFloat(zp[4])
                    cq = (tay_return(Eb, q2 + h, p2), tay_return(Eb, q2 - h, p2)); cp = (tay_return(Eb, q2, p2 + h), tay_return(Eb, q2, p2 - h))
                    Jref = [Float64((cq[1][1] - cq[2][1]) / (2h)) Float64((cp[1][1] - cp[2][1]) / (2h)); Float64((cq[1][2] - cq[2][2]) / (2h)) Float64((cp[1][2] - cp[2][2]) / (2h))]
                    wj = max(wj, maximum(abs.(Matrix(J) .- Jref) ./ max.(1, abs.(Jref))))
                    nsteps = round(Int, (rec.t[i] - rec.t[i - 1]) / tau) + 1
                    M = max(1.0, maximum(abs, tangent(im, zp)), maximum(abs, J))      # Amendment D117
                    bound = dtol * et * (1 + M) * nsteps + 4 * et
                    wd = max(wd, abs(det(J) - 1) / bound); abs(det(J) - 1) <= bound || (ok_det = false)
                    nl += 1
                end
            end
            Es = @sprintf("%.6g", E)
            emit("8", "O3b: E = $(Es), $(nl) legs: max relative |J - Taylor finite-difference Jacobian| / max(1, |entry|)", wj, "<= 1e-6", pf(wj <= 1.0e-6))
            emit("8", "O3b: E = $(Es): worst |det J - 1| / (det_residual_tol eps (1 + M) n_steps + 4 eps)", wd, "<= 1", pf(ok_det))
        end
    end
end

# ------------------------------------------------------------------------------------------------ O3c: escape times (criterion 12b)
section("O3c") do
    flow = HenonHeilesFlow(); m = KahanLi(Float64, 8); tau = 0.05; R = 10.0; n = round(Int, 200 / tau); NREF = 100
    say("# O3c: reference first passage r = R_esc = 10 by the dense Taylor polynomials (order 40, abstol 2^-200, 256 bits), horizon = module time + 0.05; $(NREF) escapes before t = 100 per energy (starts as AC12b of the package suite)")
    setprecision(BigFloat, BITS) do
        for (k, E) in enumerate((0.2, 0.25))
            rng = StableRNG(hh_seed(12, 10 + k)); done = 0; worst = 0.0; nmis = 0; tmax = 0.0; ntry = 0
            rows = Tuple{Float64, Float64, Float64}[]                       # (t_escape, |t - ref|, |t(tau/8) - ref|)
            while done < NREF
                z0, _ = random_energy_start(flow, E, (-1.0, 1.0), (-1.0, 1.0), rng); ntry += 1
                res = escape_orbit(m, flow, z0, tau, n; R_esc = R)
                (res.status === :escaped && res.t < 100) || continue
                ref = tay_escape(z0, R, res.t + 0.05)
                if ref === nothing
                    nmis += 1
                else
                    r8 = escape_orbit(m, flow, z0, tau / 8, 8n; R_esc = R)
                    push!(rows, (res.t, abs(res.t - Float64(ref)), abs(r8.t - Float64(ref)))); res.t < 50 && (worst = max(worst, rows[end][2]))   # Amendment D117: scored before t = 50
                end
                tmax = max(tmax, res.t); done += 1
            end
            emit("12", @sprintf("O3c: E = %.2f, %d escapes before t = 100 (of %d starts, latest escape t = %.1f): reference found / max |t_escape - reference| for t < 50 (D117; t >= 50 reported below)", E, NREF, ntry, tmax), string(NREF - nmis, "/", NREF, " ; ", fmt(worst)), "all found; <= 1e-6", pf(nmis == 0 && worst <= 1.0e-6))
            nover = count(r -> r[2] > 1.0e-6, rows)
            emit("12", @sprintf("O3c: E = %.2f: escapes with |t - reference| > 1e-6 at tau = 1/20 / the same count at tau/8 / max difference at tau/8", E), string(nover, " ; ", count(r -> r[3] > 1.0e-6, rows), " ; ", fmt(maximum(r[3] for r in rows))), "info", "INFO")
            for (lo, hi) in ((0.0, 25.0), (25.0, 50.0), (50.0, 100.0))
                sub = filter(r -> lo <= r[1] < hi, rows)
                isempty(sub) || emit("12", @sprintf("O3c: E = %.2f: escapes with t in [%g, %g): count, max |t - reference| at tau, at tau/8", E, lo, hi), string(length(sub), " ; ", fmt(maximum(r[2] for r in sub)), " ; ", fmt(maximum(r[3] for r in sub))), "info", "INFO")
            end
        end
    end
end

# ------------------------------------------------------------------------------------------------ O3d: finite-time spectrum at T = 100 against the Taylor variational reference (criterion 13d)
function tay_lyap(z0, W0, T)
    x = vcat(BigFloat.(collect(z0)), BigFloat.(vec(Matrix(W0)))); lam = zeros(BigFloat, 4)
    for _ in 1:T
        sol = tay(var_f!, x, 1); x = copy(sol.x[end, :])
        W = reshape(x[5:20], 4, 4); Q = similar(W)
        for j in 1:4
            v = W[:, j]
            for i in 1:(j - 1)
                v = v - dot(Q[:, i], v) * Q[:, i]
            end
            nv = sqrt(sum(abs2, v)); lam[j] += log(nv); Q[:, j] = v / nv
        end
        x[5:20] = vec(Q)
    end
    return lam ./ T
end
function hh_chaotic_start(flow, m, E, rng; tau = 0.05)
    a, b = section_interval(flow, E); pm = sqrt(2E)
    for _ in 1:200
        z, st = random_energy_start(flow, E, (a, b), (-pm, pm), rng)
        st === :ok && hh_lambda1(m, flow, z, tau, 500.0) > 0.02 && return z
    end
    return nothing
end
section("O3d") do
    flow = HenonHeilesFlow(); m = KahanLi(Float64, 8); tau = 0.05
    zC, _ = energy_start(flow, 1 / 8, -0.015, 0.0)
    z6 = hh_chaotic_start(flow, m, 1 / 6, StableRNG(hh_seed(13, 1)))
    setprecision(BigFloat, BITS) do
        for (name, z0, sd) in (("C (E = 1/8)", zC, 2), ("chaotic start at E = 1/6", z6, 5))
            W0 = initial_deviations(hh_seed(13, sd), 1, Val(4), Val(4), Float64)
            r = hh_lyap(m, flow, z0, tau, 100.0; W0 = W0)
            t1 = @elapsed lref = tay_lyap(z0, W0, 100)
            d = maximum(abs.(r.lambda .- Float64.(lref)))
            say(@sprintf("# O3d %s: module (%s), Taylor (%s), %.1f s", name, join(fmt.(r.lambda), ", "), join(fmt.(Float64.(lref)), ", "), t1))
            emit("13", "O3d: $(name), T = 100, seeded Haar frame, renormalization every time unit, Taylor variational equations (MGS in 256 bits): max |lambda_i - reference|", d, "<= 1e-6", pf(d <= 1.0e-6))
        end
    end
end

# ------------------------------------------------------------------------------------------------ O4: ChaosTools (criterion 13d)
section("O4") do
    flow = HenonHeilesFlow(); m = KahanLi(Float64, 8); tau = 0.05
    f(u, p, t) = SVector(u[3], u[4], -(u[1] + 2 * u[1] * u[2]), -(u[2] + u[1]^2 - u[2]^2))
    jac(u, p, t) = SMatrix{4, 4, Float64}(0, 0, -(1 + 2 * u[2]), -2 * u[1], 0, 0, -2 * u[1], -(1 - 2 * u[2]), 1, 0, 0, 0, 0, 1, 0, 0)
    zC, _ = energy_start(flow, 1 / 8, -0.015, 0.0)
    z6 = hh_chaotic_start(flow, m, 1 / 6, StableRNG(hh_seed(13, 1)))
    say("# O4: ChaosTools $(pkgversion(CT)), CoupledODEs with the default solver of DynamicalSystemsBase (Tsit5, order 5; no higher-order explicit solver is in the verify environment, none added), abstol = reltol = 1e-14; hand-coded Jacobian; renormalization every time unit")
    for (name, z0, sd) in (("C (E = 1/8)", zC, 2), ("chaotic start at E = 1/6", z6, 5))
        u0 = SVector{4, Float64}(z0)
        ode = CT.CoupledODEs(f, u0, nothing; diffeq = (abstol = 1.0e-14, reltol = 1.0e-14)); tds = CT.TangentDynamicalSystem(ode; J = jac, k = 4)
        # the stock lyapunovspectrum steps with step!(tds, 1) without stop_at_Δt: each step overshoots the time unit, so it averages over a slightly longer time with shifted renormalization instants
        lc = collect(CT.lyapunovspectrum(tds, 100; u0 = u0)); tel = CT.current_time(tds)
        li = hh_lyap(m, flow, z0, tau, 100.0).lambda                                      # identity frame (default W0), the frame of lyapunovspectrum
        say(@sprintf("# O4 %s: stock lyapunovspectrum (identity frame) (%s) over elapsed time %.4f (not 100); Chaoslab identity frame T = 100 (%s)", name, join(fmt.(lc), ", "), tel, join(fmt.(li), ", ")))
        emit("13", "O4: $(name): stock ChaosTools.lyapunovspectrum (identity frame, elapsed time $(round(tel; digits = 3)) instead of 100) against Chaoslab at T = 100, max |lambda_i - reference| (informational: not the same time window)", maximum(abs.(lc .- li)), "info", "INFO")
        function ct_spectrum(Q0)
            CT.reinit!(tds, u0; Q0 = Q0); lam = zeros(4)
            for _ in 1:100
                CT.step!(tds, 1.0, true); Q, Rm = qr(CT.current_deviations(tds))      # stop_at_Δt = true: exactly one time unit per step
                for j in 1:4
                    lam[j] += log(abs(Rm[j, j]))
                end
                CT.set_Q_as_deviations!(tds, Q)
            end
            return lam ./ 100, CT.current_time(tds)
        end
        lci, tI = ct_spectrum(one(SMatrix{4, 4, Float64}))
        d1 = maximum(abs.(lci .- li))
        emit("13", "O4: $(name), T = 100 exactly (stop_at_Δt), ChaosTools identity frame against Chaoslab identity frame: max |lambda_i - reference|", d1, "<= 1e-6", pf(d1 <= 1.0e-6 && tI == 100.0))
        W0 = initial_deviations(hh_seed(13, sd), 1, Val(4), Val(4), Float64)
        lam, tH = ct_spectrum(SMatrix{4, 4, Float64}(W0))
        lh = hh_lyap(m, flow, z0, tau, 100.0; W0 = W0).lambda
        d2 = maximum(abs.(lam .- lh))
        say(@sprintf("# O4 %s: ChaosTools Haar frame exact time (%s), Chaoslab (%s)", name, join(fmt.(lam), ", "), join(fmt.(lh), ", ")))
        emit("13", "O4: $(name), T = 100 exactly, seeded Haar frame (same as O3d): ChaosTools against Chaoslab, max |lambda_i - reference|", d2, "<= 1e-6", pf(d2 <= 1.0e-6 && tH == 100.0))
    end
end

# ------------------------------------------------------------------------------------------------ O5: portrait at E = 1/8 (criterion 10)
section("O5") do
    flow = HenonHeilesFlow(); m = KahanLi(Float64, 8); E = 1 / 8; N = 200
    t1 = @elapsed p = hh_portrait(m, E, N, 2000.0)
    say("# O5: $(N) x $(N) cell-centre grid over [a, b] x [-sqrt(2E), sqrt(2E)] kept inside D(E) (energy_grid), $(p.n) points, T_f = 2000, order-8 table, tau = 1/20, SALI frame (1,0,0,0), (0,0,1,0); $(round(t1; digits = 1)) s (the paper: 500 x 500 preregistered full size, not run here)")
    serialize(joinpath(SCR, "portrait-e8.jls"), (q2 = [z[2] for z in p.starts], p2 = [z[4] for z in p.starts], sali = p.sali, n = p.n))
    emit("10", "O5: chaotic fraction at E = 1/8 ($(N) x $(N), $(p.n) points, SALI(2000) < 1e-8), against the paper's about 55 percent", p.fraction, "within 0.10 of 0.55", pf(abs(p.fraction - 0.55) <= 0.10))
    for thr in (1.0e-6, 1.0e-10, 1.0e-12)
        emit("10", @sprintf("O5: chaotic fraction at E = 1/8 for the threshold %.0e (sensitivity, same data)", thr), count(<(thr), p.sali) / p.n, "info", "INFO")
    end
    emit("10", "O5: share of points with SALI(2000) in [1e-8, 1] (not classified by a clean gap, sticky and slow-mixing) ", count(s -> 1.0e-8 <= s < 1.0e-3, p.sali) / p.n, "info", "INFO")
    zR, _ = energy_start(flow, E, 0.1, 0.0); zS, _ = energy_start(flow, E, -0.0135, 0.0); zC, _ = energy_start(flow, E, -0.015, 0.0)
    emit("10", "O5: orbit R, p1 rebuilt against printed 0.490578 (difference)", abs(zR[3] - 0.490578), "<= 1e-6", pf(abs(zR[3] - 0.490578) <= 1.0e-6))
    emit("10", "O5: orbit S, p1 rebuilt against printed 0.499816 (difference)", abs(zS[3] - 0.499816), "<= 1e-6", pf(abs(zS[3] - 0.499816) <= 1.0e-6))
    emit("10", "O5: orbit C, p1 rebuilt against printed 0.499773 (difference)", abs(zC[3] - 0.499773), "<= 1e-6", pf(abs(zC[3] - 0.499773) <= 1.0e-6))
    sR = hh_sali(m, flow, zR, 2000.0)
    emit("10", "O5: orbit R, SALI(2000) (printed about 1.25)", sR, "in [0.5, 2]", pf(0.5 <= sR <= 2))
    sC = [hh_sali(m, flow, zC, T) for T in (1000.0, 2000.0, 4000.0, 10000.0)]
    emit("10", "O5: orbit C, SALI(2000)", sC[2], "< 1e-8", pf(sC[2] < 1.0e-8))
    emit("10", "O5: orbit C, SALI at T = 1000, 2000, 4000, 10^4: minimum", minimum(sC), "<= 1e-12", pf(minimum(sC) <= 1.0e-12))
    emit("10", "O5: orbit S (reported only), SALI at T = 2000, 5000, 10000 (printed: sticky, chaotic after about 5000 steps)", join(fmt.([hh_sali(m, flow, zS, T) for T in (2000.0, 5000.0, 10000.0)]), ", "), "report", "INFO")
    for (nm, z) in (("R", zR), ("S", zS), ("C", zC))
        res = [hh_sali(m, flow, z, T) for T in (100.0, 500.0, 1000.0, 2000.0)]
        say("# O5 SALI of orbit $(nm) at T = 100, 500, 1000, 2000: " * join(fmt.(res), ", "))
    end
end

function hh_sector_ref(z)
    th = atan(z[2], z[1]); c = (pi / 2, 7pi / 6, 11pi / 6)
    return argmin([abs(mod(th - ci + pi, 2pi) - pi) for ci in c]) - 1
end

# ------------------------------------------------------------------------------------------------ O6: escape continuation and tau/8 (criterion 12a)
section("O6") do
    flow = HenonHeilesFlow(); m = KahanLi(Float64, 8); tau = 0.05; R = 10.0; n = round(Int, 200 / tau); NST = 1000
    et = hh_eps(Float64)
    setprecision(BigFloat, BITS) do
        for (k, E) in enumerate((0.2, 0.25))
            rng = StableRNG(hh_seed(12, k)); nesc = 0; nret = 0; nrad = 0; nsec = 0; nstuck = 0; rminw = Inf
            for _ in 1:NST
                z0, _ = random_energy_start(flow, E, (-1.0, 1.0), (-1.0, 1.0), rng)
                res = escape_orbit(m, flow, z0, tau, n; R_esc = R)
                res.status === :escaped || continue
                nesc += 1
                hh_rad(res.z) >= R * (1 - 16 * et) || (nrad += 1)
                # continuation by the Taylor reference from the declared state (exact in 256 bits) to r = 2 R_esc: minimum radius over the nodes
                x = BigFloat.(collect(res.z)); rmin = BigFloat(Inf); reached = false
                for _ in 1:60
                    sol = tay(hh_f!, x, 0.25; maxsteps = 20000)
                    for i in 1:size(sol.x, 1)
                        r = sqrt(sol.x[i, 1]^2 + sol.x[i, 2]^2); rmin = min(rmin, r)
                        r > 2R && (reached = true; break)
                    end
                    reached && break
                    x = copy(sol.x[end, :])
                end
                reached || (nstuck += 1)
                rminw = min(rminw, Float64(rmin))
                rmin >= R / 2 || (nret += 1)
                (res.sector == hh_sector_ref(res.z)) || (nsec += 1)
            end
            emit("12", @sprintf("O6: E = %.2f, %d starts, %d escaped: r at the flag below R_esc (1 - 16 eps)", E, NST, nesc), nrad, "0", pf(nrad == 0))
            emit("12", @sprintf("O6: E = %.2f: Taylor continuation to r = 2 R_esc: orbits that reach 2 R_esc / orbits returning below R_esc / 2 (minimum r seen %.4g)", E, rminw), string(nesc - nstuck, "/", nesc, " ; ", nret), "all reach; 0 returns", pf(nstuck == 0 && nret == 0))
            emit("12", @sprintf("O6: E = %.2f: sector against the polar-angle rule (atan2, independent)", E), nsec, "0", pf(nsec == 0))
        end
    end
    # tau/8 self-convergence of the escape time and the label
    for (k, E) in enumerate((0.2, 0.25))
        rng = StableRNG(hh_seed(12, k)); worst = 0.0; nlab = 0; ncomp = 0; nstat = 0; rows = Tuple{Float64, Float64}[]
        for _ in 1:NST
            z0, _ = random_energy_start(flow, E, (-1.0, 1.0), (-1.0, 1.0), rng)
            r1 = escape_orbit(m, flow, z0, tau, n; R_esc = R); r8 = escape_orbit(m, flow, z0, tau / 8, 8n; R_esc = R)
            r1.status === r8.status || (r1.t < 100 || r8.t < 100) && (nstat += 1)
            (r1.status === :escaped && r1.t < 100 && r8.status === :escaped) || continue
            ncomp += 1; worst = max(worst, abs(r1.t - r8.t)); push!(rows, (r1.t, abs(r1.t - r8.t))); r1.sector == r8.sector || (nlab += 1)
        end
        emit("12", @sprintf("O6: E = %.2f: tau = 1/20 against tau/8, %d escapes before t = 100: max |t_escape(tau) - t_escape(tau/8)|", E, ncomp), worst, "<= 1e-6", pf(worst <= 1.0e-6))
        for (lo, hi) in ((0.0, 25.0), (25.0, 50.0), (50.0, 100.0))
            sub = filter(r -> lo <= r[1] < hi, rows)
            isempty(sub) || emit("12", @sprintf("O6: E = %.2f: escapes with t in [%g, %g): count, max |t(tau) - t(tau/8)|", E, lo, hi), string(length(sub), " ; ", fmt(maximum(r[2] for r in sub))), "<= 1e-6 (below t = 50: spec point U1 reading)", pf(maximum(r[2] for r in sub) <= 1.0e-6))
        end
        emit("12", @sprintf("O6: E = %.2f: label differences, status differences (tau against tau/8)", E), string(nlab, " ; ", nstat), "0 ; 0", pf(nlab == 0 && nstat == 0))
    end
end

# ------------------------------------------------------------------------------------------------ AC16T: 1 thread against N threads, bitwise
function portrait_data()
    p = hh_portrait(KahanLi(Float64, 8), 1 / 8, 200, 2000.0)
    return (x = p.x, W = p.W, sali = p.sali, chaotic = p.chaotic, n = p.n)
end
section("AC16child") do
    serialize(get(ENV, "HH_CHILD_OUT", joinpath(SCR, "child.jls")), portrait_data())
end
section("AC16T") do
    child = joinpath(SCR, "child-t1.jls")
    cmd = addenv(`$(Base.julia_cmd()) -t 1 --project=$(joinpath(ROOT, "verify")) $(@__FILE__)`, "HH_ONLY" => "AC16child", "HH_CHILD_OUT" => child, "HH_LOG" => joinpath(SCR, "child-t1.log"), "HH_SCRATCH" => SCR)
    t1 = @elapsed mine = portrait_data()
    t2 = @elapsed run(cmd)
    theirs = deserialize(child)
    say(@sprintf("# AC16T: 200 x 200 portrait at E = 1/8, T_f = 2000: %d threads %.1f s, 1 thread (child) %.1f s, %d points", Threads.nthreads(), t1, t2, mine.n))
    emit("16", "AC16T: final states bitwise identical with $(Threads.nthreads()) and 1 CPU threads ($(mine.n) starts)", isequal(mine.x, theirs.x), "bitwise", pf(isequal(mine.x, theirs.x)))
    emit("16", "AC16T: frames bitwise identical", isequal(mine.W, theirs.W), "bitwise", pf(isequal(mine.W, theirs.W)))
    emit("16", "AC16T: SALI values bitwise identical", isequal(mine.sali, theirs.sali), "bitwise", pf(isequal(mine.sali, theirs.sali)))
    emit("16", "AC16T: class labels identical (disagreeing points)", count(mine.chaotic .!= theirs.chaotic), "0", pf(mine.chaotic == theirs.chaotic))
    emit("16", "AC16T: GPU against CPU (chaotic fraction within 1 point, class disagreement at most 2 percent)", "not run", "-", "NOT RUN: needs verify/gpu.jl (writes data/); owner decision")
end

# ------------------------------------------------------------------------------------------------ FIG
section("FIG") do
    flow = HenonHeilesFlow(); E = 1 / 8; a, b = section_interval(flow, E); pm = sqrt(2E); T = 1000; nst = 100
    rng = StableRNG(VSEED + 6)
    starts = [random_energy_start(flow, E, (a, b), (-pm, pm), rng)[1] for _ in 1:nst]
    taus = [1 / 2, 1 / 3, 1 / 4, 1 / 5, 1 / 10, 1 / 20, 1 / 40, 1 / 80]; res = Dict{String, Vector{Float64}}()
    for (name, m) in (("order-8 (KahanLi)", KahanLi(Float64, 8)), ("order-4 (Forest-Ruth)", ForestRuth(Float64)))
        w = Float64[]
        for tau in taus
            ws = zeros(length(starts))
            Threads.@threads for i in eachindex(starts)
                z = starts[i]; H0 = hh_H(z); mx = 0.0
                for _ in 1:round(Int, T / tau)
                    z = advance(m, flow, z, tau); mx = max(mx, abs(hh_H(z) - H0))
                end
                ws[i] = mx
            end
            push!(w, maximum(ws))
        end
        res[name] = w
        emit("6", "FIG: $(name), E = 1/8, $(nst) seeded starts, T = 1000: max |H - H0| at tau = 1/2, 1/3, 1/4, 1/5, 1/10, 1/20, 1/40, 1/80", join(fmt.(w), ", "), "tau = 1/20: <= $(name[1:7] == "order-8" ? "1e-8" : "1e-4")", pf(w[6] <= (name[1:7] == "order-8" ? 1.0e-8 : 1.0e-4)))
        for i in 1:(length(taus) - 1)
            w[i + 1] > 1.0e-12 && emit("6", "FIG: $(name): observed order log(err_i / err_(i+1)) / log(tau_i / tau_(i+1)) between tau = 1/$(round(Int, 1 / taus[i])) and 1/$(round(Int, 1 / taus[i + 1]))", log(w[i] / w[i + 1]) / log(taus[i] / taus[i + 1]), "about $(name[1:7] == "order-8" ? 8 : 4) (pairs above the 1e-12 round-off floor only)", "INFO")
        end
    end
    fig = CM.Figure(size = (820, 560))
    ax = CM.Axis(fig[1, 1]; xscale = log10, yscale = log10, xlabel = "step size tau", ylabel = "max |H - H0|, T = 1000, E = 1/8", title = "Henon-Heiles energy error (100 seeded starts)")
    CM.scatterlines!(ax, taus, res["order-8 (KahanLi)"]; label = "order-8 table", marker = :circle)
    CM.scatterlines!(ax, taus, res["order-4 (Forest-Ruth)"]; label = "order-4 table", marker = :rect)
    CM.lines!(ax, taus, res["order-8 (KahanLi)"][1] .* (taus ./ taus[1]) .^ 8; linestyle = :dash, color = :gray, label = "slope 8")
    CM.lines!(ax, taus, res["order-4 (Forest-Ruth)"][1] .* (taus ./ taus[1]) .^ 4; linestyle = :dot, color = :gray, label = "slope 4")
    CM.hlines!(ax, [1.0e-8, 1.0e-4]; color = :red, linestyle = :dashdot, label = "bounds 1e-8, 1e-4 (tau = 1/20)")
    CM.axislegend(ax; position = :rb)
    CM.save(joinpath(OUT, "henon-heiles-energy.png"), fig)
    pdat = deserialize(joinpath(SCR, "portrait-e8.jls"))
    fig2 = CM.Figure(size = (780, 720))
    ax2 = CM.Axis(fig2[1, 1]; xlabel = "q2", ylabel = "p2", title = "Henon-Heiles E = 1/8: log10 SALI(2000) on the section q1 = 0, p1 > 0 ($(pdat.n) starts, order-8 table)", titlesize = 12)
    sc = CM.scatter!(ax2, pdat.q2, pdat.p2; color = log10.(max.(pdat.sali, 1.0e-16)), colormap = :viridis, colorrange = (-16, 0.3), markersize = 3)
    CM.Colorbar(fig2[1, 2], sc; label = "log10 SALI(2000); chaotic: below -8")
    zR, _ = energy_start(flow, E, 0.1, 0.0); zS, _ = energy_start(flow, E, -0.0135, 0.0); zC, _ = energy_start(flow, E, -0.015, 0.0)
    for (nm, z, col) in (("R", zR, :red), ("S", zS, :orange), ("C", zC, :magenta))
        CM.scatter!(ax2, [z[2]], [z[4]]; color = col, marker = :xcross, markersize = 14, label = "orbit $(nm)")
    end
    CM.axislegend(ax2; position = :rt)
    CM.save(joinpath(OUT, "henon-heiles-portrait.png"), fig2)
    say("# FIG: wrote henon-heiles-energy.png and henon-heiles-portrait.png in the stage 04 output folder")
end
