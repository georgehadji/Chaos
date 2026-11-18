# Stage 04 verification of module ofli (catalog row I6: Chaoslab.Indicators ofli.jl; spec D118).
# Run from the repo root, one process per group (OF_ONLY; OF_LOG names the log, OF_APPEND=1 appends):
#   OF_ONLY=FULL1 .. FULL4  julia --project=verify -t 8 verify/ofli.jl   (the package test files of test/ofli/ at the FULL sizes, OF_FULL_RUN = true; one file each; FULL = all four)
#   OF_ONLY=Oa      (BigFloat 256-bit closed forms: free particle, oscillator and pendulum identities, period-1 orbit and monodromy: criteria 1, 2, 3a, 4a, 5)
#   OF_ONLY=Ob      (same-orbit BigFloat reference, singular values of P Y, sums against the FLI and SALI reducers: criteria 9, 10, 12, 15)
#   OF_ONLY=Oc      (TaylorIntegration BigFloat references of F3, orbits R and C: criterion 17)
#   OF_ONLY=PLATEAU (plateau-loss time on a pendulum start that is not a turning point and on the F3 period-1 orbit: R1, criteria 4c, 5c, 11c)
#   OF_ONLY=Od      (ChaosTools Lyapunov cross-check of the sums, INFO)
#   OF_ONLY=R       (R2, R4, R5 recorded lines)
#   OF_ONLY=AC16T   (1 against N CPU threads, bitwise; spawns a one-thread child)
#   OF_ONLY=FIG     (ofli-plateau.png, ofli-classes.png)
# Not run: criterion 18 device part (verify/gpu.jl writes data/; owner decision). Nothing here writes under data/.
# Spec: workspaces/module-build/stages/01-spec/output/ofli-spec.md; test plan points U1 to U13 are the operative reading. Verify-only seeds 20261012 + 43900 ... (package seeds 43000 + 50 ac + sub).
# One printed line per measurement:  AC<n> | quantity | value | bound | PASS/FAIL/INFO
const ONLY = filter(!isempty, split(get(ENV, "OF_ONLY", ""), ","))
using Printf, Random, LinearAlgebra, Statistics, Serialization, Logging
using Test
const OF_FULL_RUN = true
using Chaoslab, StaticArrays, DoubleFloats, StableRNGs, TOML
const NEEDTAY = any(o -> startswith(o, "Oc") || startswith(o, "Oa") || startswith(o, "PLATEAU") || startswith(o, "FIG") || startswith(o, "R"), ONLY)
NEEDTAY && (import TaylorIntegration as TI)
"Od" in ONLY && (import ChaosTools as CT)
"FIG" in ONLY && (import CairoMakie as CM)

const ROOT = dirname(@__DIR__)
const OUT = joinpath(ROOT, "workspaces", "module-build", "stages", "04-verify", "output")
const SCR = get(ENV, "OF_SCRATCH", joinpath(tempdir(), "of-verify-scratch"))
mkpath(SCR)
const LOG = open(get(ENV, "OF_LOG", joinpath(OUT, "ofli-verify-log.txt")), get(ENV, "OF_APPEND", "") == "1" ? "a" : "w")

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

say("# ofli verification, julia $(VERSION), threads $(Threads.nthreads()), base seed 20261012 + 43000, OF_FULL_RUN = true, groups $(isempty(ONLY) ? "all" : join(ONLY, ",")), time $(round(time())) (unix)")
say("# timings in this log were taken on a machine that may be loaded by other processes")

const TD = joinpath(ROOT, "test", "ofli")
include(joinpath(TD, "of_helpers.jl"))
const VSEED = 20261012 + 43900
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
suite_section("FULL1", "of_ac01_05.jl")      # criteria 1 to 5
suite_section("FULL2", "of_ac06_10.jl")      # criteria 6 to 10
suite_section("FULL3", "of_ac11_15.jl")      # criteria 11 to 15
suite_section("FULL4", "of_ac16_18.jl")      # criteria 16 to 18 (CPU part)

# ---------------------------------------------------------------------------------------------------
# Taylor helpers (hand-written fields, independent of Chaoslab's): F3 state + 4 variational columns, 1 dof state + 1 variational column
bV(q1, q2) = (q1^2 + q2^2) / 2 + q1^2 * q2 - q2^3 / 3
bu(q2) = q2^2 / 2 - q2^3 / 3
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
tay(f!, x0, T; order = 40, maxsteps = 10^6, kw...) = TI.taylorinteg(f!, x0, zero(BigFloat), BigFloat(T), order, ldexp(BigFloat(1), -200); maxsteps = maxsteps, kw...)
# the dense Taylor solution at time t (BigFloat)
function tay_eval(sol, t)
    i = clamp(searchsortedlast(sol.t, t), 1, length(sol.t) - 1)
    s = t - sol.t[i]
    return [sol.p[i, k](s) for k in 1:size(sol.x, 2)]
end
function poly_root(g, dt)
    lo = zero(dt); hi = dt
    for _ in 1:250
        mid = (lo + hi) / 2
        g(mid) < 0 ? (lo = mid) : (hi = mid)
    end
    return (lo + hi) / 2
end
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
# first return to the section q1 = 0, p1 > 0 at energy E from (q2, p2) on it: (elapsed time, q2', p2'), 256 bits
function tay_ret(E, q2, p2)
    p1 = sqrt(2 * (E - p2^2 / 2 - bu(q2))); x = [BigFloat(0), q2, p1, p2]; tacc = BigFloat(0)
    for _ in 1:40
        sol = tay(hh_f!, x, 8)
        c = tay_section(sol; start = false)
        isempty(c.t) || return (tacc + c.t[1], c.q2[1], c.p2[1])
        tacc += sol.t[end]; x = copy(sol.x[end, :])
    end
    return (BigFloat(NaN), BigFloat(NaN), BigFloat(NaN))
end
# Taylor series of the exact flow at integer times t = 0 .. T for the frame W0 (4 columns): GB, G1 (projected running maxima incl. the s = 0 term), FB (unprojected running max of ln ||w_i||, s = 0 term 0)
function tay_series(z0, W0, T::Integer; sol = nothing)
    return setprecision(BigFloat, BITS) do
        sol === nothing && (sol = tay(var_f!, vcat(BigFloat.(collect(z0)), BigFloat.(vec(Matrix(W0)))), T))
        GB = zeros(T + 1); G1 = zeros(T + 1); FB = zeros(T + 1); gB = -BigFloat(Inf); g1 = -BigFloat(Inf); fB = BigFloat(0)
        for t in 0:T
            x = tay_eval(sol, BigFloat(t)); z = x[1:4]
            f = [z[3], z[4], -(z[1] + 2 * z[1] * z[2]), -(z[2] + z[1]^2 - z[2]^2)]; ff = dot(f, f)
            for i in 1:4
                w = x[(4i + 1):(4i + 4)]; pw = w - (dot(f, w) / ff) * f
                g = log(norm(pw)); S = log(norm(w))
                g > gB && (gB = g); i == 1 && g > g1 && (g1 = g); S > fB && (fB = S)
            end
            GB[t + 1] = Float64(gB); G1[t + 1] = Float64(g1); FB[t + 1] = Float64(fB)
        end
        (GB = GB, G1 = G1, FB = FB, nsteps = length(sol.t) - 1)
    end
end
# the period-1 orbit of F3 at E = 1/8 (Newton on the first return map of the package's section_map, Float64): z0, T, reduced monodromy A
function vperiodic(E, q2, p2; h = 1 / 80)
    flow = HenonHeilesFlow(); m = KahanLi(Float64, 8)
    ret(q2, p2) = (z = energy_start(flow, E, q2, p2)[1]; z === nothing ? nothing : section_map(m, flow, z, h, 4000; jacobian = true, max_crossings = 2))
    x = SVector(q2, p2)
    for _ in 1:30
        r = ret(x[1], x[2]); (r === nothing || length(r.z) < 2) && return nothing
        y = SVector(r.z[2][2], r.z[2][4]); F = y - x
        if maximum(abs, F) <= 1.0e-13
            return (z0 = energy_start(flow, E, x[1], x[2])[1], T = r.t[2], A = r.J[2], closure = maximum(abs, F), x = x)
        end
        x = x - (r.J[2] - one(r.J[2])) \ F
    end
    return nothing
end
const PO = Ref{Any}(nothing)
periodic_orbit() = (PO[] === nothing && (PO[] = vperiodic(1 / 8, 0.1, 0.0)); PO[])
tstar(dep, thr = 0.05) = (k = findfirst(>(thr), dep); k === nothing ? Inf : Float64(k - 1))     # dep[t + 1] at integer time t

# ------------------------------------------------------------------------------------------------ Oa: BigFloat(256) closed forms
section("Oa") do
    # (1) free particle d = 2, 4 in BigFloat(256) (exact dyadic step): criteria 1, 2
    setprecision(BigFloat, BITS) do
        eT = Float64(effective_eps(BigFloat)); N = 10_000
        say("# Oa: free particle in BigFloat(256), h = 1/16, N = $N steps, p0 = (0.6, 0.8) for d = 4 and 0.5 for d = 2, canonical frame, one single-column run per column (OFLI_1 of each column)")
        for (Nd, p0f) in ((2, (0.6, 0.8)), (1, (0.5,))), tau in (1, 10)
            flow = of_FreeFlow{Nd}(); d = 2Nd; h = BigFloat(1) / 16; sys = of_FreeMap(flow, h)
            z0 = SVector{d, BigFloat}(vcat(zeros(BigFloat, Nd), of_big.(collect(p0f))))
            f0 = vcat(of_big.(collect(p0f)), zeros(BigFloat, Nd)); ff = dot(f0, f0)
            worst = 0.0; nconst = 0; ngrow = 0
            for i in 1:d
                W = SMatrix{d, 1, BigFloat}(ntuple(j -> j == i ? BigFloat(1) : BigFloat(0), d))
                snaps = of_snaps(OFLI(flow), sys, z0, W, N, tau); g0 = ofli_initial(flow, z0, W)[1]
                cf(n) = (t = n * h; j = i - Nd;
                         w = i <= Nd ? Vector(W[:, 1]) : vcat([(l == j ? t : BigFloat(0)) for l in 1:Nd], [(l == j ? BigFloat(1) : BigFloat(0)) for l in 1:Nd]);
                         pw = w - (dot(f0, w) / ff) * f0; nn = norm(pw); nn > 0 ? log(nn) : -BigFloat(Inf))
                gr = cf(0)
                (isfinite(gr) == isfinite(g0)) || (worst = Inf)
                for s in snaps
                    gr = max(gr, cf(s.n))
                    num = max(g0, s.ofli1)
                    if isfinite(gr) && isfinite(num)
                        worst = max(worst, Float64(abs(num - gr)) / (32 * eT * max(1.0, abs(Float64(gr)))))
                        s.n in (100, 1000, 10000) && say(@sprintf("# Oa:   d = %d tau = %d column %d n = %d: |error| / eps = %.3g", d, tau, i, s.n, Float64(abs(num - gr)) / eT))
                    else
                        isfinite(gr) == isfinite(num) || (worst = Inf)
                    end
                end
                last(snaps).ofli1 > g0 + 1 ? (ngrow += 1) : (nconst += 1)
            end
            emit("1", "Oa: free particle d = $d, BigFloat(256), tau = $tau, N = $N: worst error / (32 eps max(1, |g|), Amendment D119) over $d single columns ($nconst constant, $ngrow growing like ln t; minus-infinity column agrees with the closed form)", worst, "<= 1", pf(worst <= 1))
        end
    end
    # (2) oscillator F1 (omega = 1.3) and pendulum F4: the exact-flow identity ||w_perp|| = |gradH(z0).w(0)| / ||f(z(t))|| on a TaylorIntegration BigFloat(256) reference, 100 sample times
    setprecision(BigFloat, BITS) do
        w0 = (BigFloat(6) / 10, BigFloat(8) / 10)
        function ident(name, V1, V2, z0, ac)
            fv!(dx, x, p, s) = (dx[1] = x[2]; dx[2] = -V1(x[1]); dx[3] = x[4]; dx[4] = -V2(x[1]) * x[3]; nothing)
            sol = tay(fv!, [z0[1], z0[2], w0[1], w0[2]], 100; parse_eqs = false)
            gH = [V1(z0[1]), z0[2]]; c = abs(gH[1] * w0[1] + gH[2] * w0[2]); worst = BigFloat(0)
            for t in range(0.37, 99.0, length = 100)
                x = tay_eval(sol, BigFloat(t)); f = [x[2], -V1(x[1])]; w = [x[3], x[4]]
                pw = w - (dot(f, w) / dot(f, f)) * f
                worst = max(worst, abs(norm(pw) - c / norm(f)) / (c / norm(f)))
            end
            emit(ac, "Oa: $name, Taylor order 40, abstol 2^-200, 256 bits, 100 sample times in [0.37, 99]: max relative |norm(P w) - |gradH0.w0| / norm(f)|", Float64(worst), "<= 1e-30", pf(worst <= BigFloat(10)^-30))
        end
        om2 = BigFloat(1.3)^2
        ident("oscillator F1 omega = 1.3, (q, p) = (1, 0.3)", q -> om2 * q, q -> om2, (BigFloat(1), BigFloat(3) / 10), "3")
        for A in (2.0, 3.0)
            ident("pendulum F4 amplitude $(Int(A)) (turning point start)", q -> sin(q), q -> cos(q), (BigFloat(A), BigFloat(0)), "4")
            ident("pendulum F4 amplitude $(Int(A)) from q = 0 (not a turning point)", q -> sin(q), q -> cos(q), (BigFloat(0), sqrt(2 * (1 - cos(BigFloat(A))))), "4")
        end
    end
    # (3) the period-1 orbit of F3 at E = 1/8 refined at 256 bits; reduced monodromy, T'(E); exact (Taylor) series for criterion 5(a)
    po = periodic_orbit()
    setprecision(BigFloat, BITS) do
        E = BigFloat(1) / 8; hd = ldexp(BigFloat(1), -100)
        R(E, x) = (r = tay_ret(E, x[1], x[2]); (r[1], [r[2], r[3]]))
        x = BigFloat.(collect(po.x)); J = zeros(BigFloat, 2, 2); F = x; T = BigFloat(0)
        for it in 1:7
            T, y = R(E, x); F = y - x
            say(@sprintf("# Oa: Newton iteration %d: closure |R(x) - x| = %.3g", it, Float64(maximum(abs, F))))
            maximum(abs, F) < BigFloat(10)^-70 && break
            cq = (R(E, x .+ [hd, 0])[2], R(E, x .- [hd, 0])[2]); cp = (R(E, x .+ [0, hd])[2], R(E, x .- [0, hd])[2])
            J = hcat((cq[1] - cq[2]) / (2hd), (cp[1] - cp[2]) / (2hd))
            x = x - (J - I) \ F
        end
        T, y = R(E, x); F = y - x
        trA = J[1, 1] + J[2, 2]; detA = J[1, 1] * J[2, 2] - J[1, 2] * J[2, 1]
        emit("5", "Oa: period-1 orbit of F3 at E = 1/8, 256 bits: closure |R(x) - x| after Newton", Float64(maximum(abs, F)), "<= 1e-30", pf(maximum(abs, F) <= BigFloat(10)^-30))
        emit("5", "Oa: period T (section time of the orbit), 256 bits", Float64(T), "report", "INFO")
        emit("5", "Oa: reduced monodromy A: trace, det - 1 (A2: |tr A| < 2)", string(fmt(Float64(trA)), " ; ", fmt(Float64(detA - 1))), "|tr A| < 2, |det - 1| < 1e-6", pf(abs(trA) < 2 && abs(detA - 1) < 1.0e-6))
        # T'(E) by central difference with the fixed point re-solved at E +- dE (chord Newton with the J of E)
        dE = ldexp(BigFloat(1), -30); Tp = BigFloat[]
        for s in (+1, -1)
            Es = E + s * dE; xs = copy(x); Ts = BigFloat(0)
            for _ in 1:8
                Ts, y = R(Es, xs); xs = xs - (J - I) \ (y - xs)
            end
            push!(Tp, R(Es, xs)[1])
        end
        Tprime = (Tp[1] - Tp[2]) / (2dE)
        emit("5", "Oa: T'(E) at E = 1/8 (central difference, dE = 2^-30, fixed point re-solved at 256 bits)", Float64(Tprime), "not zero (O6 e); package AC5: -1.347", pf(abs(Tprime) > 1.0e-3))
        # exact series: 100 T, canonical and seeded frame (criterion 5a)
        z0b = begin
            p2 = x[2]; q2 = x[1]; p1 = sqrt(2 * (E - p2^2 / 2 - bu(q2))); [BigFloat(0), q2, p1, p2]
        end
        Tf = Float64(T); n10, n100 = round(Int, 10Tf), round(Int, 100Tf)
        for (fname, W) in (("canonical", SMatrix{4, 4, Float64}(I)), ("seeded (of_frame)", of_frame(VSEED + 5, 1, 4, 4, Float64)))
            t0 = time()
            r = tay_series(z0b, W, n100)
            dG = r.GB[n100 + 1] - r.GB[n10 + 1]; dF = r.FB[n100 + 1] - r.FB[n10 + 1]
            say(@sprintf("# Oa: exact series, %s frame, %d Taylor steps, %.1f s", fname, r.nsteps, time() - t0))
            emit("5a", "Oa: exact (Taylor 256 bits) series on the period-1 orbit, $fname frame, OFLI_B(100 T) - OFLI_B(10 T) with 10 T -> $n10, 100 T -> $n100", dG, "<= 0.1", pf(dG <= 0.1))
            emit("5a", "Oa: exact series, $fname frame, FLI_B(100 T) - FLI_B(10 T)", dF, ">= 1.5 (ln 10 = 2.30 expected)", pf(dF >= 1.5))
        end
    end
end

# ------------------------------------------------------------------------------------------------ Ob: same-orbit BigFloat reference, basis-free value from the singular values of P Y
# ln of the largest singular value of a 4 x 4 matrix B in BigFloat: cyclic Jacobi eigenvalues of B'B (test-side, independent of the reducer's per-column norms)
function lnsigma1(B::Matrix{BigFloat})
    A = B' * B; n = size(A, 1); scale = maximum(abs, diag(A)); tol = scale * ldexp(BigFloat(1), -(BITS - 24))
    for _ in 1:60
        sum(abs2, (A[i, j] for i in 1:n for j in (i + 1):n)) <= tol^2 && break
        for p in 1:(n - 1), q in (p + 1):n
            A[p, q] == 0 && continue
            th = (A[q, q] - A[p, p]) / (2 * A[p, q]); t = (th >= 0 ? 1 : -1) / (abs(th) + sqrt(th^2 + 1)); c = 1 / sqrt(t^2 + 1); s = t * c
            G = Matrix{BigFloat}(I, n, n); G[p, p] = c; G[q, q] = c; G[p, q] = s; G[q, p] = -s
            A = G' * A * G
        end
    end
    return log(maximum(diag(A))) / 2
end
function ob_fstar(flow, c, tau::Integer, N::Integer)       # F*(m tau) = running max of ln sigma_1(P Y) over the handoffs incl. m = 0 (frame free)
    return setprecision(BigFloat, BITS) do
        out = zeros(N ÷ tau + 1); run = -Inf
        for m in 0:(N ÷ tau)
            n = m * tau; f = collect(vector_field(flow, SVector{4, BigFloat}(map(of_big, c.pts[n + 1])))); P = Matrix{BigFloat}(I, 4, 4) - (f * f') / dot(f, f)
            run = max(run, Float64(lnsigma1(P * Matrix(c.Ys[n + 1])))); out[m + 1] = run
        end
        out
    end
end
function ob_run(keys_, N, tag)
    flow = HenonHeilesFlow(); tau = 10; cps = collect(100:100:N); ln2 = log(2.0)
    nbad = 0; nrefbad = 0; nrow = 0; maxratio = 0.0; ndiff = 0; fr = Float64[]; t0 = time()
    for key in keys_
        cc = of_case(key, Float64, :o8, N; frame = Matrix{Float64}(I, 4, 4)); cs = of_case(key, Float64, :o8, N)
        Fs = ob_fstar(flow, cc, tau, N)
        for c in (cc, cs)
            snaps = of_snaps(OFLI(flow), c.im, c.x0, c.W0, N, tau); g0 = ofli_initial(flow, c.x0, c.W0)
            rows = of_rows(c, Float64, tau, snaps, g0, cps); mx = of_refmax(c.ref.g, tau, N)
            for r in rows
                r.form === :B || continue
                m = r.n ÷ tau; G = max(maximum(g0), snaps[m].ofli); F = Fs[m + 1]
                (F - ln2 - r.bound <= G <= F + r.bound) || (nbad += 1)
                rg = Float64(mx.GB[m + 1]); (F - ln2 - 1.0e-12 <= rg <= F + 1.0e-12) || (nrefbad += 1)
                nrow += 1; maxratio = max(maxratio, r.ratio); r.maxdiff && (ndiff += 1)
                r.n == N && push!(fr, (F - G) / ln2)
            end
            for r in rows
                r.form === :one || continue
                nrow += 1; maxratio = max(maxratio, r.ratio); r.maxdiff && (ndiff += 1)
            end
        end
    end
    say(@sprintf("# Ob %s: %d orbits x 2 frames, N = %d steps, tau = %d, checkpoints every 100 steps, %.1f s", tag, length(keys_), N, tau, time() - t0))
    emit("9", "Ob $tag: sandwich F* - ln 2 - B_G <= OFLI_B <= F* + B_G, F* from the BigFloat(256) singular values of P Y (Jacobi, independent): violations over $nrow rows", nbad, "0", pf(nbad == 0))
    emit("9", "Ob $tag: the same sandwich on the BigFloat reference series itself (slack 1e-12): violations", nrefbad, "0", pf(nrefbad == 0))
    emit("9", "Ob $tag: (F* - OFLI_B) / ln 2 at the last step: min, median, max (must lie in [0, 1])", string(fmt(minimum(fr)), " ; ", fmt(of_median(fr)), " ; ", fmt(maximum(fr))), "[0, 1]", pf(all(x -> -1.0e-9 <= x <= 1 + 1.0e-9, fr)))
    emit("12", "Ob $tag: Float64, order-8, same-orbit BigFloat reference, both forms, $nrow rows: max error / B_G (R3)", maxratio, "<= 1", pf(maxratio <= 1))
    emit("12", "Ob $tag: rows whose maximiser differs between computed and reference series", ndiff, "report", "INFO")
end
section("Ob1") do
    say("# Ob1 REDUCTION: orbits R and C at N = 10000 steps (t = 500) with checkpoints every 100 steps (spec criterion 9 horizon list reaches 1e4 steps; the BigFloat product at N = 1e4 is the cost)")
    ob_run([:R, :C], 10_000, "R, C")
end
section("Ob2") do
    say("# Ob2 REDUCTION: 20 regular-class starts at N = 2000 steps (t = 100), tau = 10 (spec: 20 starts per class; horizon reduced from 1e4 steps for the cost of the BigFloat products)")
    ob_run([(:reg, i) for i in 1:20], 2000, "regular class")
end
section("Ob3") do
    say("# Ob3 REDUCTION: 20 chaotic-class starts at N = 2000 steps (t = 100), tau = 10 (spec: 20 starts per class; horizon reduced from 1e4 steps for the cost of the BigFloat products)")
    ob_run([(:cha, i) for i in 1:20], 2000, "chaotic class")
end

# ------------------------------------------------------------------------------------------------ Oc: TaylorIntegration BigFloat(256) reference of F3, orbits R and C (criterion 17)
section("Oc") do
    flow = HenonHeilesFlow(); delta8 = 1.0e-6
    say("# Oc: reference = Taylor order 40, abstol 2^-200, 256 bits, state + 4 variational columns (hand-written field), OFLI_B from O1 with the exact field at the reference points at integer times; frame = seeded of_frame; numerics: order-8 (KahanLi) and order-4 (ForestRuth) tables, h = 1/20, tau = 20")
    for (name, z0, seed, T) in (("R", of_orbit_R(), VSEED + 1, 1000), ("C", of_orbit_C(), VSEED + 2, 60))
        W0 = of_frame(seed, 1, 4, 4, Float64); t0 = time()
        ref = tay_series(z0, W0, T)
        say(@sprintf("# Oc: orbit %s, reference to t = %d: %d Taylor steps, %.1f s", name, T, ref.nsteps, time() - t0))
        g0 = ofli_initial(flow, z0, W0); err = Dict{Symbol, Function}()
        for key in (:o8, :o4)
            snaps = of_snaps(OFLI(flow), of_im(Float64, key), z0, W0, 20 * T, 20)
            err[key] = t -> abs(max(maximum(g0), snaps[t].ofli) - ref.GB[t + 1])
        end
        for t in (10, 30, 60)
            e8 = err[:o8](t); e4 = err[:o4](t)
            emit("17", "Oc: orbit $name, t = $t, order-8 table h = 1/20: |OFLI_B,num - OFLI_B,Taylor| (OFLI_B,Taylor = $(fmt(ref.GB[t + 1])))", e8, "<= delta_8 = 1e-6", pf(e8 <= delta8))
            emit("17", "Oc: orbit $name, t = $t, order-4 table h = 1/20: error (reported), ratio to the order-8 error", string(fmt(e4), " ; ", fmt(e4 / e8)), "report", "INFO")
        end
        if name == "R"
            for t in (100, 300, 1000)
                emit("17", "Oc: orbit R, t = $t: order-8 error ; order-4 error (reported)", string(fmt(err[:o8](t)), " ; ", fmt(err[:o4](t))), "report", "INFO")
            end
            emit("17", "Oc: orbit R to t = 1000: max over integer t of the order-8 error", maximum(err[:o8](t) for t in 1:T), "report", "INFO")
        end
    end
end

# ------------------------------------------------------------------------------------------------ PLATEAU: plateau-loss time t* (departure of the numerical OFLI from the exact one above 0.05)
# Reader storing one Float64 field per handoff only (memory for very long series); CPU, use with propagate_orbit
struct Probe1{R}
    inner::R
    field::Symbol
end
Chaoslab.reducer_init(r::Probe1, ::Type{T}, v::Val, k::Val) where {T} = (Chaoslab.reducer_init(r.inner, T, v, k), Float64[])
function Chaoslab.reducer_step(r::Probe1, st, n::Int, x, W)
    inner2, W2 = Chaoslab.reducer_step(r.inner, st[1], n, x, W)
    push!(st[2], Float64(getfield(Chaoslab.reducer_final(r.inner, inner2), r.field)))
    return (inner2, st[2]), W2
end
Chaoslab.reducer_failed(r::Probe1, st) = Chaoslab.reducer_failed(r.inner, st[1])
Chaoslab.reducer_error(r::Probe1, st, orbit, backend) = Chaoslab.reducer_error(r.inner, st[1], orbit, backend)
Chaoslab.reducer_final(r::Probe1, st) = (series = st[2],)
function light_series(r, field, sys, x0, W0, N, tau)
    p = Probe1(r, field); _, _, st = propagate_orbit(sys, x0, W0, N, tau, (p,))
    return Chaoslab.reducer_final(p, st[1]).series
end
function pend_start(A, kind)                              # kind :zero: q = 0, p from the energy (not a turning point); :turn: the turning point (A, 0)
    return kind === :zero ? SVector(0.0, sqrt(2 * (1 - cos(A)))) : SVector(A, 0.0)
end
function pend_lnf(z0, Tmax::Integer; per = 20)            # ln ||f|| at the integer times, Double64 order-10 reference at h = 1/per (state only, O5 needs no tangent)
    D = Double64; m = SofroniouSpaletta10(D); flow = Pendulum(D(1)); h = D(1) / per
    z = SVector{2, D}(D(z0[1]), D(z0[2])); out = zeros(Tmax + 1)
    for t in 0:Tmax
        out[t + 1] = Float64(log(norm(vector_field(flow, z))))
        t < Tmax && for _ in 1:per
            z = advance(m, flow, z, h)
        end
    end
    return out
end
function decimate(t, g; n = 1500)
    idx = unique(round.(Int, exp.(range(0, log(length(t)), length = n))))
    return (t = t[idx], g = g[idx])
end
section("PLATEAUp") do
    Tmax = 1_000_000; W = SMatrix{2, 1, Float64}(0.6, 0.8); flow = Pendulum(1.0); data = Dict{Any, Any}()
    say("# PLATEAUp REDUCTION: reference horizon Tmax = $Tmax time units, 10^6/(2 pi) periods or more (a t* above Tmax is reported as > Tmax); reference = Double64 order-10 table at h = 1/20 (state only), exact OFLI_1 from O5: G(t) = ln |gradH0.w0| - running min ln ||f(z(t))||")
    for A in (2.0, 3.0)
        z0 = pend_start(A, :zero); P = 2pi / (begin a, b = 1.0, cos(A / 2); for _ in 1:40; a, b = (a + b) / 2, sqrt(a * b); end; a end)
        t0 = time(); lnf = pend_lnf(z0, Tmax); say(@sprintf("# PLATEAUp: amplitude %g, start q = 0, p = %.6f, period P = %.6f, reference %.1f s", A, z0[2], P, time() - t0))
        gH = SVector(sin(z0[1]), z0[2]); c = log(abs(dot(gH, W[:, 1])))
        Gex = c .- accumulate(min, lnf)
        data[(:pend, A, :ref)] = decimate(collect(0.0:Tmax), Gex)
        res = Dict{Any, Float64}()
        for (key, hs) in ((:o4, (1 / 2, 1 / 4, 1 / 8, 1 / 16, 1 / 20)), (:o8, (1.0, 1 / 2, 1 / 4, 1 / 20)))
            p = key === :o4 ? 4 : 8
            for h in hs
                tau = round(Int, 1 / h); im = IntegratorMap(of_table(Float64, key), flow, h)
                ser = light_series(OFLI(flow), :ofli1, im, z0, W, Tmax * tau, tau); g0 = ofli_initial(flow, z0, W)[1]
                num = vcat(g0, [max(g0, v) for v in ser]); dep = abs.(num .- Gex)
                ts = tstar(dep); res[(key, h)] = ts
                data[(:ts, :pend, A, key, h)] = ts
                imax10 = min(Tmax, floor(Int, 10P))
                data[(:pend, A, key, h)] = decimate(collect(0.0:Tmax), num)
                emit("4c", "PLATEAUp: pendulum amplitude $(Int(A)), start q = 0, order-$p table, h = 1/$(round(Int, 1 / h)): t* (first t with |OFLI_1,num - OFLI_1,exact| > 0.05); max departure for t <= 10 P = $(fmt(imax10)); departure at Tmax", string(ts == Inf ? "> $Tmax" : fmt(ts), " ; ", fmt(maximum(dep[1:(imax10 + 1)])), " ; ", fmt(dep[end])), "report (R1)", "INFO")
                if key === :o4 && h == 1 / 20
                    e = maximum(dep[1:(imax10 + 1)])
                    emit("4c", "PLATEAUp: pendulum amplitude $(Int(A)), order-4 table h = 1/20, start q = 0: max |OFLI_num - OFLI_exact| for t <= 10 P (= $imax10)", e, "<= delta = 1e-4 (Amendment D119)", pf(e <= 1.0e-4))
                end
            end
            for h in hs[1:(end - 1)]
                h2 = h / 2; (haskey(res, (key, h)) && haskey(res, (key, h2))) || continue
                a, b = res[(key, h)], res[(key, h2)]
                emit("4c", "PLATEAUp: pendulum amplitude $(Int(A)), order-$p: t*(h/2) / t*(h = 1/$(round(Int, 1 / h))) against 2^$p = $(2^p) (R1, prediction O11 ii: t* ~ h^-p)", (isfinite(a) && isfinite(b)) ? string(fmt(b / a), " (t* ", fmt(a), ", ", fmt(b), ")") : string("not both below Tmax (", fmt(a), ", ", fmt(b), ")"), "report, within a factor 2 of 2^p", "INFO")
            end
        end
    end
    serialize(joinpath(SCR, "plateau-pend.jls"), data)
end
section("PLATEAUf") do
    Tmax = 200_000; flow = HenonHeilesFlow(); po = periodic_orbit(); z0 = po.z0; W = of_frame(VSEED + 7, 1, 4, 4, Float64); data = Dict{Any, Any}()
    say("# PLATEAUf REDUCTION: reference horizon Tmax = $Tmax time units = $(round(Int, Tmax / po.T)) periods (the spec checks 100 T; t* above Tmax is reported as > Tmax); reference = Double64 order-10 table at h = 1/20 with tangent, seeded frame, k = 4; closure of the Float64 start $(fmt(po.closure))")
    D = Double64; m = SofroniouSpaletta10(D); h = D(1) / 20; z = SVector{4, D}(map(D, z0)); Wd = SMatrix{4, 4, D}(map(D, W)); t0 = time()
    Gref = zeros(Tmax + 1); G = -Inf
    for t in 0:Tmax
        f = vector_field(HenonHeilesFlow(), z); ff = dot(f, f)
        for i in 1:4
            w = Wd[:, i]; G = max(G, Float64(log(norm(w - (dot(f, w) / ff) * f))))
        end
        Gref[t + 1] = G
        t < Tmax && for _ in 1:20
            z, Wd = tangent_advance(m, flow, z, Wd, h)
        end
    end
    say(@sprintf("# PLATEAUf: reference done, %.1f s; reference OFLI_B(10 T) = %.6f, OFLI_B(100 T) = %.6f, OFLI_B(Tmax) = %.6f", time() - t0, Gref[round(Int, 10po.T) + 1], Gref[round(Int, 100po.T) + 1], Gref[end]))
    data[(:f3, :ref)] = decimate(collect(0.0:Tmax), Gref)
    g0 = ofli_initial(flow, z0, W); res = Dict{Any, Float64}()
    for (key, hs) in ((:o4, (1 / 4, 1 / 8, 1 / 16, 1 / 20)), (:o8, (1 / 2, 1 / 4, 1 / 20)))
        p = key === :o4 ? 4 : 8
        for hh in hs
            tau = round(Int, 1 / hh); im = IntegratorMap(of_table(Float64, key), flow, hh)
            snaps = of_snaps(OFLI(flow), im, z0, W, Tmax * tau, tau)
            if !all(s -> isfinite(s.ofli), snaps)
                emit("5c", "PLATEAUf: order-$p table, h = 1/$tau: non-finite series", 1, "report", "INFO"); continue
            end
            num = vcat(maximum(g0), [max(maximum(g0), s.ofli) for s in snaps]); dep = abs.(num .- Gref); ts = tstar(dep); res[(key, hh)] = ts
            data[(:ts, :f3, key, hh)] = ts
            data[(:f3, key, hh)] = decimate(collect(0.0:Tmax), num)
            n10, n100 = round(Int, 10po.T), round(Int, 100po.T)
            emit("5c", "PLATEAUf: F3 period-1 orbit, order-$p table, h = 1/$tau: t* (departure of OFLI_B from the reference above 0.05); max departure for t <= 100 T; OFLI_B(100 T) - OFLI_B(10 T) numerical; FLI_B(100 T) - FLI_B(10 T) numerical",
                 string(ts == Inf ? "> $Tmax" : fmt(ts), " ; ", fmt(maximum(dep[1:(n100 + 1)])), " ; ", fmt(num[n100 + 1] - num[n10 + 1]), " ; ", fmt(snaps[n100].fli - snaps[n10].fli)), "report (R1); dG <= 0.1, dFLI >= 1.5 for order 8", "INFO")
            if key === :o8 && hh == 1 / 20
                dG = num[n100 + 1] - num[n10 + 1]; dF = snaps[n100].fli - snaps[n10].fli
                emit("5b", "PLATEAUf: order-8 table h = 1/20 seeded frame: OFLI_B(100 T) - OFLI_B(10 T)", dG, "<= 0.1", pf(dG <= 0.1))
                emit("5b", "PLATEAUf: order-8 table h = 1/20 seeded frame: FLI_B(100 T) - FLI_B(10 T)", dF, ">= 1.5", pf(dF >= 1.5))
            end
        end
        for hh in hs[1:(end - 1)]
            hh == 1 / 16 && continue
            h2 = hh / 2; (haskey(res, (key, hh)) && haskey(res, (key, h2))) || continue
            a, b = res[(key, hh)], res[(key, h2)]
            emit("5c", "PLATEAUf: order-$p: t*(h/2) / t*(h = 1/$(round(Int, 1 / hh))) against 2^$p = $(2^p) (R1)", (isfinite(a) && isfinite(b)) ? string(fmt(b / a), " (t* ", fmt(a), ", ", fmt(b), ")") : string("not both below Tmax (", fmt(a), ", ", fmt(b), ")"), "report, within a factor 2 of 2^p", "INFO")
        end
    end
    serialize(joinpath(SCR, "plateau-f3.jls"), data)
end

# ------------------------------------------------------------------------------------------------ Od: ChaosTools sums cross-check (INFO)
section("Od") do
    flow = HenonHeilesFlow(); m = KahanLi(Float64, 8)
    f(u, p, t) = SVector(u[3], u[4], -(u[1] + 2 * u[1] * u[2]), -(u[2] + u[1]^2 - u[2]^2))
    jac(u, p, t) = SMatrix{4, 4, Float64}(0, 0, -(1 + 2 * u[2]), -2 * u[1], 0, 0, -2 * u[1], -(1 - 2 * u[2]), 1, 0, 0, 0, 0, 1, 0, 0)
    say("# Od: ChaosTools $(pkgversion(CT)), CoupledODEs Tsit5 abstol = reltol = 1e-14, hand-coded Jacobian, one renormalization per time unit (stop_at_Δt); loose sanity only (the ChaosTools spectrum is Gram-Schmidt, the OFLI sum is one unnormalized column), T = 100")
    for (name, z0) in (("R", of_orbit_R()), ("C", of_orbit_C()))
        u0 = SVector{4, Float64}(z0)
        ode = CT.CoupledODEs(f, u0, nothing; diffeq = (abstol = 1.0e-14, reltol = 1.0e-14)); tds = CT.TangentDynamicalSystem(ode; J = jac, k = 4)
        CT.reinit!(tds, u0; Q0 = one(SMatrix{4, 4, Float64})); lam = zeros(4); T = 100
        for _ in 1:T
            CT.step!(tds, 1.0, true); Q, Rm = qr(CT.current_deviations(tds))
            for j in 1:4
                lam[j] += log(abs(Rm[j, j]))
            end
            CT.set_Q_as_deviations!(tds, Q)
        end
        lam ./= T
        W0 = SMatrix{4, 4, Float64}(I); im = of_im(Float64); fin = propagate(im, [z0], 20T; k = 4, tau = 20, seed = 1, reducers = (OFLI(flow),), W0 = [W0]).finals[1][1]
        s1 = (fin.sum[1] + fin.comp[1]) / T; sB = maximum(fin.sum .+ fin.comp) / T
        emit("7", "Od: orbit $name, T = 100: ChaosTools lambda_1 (Gram-Schmidt); OFLI-reducer column sum S_1 / T (identity frame, one column); max_i S_i / T; OFLI_B / T", string(fmt(lam[1]), " ; ", fmt(s1), " ; ", fmt(sB), " ; ", fmt(fin.ofli / T)), "info (loose: window length 100, O(1 / T) end terms)", "INFO")
    end
end

# ------------------------------------------------------------------------------------------------ R lines (no threshold)
quant(v) = string(fmt(minimum(v)), " ; ", fmt(quantile(v, 0.25)), " ; ", fmt(median(v)), " ; ", fmt(quantile(v, 0.75)), " ; ", fmt(maximum(v)))
section("R2") do
    flow = HenonHeilesFlow(); im = of_im(Float64); cl = of_classes(OF_NCLASS, 1); h = OF_H
    say("# R2: the $(OF_NCLASS) starts per class of the package (random_energy_start, E = 1/8, classed by SALI at 2000 below 1e-8), k = 4, seeded frames; FLI_B - OFLI_B with the s = 0 term in OFLI_B; quantiles reported as min ; q25 ; median ; q75 ; max")
    for (cls, xs) in (("regular", cl.regular), ("chaotic", cl.chaotic)), T in (1000, 10_000)
        Ws = [of_frame(VSEED + 2, j, 4, 4, Float64) for j in eachindex(xs)]
        res = propagate(im, xs, 20T; k = 4, tau = 20, seed = 1, reducers = (OFLI(flow),), W0 = Ws)
        d = [res.finals[j][1].fli - max(maximum(ofli_initial(flow, xs[j], Ws[j])), res.finals[j][1].ofli) for j in eachindex(xs)]
        G = [max(maximum(ofli_initial(flow, xs[j], Ws[j])), res.finals[j][1].ofli) for j in eachindex(xs)]
        emit("R2", "$cls class ($(length(xs)) starts), t = $T: FLI_B - OFLI_B quantiles; negative count (must be 0 up to rounding)", string(quant(d), " ; neg ", count(<(-1.0e-12), d)), "no threshold", "INFO")
        emit("R2", "$cls class, t = $T: OFLI_B quantiles ; FLI_B quantiles", string(quant(G), " || ", quant(G .+ d)), "no threshold", "INFO")
    end
    zC = of_orbit_C(); _, _, st = propagate_orbit(im, zC, of_frame(VSEED + 3, 1, 4, 1, Float64), 20_000, 20, (Benettin(),)); lam = reducer_final(Benettin(), st[1]).lambda[1] / h
    W = of_frame(VSEED + 3, 1, 4, 4, Float64); fin = propagate(im, [zC], 20_000; k = 4, tau = 20, seed = 1, reducers = (OFLI(flow),), W0 = [W]).finals[1][1]
    GB = max(maximum(ofli_initial(flow, zC, W)), fin.ofli)
    emit("R2", "loose sanity, orbit C, t = 1000: OFLI_B / t ; L1 finite-time exponent of the identical orbit ; FLI_B - OFLI_B", string(fmt(GB / 1000), " ; ", fmt(lam), " ; ", fmt(fin.fli - GB)), "no threshold", "INFO")
end
section("R4") do
    flow = HenonHeilesFlow(); po = periodic_orbit(); z0 = po.z0; W = SMatrix{4, 4, Float64}(I); n10, n100 = round(Int, 10po.T), round(Int, 100po.T); g0 = maximum(ofli_initial(flow, z0, W))
    for (key, h) in ((:o8, 1 / 20), (:o4, 1 / 20))
        snaps = of_snaps(OFLI(flow), IntegratorMap(of_table(Float64, key), flow, h), z0, W, n100 * 20, 20)
        G = [max(g0, s.ofli) for s in snaps]; inst = [s.inst for s in snaps]
        win = (n10 + 1):n100
        nrise = count(i -> G[i] > G[i - 1], win)
        emit("R4", "periodic orbit of F3, order-$(key === :o8 ? 8 : 4) table, canonical frame, t in [10 T, 100 T] = [$n10, $n100]: running supremum OFLI_B at 10 T, 100 T, and increase", string(fmt(G[n10]), " ; ", fmt(G[n100]), " ; ", fmt(G[n100] - G[n10])), "report (F2)", "INFO")
        emit("R4", "same: instantaneous value max_i g_i(n): min ; max ; max - min over the window, handoffs where the supremum rose, handoffs with inst < sup - 0.01 (oscillates below a converging supremum if the last two are large and small)", string(fmt(minimum(inst[win])), " ; ", fmt(maximum(inst[win])), " ; ", fmt(maximum(inst[win]) - minimum(inst[win])), " ; ", nrise, " ; ", count(i -> inst[i] < G[i] - 0.01, win), " of ", length(win)), "report (F2)", "INFO")
    end
end
section("R5") do
    flow = HenonHeilesFlow(); im = of_im(Float64); xs = of_classes(OF_NCLASS, 1).chaotic; cap = log(1.0e20); Tm = 10_000
    say("# R5: $(length(xs)) chaotic-class starts, seeded frames, horizon t = $Tm (the sticky starts that do not reach ln 1e20 = $(fmt(cap)) by then are counted); first handoff with OFLI_B >= ln 1e20")
    tr = fill(NaN, length(xs))
    Threads.@threads for j in eachindex(xs)
        W = of_frame(VSEED + 5, j, 4, 4, Float64); g0 = maximum(ofli_initial(flow, xs[j], W))
        snaps = of_snaps(OFLI(flow), im, xs[j], W, 20Tm, 20); k = findfirst(s -> s.ofli >= cap, snaps)
        g0 >= cap && (k = 0)
        k === nothing || (tr[j] = Float64(k))
    end
    ok = filter(!isnan, tr)
    emit("R5", "chaotic class: starts reaching OFLI_B = ln 1e20 = $(fmt(cap)) within t = $Tm: count of $(length(xs)); time quantiles min ; q25 ; median ; q75 ; max", string(length(ok), " ; ", isempty(ok) ? "-" : quant(ok)), "no threshold", "INFO")
end

# ------------------------------------------------------------------------------------------------ AC16T: 1 thread against N threads, bitwise
function ac16_data()
    flow = HenonHeilesFlow(); im = of_im(Float64); n = 4096; a, b = section_interval(flow, 1 / 6); w = sqrt(2 / 6); rng = StableRNG(VSEED + 16)
    xs = [random_energy_start(flow, 1 / 6, (a, b), (-w, w), rng)[1] for _ in 1:n]
    return [(tau = tau, res = (res = propagate(im, xs, 20_000; k = 4, tau = tau, seed = VSEED + 16, reducers = (OFLI(flow),)); (finals = res.finals, x = res.x, W = res.W))) for tau in (1, 20)]
end
section("AC16child") do
    serialize(get(ENV, "OF_CHILD_OUT", joinpath(SCR, "child.jls")), ac16_data())
end
section("AC16T") do
    child = joinpath(SCR, "child-t1.jls")
    cmd = addenv(`$(Base.julia_cmd()) -t 1 --project=$(joinpath(ROOT, "verify")) $(@__FILE__)`, "OF_ONLY" => "AC16child", "OF_CHILD_OUT" => child, "OF_LOG" => joinpath(SCR, "child-t1.log"), "OF_SCRATCH" => SCR)
    t1 = @elapsed mine = ac16_data()
    t2 = @elapsed run(cmd)
    theirs = deserialize(child)
    say(@sprintf("# AC16T: 4096 starts of F3 at E = 1/6 (random_energy_start, as AC16a), t = 1e3 (20000 steps), tau in {1, 20}: %d threads %.1f s, 1 thread (child) %.1f s; compared: the final OFLI record (all fields: maxima, sums, compensation terms, maximisers, diagnostics), final state, final frame", Threads.nthreads(), t1, t2))
    for (a, b) in zip(mine, theirs)
        emit("16", "AC16T: tau = $(a.tau): OFLI final records bitwise identical, $(Threads.nthreads()) against 1 CPU threads (4096 starts)", isequal(a.res.finals, b.res.finals), "bitwise", pf(isequal(a.res.finals, b.res.finals)))
        emit("16", "AC16T: tau = $(a.tau): final states and frames bitwise identical", isequal(a.res.x, b.res.x) && isequal(a.res.W, b.res.W), "bitwise", pf(isequal(a.res.x, b.res.x) && isequal(a.res.W, b.res.W)))
    end
    emit("16", "AC16T: GPU against CPU", "not run", "-", "NOT RUN: needs verify/gpu.jl (writes data/); owner decision")
end

# AC18r: the CPU part of criterion 18 with a reduced start count (24 per class instead of 1024), from a patched scratch copy of the AC18 testset (the test file itself is not changed); the full-size FULL4 run is the primary
section("AC18r") do
    txt = read(joinpath(TD, "of_ac16_18.jl"), String); i = findfirst("@testset \"AC18", txt)
    body = replace(txt[first(i):end], "m = min(OF_NCLASS, of_size(2, 1024))" => "m = 24")
    f = joinpath(SCR, "ac18-reduced.jl"); write(f, body)
    say("# AC18r REDUCTION: 24 starts per class (spec: 1024), t = 400 (order-8 table), CPU backends Float64 and Double64 against their own BigFloat(256) references; patched copy of the AC18 testset in the scratch folder")
    Test.@testset "outer AC18r" begin
        ts = Test.@testset "AC18r" begin
            include(f)
        end
        report_ts("AC18r", ts)
    end
end

# ------------------------------------------------------------------------------------------------ FIG
section("FIG") do
    pend = deserialize(joinpath(SCR, "plateau-pend.jls")); f3 = deserialize(joinpath(SCR, "plateau-f3.jls"))
    fig = CM.Figure(size = (1500, 900)); cols = Dict((:o4, 1 / 4) => :red, (:o4, 1 / 8) => :orange, (:o8, 1 / 2) => :blue, (:o8, 1 / 4) => :teal)
    names = Dict((:o4, 1 / 4) => "order 4, h = 1/4", (:o4, 1 / 8) => "order 4, h = 1/8", (:o8, 1 / 2) => "order 8, h = 1/2", (:o8, 1 / 4) => "order 8, h = 1/4")
    panels = ((pend, (:pend, 2.0), "pendulum amplitude 2, start q = 0, OFLI_1"), (pend, (:pend, 3.0), "pendulum amplitude 3, start q = 0, OFLI_1"), (f3, (:f3,), "F3 period-1 orbit, E = 1/8, OFLI_B (k = 4)"))
    for (j, (D, pre, title)) in enumerate(panels)
        ax1 = CM.Axis(fig[1, j]; xscale = log10, xlabel = "t", ylabel = "OFLI (log form)", title = title, titlesize = 13)
        ax2 = CM.Axis(fig[2, j]; xscale = log10, yscale = log10, xlabel = "t", ylabel = "|OFLI_num - OFLI_ref|")
        ref = D[(pre..., :ref)]; m = ref.t .> 0
        CM.lines!(ax1, ref.t[m], ref.g[m]; color = :black, linestyle = :dash, label = "reference (exact flow)")
        for (key, h) in ((:o4, 1 / 4), (:o4, 1 / 8), (:o8, 1 / 2), (:o8, 1 / 4))
            haskey(D, (pre..., key, h)) || continue
            s = D[(pre..., key, h)]; col = cols[(key, h)]
            CM.lines!(ax1, s.t[m], s.g[m]; color = col, label = names[(key, h)])
            dep = max.(abs.(s.g .- ref.g), 1.0e-18); CM.lines!(ax2, s.t[m], dep[m]; color = col)
            ts = D[(:ts, pre..., key, h)]
            isfinite(ts) && (CM.vlines!(ax1, [max(ts, 1.0)]; color = col, linestyle = :dot); CM.vlines!(ax2, [max(ts, 1.0)]; color = col, linestyle = :dot))
        end
        CM.hlines!(ax2, [0.05]; color = :gray, linestyle = :dashdot)
        j == 1 && CM.axislegend(ax1; position = :lt, labelsize = 11)
    end
    CM.Label(fig[0, 1:3], "OFLI plateau and its loss on periodic orbits (dotted verticals: t*, departure above 0.05; dash-dot: 0.05)"; fontsize = 15)
    CM.save(joinpath(OUT, "ofli-plateau.png"), fig)
    # classes: OFLI_B and FLI_B against t for R, C and the periodic orbit, order-8 table, h = 1/20, k = 4, canonical frame
    flow = HenonHeilesFlow(); po = periodic_orbit(); W = SMatrix{4, 4, Float64}(I); T = 2000; im = of_im(Float64)
    fig2 = CM.Figure(size = (900, 650)); ax = CM.Axis(fig2[1, 1]; xscale = log10, xlabel = "t", ylabel = "OFLI_B, FLI_B (log form)", title = "OFLI_B (solid) and FLI_B (dashed), E = 1/8, order-8 table, h = 1/20")
    for (nm, z0, col) in (("orbit R (quasi-periodic)", of_orbit_R(), :blue), ("orbit C (chaotic)", of_orbit_C(), :red), ("period-1 orbit", po.z0, :green))
        snaps = of_snaps(OFLI(flow), im, z0, W, 20T, 20); g0 = maximum(ofli_initial(flow, z0, W)); t = 1:T
        CM.lines!(ax, t, [max(g0, snaps[i].ofli) for i in t]; color = col, label = nm)
        CM.lines!(ax, t, [snaps[i].fli for i in t]; color = col, linestyle = :dash)
    end
    CM.axislegend(ax; position = :lt)
    CM.save(joinpath(OUT, "ofli-classes.png"), fig2)
    say("# FIG: wrote ofli-plateau.png and ofli-classes.png in the stage 04 output folder")
end
