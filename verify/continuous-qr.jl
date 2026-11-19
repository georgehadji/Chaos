# Stage 04 verification of module continuous-qr (catalog row L8: Chaoslab.Lyapunov continuous_qr.jl; spec D120).
# Run from the repo root, one process per group (CQ_ONLY, comma list; CQ_LOG names the log, CQ_APPEND=1 appends):
#   CQ_ONLY=FULL1 .. FULL4  julia --project=verify -t 8 verify/continuous-qr.jl   (the package test files of test/continuous-qr/ at the FULL sizes, CQ_FULL_RUN = true; one file each)
#   CQ_ONLY=Oa      (BigFloat 256-bit closed forms of the linear flows: criteria 1 to 5, exact and degenerate exponents, order slopes)
#   CQ_ONLY=Oid     (flow-direction and energy-covector identities on a TaylorIntegration BigFloat(256) reference: criteria 7, 8)
#   CQ_ONLY=Od20    (TaylorIntegration BigFloat(256) variational reference of F3 orbits R, C and the pendulum: criterion 20, R2, R4, convergence data)
#   CQ_ONLY=Od20R   (orbit R against the Taylor reference to t = 1000, reported)
#   CQ_ONLY=Oc9     (same-orbit comparison with the discrete QR / Benettin at T = 1000, criterion 9, R6)
#   CQ_ONLY=Oord    (discretization order of the continuous against the discrete spectrum, slopes at h and h/2)
#   CQ_ONLY=Od21    (ChaosTools Lyapunov spectrum cross-check, criterion 21)
#   CQ_ONLY=AC17T   (1 against N CPU threads, bitwise; spawns a one-thread child), AC17alloc (allocation at t = 1e3)
#   CQ_ONLY=AC22    (CPU fallback record; the device part is NOT RUN), R1 (wall-clock ratios), FIG (the two figures)
# Not run: criterion 22 device part (verify/gpu.jl writes data/gpu-benchmark.toml; owner decision). Nothing here writes under data/.
# Spec: workspaces/module-build/stages/01-spec/output/continuous-qr-spec.md; test plan points U1 to U13 and E_max = 1e-2 are the operative reading. Verify-only seeds 20261012 + 44900 ... (package seeds 44000 + 50 ac + sub).
# One printed line per measurement:  AC<n> | quantity | value | bound | PASS/FAIL/INFO
const ONLY = filter(!isempty, split(get(ENV, "CQ_ONLY", ""), ","))
using Printf, Random, LinearAlgebra, Statistics, Serialization, Logging
using Test
const CQ_FULL_RUN = true
using Chaoslab, StaticArrays, DoubleFloats, StableRNGs, TOML
want(name) = name in ONLY || (isempty(ONLY) && !startswith(name, "AC17child"))
const NEEDTAY = isempty(ONLY) || any(o -> o in ("Oid", "Od20", "Od20R"), ONLY)
NEEDTAY && (import TaylorIntegration as TI)
(isempty(ONLY) || "Od21" in ONLY) && (import ChaosTools as CT)
(isempty(ONLY) || "FIG" in ONLY) && (import CairoMakie as CM)

const ROOT = dirname(@__DIR__)
const OUT = joinpath(ROOT, "workspaces", "module-build", "stages", "04-verify", "output")
const SCR = get(ENV, "CQ_SCRATCH", joinpath(tempdir(), "cq-verify-scratch"))
mkpath(SCR)
const LOG = open(get(ENV, "CQ_LOG", joinpath(OUT, "continuous-qr-verify-log.txt")), get(ENV, "CQ_APPEND", "") == "1" ? "a" : "w")

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

say("# continuous-qr verification, julia $(VERSION), threads $(Threads.nthreads()), base seed 20261012 + 44000, CQ_FULL_RUN = true, groups $(isempty(ONLY) ? "all" : join(ONLY, ",")), time $(round(time())) (unix)")
say("# timings in this log were taken on a machine that may be loaded by other processes")

const TD = joinpath(ROOT, "test", "continuous-qr")
include(joinpath(TD, "cq_helpers.jl"))
# Definitions of the test files without their testsets (the files depend on each other's helpers: cq_refl, cq_AFmax, cq_ref_run, cq_hyp_error ...).
function load_defs(file)
    for ex in Meta.parseall(read(joinpath(TD, file), String)).args
        ex isa Expr || continue
        (ex.head === :macrocall && ex.args[1] === Symbol("@testset")) && continue
        Core.eval(Main, ex)
    end
end
for f in ("cq_ac01_05.jl", "cq_ac06_10.jl", "cq_ac11_15.jl", "cq_ac16_22.jl")
    load_defs(f)
end
const VSEED = 20261012 + 44900
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
function npass(ts)
    n = ts.n_passed
    for r in ts.results
        r isa Test.AbstractTestSet && (n += npass(r))
    end
    return n
end
function nbroken(ts)
    n = 0
    for r in ts.results
        r isa Test.AbstractTestSet ? (n += nbroken(r)) : (r isa Test.Broken && (n += 1))
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
        emit(ac, string("testset ", ch.description, " (", npass(ch), " pass, ", nf, " fail or error, ", nbroken(ch), " broken, ", @sprintf("%.1f", ch.time_end - ch.time_start), " s)"), nf, "0", pf(nf == 0))
    end
end
# `patches` replace text in a scratch copy of the file (the test file itself is never changed); `note` is the INFO line naming the size reduction.
function suite_section(name, file; patches = Pair{String, String}[], note = "")
    section(name) do
        path = joinpath(TD, file)
        if !isempty(patches)
            txt = read(path, String)
            for (a, b) in patches
                occursin(a, txt) || error("patch target not found: $a")
                txt = replace(txt, a => b)
            end
            path = joinpath(SCR, name * "-" * file); write(path, txt)
            say("# $name REDUCTION: $note (patched copy in the scratch folder)")
        end
        try
            Test.@testset "outer $name" begin
                ts = Test.@testset "$name" begin
                    include(path)
                end
                report_ts(name, ts)
            end
        catch e
            e isa Test.TestSetException || rethrow()
        end
    end
end
suite_section("FULL1", "cq_ac01_05.jl")      # criteria 1 to 5
suite_section("FULL2", "cq_ac06_10.jl")      # criteria 6 to 10
suite_section("FULL3", "cq_ac11_15.jl")      # criteria 11 to 15
suite_section("FULL4", "cq_ac16_22.jl")      # criteria 16 to 22 (CPU part), R1, R2, R5
# Diagnostic re-runs of FULL2 at the horizons the stage 02 calibration used (AC6b t = 60, AC8 Tmax = 100), to separate the horizon from the sample size
suite_section("FULL2h", "cq_ac06_10.jl"; patches = ["t = cq_size(60.0, 100.0); n = round" => "t = 60.0; n = round", "Tmax = cq_size(100, 300); n = Tmax * 20" => "Tmax = 100; n = Tmax * 20"],
              note = "DIAGNOSTIC, not a reduction: AC6b at t = 60 (calibration horizon, spec 100) and AC8 at Tmax = 100 (calibration horizon, spec 300); every other size stays full")

# ---------------------------------------------------------------------------------------------------
# Taylor helpers (hand-written fields, independent of Chaoslab's): F3 state + 4 variational columns, pendulum and oscillator with k columns
function hh_field(z)
    return [z[3], z[4], -(z[1] + 2 * z[1] * z[2]), -(z[2] + z[1]^2 - z[2]^2)]
end
hh_grad(z) = [z[1] + 2 * z[1] * z[2], z[2] + z[1]^2 - z[2]^2, z[3], z[4]]
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
function pend_var_f!(dx, x, p, s)               # pendulum H = p^2/2 - cos q, two variational columns
    dx[1] = x[2]; dx[2] = -sin(x[1]); c = cos(x[1])
    for j in 1:2
        o = 2 * j
        dx[o + 1] = x[o + 2]; dx[o + 2] = -c * x[o + 1]
    end
    return nothing
end
tay(f!, x0, T; order = 40, maxsteps = 10^6, kw...) = TI.taylorinteg(f!, x0, zero(BigFloat), BigFloat(T), order, ldexp(BigFloat(1), -200); maxsteps = maxsteps, kw...)
function tay_eval(sol, t)
    i = clamp(searchsortedlast(sol.t, t), 1, length(sol.t) - 1)
    s = t - sol.t[i]
    return [sol.p[i, k](s) for k in 1:size(sol.x, 2)]
end
# l_i of the evolved frame W (BigFloat matrix) against the stored frame X in 256 bits (C16), kept in BigFloat
function lsplit_big(W, X)
    k = size(W, 2); out = zeros(BigFloat, k); prev = zero(BigFloat)
    for m in 1:k
        lv = (log(det(W[:, 1:m]' * W[:, 1:m])) - log(det(X[:, 1:m]' * X[:, 1:m]))) / 2
        out[m] = lv - prev; prev = lv
    end
    return out
end
# the Taylor reference of one case: state z0 (Float64), frame Q0 (Float64 d x k), horizon T; returns a function t -> (z(t), W(t)) in BigFloat
function tay_case(kind::Symbol, z0, Q0, T)
    d = length(z0); k = size(Q0, 2)
    x0 = vcat(BigFloat.(collect(z0)), BigFloat.(vec(Matrix(Q0))))
    sol = kind === :hh ? tay(var_f!, x0, T) : tay(pend_var_f!, x0, T; parse_eqs = false)
    return t -> (x = tay_eval(sol, BigFloat(t)); (z = x[1:d], W = reshape(x[(d + 1):(d + d * k)], d, k)))
end
const HH_ORBITS = (("R", cq_orbit_R()), ("C", cq_orbit_C()))

# ------------------------------------------------------------------------------------------------ Oa: BigFloat(256) closed forms (criteria 1 to 5)
section("Oa") do
    T = BigFloat; eT = cq_eps(T); Tf = Float64
    say("# Oa: element type BigFloat (256 bits), the closed forms of C19 on constant matrices (no orbit); references: exact matrix exponential and Gram determinants in 256 bits (C16)")
    # criterion 1: fixed point, n = 2000 at h = 0.05, both schemes
    n = 2000; h = 0.05
    for (name, A) in (("diag(0.7, -0.7)", cq_diag2()), ("diag(0.7, 0.3, -0.3, -0.7)", cq_diag4()), ("triangular [[0.5, 2], [0, -0.5]]", cq_tri2()), ("direct sum d = 4", cq_tri4())), q in (2, 4)
        d = size(A, 1); cq = cq_est(A, T, h, q)
        _, Q, fin, ser = propagate_cqr_series(cq, cq_x0(d, T), cq_eye(d, d, T), n; every = 1)
        qok = all(s -> s.Q == cq_eye(d, d, T), ser) && Q == cq_eye(d, d, T)
        w = 0.0
        for s in ser, i in 1:d
            a = T(A[i, i]); t = s.n * T(h)
            w = max(w, Float64(abs(s.l[i] - a * t)) / (16 * eT * max(1.0, Float64(abs(a * t)))))
        end
        emit("1", "Oa: $name, q = $q, n = $n: Q = I at every Q-step (bitwise); E_max seen = $(fmt(Float64(fin.Emax)))", qok && fin.Emax == 0, "true", pf(qok && fin.Emax == 0))
        emit("1", "Oa: $name, q = $q: worst |l_i - A_ii t| / (16 eps max(1, |A_ii| t)) over every Q-step", w, "<= 1", pf(w <= 1))
    end
    # criterion 2: skew A = [[0, I], [-I, 0]], Q(t) = exp(A t) Q0 = [[cos t, sin t], [-sin t, cos t]] Q0 in 256 bits, t up to 100
    h = 0.05; n = 2000
    for N in (1, 2), q in (2, 4)
        d = 2N; A = cq_skew(N); cq = cq_est(A, T, h, q); rng = StableRNG(VSEED + N); AF = cq_AF(A)
        frames = vcat([cq_eye(d, d, T)], [cq_frame(rng, d, d, T) for _ in 1:3])
        wl = 0.0; wq = 0.0; wEc = 0.0; wl_over = 0.0
        for Q0 in frames
            every = max(1, n ÷ (q ÷ 2) ÷ 100)
            _, Q, fin, ser = propagate_cqr_series(cq, cq_x0(d, T), Q0, n; every = every)
            wEc = max(wEc, Float64(fin.Emax))
            for s in ser
                t = s.n * T(h); c = cos(t); sn = sin(t); Id = Matrix{T}(I, N, N)
                Y = [c * Id sn * Id; -sn * Id c * Id]
                wl = max(wl, maximum(i -> Float64(abs(s.l[i])) / (8 * eT * Float64(t)), 1:d))
                B = min(cq_Blin(Float64(t), h, q, AF, T), 0.25)
                wq = max(wq, Float64(maximum(abs, Matrix(s.Q) - Y * Matrix(Q0))) / B)
            end
        end
        emit("2", "Oa: skew F1 omega = 1, N = $N (d = $d), q = $q, 4 frames, t to 100: worst |l_i| / (8 eps t)", wl, "<= 1", pf(wl <= 1))
        emit("2", "Oa: same: worst max|Q(t) - exp(A t) Q0| / min(B_lin, 0.25) (cap of plan U11); E_max seen $(fmt(wEc))", wq, "<= 1 (E_max <= $(CQ_EMAX))", pf(wq <= 1 && wEc <= CQ_EMAX))
    end
    # criterion 3: hyperbolic closed forms, individually for the degenerate pair
    rng = StableRNG(VSEED + 3)
    for (name, A) in cq_hyp_fixtures(), q in (2, 4)
        d = size(A, 1); AF = cq_AF(A); worst = 0.0; wE = 0.0; spread = 0.0
        frames = vcat([cq_eye(d, d, T)], [cq_frame(rng, d, d, T) for _ in 1:3]); l2l1 = Float64[]
        for Q0 in frames, h in (0.05, 0.025), t in (2.0, 10.0)
            r = cq_hyp_error(A, T, Q0, h, q, t); B = cq_Blin(t, h, q, AF, T)
            worst = max(worst, maximum(r.err) / B); wE = max(wE, Float64(r.fin.Emax))
            h == 0.025 && t == 10.0 && push!(l2l1, d == 4 ? r.ref[2] - r.ref[1] : 0.0)
        end
        what = name == "a" ? "inverted oscillator lam = 0.8 (d = 2)" : name == "b" ? "coupled hyperbolic lam = (0.8, 0.35), d = 4" : "degenerate lam = (0.5, 0.5), d = 4 (every l_i individually, not their sum)"
        emit("3", "Oa: $what, q = $q, 4 frames x 2 h x 2 t: worst max_i |l_i - l_i,ref| / B_lin (E_max seen $(fmt(wE)))", worst, "<= 1 (E_max <= $(CQ_EMAX))", pf(worst <= 1 && wE <= CQ_EMAX))
        name == "c" && q == 4 && emit("3", "Oa: degenerate case, l_2,ref - l_1,ref over the 4 frames at t = 10 (frame dependence of the individual values; the reference is the same Gram determinant)", string(map(fmt, l2l1)), "report", "INFO")
    end
    # criterion 4: order of the Q-scheme in BigFloat (no rounding floor): ratios and fitted slopes
    t = 10.0; hs = (0.1, 0.05, 0.025, 0.0125, 0.00625); conv = Dict{Any, Any}(); rng = StableRNG(VSEED + 4)
    for (name, A) in cq_hyp_fixtures()[1:2], q in (2, 4)
        d = size(A, 1); Q0 = cq_frame(rng, d, d, T)
        errs = [maximum(cq_hyp_error(A, T, Q0, h, q, t).err) for h in hs]
        ratios = errs[1:(end - 1)] ./ errs[2:end]; sl = cq_slope(log.(collect(hs)), log.(errs))
        # D121 (C4): the ratio gate holds for the pairs whose coarse step has Delta ||A||_F <= 0.35 (Delta = q h / 2); the coarser pairs are pre-asymptotic and reported
        gated = [q / 2 * hs[j] * cq_AF(A) <= 0.35 * (1 + 1.0e-9) for j in 1:(length(hs) - 1)]; inr(r) = 2^q / 4 <= r <= 2 * 2^q; ok = all(inr, ratios[gated])
        emit("4", "Oa: fixture $name, q = $q, t = 10, h = 0.1 .. 0.00625: error ratios of consecutive h (errors $(join(map(fmt, errs), ", "))); pairs with Delta ||A||_F <= 0.35 gated, Delta ||A||_F of the coarse steps $(join(map(fmt, [q / 2 * hs[j] * cq_AF(A) for j in 1:(length(hs) - 1)]), ", "))", join(map(fmt, ratios), ", "), "gated pairs ($(join(findall(gated), ", "))) each in [$(2^q / 4), $(2 * 2^q)]", pf(ok))
        any(.!gated) && emit("4", "Oa: fixture $name, q = $q: ratios of the pre-asymptotic pairs (Delta ||A||_F > 0.35), reported", join(map(fmt, ratios[.!gated]), ", "), "report", "INFO")
        emit("4", "Oa: fixture $name, q = $q: least-squares slope of ln error against ln h", sl, "report (expected $q)", "INFO")
        conv[(name, q)] = (collect(hs), errs, sl)
    end
    # the same ladder at h up to 0.4 for the figure (coupled fixture, stability limit nu Delta <= 1.39)
    hf = [0.2, 0.1, 0.05, 0.025, 0.0125, 0.00625, 0.003125]; A = cq_hyp_fixtures()[2][2]; d = 4; Q0 = cq_frame(StableRNG(VSEED + 44), d, d, T)
    figd = Dict{Int, Any}()
    for q in (2, 4)
        e = Float64[]
        for h in hf
            v = try
                maximum(cq_hyp_error(A, T, Q0, h, q, 10.0).err)
            catch
                NaN
            end
            push!(e, v)
        end
        figd[q] = e
    end
    serialize(joinpath(SCR, "conv-closed.jls"), (hs = hf, errs = figd))
    # criterion 5: rotated frame, closed form in 256 bits
    lam = T(0.8); A = [0.8 0.0; 0.0 -0.8]; h = 0.05; n = 400; AF = cq_AF(A)
    for theta in (0.3, 1.0, 1.4), q in (2, 4)
        cq = cq_est(A, T, h, q); th = T(theta); Q0 = cq_rot2(th, T)
        _, _, fin, ser = propagate_cqr_series(cq, cq_x0(2, T), Q0, n; every = max(1, n ÷ (q ÷ 2) ÷ 40))
        w1 = 0.0; w2 = 0.0; wQ = 0.0
        for s in ser
            t = s.n * T(h); B = min(cq_Blin(Float64(t), h, q, AF, T), 0.25)
            w1 = max(w1, Float64(abs(s.l[1] - cq_l1_rot(lam, th, t))) / B); w2 = max(w2, Float64(abs(s.l[1] + s.l[2])) / B)
            wQ = max(wQ, Float64(maximum(abs, s.Q - cq_rot2(cq_phi_rot(lam, th, t), T))) / (B / max(AF * Float64(t), 1)))
        end
        emit("5", "Oa: A = diag(0.8, -0.8), theta = $theta, q = $q, t <= 20: worst error / bound for l_1 ; l_1 + l_2 ; Q(t) against the rotation by phi(t)", string(fmt(w1), " ; ", fmt(w2), " ; ", fmt(wQ)), "<= 1 each", pf(max(w1, w2, wQ) <= 1))
    end
end

# ------------------------------------------------------------------------------------------------ Oid: the two exact identities on the Taylor BigFloat(256) reference (criteria 7, 8)
section("Oid") do
    say("# Oid: reference = Taylor order 40, abstol 2^-200, 256 bits (integrator-catalog O7, henon-heiles O3 settings), hand-written fields; frames in 256 bits")
    setprecision(BigFloat, BITS) do
        rng = StableRNG(VSEED + 7); ts = (10, 30, 60)
        # (a) exactness of the identities on the exact flow (no module): F3 orbits R, C; pendulum; oscillator
        for (name, z0) in HH_ORBITS
            z0b = BigFloat.(collect(z0)); f0 = hh_field(z0b)
            Qd = Matrix(cq_orth(hcat(f0 / norm(f0), BigFloat.(randn(rng, 4, 3)))))
            x0 = vcat(z0b, vec(Qd)); solb = tay(var_f!, x0, 60)
            Qe = cq_frame(rng, 4, 4, Float64); g0 = hh_grad(z0b)
            while abs(dot(g0, BigFloat.(Qe[:, 1]))) < 0.05 * norm(g0); Qe = cq_frame(rng, 4, 4, Float64); end
            sole = tay_case(:hh, z0, Qe, 60)
            for t in ts
                x = tay_eval(solb, BigFloat(t)); zt = x[1:4]; W = reshape(x[5:20], 4, 4)
                e17 = abs(log(norm(W[:, 1]) / norm(Qd[:, 1])) - log(norm(hh_field(zt)) / norm(f0)))
                r = sole(t); q10 = BigFloat.(Qe[:, 1]); q10 = q10 / norm(q10); Yq = r.W[:, 1] / norm(BigFloat.(Qe[:, 1])); l1 = log(norm(Yq)); qt = Yq / norm(Yq)
                e18 = abs(l1 - (log(abs(dot(g0, q10))) - log(abs(dot(hh_grad(r.z), qt)))))
                emit("7", "Oid: F3 orbit $name, t = $t, exact flow: |ln ||Y q_1|| - ln(||f(z_t)|| / ||f(z_0)||)| with q_1(0) = f(z_0)/||f(z_0)|| (256 bits)", Float64(e17), "<= 1e-30", pf(e17 <= BigFloat(10)^-30))
                emit("8", "Oid: F3 orbit $name, t = $t, exact flow, seeded frame: |l_1 - (ln|gradH(z_0).q_1(0)| - ln|gradH(z_t).q_1(t)|)| (256 bits)", Float64(e18), "<= 1e-30", pf(e18 <= BigFloat(10)^-30))
            end
        end
        # (b) module against the Taylor orbit: direction identity on F3, pendulum, F1; the same plan bounds as AC7 (U7)
        bdir(name, q) = q == 4 ? (name == "pendulum A=3" ? 1.0e-1 : 5.0e-3) : 1.0
        cases = (("pendulum A=2", Pendulum(1.0), cq_pendulum_start(2.0), :o4, :pend), ("pendulum A=3", Pendulum(1.0), cq_pendulum_start(3.0), :o4, :pend),
                 ("F1 w=1.3", HarmonicOscillator([1.3]), SVector(0.5, 0.3), :o4, :osc), ("F3 R", HenonHeilesFlow(), cq_orbit_R(), :o8, :hh), ("F3 C", HenonHeilesFlow(), cq_orbit_C(), :o8, :hh))
        for (name, flow, z0, key, kind) in cases
            z0b = BigFloat.(collect(z0)); d = length(z0)
            reft = if kind === :hh
                sol = tay(hh_f!, z0b, 60)
                t -> (x = tay_eval(sol, BigFloat(t)); (z = x, f = hh_field(x)))
            else
                om2 = kind === :osc ? BigFloat(1.3)^2 : zero(BigFloat)
                fv! = kind === :osc ? ((dx, x, p, s) -> (dx[1] = x[2]; dx[2] = -om2 * x[1]; nothing)) : ((dx, x, p, s) -> (dx[1] = x[2]; dx[2] = -sin(x[1]); nothing))
                sol = tay(fv!, z0b, 60; parse_eqs = false)
                t -> (x = tay_eval(sol, BigFloat(t)); (z = x, f = kind === :osc ? [x[2], -om2 * x[1]] : [x[2], -sin(x[1])]))
            end
            f0 = vector_field(flow, z0); k = min(d, 4)
            for q in (2, 4), h in (1 / 20, 1 / 40)
                (q == 2 && name == "pendulum A=3") && continue
                Q0 = cq_dir_frame(flow, z0, k, StableRNG(VSEED + 70)); cq = cq_cqr(cq_im(flow, Float64, key, h), flow; order = q)
                for t in ts
                    r = cq_final(cq, z0, Q0, round(Int, t / h)); rf = reft(t)
                    ref = Float64(log(norm(rf.f) / norm(BigFloat.(collect(f0))))); e = abs(r.fin.l[1] - ref)
                    emit("7", "Oid: $name, q = $q, h = 1/$(round(Int, 1 / h)), t = $t: |l_1 - ln(||f(z_t)|| / ||f(z_0)||)| with z_t from the Taylor reference", e, "<= $(bdir(name, q)) (plan U7; spec 1e-3 / 1e-1)", pf(e <= bdir(name, q)))
                end
            end
            name == "pendulum A=3" && emit("7", "Oid: pendulum A=3, order-2 scheme", "skipped", "plan U7 (realized error 3.4 at t = 60)", "INFO")
        end
        # (c) module against the Taylor orbit: energy covector identity on F3 R, C (grad H at the Taylor point, q_1(t) of the module), plan bound B_en = 3e-2 + 8 eps / c (U8)
        flow = HenonHeilesFlow(); im = cq_im(flow, Float64, :o8, 1 / 20); cq = cq_cqr(im, flow; order = 4)
        for (name, z0) in HH_ORBITS
            Qe = cq_frame(rng, 4, 4, Float64); g0 = hh_grad(BigFloat.(collect(z0)))
            while abs(dot(g0, BigFloat.(Qe[:, 1]))) < 0.05 * norm(g0); Qe = cq_frame(rng, 4, 4, Float64); end
            tc = tay_case(:hh, z0, Qe, 60)
            for t in ts
                r = cq_final(cq, z0, Qe, round(Int, t * 20)); zt = tc(t).z; g = hh_grad(zt); c = Float64(abs(dot(g, BigFloat.(r.Q[:, 1]))) / norm(g))
                rhs = Float64(log(abs(dot(g0, BigFloat.(Qe[:, 1])))) - log(abs(dot(g, BigFloat.(r.Q[:, 1]))))); e = abs(r.fin.l[1] - rhs); B = 3.0e-2 + 8 * cq_eps(Float64) / c
                emit("8", "Oid: F3 orbit $name, t = $t, q = 4, h = 1/20: |l_1 - (ln|gradH(z_0).q_1(0)| - ln|gradH(z_t^Taylor).q_1(t)|)| (c = $(fmt(c)))", e, "<= $(fmt(B)) (plan U8; spec 1e-3 + cancellation: ratio $(fmt(e / (1.0e-3 + 8 * cq_eps(Float64) / c))))", pf(e <= B))
            end
        end
    end
end

# ------------------------------------------------------------------------------------------------ Od20: Taylor BigFloat(256) variational reference (criterion 20), R2, R4, convergence data
section("Od20") do
    say("# Od20: reference = Taylor order 40, abstol 2^-200, 256 bits, state + k variational columns of the seeded frame Q0 (hand-written field), l_ref from the Gram determinants of Y Q0 (C16, 256 bits); order-8 table on F3, order-4 table on the pendulum; bounds: plan AC20 (q = 4: 5e-3, pendulum A = 3: 1e-1; q = 2: 1.0, pendulum A = 3 skipped) with the spec delta_ext (1e-3, 1e-1) reported next to them")
    setprecision(BigFloat, BITS) do
        rng = StableRNG(VSEED + 20); h0 = CQ_H; keep = Dict{Any, Any}()
        cases = (("F3 R", HenonHeilesFlow(), cq_orbit_R(), :o8, 4, (10, 30, 60), :hh), ("F3 C", HenonHeilesFlow(), cq_orbit_C(), :o8, 4, (10, 30, 60), :hh),
                 ("pendulum A=2", Pendulum(1.0), cq_pendulum_start(2.0), :o4, 2, (10, 30, 60, 100), :pend), ("pendulum A=3", Pendulum(1.0), cq_pendulum_start(3.0), :o4, 2, (10, 30, 60, 100), :pend))
        for (name, flow, z0, key, k, ts, kind) in cases
            d = dimension(flow); Q0 = cq_frame(rng, d, k, Float64); t0 = time(); tc = tay_case(kind, z0, Q0, maximum(ts))
            ref = Dict(t => (r = tc(t); Float64.(lsplit_big(r.W, BigFloat.(Matrix(Q0))))) for t in ts)
            say(@sprintf("# Od20: %s Taylor reference to t = %d, %.1f s", name, maximum(ts), time() - t0))
            keep[name] = (z0 = z0, Q0 = Q0, ref = ref, key = key, flow = flow)
            for h in (1 / 20, 1 / 40), t in ts
                im = cq_im(flow, Float64, key, h); N = round(Int, t / h)
                ld = cq_discrete(im, z0, Q0, N, h) * t; ed = maximum(abs, ld - ref[t]); errs = Dict{Int, Float64}()
                for q in (2, 4)
                    cq = cq_cqr(im, flow; order = q); e = maximum(abs, cq_final(cq, z0, Q0, N).fin.l - ref[t]); errs[q] = e
                    B = q == 4 ? (name == "pendulum A=3" ? 1.0e-1 : 5.0e-3) : 1.0
                    if q == 2 && name == "pendulum A=3"
                        emit("20", "Od20: $name, q = 2, h = 1/$(round(Int, 1 / h)), t = $t: max_i |l_i - l_i,ref| (skipped by plan U7)", e, "report", "INFO")
                    else
                        emit("20", "Od20: $name, q = $q, h = 1/$(round(Int, 1 / h)), t = $t: max_i |l_i - l_i,ref| (Taylor 256 bits)", e, "<= $B (plan); h = 1/20 spec delta_ext $(q == 4 ? "1e-3" : "1e-1"): ratio $(fmt(e / (q == 4 ? 1.0e-3 : 1.0e-1)))", pf(e <= B))
                    end
                end
                emit("20", "Od20: $name, h = 1/$(round(Int, 1 / h)), t = $t: discrete RGS spectrum error against the same reference ; continuous q = 4 ; q = 2", string(fmt(ed), " ; ", fmt(errs[4]), " ; ", fmt(errs[2])), "report", "INFO")
            end
        end
        # R4 and the convergence data: orbit C at t = 60, order-8 table, h = 1/10 .. 1/320
        c = keep["F3 C"]; hs = [1 / 10, 1 / 20, 1 / 40, 1 / 80, 1 / 160, 1 / 320]; flow = c.flow; conv = Dict{Int, Vector{Float64}}(2 => Float64[], 4 => Float64[]); emax = Dict{Int, Vector{Float64}}(2 => Float64[], 4 => Float64[])
        for q in (2, 4), h in hs
            im = cq_im(flow, Float64, :o8, h); cq = cq_cqr(im, flow; order = q, orthogonality_tol = 0.5); r = cq_final(cq, c.z0, c.Q0, round(Int, 60 / h))
            push!(conv[q], maximum(abs, r.fin.l - c.ref[60])); push!(emax[q], r.fin.Emax)
        end
        say("# Od20: convergence sweep uses orthogonality_tol = 0.5 (a guard, not a bound; the default 1e-2 / 5e-2 would stop the coarse order-2 runs)")
        for q in (2, 4)
            ratios = conv[q][1:(end - 1)] ./ conv[q][2:end]; sl = cq_slope(log.(hs), log.(conv[q]))
            sl_ok = cq_slope(log.(hs[2:4]), log.(conv[q][2:4]))
            emit("20", "Od20: orbit C, t = 60, order-8 table, q = $q, h = 1/10 .. 1/320: errors $(join(map(fmt, conv[q]), ", ")); E_max seen $(join(map(fmt, emax[q]), ", "))", string("ratios ", join(map(fmt, ratios), ", ")), "report; expected about $(2^q)", "INFO")
            emit("20", "Od20: orbit C, t = 60, q = $q: slope of ln error against ln h over all h ; over h = 1/20 .. 1/80", string(fmt(sl), " ; ", fmt(sl_ok)), "report (expected $q)", "INFO")
        end
        e4 = conv[4][2]; j = findfirst(i -> conv[2][i] <= e4, eachindex(hs))
        hreach = j === nothing ? "not reached down to h = 1/320" : j == 1 ? "already at h = 1/10" : string("about h = ", fmt(exp(log(hs[j - 1]) + (log(e4) - log(conv[2][j - 1])) / (log(conv[2][j]) - log(conv[2][j - 1])) * (log(hs[j]) - log(hs[j - 1])))))
        emit("R4", "order-2 against order-4 error at h = 1/20, orbit C, t = 60: q = 2 ; q = 4 ; h at which the order-2 error reaches the order-4 error at h = 1/20", string(fmt(conv[2][2]), " ; ", fmt(e4), " ; ", hreach), "no threshold", "INFO")
        serialize(joinpath(SCR, "conv-c.jls"), (hs = hs, errs = conv))
        # R2: equilibrium defect and exponent bias against h, c_nu and the form, orbit C, t = 60, order-4 scheme
        for (k, form) in ((4, :f84), (2, :f92)), h in (1 / 20, 1 / 40), cnu in (1.0, 2.0, 4.0)
            im = cq_im(flow, Float64, :o8, h); cq = cq_cqr(im, flow; order = 4, penalty_factor = cnu, form = form, orthogonality_tol = 0.1)
            r = cq_final(cq, c.z0, SMatrix{4, k, Float64}(c.Q0[:, 1:k]), round(Int, 60 / h)).fin; bias = maximum(abs, r.l - c.ref[60][1:k])
            emit("R2", "orbit C, t = 60, q = 4, k = $k, form $(form === :f84 ? 84 : 92), h = 1/$(round(Int, 1 / h)), c_nu = $cnu: E_eq (E_max seen) ; exponent error against the Taylor reference ; error / E_eq", string(fmt(r.Emax), " ; ", fmt(bias), " ; ", fmt(bias / r.Emax)), "no threshold", "INFO")
        end
    end
end

section("Od20R") do
    setprecision(BigFloat, BITS) do
        rng = StableRNG(VSEED + 21); z0 = cq_orbit_R(); Q0 = cq_frame(rng, 4, 4, Float64); T = 1000; t0 = time(); tc = tay_case(:hh, z0, Q0, T)
        ts = (60, 100, 300, 1000); ref = Dict(t => (r = tc(t); Float64.(lsplit_big(r.W, BigFloat.(Matrix(Q0))))) for t in ts)
        say(@sprintf("# Od20R: orbit R Taylor reference to t = %d (4 columns), %.1f s", T, time() - t0))
        flow = HenonHeilesFlow()
        for (h, q) in ((1 / 20, 4), (1 / 40, 4), (1 / 20, 2)), t in ts
            im = cq_im(flow, Float64, :o8, h); cq = cq_cqr(im, flow; order = q); r = cq_final(cq, z0, Q0, round(Int, t / h)); e = maximum(abs, r.fin.l - ref[t])
            ld = maximum(abs, cq_discrete(im, z0, Q0, round(Int, t / h), h) * t - ref[t])
            emit("20", "Od20R: orbit R, q = $q, h = 1/$(round(Int, 1 / h)), t = $t: max_i |l_i - l_i,ref| ; discrete RGS error (reported)", string(fmt(e), " ; ", fmt(ld)), "report (to t = 1000)", "INFO")
        end
    end
end

# ------------------------------------------------------------------------------------------------ Oc9: same-orbit comparison at T = 1000 (criterion 9, R6)
const C9_CACHE = Ref{Any}(nothing)
section("Oc9") do
    flow = HenonHeilesFlow(); cls = cq_classes(20, 1)
    starts = vcat([("R", cq_orbit_R())], [("C", cq_orbit_C())], [("reg$i", z) for (i, z) in enumerate(cls.regular[1:20])], [("cha$i", z) for (i, z) in enumerate(cls.chaotic[1:20])])
    say("# Oc9: orbits R, C and 20 starts per class (cq_classes, SALI at 2000 below 1e-8), order-8 table, q = 4, T = 1000; lambda form |lambda^C - lambda^D| against c_d h^4 with c_d = 5 (spec) and the l form |l^C - l^D| against 5e4 h^4 (plan U9, calibrated at T = 200); discrete = repeated Gram-Schmidt k = 4 at tau = 1 divided by h")
    T = 1000.0; agg = Dict{Any, Any}()
    for (j, (name, z0)) in enumerate(starts), h in (1 / 20, 1 / 40)
        im = cq_im(flow, Float64, :o8, h); Q0 = initial_deviations(VSEED + 90, j, Val(4), Val(4), Float64); N = round(Int, T / h)
        lc = propagate_cqr_orbit(cq_cqr(im, flow; order = 4), z0, Q0, N)[3].lambda
        ld = cq_discrete(im, z0, Q0, N, h; alg = RGS()); lh = cq_discrete(im, z0, Q0, N, h; alg = HouseholderQR()); lb = cq_benettin(im, z0, Q0[:, 1], N, h)
        e = maximum(abs, lc - ld); eh = maximum(abs, lc - lh); eb = abs(lc[1] - lb); sp = maximum(abs, ld - lh)
        cl = startswith(name, "reg") ? "regular" : startswith(name, "cha") ? "chaotic" : name
        a = get!(agg, (cl, h), (lam = 0.0, l = 0.0, eh = 0.0, eb = 0.0, sp = 0.0, c = 0.0))
        agg[(cl, h)] = (lam = max(a.lam, e / (5 * h^4)), l = max(a.l, e * T / (5.0e4 * h^4)), eh = max(a.eh, eh / (5 * h^4)), eb = max(a.eb, eb / (5 * h^4)), sp = max(a.sp, sp), c = max(a.c, e / h^4))
    end
    for ((cl, h), a) in sort!(collect(agg); by = x -> (string(x[1][1]), -x[1][2]))
        hs = "h = 1/$(round(Int, 1 / h))"
        emit("9", "Oc9: $cl, $hs, T = 1000, q = 4: worst |lambda^C - lambda^D(RGS)| / (5 h^4) (lambda form of the spec; realized c_d = $(fmt(a.c))); D121: reported at T = 1000, scored at T = 60 (FULL2)", a.lam, "report (spec c_d = 5; D121 scores T = 60 only)", "INFO")
        emit("9", "Oc9: $cl, $hs, T = 1000: worst |l^C - l^D| / (5e4 h^4) (l form of plan U9); D121: reported", a.l, "report", "INFO")
        emit("9", "Oc9: $cl, $hs, T = 1000: worst vs Householder / (5 h^4) ; vs Benettin column 1 / (5 h^4) ; R6 spread among the discrete variants max|RGS - Householder| (compare with the difference to the continuous value)", string(fmt(a.eh), " ; ", fmt(a.eb), " ; ", fmt(a.sp)), "report (R6)", "INFO")
    end
end

# ------------------------------------------------------------------------------------------------ AC11d: diagnosis of the AC11 failure at the full size (which assertion, by how much)
section("AC11d") do
    flow, im = cq_hh(); cq = cq_cqr(im, flow; order = 4); T = 2000.0; n = round(Int, T / CQ_H); A = cq_classes(CQ_NCLASS, 1).chaotic; allow = 5.0e4 * CQ_H^4 / T
    npos = 0; ndisc = 0; ratios = Float64[]; lam1 = Float64[]
    for (j, z0) in enumerate(A)
        Q0 = initial_deviations(cq_seed(11, 1), j, Val(4), Val(4), Float64)
        lam = propagate_cqr_orbit(cq, z0, Q0, n)[3].lambda; dis = cq_discrete(im, z0, Q0, n, CQ_H)
        lam[1] > 0 && (npos += 1); push!(ratios, abs(lam[1] - dis[1]) / allow); push!(lam1, lam[1])
        abs(lam[1] - dis[1]) <= allow && (ndisc += 1)
    end
    emit("11", "AC11d: $(length(A)) chaotic starts, T = 2000: starts with lambda_1 > 0 ; starts within the allowance 5e4 h^4 / T of the discrete lambda_1 (D121: at least 95 percent)", string(npos, " ; ", ndisc), "all $(length(A)) ; >= $(ceil(Int, 0.95 * length(A)))", pf(npos == length(A) && ndisc >= 0.95 * length(A)))
    emit("11", "AC11d: |lambda_1^C - lambda_1^D| / allowance: median ; q90 ; max (reported, D121) ; smallest lambda_1^C", string(fmt(median(ratios)), " ; ", fmt(quantile(ratios, 0.9)), " ; ", fmt(maximum(ratios)), " ; ", fmt(minimum(lam1))), "report", "INFO")
end

# ------------------------------------------------------------------------------------------------ Oord: discretization order against the discrete spectrum
section("Oord") do
    flow = HenonHeilesFlow(); hs = [1 / 10, 1 / 20, 1 / 40, 1 / 80]; T = 60.0
    say("# Oord: error |lambda^C - lambda^D| (max over i, k = 4, RGS at tau = 1 over h) against h on orbits R and C, T = 60; allowance c h^m with m = min(p, q); tables order 8 (KahanLi) and order 4 (ForestRuth)")
    for (name, z0) in HH_ORBITS, (key, p) in ((:o8, 8), (:o4, 4)), q in (2, 4)
        m = min(p, q); Q0 = initial_deviations(VSEED + 91, 1, Val(4), Val(4), Float64); errs = Float64[]
        for h in hs
            im = cq_im(flow, Float64, key, h); N = round(Int, T / h)
            lc = propagate_cqr_orbit(cq_cqr(im, flow; order = q, orthogonality_tol = 0.5), z0, Q0, N)[3].lambda; push!(errs, maximum(abs, lc - cq_discrete(im, z0, Q0, N, h)))
        end
        ratios = errs[1:(end - 1)] ./ errs[2:end]; floor_ = 1.0e3 * cq_eps(Float64); sl = cq_slope(log.(hs), log.(errs))
        ok = all(i -> !(errs[i] > floor_ && errs[i + 1] > floor_) || (2^m / 4 <= ratios[i] <= 4 * 2^m), 1:3)
        gate = name == "C" && key === :o8 && q == 4
        emit("9", "Oord: orbit $name, table order $p, q = $q (m = $m), T = 60, h = 1/10 .. 1/80: errors $(join(map(fmt, errs), ", ")); fitted slope $(fmt(sl))", string("ratios ", join(map(fmt, ratios), ", ")), "each in [$(2^m / 4), $(4 * 2^m)] where above the floor $(fmt(floor_))$(gate ? " (spec gate: orbit C, q = 4)" : " (reported)")", gate ? pf(ok) : "INFO")
    end
end

# ------------------------------------------------------------------------------------------------ Od21: ChaosTools cross-check (criterion 21)
section("Od21") do
    flow = HenonHeilesFlow(); t_end = 500.0; n = 100
    f(u, p, t) = SVector(u[3], u[4], -(u[1] + 2 * u[1] * u[2]), -(u[2] + u[1]^2 - u[2]^2))
    jac(u, p, t) = SMatrix{4, 4, Float64}(0, 0, -(1 + 2 * u[2]), -2 * u[1], 0, 0, -2 * u[1], -(1 - 2 * u[2]), 1, 0, 0, 0, 0, 1, 0, 0)
    function chaotic_starts(seed, m)
        a, b = section_interval(flow, CQ_E); w = sqrt(2 * CQ_E); rng = StableRNG(seed); out = SVector{4, Float64}[]
        while length(out) < m
            z, st = random_energy_start(flow, CQ_E, (a, b), (-w, w), rng); st === :ok || continue
            cq_sali(z, 2000.0) < 1.0e-8 && push!(out, SVector{4, Float64}(z))
        end
        return out
    end
    A = chaotic_starts(VSEED + 210, n); B = chaotic_starts(VSEED + 211, n)
    say("# Od21: ChaosTools $(pkgversion(CT)) lyapunovspectrum (method H2 of Geist1990, hand-coded Jacobian, CoupledODEs Tsit5 abstol = reltol = 1e-12, Delta t = 1, identity frame), $n independent chaotic-class starts (SALI at 2000 below 1e-8) on each side, t = $t_end; module: q = 4, order-8 table, h = 1/20, identity frame")
    im = cq_im(flow, Float64, :o8, CQ_H); cq = cq_cqr(im, flow; order = 4)
    lc = [propagate_cqr_orbit(cq, A[j], cq_eye(4, 4, Float64), round(Int, t_end / CQ_H))[3].lambda for j in 1:n]
    ode = CT.CoupledODEs(f, B[1], nothing; diffeq = (abstol = 1.0e-12, reltol = 1.0e-12)); tds = CT.TangentDynamicalSystem(ode; J = jac, k = 4)
    lt = Vector{Vector{Float64}}(); t0 = time()
    for j in 1:n
        push!(lt, CT.lyapunovspectrum(tds, round(Int, t_end); Δt = 1.0, u0 = B[j]))
    end
    say(@sprintf("# Od21: ChaosTools side %.1f s", time() - t0))
    for i in 1:2
        a = [l[i] for l in lc]; b = [l[i] for l in lt]; se = sqrt(cq_var(a) / n + cq_var(b) / n); allow = 4 / t_end
        d = abs(cq_mean(a) - cq_mean(b))
        emit("21", "Od21: lambda_$i ensemble means: module $(fmt(cq_mean(a))) (sd $(fmt(sqrt(cq_var(a))))) ; ChaosTools $(fmt(cq_mean(b))) (sd $(fmt(sqrt(cq_var(b))))): |difference| ; 4 combined SE + 4 / t", string(fmt(d), " ; ", fmt(4 * se + allow)), "<= 4 SE + 4 / t = $(fmt(4 * se + allow)) (difference / SE = $(fmt(d / se)))", pf(d <= 4 * se + allow))
    end
    emit("21", "Od21: lambda_1 + lambda_4 and lambda_2 + lambda_3 means, module ; ChaosTools (pairing, reported)", string(fmt(cq_mean([l[1] + l[4] for l in lc])), ", ", fmt(cq_mean([l[2] + l[3] for l in lc])), " ; ", fmt(cq_mean([l[1] + l[4] for l in lt])), ", ", fmt(cq_mean([l[2] + l[3] for l in lt]))), "report", "INFO")
end

# ------------------------------------------------------------------------------------------------ AC17: 1 against N threads, bitwise; allocation at t = 1e3
function ac17_xs()
    flow = HenonHeilesFlow(); a, b = section_interval(flow, CQ_E); w = sqrt(2 * CQ_E); rng = StableRNG(VSEED + 17)
    return [SVector{4, Float64}(random_energy_start(flow, CQ_E, (a, b), (-w, w), rng)[1]) for _ in 1:4096]
end
function ac17_serial(q, chunk)
    flow = HenonHeilesFlow(); im = cq_im(flow, Float64); cq = cq_cqr(im, flow; order = q); xs = ac17_xs()
    res = propagate_cqr(cq, xs, 20_000; k = 4, seed = VSEED + 17, chunk = chunk)
    return (l = [f.l for f in res.finals], fin = res.finals, x = res.x, W = res.W)
end
section("AC17child") do
    serialize(get(ENV, "CQ_CHILD_OUT", joinpath(SCR, "child.jls")), [(q = q, r = ac17_serial(q, 64)) for q in (2, 4)])
end
section("AC17T") do
    flow = HenonHeilesFlow(); im = cq_im(flow, Float64); xs = ac17_xs(); n = length(xs)
    child = joinpath(SCR, "child-t1.jls")
    cmd = addenv(`$(Base.julia_cmd()) -t 1 --project=$(joinpath(ROOT, "verify")) $(@__FILE__)`, "CQ_ONLY" => "AC17child", "CQ_CHILD_OUT" => child, "CQ_LOG" => joinpath(SCR, "child-t1.log"), "CQ_SCRATCH" => SCR)
    mine = Dict{Int, Any}(); t1 = @elapsed for q in (2, 4)
        cq = cq_cqr(im, flow; order = q); fin = Vector{Any}(undef, n); xo = Vector{Any}(undef, n); Wo = Vector{Any}(undef, n)
        Threads.@threads for j in 1:n
            xo[j], Wo[j], fin[j] = propagate_cqr_orbit(cq, xs[j], initial_deviations(VSEED + 17, j, Val(4), Val(4), Float64), 20_000; orbit = j)
        end
        mine[q] = (fin = fin, x = xo, W = Wo)
    end
    t2 = @elapsed run(cmd)
    theirs = deserialize(child)
    say(@sprintf("# AC17T: 4096 random starts of F3 at E = 1/8, t = 1e3 (20000 steps), both schemes: %d threads over orbits (propagate_cqr_orbit in Threads.@threads) %.1f s; 1 thread child (propagate_cqr, chunk 64) %.1f s; propagate_cqr itself is a serial loop over orbits (impl notes), so the threaded side is the orbit-level use", Threads.nthreads(), t1, t2))
    for q in (2, 4)
        a = mine[q]; b = [x for x in theirs if x.q == q][1].r
        emit("17", "AC17T: q = $q: final records (l, lambda, E, Emax, counts) bitwise identical, $(Threads.nthreads()) threads against 1 thread child (4096 starts)", isequal(a.fin, b.fin), "bitwise", pf(isequal(a.fin, b.fin)))
        emit("17", "AC17T: q = $q: final states and final frames bitwise identical", isequal(a.x, b.x) && isequal(a.W, b.W), "bitwise", pf(isequal(a.x, b.x) && isequal(a.W, b.W)))
        c = ac17_serial(q, 1000)
        emit("17", "AC17T: q = $q: chunk 64 (child) against chunk 1000 (this process, serial): bitwise identical records, states, frames", isequal(c.fin, b.fin) && isequal(c.x, b.x) && isequal(c.W, b.W), "bitwise", pf(isequal(c.fin, b.fin) && isequal(c.x, b.x) && isequal(c.W, b.W)))
    end
end
section("AC17alloc") do
    flow = HenonHeilesFlow(); z0 = cq_orbit_C(); worst = 0
    for T in (Float64, Double64), q in (2, 4), k in (1, 2, 4)
        im = cq_im(flow, T); cq = cq_cqr(im, flow; order = q); x0 = SVector{4, T}(map(T, z0)); Q0 = SMatrix{4, k, T}(initial_deviations(VSEED + 18, 1, Val(4), Val(4), T)[:, 1:k])
        b = cq_alloc_bytes(propagate_cqr_orbit, cq, x0, Q0, 20_000); worst = max(worst, b)
        emit("17", "AC17alloc: $T, q = $q, k = $k, t = 1e3 (20000 steps): bytes allocated after warm-up", b, "0", pf(b == 0))
    end
end

# ------------------------------------------------------------------------------------------------ AC22: device part not run; CPU fallback record
section("AC22") do
    flow = HenonHeilesFlow(); z0 = cq_orbit_C()
    for T in (Float64, Double64), bk in (:gpu, :auto, :cpu)
        im = cq_im(flow, T); cq = cq_cqr(im, flow; order = 4); xs = [SVector{4, T}(map(T, z0))]
        r = propagate_cqr(cq, xs, 400; k = 4, seed = 1, backend = bk); p = r.provenance
        emit("22", "AC22: $T, request :$bk -> backend $(p.backend), device $(p.device), fallback reason: $(p.backend_fallback_reason)", string(p.backend), "backend :cpu", pf(p.backend === :cpu))
    end
    emit("22", "AC22: CPU against GPU, Float64: device part and 1024 + 1024 starts at t = 400", "not run", "-", "NOT RUN: needs verify/gpu.jl (writes data/gpu-benchmark.toml); owner decision. The module has no kernel path (impl notes): a :gpu request falls back to the CPU with the reason of a CPU-only reducer, not the Double64 kernel reason of the spec Compute section")
end

# ------------------------------------------------------------------------------------------------ R1: wall-clock ratios to the discrete QR (reported)
section("R1") do
    flow = HenonHeilesFlow(); z0 = cq_orbit_C(); N = 20_000
    say("# R1: orbit C, N = $N steps (t = 1000), h = 1/20, order-8 table, k = d = 4, tau = 1, minimum of 5 repeats, 1 thread; the machine may be loaded by other processes")
    for T in (Float64, Double64)
        im = cq_im(flow, T); x0 = SVector{4, T}(map(T, z0)); Q0 = initial_deviations(VSEED + 1, 1, Val(4), Val(4), T)
        tmin(f) = (f(); minimum(@elapsed(f()) for _ in 1:5))
        td = Dict{Any, Float64}()
        td[:rgs] = tmin(() -> propagate(im, [x0], N; k = 4, tau = 1, seed = 1, reducers = (QRSpectrum(RGS()),), W0 = [Q0]))
        td[:hh] = tmin(() -> propagate(im, [x0], N; k = 4, tau = 1, seed = 1, reducers = (QRSpectrum(HouseholderQR()),), W0 = [Q0]))
        for q in (2, 4)
            cq = cq_cqr(im, flow; order = q); tq = tmin(() -> propagate_cqr_orbit(cq, x0, Q0, N))
            emit("R1", "$T, q = $q: continuous QR seconds ; discrete RGS ; discrete Householder ; ratios to RGS and to Householder (Geist1990 Table III p.886 reports 40 against Gram-Schmidt, a reproduction target only)", string(fmt(tq), " ; ", fmt(td[:rgs]), " ; ", fmt(td[:hh]), " ; ", fmt(tq / td[:rgs]), " ; ", fmt(tq / td[:hh])), "no threshold", "INFO")
        end
    end
end

# ------------------------------------------------------------------------------------------------ FIG
section("FIG") do
    ec = deserialize(joinpath(SCR, "conv-closed.jls")); cc = deserialize(joinpath(SCR, "conv-c.jls"))
    fig = CM.Figure(size = (1300, 520)); col = Dict(2 => :royalblue, 4 => :crimson)
    for (j, (title, hs, errs)) in enumerate((("coupled hyperbolic fixture (C19 e), BigFloat(256), t = 10", ec.hs, ec.errs), ("Henon-Heiles orbit C, t = 60, order-8 table, Taylor 256-bit reference", cc.hs, cc.errs)))
        ax = CM.Axis(fig[1, j]; xscale = log10, yscale = log10, xlabel = "h (orbit step)", ylabel = "max_i |l_i - l_i,ref|", title = title, titlesize = 13)
        for q in (2, 4)
            e = errs[q]; ok = findall(x -> isfinite(x) && x > 0, e); hh = collect(hs)[ok]; ee = e[ok]
            sl = cq_slope(log.(hh), log.(ee)); CM.scatterlines!(ax, hh, ee; color = col[q], label = "order $q scheme, fitted slope $(round(sl, digits = 2))")
            CM.lines!(ax, hh, ee[end] .* (hh ./ hh[end]) .^ q; color = col[q], linestyle = :dash, label = "h^$q reference")
        end
        CM.axislegend(ax; position = :rb, labelsize = 11)
    end
    CM.save(joinpath(OUT, "continuous-qr-convergence.png"), fig)
    # orthogonality: E(t) = ||Q'Q - I||_2 on F3 orbit C, with and without the penalty; left: undisturbed run, right: initial defect 1e-3 (form 84, k = 4)
    flow, im = cq_hh(); z0 = cq_orbit_C(); Q4 = initial_deviations(VSEED + 5, 1, Val(4), Val(4), Float64); Q2 = SMatrix{4, 2, Float64}(Q4[:, 1:2]); n = 12_000
    En(Q) = opnorm(Matrix(Q' * Q - I))
    function run_E(Q0, form, cnu; n = n)
        cq = cq_cqr(im, flow; order = 4, form = form, penalty_factor = cnu, orthogonality_tol = 1.0e300)
        ser = try
            propagate_cqr_series(cq, z0, Q0, n; every = 20)[4]
        catch e
            say("# FIG: run $form c_nu = $cnu stopped: $(sprint(showerror, e)[1:min(end, 160)])"); NamedTuple[]
        end
        return [s.n * CQ_H for s in ser], [max(En(s.Q), 1.0e-17) for s in ser]
    end
    fig2 = CM.Figure(size = (1300, 520))
    ax1 = CM.Axis(fig2[1, 1]; yscale = log10, xlabel = "t", ylabel = "E(t) = ||Q'Q - I||_2", title = "undisturbed start, orbit C, order-4 scheme, h = 1/20")
    ax2 = CM.Axis(fig2[1, 2]; yscale = log10, xlabel = "t", ylabel = "E(t) = ||Q'Q - I||_2", title = "initial defect delta = 1e-3 (form 84, k = 4)")
    for (Q0, form, cnu, c, lab) in ((Q4, :f84, 1.0, :royalblue, "form 84, k = 4, c_nu = 1"), (Q4, :f84, 0.0, :navy, "form 84, k = 4, no penalty"), (Q2, :f92, 2.0, :crimson, "form 92, k = 2, c_nu = 2"), (Q2, :f92, 0.0, :darkorange, "form 92, k = 2, no penalty"))
        t, E = run_E(Q0, form, cnu; n = (form === :f92 && cnu == 0.0) ? 3700 : n); isempty(t) || CM.lines!(ax1, t, E; color = c, label = lab)   # form 92 without penalty overflows at step 3878: stop earlier
    end
    CM.axislegend(ax1; position = :lt, labelsize = 11)
    rngd = StableRNG(VSEED + 55); dl = cq_sym_delta(rngd, 4, 1.0e-3); Qi = SMatrix{4, 4, Float64}(Q4 * (I + SMatrix{4, 4, Float64}(dl)))
    for (cnu, c, lab) in ((1.0, :royalblue, "c_nu = 1"), (0.0, :navy, "no penalty"))
        t, E = run_E(Qi, :f84, cnu); CM.lines!(ax2, t, E; color = c, label = lab)
    end
    CM.axislegend(ax2; position = :rt, labelsize = 11)
    CM.save(joinpath(OUT, "continuous-qr-orthogonality.png"), fig2)
    say("# FIG: wrote continuous-qr-convergence.png and continuous-qr-orthogonality.png in the stage 04 output folder")
end
