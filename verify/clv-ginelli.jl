# Stage 04 verification of module clv-ginelli (catalog row L12: Chaoslab.Lyapunov clv_ginelli.jl, clv_backward.jl; spec D122, test amendment D123).
# Run from the repo root, one process per group (CLV_ONLY, comma list; CLV_LOG names the log, CLV_APPEND=1 appends):
#   CLV_ONLY=FULL1 .. FULL3  julia --project=verify -t 8 verify/clv-ginelli.jl   (the three package test files of test/clv-ginelli/ at the FULL sizes, CLV_FULL_RUN = true; one file each)
#   CLV_ONLY=Od13a  (criterion 13, standard map K = 5: 1 start + 10 chaotic starts, Float64 and Double64 against the 256-bit static null-space algorithm; 256-bit Ginelli within 1e-30)
#   CLV_ONLY=Od13b  (criterion 13, coupled maps d = 4)
#   CLV_ONLY=Od13c, Od13d  (criterion 13, F3 orbit C, Float64 and Double64 module against the static algorithm; zero block as the span of (v_2, v_3))
#   CLV_ONLY=Od13e  (criterion 13, F3 orbit C, the module in BigFloat(256) against the static algorithm within 1e-30, long transients)
#   CLV_ONLY=Oe     (O-e and R6: KuptsovParlitz2012 3D Henon map exponents, reported only)
#   CLV_ONLY=R1     (wall-clock of the forward and backward phases and the checkpointed overhead; the static-algorithm times are in the Od13 groups)
#   CLV_ONLY=R2     (angle distributions over 1e5 samples; writes the data of the angles figure), R3idx (worst-conditioned vector index), R4 (transient actually needed against 40 and 25 / (g Dt)), R5 (orthonormalizer and tau dependence)
#   CLV_ONLY=THR    (1 against N threads, bitwise; spawns a one-thread child, group THRchild), FIG (the two figures)
# Not run: no device part (the module is CPU only, P9). Nothing here writes under data/.
# Spec: workspaces/module-build/stages/01-spec/output/clv-ginelli-spec.md; test plan points U1 to U15 and Amendment D123 are the operative reading. Verify-only seeds 20261012 + 45900 ... (package seeds 45000 + 50 ac + sub).
# One printed line per measurement:  <tag> | quantity | value | bound | PASS/FAIL/INFO
const ONLY = filter(!isempty, split(get(ENV, "CLV_ONLY", ""), ","))
using Printf, Random, LinearAlgebra, Statistics, Serialization, Logging
using Test
const CLV_FULL_RUN = true
using Chaoslab, StaticArrays, DoubleFloats, StableRNGs
want(name) = name in ONLY || (isempty(ONLY) && !startswith(name, "THRchild")) || (startswith(name, "Od13") && "Od13" in ONLY)
(isempty(ONLY) || "FIG" in ONLY) && (import CairoMakie as CM)

const ROOT = dirname(@__DIR__)
const OUT = joinpath(ROOT, "workspaces", "module-build", "stages", "04-verify", "output")
const SCR = get(ENV, "CLV_SCRATCH", joinpath(tempdir(), "clv-verify-scratch"))
mkpath(SCR)
const LOG = open(get(ENV, "CLV_LOG", joinpath(OUT, "clv-ginelli-verify-log.txt")), get(ENV, "CLV_APPEND", "") == "1" ? "a" : "w")

fmt(x::AbstractFloat) = @sprintf("%.4g", x)
fmt(x) = string(x)
const LOGLOCK = ReentrantLock()
function say(s)
    lock(LOGLOCK) do
        println(s); println(LOG, s); flush(stdout); flush(LOG)
    end
end
emit(tag, quantity, value, bound, verdict) = say(string(tag, " | ", quantity, " | ", fmt(value), " | ", bound, " | ", verdict))
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

say("# clv-ginelli verification, julia $(VERSION), threads $(Threads.nthreads()), base seed 20261012 + 45000, CLV_FULL_RUN = true, groups $(isempty(ONLY) ? "all" : join(ONLY, ",")), time $(round(time())) (unix)")
say("# timings in this log were taken on a machine that may be loaded by other processes (several groups run in parallel processes)")

const TD = joinpath(ROOT, "test", "clv-ginelli")
include(joinpath(TD, "clv_helpers.jl"))
const VSEED = 20261012 + 45900
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
function report_ts(tag, ts)
    el = ts.time_end - ts.time_start
    c = Test.get_test_counts(ts)
    ok = c.cumulative_fails == 0 && c.cumulative_errors == 0
    emit(tag, string("file total (", c.cumulative_passes, " pass, ", c.cumulative_fails, " fail, ", c.cumulative_errors, " error, ", c.cumulative_broken, " broken, ", @sprintf("%.1f", el), " s)"),
         c.cumulative_passes, "0 fail, 0 error", pf(ok))
    for ch in ts.results
        ch isa Test.AbstractTestSet || continue
        nf = nfail(ch)
        emit(tag, string("testset ", first(ch.description, 70), " (", npass(ch), " pass, ", nf, " fail or error, ", nbroken(ch), " broken, ", @sprintf("%.1f", ch.time_end - ch.time_start), " s)"), nf, "0", pf(nf == 0))
    end
end
function suite_section(name, file)
    section(name) do
        try
            Test.@testset "outer $name" begin
                ts = Test.@testset "$name" begin
                    include(joinpath(TD, file))
                end
                report_ts(name, ts)
            end
        catch e
            e isa Test.TestSetException || rethrow()
        end
    end
end
suite_section("FULL1", "clv_ac01_05.jl")      # criteria 1 to 5
suite_section("FULL2", "clv_ac06_10.jl")      # criteria 6 to 10
suite_section("FULL3", "clv_ac11_14.jl")      # criteria 11 to 14
# DIAGNOSTIC, not a size reduction and not a package-test change: at the full size (20 F3 chaotic starts) the test-side pilot of n* (clv_pilot, 300 handoffs) returns lambda_1 < lambda_2 for start 4
# (a sticky start), n* = -6548 and the module rightly throws ArgumentError (n_f < 0), which ends AC3b, AC5 and AC7 early (test-side reading U3). FULL1p and FULL2p rerun the two files
# with the pilot lengthened to 1000 handoffs (the helper redefined in this process only; the files themselves are untouched) to obtain the evidence the pilot error withheld.
function lengthen_pilot()
    Core.eval(Main, :(function clv_pilot(map, x0, tau::Integer, h; ns::Integer = 1000)
        d = length(x0); pts = clv_orbit(map, x0, tau * (ns + 1))
        G0 = clv_frame(1, 1, d, d, Float64); C0 = clv_eye(d, d, Float64)
        return clv_ref_run(map, pts, tau, 0, ns, 1, G0, C0, Float64).lnR ./ (ns * tau * h)
    end))
    say("# DIAGNOSTIC: clv_pilot redefined with ns = 1000 handoffs (default was 300)")
end
want("FULL1p") && lengthen_pilot()
want("FULL2p") && lengthen_pilot()
suite_section("FULL1p", "clv_ac01_05.jl")
suite_section("FULL2p", "clv_ac06_10.jl")

# ---------------------------------------------------------------------------------------------------
# AC3bd: diagnosis of the AC3b failures of FULL1p: per F3 chaotic start, the distance of criterion 3 (b) at the pilot transient N_b (pilot 1000 handoffs, cap 6000) and at 2, 3 N_b (no cap) and the finite-time
# exponents of the 1000-handoff pilot against those of a 6000-handoff run (the sticky starts change their spectrum with the horizon)
section("AC3bd") do
    T = Float64; im = clv_im(T)
    for (is, z0) in enumerate(clv_chaotic_starts(20, 1))
        lam = clv_pilot(im, z0, 20, CLV_H; ns = 1000); lam6 = clv_pilot(im, z0, 20, CLV_H; ns = 6000); Nb = clv_nstar(lam, 1.0; skip23 = true); nsb = 6; out = Float64[]
        for mult in (1, 2, 3)
            runs = []
            for (sg, sc) in ((1, 1), (2, 2))
                W0 = clv_frame(clv_seed(3, 10 + sg), 1, 4, 4, T); C0 = clv_rand_C(StableRNG(clv_seed(3, 20 + sc)), 4, T)
                push!(runs, clv_go(im, z0, mult * Nb, nsb, mult * Nb; tau = 20, h = CLV_H, W0 = W0, C0 = C0, outputs = (:V,)),
                      clv_go(im, z0, 2 * mult * Nb, 2 * nsb, 2 * mult * Nb; tau = 10, h = CLV_H, W0 = W0, C0 = C0, outputs = (:V,)))
            end
            w = 0.0
            for q in 2:4, m in 0:nsb
                Vq = q % 2 == 1 ? runs[q].V[:, :, m + 1] : runs[q].V[:, :, 2m + 1]
                w = max(w, clv_vdist(Vq, runs[1].V[:, :, m + 1]; singles = [1, 4], block = 2:3))
            end
            push!(out, w)
        end
        emit("AC3b-diag", string("F3 chaotic start ", is, " max distance at 1, 2, 3 x N_b (N_b = ", Nb, ")"), out[1], string(join(fmt.(out), ", "), " ; pilot-1000 lambda ", join(fmt.(lam), ","), " ; pilot-6000 lambda ", join(fmt.(lam6), ",")), out[1] <= 1.0e-10 ? "PASS" : "FAIL")
    end
end

# ---------------------------------------------------------------------------------------------------
# Criterion 13: the static null-space algorithm (spec V19, KuptsovParlitz2012 Eq. (55) p.13) in BigFloat(256) on the realized orbit points.
# Forward Gram-Schmidt frame G+_n (from a random orthonormal start), backward Gram-Schmidt frame G-_n by the transposed tangent propagation run backward in time
# (Q-_n R = M_n^T Q-_{n+1}, from a different random start), P = (G-)^T G+ at the same point, v_j = G+(:, 1:j) a with P(1:j-1, 1:j) a = 0 (null vector by the SVD at 256 bits).
# It uses no backward recursion V6, no back-substitution with the stored R and no Ginelli C matrix.
function nullvec(P::AbstractMatrix{BigFloat})
    j = size(P, 2)
    a = try
        svd(P; full = true).V[:, j]
    catch
        Matrix(qr(Matrix(P')).Q * Matrix{BigFloat}(I, j, j))[:, j]         # shortcut: QR of the transposed constraint, same null vector for a full-rank (j-1) x j matrix
    end
    return a / norm(a)
end
function static_clvs(Ms::Vector{Matrix{BigFloat}}, nf::Integer, ns::Integer, nb::Integer, Gf::Matrix{BigFloat}, Gb::Matrix{BigFloat})
    d = size(Gf, 1); N = nf + ns + nb
    G = Vector{Matrix{BigFloat}}(undef, nf + ns + 1); G[1] = Gf
    for n in 1:(nf + ns)
        G[n + 1] = clv_qr(Ms[n] * G[n])[1]
    end
    Gm = Vector{Matrix{BigFloat}}(undef, N + 1); Gm[N + 1] = Gb
    for n in N:-1:(nf + 1)
        Gm[n] = clv_qr(Ms[n]' * Gm[n + 1])[1]
    end
    Vs = Vector{Matrix{BigFloat}}(undef, ns + 1)
    for s in 0:ns
        Gp = G[nf + s + 1]; P = Gm[nf + s + 1]' * Gp; V = zeros(BigFloat, d, d)
        V[:, 1] = Gp[:, 1]
        for j in 2:d
            a = nullvec(P[1:(j - 1), 1:j]); v = Gp[:, 1:j] * a; V[:, j] = v / norm(v)
        end
        Vs[s + 1] = V
    end
    return Vs
end
widen_mat(M, S = BigFloat) = map(x -> clv_conv(S, clv_big(x)), Matrix(M))
bigfr(seed, d) = Matrix{BigFloat}(map(clv_big, Matrix(clv_frame(seed, 1, d, d, Float64))))
# One case: module in element type of x0 (n_f = n_b = nf, window ns), the 256-bit tangent products along the realized orbit, the static algorithm and the 256-bit reference Ginelli (clv_ref) on the same Ms.
function od13_run(map, x0, tau, h; nf, ns = 50, singles, block, wseed)
    T = eltype(x0); d = dimension(map); N = 2nf + ns
    W0 = clv_frame(wseed, 1, d, d, T); C0 = clv_rand_C(StableRNG(wseed + 1), d, T)
    t_mod = @elapsed res = clv_go(map, x0, nf, ns, nf; tau = tau, h = h, W0 = W0, C0 = C0, outputs = (:V,))
    pts = clv_orbit(map, x0, N * tau)
    t_ms = @elapsed Ms = clv_interval_tangents(clv_bigmap(map), pts, tau, N, BigFloat)
    t_ref = @elapsed ref = clv_ref(Ms, widen_mat(W0), widen_mat(C0), nf, ns, nf)
    t_st = @elapsed Vs = static_clvs(Ms, nf, ns, nf, bigfr(wseed + 5, d), bigfr(wseed + 6, d))
    dm = maximum(n -> clv_vdist(res.V[:, :, n + 1], Vs[n + 1]; singles = singles, block = block), 0:ns)
    dr = maximum(n -> clv_vdist(clv_refV(ref, n), Vs[n + 1]; singles = singles, block = block), 0:ns)
    return (dmod = dm, dref = dr, t_mod = t_mod, t_ms = t_ms, t_ref = t_ref, t_st = t_st, nf = nf)
end
function od13_emit(label, T, r, tol, tol256; info = false)
    vd(ok) = info ? "INFO" : pf(ok)
    if T === BigFloat
        emit("AC13", string(label, " module BigFloat(256) vs static, max delta (n_f = n_b = ", r.nf, ")"), r.dmod, string("<= ", tol256), vd(r.dmod <= tol256))
        emit("AC13", string(label, " 256-bit Ginelli (clv_ref) vs static, max delta"), r.dref, string("<= ", tol256), vd(r.dref <= tol256))
    else
        emit("AC13", string(label, " module ", T, " vs static, max delta (n_f = n_b = ", r.nf, ")"), r.dmod, string("<= ", tol), vd(r.dmod <= tol))
        emit("AC13", string(label, " 256-bit Ginelli (clv_ref) vs static, max delta (n_f = n_b = ", r.nf, ")"), r.dref, "reported (transient n*, not 1e-30)", "INFO")
    end
    emit("R1", string(label, " ", T, " seconds: module / tangents 256 / Ginelli 256 (clv_ref) / static"), r.t_st / max(r.t_ref, 1.0e-9), string(@sprintf("%.3g", r.t_mod), " / ", @sprintf("%.3g", r.t_ms), " / ", @sprintf("%.3g", r.t_ref), " / ", @sprintf("%.3g", r.t_st), " (value = static over Ginelli 256 time)"), "INFO")
end
# n* = ceil(c / (g_min Dt)) from a pilot spectrum; the transient of the module comparison (c = 120, Amendment D124; c = 40 of U3 left the coupled maps at 6.0e-5 in Double64)
# and of the 1e-30 comparison (c = 180, D124; c = 90 left 1.25e-20 on the coupled maps, 1.4e-43 at twice the transient)
od13_pilot(map, x0, tau, h; ns = 1000) = clv_pilot(map, x0, tau, h; ns = ns)   # D124: 1000 handoffs
function od13_maps(sysname)
    if sysname === :std
        return (T -> clv_std(T), 1, 1.0, 2, 1:2, nothing, false)
    elseif sysname === :cm
        return (T -> clv_cm(T), 1, 1.0, 4, 1:4, nothing, false)
    else
        return (T -> clv_im(T), 20, CLV_H, 4, [1, 4], 2:3, true)
    end
end
function od13_system(sysname, label, x0_64; Ts, ctag, big256 = true, ns = 50, nf_cap = 6000, wseed = VSEED + 100)
    mk, tau, h, d, singles, block, skip = od13_maps(sysname)
    m64 = mk(Float64); lam = od13_pilot(m64, x0_64, tau, h; ns = 1000)
    sys_ok = lam[1] > (sysname === :std ? 0.3 : 0.0)
    nf = max(clv_nstar(lam, tau * h; skip23 = skip, nmax = nf_cap), 40)
    emit("AC13", string(label, " pilot exponents (per unit time, Float64, 1000 handoffs)"), lam[1], string(join(fmt.(lam), ", "), " ; n* = ", nf), "INFO")
    sys_ok || (say(string("# ", label, " skipped: pilot lambda_1 below 0.3 (regular or island start)")); return false)
    for T in Ts
        x0 = T === Float64 ? x0_64 : map_to(T, x0_64)
        r = od13_run(mk(T), x0, tau, h; nf = nf, ns = ns, singles = singles, block = block, wseed = wseed)
        od13_emit(label, T, r, 1.0e-8, 1.0e-30)
        if r.dmod > 1.0e-8        # diagnosis: does the distance fall with longer transients (convergence-limited) or stay (a real disagreement)?
            for mult in (2, 3, 4, 6)
                r2 = od13_run(mk(T), x0, tau, h; nf = mult * nf, ns = ns, singles = singles, block = block, wseed = wseed)
                od13_emit(string(label, " [diagnostic, ", mult, " x n*]"), T, r2, 1.0e-8, 1.0e-30; info = true)
                r2.dmod <= 1.0e-8 && break
            end
        end
    end
    if big256
        gaps = [lam[i] - lam[i + 1] for i in 1:(d - 1)]; skip && deleteat!(gaps, 2)
        n256 = max(ceil(Int, 180 / (minimum(gaps) * tau * h)), 100)
        xb = map_to(BigFloat, x0_64)
        r = od13_run(mk(BigFloat), xb, tau, h; nf = n256, ns = ns, singles = singles, block = block, wseed = wseed + 20)
        od13_emit(label, BigFloat, r, 1.0e-8, 1.0e-30)
        if r.dmod > 1.0e-30
            for mult in (2, 3, 4, 6)
                r2 = od13_run(mk(BigFloat), xb, tau, h; nf = mult * n256, ns = ns, singles = singles, block = block, wseed = wseed + 20)
                od13_emit(string(label, " [diagnostic, ", mult, " x n256]"), BigFloat, r2, 1.0e-8, 1.0e-30; info = true)
                r2.dmod <= 1.0e-30 && break
            end
        end
    end
    return true
end
section("Od13a") do
    rng = StableRNG(VSEED + 1)
    starts = clv_starts(rng, 2, 16, Float64); kept = 0
    for (i, x0) in enumerate(starts)
        kept >= 11 && break
        label = i == 1 ? "std start 1 (the single start)" : string("std chaotic start ", i)
        od13_system(:std, label, x0; Ts = (Float64, Double64), ctag = "std", big256 = kept < 3, wseed = VSEED + 100 + 10 * i) && (kept += 1)
    end
    emit("AC13", "standard map starts compared (1 + 10 chaotic)", kept, "11", pf(kept >= 11))
end
section("Od13b") do
    x0 = clv_starts(StableRNG(VSEED + 2), 4, 1, Float64)[1]
    od13_system(:cm, "coupled maps d = 4", x0; Ts = (Float64, Double64), ctag = "cm", wseed = VSEED + 300)
end
f3_x0 = () -> clv_orbit_C()
section("Od13c") do
    od13_system(:f3, "F3 orbit C", f3_x0(); Ts = (Float64,), ctag = "f3", big256 = false, wseed = VSEED + 400)
end
section("Od13d") do
    od13_system(:f3, "F3 orbit C", f3_x0(); Ts = (Double64,), ctag = "f3", big256 = false, wseed = VSEED + 410)
end
section("Od13e") do
    mk, tau, h, d, singles, block, skip = od13_maps(:f3); x0 = f3_x0()
    lam = od13_pilot(mk(Float64), x0, tau, h; ns = 1000); gaps = [lam[i] - lam[i + 1] for i in 1:3]; deleteat!(gaps, 2)
    n256 = ceil(Int, 180 / (minimum(gaps) * tau * h))
    emit("AC13", "F3 orbit C pilot exponents (per unit time, Float64, 1000 handoffs)", lam[1], string(join(fmt.(lam), ", "), " ; n256 = ", n256), "INFO")
    r = od13_run(mk(BigFloat), map_to(BigFloat, x0), tau, h; nf = n256, ns = 50, singles = singles, block = block, wseed = VSEED + 420)
    od13_emit("F3 orbit C", BigFloat, r, 1.0e-8, 1.0e-30)
end

# ---------------------------------------------------------------------------------------------------
# O-e and R6: KuptsovParlitz2012 generalized three-dimensional Henon map (Eq. (76)), a = 1.76, b = 0.1; test double of the map (no change to src), reported only
struct Henon3D{T} <: Chaoslab.AbstractMap
    a::T
    b::T
end
Chaoslab.dimension(::Henon3D) = 3
Chaoslab.forward(m::Henon3D, z::SVector{3}) = SVector(m.a - z[2]^2 - m.b * z[3], z[1], z[2])
Chaoslab.tangent(m::Henon3D, z::SVector{3, T}) where {T} = SMatrix{3, 3, T}(0, 1, 0, -2 * z[2], 0, 1, -m.b, 0, 0)     # rows (0, -2 x2, -b), (1, 0, 0), (0, 1, 0)
section("Oe") do
    m = Henon3D(1.76, 0.1); target = [0.225, 0.188, -2.716]; lnb = log(0.1)
    for (is, x0) in enumerate((SVector(0.1, 0.1, 0.1), SVector(-0.3, 0.2, 0.5)))
        clv = ClvGinelli(m; tau = 1, h = 1.0)
        t = @elapsed res = propagate_clv_orbit(clv, x0, 3000, 1_000_000, 3000; k = 3, seed = VSEED + 500 + is, outputs = ())
        ns = 1_000_000; lam = [mean(@view res.lambda[i, 1:ns]) for i in 1:3]
        bterm = (log(abs(det(Matrix(res.C[:, :, 1])))) - log(abs(det(Matrix(res.C[:, :, ns + 1]))))) / ns       # telescoping boundary term of V14 (c) divided by n_s
        for i in 1:3
            emit("R6", string("start ", is, " window mean lambda_", i, " (1e6 steps, per step) against KuptsovParlitz2012 printed value ", target[i]), lam[i], string("printed 3 digits: ", target[i]), "INFO")
        end
        emit("R6", string("start ", is, " forward exponents (lambda_fwd) sum minus ln b"), sum(res.lambda_fwd) - lnb, "reported (spec: 1e-12)", "INFO")
        emit("R6", string("start ", is, " CLV window mean sum minus ln b (raw)"), sum(lam) - lnb, "reported (carries the boundary term of V14 (c))", "INFO")
        emit("R6", string("start ", is, " CLV window mean sum minus ln b minus boundary term (ln det C_0 - ln det C_ns)/ns"), sum(lam) - lnb - bterm, "reported (spec: 1e-12)", "INFO")
        emit("R6", string("start ", is, " forward lambda_fwd"), res.lambda_fwd[1], join(fmt.(res.lambda_fwd), ", "), "INFO")
        emit("R6", string("start ", is, " seconds for 1e6 steps"), t, "-", "INFO")
        if is == 1        # independent check on the same window: plain LinearAlgebra QR of the tangent product, no Chaoslab estimator (frame start I: boundary effect of order 1e-5)
            z = res.x0; Q = Matrix{Float64}(I, 3, 3); acc = zeros(3)
            for _ in 1:ns
                F = qr(Matrix(tangent(m, z)) * Q); z = forward(m, z); Q = Matrix(F.Q)
                for i in 1:3
                    acc[i] += log(abs(F.R[i, i]))
                end
            end
            emit("R6", "start 1 independent plain-QR exponents of the same 1e6-step window, max difference to lambda_fwd", maximum(abs.(acc ./ ns .- res.lambda_fwd)), string("plain QR: ", join(fmt.(acc ./ ns), ", ")), "INFO")
        end
    end
    ws = Vector{Vector{Float64}}()
    for is in 1:8
        x0 = SVector(clv_starts(StableRNG(VSEED + 520 + is), 1, 1, Float64)[1][1] / 10, 0.2, 0.5)
        res = propagate_clv_orbit(ClvGinelli(m; tau = 1, h = 1.0), x0, 3000, 1_000_000, 3000; k = 3, seed = VSEED + 530 + is, outputs = ())
        all(isfinite, res.lambda_fwd) && push!(ws, [mean(@view res.lambda[i, 1:1_000_000]) for i in 1:3])
    end
    for i in 1:3
        v = [w[i] for w in ws]
        emit("R6", string("mean over ", length(ws), " random starts of the 1e6-step window means, lambda_", i, " (standard error in the bound column)"), mean(v), string("printed ", target[i], " ; standard error ", fmt(std(v) / sqrt(length(v)))), "INFO")
    end
    emit("O-e", "Ginelli2013 pp.25 to 26 (minimum angle scales as 1/N in the FPU chain and the Henon chain): system not in the catalog, not run", 0, "reported only", "INFO")
end

# ---------------------------------------------------------------------------------------------------
# R1: wall-clock of the phases (forward reference = delivered discrete QR spectrum over the same handoffs; backward inferred by difference), checkpointed against store-all
timeit(f; reps = 3) = (f(); minimum(@elapsed(f()) for _ in 1:reps))
section("R1") do
    for T in (Float64, Double64)
        cases = [("standard map k = 2", clv_std(T), SVector{2, T}(T(1.3), T(0.7)), 1, 1.0, 2, 3000),
                 ("coupled maps k = 4", clv_cm(T), SVector{4, T}(map_to(T, SVector(1.1, 2.3, 0.4, 5.1))), 1, 1.0, 4, 3000),
                 ("F3 orbit C k = 4", clv_im(T), map_to(T, clv_orbit_C()), 20, CLV_H, 4, 300)]
        for (lab, map, x0, tau, h, k, n) in cases
            N = 3n; W0 = clv_frame(VSEED + 600, 1, dimension(map), k, T)
            tfwd = timeit(() -> propagate(map, [x0], N * tau; k = k, tau = tau, seed = 1, reducers = (QRSpectrum(HouseholderQR()),), W0 = [W0]))
            clv = ClvGinelli(map; tau = tau, h = h); clvc = ClvGinelli(map; tau = tau, h = h, storage = :checkpoint, block = max(1, N ÷ 10))
            ttot = timeit(() -> propagate_clv_orbit(clv, x0, n, n, n; k = k, seed = 1, W0 = W0))
            tV = timeit(() -> propagate_clv_orbit(clv, x0, n, n, n; k = k, seed = 1, W0 = W0, outputs = (:V,)))
            tck = timeit(() -> propagate_clv_orbit(clvc, x0, n, n, n; k = k, seed = 1, W0 = W0))
            tbwd = max(ttot - tfwd, 0.0)
            emit("R1", string(lab, " ", T, " delivered spectrum (forward reference), ", N, " handoffs, seconds"), tfwd, "-", "INFO")
            emit("R1", string(lab, " ", T, " module store-all no outputs, seconds"), ttot, "-", "INFO")
            emit("R1", string(lab, " ", T, " module store-all outputs (:V), seconds"), tV, "-", "INFO")
            emit("R1", string(lab, " ", T, " module checkpointed (block ", max(1, N ÷ 10), ") seconds"), tck, "-", "INFO")
            emit("R1", string(lab, " ", T, " forward / backward ratio (backward = store-all minus forward reference)"), tbwd > 0 ? tfwd / tbwd : Inf, "reported 2.4 (6 without the phase-space expansion), Ginelli2013 p.15 (V18)", "INFO")
            emit("R1", string(lab, " ", T, " checkpointed / store-all"), tck / ttot, "reported about 2 (V17)", "INFO")
        end
    end
end

# ---------------------------------------------------------------------------------------------------
# R2: angle distributions over 1e5 samples (minimum principal angle between the unstable and stable sets, and the pairwise angles)
cdf_slope(a, xs) = (ys = [count(<(x), a) / length(a) for x in xs]; ok = findall(>(10 / length(a)), ys); length(ok) < 3 ? NaN : (u = log.(xs[ok]); v = log.(ys[ok]); sum((u .- mean(u)) .* (v .- mean(v))) / sum((u .- mean(u)) .^ 2)))
section("R2") do
    NS = 100_000; data = Dict{String, Any}()
    cases = [("standard map K = 5", clv_std(Float64), SVector(1.3, 0.7), 1, 1.0, 2, (1:1, 2:2), 400),
             ("coupled maps d = 4", clv_cm(Float64), SVector(1.1, 2.3, 0.4, 5.1), 1, 1.0, 4, (1:2, 3:4), 400),
             ("F3 orbit C", clv_im(Float64), clv_orbit_C(), 20, CLV_H, 4, (1:1, 4:4), 5000)]
    for (lab, map, x0, tau, h, k, sets, nt) in cases
        clv = ClvGinelli(map; tau = tau, h = h)
        t = @elapsed res = propagate_clv_orbit(clv, x0, nt, NS, nt; k = k, seed = VSEED + 700, angle_sets = sets)
        th = vec(res.angles[1, :]); data[lab] = sort(th)
        emit("R2", string(lab, " minimum principal angle, samples, seconds"), length(th), string(@sprintf("%.1f", t), " s (", nt, " transient handoffs each)"), "INFO")
        emit("R2", string(lab, " minimum principal angle: min, median, max"), minimum(th), string(fmt(median(th)), ", ", fmt(maximum(th))), "INFO")
        for x in (1.0e-1, 1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5)
            emit("R2", string(lab, " fraction of samples with minimum principal angle below ", x), count(<(x), th) / length(th), "-", "INFO")
        end
        emit("R2", string(lab, " log-log slope of the empirical CDF of the minimum angle on 1e-4 .. 1e-1 (P(theta < x) ~ x^slope)"), cdf_slope(th, 10 .^ range(-4, -1; length = 13)), "-", "INFO")
        for i in 1:(k - 1), j in (i + 1):k
            a = acos.(clamp.(vec(res.cosines[i, j, :]), 0.0, 1.0))
            emit("R2", string(lab, " pair (", i, ",", j, ") angle: min, fraction below 1e-2, fraction below 1e-3"), minimum(a), string(fmt(count(<(1.0e-2), a) / length(a)), ", ", fmt(count(<(1.0e-3), a) / length(a))), "INFO")
        end
    end
    serialize(joinpath(SCR, "r2-angles.jls"), data)
end

# ---------------------------------------------------------------------------------------------------
# R3 index: the worst-conditioned vector index (error of each CLV against the 256-bit reference, F3 orbit C and standard-map starts, B_v with c_v = 1e4 as the package test)
section("R3idx") do
    for T in (Float64, Double64)
        cases = Any[(label = "std $j", map = clv_std(T), x0 = x0, tau = 1, h = 1.0, nf = 60, ns = 60, nb = 60) for (j, x0) in enumerate(clv_starts(StableRNG(VSEED + 800), 2, 20, T))]
        push!(cases, (label = "F3 orbit C", map = clv_im(T), x0 = map_to(T, clv_orbit_C()), tau = 20, h = CLV_H, nf = 60, ns = 100, nb = 60))
        wv = Dict{String, Vector{Float64}}(); wl = Dict{String, Vector{Float64}}()
        for (ic, c) in enumerate(cases)
            d = dimension(c.map); Dt = c.tau * c.h; W0 = clv_frame(VSEED + 810 + ic, 1, d, d, T); C0 = clv_rand_C(StableRNG(VSEED + 850 + ic), d, T)
            res = clv_go(c.map, c.x0, c.nf, c.ns, c.nb; tau = c.tau, h = c.h, W0 = W0, C0 = C0)
            pts = clv_orbit(c.map, c.x0, (c.nf + c.ns + c.nb) * c.tau)
            ref = clv_ref_run(clv_bigmap(c.map), pts, c.tau, c.nf, c.ns, c.nb, W0, C0, BigFloat)
            Bv = 1.0e4 * clv_eps(T) * (1 + clv_Gmax(res.R, c.ns)); key = startswith(c.label, "std") ? "standard map (20 starts)" : c.label
            ev = zeros(d); el = zeros(d)
            for n in 0:c.ns
                Vr = clv_refV(ref, n)
                for i in 1:d
                    ev[i] = max(ev[i], clv_delta(res.V[:, i, n + 1], Vr[:, i]) / Bv)
                    el[i] = max(el[i], Float64(abs(clv_big(res.lambda[i, n + 1]) * Dt - log(ref.Lambda[i, n + 1]))) / Bv)
                end
            end
            wv[key] = max.(get(wv, key, zeros(d)), ev); wl[key] = max.(get(wl, key, zeros(d)), el)
        end
        for key in sort(collect(keys(wv)))
            emit("R3", string(T, " ", key, " per vector index i: max delta(v_i) / B_v"), maximum(wv[key]), string(join(fmt.(wv[key]), ", "), " ; worst index ", argmax(wv[key])), "INFO")
            emit("R3", string(T, " ", key, " per index i: max |ln Lambda_i error| / B_v"), maximum(wl[key]), string(join(fmt.(wl[key]), ", "), " ; worst index ", argmax(wl[key])), "INFO")
        end
    end
end

# ---------------------------------------------------------------------------------------------------
# R4: the transient actually needed for 1e-10 (two seeds, same orbit window, distance between the CLVs) against 40 / (g_min Dt) (plan U3) and 25 / (g_min Dt) (spec); curves saved for the figure
section("R4") do
    function scan(label, map, x0, tau, h; singles, block, skip, ladder)
        d = dimension(map); lam = clv_pilot(map, x0, tau, h; ns = 1000); gaps = [lam[i] - lam[i + 1] for i in 1:(d - 1)]; skip && deleteat!(gaps, 2)
        g = minimum(gaps); Dt = tau * h; dist = Float64[]
        for m in ladder
            runs = [clv_go(map, x0, m, 4, m; tau = tau, h = h, W0 = clv_frame(VSEED + 900 + s, 1, d, d, Float64), C0 = clv_rand_C(StableRNG(VSEED + 920 + s), d, Float64), outputs = (:V,)) for s in 1:2]
            push!(dist, maximum(n -> clv_vdist(runs[1].V[:, :, n + 1], runs[2].V[:, :, n + 1]; singles = singles, block = block), 0:4))
        end
        bad = findall(>(1.0e-10), dist); need = isempty(bad) ? ladder[1] : (bad[end] == length(ladder) ? NaN : ladder[bad[end] + 1])
        p40 = 40 / (g * Dt); p25 = 25 / (g * Dt)
        emit("R4", string(label, " g_min Dt = ", fmt(g * Dt), ": n needed for 1e-10 (ladder resolution about 15 percent)"), need, string("40/(g Dt) = ", fmt(p40), " ; 25/(g Dt) = ", fmt(p25), " ; needed/25-prediction = ", fmt(need / p25), " ; needed/40-prediction = ", fmt(need / p40)), "INFO")
        return (label = label, ms = collect(ladder), dist = dist, need = need, p25 = p25, p40 = p40, g = g * Dt)
    end
    curves = Any[]
    ladder_small = unique(round.(Int, 10 .^ range(log10(2), log10(300); length = 40)))
    ladder_f3 = unique(round.(Int, 10 .^ range(log10(20), log10(7000); length = 45)))
    for (j, x0) in enumerate(clv_starts(StableRNG(VSEED + 3), 2, 4, Float64))
        push!(curves, scan("standard map start $j", clv_std(Float64), x0, 1, 1.0; singles = 1:2, block = nothing, skip = false, ladder = ladder_small))
    end
    for (j, x0) in enumerate(clv_starts(StableRNG(VSEED + 4), 4, 3, Float64))
        push!(curves, scan("coupled maps start $j", clv_cm(Float64), x0, 1, 1.0; singles = 1:4, block = nothing, skip = false, ladder = ladder_small))
    end
    for (j, x0) in enumerate(clv_chaotic_starts(5, 1))
        push!(curves, scan("F3 chaotic start $j", clv_im(Float64), x0, 20, CLV_H; singles = [1, 4], block = 2:3, skip = true, ladder = ladder_f3))
    end
    push!(curves, scan("F3 orbit C", clv_im(Float64), clv_orbit_C(), 20, CLV_H; singles = [1, 4], block = 2:3, skip = true, ladder = ladder_f3))
    serialize(joinpath(SCR, "r4-curves.jls"), curves)
end

# ---------------------------------------------------------------------------------------------------
# R5: orthonormalizer and tau dependence of the local exponents at common samples
section("R5") do
    function lamdiff(map, x0, tau_s, tau_b, h, nf, nsb; d = dimension(map))
        r = tau_b ÷ tau_s; W0 = clv_frame(VSEED + 1000, 1, d, d, Float64); C0 = clv_rand_C(StableRNG(VSEED + 1001), d, Float64)
        big = clv_go(map, x0, nf, nsb, nf; tau = tau_b, h = h, W0 = W0, C0 = C0)
        small = clv_go(map, x0, r * nf, r * nsb, r * nf; tau = tau_s, h = h, W0 = W0, C0 = C0)
        w = 0.0
        for m in 0:(nsb - 1), i in 1:d
            ssum = sum(small.lambda[i, r * m + l + 1] * tau_s * h for l in 0:(r - 1))
            w = max(w, abs(big.lambda[i, m + 1] * tau_b * h - ssum))
        end
        return w
    end
    emit("R5", "standard map: max |ln Lambda_i(tau = 5) - sum of five tau = 1 sub-interval values| over samples and i (n_f = n_b = 200)", lamdiff(clv_std(Float64), SVector(1.3, 0.7), 1, 5, 1.0, 200, 40), "reported (intrinsic: rounding only)", "INFO")
    emit("R5", "coupled maps: same, tau 1 against 5 (n_f = n_b = 200)", lamdiff(clv_cm(Float64), SVector(1.1, 2.3, 0.4, 5.1), 1, 5, 1.0, 200, 40), "reported", "INFO")
    z = clv_chaotic_starts(1, 1)[1]
    emit("R5", "F3 chaotic start: max |ln Lambda_i(tau = 20) - sum of two tau = 10 values| (n_f = n_b = 3000)", lamdiff(clv_im(Float64), z, 10, 20, CLV_H, 3000, 10), "reported", "INFO")
    for (lab, map, x0, tau, h, nf) in (("standard map", clv_std(Float64), SVector(1.3, 0.7), 1, 1.0, 200), ("F3 chaotic start", clv_im(Float64), z, 20, CLV_H, 3000))
        d = dimension(map); W0 = clv_frame(VSEED + 1010, 1, d, d, Float64); C0 = clv_rand_C(StableRNG(VSEED + 1011), d, Float64)
        a = clv_go(map, x0, nf, 20, nf; tau = tau, h = h, W0 = W0, C0 = C0, alg = HouseholderQR()); b = clv_go(map, x0, nf, 20, nf; tau = tau, h = h, W0 = W0, C0 = C0, alg = RGS())
        emit("R5", string(lab, ": max |lambda_i(n) Householder - lambda_i(n) RGS| over samples (n_f = n_b = ", nf, ")"), maximum(abs.(a.lambda .- b.lambda)), "reported (intrinsic: rounding only)", "INFO")
        emit("R5", string(lab, ": max delta CLV Householder against RGS"), maximum(n -> clv_vdist(a.V[:, :, n + 1], b.V[:, :, n + 1]; singles = (d == 4 ? [1, 4] : 1:d), block = (d == 4 ? (2:3) : nothing)), 0:20), "reported", "INFO")
    end
end

# ---------------------------------------------------------------------------------------------------
# THR: execution invariance, 1 against N threads (criterion 10 / I6 beyond the package test, which varies tasks): tasks and chunks in one process, then a one-thread child process
function thr_run(; ntasks, chunk)
    xs = clv_starts(StableRNG(VSEED + 300), 2, 1024, Float64)
    return propagate_clv(ClvGinelli(clv_std(Float64)), xs, 150, 200, 200; k = 2, seed = 7, chunk = chunk, ntasks = ntasks, outputs = (:V,))
end
thr_pack(r) = [(C = x.C, V = x.V, lambda = x.lambda, x = x.x, lf = x.lambda_fwd) for x in r.results]
same(a, b) = all(isequal(p.C, q.C) && isequal(p.V, q.V) && isequal(p.lambda, q.lambda) && isequal(p.x, q.x) && isequal(p.lf, q.lf) for (p, q) in zip(a, b))
section("THRchild") do
    serialize(joinpath(SCR, "thr-child.jls"), thr_pack(thr_run(ntasks = 1, chunk = 64)))
    say("# THRchild wrote the one-thread results")
end
section("THR") do
    nt = Threads.nthreads(); emit("AC10", "threads of this process", nt, ">= 2", pf(nt >= 2))
    r1 = thr_pack(thr_run(ntasks = 1, chunk = 64)); rN = thr_pack(thr_run(ntasks = nt, chunk = 1000)); rM = thr_pack(thr_run(ntasks = nt, chunk = 64))
    emit("AC10", "1024 orbits: ntasks 1 against ntasks N, bitwise (C, V, lambda, final state, forward exponents), chunk 64 against 1000", same(r1, rN) ? 0 : 1, "0 differences", pf(same(r1, rN)))
    emit("AC10", "1024 orbits: chunk 64 against 1000 with N tasks, bitwise", same(rM, rN) ? 0 : 1, "0 differences", pf(same(rM, rN)))
    cmd = addenv(`$(Base.julia_cmd()) --project=$(joinpath(ROOT, "verify")) -t 1 $(@__FILE__)`, "CLV_ONLY" => "THRchild", "CLV_LOG" => joinpath(SCR, "thr-child-log.txt"), "CLV_SCRATCH" => SCR)
    run(cmd)
    child = deserialize(joinpath(SCR, "thr-child.jls"))
    emit("AC10", "one-thread process (-t 1) against the N-thread process, bitwise over 1024 orbits", same(child, rN) ? 0 : 1, "0 differences", pf(same(child, rN)))
end

# ---------------------------------------------------------------------------------------------------
# FIG: convergence of the backward and forward transients (criterion 3 (a) and (b)), and the angle distributions (R2)
section("FIG") do
    curves = deserialize(joinpath(SCR, "r4-curves.jls"))
    fig = CM.Figure(size = (1400, 560))
    ax1 = CM.Axis(fig[1, 1]; xscale = log10, yscale = log10, xlabel = "transient n_f = n_b (handoffs)", ylabel = "max delta between two seeds of (frame, C)", title = "independence of the start (criterion 3 b, R4), 1e-10 level dashed", titlesize = 13)
    cols = Dict("standard map" => :royalblue, "coupled maps" => :seagreen, "F3 chaotic" => :crimson, "F3 orbit C" => :black)
    done = Set{String}()
    for c in curves
        key = first(k for k in keys(cols) if startswith(c.label, k)); ok = findall(>(0), c.dist)
        CM.lines!(ax1, c.ms[ok], max.(c.dist[ok], 1.0e-17); color = cols[key], label = key in done ? nothing : key); push!(done, key)
    end
    CM.hlines!(ax1, [1.0e-10]; color = :gray, linestyle = :dash); CM.axislegend(ax1; position = :lb, labelsize = 11)
    # panel 2: criterion 3 (a), column j of the first delivered sample against the backward transient m, fixture b, exact vectors, with the predicted slope
    ax2 = CM.Axis(fig[1, 2]; yscale = log10, xlabel = "backward transient m (handoffs)", ylabel = "delta(v_j computed, v_j exact)", title = "criterion 3 a: fixture b, columns j = 2 .. 4, slope -(lambda_(j-1) - lambda_j)", titlesize = 13)
    fx = clv_fixture(:b, Float64); d = fx.d; W0 = clv_frame(VSEED + 1100, 1, d, d, Float64); C0 = clv_rand_C(StableRNG(VSEED + 1101), d, Float64)
    errs = zeros(d, 40)
    for m in 2:40
        res = clv_go(fx.map, zero(SVector{d, Float64}), 100, 1, m; W0 = W0, C0 = C0)
        for j in 1:d
            errs[j, m] = clv_delta(res.V[:, j, 1], fx.V[:, j])
        end
    end
    for (j, col) in zip(2:d, (:royalblue, :seagreen, :crimson))
        ms = [m for m in 2:40 if errs[j, m] > 1.0e3 * clv_eps(Float64)]
        CM.scatterlines!(ax2, ms, errs[j, ms]; color = col, label = "column $j")
        sl = -Float64(fx.lam[j - 1] - fx.lam[j]); CM.lines!(ax2, ms, errs[j, ms[1]] .* exp.(sl .* (ms .- ms[1])); color = col, linestyle = :dash)
    end
    CM.axislegend(ax2; position = :rt, labelsize = 11)
    CM.save(joinpath(OUT, "clv-ginelli-convergence.png"), fig)
    if isfile(joinpath(SCR, "r2-angles.jls"))
        data = deserialize(joinpath(SCR, "r2-angles.jls")); fig2 = CM.Figure(size = (900, 560))
        ax = CM.Axis(fig2[1, 1]; xscale = log10, yscale = log10, xlabel = "theta (rad), minimum principal angle between unstable and stable sets", ylabel = "P(angle < theta), 1e5 samples", title = "angle distributions, R2", titlesize = 13)
        for (lab, col) in (("standard map K = 5", :royalblue), ("coupled maps d = 4", :seagreen), ("F3 orbit C", :crimson))
            haskey(data, lab) || continue
            a = data[lab]; n = length(a); idx = unique(round.(Int, 10 .^ range(0, log10(n); length = 400)))
            CM.lines!(ax, max.(a[idx], 1.0e-9), idx ./ n; color = col, label = lab)
        end
        CM.axislegend(ax; position = :lt, labelsize = 11); CM.save(joinpath(OUT, "clv-ginelli-angles.png"), fig2)
    end
    say("# FIG: wrote clv-ginelli-convergence.png and clv-ginelli-angles.png in the stage 04 output folder")
end
