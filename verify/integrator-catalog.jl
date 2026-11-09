# Stage 04 verification of module integrator-catalog (Chaoslab.Integrators, flows in Chaoslab.Systems; catalog rows S1 to S9, S11 to S13, P1 to P3, G1 to G3, R1).
# Run from the repo root, one process per group (the groups are independent; IC_LOG names the log, IC_APPEND=1 appends):
#   IC_ONLY=AC01-03   julia --project=verify -t 4 verify/integrator-catalog.jl     (also AC04-08, AC09-12, AC13-16, AC17-23 = the package test files at the FULL sizes, IC_FULL_RUN = true)
#   IC_ONLY=O6        (OrdinaryDiffEqSymplecticRK one-step maps and tables, criteria 3(c), 3(b4), 3(b5))
#   IC_ONLY=O7        (TaylorIntegration reference orbits R1 of F2, F4, T1 and a reference tangent)
#   IC_ONLY=O8        (REBOUND integrator_saba.c / integrator_leapfrog.c arrays, criterion 3(d), the ABA864 b4 verdict)
#   IC_ONLY=O9        (OrdinaryDiffEq ImplicitMidpoint against G1)
#   IC_ONLY=AC15D     (criterion 15(d): Gauss s = 4, Brouwer's law, standard implementation as negative control; IC_D_STARTS, IC_D_NMAX)
#   IC_ONLY=AC18T     (criterion 18: 1 against N CPU threads, bitwise; spawns a one-thread child)
#   IC_ONLY=MARGINS   (worst measured ratios of criteria 1, 9, 10, 11, 13, 14, 16, 17)
#   IC_ONLY=GPUHOOK   (D113: the device hook of the ensemble kernel exists in ext/ChaoslabCUDAExt.jl; no device code runs)
#   IC_ONLY=FIG       (integrator-catalog-order.png, integrator-catalog-brouwer.png)
# Not run here: the GPU part of criterion 18 (needs verify/gpu.jl, which writes data/gpu-benchmark.toml; owner decision), AC3(b2) (McLachlan1995 prints, text not in the repo, open point P15).
# Spec: workspaces/module-build/stages/01-spec/output/integrator-catalog-spec.md (criteria 1 to 23, Amendment D112). Seeds 20261012 + 41000 .. 41999 (ic_seed); verify-only seeds 41900 .. 41999 (VSEED).
# One printed line per measurement:  AC<n> | quantity | value | bound | PASS/FAIL/INFO
const ONLY = filter(!isempty, split(get(ENV, "IC_ONLY", ""), ","))
using Printf, Random, LinearAlgebra, Statistics, Serialization, Logging
using Test
const IC_FULL_RUN = true
using Chaoslab, StaticArrays, DoubleFloats, StableRNGs, ForwardDiff
import TaylorIntegration                      # loads the package extension ChaoslabTaylorExt (qualified import: `using` makes `inverse` ambiguous)
import OrdinaryDiffEqSymplecticRK as SRK
import OrdinaryDiffEqSDIRK as SD
import Pkg
if "FIG" in ONLY
    import CairoMakie as CM
end
if "GPUHOOK" in ONLY
    import CUDA                               # loads the extension ChaoslabCUDAExt; nothing runs on a device and nothing is written under data/
end

const ROOT = dirname(@__DIR__)
const OUT = joinpath(ROOT, "workspaces", "module-build", "stages", "04-verify", "output")
const SCR = get(ENV, "IC_SCRATCH", joinpath(tempdir(), "ic-scratch"))
mkpath(SCR)
const LOG = open(get(ENV, "IC_LOG", joinpath(OUT, "integrator-catalog-verify-log.txt")), get(ENV, "IC_APPEND", "") == "1" ? "a" : "w")

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

# @info of the tests -> the log; the AC15 payload is kept for the figure
struct TeeLog <: Logging.AbstractLogger end
Logging.min_enabled_level(::TeeLog) = Logging.Info
Logging.shouldlog(::TeeLog, args...) = true
Logging.catch_exceptions(::TeeLog) = true
function Logging.handle_message(::TeeLog, level, message, _mod, group, id, file, line; kwargs...)
    kv = join(("$k = $(first(repr(v), 700))" for (k, v) in kwargs), "; ")
    say(string("# ", level, " ", message, " | ", kv))
    if startswith(string(message), "AC15 reported")
        serialize(joinpath(SCR, "ac15.jls"), Dict(string(k) => v for (k, v) in kwargs))
    end
end
global_logger(TeeLog())

function section(f, name)
    isempty(ONLY) || name in ONLY || return
    t0 = time()
    try
        f()
    catch e
        say(string(name, " | ERROR | ", sprint(showerror, e, catch_backtrace())[1:min(end, 1500)], " | - | FAIL"))
    end
    say(@sprintf("# section %s took %.1f s", name, time() - t0))
end

say("# integrator-catalog verification, julia $(VERSION), threads $(Threads.nthreads()), base seed 20261012 + 41000, IC_FULL_RUN = true, time $(round(time())) (unix)")
say("# timings in this log were taken on a machine that may be loaded by other processes")
for (_, info) in Pkg.dependencies()
    info.name in ("OrdinaryDiffEqSymplecticRK", "OrdinaryDiffEqSDIRK", "TaylorIntegration", "OrdinaryDiffEqCore") && say("# package $(info.name) $(info.version)")
end
say("# extension ChaoslabTaylorExt loaded: $(Base.get_extension(Chaoslab, :ChaoslabTaylorExt) !== nothing)")

# ---------------------------------------------------------------------------------------------------
# The package suite at full size (IC_FULL_RUN = true: AC11 all 41 methods, AC12 and AC13 and AC15 at spec size, AC17 200 orbits, AC18 1e5 orbits).
const TD = joinpath(ROOT, "test", "integrator-catalog")
include(joinpath(TD, "ic_helpers.jl"))
say("# IC_FULL = $IC_FULL")

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
        nf = nfail(ch); c = Test.get_test_counts(ch)   # per-testset pass counts are not available (rolled into the file total)
        emit(ac, string("testset ", ch.description, " (", nf, " fail or error, ", c.cumulative_broken, " broken, ", @sprintf("%.1f", ch.time_end - ch.time_start), " s)"), nf, "0", pf(nf == 0))
        report_ts(ac, ch; top = false)
    end
end
# Definitions (function, const, assignment) of an earlier test file, without running its testsets.
function include_defs(file)
    for e in Meta.parseall(read(file, String)).args
        e isa Expr && e.head in (:function, :const, :(=)) && Core.eval(Main, e)
    end
end
function suite_section(name, files; defs = String[])
    section(name) do
        for f in defs
            include_defs(joinpath(TD, f))
        end
        Test.@testset "outer $name" begin
            ts = Test.@testset "$name" begin
                for f in files
                    include(joinpath(TD, f))
                end
            end
            report_ts(name, ts)
        end
    end
end
suite_section("AC01-03", ["ic_flows.jl", "ic_ac01_03.jl"])
suite_section("AC04-08", ["ic_ac04_08.jl"]; defs = ["ic_flows.jl", "ic_ac01_03.jl"])
suite_section("AC09-12", ["ic_ac09_12.jl"]; defs = ["ic_flows.jl", "ic_ac01_03.jl", "ic_ac04_08.jl"])
suite_section("AC13-16", ["ic_ac13_16.jl"]; defs = ["ic_flows.jl", "ic_ac01_03.jl", "ic_ac04_08.jl", "ic_ac09_12.jl"])
suite_section("AC17-23", ["ic_ac17_23.jl"]; defs = ["ic_flows.jl", "ic_ac01_03.jl", "ic_ac04_08.jl", "ic_ac09_12.jl", "ic_ac13_16.jl"])

# ---------------------------------------------------------------------------------------------------
# From here: sections that are not the package testsets.
const VSEED = 20261012 + 41900
need_defs() = for f in ("ic_flows.jl", "ic_ac01_03.jl", "ic_ac04_08.jl", "ic_ac09_12.jl", "ic_ac13_16.jl", "ic_ac17_23.jl")
    include_defs(joinpath(TD, f))
end
any(x -> x in ("O6", "O7", "O8", "O9", "AC15D", "AC18T", "AC18child", "FIG", "MARGINS"), ONLY) && need_defs()

# ------------------------------------------------------------------------------------------------ O6 (OrdinaryDiffEqSymplecticRK)
# The package stage is (drift b_i, kick a_i) with the drift first (symplectic_perform_step.jl, generic tableau): our table (c, d) = (b, a).
# `pkg_map` runs one step of the package solver on F1 or F2; `ops` lists the drift/kick operations of an ours-method, `swap_ops` exchanges the roles of drift and kick.
const W1 = IC_F1W
function pkg_map(alg, kind::Symbol, z, dt)
    n = length(z) ÷ 2
    f1 = kind === :f1 ? (v, u, p, t) -> -(collect(W1[1:n]) .^ 2) .* u : (v, u, p, t) -> -u ./ (u[1]^2 + u[2]^2)^1.5
    prob = SRK.DynamicalODEProblem(f1, (v, u, p, t) -> v, collect(z[(n + 1):end]), collect(z[1:n]), (0.0, dt))
    sol = SRK.solve(prob, alg; dt = dt, adaptive = false, save_everystep = false)
    u = sol.u[end]
    return SVector{2n}(vcat(u.x[2], u.x[1]))
end
function ops(m)
    c, d = coefficients(m); o = Tuple{Symbol, Float64}[]
    for i in eachindex(c)
        iszero(c[i]) || push!(o, (:D, Float64(c[i]))); iszero(d[i]) || push!(o, (:K, Float64(d[i])))
    end
    return o
end
function from_ops(o)
    c = Float64[]; d = Float64[]; i = 1
    while i <= length(o)
        if o[i][1] === :D
            if i < length(o) && o[i + 1][1] === :K
                push!(c, o[i][2]); push!(d, o[i + 1][2]); i += 2
            else
                push!(c, o[i][2]); push!(d, 0.0); i += 1
            end
        else
            push!(c, 0.0); push!(d, o[i][2]); i += 1
        end
    end
    return SplittingMethod(c, d; order = 1, palindromic = false, source = "swap", name = "swapped")
end
swap_ops(m) = from_ops([(t === :D ? :K : :D, w) for (t, w) in ops(m)])

section("O6") do
    dcoef = Dict{String, Vector{Tuple{String, Float64}}}()
    # coefficients: the package tableaux in BigFloat against ours (operation lists, drift/kick roles as is and exchanged); computed first: D113 (F3) the bound of the one-step maps scales with it
    setprecision(BigFloat, 256) do
        for (nm, lab, mk) in (("Ruth3", "S4", () -> Ruth3(BigFloat)), ("CandyRoz4", "S5", () -> ForestRuth(BigFloat)), ("CalvoSanz4", "S9", () -> CalvoSanzSerna4(BigFloat)),
                              ("Yoshida6", "S11-6A", () -> Yoshida(BigFloat, 6, :A)), ("Yoshida6", "S12 s7odr6 (20 digits)", () -> KahanLi(BigFloat, :s7odr6)), ("KahanLi6", "S12-6", () -> KahanLi(BigFloat, 6)), ("KahanLi8", "S12-8", () -> KahanLi(BigFloat, 8)),
                              ("SofSpa10", "S13", () -> SofroniouSpaletta10(BigFloat)))
            tab = getfield(SRK, Symbol(nm, "ConstantCache"))(BigFloat, BigFloat)
            pops = Tuple{Symbol, BigFloat}[]
            for i in eachindex(tab.a)
                iszero(tab.b[i]) || push!(pops, (:D, BigFloat(tab.b[i]))); iszero(tab.a[i]) || push!(pops, (:K, BigFloat(tab.a[i])))
            end
            m = mk(); c, d = coefficients(m); mops = Tuple{Symbol, BigFloat}[]
            for i in eachindex(c)
                iszero(c[i]) || push!(mops, (:D, c[i])); iszero(d[i]) || push!(mops, (:K, d[i]))
            end
            sw = [(t === :D ? :K : :D, w) for (t, w) in mops]
            cmp(a, b) = length(a) == length(b) && all(x[1] === y[1] for (x, y) in zip(a, b)) ? Float64(maximum(abs(x[2] - y[2]) for (x, y) in zip(a, b))) : NaN
            e_as = cmp(pops, mops); e_sw = cmp(pops, sw)
            bound = nm == "SofSpa10" ? 1e-26 : nm == "CalvoSanz4" ? 1e-14 : 1e-13
            best = min(isnan(e_as) ? Inf : e_as, isnan(e_sw) ? Inf : e_sw)
            push!(get!(dcoef, nm, Tuple{String, Float64}[]), (lab, best))
            emit("3(c)", string("O6 table ", nm, " (BigFloat tableau) vs ", lab, ": max |coefficient difference| as is / roles exchanged (", length(pops), " vs ", length(mops), " operations)"),
                 string(fmt(e_as), " / ", fmt(e_sw)), string(bound), best <= bound ? (isnan(e_as) || e_as > bound ? "PASS (roles exchanged)" : "PASS") : "FAIL")
        end
    end
    rng = StableRNG(VSEED + 1)
    npts = 1000
    pts = Dict(:f1 => [ic_draw(rng, :f1_3) for _ in 1:npts], :f2 => [ic_draw(rng, :f2) for _ in 1:npts])
    dts = Dict(:f1 => 0.1, :f2 => 0.02)
    flows = Dict(:f1 => HarmonicOscillator(collect(W1)), :f2 => Kepler())
    rows = (  # (package solver name, our candidates (label, constructor))
        ("SymplecticEuler", () -> SRK.SymplecticEuler(), [("S1 DK", SymplecticEuler(:DK)), ("S1 KD", SymplecticEuler(:KD))]),
        ("VelocityVerlet", () -> SRK.VelocityVerlet(), [("S2 KDK", Leapfrog(:KDK)), ("S2 DKD", Leapfrog(:DKD))]),
        ("VerletLeapfrog", () -> SRK.VerletLeapfrog(), [("S2 KDK", Leapfrog(:KDK)), ("S2 DKD", Leapfrog(:DKD))]),
        ("LeapfrogDriftKickDrift", () -> SRK.LeapfrogDriftKickDrift(), [("S2 DKD", Leapfrog(:DKD)), ("S2 KDK", Leapfrog(:KDK))]),
        ("McAte2", () -> SRK.McAte2(), [("S3 McLachlan1995", McLachlan2(Float64))]),
        ("Ruth3", () -> SRK.Ruth3(), [("S4", Ruth3(Float64))]),
        ("CandyRoz4", () -> SRK.CandyRoz4(), [("S5 Forest-Ruth", ForestRuth(Float64))]),
        ("CalvoSanz4", () -> SRK.CalvoSanz4(), [("S9", CalvoSanzSerna4(Float64))]),
        ("Yoshida6", () -> SRK.Yoshida6(), [("S11-6A", Yoshida(6, :A)), ("S12 s7odr6", KahanLi(:s7odr6))]),
        ("KahanLi6", () -> SRK.KahanLi6(), [("S12-6 s9odr6a", KahanLi(6))]),
        ("KahanLi8", () -> SRK.KahanLi8(), [("S12-8 s17odr8a", KahanLi(8))]),
        ("SofSpa10", () -> SRK.SofSpa10(), [("S13", SofroniouSpaletta10(Float64))]))
    for (nm, mk, cands) in rows
        alg = mk()
        for kind in (:f1, :f2)
            dt = dts[kind]; fl = flows[kind]
            pk = [pkg_map(alg, kind, p, dt) for p in pts[kind]]
            best = (Inf, "", "")
            for (lab, m) in cands
                for (variant, mm) in (("as is", m), ("roles exchanged", swap_ops(m)))
                    mx = maximum(maximum(abs, collect(advance(mm, fl, SVector{length(p)}(p), dt)) .- collect(q)) for (p, q) in zip(pts[kind], pk))
                    mx < best[1] && (best = (mx, lab, variant))
                    emit("3(c)", string("O6 ", nm, " vs ", lab, " (", variant, "), ", kind === :f1 ? "F1 n=3" : "F2 e in {0.1,0.5,0.9}", ", ", npts, " points, max |difference|"), mx, "max(16 eps, 1e-13)", "INFO")
                end
            end
            # D113 (F3): the printed tableau carries 14 to 26 digits. A one-step map is first order in the coefficients; measured, the map difference is 3.5 times the largest
            # coefficient difference (Yoshida6 on F2: 1.03e-13 against 3.0e-14). The bound is max(16 eps, 1e-13, 8 x the largest coefficient difference of the matched table).
            dc = something(findfirst(t -> startswith(best[2], t[1]) || startswith(t[1], best[2]), get(dcoef, nm, Tuple{String, Float64}[])), 0)
            dcv = dc == 0 ? 0.0 : dcoef[nm][dc][2]
            bound = max(16 * eps(), 1e-13, 8 * dcv)
            bstr = "max(16 eps, 1e-13, 8 dc), dc = " * fmt(dcv)
            verdict = nm == "McAte2" ? "INFO (different method from S3, Inconsistencies (f))" : best[1] <= bound ? (best[3] == "as is" ? "PASS" : "PASS (roles of drift and kick exchanged)") : "FAIL"
            emit("3(c)", string("O6 ", nm, " best match ", best[2], " (", best[3], "), ", kind === :f1 ? "F1" : "F2", ", max |difference|"), best[1], bstr, verdict)
        end
    end
    # the McAte rows are not catalog rows: only listed; McAte42 is "Broken" in the source and not used
    emit("3(c)", "O6 McAte3, McAte4, McAte5, McAte8 (McLachlan-type, not catalog rows)", "not compared", "-", "INFO")
end

# ------------------------------------------------------------------------------------------------ O7 (TaylorIntegration, R1)
section("O7") do
    ext = Base.get_extension(Chaoslab, :ChaoslabTaylorExt)
    emit("O7", "extension ChaoslabTaylorExt loaded", ext !== nothing, "true", pf(ext !== nothing))
    setprecision(BigFloat, IC_BITS) do
        for (kind, z0, T) in ((:f4, IC_START[:f4], 8), (:f2, IC_START[:f2], 8), (:t1, IC_START[:t1], 8))
            flow = ic_flow(BigFloat, kind); zb = ic_to(BigFloat, z0)
            t1 = @elapsed r1 = taylor_reference(flow, zb, BigFloat(T))
            g = Gauss(BigFloat, 8; maxiter = 500, tol = 1e-40)
            e64 = ic_err(advance_n(g, flow, zb, BigFloat(1) / 64, 64 * T), r1)
            e128 = ic_err(advance_n(g, flow, zb, BigFloat(1) / 128, 128 * T), r1)
            dH = Float64(abs(ic_H(kind, collect(r1)) - ic_H(kind, collect(zb))))
            emit("O7", string("R1 orbit of ", kind, " T = 8 (", @sprintf("%.1f", t1), " s): max difference to Gauss s = 8 at tau = 1/64 and 1/128"), string(fmt(e64), " / ", fmt(e128)), "1e-30 at 1/128", pf(e128 <= 1e-30))
            emit("O7", string("R1 orbit of ", kind, " T = 8: |H(T) - H(0)|"), dH, "2^-180 = 6.5e-55", pf(dH <= 2.0^-180))
            kind === :f2 && emit("O7", "R1 orbit of f2 T = 8: |L(T) - L(0)|", Float64(abs(ic_L(r1) - ic_L(zb))), "2^-180", pf(Float64(abs(ic_L(r1) - ic_L(zb))) <= 2.0^-180))
        end
        # reference tangent: central difference of R1 orbits (step 2^-60) against the tangent of S12-8 (Float64, tau = 1/32) and S12-6 (Double64, tau = 1/128)
        for (kind, z0, T) in ((:f2, IC_START[:f2], 8), (:f4, IC_START[:f4], 8), (:t1, IC_START[:t1], 8))
            D = length(z0); flowb = ic_flow(BigFloat, kind); h = ldexp(BigFloat(1), -60)
            Wref = zeros(BigFloat, D, D)
            for j in 1:D
                zp = ic_to(BigFloat, z0); zm = ic_to(BigFloat, z0)
                zp = setindex(zp, zp[j] + h, j); zm = setindex(zm, zm[j] - h, j)
                Wref[:, j] = (collect(taylor_reference(flowb, zp, BigFloat(T))) .- collect(taylor_reference(flowb, zm, BigFloat(T)))) ./ (2h)
            end
            if kind !== :t1
                flow = ic_flow(Float64, kind); m = KahanLi(Float64, 8); tau = 1 / 32; W = SMatrix{D, D, Float64}(I); z = ic_to(Float64, z0)
                for _ in 1:round(Int, T / tau)
                    z, W = tangent_advance(m, flow, z, W, tau)
                end
                rel = Float64(maximum(abs(Float64(Wref[i, j]) - W[i, j]) for i in 1:D, j in 1:D) / maximum(abs.(Wref)))
                emit("O7", string("R1 tangent (central difference of orbits, 2^-60) of ", kind, " T = 8 vs tangent of S12-8 (Float64, tau = 1/32), relative max difference"), rel, "2^-30 = 9.3e-10", pf(rel <= 2.0^-30))
            else
                emit("O7", "R1 tangent of t1: no tangent step for the nonseparable T1 outside the extended space of G3 (open point P11)", "not compared", "-", "INFO")
            end
        end
    end
end

# ------------------------------------------------------------------------------------------------ O8 (REBOUND arrays)
section("O8") do
    srcdir = joinpath(SCR, "rebound-src")
    saba = read(joinpath(srcdir, "integrator_saba.c"), String); leap = read(joinpath(srcdir, "integrator_leapfrog.c"), String)
    say("# O8 REBOUND files read (not executed): integrator_saba.c $(length(saba)) bytes, integrator_leapfrog.c $(length(leap)) bytes; raw.githubusercontent.com/hannorein/rebound/main/src (HEAD fd484cd6f6ed3d538b8e9606cdc2ddf251e646e9, latest tag 5.2.2 c34134716cd484ed9cdb9f7d54ef8d859bdf547d, fetched 2026-11-06)")
    function arrays(name)
        m = match(Regex("$name\\[10\\]\\[5\\] = \\{(.*?)\\n\\}\\s*;", "s"), saba); body = m.captures[1]
        rows = Vector{Vector{String}}()
        for r in eachmatch(r"\{([^{}]*)\}", body)
            push!(rows, [strip(String(x)) for x in split(replace(replace(r.captures[1], r"//[^\n]*" => ""), "\n" => " "), ",") if !isempty(strip(String(x)))])
        end
        return rows
    end
    rc = arrays("reb_saba_c"); rd = arrays("reb_saba_d")
    emit("3(d)", "REBOUND arrays parsed: rows of reb_saba_c / reb_saba_d", string(length(rc), " / ", length(rd)), "10 / 10", pf(length(rc) == 10 && length(rd) == 10))
    names = ["SABA1", "SABA2", "SABA3", "SABA4", "ABA(10,4)", "ABA(8,6,4)", "ABA(10,6,4)", "ABAH(8,4,4)", "ABAH(8,6,4)", "ABAH(10,6,4)"]
    setprecision(BigFloat, 400) do
        P(s) = parse(BigFloat, s)
        # SABA_n: REBOUND stores the first halves of the drifts (n + 1 entries, the half rounded up) and of the kicks (n entries, the half rounded up)
        for n in 1:4
            m = SABA(BigFloat, n); c, d = coefficients(m); k = n + 1
            hc = (k + 1) ÷ 2; hd = (n + 1) ÷ 2
            ec = maximum(abs(c[i] - P(rc[n][i])) for i in 1:hc); ed = maximum(abs(d[i] - P(rd[n][i])) for i in 1:hd)
            e = Float64(max(ec, ed))
            emit("3(d)", string("REBOUND ", names[n], " drifts and kicks (first halves) vs SABA", n, " (BigFloat), max |difference|"), e, "1e-30", pf(e <= 1e-30))
        end
        # ABA864 and ABAH864
        a864 = [P(s) for s in rc[6]]; b864 = [P(s) for s in rd[6]]
        m = ABA864(BigFloat); c, d = coefficients(m)
        ea = Float64(maximum(abs(c[i] - a864[i]) for i in 1:4)); eb = Float64(maximum(abs(d[i] - b864[i]) for i in 1:4))
        emit("3(d)", "REBOUND ABA(8,6,4) a_1..a_4 vs ours ABA864 (BigFloat), max |difference|", ea, "1e-30", pf(ea <= 1e-30))
        emit("3(d)", "REBOUND ABA(8,6,4) b_1..b_4 vs ours ABA864 (BigFloat, b4 recomputed), max |difference|", eb, "1e-30", pf(eb <= 1e-30))
        emit("3(d)", "  of which the b4 difference (ours recomputed minus REBOUND print)", Float64(d[4] - b864[4]), "-", "INFO")
        mH = ABAH864(BigFloat); cH, dH = coefficients(mH); aH = [P(s) for s in rc[9]]; bH = [P(s) for s in rd[9]]
        ea = Float64(maximum(abs(cH[i] - aH[i]) for i in 1:5)); eb = Float64(maximum(abs(dH[i] - bH[i]) for i in 1:4))
        emit("3(d)", "REBOUND ABAH(8,6,4) a_1..a_5 and b_1..b_4 vs ours ABAH864, max |difference|", max(ea, eb), "1e-30", pf(max(ea, eb) <= 1e-30))
        # the printed sums
        for (nm, row, a, b, kb) in (("ABA(8,6,4)", 6, a864, b864, 4), ("ABAH(8,6,4)", 9, aH, bH, 4))
            # layout: a_1 .. a_k [a_(k+1)] b_k ... : ABA864 a1..a4 (sum of all drifts 2(a1+a2+a3)+2 a4 = 1), b: 2(b1+b2+b3)+b4 = 1
            sa = nm == "ABAH(8,6,4)" ? 2 * (sum(a[1:4])) + a[5] : 2 * sum(a[1:3]) + 2 * a[4]
            sb = nm == "ABAH(8,6,4)" ? 2 * sum(b[1:4]) : 2 * sum(b[1:3]) + b[4]
            emit("3(d)", string("REBOUND ", nm, " printed digits: sum of drifts - 1 / sum of kicks - 1"), string(fmt(Float64(sa - 1)), " / ", fmt(Float64(sb - 1))), "-", "INFO")
        end
        # ABA864 b4 verdict: the kick positions and the quadrature conditions of order eps tau^8 (moments 0, 2, 4, 6 about the midpoint); solve for b1..b4 from the printed a
        function solve_b(a, kind)
            # kick positions: ABA864 a1 b1 a2 b2 a3 b3 a4 b4 a4 ...: s_i = a1 + ... + a_i; the mirror points 1 - s_i; ABAH: b4 at s4 and again at s4 + a5
            s = cumsum(a[1:4])
            xs = kind === :aba ? vcat(s[1:3], s[4], reverse(1 .- s[1:3])) : vcat(s[1:4], reverse(1 .- s[1:4]))
            # unknown weights (b1, b2, b3, b4): symmetric pairs (and the central b4 once for ABA, as a pair for ABAH)
            mom(k) = (BigFloat(1) / 2)^(k) / (k + 1)           # int over [-1/2, 1/2] of y^k dy = 2 (1/2)^(k+1)/(k+1) = (1/2)^k/(k+1)
            A = zeros(BigFloat, 4, 4); r = zeros(BigFloat, 4)
            for (row, k) in enumerate((0, 2, 4, 6))
                y(x) = x - BigFloat(1) / 2
                for j in 1:3
                    A[row, j] = y(xs[j])^k + y(xs[end + 1 - j])^k
                end
                A[row, 4] = kind === :aba ? y(xs[4])^k : y(xs[4])^k + y(xs[5])^k
                r[row] = mom(k)
            end
            return A \ r
        end
        for (nm, a, b, kind) in (("ABA864", a864, b864, :aba), ("ABAH864", aH, bH, :abah))
            bs = solve_b(a, kind)
            emit("3(d)", string(nm, " b_1..b_4 solved from the printed a_i and the quadrature moments k = 0, 2, 4, 6 (eps tau^8 conditions): max |solved - printed|"), Float64(maximum(abs.(bs .- b))), "-", "INFO")
            for i in 1:4
                emit("3(d)", string(nm, " b", i, ": solved - printed (the print of REBOUND = the print of Blanes2013 Table 3 / 4)"), Float64(bs[i] - b[i]), "-", "INFO")
            end
            if nm == "ABA864"
                brec = 1 - 2 * (b[1] + b[2] + b[3])
                emit("3(d)", "ABA864 b4: recomputed 1 - 2(b1+b2+b3) minus solved (stage 03 choice)", Float64(brec - bs[4]), "<= 1e-36 confirms", pf(abs(brec - bs[4]) <= 1e-36))
                emit("3(d)", "ABA864 b4: printed minus solved", Float64(b[4] - bs[4]), "-", "INFO")
                # moments k = 0..7 of the kick weights at the kick positions, with the printed b, the recomputed b4 and the solved b
                s3 = cumsum(a[1:4]); xs7 = vcat(s3[1:3], s3[4], reverse(1 .- s3[1:3]))
                resid(bb) = Float64(maximum(abs(sum(bb[j] * (xs7[j] - BigFloat(1) / 2)^k for j in 1:7) - (iseven(k) ? (BigFloat(1) / 2)^k / (k + 1) : zero(BigFloat))) for k in 0:7))
                full(bb) = [bb[1], bb[2], bb[3], bb[4], bb[3], bb[2], bb[1]]
                emit("3(d)", "ABA864 first-order conditions (moments k = 0..7 about the midpoint): max residual with the printed b", resid(full(b864)), "-", "INFO")
                emit("3(d)", "ABA864 first-order conditions: max residual with b4 recomputed (stage 03)", resid(full([b864[1], b864[2], b864[3], brec])), "-", "INFO")
                emit("3(d)", "ABA864 first-order conditions: max residual with the solved b", resid(full(bs)), "-", "INFO")
            end
        end
        # corrector coefficients and the ABA(10,4) entry: sums only (ABA104 is not implemented in this module)
        cc = [P(String(m.captures[1])) for m in eachmatch(r"([0-9]\.[0-9]+),\s*//\s*SABAC", saba)]
        for n in 1:min(4, length(cc))
            e = Float64(abs(Chaoslab.Integrators._corrector_c(BigFloat, n, false) - cc[n]))
            emit("3(d)", string("REBOUND SABAC", n, " corrector coefficient vs ours (Laskar-Robutel Table II), |difference|"), e, "1e-30 (print has 34 to 43 digits)", pf(e <= 1e-30))
        end
        # leapfrog 4 (Yoshida): drift a, kick 2a, drift 1/2 - a, kick 1 - 4a, drift 1/2 - a, kick 2a, drift a
        la = P(match(r"lf4_a = ([0-9.]+);", leap).captures[1])
        fr = ForestRuth(BigFloat); c, d = coefficients(fr)
        # ours: c = (theta/2, (1-theta)/2, (1-theta)/2, theta/2), d = (theta, 1 - 2 theta, theta, 0), after merging adjacent drifts
        e1 = Float64(abs(c[1] - la)); e2 = Float64(maximum(abs.([c[2] - (BigFloat(1) / 2 - la), d[1] - 2la, d[2] - (1 - 4la)])))
        emit("3(d)", "REBOUND leapfrog order 4 (lf4_a, 30 digits) vs ours Forest-Ruth drift c1, |difference|", e1, "1e-30", pf(e1 <= 1e-30))
        emit("3(d)", "REBOUND leapfrog order 4: derived entries (1/2 - a, 2a, 1 - 4a of the C code) vs ours, max |difference| (the print rounds a at 5e-31, 1 - 4a carries 4 times that)", e2, "5e-30", pf(e2 <= 5e-30))
        theta = 1 / (2 - 2^(BigFloat(1) / 3))
        emit("3(d)", "REBOUND lf4_a minus theta/2, theta = 1/(2 - 2^(1/3)) (Yoshida1990 Eq. (2.11))", Float64(la - theta / 2), "1e-30", pf(abs(la - theta / 2) <= 1e-30))
        emit("3(d)", "REBOUND leapfrog orders 6 and 8 (Blanes and Casas 2016 arrays, not catalog rows): not compared", "-", "-", "INFO")
    end
end

# ------------------------------------------------------------------------------------------------ O9 (ImplicitMidpoint of OrdinaryDiffEq)
section("O9") do
    for (kind, z0, h, nst) in ((:f1_3, IC_START[:f1_3], 0.1, 100), (:f4, IC_START[:f4], 0.1, 100), (:f4, [2.5, 0.0], 0.05, 200))
        flow = ic_flow(Float64, kind)
        f = kind === :f1_3 ? (u, p, t) -> vcat(u[4:6], -(collect(W1) .^ 2) .* u[1:3]) : (u, p, t) -> [u[2], -sin(u[1])]
        prob = SD.ODEProblem(f, collect(z0), (0.0, h * nst))
        sol = SD.solve(prob, SD.ImplicitMidpoint(); dt = h, adaptive = false, reltol = 1e-14, abstol = 1e-14, save_everystep = true)
        ours = [advance_n(ImplicitMidpoint(Float64), flow, ic_to(Float64, z0), h, k) for k in 0:nst]
        # same ImplicitMidpoint solve compared at every step
        mx1 = maximum(maximum(abs, collect(ours[k + 1]) .- sol.u[k + 1]) for k in 0:nst)
        # and each against the 256-bit midpoint solution
        ref = setprecision(BigFloat, IC_BITS) do
            Vector(advance_n(ImplicitMidpoint(BigFloat), ic_flow(BigFloat, kind), ic_to(BigFloat, z0), BigFloat(h), nst))
        end
        eo = ic_err(ours[end], ref); ep = Float64(maximum(abs(sol.u[end][i] - ref[i]) for i in eachindex(ref)))
        emit("O9", string("ImplicitMidpoint vs G1, ", kind, ", h = ", h, ", ", nst, " steps, solver tolerance 1e-14: max over steps |difference|"), mx1, "1e-12", pf(mx1 <= 1e-12))
        emit("O9", string("  end state against the 256-bit implicit midpoint: ours / OrdinaryDiffEq"), string(fmt(eo), " / ", fmt(ep)), "-", "INFO")
    end
end

# ------------------------------------------------------------------------------------------------ AC15(d): Gauss s = 4, Brouwer
section("AC15D") do
    nst = parse(Int, get(ENV, "IC_D_STARTS", "200")); nmax = parse(Int, get(ENV, "IC_D_NMAX", "3000000"))
    ns = filter(<=(nmax), [1000, 3000, 10000, 30000, 100000, 300000, 1000000, 3000000])
    rng = StableRNG(ic_seed(15, 2)); tau = 2 / 240; flow = HenonHeilesFlow()
    starts = Vector{Vector{Float64}}()
    for _ in 1:nst
        q1 = 0.0 + (2rand(rng) - 1) * 5e-3; q2 = 0.3 + (2rand(rng) - 1) * 5e-3; p2 = 0.2 + (2rand(rng) - 1) * 5e-3
        V = (q1^2 + q2^2) / 2 + q1^2 * q2 - q2^3 / 3
        push!(starts, [q1, q2, sqrt(2 * (1 / 8 - V - p2^2 / 2)), p2])
    end
    results = Dict{String, Matrix{Float64}}(); nbad = Dict{String, Int}()
    for (label, mk, comp) in (("Hairer2008 implementation (iteration to convergence, split coefficients, compensated)", () -> Gauss(Float64, 4), true),
                              ("standard implementation (fixed threshold, rounded coefficients)", () -> Gauss(Float64, 4; standard = true), false),
                              ("diagnostic: Hairer2008 iteration and coefficients, plain summation", () -> Gauss(Float64, 4), false))
        E = fill(NaN, length(ns), nst); lk = ReentrantLock()
        t0 = time()
        Threads.@threads for j in 1:nst
            m = mk(); zs = advance_to(m, flow, ic_to(Float64, starts[j]), tau, ns; compensated = comp)
            setprecision(BigFloat, IC_BITS) do
                H0 = ic_H(:t2, ic_to(BigFloat, starts[j]))
                for (i, z) in enumerate(zs)
                    all(isfinite, z) && (E[i, j] = Float64(ic_H(:t2, ic_to(BigFloat, collect(z))) - H0))
                end
            end
        end
        good = [j for j in 1:nst if all(isfinite, E[:, j])]; nbad[label] = nst - length(good)
        Eg = E[:, good]; rms = [sqrt(sum(abs2, Eg[i, :]) / length(good)) for i in axes(Eg, 1)]
        idx = findall(>=(10000), ns)
        slope = length(idx) >= 2 ? (log(rms[idx[end]]) - log(rms[idx[1]])) / (log(ns[idx[end]]) - log(ns[idx[1]])) : NaN
        slope_all = (log(rms[end]) - log(rms[findfirst(>=(1000), ns)])) / (log(ns[end]) - log(ns[1]))
        results[label] = Eg
        emit("15(d)", string(label, ": starts used / failed, steps up to ", ns[end], " (", @sprintf("%.0f", time() - t0), " s)"), string(length(good), " / ", nst - length(good)), "-", "INFO")
        emit("15(d)", string(label, ": rms of E(n) at n = ", join(ns, ", ")), join((fmt(x) for x in rms), ", "), "-", "INFO")
        emit("15(d)", string(label, ": log-log slope of the rms over n in [1e4, ", ns[end], "]"), slope, comp ? "0.5 within 0.1" : "reported (negative control)", comp ? pf(abs(slope - 0.5) <= 0.1) : "INFO")
        emit("15(d)", string(label, ": slope over the whole range from n = 1000"), slope_all, "-", "INFO")
        # D113: the rule of criterion 15(a) for the start of the fit, the first n whose rms is at least 1.3 times the evaluation floor (the rms at n = 1000); reported next to the [1e4, 3e6] slope
        fl = findall(i -> ns[i] >= 10000 && rms[i] >= 1.3 * rms[1], eachindex(ns))
        slope_fl = length(fl) >= 2 ? (log(rms[end]) - log(rms[fl[1]])) / (log(ns[end]) - log(ns[fl[1]])) : NaN
        emit("15(d)", string(label, ": slope from the first n with rms >= 1.3 x floor (n = ", isempty(fl) ? "none" : string(ns[fl[1]]), ")"), slope_fl, "reported", "INFO")
        mu = vec(sum(Eg, dims = 2)) ./ length(good); se = [sqrt(sum(abs2, Eg[i, :] .- mu[i]) / (length(good) - 1)) / sqrt(length(good)) for i in axes(Eg, 1)]
        emit("15(d)", string(label, ": mean of E(n) over the starts in standard errors at n = ", join(ns, ", ")), join((fmt(mu[i] / se[i]) for i in eachindex(mu)), ", "), "-", "INFO")
        emit("15(d)", string(label, ": local slopes of the rms between consecutive n"), join((fmt((log(rms[i + 1]) - log(rms[i])) / (log(ns[i + 1]) - log(ns[i]))) for i in 1:(length(ns) - 1)), ", "), "-", "INFO")
    end
    serialize(joinpath(SCR, "ac15d.jls"), Dict("ns" => ns, "results" => results, "tau" => tau))
end

# ------------------------------------------------------------------------------------------------ GPUHOOK (D113, F4): the CUDA extension defines the ensemble device hook
section("GPUHOOK") do
    ext = Base.get_extension(Chaoslab, :ChaoslabCUDAExt)
    emit("18", "extension ChaoslabCUDAExt loaded (CUDA imported, no device code run)", ext !== nothing, "true", pf(ext !== nothing))
    ext === nothing && return
    has = hasmethod(Chaoslab.Integrators.ensemble_backend, Tuple{ext.CUDADevice})
    emit("18", "Integrators.ensemble_backend(::CUDADevice) is defined (a MethodError of advance_ensemble(...; backend = :gpu) otherwise)", has, "true", pf(has))
    be = has ? Chaoslab.Integrators.ensemble_backend(ext.CUDADevice()) : nothing
    emit("18", "the hook returns the CUDA backend", be isa CUDA.CUDABackend, "true", pf(be isa CUDA.CUDABackend))
    emit("18", "CUDA.functional() on this machine (informational; the GPU runs themselves are in verify/gpu.jl, not run here)", CUDA.functional(), "-", "INFO")
end

# ------------------------------------------------------------------------------------------------ AC18T: 1 thread against N threads, bitwise
const T_CASES = ((:S2, () -> Leapfrog(Float64, :DKD)), (:S5, () -> ForestRuth(Float64)), (:S8, () -> PEFRL(Float64)), (:S12_6, () -> KahanLi(Float64, 6)), (:SABA3, () -> SABA(Float64, 3)))
function t_ensembles()
    rng = StableRNG(VSEED + 18); out = Dict{String, Any}()
    for (nm, mk) in T_CASES, kind in (:f4, :f2)
        tau, n = IC_NH[kind]; n = min(n, 200)
        x0 = [ic_to(Float64, kind === :f2 ? ic_kepler_state(0.5, 2 * Float64(pi) * rand(rng)) : ic_draw(rng, kind)) for _ in 1:100000]
        out[string(nm, "_", kind)] = advance_ensemble(mk(), ic_flow(Float64, kind), x0, tau, n; backend = :cpu)[1]
    end
    out["GAUSS"] = advance_ensemble(Gauss(Float64, 2), ic_flow(Float64, :f4), [ic_to(Float64, ic_draw(rng, :f4)) for _ in 1:10000], 0.05, 20; backend = :cpu)
    return out
end
section("AC18child") do
    serialize(get(ENV, "IC_CHILD_OUT", joinpath(SCR, "child.jls")), t_ensembles())
end
section("AC18T") do
    child = joinpath(SCR, "child-t1.jls")
    cmd = addenv(`$(Base.julia_cmd()) -t 1 --project=$(joinpath(ROOT, "verify")) $(@__FILE__)`, "IC_ONLY" => "AC18child", "IC_CHILD_OUT" => child, "IC_LOG" => joinpath(SCR, "child-t1.log"), "IC_SCRATCH" => SCR)
    mine = t_ensembles()
    run(cmd)
    theirs = deserialize(child)
    for k in sort(collect(keys(mine)))
        a = mine[k]; b = theirs[k]
        same = k == "GAUSS" ? (isequal(a[1], b[1]) && isequal(a[2], b[2])) : isequal(a, b)
        emit("18", string("orbit results bitwise identical with ", Threads.nthreads(), " and 1 CPU threads: ", k, k == "GAUSS" ? " (10^4 orbits, states and statuses)" : " (10^5 orbits)"), same, "bitwise", pf(same))
    end
    emit("18", "GPU against CPU (10^5 orbits, double-double kernel compile), G1 and G2 10^4 orbits CPU against GPU", "not run", "-", "NOT RUN: needs verify/gpu.jl (writes data/gpu-benchmark.toml) and a device hook ensemble_backend(::CUDADevice) that ext/ChaoslabCUDAExt.jl does not define")
end

# ------------------------------------------------------------------------------------------------ MARGINS: worst measured ratio (measured / bound) of the criteria whose package tests only pass or fail
# Same bounds and flows as the package tests (the helpers are the test-side ones), at the spec sizes; each line is one number per criterion (at most 1 passes).
section("MARGINS") do
    ratio_line(ac, what, r) = emit(ac, what, r, "<= 1", pf(r <= 1))
    # 1: sum conditions in 256 bits
    r1 = 0.0; wm = ""
    setprecision(BigFloat, IC_BITS) do
        for M in IC_EXPL
            m = M.mk(BigFloat); c, d = coefficients(m); k = stage_count(m); et = ic_etab(BigFloat, M.dtab)
            r = Float64(max(abs(sum(c) - 1), abs(sum(d) - 1))) / (4 * k * et)
            r > r1 && (r1 = r; wm = M.id)
        end
    end
    ratio_line("1", string("sum conditions, worst |sum - 1| / (4 k eps_tab) over ", length(IC_EXPL), " explicit tables in 256 bits (worst: ", wm, ")"), r1)
    # 9: symplecticity
    rng = StableRNG(VSEED + 9)
    for T in (Float64, Double64)
        worst = 0.0; wid = ""
        for M in IC_ALLX, (kind, tau) in IC_SYMFL
            m = M.mk(T); flow = ic_flow(T, kind); d = ic_dim(kind)
            for _ in 1:200
                z = ic_to(T, ic_draw(rng, kind)); Mx = Matrix(tangent(ic_imap(m, flow, T(tau)), z))
                r = Float64(maximum(abs, Mx' * Matrix(symplectic_form(T, d)) * Mx - symplectic_form(T, d))) / (8 * d * ic_eps(T) * Float64(maximum(abs, Mx))^2)
                r > worst && (worst = r; wid = string(M.id, " ", kind))
            end
        end
        ratio_line("9", string("symplecticity, worst residual / (8 d eps ||M||^2), ", T, ", ", length(IC_ALLX), " methods x 3 flows x 200 points (worst: ", wid, ")"), worst)
    end
    # 10: reversibility of the palindromic tables on F4 (tau = 0.05, n_H = 200) and F2 (tau = 0.01, n_H = 1000), Float64
    worst = 0.0; wid = ""
    for M in IC_ALLX, kind in (:f4, :f2)
        M.sym === true || continue
        tau, n = IC_NH[kind]; m = M.mk(Float64); flow = ic_flow(Float64, kind)
        for _ in 1:100, nn in (1, n)
            z = ic_to(Float64, ic_draw(rng, kind)); y = advance_n(m, flow, advance_n(m, flow, z, tau, nn), -tau, nn)
            r = Float64(maximum(abs, collect(y) .- collect(z))) / (16 * ic_B(Float64, ic_mmax(z), nn, tau) * 4)
            r > worst && (worst = r; wid = string(M.id, " ", kind, " n = ", nn))
        end
    end
    ratio_line("10", string("reversibility, worst |S(-tau)^n S(tau)^n y - y| / (16 B(n) x 4), Float64 (worst: ", wid, ")"), worst)
    # 11: tangent step against ForwardDiff, Float64, all 41 explicit methods
    worst = 0.0; wid = ""
    for M in IC_ALLX, (kind, tau) in IC_SYMFL[2:3]
        m = M.mk(Float64); flow = ic_flow(Float64, kind); k = max(1, stage_count(m))
        for _ in 1:30
            z = ic_to(Float64, ic_draw(rng, kind)); h = tau
            A = Matrix(tangent(ic_imap(m, flow, h), z)); B = ic_jacobian(w -> advance(m, flow, w, h), z)
            r = Float64(maximum(abs, A .- B)) / (16 * ic_eps(Float64) * (1 + Float64(maximum(abs, A))) * k)
            r > worst && (worst = r; wid = string(M.id, " ", kind))
        end
    end
    ratio_line("11", string("tangent step against ForwardDiff, worst error / (16 eps (1 + |M|) k), Float64 (worst: ", wid, ")"), worst)
    # 13: bounded energy, no drift; 1e5 steps, 1000 samples, representative S2, S5, S12-6 on three flows
    w1 = 0.0; w2 = 0.0
    for M in (IC_S[3], IC_S[7], IC_S[19]), (kind, z0) in ((:f4, IC_START[:f4]), (:f2, IC_START[:f2]), (:f1_3, IC_START[:f1_3]))
        T = ic_ac13_type(M); tau = T(1) / 20
        E = ic_energy_series(M.mk(T), kind, ic_flow(T, kind), z0, tau, 100000, 1000); half = 500
        w1 = max(w1, maximum(E[(half + 1):end]) / (2 * maximum(E[1:half]) + 16 * ic_eps(T) * ic_mmax(z0)))
        w2 = max(w2, abs(ic_lsslope(E[(half + 1):end]) * half) / (0.1 * maximum(E)))
    end
    ratio_line("13", "energy: max over the second half / (2 x max over the first half + 16 eps M), S2, S5, S12-6 on F4, F2, F1 at tau = 1/20, 1e5 steps", w1)
    ratio_line("13", "energy: |least-squares slope x run length| / (0.1 max), same cases", w2)
    # 14: angular momentum of F2, 1e5 steps
    worst = 0.0; wid = ""
    for e in (0.1, 0.5), M in IC_ALLX
        z0 = ic_kepler_state(e, 0.7); z = advance_n(M.mk(Float64), ic_flow(Float64, :f2), ic_to(Float64, z0), 1 / 50, 100000)
        dL = setprecision(BigFloat, IC_BITS) do
            Float64(abs(ic_L(ic_to(BigFloat, collect(z))) - ic_L(ic_to(BigFloat, z0))))
        end
        r = dL / ic_B(Float64, ic_mmax(z0), 100000, 1 / 50)
        r > worst && (worst = r; wid = string(M.id, " e = ", e))
    end
    ratio_line("14", string("angular momentum of F2, worst |L(n) - L(0)| / B(n), n = 1e5, ", length(IC_ALLX), " methods x e in {0.1, 0.5} (worst: ", wid, ")"), worst)
    # 17: element-type agreement, 5 orbits per method and flow
    for T in (Float64, Double64)
        worst = 0.0; wid = ""
        for M in IC_ALLX, kind in (:f4, :f2)
            tau, n = IC_NH[kind]
            for _ in 1:5
                z0 = ic_draw(rng, kind)
                zb = setprecision(BigFloat, IC_BITS) do
                    Vector(advance_n(M.mk(BigFloat), ic_flow(BigFloat, kind), ic_to(BigFloat, z0), BigFloat(tau), n))
                end
                r = ic_err(advance_n(M.mk(T), ic_flow(T, kind), ic_to(T, z0), T(tau), n), zb) / ic_B(T, ic_mmax(z0), n, tau)
                r > worst && (worst = r; wid = string(M.id, " ", kind))
            end
        end
        ratio_line("17", string("element-type agreement ", T, " against 256 bits, worst error / B(n), n = n_H (worst: ", wid, ")"), worst)
    end
    emit("17", "effective epsilon of double-double", ic_eps(Double64), "<= 2^-100 = 7.9e-31", pf(ic_eps(Double64) <= 2.0^-100))
    # 16: allocations
    na = 0; nbad = 0
    for T in (Float64, Double64), M in IC_ALLX, kind in (:f1_1, :f1_3, :f4)
        m = M.mk(T); flow = ic_flow(T, kind); z = ic_to(T, IC_START[kind]); d = ic_dim(kind); W = SMatrix{d, d, T}(Matrix{T}(I, d, d))
        a1 = ic_alloc(advance, m, flow, z, T(0.05)); a2 = ic_alloc(tangent_advance, m, flow, z, W, T(0.05)); na += 2; nbad += (a1 != 0) + (a2 != 0)
    end
    emit("16", "bytes allocated by advance and tangent_advance after warm-up: calls with nonzero allocation / calls (all explicit methods, F1 n = 1, 3, F4, Float64 and Double64)", string(nbad, " / ", na), "0", pf(nbad == 0))
end

# ------------------------------------------------------------------------------------------------ FIG
mk_T(T, lab) = lab == "S1 DK" ? SymplecticEuler(T, :DK) : lab == "S2 leapfrog" ? Leapfrog(T) : lab == "S4 Ruth" ? Ruth3(T) : lab == "S5 Forest-Ruth" ? ForestRuth(T) :
               lab == "S8 PEFRL" ? PEFRL(T) : lab == "S9 CSS" ? CalvoSanzSerna4(T) : lab == "S11 Yoshida 6A" ? Yoshida(T, 6, :A) : lab == "S12 KL 6" ? KahanLi(T, 6) :
               lab == "S12 KL 8" ? KahanLi(T, 8) : SofroniouSpaletta10(T)

section("FIG") do
    # order plot: final-state error against tau on F4 (libration amplitude 1 rad, T = 8), reference O1 (Gauss s = 8, 256 bits)
    z0 = IC_START[:f4]; flow = ic_flow(Float64, :f4); ref = ic_ref(:f4, z0, 8)
    sel = [("S1 DK", () -> SymplecticEuler(Float64, :DK), 1), ("S2 leapfrog", () -> Leapfrog(Float64), 2), ("S4 Ruth", () -> Ruth3(Float64), 3), ("S5 Forest-Ruth", () -> ForestRuth(Float64), 4),
           ("S8 PEFRL", () -> PEFRL(Float64), 4), ("S9 CSS", () -> CalvoSanzSerna4(Float64), 4), ("S11 Yoshida 6A", () -> Yoshida(Float64, 6, :A), 6), ("S12 KL 6", () -> KahanLi(Float64, 6), 6),
           ("S12 KL 8", () -> KahanLi(Float64, 8), 8), ("S13 SS 10", () -> SofroniouSpaletta10(Float64), 10)]
    data = Dict{String, Any}()
    for (lab, mk, p) in sel
        T = ic_slope_type(p); m = mk_T(T, lab)
        tau0 = ic_pick_tau0(m, ic_flow(T, :f4), z0, 8, ref, T)
        taus, errs, ns = ic_scan(m, ic_flow(T, :f4), z0, 8, ref, tau0, T; J = 8)
        dtab = Dict("S8 PEFRL" => 16, "S9 CSS" => 18, "S11 Yoshida 6A" => 15, "S12 KL 6" => 20, "S12 KL 8" => 20, "S13 SS 10" => 32)
        floor_ = 1e3 * max(ic_eps(T), 10.0^(-get(dtab, lab, 99)))   # criterion 4: slopes only above 10^3 eps_tab
        sl = ic_slopes(errs, fill(floor_, length(errs)))
        data[lab] = (taus, errs, p, T, sl, floor_)
        say(string("FIG | ", lab, " (", T, ") tau0 = ", tau0, " slopes ", join((fmt(s) for s in sl), ", ")))
    end
    serialize(joinpath(SCR, "fig-order.jls"), data)
    fig = CM.Figure(size = (900, 650))
    ax = CM.Axis(fig[1, 1], xscale = log10, yscale = log10, xlabel = "tau", ylabel = "final-state error, F4 pendulum, T = 8", title = "integrator-catalog: measured order against the 256-bit reference")
    cols = CM.Makie.wong_colors()
    for (i, (lab, mk, p)) in enumerate(sel)
        taus, errs, _, T, _, floor_ = data[lab]
        good = findall(e -> e > floor_, errs)
        CM.scatterlines!(ax, taus[good], errs[good], label = string(lab, " (p = ", p, ")"), color = cols[mod1(i, length(cols))], marker = i % 2 == 0 ? :circle : :utriangle)
    end
    CM.axislegend(ax, position = :rb, labelsize = 10)
    CM.save(joinpath(OUT, "integrator-catalog-order.png"), fig)
    # Brouwer's law: rms of the energy error against n (AC15 explicit order 8, AC15(d) Gauss s = 4)
    fig2 = CM.Figure(size = (900, 600))
    ax2 = CM.Axis(fig2[1, 1], xscale = log10, yscale = log10, xlabel = "n steps (tau = 2/240)", ylabel = "rms over starts of H(y_n) - H(y_0)", title = "criterion 15: Brouwer's law, Henon-Heiles H0 = 1/8")
    f15 = joinpath(SCR, "ac15.jls")
    if isfile(f15)
        d = deserialize(f15); ns = d["ns"]
        CM.scatterlines!(ax2, ns, d["rc"], label = "S12 order 8, compensated", marker = :circle)
        CM.scatterlines!(ax2, ns, d["ru"], label = "S12 order 8, plain summation", marker = :rect)
    end
    f15d = joinpath(SCR, "ac15d.jls")
    if isfile(f15d)
        d = deserialize(f15d); ns = d["ns"]
        for (lab, E) in d["results"]
            rms = [sqrt(sum(abs2, E[i, :]) / size(E, 2)) for i in axes(E, 1)]
            CM.scatterlines!(ax2, ns, rms, label = string("Gauss s = 4, ", startswith(lab, "standard") ? "standard implementation" : startswith(lab, "diagnostic") ? "Hairer2008 iteration, plain summation" : "Hairer2008 implementation, compensated"), marker = :diamond)
        end
    end
    isfile(f15) && let ns = deserialize(f15)["ns"]
        CM.lines!(ax2, ns, 2e-17 .* sqrt.(ns ./ 1e4), color = :gray, linestyle = :dash, label = "n^(1/2)")
        CM.lines!(ax2, ns, 2e-17 .* (ns ./ 1e4), color = :lightgray, linestyle = :dot, label = "n")
    end
    CM.axislegend(ax2, position = :lt, labelsize = 10)
    CM.save(joinpath(OUT, "integrator-catalog-brouwer.png"), fig2)
    say("FIG | wrote integrator-catalog-order.png, integrator-catalog-brouwer.png")
end
close(LOG)
