# Stage 04 verification of module megno (reducer MEGNO, check_frame; layer Indicators; catalog row I8).
# Run from the repo root:  julia -t 2 --project=verify verify/megno.jl
# Spec: workspaces/module-build/stages/01-spec/output/megno-spec.md (criteria 1 to 15, Amendment D74). Criterion 15 (CPU vs GPU) is in verify/gpu.jl, tag "megno".
# Section "suite": the package suite at full size (MG_FULL_RUN = true: every multiple of 10 as a checkpoint, tau = 1 at n = 1e6 in AC2, Double64 at 1e5 in AC1(v), every n in
# AC5(b)); the @info lines of the suite (the realized ratios R2, the recorded lines R1, R3 to R6) are forwarded into the log by a logger.
# Other sections: one printed line per measurement  AC<n> | quantity | value | bound | PASS/FAIL/INFO  with the MEASURED number, reusing the helpers of test/megno/mg_helpers.jl
# (one home for B_M, Psi_max, the BigFloat reference, the M5 and M6 forms). Values are "worst ratio error / bound" (pass when <= 1) unless stated. The cases use verify seeds
# (MG_SEED + 58000 + 50 ac + sub), distinct from the package tests (MG_SEED + 35000 ..).
# MG_ONLY="AC14,figure" runs only the named sections (debugging aid) and MG_LOG names the log file; the default runs everything.
# MG_CHILD=1 is the child-process mode of the thread-invariance check of AC12(c).
using Printf, Test, Random, LinearAlgebra, Statistics, Logging
using Chaoslab, StaticArrays, DoubleFloats, StableRNGs

# ---------------------------------------------------------------------------------------------------
# Thread-count invariance (AC12c): the same outputs hashed in this process and in child processes with other thread counts.
function v_flatten!(v, x)
    if x isa Number
        push!(v, Float64(x))
    elseif x isa Union{Tuple, NamedTuple, AbstractArray}
        foreach(y -> v_flatten!(v, y), x)
    else
        error("v_flatten!: $(typeof(x))")
    end
    return v
end
function thread_hashes()
    out = Dict{String, String}()
    seed = 20261012 + 58000 + 50 * 12 + 3
    chains = (("MEGNO", (MEGNO(),)), ("MEGNO,SALI", (MEGNO(replace = false), SALI())), ("FLI,MEGNO,SALI", (FLI(replace = false), MEGNO(replace = false), SALI())))
    for (cname, reducers) in chains, tau in (1, 4), chunk in (64, 1000)
        rng = StableRNG(seed)
        x0 = [SVector(2pi * rand(rng), 2pi * rand(rng)) for _ in 1:4096]
        res = propagate(StandardMap(2pi), x0, 1000; k = 2, tau = tau, seed = seed, reducers = reducers, chunk = chunk, backend = :cpu)
        v = Float64[]
        v_flatten!(v, res.x); v_flatten!(v, res.W); v_flatten!(v, res.finals)
        out["$cname tau=$tau chunk=$chunk"] = string(hash(reinterpret(UInt64, v)))
    end
    return out
end
if get(ENV, "MG_CHILD", "") == "1"
    for (k, h) in sort!(collect(thread_hashes())); println("HASH ", k, " = ", h); end
    println("THREADS ", Threads.nthreads())
    exit(0)
end

const MG_FULL_RUN = true
const ROOT = dirname(@__DIR__)
const TESTDIR = joinpath(ROOT, "test", "megno")
include(joinpath(TESTDIR, "mg_helpers.jl"))

vseed(ac, sub = 0) = MG_SEED + 58000 + 50 * ac + sub

import ChaosTools as CT
import CairoMakie as CM

const OUT = joinpath(ROOT, "workspaces", "module-build", "stages", "04-verify", "output")
mkpath(OUT)
const ONLY = filter(!isempty, split(get(ENV, "MG_ONLY", ""), ","))
const LOG = open(get(ENV, "MG_LOG", joinpath(OUT, "megno-verify-log.txt")), "w")
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
info(ac, quantity, value) = emit(ac, quantity, value, "reported", "INFO")
rat(x, B) = x <= 0 ? 0.0 : (B == 0 ? Inf : x / B)

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

const T0 = time()
const DATA = Dict{Any, Any}()
const EPS64 = mg_eps(Float64)
say("# megno verification, julia $(VERSION), threads $(Threads.nthreads()), base seed $MG_SEED (verify offset 58000), time $(round(time())) (unix)")
say("# timings in this log were taken on a machine that was not otherwise idle")

# ---------------------------------------------------------------------------------------------------
# The package suite at full size; its @info lines are forwarded into this log (ratio keys are gated at <= 1)
struct SuiteLogger <: Logging.AbstractLogger end
Logging.min_enabled_level(::SuiteLogger) = Logging.Info
Logging.shouldlog(::SuiteLogger, args...) = true
Logging.catch_exceptions(::SuiteLogger) = true
function Logging.handle_message(::SuiteLogger, level, message, _module, group, id, file, line; kwargs...)
    kw = Dict{Symbol, Any}(kwargs)
    ratios = [k for k in (:worst, :worst_b, :worst_c) if haskey(kw, k)]
    rest = join(("$k = $(fmt(v))" for (k, v) in kw if !(k in ratios)), "; ")
    for k in ratios
        v = kw[k] isa Real ? Float64(kw[k]) : maximum(Float64.(kw[k]))
        ratio_line("SUITE", string(message, " [", k, "] ", rest), v)
    end
    isempty(ratios) && emit("SUITE", string(message, " ", rest), "-", "reported", "INFO")
end

section("suite") do
    local ts
    t = @elapsed try
        ts = with_logger(SuiteLogger()) do
            @testset verbose = true "megno verify" begin
                include(joinpath(TESTDIR, "runtests.jl"))
            end
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
    ok && foreach(r -> r isa Test.AbstractTestSet && (cc = Test.get_test_counts(r); info("SUITE", "testset '$(r.description)': assertions passed (fail + error = $(cc.fails + cc.errors + cc.cumulative_fails + cc.cumulative_errors))", cc.passes + cc.cumulative_passes)), ts.results)
    emit("SUITE", "package suite at full size (MG_FULL_RUN), assertions passed (fail + error = $nbad); wall time $(round(t; digits = 1)) s", npass, "0 fail, 0 error", pf(ok && nbad == 0 && npass > 0))
end

# ---------------------------------------------------------------------------------------------------
# AC1: algebra of the definition (ratios error / 6 B_S for the closed forms and the M5 sums, negative controls)
function v_const_ref(lam::BigFloat, tau::Int, K::Int)      # closed forms of M8(a) at the handoffs 1..K
    f = [mg_const_forms(lam, tau, k) for k in 1:K]
    return (Y = [x.Y for x in f], Ybar = [x.Ybar for x in f], Ybar20 = [x.Ybar20 for x in f])
end
worst_vs(Y, Yb, Y20, ref, B) = maximum(max(mg_err(Y[k], ref.Y[k]), mg_err(Yb[k], ref.Ybar[k]), mg_err(Y20[k], ref.Ybar20[k])) / B[k] for k in eachindex(Y))

section("AC1") do
    for T in MG_TYPES
        epsT = mg_eps(T)
        I2 = SMatrix{2, 2, T}(1, 0, 0, 1)
        x0 = SVector{2, T}(0, 0)
        snaps(kind, N, tau; W = I2) = mg_run((MEGNO(),), mg_Diag(kind), [x0], [W], N, tau; all = true)[1][1].snaps
        B6(kind, n, tau) = 6 * 8 * epsT * (n + mg_G_diag(kind, tau))
        wi = 0.0; ctl = nctl = 0
        lam0 = BigFloat(mg_c(:const, 1))
        for tau in (1, 5)
            sn = snaps(:const, 40, tau); K = 40 ÷ tau
            B = [B6(:const, k * tau, tau) for k in 1:K]
            sh = [lam0 * tau * k for k in 1:K]
            for (i, sg) in ((1, 1), (2, -1))
                ref = v_const_ref(sg * lam0, tau, K)
                wi = max(wi, worst_vs(mg_col(sn, :Y, i), mg_col(sn, :Ybar, i), mg_col(sn, :Ybar20, i), ref, B))
                Ls = [sg * lam0 for _ in 1:K]
                bads = [worst_vs(Ls, Ls, Ls, ref, B)]
                for kw in ((noweight = true,), (nofactor = true,), (mean0 = true,))
                    bad = mg_M5(sg .* sh, tau; kw...)
                    push!(bads, worst_vs(bad.Y, bad.Ybar, bad.Ybar20, ref, B))
                end
                ctl += count(<=(1), bads); nctl += length(bads)
            end
        end
        ratio_line("AC1", "(i) a_m = 0.05, j <= 40, tau in {1,5}, both columns: max |M - M8(a) closed form| / 6 B_S(n), $T", wi)
        sn = snaps(:log1, 200, 1)
        B = [B6(:log1, j, 1) for j in 1:200]
        lnfact = setprecision(() -> accumulate(+, [log(BigFloat(k)) for k in 1:200]), BigFloat, MG_BITS)
        Ycf = setprecision(() -> [2 * (log(BigFloat(j + 1)) - lnfact[j] / j) for j in 1:200], BigFloat, MG_BITS)
        sh1 = setprecision(() -> [log(BigFloat(k + 1)) for k in 1:200], BigFloat, MG_BITS)
        r1, r2 = mg_M5(sh1, 1), mg_M5(-sh1, 1)
        wY = maximum(max(mg_err(sn[j].Y[1], Ycf[j]), mg_err(sn[j].Y[2], -Ycf[j])) / B[j] for j in 1:200)
        wii = max(worst_vs(mg_col(sn, :Y, 1), mg_col(sn, :Ybar, 1), mg_col(sn, :Ybar20, 1), r1, B), worst_vs(mg_col(sn, :Y, 2), mg_col(sn, :Ybar, 2), mg_col(sn, :Ybar20, 2), r2, B))
        ratio_line("AC1", "(ii) s(k) = ln(1 + k), j <= 200: max |Y - 2 (ln(j+1) - ln(j!)/j)| / 6 B_S, $T", wY)
        ratio_line("AC1", "(ii) same, Y, Ybar, Ybar20 against the BigFloat M5 sums, both columns, $T", wii)
        sn = snaps(:alt, 200, 1)
        c = setprecision(() -> cumsum([BigFloat(mg_c(:alt, m)) for m in 1:200]), BigFloat, MG_BITS)
        B = [B6(:alt, j, 1) for j in 1:200]
        w3 = max(worst_vs(mg_col(sn, :Y, 1), mg_col(sn, :Ybar, 1), mg_col(sn, :Ybar20, 1), mg_M5(c, 1), B), worst_vs(mg_col(sn, :Y, 2), mg_col(sn, :Ybar, 2), mg_col(sn, :Ybar20, 2), mg_M5(-c, 1), B))
        ratio_line("AC1", "(iii) +0.2, -0.2 alternating, j <= 200: max |M - M5 sums| / 6 B_S, $T; |Ybar(200)|, |Ybar20(200)| max = $(fmt(maximum(abs, (sn[200].Ybar..., sn[200].Ybar20...)))) (bound 0.1)", max(w3, maximum(abs, (sn[200].Ybar..., sn[200].Ybar20...)) / 0.1))
        sA, sB = snaps(:pathA, 40, 1), snaps(:pathB, 40, 1)
        B = [B6(:pathB, j, 1) for j in 1:40]
        cA = setprecision(() -> cumsum([BigFloat(mg_c(:pathA, m)) for m in 1:40]), BigFloat, MG_BITS)
        cB = setprecision(() -> cumsum([BigFloat(mg_c(:pathB, m)) for m in 1:40]), BigFloat, MG_BITS)
        w4 = max(worst_vs(mg_col(sA, :Y, 1), mg_col(sA, :Ybar, 1), mg_col(sA, :Ybar20, 1), mg_M5(cA, 1), B), worst_vs(mg_col(sB, :Y, 1), mg_col(sB, :Ybar, 1), mg_col(sB, :Ybar20, 1), mg_M5(cB, 1), B))
        SA, SB = sA[40].sum[1] + sA[40].comp[1], sB[40].sum[1] + sB[40].comp[1]
        emit("AC1", "(iv) memory M7: Y_A(40) = $(fmt(Float64(sA[40].Y[1]))) (2.05), Y_B(40) = $(fmt(Float64(sB[40].Y[1]))) (1.05), |s_A(40) - s_B(40)| = $(fmt(Float64(abs(SA - SB)))) (<= 2 B_S), worst vs M5 sums / 6 B_S = $(fmt(w4)), $T", w4,
             "<= 1", pf(w4 <= 1 && abs(sA[40].Y[1] - 2.05) < 1e-12 && abs(sB[40].Y[1] - 1.05) < 1e-12 && abs(SA - SB) <= 2 * 8 * epsT * (40 + mg_G_diag(:pathB, 1))))
        emit("AC1", "negative controls (s/j no memory, no weight, no factor 2, mean over 0..K-1) that wrongly satisfy the predicate (of $nctl), $T", ctl, "0", pf(ctl == 0))
        wv = 0.0
        for n in (T === Float64 ? [100_000, 4_000_000] : [100_000])
            sn5 = mg_run((MEGNO(),), mg_Diag(:long), [SVector{2, T}(0, 0)], [I2], n, 1; steps = [n])[1][1].snaps[1]
            lam = BigFloat(mg_c(:long, 1)); B = 6 * 8 * epsT * (n + mg_G_diag(:long, 1))
            for (i, sg) in ((1, 1), (2, -1))
                f = mg_const_forms(sg * lam, 1, n)
                wv = max(wv, max(mg_err(sn5.Y[i], f.Y), mg_err(sn5.Ybar[i], f.Ybar), mg_err(sn5.Ybar20[i], f.Ybar20)) / B)
            end
        end
        ratio_line("AC1", "(v) long run a_m = 0.1, n in {1e5$(T === Float64 ? ",4e6" : "")}: max |M - M8(a) closed form| / 6 B_S, $T (no integer weights; K^3 at 4e6 = $(fmt(Float64(4_000_000)^3)) > 2^63 = $(fmt(Float64(2)^63)))", wv)
    end
end

vlabel(label) = label === :half ? "0.5 regular class" : label === :two_pi ? "2 pi" : label === :ten ? "10" : "5"
vlabel_k(label) = label === :two_pi ? "2 pi" : label === :five ? "5" : "10"

# ---------------------------------------------------------------------------------------------------
# Negative controls of AC6 (Float32 accumulator) and AC7 (Mestre relation with the running maximum), with their measured numbers
function v_single_Y(s::AbstractVector{BigFloat}, n::Int)       # Y_1(n), tau = 1, of a double whose accumulator is a Float32 (as the AC6 control of the suite)
    return setprecision(BigFloat, MG_BITS) do
        R1 = Float32(0); prev = zero(BigFloat)
        for k in 1:n
            R1 += Float32(k) * Float32(s[k] - prev); prev = s[k]
        end
        Float64(2 * R1 / Float32(n))
    end
end
section("AC6c") do
    for label in (:two_pi, :ten)
        c = mg_case(label, Float64, 20)
        rows = mg_rows(label, Float64, 1, 20, 10_000)
        d = [(r = only(r for r in rows if r.j == j && r.n == 1000); sv = mg_ref_s(mg_mats_std(c.K, mg_pts(c.sys, c.x0[j], 1000)), c.W0[j], 1000);
              abs(v_single_Y(sv.s1, 1000) - Float64(r.r.Y[1])) / r.B[1]) for j in 1:3]
        emit("AC6", "K = $(vlabel(label)): negative control, Float32 accumulator at n = 1e3, min |Y_1(Float32) - Y_ref| / B_M over 3 starts", minimum(d), "> 1", pf(minimum(d) > 1))
    end
end
section("AC7c") do
    T = Float64; N = 10_000; epsT = EPS64
    c = mg_case(:half, T, 20)
    run = mg_run((MEGNO(replace = false), FLI(replace = false), SALI()), c.sys, c.x0, c.W0, N, 1; steps = collect(10:10:N))
    l1 = mg_run((Benettin(),), c.sys, c.x0, [SMatrix{2, 1, T}(W[:, 1]) for W in c.W0], N, 1; all = true)
    nctl = nfail = nneg = 0; minratio = Inf; Dmin = Inf; wd = 0.0; nsum = 0
    for j in 1:20
        msn, fsn = run[j][1].snaps, run[j][2].snaps
        sL = [k * s.lambda[1] for (k, s) in enumerate(l1[j][1].snaps)]
        smax = accumulate(max, abs.(sL))
        pre = setprecision(() -> accumulate(+, BigFloat.(sL)), BigFloat, MG_BITS)
        for (ms, fs) in zip(msn, fsn)
            n = ms.n
            tol = 16 * epsT * max(1, smax[n]); tolD = 16 * epsT * max(1, smax[n], fs.fli1)
            D = fs.fli1 - sL[n]
            Dmin = min(Dmin, D / tolD); nneg += D < -tolD; nsum += 1
            wd = max(wd, mg_err(ms.Y[1], setprecision(() -> 2 * ((BigFloat(fs.fli1) - BigFloat(D)) - pre[n] / n) + 2 * (BigFloat(fs.fli1) - BigFloat(D)) / n, BigFloat, MG_BITS)) / tolD)
            if D > 10 * tolD
                bad = setprecision(() -> 2 * (BigFloat(fs.fli1) - pre[n] / n) + 2 * BigFloat(fs.fli1) / n, BigFloat, MG_BITS)
                nctl += 1; nfail += mg_err(ms.Y[1], bad) > tol; minratio = min(minratio, mg_err(ms.Y[1], bad) / tol)
            end
        end
    end
    emit("AC7", "(d) K = 0.5 regular class, 20 starts x 1000 handoffs: D = FLI_1 - n X_L1 below -tol at $nneg of $nsum points (min D / tol = $(fmt(Dmin))); max |Y_1 - A2 with FLI_1 - D| / tol = $(fmt(wd))", wd, "<= 1", pf(nneg == 0 && wd <= 1))
    emit("AC7", "negative control, A2 with the running maximum FLI_1 in place of ln norm: pairs with D > 10 tol: $nctl, of which fail the relation (error > tol): $nfail; smallest error / tol among them", minratio, "> 1 for all", pf(nctl > 0 && nfail == nctl))
end

# ---------------------------------------------------------------------------------------------------
# AC8: growth laws at the ensemble level, measured against the allowances (verify seeds, 1000 starts)
section("AC8") do
    c = mg_case(:half, Float64, 20)
    run4 = mg_run((MEGNO(),), c.sys, c.x0, c.W0, 10_000, 1; steps = [10_000])
    info("AC8", "R3 K = 0.5 regular class, 20 starts, n = 1e4: Ybar_1 min / median / max; Ybar20_1 min / median / max (levels 2 and 0.5)", string(join(fmt.((minimum(r[1].snaps[1].Ybar[1] for r in run4), mg_median([r[1].snaps[1].Ybar[1] for r in run4]), maximum(r[1].snaps[1].Ybar[1] for r in run4))), " / "), "; ", join(fmt.((minimum(r[1].snaps[1].Ybar20[1] for r in run4), mg_median([r[1].snaps[1].Ybar20[1] for r in run4]), maximum(r[1].snaps[1].Ybar20[1] for r in run4))), " / ")))
    run = mg_run((MEGNO(),), c.sys, c.x0, c.W0, 100_000, 1; steps = [100_000])
    Yb1 = [r[1].snaps[1].Ybar[1] for r in run]; Y201 = [r[1].snaps[1].Ybar20[1] for r in run]
    ok = count(j -> 1.5 <= Yb1[j] <= 2.5 && 0.4 <= Y201[j] <= 0.6, 1:20)
    emit("AC8", "(a) K = 0.5 regular class, 20 starts, n = 1e5: starts with Ybar_1 in [1.5,2.5] and Ybar20_1 in [0.4,0.6]; Ybar_1 min / median / max = $(fmt(minimum(Yb1))) / $(fmt(mg_median(Yb1))) / $(fmt(maximum(Yb1))); Ybar20_1 = $(fmt(minimum(Y201))) / $(fmt(mg_median(Y201))) / $(fmt(maximum(Y201)))", ok, ">= 16", pf(ok >= 16))
    m = 1000
    for (ci, label) in enumerate((:two_pi, :five, :ten))
        seed = vseed(8, ci)
        sys = StandardMap(mg_K(Float64, label))
        x0 = [mg_cast(Float64, s) for s in mg_starts(StableRNG(seed), m)]
        W0 = [mg_frame(seed, j, Float64) for j in 1:m]
        w = collect(20:200)
        run = mg_run((MEGNO(),), sys, x0, W0, 200, 1; steps = collect(1:200))
        l1 = mg_run((Benettin(),), sys, x0, [SMatrix{2, 1, Float64}(W[:, 1]) for W in W0], 10_000, 1; steps = [10_000])
        X = [(r[1].snaps[1].sum[1] + r[1].snaps[1].comp[1]) / 10_000 for r in l1]
        mYb = [mg_slope(Float64.(w), Float64[s.Ybar[1] for s in r[1].snaps[20:200]]) for r in run]
        mY = [mg_slope(Float64.(w), Float64[s.Y[1] for s in r[1].snaps[20:200]]) for r in run]
        m20 = [mg_slope(Float64.(w), Float64[s.Ybar20[1] for s in r[1].snaps[20:200]]) for r in run]
        for (nm, a, al, fac) in (("2 m_Ybar", 2 .* mYb, 4 / 180, "first"), ("m_Y", mY, 4 / 180, "second"), ("4 m_20", 4 .* m20, 8 / 180, "third"))
            se = sqrt(mg_var(a) / m + mg_var(X) / m)
            d = abs(mg_mean(a) - mg_mean(X))
            ratio_line("AC8", "(b) K = $(vlabel_k(label)), 1000 starts: |mean($nm) - mean X(1e4)| = $(fmt(d)) (mean X = $(fmt(mg_mean(X)))) over 4 SE (= $(fmt(4 * se))) + allowance $(fmt(al)); in units of the allowance alone: $(fmt(d / al))", d / (4 * se + al))
        end
        hit = [(h = findfirst(s -> s.Ybar[1] >= 3, r[1].snaps); h === nothing ? nothing : r[1].snaps[h].n) for r in run]
        info("AC8", "R4 K = $(vlabel_k(label)): first handoff with Ybar_1 >= 3, median over starts that reach it; never reached within 200", "$(mg_median([n for n in hit if n !== nothing])); $(count(isnothing, hit))")
    end
end

# ---------------------------------------------------------------------------------------------------
# AC9: tau dependence (measured ratios against the three bounds)
section("AC9") do
    m, N = 20, 1000
    for T in MG_TYPES, label in (:half, :two_pi)
        epsT = mg_eps(T)
        lns = log(mg_sigma(Float64(mg_K(T, label))))
        O = mg_orbit_ref(label, T, m, N)
        base = Dict((r.j, r.n) => r for r in mg_rows(label, T, 1, m, N))
        for tau in (2, 5, 10)
            ra = rb = rc = 0.0
            dYmax = dYbmax = 0.0
            for r in mg_rows(label, T, tau, m, N)
                r.n in (100, 1_000) || continue
                r1 = base[(r.j, r.n)]
                o = O.per[r.j]
                for i in 1:2
                    ra = max(ra, Float64(abs(r.S[i] - r1.S[i])) / (16 * epsT * mg_Psi(o, i, r.n)))
                    Bsum = r.B[i] + r1.B[i]
                    ra_b = Float64(abs(r.Y[i] - r1.Y[i])) / ((tau - 1) * lns + Bsum)
                    ra_c = Float64(abs(r.Ybar[i] - r1.Ybar[i])) / (3 * (tau - 1) * lns + Bsum)
                    rb, rc = max(rb, ra_b), max(rc, ra_c)
                    dYmax = max(dYmax, Float64(abs(r.Y[i] - r1.Y[i])) / ((tau - 1) * lns)); dYbmax = max(dYbmax, Float64(abs(r.Ybar[i] - r1.Ybar[i])) / (3 * (tau - 1) * lns))
                end
            end
            K = vlabel(label)
            ratio_line("AC9", "K = $K, $T, tau = $tau vs 1, n in {100,1e3}, 20 starts, both columns: (a) max |s^tau - s^1| / (16 eps Psi_n)", ra)
            ratio_line("AC9", "K = $K, $T, tau = $tau: (b) max |Y_tau - Y_1| / ((tau-1) L + B_M^tau + B_M^1); without the B_M terms the largest |dY| / ((tau-1) L) = $(fmt(dYmax))", rb)
            ratio_line("AC9", "K = $K, $T, tau = $tau: (c) max |Ybar_tau - Ybar_1| / (3 (tau-1) L + B_M^tau + B_M^1); largest |dYbar| / (3 (tau-1) L) = $(fmt(dYbmax))", rc)
        end
    end
end

# ---------------------------------------------------------------------------------------------------
# AC10: precision ladder and representation (measured z-scores and ratios)
zscore(a, b) = abs(mean(a) - mean(b)) / sqrt(var(a) / length(a) + var(b) / length(b))
section("AC10") do
    for (labels, m, n) in (((:two_pi, :five), 1000, 1_000), ((:five, :ten), 200, 10_000)), label in labels
        a, b = mg_ens(Float64, label, m, n), mg_ens(Double64, label, m, n)
        zy, zb = zscore(a.Y, b.Y), zscore(a.Ybar, b.Ybar)
        emit("AC10", "(a) Float64 vs Double64, K = $(vlabel(label)), $m starts, n = $n: |difference of means| in combined SE for Y_1/n and Ybar_1/n (means $(fmt(mean(a.Y))) / $(fmt(mean(b.Y))))", max(zy, zb), "<= 4", pf(max(zy, zb) <= 4))
    end
    setprecision(BigFloat, MG_BITS) do
        rows = mg_rows(:two_pi, BigFloat, 1, 50, 60; bits = 512)
        ratio_line("AC10", "(b) BigFloat(256) on its own orbit vs 512-bit reference, K = 2 pi, 50 starts, n up to 60 (every multiple of 10): max |M - M_ref| / B_M over $(length(rows)) rows", maximum(r.ratio for r in rows))
        for label in (:two_pi, :five)
            a, b = mg_ens(BigFloat, label, 100, 300), mg_ens(Double64, label, 100, 300)
            z = max(zscore(a.Y, b.Y), zscore(a.Ybar, b.Ybar))
            emit("AC10", "(c) BigFloat vs Double64, K = $(vlabel(label)), 100 starts, n = 300: |difference of means| in combined SE, Y_1/n and Ybar_1/n", z, "<= 4", pf(z <= 4))
        end
    end
    rows = [r for r in mg_rows(:ten, Float64, 1, 20, 10_000) if r.n == 10_000]
    Ys = [Float64(r.Y[1]) for r in rows]
    w = maximum(max(r.e.Y[1] / r.B[1], r.e.Y[2] / r.B[2]) for r in rows)
    emit("AC10", "(d) K = 10, n = 1e4, Float64, 20 starts: Y_1 min / median / max = $(fmt(minimum(Ys))) / $(fmt(mg_median(Ys))) / $(fmt(maximum(Ys))) (band 1.4e4..1.9e4), max |Y - Y_ref| / B_M = $(fmt(w)), max |Y_1 - Y_ref,1| (not clipped) = $(fmt(maximum(abs(Float64(r.Y[1]) - Float64(r.r.Y[1])) for r in rows)))", w, "<= 1", pf(w <= 1 && all(y -> 1.4e4 < y < 1.9e4, Ys)))
end

# ---------------------------------------------------------------------------------------------------
# D74 item 7: Double64 K = 10 log sums against 8 eps Psi, Psi on the full handoff set and on the two grids
function full_psi(Ms, w0, n)            # Psi_g for g = 1..n (the full handoff set), Float64
    out = Vector{Float64}(undef, n)
    Threads.@threads :dynamic for g in 1:n
        out[g] = lb_Psi(view(Ms, 1:g), w0)
    end
    return out
end
const GRID_DEC = [n for n in sort(unique(vcat((d .* (1:10) for d in (1, 10, 100, 1000))...))) if n % 10 == 0]   # the 10-per-decade grid of D38
section("D74-7") do
    T = Double64
    epsT = mg_eps(T)
    c = mg_case(:ten, T, 20)
    N = 10_000
    function d74_one(j, nmax)
        pts = mg_pts(c.sys, c.x0[j], nmax)
        sv = mg_ref_s(mg_mats_std(c.K, pts), c.W0[j], nmax)
        probe = mg_Probe(MEGNO(); all = true)
        _, _, st = propagate_orbit(c.sys, c.x0[j], c.W0[j], nmax, 1, (probe,); orbit = j)
        sn = Chaoslab.reducer_final(probe, st[1]).snaps
        Ms = lb_tangents(c.sys, c.x0[j], nmax)
        ref = [mg_M5(sv.s1, 1), mg_M5(sv.s2, 1)]
        res = NamedTuple[]
        for i in 1:2
            psi = full_psi(Ms, c.W0[j][:, i], nmax)
            pmax = accumulate(max, psi)
            pgrid_dec = [maximum((psi[g] for g in GRID_DEC if g <= n); init = psi[1]) for n in 1:nmax]
            pgrid_10 = [maximum((psi[g] for g in 10:10:n); init = psi[1]) for n in 1:nmax]
            S = [Float64(abs(mg_big(sn[n].sum[i]) + mg_big(sn[n].comp[i]) - (i == 1 ? sv.s1 : sv.s2)[n])) for n in 1:nmax]
            sref = i == 1 ? sv.s1 : sv.s2
            Sr = [mg_err(sn[n].sum[i] + sn[n].comp[i], sref[n]) for n in 1:nmax]                       # the total read as T arithmetic (sum + comp rounded to T), as the suite's R2 line does
            hu = [Float64(abs(sref[n])) * 2.0^-107 for n in 1:nmax]                                      # half an ulp of a 106-bit significand at |s|
            eY = [Float64(max(mg_err(sn[n].Y[i], ref[i].Y[n]), mg_err(sn[n].Ybar[i], ref[i].Ybar[n]), mg_err(sn[n].Ybar20[i], ref[i].Ybar20[n]))) for n in 1:nmax]
            ratio_own = S ./ (8 * epsT .* psi)
            nw = argmax(ratio_own)
            gd = [n for n in GRID_DEC if n <= nmax]; gm = collect(10:10:nmax)                                   # evaluation points of the two grids (the bound is read where its Psi_max includes g = n)
            xrun = accumulate(max, [max(abs(p[1]), abs(p[2])) for p in pts])                                      # running max of |x|, |p| of the lifted orbit (the trig argument)
            push!(res, (i = i, own = maximum(ratio_own), full = maximum(S ./ (8 * epsT .* pmax)), gmult = maximum(S[n] / (8 * epsT * pgrid_10[n]) for n in gm),
                        gdec = maximum(S[n] / (8 * epsT * pgrid_dec[n]) for n in gd), series_full = maximum(eY ./ (48 * epsT .* pmax)), series_grid = maximum(eY[n] / (48 * epsT * pgrid_dec[n]) for n in gd),
                        psi_underestimate = maximum(pmax[n] / pgrid_dec[n] for n in gd),
                        rnd_full = maximum(Sr ./ (8 * epsT .* pmax)), rnd_dec = maximum(Sr[n] / (8 * epsT * pgrid_dec[n]) for n in gd), rnd_own = maximum(Sr ./ (8 * epsT .* psi)),
                        hu_dec = maximum(hu[n] / (8 * epsT * pgrid_dec[n]) for n in gd), hu_full = maximum(hu ./ (8 * epsT .* pmax)),
                        n_worst = nw, x_worst = xrun[nw], small_x = maximum((ratio_own[n] for n in 1:nmax if xrun[n] < 1e3); init = 0.0), n_small = count(n -> xrun[n] < 1e3, 1:nmax)))
        end
        res
    end
    r1000 = [setprecision(BigFloat, MG_BITS) do; d74_one(j, 1000); end for j in 1:20]
    f4(r) = string(fmt(r.own), " / ", fmt(r.full), " / ", fmt(r.gdec), " / ", fmt(r.gmult))
    allr(rrs, f) = maximum(f(r) for rr in rrs for r in rr)
    info("D74-7", "Double64, K = 10, 20 starts, n <= 1e3, both columns, every n: max |s - s_ref| / (8 eps X) with X = Psi_n (own n) / max over the FULL handoff set g <= n / max over the 10-per-decade grid / max over every multiple of 10",
         "$(fmt(allr(r1000, r -> r.own))) / $(fmt(allr(r1000, r -> r.full))) / $(fmt(allr(r1000, r -> r.gdec))) / $(fmt(allr(r1000, r -> r.gmult)))")
    info("D74-7", "same, but the total read as sum + comp ROUNDED to Double64 (what the suite's R2 line does): max |S - s_ref| / (8 eps X), X = Psi_n / full set / 10-per-decade grid; and half an ulp of a 106-bit significand at |s_ref| over 8 eps Psi_max: full set / 10-per-decade grid",
         "$(fmt(allr(r1000, r -> r.rnd_own))) / $(fmt(allr(r1000, r -> r.rnd_full))) / $(fmt(allr(r1000, r -> r.rnd_dec))); $(fmt(allr(r1000, r -> r.hu_full))) / $(fmt(allr(r1000, r -> r.hu_dec)))")
    info("D74-7", "same, three series: max |M - M_ref| / (48 eps Psi_max) with Psi_max over the full set / the 10-per-decade grid", "$(fmt(allr(r1000, r -> r.series_full))) / $(fmt(allr(r1000, r -> r.series_grid)))")
    info("D74-7", "same, ratio Psi_max(full set) / Psi_max(10-per-decade grid) at the worst n (how much the grid underestimates)", fmt(allr(r1000, r -> r.psi_underestimate)))
    worst = sortperm([maximum(r.gdec for r in rr) for rr in r1000]; rev = true)[1:3]
    starts = collect(1:20)      # all 20 starts to n = 1e4 (the suite's R2 line covers these; worst-first order is not needed)
    r10k = [(j, setprecision(BigFloat, MG_BITS) do; d74_one(j, N); end) for j in starts]
    for (j, rr) in r10k, r in rr
        info("D74-7", "start $j column $(r.i), n <= 1e4: max |s - s_ref| / (8 eps X), X = Psi_n / Psi_max full set / Psi_max 10-per-decade grid / Psi_max every multiple of 10; three series / (48 eps Psi_max full set); / (... grid); Psi underestimate by the grid", "$(f4(r)); $(fmt(r.series_full)); $(fmt(r.series_grid)); $(fmt(r.psi_underestimate)); worst own-n ratio at n = $(r.n_worst) with max(|x|,|p|) up to $(fmt(r.x_worst)); max own-n ratio over the n with max(|x|,|p|) < 1e3 ($(r.n_small) of 1e4 n) = $(fmt(r.small_x))")
    end
    rs = [r for (_, rr) in r10k for r in rr]
    for (j, rr) in r10k, r in rr
        info("D74-7", "start $j column $(r.i), n <= 1e4, total rounded to Double64: max |S - s_ref| / (8 eps X), X = Psi_n / Psi_max full set / Psi_max 10-per-decade grid; half-ulp at |s_ref| over 8 eps Psi_max full set / grid", "$(fmt(r.rnd_own)) / $(fmt(r.rnd_full)) / $(fmt(r.rnd_dec)); $(fmt(r.hu_full)) / $(fmt(r.hu_dec))")
    end
    wfull, wdec, wown, sfull = maximum(r.full for r in rs), maximum(r.gdec for r in rs), maximum(r.own for r in rs), maximum(r.series_full for r in rs)
    emit("D74-7", "RESULT n <= 1e4, all 20 starts: log sums against 8 eps Psi_max: full handoff set $(fmt(wfull)), 10-per-decade grid $(fmt(wdec)), own-n Psi_n $(fmt(wown))", wfull, "<= 1 on the full set means the excess is the grid", "INFO")
    rfull, rdec = maximum(r.rnd_full for r in rs), maximum(r.rnd_dec for r in rs)
    emit("D74-7", "RESULT n <= 1e4, total read as sum + comp rounded to Double64: against 8 eps Psi_max full handoff set $(fmt(rfull)), 10-per-decade grid $(fmt(rdec)) (the suite's worstS 1.64 is this reading on the grid)", rfull, "reported", "INFO")
    xw = [(r.x_worst, r.n_worst, r.own) for r in rs]
    xs_over = [t for t in xw if t[3] > 1]
    emit("D74-7", "cause: columns with an own-n ratio above 1: $(length(xs_over)) of $(length(rs)); max(|x|,|p|) of the lifted orbit at their worst n (min; max) = $(isempty(xs_over) ? "-" : string(fmt(minimum(t[1] for t in xs_over)), "; ", fmt(maximum(t[1] for t in xs_over)))); largest own-n ratio over the n with max(|x|,|p|) < 1e3 on any column = $(fmt(maximum(r.small_x for r in rs))) (the Double64 trig floor beyond |x| = 1e3 of the precision spec Amendment D35, `PR-trig`)", maximum(r.small_x for r in rs), "<= 1 if the cause is the trig floor", "INFO")
    nb = 0
    for j in 1:20
        ra = mg_run((MEGNO(),), c.sys, [c.x0[j]], [c.W0[j]], N, 1; steps = [1000, N])[1][1].snaps
        for i in 1:2
            rl = mg_run((Benettin(),), c.sys, [c.x0[j]], [SMatrix{2, 1, T}(c.W0[j][:, i])], N, 1; steps = [1000, N])[1][1].snaps
            nb += !(all(q -> ra[q].sum[i] == rl[q].sum[1] && ra[q].comp[i] == rl[q].comp[1], 1:2))
        end
    end
    emit("D74-7", "the MEGNO sums (sum and comp, Double64) differ from the L1 (Benettin, k = 1) sums of the same column at n = 1e3 and 1e4: columns of 40", nb, "0 (so the excess is inherited from L1, not from the weighted sums)", pf(nb == 0))
    emit("D74-7", "three series against B_M = 48 eps Psi_max with Psi_max on the full handoff set (M11 as stated), same starts", sfull, "<= 1", pf(sfull <= 1))
end

# ---------------------------------------------------------------------------------------------------
# AC11: health and frames (direct checks with the measured values)
const VCHAINS = ((MEGNO(),), (MEGNO(replace = false), SALI()), (FLI(replace = false), MEGNO(replace = false), SALI()))
section("AC11") do
    seed = vseed(11, 1)
    x0s(flagged) = [SVector{3, Float64}(0, 0, j in flagged ? 1 : 0) for j in 1:6]
    zero_map, col2_map, nan_map = mg_Fault(:zero, 0), mg_Fault(:col2, 0), mg_Fault(:nan, 3)
    bad = n = 0
    for chain in VCHAINS
        res = nothing
        e = mg_health_error(() -> (res = propagate(zero_map, x0s((4,)), 12; k = 2, tau = 3, seed = seed, reducers = chain, chunk = 2)))
        bad += !(e isa NumericalHealthError && res === nothing && e.what == "deviation matrix" && e.orbit == 4 && e.step == 3 && e.backend === :cpu && e.column == 1 && e.value == 0.0); n += 1
        e = mg_health_error(() -> propagate(zero_map, x0s((5, 2)), 12; k = 2, tau = 3, seed = seed, reducers = chain, chunk = 2))
        bad += !(e isa NumericalHealthError && e.orbit == 2 && e.step == 3); n += 1
        eye = [SMatrix{2, 2, Float64}(1, 0, 0, 1) for _ in 1:6]
        e = mg_health_error(() -> propagate(col2_map, x0s((5,)), 12; k = 2, tau = 3, seed = seed, reducers = chain, chunk = 2, W0 = eye))
        bad += !(e isa NumericalHealthError && e.orbit == 5 && e.step == 3 && e.column == 2 && e.value == 0.0); n += 1
        e = mg_health_error(() -> propagate(nan_map, x0s((3,)), 20; k = 2, tau = 5, seed = seed, reducers = chain, chunk = 2))
        bad += !(e isa NumericalHealthError && e.what == "deviation matrix" && e.orbit == 3 && e.step == 5 && e.column == 0 && e.threshold == Inf && !isfinite(e.value)); n += 1
    end
    e = try; propagate(mg_Boom(), [SVector(0.0, 0.0)], 4; k = 1, tau = 2, seed = seed, reducers = (MEGNO(),)); nothing; catch ex; ex; end
    bad += !(e isa ArgumentError && occursin("MEGNO needs k = 2, got k = 1", e.msg)); n += 1
    emit("AC11", "(a), (b) zero tangent (orbit 4, step 3, column 1), column-2 collapse, lowest orbit id first, NaN tangent (column 0, threshold Inf), k = 1 ArgumentError; three chains: failing checks", bad, "0 of $n", pf(bad == 0))
    # (c) frame check
    epsT = EPS64
    rejected = count(W -> mg_health_error(() -> check_frame(MEGNO(), W; orbit = 7)) isa NumericalHealthError,
                     (SMatrix{2, 2, Float64}(2, 0, 0, 1), SMatrix{2, 2, Float64}(1, 0, 0, 3), SMatrix{2, 2, Float64}(1, 0, 0, 0), SMatrix{2, 2, Float64}(1 + 1000 * epsT, 0, 0, 1)))
    accepted = count(W -> mg_health_error(() -> check_frame(MEGNO(), W)) === nothing, (SMatrix{2, 2, Float64}(1, 0, 0, 1), SMatrix{2, 2, Float64}(1, 0, 1, 0), SMatrix{2, 2, Float64}(1 + epsT, 0, 0, 1)))
    dmax = maximum(check_frame(MEGNO(), mg_frame(vseed(11, 3), j, Float64); orbit = j) for j in 1:10_000)
    emit("AC11", "(c) frame check: non-unit columns rejected (of 4), unit / dependent / 1+eps accepted (of 3); seeded frames ids 1..1e4 accepted, largest defect $(fmt(dmax)) = $(fmt(dmax / epsT)) eps (threshold 8 eps)", (rejected, accepted), "(4, 3)", pf(rejected == 4 && accepted == 3 && dmax <= 8 * epsT))
    # (d) overflow
    mu = exp(3.0)
    sys = LinearSymplecticMap(SMatrix{2, 2, Float64}(mu, 0, 0, 1 / mu))
    W0 = mg_frame(vseed(11, 5), 1, Float64); x0 = SVector(0.0, 0.0)
    nd = 0; rel = 0.0
    for chain in ((MEGNO(),), (MEGNO(replace = false), SALI()))
        nd += mg_health_error(() -> propagate_orbit(sys, x0, W0, 200, 200, chain)) === nothing ? 0 : 1
        _, _, st = propagate_orbit(sys, x0, W0, 200, 200, chain)
        f = reducer_final(MEGNO(), st[1])
        for i in 1:2
            S = f.sum[i] + f.comp[i]
            rel = max(rel, abs(f.Y[i] - 2S) / abs(2S), abs(f.Ybar[i] - 2S) / abs(2S), abs(f.Ybar20[i] - 3S) / abs(3S))
        end
        e = mg_health_error(() -> propagate_orbit(sys, x0, W0, 300, 300, chain))
        nd += !(e isa NumericalHealthError && e.step == 300 && e.what == "deviation matrix")
    end
    emit("AC11", "(d) diag(e^3, e^-3): tau = 200 raises no error at n = 200 (Y = Ybar = 2 s, Ybar20 = 3 s to relative $(fmt(rel))), tau = 300 raises at handoff 300; failures", nd, "0 (relative error <= 1e-13)", pf(nd == 0 && rel <= 1e-13))
    # (e) run record
    res = propagate(StandardMap(2.0), [SVector(1.0, 2.0)], 8; k = 2, tau = 2, seed = vseed(11, 6), reducers = (MEGNO(),))
    pr, f = res.provenance, res.finals[1][1]
    okrec = pr.indicator == "megno" && pr.convention == MG_CONVENTION && pr.k == 2 && pr.tau == 2 && pr.precision == "Float64" && f.n == 8 && f.handoffs == 4
    emit("AC11", "(e) run record: indicator $(pr.indicator), convention '$(pr.convention)', k $(pr.k), tau $(pr.tau), precision $(pr.precision); final n $(f.n), handoffs $(f.handoffs)", okrec, "true", pf(okrec))
end

# ---------------------------------------------------------------------------------------------------
# AC12: chain, orbit and execution invariance (child processes for the thread count)
section("AC12") do
    seed = vseed(12)
    T = Float64
    bad = bad_b = ncmp = 0
    for tau in (1, 4)
        sys = StandardMap(mg_K(T, :two_pi))
        for (j, s) in enumerate(mg_starts(StableRNG(seed + 1), 5))
            x0 = mg_cast(T, s)
            W2 = mg_frame(seed, j, T)
            W1 = SMatrix{2, 1, T}(W2[:, 1])
            xa, _, sta = propagate_orbit(sys, x0, W1, 400, tau, (mg_Trace(),); orbit = j)
            xb, _, stb = propagate_orbit(sys, x0, W1, 400, tau, (Benettin(), mg_Trace()); orbit = j)
            xc, Wc, stc = propagate_orbit(sys, x0, W2, 400, tau, (MEGNO(), mg_Trace()); orbit = j)
            xd, Wd, std = propagate_orbit(sys, x0, W2, 400, tau, (MEGNO(replace = false), SALI(), mg_Trace()); orbit = j)
            xe, We, ste = propagate_orbit(sys, x0, W2, 400, tau, (FLI(replace = false), MEGNO(replace = false), SALI(), mg_Trace()); orbit = j)
            xf, Wf, stf = propagate_orbit(sys, x0, W2, 400, tau, (MEGNO(replace = false), FLI(replace = false), SALI(), mg_Trace()); orbit = j)
            xg, Wg, stg = propagate_orbit(sys, x0, W2, 400, tau, (SALI(), mg_Trace()); orbit = j)
            xh, _, sth = propagate_orbit(sys, x0, W2, 400, tau, (mg_Trace(),); orbit = j)
            xi, Wi, sti = propagate_orbit(sys, x0, W2, 400, tau, (FLI(), mg_Trace()); orbit = j)
            tra = reducer_final(mg_Trace(), sta[1]).xs
            trs = (reducer_final(mg_Trace(), stb[2]).xs, reducer_final(mg_Trace(), stc[2]).xs, reducer_final(mg_Trace(), std[3]).xs, reducer_final(mg_Trace(), ste[4]).xs,
                   reducer_final(mg_Trace(), stf[4]).xs, reducer_final(mg_Trace(), stg[2]).xs, reducer_final(mg_Trace(), sth[1]).xs)
            bad += !(all(tr -> isequal(tra, tr), trs) && xa === xb && xa === xc && xa === xd && xa === xe && xa === xf && xa === xg && xa === xh)
            m0 = reducer_final(MEGNO(), stc[1]); f0 = reducer_final(FLI(), sti[1]); s0 = reducer_final(SALI(), stg[1])
            ok_b = isequal(m0, reducer_final(MEGNO(replace = false), std[1])) && isequal(m0, reducer_final(MEGNO(replace = false), ste[2])) && isequal(m0, reducer_final(MEGNO(replace = false), stf[1])) &&
                   isequal(f0, reducer_final(FLI(replace = false), ste[1])) && isequal(f0, reducer_final(FLI(replace = false), stf[2])) &&
                   isequal(s0, reducer_final(SALI(), std[2])) && isequal(s0, reducer_final(SALI(), ste[3])) && isequal(s0, reducer_final(SALI(), stf[3])) &&
                   isequal(Wc, Wg) && isequal(Wd, Wg) && isequal(We, Wg) && isequal(Wf, Wg) && isequal(Wi, Wg) &&
                   m0.sum == s0.sum && m0.comp == s0.comp && m0.sum == f0.sum && m0.comp == f0.comp
            for i in 1:2
                _, _, stl = propagate_orbit(sys, x0, SMatrix{2, 1, T}(W2[:, i]), 400, tau, (Benettin(),); orbit = j)
                l = reducer_final(Benettin(), stl[1])
                ok_b &= m0.sum[i] + m0.comp[i] == l.sum[1] + l.comp[1]
            end
            bad_b += !ok_b
            ncmp += 1
        end
    end
    emit("AC12", "(a) orbit points at every handoff differ across the chains {}, {L1} (k=1), {MEGNO}, {MEGNO(f),SALI}, {FLI(f),MEGNO(f),SALI}, {MEGNO(f),FLI(f),SALI}, {SALI}, {} (k=2); cases = (5 starts x tau in {1,4})", bad, "0 of $ncmp", pf(bad == 0))
    emit("AC12", "(b) MEGNO, FLI, SALI finals differ from their stand-alone values, or final W differs among chains with a replacer, or sums differ from SALI / FLI / the L1 run of each column (bitwise)", bad_b, "0 of $ncmp", pf(bad_b == 0))
    for (reducers, taus) in (((MEGNO(replace = false), SALI()), (1, 4)), ((FLI(replace = false), MEGNO(replace = false), SALI()), (4,)), ((MEGNO(),), (4,))), tau in taus
        sys = StandardMap(mg_K(T, :two_pi))
        m, N = 4096, 1000
        x0 = [mg_cast(T, s) for s in mg_starts(StableRNG(seed + 2), m)]
        W0 = [mg_frame(seed, j, T) for j in 1:m]
        serial = [propagate_orbit(sys, x0[j], W0[j], N, tau, reducers; orbit = j) for j in 1:m]
        xs, Ws = [s[1] for s in serial], [s[2] for s in serial]
        fs = [map(reducer_final, reducers, s[3]) for s in serial]
        for chunk in (64, 1000)
            res = propagate(sys, x0, N; k = 2, tau = tau, seed = seed, reducers = reducers, chunk = chunk, W0 = W0)
            okb = isequal(res.x, xs) && isequal(res.W, Ws) && isequal(res.finals, fs)
            emit("AC12", "(b') propagate equals the serial map of propagate_orbit bitwise (x, W, finals), chain of $(length(reducers)) reducers ($(join(map(r -> r isa MEGNO ? (r.replace ? "MEGNO" : "MEGNO(f)") : r isa FLI ? "FLI(f)" : "SALI", reducers), ","))), 4096 orbits, K = 2 pi, n = 1e3, tau = $tau, chunk = $chunk, $(Threads.nthreads()) threads", okb, "true", pf(okb))
        end
    end
    # (c) 1 vs N threads in child processes
    hashes = Dict{Int, Dict{String, String}}(Threads.nthreads() => thread_hashes())
    for nt in (1, 2, 4)
        haskey(hashes, nt) && continue
        cmd = addenv(`$(Base.julia_cmd()) -t $nt --project=$(dirname(Base.active_project())) $(@__FILE__)`, "MG_CHILD" => "1")
        txt = read(pipeline(cmd; stderr = devnull), String)
        d = Dict{String, String}(m[1] => m[2] for m in eachmatch(r"HASH (.+) = (\d+)\n", txt))
        mt = match(r"THREADS (\d+)", txt)
        hashes[nt] = (mt !== nothing && parse(Int, mt.captures[1]) == nt) ? d : Dict("child failed" => first(txt, 100))
    end
    keys_ = sort(collect(keys(hashes[Threads.nthreads()])))
    ok_all = true
    for k in keys_
        hs = [get(hashes[nt], k, "missing") for nt in sort(collect(keys(hashes)))]
        ok = length(unique(hs)) == 1
        ok_all &= ok
        ok || emit("AC12", "(c) hashes differ for $k across threads $(sort(collect(keys(hashes)))): $hs", 0, "equal", "FAIL")
    end
    for cname in ("MEGNO", "MEGNO,SALI", "FLI,MEGNO,SALI"), tau in (1, 4)
        ok_all &= length(unique([hashes[Threads.nthreads()]["$cname tau=$tau chunk=$chunk"] for chunk in (64, 1000)])) == 1
    end
    emit("AC12", "(c) x, W and every final bitwise identical for thread counts $(sort(collect(keys(hashes)))) (hashes from separate processes) and chunk sizes 64 and 1000: $(length(keys_)) cases (chains {MEGNO}, {MEGNO(f),SALI}, {FLI(f),MEGNO(f),SALI} x tau in {1,4} x chunk), 4096 orbits, K = 2 pi, n = 1e3; example hash $(first(values(hashes[1])))", ok_all, "true", pf(ok_all))
    for T in MG_TYPES, tau in (1, 4), reducers in ((MEGNO(),), (MEGNO(replace = false), SALI()))
        a = mg_alloc_bytes(StandardMap(mg_K(T, :two_pi)), mg_cast(T, SVector(1.0, 2.0)), mg_frame(seed, 1, T), 10_000, tau, reducers)
        emit("AC12", "(d) bytes allocated by propagate_orbit after warm-up, n = 1e4, tau = $tau, chain length $(length(reducers)), $T", a, "0", pf(a == 0))
    end
    badx = badn = ncase = 0
    flds = (:Y, :Ybar, :Ybar20, :sum, :comp, :R1, :R2, :Q, :Q2)
    for T in MG_TYPES, label in (:two_pi, :half), tau in (1, 4)
        c = mg_case(label, T, 5)
        W0 = c.W0
        swap(W) = SMatrix{2, 2, T}(W[1, 2], W[2, 2], W[1, 1], W[2, 1])
        neg1(W) = SMatrix{2, 2, T}(-W[1, 1], -W[2, 1], W[1, 2], W[2, 2])
        neg2(W) = SMatrix{2, 2, T}(W[1, 1], W[2, 1], -W[1, 2], -W[2, 2])
        base = mg_run((MEGNO(),), c.sys, c.x0, W0, 400, tau; all = true)
        sw = mg_run((MEGNO(),), c.sys, c.x0, swap.(W0), 400, tau; all = true)
        n1 = mg_run((MEGNO(),), c.sys, c.x0, neg1.(W0), 400, tau; all = true)
        n2 = mg_run((MEGNO(),), c.sys, c.x0, neg2.(W0), 400, tau; all = true)
        for j in 1:5
            a, b = base[j][1].snaps, sw[j][1].snaps
            ncase += 1
            badx += !all(f -> isequal([reverse(getproperty(s, f)) for s in a], [getproperty(s, f) for s in b]), flds)
            for other in (n1, n2)
                badn += !isequal(a, other[j][1].snaps)
            end
        end
    end
    emit("AC12", "(e) exchange of the two columns changing a per-column output not by exchange (bitwise): $badx of $ncase; sign change of either column changing anything: $badn of $(2ncase) (K = 2 pi, 0.5; Float64, Double64; tau in {1,4}; 5 starts; 400 steps, every handoff)", badx + badn, "0", pf(badx + badn == 0))
end

# ---------------------------------------------------------------------------------------------------
# AC13: wrapped versus lifted
section("AC13") do
    seed = vseed(13)
    rng = StableRNG(seed)
    dy() = rand(rng, 0:2048) / 1024
    T = Float64
    nbad = 0
    for tau in (1, 3)
        zs = [SVector{2, T}(T(dy()), T(dy())) for _ in 1:50]
        W0 = [mg_frame(seed, j, T) for j in 1:50]
        rl = mg_run((MEGNO(),), StandardMap(exact(T, 0)), zs, W0, 300, tau; all = true)
        rw = mg_run((MEGNO(),), WrappedStandardMap(exact(T, 0)), torus_state.(zs), W0, 300, tau; all = true)
        nbad += mg_nbad(j -> isequal(rl[j][1].snaps, rw[j][1].snaps), 1:50)
    end
    emit("AC13", "(a) K = 0, 50 dyadic starts, tau in {1,3}, n = 300, every handoff, whole final: orbits with a bitwise difference between wrapped and lifted", nbad, "0 of 100", pf(nbad == 0))
    m = 50
    c = mg_case(:half, T, m)
    for n in (10, 100)
        rl = mg_run((MEGNO(),), c.sys, c.x0, c.W0, n, 1; steps = [n])
        rw = mg_run((MEGNO(),), WrappedStandardMap(0.5), torus_state.(c.x0), c.W0, n, 1; steps = [n])
        excluded = 0; w = 0.0; dmax = 0.0
        for j in 1:m, i in 1:2
            a, b = rl[j][1].snaps[1], rw[j][1].snaps[1]
            rmax, Mst = mg_rho_max(0.5, c.x0[j], c.W0[j][:, i], n)
            if rmax > 50
                excluded += 1
                continue
            end
            B = 6 * 64 * 0.5 * Mst * n^2 * eps(Float64) * rmax
            d = max(abs(a.Y[i] - b.Y[i]), abs(a.Ybar[i] - b.Ybar[i]), abs(a.Ybar20[i] - b.Ybar20[i]))
            dmax = max(dmax, d); w = max(w, d / B)
        end
        ratio_line("AC13", "(b) K = 0.5 regular class, n = $n, 50 starts, both columns, three series: max |M_w - M_l| / (6 * 64 K M_state n^2 eps rho_max) (largest difference $(fmt(dmax)))", w)
        emit("AC13", "(b) n = $n: excluded comparisons (rho_max > 50, of 100)", excluded, "<= 10", pf(excluded <= 10))
    end
end

# ---------------------------------------------------------------------------------------------------
# AC14: ChaosTools one-column tangent system (the only external line: there is no MEGNO routine)
ct_rule(u, p, n) = (pp = u[2] + p[1] * sin(u[1]); SVector(u[1] + pp, pp))
ct_jac(u, p, n) = (c = p[1] * cos(u[1]); SMatrix{2, 2}(1 + c, c, 1, 1))
ct_tds(K) = CT.TangentDynamicalSystem(CT.DeterministicIteratedMap(ct_rule, SVector(0.1, 0.1), [K]); J = ct_jac, k = 1)   # one per call (threads)
# One column w0 from z0: steps the oracle, reads the deviation after each step, takes the log norm in BigFloat(256), then renormalises the deviation (Float64).
# Returns the per-step log increments a_k (BigFloat) and the number of steps at which the oracle orbit differs bitwise from forward(StandardMap).
function ct_series(K::Float64, z0, w0, n::Int)
    tds = ct_tds(K)
    CT.reinit!(tds, SVector{2, Float64}(z0); Q0 = SMatrix{2, 1, Float64}(w0[1], w0[2]))
    sys = StandardMap(K)
    z = SVector{2, Float64}(z0)
    a = Vector{BigFloat}(undef, n)
    mism = 0
    setprecision(BigFloat, MG_BITS) do
        for i in 1:n
            CT.step!(tds)
            w = SVector{2, Float64}(CT.current_deviations(tds)[:, 1])
            a[i] = log(BigFloat(w[1])^2 + BigFloat(w[2])^2) / 2
            CT.set_deviations!(tds, SMatrix{2, 1, Float64}(w[1] / norm(w), w[2] / norm(w)))
            z = forward(sys, z)
            mism += SVector{2, Float64}(CT.current_state(tds)) == z ? 0 : 1
        end
    end
    return a, mism
end
psimax_cp(Ms, w0, cps) = accumulate(max, [lb_Psi(view(Ms, 1:g), w0) for g in cps])

section("AC14") do
    DATA[:ct] = NamedTuple[]
    for (label, m, ns) in ((:half, 20, [10, 100, 1000, 10_000]), (:two_pi, 1000, [10, 100, 1000]), (:five, 1000, [10, 100, 1000]))
        c = mg_case(label, Float64, m)
        K = Float64(c.K)
        N = last(ns)
        cps = collect(10:10:N)
        runs = mg_run((MEGNO(),), c.sys, c.x0, c.W0, N, 1; steps = cps)
        res = Vector{Any}(undef, m)
        Threads.@threads :dynamic for j in 1:m
            Ms = lb_tangents(c.sys, c.x0[j], N)
            out = []
            for i in 1:2
                a, mism = ct_series(K, c.x0[j], c.W0[j][:, i], N)
                ref = mg_M5(cumsum(a), 1)
                pm = psimax_cp(Ms, c.W0[j][:, i], cps)
                push!(out, (ref = ref, mism = mism, pm = pm))
            end
            res[j] = out
        end
        mism_total = sum(res[j][i].mism for j in 1:m for i in 1:2)
        nbit = count(j -> all(res[j][i].mism == 0 for i in 1:2), 1:m)
        emit("AC14", "K = $K ($label): oracle orbit vs forward(StandardMap) over $m starts x $N steps x 2 columns: steps with a bitwise difference; starts with a fully bitwise oracle orbit", "$mism_total; $nbit of $m", "0 (else excluded, reported)", pf(mism_total == 0))
        tot = excl = 0
        w = 0.0; wlist = Float64[]
        for j in 1:m, n in ns
            tot += 1
            if !all(res[j][i].mism == 0 for i in 1:2)
                excl += 1
                continue
            end
            snaps = runs[j][1].snaps
            s = snaps[findfirst(s -> s.n == n, snaps)]
            idx = n ÷ 10
            for i in 1:2
                BM = 48 * EPS64 * res[j][i].pm[idx]
                BO = BM                      # same orbit bitwise: the oracle's own tangents are ours; B_oracle is the M11 form with eps(Float64)
                ref = res[j][i].ref
                d = max(mg_err(s.Y[i], ref.Y[n]), mg_err(s.Ybar[i], ref.Ybar[n]), mg_err(s.Ybar20[i], ref.Ybar20[n]))
                w = max(w, d / (BM + BO))
                push!(wlist, d / (BM + BO))
                push!(DATA[:ct], (K = K, n = n, Y = Float64(s.Y[i]), Yref = Float64(ref.Y[n]), ratio = d / (BM + BO)))
            end
        end
        ratio_line("AC14", "K = $K ($label), $m starts, n in $ns, both columns, three series: max |M_Chaoslab - M_oracle| / (B_M + B_oracle), $(2 * (tot - excl)) pairs; median ratio $(fmt(mg_median(wlist)))", w)
        emit("AC14", "K = $K ($label): pairs excluded (oracle orbit not bitwise equal)", excl / tot, "<= 0.10", pf(excl / tot <= 0.10))
    end
end

# ---------------------------------------------------------------------------------------------------
# Figures
section("figure") do
    # panel 1: regular class, Ybar_1 and 4 Ybar20_1 to the levels 2 and 2
    c = mg_case(:half, Float64, 8)
    nn = unique(round.(Int, 10 .^ range(1, 5; length = 41)))
    run = mg_run((MEGNO(),), c.sys, c.x0, c.W0, 100_000, 1; steps = nn)
    fig = CM.Figure(size = (1800, 520))
    ax1 = CM.Axis(fig[1, 1]; xscale = log10, xlabel = "n", ylabel = "Ybar_1(n)", title = "regular class, K = 0.5: Ybar_1 -> 2 (M3, M8(b))")
    for r in run
        CM.lines!(ax1, Float64.(nn), Float64[s.Ybar[1] for s in r[1].snaps]; linewidth = 1)
    end
    CM.hlines!(ax1, [2.0]; color = :black, linestyle = :dash, label = "2")
    CM.axislegend(ax1; position = :rb, labelsize = 10)
    # panel 2: chaotic, Ybar_1 and lambda_1 n / 2
    ax2 = CM.Axis(fig[1, 2]; xlabel = "n", ylabel = "Ybar_1(n)", title = "chaotic: ensemble mean of Ybar_1 (200 starts) and lambda_1 n / 2")
    for (label, col) in ((:two_pi, :darkorange), (:ten, :purple))
        sys = StandardMap(mg_K(Float64, label)); seed = vseed(20, label === :two_pi ? 1 : 2)
        x0 = [mg_cast(Float64, s) for s in mg_starts(StableRNG(seed), 200)]; W0 = [mg_frame(seed, j, Float64) for j in 1:200]
        nl = collect(1:300)
        r = mg_run((MEGNO(),), sys, x0, W0, 300, 1; steps = nl)
        l1 = mg_run((Benettin(),), sys, x0, [SMatrix{2, 1, Float64}(W[:, 1]) for W in W0], 10_000, 1; steps = [10_000])
        lam = mean((q[1].snaps[1].sum[1] + q[1].snaps[1].comp[1]) / 10_000 for q in l1)
        CM.lines!(ax2, Float64.(nl), [mean(Float64(q[1].snaps[i].Ybar[1]) for q in r) for i in eachindex(nl)]; color = col, label = "K = $(fmt(Float64(mg_K(Float64, label)))): mean Ybar_1")
        CM.lines!(ax2, Float64.(nl), lam .* nl ./ 2; color = col, linestyle = :dash, label = "lambda_1 n / 2, lambda_1 = $(fmt(lam))")
    end
    CM.axislegend(ax2; position = :lt, labelsize = 9)
    # panel 3: Mestre Eq. (A2) residual of Y_1 in units of its tolerance
    ax3 = CM.Axis(fig[1, 3]; xscale = log10, yscale = log10, xlabel = "n", ylabel = "|Y_1 - (2 (FLI - FLIbar) + 2 FLI / n)| / (16 eps max(1, max |s_1|))", title = "Mestre2011 Eq. (A2) residual (bound 1)")
    worst = 0.0
    for (label, col) in ((:half, :seagreen), (:two_pi, :darkorange), (:ten, :purple))
        cc = mg_case(label, Float64, 3)
        N = 10_000
        cps = collect(10:10:N)
        rr = mg_run((MEGNO(),), cc.sys, cc.x0, cc.W0, N, 1; all = true)
        for j in 1:3
            sn = rr[j][1].snaps
            S1 = [s.sum[1] + s.comp[1] for s in sn]
            smax = accumulate(max, abs.(S1))
            pre = setprecision(() -> accumulate(+, BigFloat.(S1)), BigFloat, MG_BITS)
            res = [setprecision(BigFloat, MG_BITS) do
                       a2 = 2 * (BigFloat(S1[n]) - pre[n] / n) + 2 * BigFloat(S1[n]) / n
                       Float64(abs(BigFloat(sn[n].Y[1]) - a2)) / (16 * EPS64 * max(1, smax[n]))
                   end for n in cps]
            worst = max(worst, maximum(res))
            CM.scatter!(ax3, Float64.(cps), max.(res, 1e-4); color = (col, 0.5), markersize = 3, label = j == 1 ? "K = $(fmt(Float64(mg_K(Float64, label))))" : nothing)
        end
    end
    CM.hlines!(ax3, [1.0]; color = :black, linestyle = :dash)
    ax3.title = "Mestre Eq. (A2) residual / tolerance, worst $(fmt(worst))"
    CM.axislegend(ax3; position = :rb, labelsize = 9)
    CM.save(joinpath(OUT, "megno-growth-laws.png"), fig)
    say("# figure saved: megno-growth-laws.png (panel 3 worst ratio $(fmt(worst)))")
    haskey(DATA, :ct) || return
    pts = DATA[:ct]
    fig2 = CM.Figure(size = (1200, 480))
    ax4 = CM.Axis(fig2[1, 1]; xlabel = "Y_1 or Y_2 from the ChaosTools one-column tangent system (BigFloat M5)", ylabel = "Y (Chaoslab)", title = "AC14: Y, Chaoslab against the oracle line")
    for (K, col) in ((0.5, :seagreen), (Float64(mg_K(Float64, :two_pi)), :darkorange), (5.0, :purple))
        sel = [p for p in pts if p.K == K]
        CM.scatter!(ax4, [p.Yref for p in sel], [p.Y for p in sel]; color = (col, 0.4), markersize = 4, label = "K = $(fmt(K))")
    end
    mx = maximum(p.Yref for p in pts)
    CM.lines!(ax4, [minimum(p.Yref for p in pts), mx], [minimum(p.Yref for p in pts), mx]; color = :black)
    CM.axislegend(ax4; position = :lt, labelsize = 10)
    ax5 = CM.Axis(fig2[1, 2]; yscale = log10, xlabel = "n", ylabel = "|M_Chaoslab - M_oracle| / (B_M + B_oracle)", title = "AC14: ratio to the bound (bound 1)")
    for (K, col) in ((0.5, :seagreen), (Float64(mg_K(Float64, :two_pi)), :darkorange), (5.0, :purple))
        sel = [p for p in pts if p.K == K]
        CM.scatter!(ax5, [Float64(p.n) for p in sel], [max(p.ratio, 1e-6) for p in sel]; color = (col, 0.4), markersize = 4, label = "K = $(fmt(K))")
    end
    CM.hlines!(ax5, [1.0]; color = :black, linestyle = :dash)
    ax5.xscale = log10
    CM.axislegend(ax5; position = :rb, labelsize = 9)
    CM.save(joinpath(OUT, "megno-chaostools.png"), fig2)
    say("# figure saved: megno-chaostools.png")
end

say(@sprintf("# verdicts: %d PASS, %d FAIL, %d other", count(==("PASS"), VERDICTS), count(==("FAIL"), VERDICTS), count(v -> v != "PASS" && v != "FAIL", VERDICTS)))
say(@sprintf("# total wall time %.1f s", time() - T0))
close(LOG)
