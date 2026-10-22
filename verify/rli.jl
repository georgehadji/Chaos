# Stage 04 verification of module rli (reducer RLI, PairMap, pair_start, pair_distance, check_frame; layers Indicators and Systems; catalog row I9).
# Run from the repo root:  julia -t 2 --project=verify verify/rli.jl
# Spec: workspaces/module-build/stages/01-spec/output/rli-spec.md (criteria 1 to 15, Amendment D76). Criterion 15 (CPU vs GPU) is in verify/gpu.jl, tag "rli".
# Section "suite": the package suite at full size (RL_FULL_RUN = true: 20 starts to n = 1e4 in AC3, n = 1e4 in AC4, 1000 starts in AC5 at 1e4, the spec sizes in AC10(a), n = 1e5 in AC2
# for Double64, every multiple of 10 for the Psi grid); the @info lines of the suite (the realized ratios R2, the recorded lines R1, R3, R8) are forwarded into the log by a logger.
# Other sections: one printed line per measurement  AC<n> | quantity | value | bound | PASS/FAIL/INFO  with the MEASURED number, reusing the helpers of test/rli/rl_helpers.jl
# (one home for B_R, Psi_max, the BigFloat reference, the RL4 post-processing). Values are "worst ratio error / bound" (pass when <= 1) unless stated. The cases use verify seeds
# (RL_VSEED0 + 50 ac + sub), distinct from the package tests (RL_SEED + 36000 ..).
# RL_ONLY="AC14,figure" runs only the named sections (debugging aid) and RL_LOG names the log file; the default runs everything.
# RL_CHILD=1 is the child-process mode of the thread-invariance check of AC12(b).
using Printf, Test, Random, LinearAlgebra, Statistics, Logging
using Chaoslab, StaticArrays, DoubleFloats, StableRNGs

const RL_VSEED0 = 20261012 + 59000      # verify offset (megno used 58000)

# ---------------------------------------------------------------------------------------------------
# Thread-count invariance (AC12b): the same outputs hashed in this process and in child processes with other thread counts.
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
    seed = RL_VSEED0 + 50 * 12 + 3
    base = WrappedStandardMap(2pi)
    P = PairMap(base)
    chains = (("RLI", (RLI(),)), ("RLI(f),MEGNO", (RLI(replace = false), MEGNO())))
    m, N = 4096, 1000
    for (cname, reducers) in chains, tau in (1, 4), chunk in (64, 1000)
        rng = StableRNG(seed)
        z0 = Vector{SVector{8, Float64}}(undef, m)
        W0 = Vector{SMatrix{4, 2, Float64, 8}}(undef, m)
        for j in 1:m
            x0 = SVector(2pi * rand(rng), 2pi * rand(rng))
            xi = SVector{2, Float64}(initial_deviations(seed, j, Val(2), Val(1), Float64)[:, 1])
            z0[j] = pair_start(base, x0, xi, 1e-9; orbit = j).z
            W0[j] = SMatrix{4, 2, Float64}(xi[1], xi[2], 0.0, 0.0, 0.0, 0.0, xi[1], xi[2])
        end
        res = propagate(P, z0, N; k = 2, tau = tau, seed = seed, reducers = reducers, chunk = chunk, W0 = W0, backend = :cpu)
        v = Float64[]
        v_flatten!(v, res.x); v_flatten!(v, res.W); v_flatten!(v, res.finals)
        out["$cname tau=$tau chunk=$chunk"] = string(hash(reinterpret(UInt64, v)))
    end
    return out
end
if get(ENV, "RL_CHILD", "") == "1"
    for (k, h) in sort!(collect(thread_hashes())); println("HASH ", k, " = ", h); end
    println("THREADS ", Threads.nthreads())
    exit(0)
end

const RL_FULL_RUN = true
const ROOT = dirname(@__DIR__)
const TESTDIR = joinpath(ROOT, "test", "rli")
include(joinpath(TESTDIR, "rl_helpers.jl"))

vseed(ac, sub = 0) = RL_VSEED0 + 50 * ac + sub

import ChaosTools as CT
import CairoMakie as CM

const OUT = joinpath(ROOT, "workspaces", "module-build", "stages", "04-verify", "output")
mkpath(OUT)
const ONLY = filter(!isempty, split(get(ENV, "RL_ONLY", ""), ","))
const LOG = open(get(ENV, "RL_LOG", joinpath(OUT, "rli-verify-log.txt")), "w")
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
rmedian(v) = sort(v)[max((length(v) + 1) ÷ 2, 1)]

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
const EPS64 = rl_eps(Float64)
say("# rli verification, julia $(VERSION), threads $(Threads.nthreads()), base seed $RL_SEED (verify offset 59000), time $(round(time())) (unix)")
say("# timings in this log were taken on a machine that was not otherwise idle")

# ---------------------------------------------------------------------------------------------------
# The package suite at full size; its @info lines are forwarded into this log (the gated ratio keys are listed in GATES)
const GATES = ("AC1 realized" => ((:worst, 1.0),), "AC3 realized" => ((:worst, 1.0), (:worstD, 1.0)), "AC4(b)" => ((:worst, 1e-3),),
               "AC8(c)" => ((:worstD, 1.0), (:worstR, 1.0)), "AC10(b)" => ((:worst, 1.0),))
struct SuiteLogger <: Logging.AbstractLogger end
Logging.min_enabled_level(::SuiteLogger) = Logging.Info
Logging.shouldlog(::SuiteLogger, args...) = true
Logging.catch_exceptions(::SuiteLogger) = true
function Logging.handle_message(::SuiteLogger, level, message, _module, group, id, file, line; kwargs...)
    kw = Dict{Symbol, Any}(kwargs)
    msg = string(message)
    gates = ()
    for (pre, g) in GATES
        startswith(msg, pre) && (gates = g)
    end
    gated = [(k, b) for (k, b) in gates if haskey(kw, k)]
    rest = join(("$k = $(fmt(v))" for (k, v) in kw if !any(k == g[1] for g in gated)), "; ")
    for (k, b) in gated
        v = kw[k] isa Real ? Float64(kw[k]) : maximum(Float64.(kw[k]))
        emit("SUITE", string(msg, " [", k, "] ", rest), v, "<= $b", pf(v <= b))
    end
    isempty(gated) && emit("SUITE", string(msg, " ", rest), "-", "reported", "INFO")
end

section("suite") do
    local ts
    t = @elapsed try
        ts = with_logger(SuiteLogger()) do
            @testset verbose = true "rli verify" begin
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
    emit("SUITE", "package suite at full size (RL_FULL_RUN), assertions passed (fail + error = $nbad); wall time $(round(t; digits = 1)) s", npass, "0 fail, 0 error", pf(ok && nbad == 0 && npass > 0))
end

# ---------------------------------------------------------------------------------------------------
# AC4(a) diagnosis at the full size: every (start, n) of the BigFloat(256) own-orbit run whose O(Dx) law test or slope test of the suite fails is listed with its numbers
# (rl_own, rl_regular_law are the suite's own functions, defined in Main by the include above).
section("AC4diag") do
    DxA = (1e-14, 1e-12, 1e-10, 1e-8, 1e-6); ns = (100, 1000, 10_000)
    bf = rl_own(BigFloat, :half, 20, DxA, ns)
    nlaw = nslope = 0; minslope, maxslope = Inf, -Inf
    for j in 1:20, ni in eachindex(ns)
        r = rl_regular_law(bf, j, ni, DxA)
        minslope, maxslope = min(minslope, r.slope), max(maxslope, r.slope)
        nslope += !(0.99 <= r.slope <= 1.01)
        if !r.law_ok
            nlaw += 1
            dxr = [bf[(j, Float64(D))].Dxr for D in DxA]
            rhoT = [bf[(j, Float64(D))].rli[ni] / d for (D, d) in zip(DxA, dxr)]
            dev = [Float64(abs(rhoT[i] / rhoT[1] - 1)) for i in 2:5]
            info("AC4diag", "law test dev_i <= 2 (Dx_i / Dx_{i+1}) dev_{i+1} FAILS at start $j, n = $(ns[ni]): dev at Dx = 1e-12, 1e-10, 1e-8, 1e-6 = $(join(fmt.(dev), ", ")); slope $(fmt(r.slope)); max dev_i / Dx_i = $(fmt(r.coeff)); rho at Dx = 1e-14 = $(fmt(Float64(rhoT[1])))", "-")
        end
    end
    info("AC4diag", "AC4(a), 20 starts x n in {100, 1e3, 1e4} = 60 cases: slope tests failing $nslope (slope min $(fmt(minslope)), max $(fmt(maxslope))); O(Dx) law tests failing $nlaw", "$nslope; $nlaw")
end

# ---------------------------------------------------------------------------------------------------
# AC14: ChaosTools two one-column tangent systems (base and shadow, the same supplied xi_0): the only external line (there is no RLI routine, spec Oracle).
# Lifted standard map (the oracle's rule is the lifted one, as in the sali and megno lines), Dx = 1e-12 along xi_0, tau = 1.
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
    setprecision(BigFloat, RL_BITS) do
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

section("AC14") do
    DATA[:ct] = NamedTuple[]
    Dx = 1e-12
    for (label, m, ns) in ((:half, 20, [10, 100, 1000, 10_000]), (:two_pi, 1000, [10, 100, 1000]), (:five, 1000, [10, 100, 1000]))
        c = rl_case(label, Float64, m; wrapped = false)
        K = Float64(c.K)
        N = last(ns)
        cps = collect(10:10:N)
        prs = rl_pairs(c, Dx)
        P = PairMap(c.sys)
        res = Vector{Any}(undef, m)
        Threads.@threads :dynamic for j in 1:m
            snaps = rl_run((RLI(),), P, [prs[j].z], [rl_frame(c.xi0[j])], N, 1; steps = ns, ids = [c.ids[j]])[1][1].snaps
            aA, mA = ct_series(K, prs[j].zA, c.xi0[j], N)
            aB, mB = ct_series(K, prs[j].zB, c.xi0[j], N)
            sA, sB = cumsum(aA), cumsum(aB)
            ref = rl_ref_rli(sA, sB, 1)
            o = rl_PairRef(sA, sB, lb_tangents(c.sys, prs[j].zA, N), lb_tangents(c.sys, prs[j].zB, N), c.xi0[j], Dict{Tuple{Int, Int}, Float64}())
            rows = NamedTuple[]
            for (f, n) in zip(snaps, ns)
                f.n == n || error("AC14: snapshot step $(f.n) vs $n")
                PA, PB = rl_Psimax(o, 1, n, cps), rl_Psimax(o, 2, n, cps)
                BR = rl_BR(EPS64, PA, PB, n, n, Float64(ref.RLI[n]))
                BD = rl_BDL(EPS64, rl_Psi(o, 1, n), rl_Psi(o, 2, n), n, Float64(ref.DL[n]))
                push!(rows, (n = n, RLI = f.RLI, ref = Float64(ref.RLI[n]), rR = rat(rl_err(f.RLI, ref.RLI[n]), 2BR), rD = rat(rl_err(f.DL, ref.DL[n]), 2BD)))
            end
            res[j] = (mism = mA + mB, rows = rows)
        end
        mism_total = sum(r.mism for r in res)
        nbit = count(r -> r.mism == 0, res)
        emit("AC14", "K = $K ($label): oracle orbits (base and shadow) vs forward(StandardMap) over $m pairs x $N steps x 2 orbits: steps with a bitwise difference; pairs with both oracle orbits fully bitwise", "$mism_total; $nbit of $m", "0 (else excluded, reported)", pf(mism_total == 0))
        tot = excl = 0
        wR = wD = 0.0; wlist = Float64[]
        for r in res, row in r.rows
            tot += 1
            if r.mism != 0
                excl += 1
                continue
            end
            wR = max(wR, row.rR); wD = max(wD, row.rD); push!(wlist, row.rR)
            push!(DATA[:ct], (K = K, n = row.n, RLI = row.RLI, ref = row.ref, ratio = row.rR))
        end
        ratio_line("AC14", "K = $K ($label), $m pairs, n in $ns: max |RLI_Chaoslab - RLI_oracle| / (B_R + B_oracle), $(tot - excl) comparisons; median ratio $(fmt(rmedian(wlist)))", wR)
        ratio_line("AC14", "K = $K ($label): same for DL_n over (B_DL + B_oracle)", wD)
        emit("AC14", "K = $K ($label): comparisons excluded (oracle orbit not bitwise equal)", excl / tot, "<= 0.10", pf(excl / tot <= 0.10))
    end
end

# ---------------------------------------------------------------------------------------------------
# AC12(b): 1 vs N threads in child processes (output hashes), with the in-process chunk-size comparison
section("AC12") do
    hashes = Dict{Int, Dict{String, String}}(Threads.nthreads() => thread_hashes())
    for nt in (1, 2, 4)
        haskey(hashes, nt) && continue
        cmd = addenv(`$(Base.julia_cmd()) -t $nt --project=$(dirname(Base.active_project())) $(@__FILE__)`, "RL_CHILD" => "1")
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
        ok || emit("AC12", "(b) hashes differ for $k across threads $(sort(collect(keys(hashes)))): $hs", 0, "equal", "FAIL")
    end
    for cname in ("RLI", "RLI(f),MEGNO"), tau in (1, 4)
        ok_all &= length(unique([hashes[Threads.nthreads()]["$cname tau=$tau chunk=$chunk"] for chunk in (64, 1000)])) == 1
    end
    emit("AC12", "(b) x, W and every final bitwise identical for thread counts $(sort(collect(keys(hashes)))) (hashes from separate processes) and chunk sizes 64 and 1000: $(length(keys_)) cases (chains {RLI}, {RLI(f),MEGNO} x tau in {1,4} x chunk), 4096 pairs, K = 2 pi, wrapped, Dx = 1e-9, n = 1e3; example hash $(first(values(hashes[1])))", ok_all, "true", pf(ok_all))
end

# ---------------------------------------------------------------------------------------------------
# R4: threshold demonstration of Maffione2011a Table 1 (RLI above 1e-12 at Dx = 1e-12, scaled with Dx_r) against the labels of L1, FLI and MEGNO
section("R4") do
    m = 500
    Dx = 1e-12
    for (label, K) in ((:half, 0.5), (:two_pi, 2pi), (:ten, 10.0))
        base = WrappedStandardMap(K)
        P = PairMap(base)
        seed = vseed(30, label === :half ? 1 : label === :two_pi ? 2 : 3)
        rng = StableRNG(seed)
        x0 = [SVector(2pi * rand(rng), 2pi * rand(rng)) for _ in 1:m]
        xis = [SVector{2, Float64}(initial_deviations(seed, j, Val(2), Val(1), Float64)[:, 1]) for j in 1:m]
        ps = [pair_start(base, x0[j], xis[j], Dx; orbit = j) for j in 1:m]
        z0 = [p.z for p in ps]
        xw = [torus_state(x) for x in x0]
        W0 = [rl_frame(x) for x in xis]
        for n in (1_000, 10_000, 100_000)
            r = propagate(P, z0, n; k = 2, tau = 1, seed = seed, reducers = (RLI(),), W0 = W0, backend = :cpu)
            rli = [r.finals[j][1].RLI / ps[j].Dx_r for j in 1:m]                 # rho = RLI / Dx_r; threshold 1e-12 at Dx = 1e-12 is rho > 1 (scaled with Dx_r)
            l1 = propagate(base, xw, n; k = 1, tau = 1, seed = seed, reducers = (Benettin(),), backend = :cpu)
            l1s = [f[1].lambda[1] * n for f in l1.finals]                         # s(n) = n X(n), X(n) = lambda_1
            fm = propagate(base, xw, n; k = 2, tau = 1, seed = seed, reducers = (FLI(replace = false), MEGNO()), backend = :cpu)
            fli = [f[1].fli for f in fm.finals]; yb = [f[2].Ybar[1] for f in fm.finals]
            lab_rli = rli .> 1.0
            lab_l1 = l1s .> log(n); lab_fli = fli .> log(n); lab_mg = yb .> 3
            info("R4", "K = $(fmt(K)), $m uniform starts, n = $n, Dx = 1e-12: fraction labelled chaotic by RLI (rho = RLI / Dx_r > 1, i.e. RLI above 1e-12 scaled), L1 (n X > ln n), FLI_B (> ln n), MEGNO (Ybar_1 > 3); disagreements of RLI with L1 / FLI / MEGNO",
                 "$(fmt(mean(lab_rli))), $(fmt(mean(lab_l1))), $(fmt(mean(lab_fli))), $(fmt(mean(lab_mg))); $(count(lab_rli .!= lab_l1)) / $(count(lab_rli .!= lab_fli)) / $(count(lab_rli .!= lab_mg)); rho median $(fmt(rmedian(rli))), max $(fmt(maximum(rli)))")
        end
    end
end

# ---------------------------------------------------------------------------------------------------
# R5 at the spec size: 20 regular-class starts x 20 supplied unit vectors xi_0 x 20 directions u, n = 1e3 and 1e4 (Dx = 1e-9, wrapped, tau = 1)
section("R5") do
    c = rl_case(:half, Float64, 20)
    P = PairMap(c.sys)
    rng = StableRNG(vseed(31, 1))
    z0 = SVector{8, Float64}[]; W0 = SMatrix{4, 2, Float64, 8}[]; Dxr = Float64[]; start_of = Int[]
    nbadframe = 0
    for j in 1:20, _ in 1:20, _ in 1:20
        a, b = 2pi * rand(rng), 2pi * rand(rng)
        xi, u = SVector(cos(a), sin(a)), SVector(cos(b), sin(b))
        W = rl_frame(xi)
        nbadframe += rl_health_error(() -> check_frame(RLI(), W; orbit = j)) === nothing ? 0 : 1
        p = pair_start(c.sys, c.x0[j], u, 1e-9; orbit = j)
        push!(z0, p.z); push!(W0, W); push!(Dxr, p.Dx_r); push!(start_of, j)
    end
    emit("R5", "frames rejected by check_frame (of 8000)", nbadframe, "0", pf(nbadframe == 0))
    for n in (1_000, 10_000)
        r = propagate(P, z0, n; k = 2, tau = 1, seed = vseed(31, 2), reducers = (RLI(),), W0 = W0, backend = :cpu)
        rho = [r.finals[i][1].RLI / Dxr[i] for i in eachindex(z0)]
        per = [rho[findall(==(j), start_of)] for j in 1:20]
        spread = [maximum(v) / minimum(v) for v in per]
        info("R5", "n = $n, 20 starts x 20 xi_0 x 20 directions (8000 pairs): rho = RLI / Dx_r min / median / max over all; per-start max/min ratio (median; max over starts); finite and nonnegative all",
             "$(fmt(minimum(rho))) / $(fmt(rmedian(rho))) / $(fmt(maximum(rho))); $(fmt(rmedian(spread))); $(fmt(maximum(spread))); $(all(isfinite, rho) && all(>=(0), rho))")
    end
end

# ---------------------------------------------------------------------------------------------------
# R6: the first n at which the Float64 RLI exceeds 10 times the Float64 floor of RL11 (8 eps (Psi_max^A + Psi_max^B) H_K / n), ordered start of R1 (nu = 0.5, (q, y) = (2, 0))
section("R6") do
    N = 20_000
    wsys = WrappedStandardMap(-0.5)
    P = PairMap(wsys)
    grid = collect(100:100:N)
    for Dx in (1e-14, 1e-12, 1e-10)
        xi = rl_xi0(rl_seed(16, 1), 1, Float64)
        pr = rl_pair(wsys, SVector(2.0, 0.5 * sin(2.0)), xi, Dx)
        snaps = rl_run((RLI(),), P, [pr.z], [rl_frame(xi)], N, 1; steps = collect(10:10:N))[1][1].snaps
        sA = rl_ref_s(rl_mats_std(rl_K(Float64, :half) * -1, rl_pts(wsys, pr.zA, N)), xi, N)
        sB = rl_ref_s(rl_mats_std(rl_K(Float64, :half) * -1, rl_pts(wsys, pr.zB, N)), xi, N)
        o = rl_PairRef(sA, sB, lb_tangents(wsys, pr.zA, N), lb_tangents(wsys, pr.zB, N), xi, Dict{Tuple{Int, Int}, Float64}())
        ref = rl_ref_rli(sA, sB, 1)
        first_n = nothing; first_bad = nothing
        for f in snaps
            n = f.n
            first_bad === nothing && abs(f.RLI - Float64(ref.RLI[n])) > 0.1 * Float64(ref.RLI[n]) && (first_bad = n)
            floor_n = 8 * EPS64 * (rl_Psimax(o, 1, max(n, 100), grid) + rl_Psimax(o, 2, max(n, 100), grid)) * rl_H(n) / n
            if f.RLI > 10 * floor_n && first_n === nothing
                first_n = n
            end
        end
        n_last = last(snaps).n
        info("R6", "ordered start (2, 0) of R1, nu = 0.5, Dx = $Dx (Dx_r = $(fmt(pr.Dxr))): first n (multiples of 10) with RLI(n) > 10 x the RL11 floor 8 eps (Psi_max^A + Psi_max^B) H_K / n; RLI(2e4) = $(fmt(last(snaps).RLI)), BigFloat $(fmt(Float64(ref.RLI[N]))), relative error $(fmt(Float64(abs(last(snaps).RLI - ref.RLI[N]) / ref.RLI[N])))",
             first_n === nothing ? "never within 2e4" : first_n)
        info("R6", "same start and Dx = $Dx: first n (multiples of 10) with the measured Float64 relative error |RLI - RLI_BigFloat| / RLI_BigFloat above 0.1 (the empirical resolution horizon; the worst-case floor above is pessimistic, F5)", first_bad === nothing ? "never within 2e4" : first_bad)
    end
end

# ---------------------------------------------------------------------------------------------------
# R7 (CPU): wall time of the pair run against a single-orbit k = 1 and k = 2 run on the same starts (GPU: verify/gpu.jl)
section("R7") do
    m, N = 20_000, 1_000
    rng = StableRNG(vseed(32, 1))
    x0 = [SVector(2pi * rand(rng), 2pi * rand(rng)) for _ in 1:m]
    base = WrappedStandardMap(2pi)
    P = PairMap(base)
    xis = [SVector{2, Float64}(initial_deviations(vseed(32, 2), j, Val(2), Val(1), Float64)[:, 1]) for j in 1:m]
    z0 = [pair_start(base, x0[j], xis[j], 1e-9; orbit = j).z for j in 1:m]
    W0 = [rl_frame(x) for x in xis]
    tbest(f) = (f(); minimum(@elapsed(f()) for _ in 1:3))
    tp = tbest(() -> propagate(P, z0, N; k = 2, tau = 1, seed = 1, reducers = (RLI(),), W0 = W0, backend = :cpu))
    t1 = tbest(() -> propagate(base, [torus_state(x) for x in x0], N; k = 1, tau = 1, seed = 1, reducers = (Benettin(),), backend = :cpu))
    t2 = tbest(() -> propagate(base, [torus_state(x) for x in x0], N; k = 2, tau = 1, seed = 1, reducers = (MEGNO(),), backend = :cpu))
    st = m * N
    info("R7", "CPU ($(Threads.nthreads()) threads), $m orbits x $N steps, K = 2 pi, wrapped, tau = 1: pair RLI (k = 2, d = 4) $(fmt(tp / st * 1e9)) ns/orbit-step; single L1 (k = 1) $(fmt(t1 / st * 1e9)); single MEGNO (k = 2) $(fmt(t2 / st * 1e9)); ratio pair / L1 $(fmt(tp / t1)), pair / MEGNO $(fmt(tp / t2))", "-")
end

# ---------------------------------------------------------------------------------------------------
# Figures
section("figure") do
    ns = unique(round.(Int, 10 .^ range(1, log10(20_000); length = 61)))
    wsys = WrappedStandardMap(-0.5)
    P = PairMap(wsys)
    Dxs = (1e-12, 1e-10, 1e-8)
    cols = (:royalblue, :darkorange, :seagreen)
    runs = Dict{Tuple{Symbol, Float64}, Any}()
    for (nm, q0, id) in ((:ordered, 2.0, 1), (:chaotic, 3.0, 2)), Dx in Dxs
        xi = rl_xi0(rl_seed(16, 1), id, Float64)
        pr = rl_pair(wsys, SVector(q0, 0.5 * sin(q0)), xi, Dx)
        snaps = rl_run((RLI(),), P, [pr.z], [rl_frame(xi)], 20_000, 1; steps = ns)[1][1].snaps
        runs[(nm, Dx)] = (n = Float64.(ns), RLI = Float64[s.RLI for s in snaps], Dxr = pr.Dxr)
    end
    fig = CM.Figure(size = (1800, 560))
    ax1 = CM.Axis(fig[1, 1]; xscale = log10, yscale = log10, xlabel = "n", ylabel = "RLI(n)", title = "standard map (SandorMaffione2016 Eq. (2), nu = 0.5): ordered (2,0) solid, chaotic (3,0) dashed")
    for (Dx, col) in zip(Dxs, cols), (nm, ls) in ((:ordered, :solid), (:chaotic, :dash))
        r = runs[(nm, Dx)]
        CM.lines!(ax1, r.n, max.(r.RLI, 1e-30); color = col, linestyle = ls, label = "$nm Dx = $(fmt(Dx))")
    end
    CM.axislegend(ax1; position = :rb, labelsize = 9)
    ax2 = CM.Axis(fig[1, 2]; xscale = log10, yscale = log10, xlabel = "n", ylabel = "rho = RLI / Dx_r", title = "rho = RLI / Dx_r: ordered curves collapse (RL6), chaotic rho scales as 1/Dx")
    for (Dx, col) in zip(Dxs, cols), (nm, ls) in ((:ordered, :solid), (:chaotic, :dash))
        r = runs[(nm, Dx)]
        CM.lines!(ax2, r.n, max.(r.RLI ./ r.Dxr, 1e-30); color = col, linestyle = ls, label = "$nm Dx = $(fmt(Dx))")
    end
    CM.axislegend(ax2; position = :lb, labelsize = 9)
    # panel 3: RLI against Dx_r at n = 1e3, regular class (K = 0.5) and chaotic (K = 2 pi), Float64 wrapped
    ax3 = CM.Axis(fig[1, 3]; xscale = log10, yscale = log10, xlabel = "Dx_r", ylabel = "RLI(1e3)", title = "RLI(1e3) against Dx_r: regular class (slope 1, RL6), chaotic (flat, RL7)")
    Dsw = [1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6]
    slopes = Float64[]
    for (label, col) in ((:half, :seagreen), (:two_pi, :purple))
        c = rl_case(label, Float64, 8)
        for j in 1:8
            xs, ys = Float64[], Float64[]
            for Dx in Dsw
                pr = rl_pair(c.sys, c.x0[j], c.xi0[j], Dx)
                f = rl_run((RLI(),), PairMap(c.sys), [pr.z], [rl_frame(c.xi0[j])], 1000, 1; steps = [1000], ids = [c.ids[j]])[1][1].snaps[1]
                push!(xs, pr.Dxr); push!(ys, f.RLI)
            end
            CM.lines!(ax3, xs, max.(ys, 1e-30); color = (col, 0.6), label = j == 1 ? (label === :half ? "K = 0.5 regular class" : "K = 2 pi") : nothing)
            label === :half && push!(slopes, rl_slope(log.(xs), log.(ys)))
        end
    end
    CM.lines!(ax3, Dsw, 0.05 .* Dsw; color = :black, linestyle = :dot, label = "slope 1 guide")
    CM.axislegend(ax3; position = :rb, labelsize = 9)
    CM.save(joinpath(OUT, "rli-growth-laws.png"), fig)
    say("# figure saved: rli-growth-laws.png (panel 3: Float64 least-squares slope of log RLI against log Dx_r over Dx 1e-12..1e-6, regular class, 8 starts: min $(fmt(minimum(slopes))), max $(fmt(maximum(slopes))))")
    info("figure", "ordered start (2,0), nu = 0.5: rho = RLI / Dx_r at n = 2e4 for Dx = 1e-12, 1e-10, 1e-8; chaotic (3,0): RLI(2e4) for the same Dx", "$(join((fmt(runs[(:ordered, D)].RLI[end] / runs[(:ordered, D)].Dxr) for D in Dxs), ", ")); $(join((fmt(runs[(:chaotic, D)].RLI[end]) for D in Dxs), ", "))")
    haskey(DATA, :ct) || return
    pts = DATA[:ct]
    fig2 = CM.Figure(size = (1200, 480))
    ax4 = CM.Axis(fig2[1, 1]; xscale = log10, yscale = log10, xlabel = "RLI from the ChaosTools one-column tangent systems (BigFloat RL4)", ylabel = "RLI (Chaoslab)", title = "AC14: RLI, Chaoslab against the oracle line (Dx = 1e-12)")
    for (K, col) in ((0.5, :seagreen), (Float64(2pi), :darkorange), (5.0, :purple))
        sel = [p for p in pts if p.K == K && p.ref > 0]
        CM.scatter!(ax4, [p.ref for p in sel], [p.RLI for p in sel]; color = (col, 0.4), markersize = 4, label = "K = $(fmt(K))")
    end
    lo, hi = minimum(p.ref for p in pts if p.ref > 0), maximum(p.ref for p in pts)
    CM.lines!(ax4, [lo, hi], [lo, hi]; color = :black)
    CM.axislegend(ax4; position = :lt, labelsize = 10)
    ax5 = CM.Axis(fig2[1, 2]; yscale = log10, xscale = log10, xlabel = "n", ylabel = "|RLI_Chaoslab - RLI_oracle| / (B_R + B_oracle)", title = "AC14: ratio to the bound (bound 1)")
    for (K, col) in ((0.5, :seagreen), (Float64(2pi), :darkorange), (5.0, :purple))
        sel = [p for p in pts if p.K == K]
        CM.scatter!(ax5, [Float64(p.n) for p in sel], [max(p.ratio, 1e-6) for p in sel]; color = (col, 0.4), markersize = 4, label = "K = $(fmt(K))")
    end
    CM.hlines!(ax5, [1.0]; color = :black, linestyle = :dash)
    CM.axislegend(ax5; position = :rb, labelsize = 9)
    CM.save(joinpath(OUT, "rli-chaostools.png"), fig2)
    say("# figure saved: rli-chaostools.png")
end

say(@sprintf("# verdicts: %d PASS, %d FAIL, %d other", count(==("PASS"), VERDICTS), count(==("FAIL"), VERDICTS), count(v -> v != "PASS" && v != "FAIL", VERDICTS)))
say(@sprintf("# total wall time %.1f s", time() - T0))
close(LOG)
