# Stage 04 verification of module stretching-numbers (reducer StretchingNumbers, spectral_distance, stretching_range, check_frame; layer Indicators, workspace surface of Tangent; catalog row I10).
# Run from the repo root:  julia -t 4 --project=verify verify/stretching.jl
# Spec: workspaces/module-build/stages/01-spec/output/stretching-numbers-spec.md (criteria 1 to 16, Amendments D78, D79). Criterion 16 (CPU vs GPU) is in verify/gpu.jl, tag "stretching".
# Section "suite": the package suite at full size (SN_FULL_RUN = true: 20 starts to n = 1e4 for the reference criteria of AC4 to AC9, as the spec states them); the @info lines of the suite are
# forwarded into the log by a logger. Other sections: one printed line per measurement  AC<n> | quantity | value | bound | PASS/FAIL/INFO  with the MEASURED number, reusing the helpers of
# test/stretching/sn_helpers.jl (one home for the BigFloat reference, the edge certificate SN14, the D reference SN15, the plain-array histogram). The cases use verify seeds
# (SN_VSEED0 + 50 ac + sub), distinct from the package tests (SN_SEED + 37000 ..). The full error of a failing section is written to the log (message and backtrace), not truncated.
# SN_ONLY="AC14,R1" runs only the named sections (debugging aid) and SN_LOG names the log file; the default runs everything.
# SN_CHILD=1 is the child-process mode of the thread-invariance check of AC12(b).
using Printf, Test, Random, LinearAlgebra, Statistics, Logging
using Chaoslab, StaticArrays, DoubleFloats, StableRNGs

const SN_VSEED0 = 20261012 + 60000      # verify offset (rli used 59000)

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
    seed = SN_VSEED0 + 50 * 12 + 3
    sys = WrappedStandardMap(2pi)
    r = StretchingNumbers(; range = stretching_range(sys), bins = 1000, interval = 100)
    m, N = 4096, 1000
    for tau in (1, 4), chunk in (64, 1000)
        rng = StableRNG(seed)
        x0 = [torus_state(SVector(2pi * rand(rng), 2pi * rand(rng))) for _ in 1:m]
        res = propagate(sys, x0, N; k = 2, tau = tau, seed = seed, reducers = (r,), chunk = chunk, backend = :cpu)
        v = Float64[]
        v_flatten!(v, res.x); v_flatten!(v, res.W); v_flatten!(v, res.finals)
        out["StretchingNumbers tau=$tau chunk=$chunk"] = string(hash(reinterpret(UInt64, v)))
    end
    return out
end
if get(ENV, "SN_CHILD", "") == "1"
    for (k, h) in sort!(collect(thread_hashes())); println("HASH ", k, " = ", h); end
    println("THREADS ", Threads.nthreads())
    exit(0)
end

const SN_FULL_RUN = true
const ROOT = dirname(@__DIR__)
const TESTDIR = joinpath(ROOT, "test", "stretching")
include(joinpath(TESTDIR, "sn_helpers.jl"))

vseed(ac, sub = 0) = SN_VSEED0 + 50 * ac + sub

import ChaosTools as CT
import CairoMakie as CM

const OUT = joinpath(ROOT, "workspaces", "module-build", "stages", "04-verify", "output")
mkpath(OUT)
const ONLY = filter(!isempty, split(get(ENV, "SN_ONLY", ""), ","))
const LOG = open(get(ENV, "SN_LOG", joinpath(OUT, "stretching-numbers-verify-log.txt")), "w")
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
rat(x, B) = x <= 0 ? 0.0 : (B <= 0 ? Inf : x / B)
rmedian(v) = sort(v)[max((length(v) + 1) ÷ 2, 1)]
quants(v, qs = (0.05, 0.25, 0.5, 0.75, 0.95)) = join((fmt(quantile(v, q)) for q in qs), " / ")

function section(f, name)
    isempty(ONLY) || name in ONLY || return
    t0 = time()
    try
        f()
    catch e
        say(string(name, " | ERROR | ", sprint(showerror, e, catch_backtrace()), " | - | FAIL"))
        push!(VERDICTS, "FAIL")
    end
    say(@sprintf("# section %s took %.1f s", name, time() - t0))
end

const T0 = time()
const DATA = Dict{Any, Any}()
const EPS64 = sn_eps(Float64)
say("# stretching-numbers verification, julia $(VERSION), threads $(Threads.nthreads()), base seed $SN_SEED (verify offset 60000), time $(round(time())) (unix)")
say("# timings in this log were taken on a machine that was not otherwise idle")

# ---------------------------------------------------------------------------------------------------
# The package suite at full size; its @info lines are forwarded into this log (the gated ratio keys are listed in GATES)
const GATES = ("AC6 realized" => ((:worst_D_ratio, 1.0),), "AC6b realized" => ((:worst_D_ratio, 1.0),), "AC5(d) realized" => ((:worst_sd, 1.0),), "AC9 realized" => ((:worst_a, 1.0), (:worst_b, 1.0)))
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
            @testset verbose = true "stretching-numbers verify" begin
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
    emit("SUITE", "package suite at full size (SN_FULL_RUN), assertions passed (fail + error = $nbad); wall time $(round(t; digits = 1)) s", npass, "0 fail, 0 error", pf(ok && nbad == 0 && npass > 0))
end

# ---------------------------------------------------------------------------------------------------
# The independent tangent line of ChaosTools 3.5.4 (Oracle a): a one-column TangentDynamicalSystem with the hand-coded Jacobian of the lifted standard map. The oracle's rule is the lifted
# one (as the sali, rli and megno lines), so the module side of criteria 14 and 15 is the lifted StandardMap.
ct_rule(u, p, n) = (pp = u[2] + p[1] * sin(u[1]); SVector(u[1] + pp, pp))
ct_jac(u, p, n) = (c = p[1] * cos(u[1]); SMatrix{2, 2}(1 + c, c, 1, 1))
ct_tds(K) = CT.TangentDynamicalSystem(CT.DeterministicIteratedMap(ct_rule, SVector(0.1, 0.1), [K]); J = ct_jac, k = 1)   # one per call (threads)
# One column w0 from z0: tau oracle steps, then the deviation is read, its log norm taken in BigFloat(256), and the deviation renormalised (Float64). Returns the log increments a_k (BigFloat, one
# per handoff), the number of steps at which the oracle orbit differs bitwise from forward(StandardMap), the unit column and the orbit point at the start of each window.
function ct_run(K::Float64, z0, w0, n::Int, tau::Int)
    tds = ct_tds(K)
    CT.reinit!(tds, SVector{2, Float64}(z0); Q0 = SMatrix{2, 1, Float64}(w0[1], w0[2]))
    sys = StandardMap(K)
    z = SVector{2, Float64}(z0)
    w = SVector{2, Float64}(w0[1], w0[2])
    nk = n ÷ tau
    a = Vector{BigFloat}(undef, nk)
    ws = Vector{SVector{2, Float64}}(undef, nk)
    pts = Vector{SVector{2, Float64}}(undef, nk)
    mism = 0
    setprecision(BigFloat, SN_BITS) do
        for k in 1:nk
            ws[k], pts[k] = w, z
            for _ in 1:tau
                CT.step!(tds)
                z = forward(sys, z)
                mism += SVector{2, Float64}(CT.current_state(tds)) == z ? 0 : 1
            end
            w = SVector{2, Float64}(CT.current_deviations(tds)[:, 1])
            a[k] = log(BigFloat(w[1])^2 + BigFloat(w[2])^2) / 2
            w = w / norm(w)
            CT.set_deviations!(tds, SMatrix{2, 1, Float64}(w[1], w[2]))
        end
    end
    return a, mism, ws, pts
end

# ---------------------------------------------------------------------------------------------------
# AC14: the z series of the oracle, histogram and D by the plain formula in BigFloat, against the module's counts and D (envelope of SN14 and bound of SN15, each doubled: ours plus the oracle's own,
# which has the same width because the orbits are bitwise equal and the tangents the same).
section("AC14") do
    DATA[:ct] = NamedTuple[]
    B = 1000
    panels = ((:half, 20, [100, 1000, 10_000], 1, "K = 0.5 regular class"), (:two_pi, 1000, [100, 1000], 1, "K = 2 pi"), (:five, 1000, [100, 1000], 1, "K = 5"),
              (:half, 20, [100, 1000, 10_000], 4, "K = 0.5 regular class, tau = 4 (extra)"))
    for (label, m, ns, tau, descr) in panels
        c = sn_case(label, Float64, m; wrapped = false)
        K = Float64(c.K)
        N = last(ns)
        ks = [n ÷ tau for n in ns]
        lo, hi = stretching_range(c.sys)
        h = setprecision(BigFloat, SN_BITS) do
            (sn_big(hi) - sn_big(lo)) / B
        end
        res = Vector{Any}(undef, m)
        Threads.@threads :dynamic for j in 1:m
            ref = sn_ref_orbit(c.sys, c.K, c.x0[j], c.W0[j], N)
            snaps = sn_run((sn_red(c.sys, B, tau),), c.sys, [c.x0[j]], [c.W0[j]], N, tau; steps = ns, ids = [c.ids[j]])[1][1].snaps
            o = [ct_run(K, c.x0[j], c.W0[j][:, i], N, tau) for i in 1:2]
            mism = o[1][2] + o[2][2]
            zref = [sn_z_from_s(ref.s1, tau), sn_z_from_s(ref.s2, tau)]
            zorc = [o[i][1] ./ tau for i in 1:2]
            eps = [sn_halfwidths(ref, i, EPS64, tau, N ÷ tau, lo, hi) for i in 1:2]
            resolved = all(maximum(eps[i]) < Float64(h) / 2 for i in 1:2)
            env = [sn_envelope(zref[i], eps[i], lo, hi, B, ks) for i in 1:2]
            orc = [sn_prefix_counts(zorc[i], lo, hi, B, ks) for i in 1:2]
            rows = NamedTuple[]
            for (f, n, k) in zip(snaps, ns, ks)
                f.n == n || error("AC14: snapshot step $(f.n) vs $n")
                rN = 0.0; ndiff = 0
                for i in 1:2
                    cnt = sn_cnt(f, i)
                    width = 2 .* (env[i][k].Np .- env[i][k].Nm)
                    for cc in eachindex(cnt)
                        d = abs(cnt[cc] - orc[i][k][cc])
                        ndiff += d != 0
                        rN = max(rN, d == 0 ? 0.0 : (width[cc] == 0 ? Inf : d / width[cc]))
                    end
                end
                Dorc = sn_D_ref(orc[1][k], orc[2][k], lo, hi, B)
                Dref = sn_D_ref(env[1][k].N, env[2][k].N, lo, hi, B)
                U = env[1][k].U + env[2][k].U
                bound = 2 * sqrt(BigFloat(2)) * U / (k * sqrt(h)) + 8 * EPS64 * Dref
                eD = abs(sn_big(f.D) - Dorc)
                push!(rows, (n = n, rN = rN, ndiff = ndiff, rD = rat(Float64(eD), Float64(bound)), Dours = Float64(f.D), Dorc = Float64(Dorc), U = U))
            end
            res[j] = (mism = mism, resolved = resolved, rows = rows)
        end
        mism_total = sum(r.mism for r in res)
        nbit = count(r -> r.mism == 0, res)
        emit("AC14", "$descr: oracle orbits (one-column systems) vs forward(StandardMap) over $m starts x $N steps x 2 columns: steps with a bitwise difference; starts with both columns fully bitwise", "$mism_total; $nbit of $m", "0 (else excluded, reported)", pf(mism_total == 0))
        tot = excl = unres = 0
        wN = wD = 0.0; nd = 0; Ucount = 0
        for r in res, row in r.rows
            tot += 1
            if r.mism != 0 || !r.resolved
                excl += 1
                unres += (r.mism == 0 && !r.resolved)
                continue
            end
            wN = max(wN, row.rN); wD = max(wD, row.rD); nd += row.ndiff != 0; Ucount += row.U
            tau == 1 && push!(DATA[:ct], (K = K, n = row.n, Dours = row.Dours, Dorc = row.Dorc, rD = row.rD, rN = row.rN))
        end
        ratio_line("AC14", "$descr, $m starts, n in $ns: max over cells and columns |N_ours - N_oracle| / (2 (N+ - N-)) ($(tot - excl) comparisons)", wN)
        ratio_line("AC14", "$descr: max |D_ours - D_oracle| / (SN15 bound of ours + the oracle's own, 2 sqrt 2 U / (k sqrt h) + 8 eps D_ref)", wD)
        info("AC14", "$descr: comparisons with at least one count differing between module and oracle; sum of U (samples within a half-width of two cells) over the comparisons", "$nd of $(tot - excl); $Ucount")
        emit("AC14", "$descr: comparisons excluded (oracle orbit not bitwise equal, or unresolved; unresolved alone: $unres)", excl / tot, "<= 0.10", pf(excl / tot <= 0.10))
    end
end

# ---------------------------------------------------------------------------------------------------
# AC15: ChaosTools local_growth_rates (finite perturbation of the full nonlinear map, one window, no rescaling) with a supplied perturbation e = 1e-6 along the unit column at each window start,
# against the z_k of the module (the difference of its compensated totals), at every handoff up to 100, tau = 1, both columns.
section("AC15") do
    DATA[:lg] = NamedTuple[]
    e = 1e-6
    Nw = 100
    for (label, descr) in ((:half, "K = 0.5"), (:five, "K = 5"), (:ten, "K = 10"))
        m = 20
        c = sn_case(label, Float64, m; wrapped = false)
        K = Float64(c.K)
        sig = Float64(exp(sn_Lsigma(c.K)))
        rho = K * e / sqrt(2.0)
        term1 = -log1p(-sig * rho)
        res = Vector{Any}(undef, m)
        Threads.@threads :dynamic for j in 1:m
            lights = sn_run((sn_red(c.sys, 1000, 1),), c.sys, [c.x0[j]], [c.W0[j]], Nw, 1; light = true, ids = [c.ids[j]])[1][1].lights
            ds = CT.DeterministicIteratedMap(ct_rule, c.x0[j], [K])
            out = NamedTuple[]
            mism = 0
            for col in 1:2
                a, mi, ws, pts = ct_run(K, c.x0[j], c.W0[j][:, col], Nw, 1)
                mism += mi
                zmod = sn_z_from_lights(lights, col, 1)
                Mst = max(1.0, abs(K))
                for k in 1:Nw
                    Mst = max(Mst, maximum(abs, pts[k]), maximum(abs, forward(c.sys, pts[k])))
                    lg = CT.local_growth_rates(ds, [pts[k]]; S = 1, Δt = 1, perturbation = (ds_, u_, jj_) -> e * ws[k])[1, 1]
                    err = abs(Float64(zmod[k]) - lg)
                    err_orc = abs(Float64(a[k]) - lg)
                    b2 = 2 * (K + 3) * Mst * EPS64 * sig / e
                    push!(out, (k = k, col = col, err = err, err_orc = err_orc, bound = term1 + b2, term1 = term1, b2 = b2))
                end
            end
            res[j] = (mism = mism, rows = out)
        end
        mism_total = sum(r.mism for r in res)
        emit("AC15", "$descr: oracle orbits vs forward(StandardMap) over $m starts x $Nw steps x 2 columns: steps with a bitwise difference", mism_total, "0", pf(mism_total == 0))
        rows = [row for r in res for row in r.rows]
        worst = maximum(rat(row.err, row.bound) for row in rows)
        worst1 = maximum(rat(row.err, row.term1) for row in rows)
        worsto = maximum(rat(row.err_orc, row.bound) for row in rows)
        ratio_line("AC15", "$descr, $m starts, handoffs 1 to $Nw, both columns ($(length(rows)) comparisons), e = 1e-6: max |z_k - local growth rate| / (-ln(1 - sigma rho) + 2 (K + 3) M_state eps sigma / e); first term $(fmt(term1)), sigma $(fmt(sig))", worst)
        info("AC15", "$descr: same ratio against the first term alone (SN17 without the rounding term); against the oracle's own log norm a_k instead of the module's z_k; median |z_k - lambda_local|", "$(fmt(worst1)); $(fmt(worsto)); $(fmt(rmedian([row.err for row in rows])))")
        for row in rows
            push!(DATA[:lg], (K = K, err = max(row.err, 1e-18), bound = row.bound))
        end
    end
end

# ---------------------------------------------------------------------------------------------------
# AC12(b): 1 vs N threads in child processes (output hashes), with the in-process chunk-size comparison
section("AC12") do
    hashes = Dict{Int, Dict{String, String}}(Threads.nthreads() => thread_hashes())
    for nt in (1, 2, 4)
        haskey(hashes, nt) && continue
        cmd = addenv(`$(Base.julia_cmd()) -t $nt --project=$(dirname(Base.active_project())) $(@__FILE__)`, "SN_CHILD" => "1")
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
    for tau in (1, 4)
        ok_all &= length(unique([hashes[Threads.nthreads()]["StretchingNumbers tau=$tau chunk=$chunk"] for chunk in (64, 1000)])) == 1
    end
    emit("AC12", "(b) x, W and every final (counts, D, sums) bitwise identical for thread counts $(sort(collect(keys(hashes)))) (hashes from separate processes) and chunk sizes 64 and 1000: $(length(keys_)) cases (tau in {1,4} x chunk), 4096 orbits, K = 2 pi, wrapped, n = 1e3; example hash $(first(values(hashes[1])))", ok_all, "true", pf(ok_all))
end

# ---------------------------------------------------------------------------------------------------
# Shared by R1 and R2: the K = 0.5 regular class (L1 S7 screening: Benettin X(1e5) <= 5e-4 on the lifted map, tau = 10) of m starts out of 3 m candidates, and uniform starts of the other classes.
function regular_class(m)
    haskey(DATA, (:reg, m)) && return DATA[(:reg, m)]
    rng = StableRNG(vseed(40, 1))
    cand = [SVector(2pi * rand(rng), 2pi * rand(rng)) for _ in 1:3m]
    X = [f[1].lambda[1] for f in propagate(StandardMap(0.5), cand, 10^5; k = 1, tau = 10, seed = vseed(40, 2), reducers = (Benettin(),), backend = :cpu).finals]
    reg = findall(<=(5e-4), X)
    length(reg) >= m || error("only $(length(reg)) of $(3m) candidates are regular")
    return DATA[(:reg, m)] = (starts = cand[reg[1:m]], ncand = 3m, nreg = length(reg))
end
function class_starts(label, m)
    label === :half && return regular_class(m).starts
    rng = StableRNG(vseed(41, label === :two_pi ? 1 : label === :five ? 2 : 3))
    return [SVector(2pi * rand(rng), 2pi * rand(rng)) for _ in 1:m]
end
const R_CLASSES = ((:half, 0.5), (:two_pi, 2pi), (:five, 5.0), (:ten, 10.0))

lslope(x, y) = (xm = mean(x); ym = mean(y); sum((x .- xm) .* (y .- ym)) / sum((x .- xm) .^ 2))

# ---------------------------------------------------------------------------------------------------
# R1: class separation of D on the 2D panel (F5). Wrapped, tau = 1, B = 1000 over the derived range, interval 100, seeded frame, 1000 starts per class.
section("R1") do
    m = 1000
    Zs = (1_000, 2_000, 5_000, 10_000, 20_000, 50_000, 100_000)
    rep = (1_000, 10_000, 100_000)
    D = Dict{Symbol, Dict{Int, Vector{Float64}}}()
    seedR = vseed(42, 1)
    reg = regular_class(m)
    info("R1", "K = 0.5 regular class: candidates screened, regular (Benettin X(1e5) <= 5e-4), taken", "$(reg.ncand), $(reg.nreg), $m")
    for (label, K) in R_CLASSES
        sys = WrappedStandardMap(K)
        x0 = [torus_state(s) for s in class_starts(label, m)]
        r = StretchingNumbers(; range = stretching_range(sys), bins = 1000, interval = 100)
        D[label] = Dict{Int, Vector{Float64}}()
        for N in Zs
            res = propagate(sys, x0, N; k = 2, tau = 1, seed = seedR, reducers = (r,), backend = :cpu)
            D[label][N] = [Float64(f[1].D) for f in res.finals]
            N == 100_000 && info("R1", "K = $(fmt(K)), $m starts, n = 1e5: undersampled flag (k < 100 B = 1e5) is set for; any overflow count nonzero", "$(count(f -> f[1].undersampled, res.finals)); $(count(f -> any(f[1].below .> 0) || any(f[1].above .> 0), res.finals))")
        end
    end
    for (label, K) in R_CLASSES, n in rep
        d = D[label][n]
        info("R1", "K = $(fmt(K)) ($label), n = Z = $n, $m starts: quantiles 5 / 25 / 50 / 75 / 95 % of D; of D times Z", "$(quants(d)); $(quants(d .* n))")
    end
    slopes = Dict{Symbol, Vector{Float64}}()
    for (label, K) in R_CLASSES
        slopes[label] = [begin
            ds = [D[label][N][j] for N in Zs]
            all(>(0), ds) ? lslope(log.(Float64.(collect(Zs))), log.(ds)) : NaN
        end for j in 1:m]
        s = filter(!isnan, slopes[label])
        info("R1", "K = $(fmt(K)) ($label): log-log slope of D against Z over Z = 1e3 .. 1e5 (7 points) per orbit: quantiles 5 / 25 / 50 / 75 / 95 %; orbits with D = 0 at some Z (no slope)", "$(quants(s)); $(m - length(s))")
        info("R1", "K = $(fmt(K)) ($label): fraction of orbits whose slope is not within 0.25 of -1 (decay not close to 1/Z); not within 0.10", "$(fmt(count(x -> abs(x + 1) > 0.25, s) / length(s))); $(fmt(count(x -> abs(x + 1) > 0.10, s) / length(s)))")
    end
    for n in rep
        q = median(D[:half][n])
        info("R1", "n = $n: regular-class median of D $(fmt(q)); fraction of orbits with D above it: K = 2 pi; 5; 10; and regular class 5 % / 95 % quantile of D", "$(fmt(mean(D[:two_pi][n] .> q))); $(fmt(mean(D[:five][n] .> q))); $(fmt(mean(D[:ten][n] .> q))); $(fmt(quantile(D[:half][n], 0.05))) / $(fmt(quantile(D[:half][n], 0.95)))")
        info("R1", "n = $n: ratio of the regular-class median of D times Z to the chaotic-class medians (K = 2 pi; 5; 10), and the ratio of the medians of D", "$(fmt(median(D[:half][n] .* n) / median(D[:two_pi][n] .* n))); $(fmt(median(D[:half][n] .* n) / median(D[:five][n] .* n))); $(fmt(median(D[:half][n] .* n) / median(D[:ten][n] .* n))); $(fmt(q / median(D[:two_pi][n]))); $(fmt(q / median(D[:five][n]))); $(fmt(q / median(D[:ten][n])))")
    end
    # cross-tabulation with the labels of L1 (n X above ln n), FLI_B (above ln n), MEGNO (Ybar_1 above 3) and RLI (rho = RLI / Dx_r above 1 at Dx = 1e-12), the thresholds of Maffione2011a as used in
    # rli R4 (not pre-registered). Stretching side: D at most / above the regular-class median at the same n. Cells: (label chaotic, D low) / (label chaotic, D high) / (label regular, D low) / (label regular, D high).
    labels = Dict{Tuple{Symbol, Int}, NamedTuple}()
    for (label, K) in R_CLASSES
        base = WrappedStandardMap(K)
        P = PairMap(base)
        starts = class_starts(label, m)
        xw = [torus_state(s) for s in starts]
        xis = [SVector{2, Float64}(initial_deviations(seedR, j, Val(2), Val(1), Float64)[:, 1]) for j in 1:m]
        ps = [pair_start(base, starts[j], xis[j], 1e-12; orbit = j) for j in 1:m]
        z0 = [p.z for p in ps]
        W0 = [SMatrix{4, 2, Float64}(xi[1], xi[2], 0.0, 0.0, 0.0, 0.0, xi[1], xi[2]) for xi in xis]
        for n in rep
            r = propagate(P, z0, n; k = 2, tau = 1, seed = seedR, reducers = (RLI(),), W0 = W0, backend = :cpu)
            l1 = propagate(base, xw, n; k = 1, tau = 1, seed = seedR, reducers = (Benettin(),), backend = :cpu)
            fm = propagate(base, xw, n; k = 2, tau = 1, seed = seedR, reducers = (FLI(replace = false), MEGNO()), backend = :cpu)
            labels[(label, n)] = (RLI = [r.finals[j][1].RLI / ps[j].Dx_r > 1.0 for j in 1:m], L1 = [f[1].lambda[1] * n > log(n) for f in l1.finals],
                                  FLI = [f[1].fli > log(n) for f in fm.finals], MEGNO = [f[2].Ybar[1] > 3 for f in fm.finals])
        end
    end
    for n in rep
        q = median(D[:half][n])
        for nm in (:L1, :FLI, :MEGNO, :RLI)
            lab = vcat((getproperty(labels[(label, n)], nm) for (label, _) in R_CLASSES)...)
            low = vcat((D[label][n] .<= q for (label, _) in R_CLASSES)...)
            info("R1", "n = $n, all 4000 orbits (4 classes), label $nm against D <= regular-class median ($(fmt(q))): (chaotic, low) / (chaotic, high) / (regular, low) / (regular, high); fraction labelled chaotic per class (K = 0.5; 2 pi; 5; 10)",
                 "$(count(lab .& low)) / $(count(lab .& .!low)) / $(count(.!lab .& low)) / $(count(.!lab .& .!low)); $(join((fmt(mean(getproperty(labels[(label, n)], nm))) for (label, _) in R_CLASSES), "; "))")
        end
    end
    DATA[:R1] = D
end

# ---------------------------------------------------------------------------------------------------
# R2: sensitivity of the regular to chaotic ratio of D at n = 1e4 to B, the interval, the initial pair and tau (reported; Maffione2011a p.20: almost two orders of magnitude in the 4D map).
section("R2") do
    n = 10_000
    nreg, nch = 100, 100
    sregs = class_starts(:half, nreg)
    chaotic = ((:two_pi, 2pi), (:ten, 10.0))
    schs = Dict(l => class_starts(l, nch) for (l, _) in chaotic)
    seedR = vseed(43, 1)
    function Dat(K, starts, frames; B = 1000, tau = 1, interval = nothing)
        sys = WrappedStandardMap(K)
        r = StretchingNumbers(; range = stretching_range(sys), bins = B, interval = interval === nothing ? sn_interval(tau) : interval)
        x0 = [torus_state(s) for s in starts]
        res = propagate(sys, x0, n; k = 2, tau = tau, seed = seedR, reducers = (r,), backend = :cpu, W0 = frames)
        return [Float64(f[1].D) for f in res.finals]
    end
    function line(name, Wfun; B = 1000, tau = 1, interval = nothing, rep = 1)
        rep_(v) = repeat(v, inner = rep)
        Wm(m) = Wfun === nothing ? nothing : Wfun(m * rep)
        dr = Dat(0.5, rep_(sregs), Wm(nreg); B = B, tau = tau, interval = interval)
        out = Dict{Symbol, Float64}()
        for (l, K) in chaotic
            d = Dat(K, rep_(schs[l]), Wm(nch); B = B, tau = tau, interval = interval)
            out[l] = median(d)
        end
        mr = median(dr)
        info("R2", "$name: median D at n = 1e4, regular class; chaotic K = 2 pi; K = 10; ratio regular / chaotic (K = 2 pi; K = 10)", "$(fmt(mr)); $(fmt(out[:two_pi])); $(fmt(out[:ten])); $(fmt(out[:two_pi] > 0 ? mr / out[:two_pi] : NaN)); $(fmt(out[:ten] > 0 ? mr / out[:ten] : NaN))")
        return dr
    end
    base = line("baseline: seeded frame, B = 1000, interval 100, tau = 1", nothing)
    for B in (100, 10_000)
        line("B = $B", nothing; B = B)
    end
    # the interval does not enter D at the last handoff (Amendment D79): bitwise equal D for 10, 100, 1000
    ds = [Dat(0.5, sregs, nothing; interval = iv) for iv in (10, 100, 1000)]
    emit("R2", "interval in {10, 100, 1000}: D at n = 1e4 of the $nreg regular starts bitwise equal (the interval fixes only the evaluation handoffs of a probe or recorder)", all(d -> d == ds[2], ds), "true (expected)", pf(all(d -> d == ds[2], ds)))
    for tau in (2, 4, 8, 16)
        line("tau = $tau (interval $(sn_interval(tau)))", nothing; tau = tau)
    end
    unit(v) = v / norm(v)
    fixed(a, b) = m -> [SMatrix{2, 2, Float64}(a[1], a[2], b[1], b[2]) for _ in 1:m]
    line("canonical basis (1,0), (0,1)", fixed(SVector(1.0, 0.0), SVector(0.0, 1.0)))
    line("unit (1,1) with unit (2,1)", fixed(unit(SVector(1.0, 1.0)), unit(SVector(2.0, 1.0))))
    line("unit (1,1) with unit (1,2)", fixed(unit(SVector(1.0, 1.0)), unit(SVector(1.0, 2.0))))
    line("parallel pair (1,0), (1,0) (control, D = 0 expected)", fixed(SVector(1.0, 0.0), SVector(1.0, 0.0)))
    rng = StableRNG(seedR + 7)
    randframes(m) = [(t = 2pi * rand(rng); SMatrix{2, 2, Float64}(cos(t), sin(t), -sin(t), cos(t))) for _ in 1:m]
    dr = line("random orthonormal orientations, 20 per start (medians over all orbits)", randframes; rep = 20)
    per = [dr[(i - 1) * 20 + 1:i * 20] for i in 1:nreg]
    info("R2", "random orientations, regular class: per-start max / min of D over the 20 orientations: quantiles 5 / 25 / 50 / 75 / 95 % (start-wise spread of the same orbit with different vectors)", quants([maximum(p) / max(minimum(p), 1e-300) for p in per]))
    info("R2", "seeded baseline, regular class: median D times Z at n = 1e4 (for the orders of magnitude between classes read the baseline line above)", fmt(median(base) * n))
end

# ---------------------------------------------------------------------------------------------------
# Figures
section("figure") do
    # spectra (SSN) of a regular and a chaotic orbit of the standard map, from the module
    n = 100_000
    reg = regular_class(100).starts[1]
    chs = SVector(2pi * 0.3141592653589793, 2pi * 0.7182818284590452)
    cases = ((WrappedStandardMap(0.5), reg, "regular orbit, K = 0.5, start ($(fmt(reg[1])), $(fmt(reg[2])))", :seagreen),
             (WrappedStandardMap(2pi), chs, "chaotic orbit, K = 2 pi, start ($(fmt(chs[1])), $(fmt(chs[2])))", :darkorange))
    fig = CM.Figure(size = (1800, 560))
    axs = [CM.Axis(fig[1, i]; xlabel = "stretching number z", ylabel = "S(z) = N_c / (k h)", yscale = log10, title = cases[i][3]) for i in 1:2]
    axD = CM.Axis(fig[1, 3]; xscale = log10, yscale = log10, xlabel = "Z = number of stretching numbers (n, tau = 1)", ylabel = "spectral distance D", title = "D(Z), B = 1000, every 100 steps (dotted: 1/Z)")
    for (i, (sys, s, _, col)) in enumerate(cases)
        x0 = [torus_state(s)]
        W0 = [sn_frame(vseed(44, 1), 1, Float64)]
        lo, hi = stretching_range(sys)
        B = 100
        snap = sn_run((sn_red(sys, B, 1),), sys, x0, W0, n, 1; steps = [n])[1][1].snaps[end]
        h = (hi - lo) / B
        centres = [lo + (c - 0.5) * h for c in 1:B]
        for (j, ls) in ((1, :solid), (2, :dash))
            cnt = sn_cnt(snap, j)
            dens = [cnt[c + 1] / (n * h) for c in 1:B]
            CM.stairs!(axs[i], centres, max.(dens, 1e-6); color = col, linestyle = ls, step = :center, label = "column $j")
        end
        CM.ylims!(axs[i], 1e-4, nothing)
        CM.axislegend(axs[i]; position = :ct, labelsize = 10)
        evs = collect(100:100:n)
        sn = sn_run((sn_red(sys, 1000, 1),), sys, x0, W0, n, 1; steps = evs)[1][1].snaps
        CM.lines!(axD, Float64.(evs), max.(Float64[f.D for f in sn], 1e-12); color = col, label = i == 1 ? "regular K = 0.5" : "chaotic K = 2 pi")
    end
    CM.lines!(axD, [1e2, 1e5], [3.0, 3e-3]; color = :black, linestyle = :dot, label = "1/Z guide")
    CM.axislegend(axD; position = :lb, labelsize = 10)
    CM.save(joinpath(OUT, "stretching-numbers-spectra.png"), fig)
    say("# figure saved: stretching-numbers-spectra.png")
    haskey(DATA, :ct) || return
    pts = DATA[:ct]
    fig2 = CM.Figure(size = (1800, 560))
    ax1 = CM.Axis(fig2[1, 1]; xscale = log10, yscale = log10, xlabel = "D from the ChaosTools one-column tangent systems (plain formula, BigFloat)", ylabel = "D (Chaoslab)", title = "AC14: spectral distance, module against the oracle line")
    cols = ((0.5, :seagreen), (Float64(2pi), :darkorange), (5.0, :purple))
    for (K, col) in cols
        sel = [p for p in pts if p.K == K && p.Dorc > 0]
        CM.scatter!(ax1, [p.Dorc for p in sel], [p.Dours for p in sel]; color = (col, 0.4), markersize = 4, label = "K = $(fmt(K))")
    end
    lo_, hi_ = minimum(p.Dorc for p in pts if p.Dorc > 0), maximum(p.Dorc for p in pts)
    CM.lines!(ax1, [lo_, hi_], [lo_, hi_]; color = :black)
    CM.axislegend(ax1; position = :lt, labelsize = 10)
    ax2 = CM.Axis(fig2[1, 2]; xscale = log10, yscale = log10, xlabel = "n", ylabel = "|D_ours - D_oracle| / bound", title = "AC14: D difference over its bound (bound 1)")
    for (K, col) in cols
        sel = [p for p in pts if p.K == K]
        CM.scatter!(ax2, [Float64(p.n) for p in sel], [max(p.rD, 1e-6) for p in sel]; color = (col, 0.4), markersize = 4, label = "K = $(fmt(K))")
    end
    CM.hlines!(ax2, [1.0]; color = :black, linestyle = :dash)
    CM.axislegend(ax2; position = :rb, labelsize = 9)
    haskey(DATA, :lg) || (CM.save(joinpath(OUT, "stretching-numbers-chaostools.png"), fig2); return)
    ax3 = CM.Axis(fig2[1, 3]; xscale = log10, yscale = log10, xlabel = "bound of SN17 (first term + rounding term)", ylabel = "|z_k - local growth rate|", title = "AC15: local growth rate, e = 1e-6 (dashed: y = x)")
    for (K, col) in ((0.5, :seagreen), (5.0, :purple), (10.0, :firebrick))
        sel = [p for p in DATA[:lg] if p.K == K]
        CM.scatter!(ax3, [p.bound for p in sel], [p.err for p in sel]; color = (col, 0.3), markersize = 3, label = "K = $(fmt(K))")
    end
    bl = extrema(p.bound for p in DATA[:lg])
    CM.lines!(ax3, [bl[1], bl[2]], [bl[1], bl[2]]; color = :black, linestyle = :dash)
    CM.axislegend(ax3; position = :lt, labelsize = 10)
    CM.save(joinpath(OUT, "stretching-numbers-chaostools.png"), fig2)
    say("# figure saved: stretching-numbers-chaostools.png")
end

say(@sprintf("# verdicts: %d PASS, %d FAIL, %d other", count(==("PASS"), VERDICTS), count(==("FAIL"), VERDICTS), count(v -> v != "PASS" && v != "FAIL", VERDICTS)))
say(@sprintf("# total wall time %.1f s", time() - T0))
close(LOG)
