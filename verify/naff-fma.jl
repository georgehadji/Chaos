# Stage 04 verification of module naff-fma (E1: window_frequency, frequency_map; E2: WeightedPhaseAdvance, weighted_phase_frequency; layer Indicators; catalog rows I12, I12b).
# Run from the repo root:  julia -t 16 --project=verify verify/naff-fma.jl
# Spec: workspaces/module-build/stages/01-spec/output/naff-fma-spec.md (criteria 1 to 16, Amendments D80 and D81). Criterion 15 (CPU vs GPU) is in verify/gpu.jl, tag "naff".
# Section "suite": the package suite at the spec sizes (NF_FULL_RUN = true: 20 starts, N up to 4096 and 10001, N = 16000 in AC8(c)); its @info lines are forwarded into the log. Section "AC14": the
# independent BigFloat orbit and the loose external line (ChaosTools estimate_period); AC12: 1 versus N threads in child processes; R1 to R6 and the Eq. (79) part of R2: reported numbers, no gate
# (the rli D77 form). Output lines:  AC<n> | quantity | value | bound | PASS/FAIL/INFO  with the MEASURED number. Cases use verify seeds (NV_SEED0 + 50 ac + sub), distinct from the package tests.
# The full error of a failing section is written to the log (message and backtrace), not truncated. NF_ONLY="AC14,R1" runs only the named sections, NF_LOG names the log file.
# NF_CHILD=1 is the child-process mode of the thread-invariance check of AC12.
using Printf, Test, Random, LinearAlgebra, Statistics, Logging
using Chaoslab, StaticArrays, DoubleFloats, StableRNGs

const NV_SEED0 = 20261012 + 61000      # verify offset (stretching used 60000)

# One orbit's states at the handoffs 1..n (tau = 1): the samples that the stream reducers see.
function orbit_states(sys, z0, n::Integer)
    out = Vector{typeof(z0)}(undef, n)
    z = z0
    for i in 1:n
        z = forward(sys, z)
        out[i] = z
    end
    return out
end
obs_cis(s) = complex(cos(s[1]), sin(s[1]))       # exp(i q) from the wrapped coordinate (spec P4)
obs_planar(s) = complex(s[1], s[2])              # q + i p, wrapped pair

# ---------------------------------------------------------------------------------------------------
# Thread-count invariance (AC12): the same outputs hashed in this process and in child processes with other thread counts.
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
e1_item(fm) = (fm.w1[1].nu, fm.w1[1].amplitude, fm.w1[1].mirror_amplitude, fm.w1[1].second_ratio, Float64(fm.w1[1].bin), fm.w1[1].tie, fm.w1[1].near_zero, fm.w1[1].mirror,
               fm.w2[1].nu, fm.w2[1].amplitude, fm.w2[1].mirror_amplitude, fm.w2[1].second_ratio, Float64(fm.w2[1].bin), fm.w2[1].tie, fm.w2[1].near_zero, fm.w2[1].mirror,
               fm.delta, fm.D, fm.below_floor)
function thread_hashes()
    out = Dict{String, String}()
    seed = NV_SEED0 + 50 * 12 + 3
    N, m = 1000, 4096
    for (K, tag) in ((0.5, "K=0.5"), (2pi, "K=2pi"))
        sys = WrappedStandardMap(K)
        rng = StableRNG(seed + (K == 0.5 ? 0 : 1))
        x0 = [torus_state(SVector(2pi * rand(rng), 2pi * rand(rng))) for _ in 1:m]
        r = WeightedPhaseAdvance(; window = N, observable = LiftedAngleObservable(1, 3))
        for chunk in (64, 1000)
            res = propagate(sys, x0, 2N + 1; k = 1, tau = 1, seed = seed, reducers = (r, Benettin()), chunk = chunk, backend = :cpu)
            v = Float64[]
            v_flatten!(v, res.x); foreach(f -> v_flatten!(v, f[1]), res.finals)
            out["E2 $tag chunk=$chunk"] = string(hash(reinterpret(UInt64, v)))
        end
        for chunk in (64, 1000)
            fms = Vector{Any}(undef, m)
            for c0 in 1:chunk:m
                idx = c0:min(c0 + chunk - 1, m)
                Threads.@threads :dynamic for j in idx
                    st = orbit_states(sys, x0[j], 2N)
                    fms[j] = e1_item(frequency_map(ComplexF64[obs_cis(s) for s in st], N; mean_removal = true))
                end
            end
            v = Float64[]
            v_flatten!(v, fms)
            out["E1 $tag chunk=$chunk"] = string(hash(reinterpret(UInt64, v)))
        end
    end
    return out
end
if get(ENV, "NF_CHILD", "") == "1"
    for (k, h) in sort!(collect(thread_hashes())); println("HASH ", k, " = ", h); end
    println("THREADS ", Threads.nthreads())
    exit(0)
end

const NF_FULL_RUN = true
const ROOT = dirname(@__DIR__)
const TESTDIR = joinpath(ROOT, "test", "naff-fma")
include(joinpath(TESTDIR, "nf_helpers.jl"))

vseed(ac, sub = 0) = NV_SEED0 + 50 * ac + sub

import ChaosTools as CT
import CairoMakie as CM

const OUT = joinpath(ROOT, "workspaces", "module-build", "stages", "04-verify", "output")
mkpath(OUT)
const ONLY = filter(!isempty, split(get(ENV, "NF_ONLY", ""), ","))
const LOG = open(get(ENV, "NF_LOG", joinpath(OUT, "naff-fma-verify-log.txt")), "w")
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
quants(v, qs = (0.05, 0.25, 0.5, 0.75, 0.95)) = isempty(v) ? "-" : join((fmt(quantile(v, q)) for q in qs), " / ")
frac(v) = isempty(v) ? NaN : count(identity, v) / length(v)

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
const EPS64 = nf_eps(Float64)
say("# naff-fma verification, julia $(VERSION), threads $(Threads.nthreads()), base seed $NF_SEED (verify offset 61000), time $(round(time())) (unix)")
say("# timings in this log were taken on a machine that was not otherwise idle")

# ---------------------------------------------------------------------------------------------------
# The package suite at full size; its @info lines are forwarded into this log (the gated ratio keys are listed in GATES)
const GATES = ("AC10d" => ((:worst, 1.0),), "AC7b" => ((:worst_ratio, 4.0),), "D81a" => ((:worst, 1.0),))
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
    rest = join(("$k = $(kw[k] isa AbstractFloat ? fmt(kw[k]) : kw[k])" for (k, v) in kw if !any(k == g[1] for g in gated)), "; ")
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
            @testset verbose = true "naff-fma verify" begin
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
    emit("SUITE", "package suite at full size (NF_FULL_RUN), assertions passed (fail + error = $nbad); wall time $(round(t; digits = 1)) s", npass, "0 fail, 0 error", pf(ok && nbad == 0 && npass > 0))
end

# ---------------------------------------------------------------------------------------------------
# AC14 (a): the independent BigFloat(256) orbit of the standard map in the project form (p' = p + K sin q, q' = q + p'), from the same exactly widened start as the module's Float64 orbit. The module's
# realized series (Float64 wrapped orbit) against the series of the BigFloat orbit, both through the module's estimators: E1 on exp(i q) (mean removed, 2N samples) and E2 on the lifted increments
# (the stream result of the module against NF6 evaluated on the exact increments of the BigFloat orbit). (b): ChaosTools estimate_period (periodogram) on cos(q_j) of the same regular starts.
function big_orbit(K::Float64, s0::SVector{4, Float64}, n::Int)
    return setprecision(BigFloat, NF_BITS) do
        tp = 2 * BigFloat(pi)
        Kb = BigFloat(K)
        X = BigFloat(s0[1]) + tp * BigFloat(s0[3])
        P = BigFloat(s0[2]) + tp * BigFloat(s0[4])
        Xs = Vector{BigFloat}(undef, n)
        for i in 1:n
            P = P + Kb * sin(X)
            X = X + P
            Xs[i] = X
        end
        Xs
    end
end

section("AC14") do
    N, m = 1000, 20
    c = nf_case(:half, Float64, m, 0)
    sys = c.sys
    r = WeightedPhaseAdvance(; window = N, observable = LiftedAngleObservable(1, 3))
    res = propagate(sys, c.x0, 2N + 1; k = 1, tau = 1, seed = vseed(14, 1), reducers = (r, Benettin()), backend = :cpu)
    rows = NamedTuple[]
    for j in 1:m
        st = orbit_states(sys, c.x0[j], 2N + 1)
        zmod = ComplexF64[obs_cis(s) for s in st[1:2N]]
        Xs = big_orbit(0.5, c.x0[j], 2N + 1)
        zbig = setprecision(BigFloat, NF_BITS) do
            ComplexF64[complex(Float64(cos(x)), Float64(sin(x))) for x in Xs[1:2N]]
        end
        fm, fb = frequency_map(zmod, N; mean_removal = true), frequency_map(zbig, N; mean_removal = true)
        fm0, fb0 = frequency_map(zmod, N; mean_removal = false), frequency_map(zbig, N; mean_removal = false)
        incr = setprecision(BigFloat, NF_BITS) do
            [Xs[i + 1] - Xs[i] for i in 1:2N]
        end
        ref = nf_E2_ref(incr, N)
        f = res.finals[j][1]
        e1 = [nf_cd(fm.w1[1].nu, fb.w1[1].nu), nf_cd(fm.w2[1].nu, fb.w2[1].nu)]
        e1n = [nf_cd(fm0.w1[1].nu, fb0.w1[1].nu), nf_cd(fm0.w2[1].nu, fb0.w2[1].nu)]
        e2 = [Float64(abs(nf_big(f.nu1[1]) - ref[1])), Float64(abs(nf_big(f.nu2[1]) - ref[2]))]
        # (b) ChaosTools period of cos(q_j) over the first window
        v = Float64[cos(s[1]) for s in st[1:N]]
        per = try
            CT.estimate_period(v .- mean(v), :periodogram)
        catch e
            NaN
        end
        push!(rows, (j = j, e1 = Float64.(e1), e1n = Float64.(e1n), e2 = e2, nu = fm.w1[1].nu, period = per, flagz = fm.w1[1].near_zero, flagm = fm.w1[1].mirror, delta_mod = fm.delta, delta_big = fb.delta))
    end
    DATA[:ac14] = rows
    w1, w2 = maximum(maximum(r.e1) for r in rows), maximum(maximum(r.e2) for r in rows)
    emit("AC14a", "regular class K = 0.5, $m starts, N = $N, Float64 (mean removal on): max over starts and both windows of the circle distance |nu_E1(module series) - nu_E1(BigFloat orbit)|", w1, "<= 1e-12", pf(w1 <= 1e-12))
    emit("AC14a", "same, E2 (stream, lifted increments, exponential weight) against NF6 on the exact increments of the BigFloat orbit: max |nu_E2 - nu_E2,ref|", w2, "<= 1e-12", pf(w2 <= 1e-12))
    info("AC14a", "realized ratios to the 1e-12 bound: E1 (mean removal on); E1 (mean removal off, max); E2", "$(fmt(w1 / 1e-12)); $(fmt(maximum(maximum(r.e1n) for r in rows) / 1e-12)); $(fmt(w2 / 1e-12))")
    info("AC14a", "median over starts of the E1 / E2 difference (both windows)", "$(fmt(rmedian(vcat([r.e1 for r in rows]...)))); $(fmt(rmedian(vcat([r.e2 for r in rows]...))))")
    info("AC14a", "starts with the near-zero flag or the mirror flag set in window 1 (E1, mean removal on); frequency map delta of module series / of BigFloat series: max", "$(count(r -> r.flagz, rows)); $(count(r -> r.flagm, rows)); $(fmt(maximum(r.delta_mod for r in rows))) / $(fmt(maximum(r.delta_big for r in rows)))")
    # (b)
    ok, devs, bad = 0, Float64[], String[]
    for r in rows
        isfinite(r.period) || (push!(bad, "start $(r.j): estimate_period gave $(r.period)"); continue)
        dev = abs(1 / r.period - abs(r.nu)) * N      # in Fourier bins of the window (1/N)
        push!(devs, dev)
        dev <= 1 ? (ok += 1) : push!(bad, "start $(r.j): period $(fmt(r.period)) -> 1/P = $(fmt(1 / r.period)), |nu_E1| = $(fmt(abs(r.nu))), deviation $(fmt(dev)) bins")
    end
    info("AC14b", "ChaosTools 3.5.4 estimate_period(cos(q_j) minus mean, :periodogram) on window 1 (N = 1000, $m regular starts): starts with | 1/P - |nu_E1| | within one Fourier bin (1/N); listed below if not", "$ok of $m")
    # the periodogram returns the strongest harmonic of cos(q_j), the module the fundamental of exp(i q_j): starts where 1/P is within one bin of k |nu_E1| (folded into [0, 1/2]) for k = 1..4
    kmatch = [begin
        f = 1 / r.period
        ks = [k for k in 1:4 if abs(f - abs(nf_red(k * abs(r.nu)))) * N <= 1]
        isempty(ks) ? 0 : first(ks)
    end for r in rows if isfinite(r.period)]
    info("AC14b", "starts where 1/P is within one bin of k |nu_E1| (folded into [0, 1/2]): k = 1; 2; 3; 4; none", "$(count(==(1), kmatch)); $(count(==(2), kmatch)); $(count(==(3), kmatch)); $(count(==(4), kmatch)); $(count(==(0), kmatch))")
    isempty(devs) || info("AC14b", "deviation in bins | 1/P - |nu_E1| | N: quantiles 5 / 50 / 95 % and max", "$(quants(devs, (0.05, 0.5, 0.95))); $(fmt(maximum(devs)))")
    for b in bad
        info("AC14b", "outside one bin / not run", b)
    end
    # (c) Python packages
    pyok = [m_ for m_ in ("nafflib", "PyNAFF") if success(pipeline(`python -c "import $m_"`; stderr = devnull))]
    info("AC14c", "Python NAFF packages (nafflib, PyNAFF) importable in the default Python (3.12): installed ones", isempty(pyok) ? "none (not run; spec: only if the owner installs them)" : join(pyok, ", "))
end

# ---------------------------------------------------------------------------------------------------
# AC12(b): 1 versus N threads in child processes (output hashes), chunk sizes 64 and 1000, E1 (series path over orbits) and E2 (stream)
section("AC12") do
    hashes = Dict{Int, Dict{String, String}}(Threads.nthreads() => thread_hashes())
    for nt in (1, 2, 4)
        haskey(hashes, nt) && continue
        cmd = addenv(`$(Base.julia_cmd()) -t $nt --project=$(dirname(Base.active_project())) $(@__FILE__)`, "NF_CHILD" => "1")
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
    for tag in ("K=0.5", "K=2pi"), est in ("E1", "E2")
        ok_all &= length(unique([hashes[Threads.nthreads()]["$est $tag chunk=$chunk"] for chunk in (64, 1000)])) == 1
    end
    emit("AC12", "(b) E1 (frequency_map over 4096 orbits, threads over orbits, dynamic) and E2 (propagate, chain {E2, Benettin}) outputs bitwise identical for thread counts $(sort(collect(keys(hashes)))) (hashes from separate processes) and chunk sizes 64 and 1000: $(length(keys_)) cases, K in {0.5, 2 pi}, N = 1000; the transform is the in-house single-threaded radix-2 FFT, so the library thread count is moot", ok_all, "all equal", pf(ok_all))
end

# ---------------------------------------------------------------------------------------------------
# Shared by R1, R3, R4: the K = 0.5 regular class (L1 S7 screening: Benettin X(1e5) <= 5e-4 on the lifted map, tau = 10) of m starts out of 3 m candidates, and uniform starts of the other classes.
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

# E1 on one class: threads over orbits, the orbit generated once to 2 Nmax samples, both observables (exp(i q) and q + i p, wrapped, window means removed), the pair of windows per N.
function e1_panel(sys, x0, Ns)
    m, Nmax = length(x0), maximum(Ns)
    mk() = Dict(N => (D = fill(NaN, m), nu1 = fill(NaN, m), mir = fill(false, m), nz = fill(false, m), floor_ = fill(false, m), tie = fill(false, m), ok = fill(false, m)) for N in Ns)
    out = (cis = mk(), planar = mk())
    Threads.@threads :dynamic for j in 1:m
        st = orbit_states(sys, x0[j], 2Nmax)
        for (name, obs) in ((:cis, obs_cis), (:planar, obs_planar))
            z = ComplexF64[obs(s) for s in st]
            for N in Ns
                fm = frequency_map(@view(z[1:2N]), N; mean_removal = true)
                o = getproperty(out, name)[N]
                o.D[j] = fm.D; o.nu1[j] = fm.w1[1].nu; o.ok[j] = fm.status[1] === :ok
                o.mir[j] = fm.w1[1].mirror | fm.w2[1].mirror; o.nz[j] = fm.w1[1].near_zero | fm.w2[1].near_zero
                o.floor_[j] = fm.below_floor; o.tie[j] = fm.w1[1].tie | fm.w2[1].tie
            end
        end
    end
    return out
end
# E2 on one class: one stream run with several horizons (windows), lifted and planar observables, chain {E2 lifted, E2 planar, Benettin} (the Benettin reducer renormalises the single column)
function e2_panel(sys, x0, Ns; seed)
    rl = WeightedPhaseAdvance(; window = Tuple(Ns), observable = LiftedAngleObservable(1, 3))
    rp = WeightedPhaseAdvance(; window = Tuple(Ns), observable = PlanarObservable(1, 2))
    res = propagate(sys, x0, 2 * maximum(Ns) + 1; k = 1, tau = 1, seed = seed, reducers = (rl, rp, Benettin()), backend = :cpu)
    get(i, h) = (D = [Float64(f[i].D[h]) for f in res.finals], floor_ = [f[i].below_floor[h] for f in res.finals], nu1 = [Float64(f[i].nu1[h]) for f in res.finals])
    return (lifted = Dict(N => get(1, h) for (h, N) in enumerate(Ns)), planar = Dict(N => get(2, h) for (h, N) in enumerate(Ns)))
end

# ---------------------------------------------------------------------------------------------------
# R1 (and R3, observables and flags): class separation on the 2D panel, E1 and E2, both observables, 1000 starts per class, N in {250, 1000, 4000, 16000, 64000}.
section("R1") do
    m = 1000
    Ns = (250, 1000, 4000, 16000, 64000)
    reg = regular_class(m)
    info("R1", "K = 0.5 regular class: candidates screened, regular (Benettin X(1e5) <= 5e-4), taken", "$(reg.ncand), $(reg.nreg), $m")
    R = Dict{Tuple{Symbol, Symbol}, Any}()      # (class, estimator) => per-N record
    for (label, K) in R_CLASSES
        sys = WrappedStandardMap(K)
        x0 = [torus_state(s) for s in class_starts(label, m)]
        t1 = @elapsed e1 = e1_panel(sys, x0, Ns)
        t2 = @elapsed e2 = e2_panel(sys, x0, Ns; seed = vseed(42, 1))
        R[(label, :e1c)] = e1.cis; R[(label, :e1p)] = e1.planar; R[(label, :e2l)] = e2.lifted; R[(label, :e2p)] = e2.planar
        say(@sprintf("# R1 class %s: E1 panel %.1f s, E2 panel %.1f s", label, t1, t2))
    end
    DATA[:R1] = R
    names = (:e1c => "E1 exp(iq)", :e1p => "E1 q+ip", :e2l => "E2 lifted", :e2p => "E2 planar")
    for (est, nm) in names, N in Ns
        for (label, K) in R_CLASSES
            d = filter(isfinite, R[(label, est)][N].D)
            info("R1", "$nm, N = $N, K = $(fmt(K)) ($label), $m starts: quantiles 5 / 25 / 50 / 75 / 95 % of D (non-finite excluded: $(m - length(d)))", quants(d))
        end
        q = median(filter(isfinite, R[(:half, est)][N].D))
        regD = filter(isfinite, R[(:half, est)][N].D)
        fl_reg = frac(R[(:half, est)][N].floor_)
        line = String[]
        for (label, K) in R_CLASSES[2:end]
            d = R[(label, est)][N].D
            fl = R[(label, est)][N].floor_
            ok = isfinite.(d)
            push!(line, "K=$(fmt(K)): below-median-or-floor $(fmt(mean((d[ok] .< q) .| fl[ok]))) (floor flag alone $(fmt(frac(fl))))")
        end
        info("R1", "$nm, N = $N: regular-class median of D $(fmt(q)); fraction of chaotic starts with D below it or with the below-floor flag; regular class: fraction above D = -6 $(fmt(mean(regD .> -6))), floor-flag share $(fmt(fl_reg))", join(line, "; "))
    end
    # modes of D (Bazzani2023 pp.13, 20): histogram of D over all 4000 starts, unit bins from -16 to 2, E1 exp(iq) and E2 lifted
    edges = collect(-16.0:1.0:2.0)
    for (est, nm) in (:e1c => "E1 exp(iq)", :e2l => "E2 lifted"), N in (1000, 4000, 16000)
        pooled = vcat((filter(isfinite, R[(label, est)][N].D) for (label, _) in R_CLASSES)...)
        h = [count(x -> lo <= x < lo + 1, pooled) for lo in edges[1:end-1]]
        hs = [(get(h, i - 1, 0) + 2 * h[i] + get(h, i + 1, 0)) / 4 for i in eachindex(h)]       # 3-point smoothing
        peaks = [i for i in eachindex(hs) if hs[i] >= 0.03 * maximum(hs) && hs[i] > get(hs, i - 1, 0) && hs[i] >= get(hs, i + 1, 0)]
        per = join(("K=$(fmt(K)):" * join([count(x -> lo <= x < lo + 1, filter(isfinite, R[(label, est)][N].D)) for lo in edges[1:end-1]], ",") for (label, K) in R_CLASSES), " ")
        info("R1", "$nm, N = $N: histogram of D over the 4000 starts, unit bins from -16 to 2 (counts; smoothed local maxima at bin lower edges: $(join(edges[peaks], ", ")) -> $(length(peaks)) modes)", join(h, ","))
        info("R1", "$nm, N = $N: the same per class (K = 0.5; 2 pi; 5; 10)", per)
    end
    # R3: flags
    for (est, nm) in (:e1c => "E1 exp(iq)", :e1p => "E1 q+ip"), N in (1000, 4000)
        info("R3", "$nm (mean removed), N = $N: share of starts with the mirror flag; with the near-zero flag; with a tie (any window), per class K = 0.5; 2 pi; 5; 10", join(("$(fmt(frac(R[(l, est)][N].mir))) / $(fmt(frac(R[(l, est)][N].nz))) / $(fmt(frac(R[(l, est)][N].tie)))" for (l, _) in R_CLASSES), "; "))
    end
    for N in (1000, 4000)
        for (l, K) in R_CLASSES
            a, b, c = (filter(isfinite, R[(l, e)][N].D) for e in (:e1c, :e1p, :e2l))
            info("R3", "N = $N, K = $(fmt(K)): median D under exp(iq) (E1) / q+ip (E1) / lifted increments (E2) / planar (E2)", "$(fmt(median(a))) / $(fmt(median(b))) / $(fmt(median(c))) / $(fmt(median(filter(isfinite, R[(l, :e2p)][N].D))))")
        end
    end
    # cross-tabulation with the labels of L1 (n X above ln n), FLI_B (above ln n), MEGNO (Ybar_1 above 3), RLI (RLI / Dx_r above 1 at Dx = 1e-12) and the stretching numbers D (above the regular-class
    # median at the same n), at n = 2N for N = 1000 and 4000 (thresholds of Maffione2011a as used in the stretching and rli runs, not pre-registered). Cells: (label chaotic, D low) / (chaotic, high) /
    # (regular, low) / (regular, high), D low = at most the regular-class median at the same N.
    labels = Dict{Tuple{Symbol, Int}, NamedTuple}()
    seedR = vseed(42, 1)
    for (label, K) in R_CLASSES
        base = WrappedStandardMap(K)
        P = PairMap(base)
        starts = class_starts(label, m)
        xw = [torus_state(s) for s in starts]
        xis = [SVector{2, Float64}(initial_deviations(seedR, j, Val(2), Val(1), Float64)[:, 1]) for j in 1:m]
        ps = [pair_start(base, starts[j], xis[j], 1e-12; orbit = j) for j in 1:m]
        z0 = [p.z for p in ps]
        W0 = [SMatrix{4, 2, Float64}(xi[1], xi[2], 0.0, 0.0, 0.0, 0.0, xi[1], xi[2]) for xi in xis]
        for N in (1000, 4000)
            n = 2N
            r = propagate(P, z0, n; k = 2, tau = 1, seed = seedR, reducers = (RLI(),), W0 = W0, backend = :cpu)
            l1 = propagate(base, xw, n; k = 1, tau = 1, seed = seedR, reducers = (Benettin(),), backend = :cpu)
            fm = propagate(base, xw, n; k = 2, tau = 1, seed = seedR, reducers = (FLI(replace = false), MEGNO()), backend = :cpu)
            sn = propagate(base, xw, n; k = 2, tau = 1, seed = seedR, reducers = (StretchingNumbers(; range = stretching_range(base), bins = 1000, interval = 100),), backend = :cpu)
            labels[(label, N)] = (RLI = [r.finals[j][1].RLI / ps[j].Dx_r > 1.0 for j in 1:m], L1 = [f[1].lambda[1] * n > log(n) for f in l1.finals],
                                  FLI = [f[1].fli > log(n) for f in fm.finals], MEGNO = [f[2].Ybar[1] > 3 for f in fm.finals], SNdist = [Float64(f[1].D) for f in sn.finals])
        end
    end
    for N in (1000, 4000), (est, nm) in (:e1c => "E1 exp(iq)", :e2l => "E2 lifted")
        regD = filter(isfinite, R[(:half, est)][N].D)
        q = median(regD)
        for lab in (:L1, :FLI, :MEGNO, :RLI, :SNdist)
            lv = vcat((lab === :SNdist ? labels[(l, N)].SNdist .> median(labels[(:half, N)].SNdist) : getproperty(labels[(l, N)], lab) for (l, _) in R_CLASSES)...)
            dv = vcat((R[(l, est)][N].D for (l, _) in R_CLASSES)...)
            ok = isfinite.(dv)
            low = (dv .<= q)
            lv, low = lv[ok], low[ok]
            info("R1", "n = $(2N) (N = $N), $nm D against label $lab ($(lab === :SNdist ? "stretching D above its regular-class median" : "chaotic")), all $(count(ok)) finite starts of 4 classes, regular-class median of D $(fmt(q)): (label chaotic, D low) / (chaotic, high) / (regular, low) / (regular, high); fraction labelled chaotic per class (K = 0.5; 2 pi; 5; 10)",
                 "$(count(lv .& low)) / $(count(lv .& .!low)) / $(count(.!lv .& low)) / $(count(.!lv .& .!low)); $(join((fmt(frac(lab === :SNdist ? labels[(l, N)].SNdist .> median(labels[(:half, N)].SNdist) : getproperty(labels[(l, N)], lab))) for (l, _) in R_CLASSES), "; "))")
        end
    end
end

# ---------------------------------------------------------------------------------------------------
# R2, the Eq. (79) part: Laskar2003 Eq. (79) test function F_1(t) = 1/(1 + (1/2) exp(i t) + (1/4) exp(-i omega t)), omega = 2.02, and its truncation F_2 to 50 terms, E1 of window order p = 0..5 at the
# spans log10(T / (2 pi)) from 2.2 to 4 (window [-T, T], sampling interval h = 1/4 time units, samples at the centres t_j = -T + (j + 1/2) h). The first frequency is the largest-amplitude term, the constant term
# 1 (reader inference of the spec): the reference is the exact value 0 and, for the rounding part, the same procedure in Double64. Errors in rad per unit time: 2 pi |nu| / h.
section("R2eq79") do
    h = 0.25
    om = Double64(202) / 100
    xs = [2.2, 2.5, 2.8, 3.1, 3.4, 3.7, 4.0]
    function series(kind, x)
        T = 2 * Double64(pi) * Double64(10)^x
        N = round(Int, 2 * Float64(T) / h)
        N += isodd(N) ? 1 : 0
        Tt = Double64(N) * Double64(h) / 2
        out = Vector{Complex{Double64}}(undef, N)
        for j in 1:N
            t = -Tt + (Double64(j - 1) + Double64(1) / 2) * Double64(h)
            u = Complex(cos(t), sin(t)) / 2
            v = Complex(cos(om * t), -sin(om * t)) / 4
            s = u + v
            out[j] = kind === :F1 ? 1 / (1 + s) : (q = Complex{Double64}(0); pw = Complex{Double64}(1); for _ in 0:49; q += pw; pw *= -s; end; q)
        end
        return N, out
    end
    cases = [(kind, x, p) for kind in (:F1, :F2) for x in xs for p in 0:5]
    res = Vector{Any}(undef, length(cases))
    cache = Dict{Tuple{Symbol, Float64}, Any}()
    for kind in (:F1, :F2), x in xs
        cache[(kind, x)] = series(kind, x)
    end
    Threads.@threads :dynamic for i in eachindex(cases)
        kind, x, p = cases[i]
        N, z = cache[(kind, x)]
        z64 = ComplexF64[complex(Float64(real(a)), Float64(imag(a))) for a in z]
        w64 = window_frequency(z64; order = p)
        wd = window_frequency(z; order = p)
        res[i] = (kind = kind, x = x, p = p, N = N, e64 = 2pi * abs(w64.nu) / h, ed = 2pi * abs(Float64(wd.nu)) / h, rnd = 2pi * abs(w64.nu - Float64(wd.nu)) / h, nu = Float64(wd.nu), ok = (w64.status === :ok, wd.status === :ok, wd.refined, w64.refined))
    end
    DATA[:eq79] = res
    a0 = Dict(:F1 => (1.98, 3.94, 5.82, 7.64, 8.89, 9.15), :F2 => (1.91, 3.98, 5.89, 8.69, 11.42, 13.50))
    for kind in (:F1, :F2), p in 0:5
        sel = [r for r in res if r.kind === kind && r.p == p]
        sort!(sel, by = r -> r.x)
        lx = [r.x + log10(2pi) for r in sel]          # log10 T
        function fit(key, floor_)
            idx = [i for i in eachindex(sel) if getproperty(sel[i], key) > floor_]
            length(idx) >= 3 ? -lslope(lx[idx], log10.([getproperty(sel[i], key) for i in idx])) : NaN
        end
        sl_d, sl_64 = fit(:ed, 1e-24), fit(:e64, 1e-13)
        info("R2eq79", "$kind, p = $p (printed theory 2p+2 = $(2p + 2); printed fitted a_0 / b_0 = $(a0[kind][p + 1])): error of the first frequency in rad per time at log10(T/2pi) = $(join(xs, ", ")): Double64 | Float64", "$(join((fmt(r.ed) for r in sel), " ")) | $(join((fmt(r.e64) for r in sel), " "))")
        info("R2eq79", "$kind, p = $p: fitted slope -d log10(err) / d log10(T), Double64 (points above 1e-24) | Float64 (points above 1e-13); Float64 minus Double64 rounding part: max | all windows ok / refined", "$(fmt(sl_d)) | $(fmt(sl_64)); $(fmt(maximum(r.rnd for r in sel))); $(all(r -> all(r.ok), sel))")
    end
end

# ---------------------------------------------------------------------------------------------------
# R4: sticky orbits. The sticky ensemble of the study does not exist yet (the P2 study is not designed); an operational one is built here: K = 0.9, uniform candidates, labelled by the L1 rule n X above ln n at two
# horizons: regular at n = 1e4 (n X <= ln n) and chaotic at n = 1e6 (n X > ln n). D (E1 exp(iq), E2 lifted) at N in {250, 1000, 4000, 16000} for the first window pair from the start; "leaves the regular range" =
# D above the 99th percentile of the K = 0.5 regular-class D at the same N.
section("R4") do
    K = 0.9
    ncand = 6000
    rng = StableRNG(vseed(44, 1))
    cand = [SVector(2pi * rand(rng), 2pi * rand(rng)) for _ in 1:ncand]
    xs_ = [torus_state(s) for s in cand]
    sys = WrappedStandardMap(K)
    r4 = propagate(sys, xs_, 10^4; k = 1, tau = 10, seed = vseed(44, 2), reducers = (Benettin(),), backend = :cpu)
    r6 = propagate(sys, xs_, 10^6; k = 1, tau = 10, seed = vseed(44, 2), reducers = (Benettin(),), backend = :cpu)
    lab4 = [f[1].lambda[1] * 10^4 > log(10^4) for f in r4.finals]
    lab6 = [f[1].lambda[1] * 10^6 > log(10^6) for f in r6.finals]
    sticky = findall(.!lab4 .& lab6)
    info("R4", "K = 0.9, $ncand uniform candidates: label chaotic at n = 1e4 (n X > ln n); chaotic at n = 1e6; sticky = regular at 1e4 and chaotic at 1e6; the other chaotic at 1e4 and 1e6", "$(count(lab4)); $(count(lab6)); $(length(sticky)); $(count(lab4 .& lab6))")
    length(sticky) >= 1 || return
    Ns = (250, 1000, 4000, 16000)
    st = e1_panel(sys, xs_[sticky], Ns)
    e2 = e2_panel(sys, xs_[sticky], Ns; seed = vseed(44, 3))
    DATA[:R4] = (n = length(sticky), cis = st.cis, e2 = e2.lifted)
    R = DATA[:R1]
    m = length(sticky)
    for (est, nm, rec) in ((:e1c, "E1 exp(iq)", st.cis), (:e2l, "E2 lifted", e2.lifted))
        first_out = fill(0, m)
        for N in Ns
            regD = filter(isfinite, R[(:half, est)][N].D)
            q99, q50 = quantile(regD, 0.99), median(regD)
            d = rec[N].D
            ok = isfinite.(d)
            info("R4", "$nm, N = $N: $(count(ok)) finite of $m sticky orbits: quantiles 5 / 25 / 50 / 75 / 95 % of D; regular-class median / 99 % quantile $(fmt(q50)) / $(fmt(q99)); share with D below the regular median; share above the 99 % quantile", "$(quants(d[ok])); $(fmt(mean(d[ok] .< q50))); $(fmt(mean(d[ok] .> q99)))")
            for j in 1:m
                first_out[j] == 0 && ok[j] && d[j] > q99 && (first_out[j] = N)
            end
        end
        info("R4", "$nm: window N at which D first leaves the regular range (above the regular-class 99 % quantile), counts over the $m sticky orbits for N = 250; 1000; 4000; 16000; never in the list", "$(join((count(==(N), first_out) for N in Ns), "; ")); $(count(==(0), first_out))")
    end
end

# ---------------------------------------------------------------------------------------------------
# R5: cost on the CPU (the device part is in verify/gpu.jl, tag naff). Wall seconds per orbit and bytes per orbit: the stream alone (Benettin), the stream with the E2 reader, and for E1 the orbit pass (recording
# the observable at every handoff, per orbit) and the transform (frequency_map of the pair) separately, at N in {1e3, 1e4, 1e5}.
section("R5") do
    sys = WrappedStandardMap(2pi)
    for (N, norb) in ((1000, 256), (10_000, 64), (100_000, 8))
        rng = StableRNG(vseed(45, 1))
        x0 = [torus_state(SVector(2pi * rand(rng), 2pi * rand(rng))) for _ in 1:norb]
        n = 2N + 1
        r = WeightedPhaseAdvance(; window = N, observable = LiftedAngleObservable(1, 3))
        run_(rs) = propagate(sys, x0, n; k = 1, tau = 1, seed = vseed(45, 2), reducers = rs, backend = :cpu)
        run_((Benettin(),)); run_((r, Benettin()))
        t_stream = minimum(@elapsed(run_((Benettin(),))) for _ in 1:3)
        t_e2 = minimum(@elapsed(run_((r, Benettin()))) for _ in 1:3)
        info("R5", "CPU ($(Threads.nthreads()) threads), N = $N ($n steps), $norb orbits, K = 2 pi: stream alone (Benettin) / stream with E2 reader: wall seconds per orbit; ns per orbit-step; overhead of E2", "$(fmt(t_stream / norb)) / $(fmt(t_e2 / norb)); $(fmt(t_stream / (norb * n) * 1e9)) / $(fmt(t_e2 / (norb * n) * 1e9)); $(fmt(t_e2 / t_stream - 1))")
        # E1: orbit pass (single thread, one orbit at a time), transform (single thread)
        z = Vector{ComplexF64}(undef, 2N)
        function pass!(z, x)
            s = x
            @inbounds for i in 1:2N
                s = forward(sys, s)
                z[i] = obs_cis(s)
            end
            return z
        end
        pass!(z, x0[1])
        t_pass = minimum(@elapsed(foreach(x -> pass!(z, x), x0[1:min(norb, 8)])) for _ in 1:3) / min(norb, 8)
        zs = [copy(pass!(z, x)) for x in x0[1:min(norb, 8)]]
        frequency_map(zs[1], N; mean_removal = true)
        t_tr = minimum(@elapsed(foreach(zz -> frequency_map(zz, N; mean_removal = true), zs)) for _ in 1:3) / length(zs)
        alloc = @allocated frequency_map(zs[1], N; mean_removal = true)
        info("R5", "CPU single thread, N = $N: E1 orbit pass (2N samples of exp(i q)) / transform (frequency_map of the pair, mean removal on): seconds per orbit; ratio transform / orbit pass; series bytes 32 N (Float64: 2N complex) and bytes allocated by the transform", "$(fmt(t_pass)) / $(fmt(t_tr)); $(fmt(t_tr / t_pass)); $(32N) B; $(alloc) B")
    end
    info("R5", "reducer state bytes (Float64 / Double64, one horizon; per extra horizon)", "$(sizeof(typeof(Chaoslab.reducer_init(WeightedPhaseAdvance(; window = 1000, observable = LiftedAngleObservable(1, 3)), Float64, Val(2), Val(1))))) / $(sizeof(typeof(Chaoslab.reducer_init(WeightedPhaseAdvance(; window = 1000, observable = LiftedAngleObservable(1, 3)), Double64, Val(2), Val(1))))); +64 / +128")
end

# ---------------------------------------------------------------------------------------------------
# R6: a 2D map with a planar elliptic point. Y2 does not exist in the library; the substitute is the elliptic fixed point (pi, 0) of the standard map (stable for 0 < K < 4). Observable (q - pi) + i p with p reduced to
# (-pi, pi]; starts uniform in the box |q - pi| < 1, |p| < 1; K in {0.5, 1}, 1000 starts each, N in {1000, 4000, 16000}; E1 (mean removal off, the point is the origin) and E2 (pure function on the same series).
section("R6") do
    m = 1000
    Ns = (1000, 4000, 16000)
    for K in (0.5, 1.0)
        sys = WrappedStandardMap(K)
        rng = StableRNG(vseed(46, K == 0.5 ? 1 : 2))
        x0 = [torus_state(SVector(pi + 2 * rand(rng) - 1, 2 * rand(rng) - 1)) for _ in 1:m]
        Nmax = maximum(Ns)
        D1 = Dict(N => fill(NaN, m) for N in Ns); D2 = Dict(N => fill(NaN, m) for N in Ns); F1 = Dict(N => fill(false, m) for N in Ns)
        Threads.@threads :dynamic for j in 1:m
            st = orbit_states(sys, x0[j], 2Nmax + 1)
            z = ComplexF64[complex(s[1] - pi, (s[2] > pi ? s[2] - 2pi : s[2])) for s in st]
            for N in Ns
                fm = frequency_map(@view(z[1:2N]), N)
                D1[N][j] = fm.D; F1[N][j] = fm.below_floor
                nu = weighted_phase_frequency(phase_increments(z[1:(2N + 1)]), N)
                D2[N][j] = frequency_index(circle_distance(nu[1], nu[2]))[1]
            end
        end
        for N in Ns
            for (nm, D) in (("E1", D1[N]), ("E2 (pure)", D2[N]))
                d = filter(isfinite, D)
                info("R6", "standard map K = $K, elliptic point (pi, 0), $nm, N = $N, $m starts in |q - pi|, |p| < 1: quantiles 5 / 25 / 50 / 75 / 95 % of D; share of starts with D above 1e-5 (D > -5); share at the floor (D <= log10 32 eps = $(fmt(log10(32 * EPS64))) + 1e-9)", "$(quants(d)); $(fmt(mean(d .> -5))); $(fmt(mean(d .<= log10(32 * EPS64) + 1e-9)))")
            end
        end
        info("R6", "standard map K = $K: the E1 below-floor flag share at N = 1000; 4000; 16000", "$(fmt(frac(F1[1000]))); $(fmt(frac(F1[4000]))); $(fmt(frac(F1[16000])))")
    end
end

# ---------------------------------------------------------------------------------------------------
# Figures
section("figure") do
    # spectra of a regular and a chaotic orbit (window 1, N = 1000, exp(i q) with the window mean removed), the found frequency of the module, and the map of D over the panel of the standard map at K = 0.9
    N = 1000
    reg = regular_class(100).starts[1]
    chs = SVector(2pi * 0.3141592653589793, 2pi * 0.7182818284590452)
    cases = ((WrappedStandardMap(0.5), reg, "regular orbit, K = 0.5, start ($(fmt(reg[1])), $(fmt(reg[2])))", :seagreen),
             (WrappedStandardMap(2pi), chs, "chaotic orbit, K = 2 pi, start ($(fmt(chs[1])), $(fmt(chs[2])))", :darkorange))
    fig = CM.Figure(size = (1800, 560))
    chi = window_weights(Float64, N)
    nus = collect(range(-0.5, 0.5; length = 4001))
    for (i, (sys, s, ttl, col)) in enumerate(cases)
        st = orbit_states(sys, torus_state(s), 2N)
        z = ComplexF64[obs_cis(a) for a in st]
        fm = frequency_map(z, N; mean_removal = true)
        zw = z[1:N] .- mean(z[1:N])
        phi = [abs(sum(chi[j] * zw[j] * cis(-2pi * nu * (j - 1)) for j in 1:N)) / sum(chi) for nu in nus]
        ax = CM.Axis(fig[1, i]; xlabel = "frequency (cycles per sample)", ylabel = "|phi(nu)| (Hanning window, N = 1000)", yscale = log10, title = ttl)
        CM.lines!(ax, nus, max.(phi, 1e-9); color = col)
        CM.vlines!(ax, [fm.w1[1].nu]; color = :black, linestyle = :dash, label = "found nu = $(fmt(fm.w1[1].nu)), window 2: $(fmt(fm.w2[1].nu)), D = $(fmt(fm.D))")
        CM.ylims!(ax, 1e-7, nothing)
        CM.axislegend(ax; position = :ct, labelsize = 11)
    end
    # D map over the panel, K = 0.9
    K = 0.9
    g = 160
    sys = WrappedStandardMap(K)
    qs = range(0, 2pi; length = g + 1)[1:g]
    D = fill(NaN, g, g)
    Threads.@threads :dynamic for a in 1:g
        for b in 1:g
            x0 = torus_state(SVector(qs[a], qs[b]))
            st = orbit_states(sys, x0, 2N)
            fm = frequency_map(ComplexF64[obs_cis(s) for s in st], N; mean_removal = true)
            D[a, b] = fm.D
        end
    end
    axm = CM.Axis(fig[1, 3]; xlabel = "q", ylabel = "p", title = "D = log10 delta (E1, exp(iq), N = 1000) over the panel, K = 0.9 ($(g)x$(g) starts)", aspect = 1)
    hm = CM.heatmap!(axm, collect(qs), collect(qs), D; colormap = :viridis, colorrange = (-14.5, 0.0))
    CM.Colorbar(fig[1, 4], hm; label = "D")
    CM.save(joinpath(OUT, "naff-fma-spectra.png"), fig)
    say("# figure saved: naff-fma-spectra.png")
    # the oracle figure: AC14 (a) and (b), and R2 Eq. (79)
    fig2 = CM.Figure(size = (1800, 560))
    if haskey(DATA, :ac14)
        rows = DATA[:ac14]
        ax1 = CM.Axis(fig2[1, 1]; yscale = log10, xlabel = "regular-class start", ylabel = "difference of the frequency (cycles per sample)", title = "AC14(a): module series against the BigFloat(256) orbit (bound 1e-12)")
        CM.scatter!(ax1, 1:length(rows), [max(maximum(r.e1), 1e-19) for r in rows]; color = :seagreen, marker = :circle, label = "E1 exp(iq), max of 2 windows")
        CM.scatter!(ax1, 1:length(rows), [max(maximum(r.e2), 1e-19) for r in rows]; color = :purple, marker = :utriangle, label = "E2 lifted, max of 2 windows")
        CM.hlines!(ax1, [1e-12]; color = :black, linestyle = :dash)
        CM.axislegend(ax1; position = :rb, labelsize = 11)
        ax2 = CM.Axis(fig2[1, 2]; xscale = log10, yscale = log10, xlabel = "|nu_E1| (module)", ylabel = "1 / estimate_period(cos q)", title = "AC14(b): ChaosTools period against the module (dotted: one bin, 1/N)")
        ok = [r for r in rows if isfinite(r.period)]
        CM.scatter!(ax2, [abs(r.nu) for r in ok], [1 / r.period for r in ok]; color = :darkorange)
        lo_, hi_ = extrema(abs(r.nu) for r in ok)
        CM.lines!(ax2, [lo_, hi_], [lo_, hi_]; color = :black)
        CM.lines!(ax2, [lo_, hi_], [lo_ + 1e-3, hi_ + 1e-3]; color = :black, linestyle = :dot)
        CM.lines!(ax2, [lo_, hi_], [max(lo_ - 1e-3, 1e-4), hi_ - 1e-3]; color = :black, linestyle = :dot)
    end
    if haskey(DATA, :eq79)
        ax3 = CM.Axis(fig2[1, 3]; xscale = log10, yscale = log10, xlabel = "T / (2 pi)", ylabel = "|nu_est| (rad per time), F_1 of Eq. (79), Double64", title = "R2: E1 on the Laskar2003 test function F_1, window order p = 0..5")
        for p in 0:5
            sel = sort([r for r in DATA[:eq79] if r.kind === :F1 && r.p == p], by = r -> r.x)
            CM.lines!(ax3, [10.0^r.x for r in sel], [max(r.ed, 1e-30) for r in sel]; label = "p = $p (theory 2p+2 = $(2p + 2))")
        end
        CM.axislegend(ax3; position = :lb, labelsize = 10)
    end
    CM.save(joinpath(OUT, "naff-fma-oracle.png"), fig2)
    say("# figure saved: naff-fma-oracle.png")
end

say(@sprintf("# verdicts: %d PASS, %d FAIL, %d other", count(==("PASS"), VERDICTS), count(==("FAIL"), VERDICTS), count(v -> v != "PASS" && v != "FAIL", VERDICTS)))
say(@sprintf("# total wall time %.1f s", time() - T0))
close(LOG)
