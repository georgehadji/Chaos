# Stage 04 verification of the analysis module (decision D60). Run from the repo root:
#   julia --project=verify -t 23 verify/analysis.jl
# Spec: workspaces/module-build/stages/01-spec/output/analysis-spec.md (acceptance criteria 3, 16, 17, 18, 19, 21, 22; oracles O1 and O3; parameter table, open point P1).
# Sections (ANV_ONLY="RUN,K0" runs only those): TRANSCRIPTION (parameter file against the preregistration), RUN (analyze on the real store with figures: time,
# memory, outputs, read-only, determinism between 1 and N threads), K0 (closed form on the real store), COMPARE (independent recomputation of
# verify/analysis-independent.jl, run beforehand: ANV_IND), CLI (the command, exit code 0).
# One printed line per check:  <id> | quantity | value | bound | PASS/FAIL/INFO
# The real run uses a manifest that is on branch study/p1-le-benchmark-run (git show e35ee5a:...): ANV_MANIFEST. Nothing in the study folders is written.
using Printf, SHA, TOML, Dates
using Chaoslab, CUDA
import CairoMakie as CM

const AN = Chaoslab.Analysis
const ROOT = dirname(@__DIR__)
const OUT = joinpath(ROOT, "workspaces", "module-build", "stages", "04-verify", "output")
const PREREG = joinpath(ROOT, "workspaces", "study", "stages", "01-design", "output", "p1-le-benchmark-preregistration.md")
const PANEL = joinpath(ROOT, "workspaces", "study", "stages", "02-panel", "output", "p1-le-benchmark-panel.toml")
const PARAMS = joinpath(ROOT, "workspaces", "study", "stages", "04-analyze", "references", "p1-le-benchmark-analysis-parameters.toml")
const SCR = raw"C:\Users\tesse\AppData\Local\Temp\claude\C--Users-tesse-Documents-Vibe-Coding-Chaos\5df15ee1-47df-4875-85bd-1bd4d7c665a6\scratchpad\real"
const MANIFEST = get(ENV, "ANV_MANIFEST", joinpath(SCR, "p1-le-benchmark-manifest.json"))
const STORE = get(ENV, "ANV_STORE", joinpath(ROOT, "data", "p1-le-benchmark"))
const IND = get(ENV, "ANV_IND", joinpath(SCR, "independent-results.txt"))
const ONLY = filter(!isempty, split(get(ENV, "ANV_ONLY", ""), ","))
const NT = Threads.nthreads()
mkpath(OUT); mkpath(SCR)
const LOG = open(get(ENV, "ANV_LOG", joinpath(OUT, "analysis-verify-log.txt")), "w")

fmt(x::AbstractFloat) = @sprintf("%.6g", x)
fmt(x) = string(x)
const COUNTS = Dict("PASS" => 0, "FAIL" => 0, "INFO" => 0)
function say(args...)
    s = string(args...)
    println(s); println(LOG, s); flush(stdout); flush(LOG)
end
function emit(id, quantity, value, bound, verdict)
    v = startswith(verdict, "FAIL") ? "FAIL" : startswith(verdict, "PASS") ? "PASS" : "INFO"
    COUNTS[v] += 1
    say(string(id, " | ", quantity, " | ", fmt(value), " | ", bound, " | ", verdict))
end
pf(ok::Bool) = ok ? "PASS" : "FAIL"
want(sec) = isempty(ONLY) || sec in ONLY
function section(f, name)
    want(name) || return
    say("== ", name, " ==")
    try
        f()
    catch e
        emit(name, "ERROR", sprint(showerror, e, catch_backtrace())[1:min(end, 900)], "-", "FAIL")
    end
end

const RES = Ref{Any}(nothing)
const OUT_RUN = joinpath(SCR, "analysis-out-run")
fingerprint() = Dict(f => bytes2hex(sha256(read(f))) for f in vcat([MANIFEST, PREREG, PANEL, PARAMS], sort(filter(f -> occursin(r"chunk-\d{5}\.h5$", f), readdir(STORE; join = true)))))

# ---------------------------------------------------------------------------------------------------------------- TRANSCRIPTION (P1)
const TOKENS = [  # key => (preregistration section, literal text that states the value there)
    :N_1 => ("Metrics", "N_1 = 1e7"), :N_long => ("Metrics", "N_long = 1e8"), :N_B => ("Metrics", "N_B = 1e6"), :n_H2 => ("H2", "n_H2 = 1e4"),
    :h1a_bound => ("H1", "<= 1e-8"), :h1a_fraction => ("H1", ">= 0.98"), :h1b_bound => ("H1", "<= 1e-6"), :h1b_fraction => ("H1", ">= 0.98"),
    :h1c_bound => ("H1", "<= 5e-3"), :h1c_fraction => ("H1", ">= 0.90"), :h1c_eps => ("H1", "eps_0 = 1e-8"),
    :control_bound => ("Amendment A4", "control_bound = 1e-8"), :control_bound_L9 => ("Amendment A6", "L9 = 5e-8"),
    :h2_median => ("H2", ">= 0.10"), :h2_interval => ("H2", ">= 0.05"), :bootstrap_resamples => ("H2", "10^4 resamples"), :bootstrap_level => ("H2", "95 percent"),
    :h3_band => ("H3", "0.05 lambda_ref"), :h3_ratio => ("H3", ">= 10"), :h3_interval => ("H3", "> 3"),
    :h4_ref_halfpi => ("H4", "0.298"), :h4_tol_halfpi => ("H4", "0.01"), :h4_ref_2pi => ("H4", "1.172"), :h4_tol_2pi => ("H4", "0.02"),
    :h4_T => ("H4", "T = 1e9"), :h4_T_cross => ("H4", "T = 1e8"), :h5a_min => ("H5", ">= 0.90"), :h5b_max => ("H5", "<= 0.30"), :h5c_max => ("H5", "<= 0.30"),
    :lambda_threshold => ("Label rules", "1e-4"), :growth_ratio => ("Label rules", "r(N) >= 5"), :rotation_tol => ("Label rules", "1e-7"),
    :return_tol => ("Label rules", "1e-10"), :ladder_tol => ("Label rules", "delta_lad = 5e-3"), :drift_bound => ("Eligible sets", "1e6"),
    :min_set => ("Analysis plan", "fewer than 20"), :K_c => ("K values", "0.971635406")]

function prereg_section(text, name)
    m = match(Regex("(?m)^#{2,3} " * replace(name, r"([\\^$.|?*+()\[\]{}])" => s"\\\1") * "\\b.*\\n"), text)
    m === nothing && return ""
    rest = text[m.offset + length(m.match):end]
    n = match(r"(?m)^#{2,3} ", rest)
    return n === nothing ? rest : rest[1:n.offset - 1]
end

section("TRANSCRIPTION") do
    text = replace(read(PREREG, String), "\r\n" => "\n")
    p = TOML.parsefile(PARAMS)
    mod = AN.read_params(PARAMS)
    emit("P1", "parameter file read by the module: keys, unknown or missing", length(mod), "37 keys", pf(length(mod) == 37 && length(p) == 37))
    bad = String[]
    for (k, (sec, tok)) in TOKENS
        body = prereg_section(text, sec)
        occursin(tok, body) || push!(bad, "$k: \"$tok\" not in section $sec")
        haskey(p, string(k)) || push!(bad, "$k: not in the parameter file")
    end
    emit("P1", "every value has its literal text in the preregistration section named for it (list of misses)", isempty(bad) ? "none" : join(bad, "; "), "none", pf(isempty(bad) && length(TOKENS) == 37))
    num = Dict(:N_1 => 1e7, :N_long => 1e8, :N_B => 1e6, :n_H2 => 1e4, :h1a_bound => 1e-8, :h1b_bound => 1e-6, :h1c_bound => 5e-3, :h1c_eps => 1e-8,
               :control_bound => 1e-8, :control_bound_L9 => 5e-8, :h4_T => 1e9, :h4_T_cross => 1e8, :K_c => 0.971635406, :h3_band => 0.05, :growth_ratio => 5.0)
    miss = [string(k) for (k, v) in num if Float64(getfield(mod, k)) != v]
    emit("P1", "numeric values equal the numbers of the cited text (mismatches)", isempty(miss) ? "none" : join(miss, ", "), "none", pf(isempty(miss)))
end

# ---------------------------------------------------------------------------------------------------------------- RUN (criteria 16, 17, 18, 21)
section("RUN") do
    rm(OUT_RUN; recursive = true, force = true)
    before = fingerprint()
    emit("RUN", "inputs fingerprinted before the run: manifest, panel, parameter file, preregistration, chunk files (count)", length(before), "1363 chunk files + 4", pf(length(before) == 1367))
    t = @timed AN.analyze(; prereg = PREREG, manifest = MANIFEST, panel = PANEL, params = PARAMS, store = STORE, out = OUT_RUN, figures = true, threads = NT)
    res = t.value
    RES[] = res
    rss = Sys.maxrss() / 2^30
    emit("O21", "wall time of analyze on the real store with figures (s), $NT threads, first call included", t.time, "<= 600 s", pf(t.time <= 600))
    emit("O21", "peak resident memory of the process after the run (GiB)", rss, "<= 8 GiB", pf(rss <= 8))
    emit("O21", "allocated by analyze (GiB), garbage-collector time (s)", @sprintf("%.1f, %.1f", t.bytes / 2^30, t.gctime), "-", "INFO")
    emit("RUN", "status", res.status, ":ok", pf(res.status === :ok))
    emit("O17", "inputs unchanged after the run (files compared by SHA-256)", fingerprint() == before, "true", pf(fingerprint() == before))
    pdfs = [joinpath(OUT_RUN, "p1-le-benchmark-fig-F$i.pdf") for i in 1:7]
    ok = all(f -> isfile(f) && filesize(f) > 1000 && String(read(f, 5)) == "%PDF-", pdfs)
    emit("O18", "seven figure files, each non-empty with the PDF signature (sizes kB)", join([isfile(f) ? round(filesize(f) / 1024; digits = 1) : 0 for f in pdfs], " "), "7 files", pf(ok))
    text = read(res.tables, String)
    dg = [match(Regex("## F" * string(i) * " data[^\\n]*\\n(?:.*\\n)*?digest: ([0-9a-f]{64})"), text).captures[1] for i in 1:7]
    emit("O18", "digest lines of the tables file equal the digests of the figure data", dg == [res.digests["F$i"] for i in 1:7], "true", pf(dg == [res.digests["F$i"] for i in 1:7]))
    OUT1 = joinpath(SCR, "analysis-out-1thread")
    rm(OUT1; recursive = true, force = true)
    t1 = @timed AN.analyze(; prereg = PREREG, manifest = MANIFEST, panel = PANEL, params = PARAMS, store = STORE, out = OUT1, figures = false, threads = 1)
    same = read(res.tables) == read(t1.value.tables)
    emit("O16", "tables file at $NT threads against 1 thread: byte-identical (1-thread run $(round(t1.time; digits = 1)) s)", same, "true", pf(same && t1.value.digests == res.digests))
    emit("O17", "inputs unchanged after the second run", fingerprint() == before, "true", pf(fingerprint() == before))
end

# ---------------------------------------------------------------------------------------------------------------- K0 (criterion 3, O1)
section("K0") do
    items = AN.read_store(MANIFEST, STORE)
    panel = TOML.parsefile(PANEL)["orbit"]
    ts = [Int(o["id"]) for o in panel if !o["control"] && o["class"] == "T" && o["K"] == 0.0]
    emit("O1", "class T orbits at K = 0 in the panel", length(ts), "40", pf(length(ts) == 40))
    a = pi / 4
    closed(n) = log(sqrt((cos(a) + n * sin(a))^2 + sin(a)^2)) / n          # deviation-vector growth (preregistration A1 item 2): L1 and the QR family
    closed_svd(n) = log((n + sqrt(n^2 + 4)) / 2) / n                       # L9 stores ln(largest singular value of the product [[1, n], [0, 1]]) / n
    idx = AN.index_items(items)
    N1 = 10_000_000
    steps_of(idx, id) = idx[AN.item_key(id, "L5", "Float64", N1, 1, nothing, "pi/4", "lifted")].t
    worst = Dict{String, Float64}()
    curve = Dict{String, Vector{Float64}}()
    for m in AN.V_METHODS
        dev = Float64[]
        for id in ts
            it = idx[AN.item_key(id, m, "Float64", N1, 1, nothing, "pi/4", "lifted")]
            d = abs.(it.lam .- (m == "L9" ? closed_svd : closed).(Float64.(it.t)))
            isempty(dev) ? append!(dev, d) : (dev .= max.(dev, d))
        end
        curve[m] = dev
        worst[m] = maximum(dev)
    end
    emit("O1", "K = 0: largest |lambda_m(n) - closed form| over 40 orbits, every recorded n, per method (L9 against the singular-value closed form)", join(["$m $(fmt(worst[m]))" for m in AN.V_METHODS], ", "), "<= 1e-9", pf(maximum(values(worst)) <= 1e-9))
    l9_wrong = maximum(maximum(abs.(idx[AN.item_key(id, "L9", "Float64", N1, 1, nothing, "pi/4", "lifted")].lam .- closed.(Float64.(steps_of(idx, id))))) for id in ts)
    emit("O1", "finding: L9 against the deviation-vector closed form (spec O1 as written), largest over 40 orbits", l9_wrong, "-", "INFO")
    steps = idx[AN.item_key(ts[1], "L5", "Float64", N1, 1, nothing, "pi/4", "lifted")].t
    fig = CM.Figure(size = (720, 460))
    ax = CM.Axis(fig[1, 1]; xscale = log10, yscale = log10, xlabel = "n", ylabel = "max over 40 orbits |lambda_m(n) - closed form|", title = "K = 0: stored exponents against ln sqrt((cos a + n sin a)^2 + sin^2 a) / n, a = pi/4")
    for m in AN.V_METHODS
        CM.lines!(ax, steps, max.(curve[m], 1e-20); label = m)
    end
    CM.hlines!(ax, [1e-9]; color = :black, linestyle = :dash)
    CM.axislegend(ax; position = :rt, nbanks = 2)
    CM.save(joinpath(OUT, "analysis-k0-closed-form.png"), fig)
end

# ---------------------------------------------------------------------------------------------------------------- COMPARE (criterion 19, O3)
section("COMPARE") do
    res = RES[] === nothing ? AN.analyze(; prereg = PREREG, manifest = MANIFEST, panel = PANEL, params = PARAMS, store = STORE, out = joinpath(SCR, "analysis-out-compare"), threads = NT) : RES[]
    isfile(IND) || error("independent results file $IND not found")
    ind = Dict{String, Vector{Vector{String}}}()
    for l in eachline(IND)
        isempty(l) && continue
        f = split(l, '|')
        push!(get!(ind, String(f[1]), Vector{String}[]), String.(f[2:end]))
    end
    toint(s) = Int(round(parse(Float64, s)))
    num(s) = s == "NaN" ? NaN : parse(Float64, s)
    near(a, b; rtol = 1e-12) = (isnan(a) && isnan(b)) || a == b || abs(a - b) <= rtol * max(abs(a), abs(b), 1e-300)
    listids(s) = s == "-" ? Int[] : sort!([toint(x) for x in split(s, ',')])
    byid = Dict(r.id => r for r in res.records)
    rows = ind["ORBIT"]
    emit("O3", "orbits in the independent file and in the module", "$(length(rows)) and $(length(res.records))", "equal", pf(length(rows) == length(res.records)))
    mism = Dict(k => 0 for k in ("lam_ref", "delta_qb", "diff_l9", "diff_l2", "R_nH2", "nconv", "defined", "horizon", "raw label", "final label", "forced", "ladder_same", "flags"))
    for r in rows
        id = toint(r[1])
        m = byid[id]
        # fields after the tag: 1 id, 2 class, 3 K, 4 N, 5 lam_ref, 6 delta_qb, 7 diff_l9, 8 diff_l2, 9 R_nH2, 10 nconv, 11 defined, 12-13 raw label, 14-15 label, 16 forced, 17 ladder_same, 18 flags
        near(num(r[5]), m.lam_ref) || (mism["lam_ref"] += 1)
        near(num(r[6]), m.h1.delta_qb) || (mism["delta_qb"] += 1)
        near(num(r[7]), m.h1.diff_l9) || (mism["diff_l9"] += 1)
        near(num(r[8]), m.h1.diff_l2) || (mism["diff_l2"] += 1)
        near(num(r[9]), m.R_nH2) || (mism["R_nH2"] += 1)
        toint(r[10]) == m.conv.n || (mism["nconv"] += 1)
        (r[11] == "true") == m.conv.defined || (mism["defined"] += 1)
        toint(r[4]) == m.horizon || (mism["horizon"] += 1)
        (r[12], r[13]) == (m.label.raw.instability, m.label.raw.structure) || (mism["raw label"] += 1)
        (r[14], r[15]) == (m.label.instability, m.label.structure) || (mism["final label"] += 1)
        Set(r[16] == "-" ? String[] : String.(split(r[16], ','))) == Set(string.(m.label.forced)) || (mism["forced"] += 1)
        r[17] == (m.label.ladder_same === nothing ? "-" : string(m.label.ladder_same)) || (mism["ladder_same"] += 1)
        Set(r[18] == "-" ? String[] : String.(split(r[18], ','))) == Set(m.flags) || (mism["flags"] += 1)
    end
    for k in sort!(collect(keys(mism)))
        emit("O3", "per-orbit quantity $k: orbits that differ from the independent value (of $(length(rows)))", mism[k], "0", pf(mism[k] == 0))
    end
    S = Dict(r[1] => r[2] for r in ind["SET"])
    sets = res.sets
    same_sets = listids(S["E_ch"]) == sort(sets.E_ch) && listids(S["W"]) == sort(sets.W) && listids(S["S"]) == sort(sets.S)
    emit("O3", "E_ch, W, S: independent ids equal the module's (sizes)", "$(length(sets.E_ch)) $(length(sets.W)) $(length(sets.S))", "equal", pf(same_sets))
    dr = Set(S["dropped"] == "-" ? String[] : String.(split(S["dropped"], ',')))
    emit("O3", "dropped orbits and S orbits without positive growth", "$(length(dr)), $(sets.S_nonpositive)", "equal", pf(dr == Set("$(d.id):$(d.reason)" for d in sets.dropped) && toint(S["S_notpositive_label"]) == sets.S_nonpositive))
    d = res.decisions
    h1 = ind["H1"][1]
    emit("O3", "H1: eligible orbits and passes a, b, c", "$(d.H1.n) $(d.H1.passes.a) $(d.H1.passes.b) $(d.H1.passes.c)", "equal", pf([toint(x) for x in h1[1:4]] == [d.H1.n, d.H1.passes.a, d.H1.passes.b, d.H1.passes.c] && h1[5] == string(d.H1.verdict)))
    wid2 = d.H2.hi - d.H2.lo
    h2 = ind["H2"][1]
    emit("O3", "H2: n, median (exact), interval ends within 0.05 of the width, verdict", "$(d.H2.n) $(fmt(d.H2.median)) [$(fmt(d.H2.lo)), $(fmt(d.H2.hi))]", "equal / 0.05 w",
         pf(toint(h2[1]) == d.H2.n && near(num(h2[2]), d.H2.median) && abs(num(h2[3]) - d.H2.lo) <= 0.05 * wid2 + 1e-15 && abs(num(h2[4]) - d.H2.hi) <= 0.05 * wid2 + 1e-15 && h2[5] == string(d.H2.verdict)))
    wid3 = d.H3.hi - d.H3.lo
    h3 = ind["H3"][1]
    mw, ms = AN.median_of([byid[i].conv.n for i in sets.W]), AN.median_of([byid[i].conv.n for i in sets.S])
    emit("O3", "H3: sizes, medians, ratio (exact), interval ends within 0.05 of the width, verdict", "$(d.H3.n_w) $(d.H3.n_s) $(fmt(mw)) $(fmt(ms)) $(fmt(d.H3.ratio)) [$(fmt(d.H3.lo)), $(fmt(d.H3.hi))]", "equal / 0.05 w",
         pf([toint(h3[1]), toint(h3[2])] == [d.H3.n_w, d.H3.n_s] && near(num(h3[3]), mw) && near(num(h3[4]), ms) && near(num(h3[5]), d.H3.ratio) &&
            abs(num(h3[6]) - d.H3.lo) <= 0.05 * wid3 + 1e-15 && abs(num(h3[7]) - d.H3.hi) <= 0.05 * wid3 + 1e-15 && h3[8] == string(d.H3.verdict)))
    h4 = ind["H4"][1]
    mine = [d.H4.halfpi.L1, d.H4.halfpi.L5, d.H4.twopi.L1, d.H4.twopi.L5, d.H4.crosscheck[1].float64, d.H4.crosscheck[1].double64, d.H4.crosscheck[2].float64, d.H4.crosscheck[2].double64]
    emit("O3", "H4: eight values equal, verdict", join(fmt.(mine), " "), "equal", pf(all(near(num(h4[i]), mine[i]) for i in 1:8) && h4[9] == string(d.H4.verdict)))
    h5 = ind["H5"][1]
    k(s) = round(Int, s.fraction * s.n)
    emit("O3", "H5: sizes and counts of the three sub-rules, verdict", "$(d.H5.h5a.n)/$(k(d.H5.h5a)) $(d.H5.h5b.n)/$(k(d.H5.h5b)) $(d.H5.h5c.n)/$(k(d.H5.h5c))", "equal",
         pf([toint(x) for x in h5[1:6]] == [d.H5.h5a.n, k(d.H5.h5a), d.H5.h5b.n, k(d.H5.h5b), d.H5.h5c.n, k(d.H5.h5c)] && h5[7] == string(d.H5.verdict)))
    mc = Dict((c.k_index, c.class, l[1], l[2]) => n for c in res.cells for (l, n) in c.labels)
    ic = Dict((toint(r[1]), r[2], r[3], r[4]) => toint(r[5]) for r in ind["CELL"])
    emit("O3", "label counts per cell and label pair: entries / differing entries", "$(length(mc)) / $(count(k -> get(ic, k, -1) != get(mc, k, -2), union(keys(mc), keys(ic))))", "0 differ", pf(mc == ic))
    kx = get(ind, "K0X", Vector{String}[])
    emit("O3", "independent K = 0 closed forms: largest deviation of the seven deviation-vector methods, and of L9 against the singular-value form, over $(length(kx)) orbits",
         isempty(kx) ? "none" : "$(fmt(maximum(num(r[2]) for r in kx))), $(fmt(maximum(num(r[3]) for r in kx)))", "<= 1e-9", pf(!isempty(kx) && maximum(num(r[2]) for r in kx) <= 1e-9 && maximum(num(r[3]) for r in kx) <= 1e-9))
end

# ---------------------------------------------------------------------------------------------------------------- CLI (criterion 22)
section("CLI") do
    out = joinpath(SCR, "analysis-out-cli")
    rm(out; recursive = true, force = true)
    jl = Base.julia_cmd()
    cmd = `$jl --project=$(joinpath(ROOT, "verify")) -t $NT $(joinpath(ROOT, "bin", "chaoslab.jl")) analyze $PREREG $MANIFEST $PANEL $PARAMS --store $STORE --out $out`
    p = run(ignorestatus(cmd))
    emit("O22", "bin/chaoslab.jl analyze on the real inputs: exit code", p.exitcode, "0", pf(p.exitcode == 0))
    n = length(filter(f -> endswith(f, ".pdf") || endswith(f, "-tables.md"), isdir(out) ? readdir(out) : String[]))
    emit("O22", "tables file and seven figure files written", n, "8", pf(n == 8))
    bad = run(ignorestatus(`$jl --project=$(joinpath(ROOT, "verify")) $(joinpath(ROOT, "bin", "chaoslab.jl")) analyze $PREREG --bogus x`))
    emit("O22", "unknown option and too few arguments: exit code", bad.exitcode, "1", pf(bad.exitcode == 1))
end

say(@sprintf("TOTAL | PASS %d | FAIL %d | INFO %d | %s", COUNTS["PASS"], COUNTS["FAIL"], COUNTS["INFO"], Dates.format(now(), "yyyy-mm-dd HH:MM")))
close(LOG)
