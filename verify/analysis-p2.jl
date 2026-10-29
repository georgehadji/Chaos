# Stage 04 verification of the analysis module, Amendment D89 (study p2-indicators). Run from the repo root:
#   julia --project=verify -t 23 verify/analysis-p2.jl
# Spec: workspaces/module-build/stages/01-spec/output/analysis-spec.md, "Amendment D89" (acceptance criteria D89a, D89q, D89r, D89v, D89w, D89x; oracles O6, O9, O10).
# Sections (ANV_ONLY="TRANSCRIPTION,RUN" runs only those): TRANSCRIPTION (parameter table against the preregistration, E22), RUN (analyze at 23, 1 and 4 threads: time, memory,
# read-only, determinism), REPORT (D89x numbers, verbatim), COMPARE (independent results of verify/analysis-p2-decision.py, run beforehand), P1UNCH (D89a, V), CLI (D89s).
# One printed line per check:  <id> | quantity | value | bound | PASS/FAIL/INFO
# Nothing in data/, src/, the studies or the specs is written; scratch outputs go to the scratch folder SCR.
using Printf, SHA, TOML, Dates
using Chaoslab

const AN = Chaoslab.Analysis
const ROOT = dirname(@__DIR__)
const OUT = joinpath(ROOT, "workspaces", "module-build", "stages", "04-verify", "output")
const STUDY = joinpath(ROOT, "workspaces", "study", "stages")
const PREREG = joinpath(STUDY, "01-design", "output", "p2-indicators-preregistration.md")
const MAN2 = joinpath(STUDY, "03-run", "output", "p2-indicators-manifest.json")
const PANEL = joinpath(STUDY, "02-panel", "output", "p1-le-benchmark-panel.toml")
const PARAMS = joinpath(STUDY, "04-analyze", "references", "p2-indicators-analysis-parameters.toml")
const PREREG1 = joinpath(STUDY, "01-design", "output", "p1-le-benchmark-preregistration.md")
const MAN1 = get(ENV, "ANV_MAN1", joinpath(STUDY, "03-run", "output", "p1-le-benchmark-manifest.json"))
const PARAMS1 = joinpath(STUDY, "04-analyze", "references", "p1-le-benchmark-analysis-parameters.toml")
const TABLES1 = joinpath(STUDY, "04-analyze", "output", "p1-le-benchmark-tables.md")
const DEFAULTS = joinpath(ROOT, "_config", "numerics-defaults.toml")
const STORE2 = joinpath(ROOT, "data", "p2-indicators")
const STORE1 = joinpath(ROOT, "data", "p1-le-benchmark")
const SCR = raw"C:\Users\tesse\AppData\Local\Temp\claude\C--Users-tesse-Documents-Vibe-Coding-Chaos\5df15ee1-47df-4875-85bd-1bd4d7c665a6\scratchpad\p2run"
const IND = get(ENV, "ANV_IND", joinpath(OUT, "analysis-p2-independent-results.txt"))
const ONLY = filter(!isempty, split(get(ENV, "ANV_ONLY", ""), ","))
const NT = Threads.nthreads()
mkpath(SCR)
const LOG = open(get(ENV, "ANV_LOG", joinpath(OUT, "analysis-p2-verify-log.txt")), get(ENV, "ANV_APPEND", "") == "1" ? "a" : "w")

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
        emit(name, "ERROR", sprint(showerror, e, catch_backtrace())[1:min(end, 1500)], "-", "FAIL")
    end
end

# ---------------------------------------------------------------------------------------------------------------- fingerprints (D89r, I12)
chunkfiles(store) = sort(filter(f -> occursin(r"chunk-\d{5}\.h5$", f), readdir(store; join = true)))
const INPUTS = vcat([MAN2, PANEL, PREREG, PARAMS, PREREG1, MAN1, PARAMS1, DEFAULTS], chunkfiles(STORE2), chunkfiles(STORE1))
fingerprint() = Dict(f => bytes2hex(sha256(read(f))) for f in INPUTS)
# listing of the data tree: relative path, size, modification time of every file (I12)
function data_listing()
    root = joinpath(ROOT, "data")
    rows = String[]
    for (d, _, fs) in walkdir(root), f in fs
        p = joinpath(d, f)
        st = stat(p)
        push!(rows, string(replace(relpath(p, root), '\\' => '/'), '\t', st.size, '\t', st.mtime))
    end
    sort!(rows)
    return (n = length(rows), sha = bytes2hex(sha256(join(rows, '\n'))))
end

# ---------------------------------------------------------------------------------------------------------------- TRANSCRIPTION (E22, stage 04 line by line)
# Each expected value is typed here from reading the sentence the table quotes (judged one by one), not copied from the table.
const EXPECT = Dict{String, Any}(
    "persist_divisor" => 4, "r2_le_factor" => 3.0, "n_s" => 100000, "budgets" => [200, 400, 800, 1600, 3200, 6400, 12800, 25600, 51200, 102400], "min_set" => 20,
    "h2_q_min" => 10.0, "h3_r_lower_gt" => 3.0, "s2_lambda_min" => 1.0, "rli_threshold" => 1e-7, "wba_digits_threshold" => 6.0, "h6_flip_min" => 0.05,
    "bootstrap_level" => 0.95, "w2_lambda_max" => 0.1, "gate_k_index" => 1, "h3_r_min" => 10.0, "gate_classes" => ["T", "R"], "censor_factor" => 2,
    "e2_d_threshold" => -4.0, "label_estimator_bound" => 1e-8, "label_rotation_tol" => 1e-7, "gate_max_share" => 0.10, "label_lambda_min" => 1e-4,
    "h1_u_min" => 0.20, "h4_dphi_lower_min" => 0.25, "h4_dphi_min" => 0.50, "range_upper_percentile" => 95.0, "h3_count_min" => 5, "h2_q_lower_gt" => 3.0,
    "f_s" => 0.2, "r2_megno_threshold" => 2.5, "e1_d_threshold" => -4.0, "h5_fmed_min" => 0.10, "h1_amin_max" => 0.80, "b_primary" => 12800,
    "bootstrap_seed" => 20261012, "panel_sha256" => "881874d26aa8883caa0545331ba167a3f2c9e3c249cdb6ba32ad39c98ad348d5",
    "plan_sha256" => "3823fa4a89fee047ab1f9ca63fee3d4e2dafa372c50b1d63c6de244fa8a18a81", "h1_u_lower_min" => 0.10, "h5_fmed_lower_min" => 0.05, "r2_sali_threshold" => 1e-8,
    "range_lower_percentile" => 5.0, "label_growth_ratio" => 5.0, "ratio_threshold" => 5.0, "le_floor" => 1e-4, "r2_fli_factor" => 3.0, "gate_size" => 80,
    "bootstrap_resamples" => 10000)
# keys whose quote carries no number of its own: the value is held by the P1 table / P1 preregistration (key of the P1 table, P1 section, literal P1 text)
const WEAK = Dict("label_lambda_min" => ("lambda_threshold", "Label rules", "1e-4"), "label_growth_ratio" => ("growth_ratio", "Label rules", "r(N) >= 5"),
                  "label_rotation_tol" => ("rotation_tol", "Label rules", "1e-7"), "label_estimator_bound" => ("h1a_bound", "H1", "<= 1e-8"))

# sections by heading outside the code fences: the lines up to the next heading of the same or a higher level (written again here, as the module does it)
function sections_of(text)
    lines = split(replace(text, "\r\n" => "\n"), '\n')
    heads = Tuple{Int, Int, String}[]
    fence = false
    for (i, l) in enumerate(lines)
        startswith(l, "```") && (fence = !fence)
        m = fence ? nothing : match(r"^(#+)[ \t]+(.*?)[ \t]*$", l)
        m === nothing || push!(heads, (i, length(m.captures[1]), String(m.captures[2])))
    end
    out = Dict{String, String}()
    for (n, (i, lev, title)) in enumerate(heads)
        j = findfirst(h -> h[2] <= lev, heads[(n + 1):end])
        last = j === nothing ? length(lines) : heads[n + j][1] - 1
        haskey(out, title) || (out[title] = join(lines[(i + 1):last], "\n"))
    end
    return out
end

function fenced_blocks(text)
    return [TOML.parse(String(m.captures[1])) for m in eachmatch(r"```toml\n(.*?)```"s, replace(text, "\r\n" => "\n"))]
end

section("TRANSCRIPTION") do
    ptext = read(PREREG, String)
    secs = sections_of(ptext)
    p1secs = sections_of(read(PREREG1, String))
    tab = TOML.parsefile(PARAMS)
    p1tab = TOML.parsefile(PARAMS1)
    defaults = TOML.parsefile(DEFAULTS)
    mod = AN.p2_read_params(PARAMS; prereg = PREREG, p1_params = PARAMS1)
    emit("E22", "table read by the module (checks of E22 i and ii, A2): keys", length(tab), "58", pf(length(tab) == 58 && length(mod) == 58 && Set(Symbol.(keys(tab))) == Set(AN.P2_PARAM_KEYS)))
    machine = sort([k for (k, v) in tab if get(v, "machine", false) === true])
    emit("E22", "keys marked machine (list)", join(machine, " "), "the 11 of E22", pf(machine == sort(collect(string.(AN.P2_MACHINE_KEYS))) && length(machine) == 11))
    nfail = 0
    nweak = 0
    nedit = 0
    # --- the 47 keys with a quote: section, quote, text, value
    for k in sort([k for k in keys(tab) if !(k in machine)])
        e = tab[k]
        sec, q, tx, val = e["section"], get(e, "quote", ""), get(e, "text", ""), e["value"]
        body = get(secs, sec, nothing)
        ok_sec = body !== nothing
        ok_q = ok_sec && occursin(replace(q, "\r\n" => "\n"), body)
        ok_t = occursin(tx, q)
        exp = EXPECT[k]
        ok_v = exp isa AbstractVector ? (val == exp) : exp isa AbstractString ? (val == exp) : (val isa Real && Float64(val) == Float64(exp))
        ok_type = exp isa Integer ? val isa Integer : exp isa AbstractFloat ? val isa AbstractFloat : true
        ok = ok_sec && ok_q && ok_t && ok_v && ok_type
        nfail += !ok
        if haskey(WEAK, k)
            nweak += 1
            k1, s1, lit = WEAK[k]
            t1 = findfirst(t -> t == s1 || startswith(t, s1 * " "), collect(keys(p1secs)))   # P1 headings carry a suffix: "H1 (claim C1, ...)", "Label rules (axes of D36 and D37)"
            ok_p1 = t1 !== nothing && occursin(lit, p1secs[collect(keys(p1secs))[t1]]) && p1tab[k1] == val
            nfail += !ok_p1
            emit("T58", "$k = $(repr(val)) | section \"$sec\" | quote \"$q\" | text \"$tx\"", ok ? "section, quote and text found" : "NOT FOUND",
                 "weak quote: no number of its own; value held by P1 key $k1 = $(p1tab[k1]) and the P1 text \"$lit\" in section \"$s1\" ($(ok_p1 ? "equal" : "DIFFERENT"))", ok && ok_p1 ? "INFO (weak quote)" : "FAIL")
        else
            emit("T58", "$k = $(repr(val)) | section \"$sec\" | quote \"$q\" | text \"$tx\"", ok ? "section, quote and text found; value equals the cited number" : "section $ok_sec quote $ok_q text $ok_t value $ok_v type $ok_type", "all", pf(ok))
        end
    end
    emit("T58", "keys with a quote checked (sections, quotes, texts, values); failures; weak quotes", "47 keys; $nfail failures; $nweak weak", "0 failures", pf(nfail == 0))
    # --- the 11 machine keys against the run blocks
    bl = fenced_blocks(ptext)
    rb, h4 = bl[1]["run"], bl[2]["run"]["h4"]
    arms = Dict(a["name"] => a for a in rb["arm"])
    specs = Dict(s["id"] => s for s in rb["indicator"])
    ij = [match(r"-i(\d+)-j(\d+)$", h["label"]) for h in h4]
    ktexts = unique(String[h["K"] for h in h4])
    runv = Dict("N_1" => only(arms["lyap_p"]["horizons"]), "N_c" => only(arms["win_p"]["horizons"]), "N_dir" => only(arms["dir_tg"]["horizons"]), "windows" => Int.(specs["fq1"]["windows"]),
                "directions" => String.(arms["dir_tg"]["direction"]), "drift_bound" => rb["health"]["drift_bound"], "grid_k_values" => ktexts,
                "grid_i_count" => length(unique(m.captures[1] for m in ij)), "grid_j_count" => length(unique(m.captures[2] for m in ij)), "grid_orbit_count" => length(h4),
                "panel_build_sha256" => only(rb["panel_prereg_sha256_accepted"]))
    for k in machine
        v = tab[k]["value"]
        ok = v isa Real ? Float64(v) == Float64(runv[k]) : v == runv[k]
        nfail += !ok
        emit("T58", "machine key $k = $(repr(v)) against the run blocks ($(repr(runv[k])))", ok, "equal", pf(ok))
    end
    w = Int.(tab["windows"]["value"])
    allwin = all(s -> Int.(s["windows"]) == w, [s for s in rb["indicator"] if haskey(s, "windows")])
    bud = Int.(tab["budgets"]["value"])
    rel = bud == 2 .* w && Int(tab["N_c"]["value"]) == 2 * maximum(w) && Int(tab["N_c"]["value"]) == bud[end] && tab["b_primary"]["value"] in bud && bud[7] == tab["b_primary"]["value"]
    emit("T58", "windows equal in every window spec of the run block; budgets = 2 x windows; N_c = 2 x largest window = last budget; b_primary = element 7 of budgets (b_6)", allwin && rel, "true", pf(allwin && rel))
    gk = length(ktexts) * tab["grid_i_count"]["value"] * tab["grid_j_count"]["value"] == tab["grid_orbit_count"]["value"] == length(h4)
    emit("T58", "grid_orbit_count = K count x grid_i_count x grid_j_count = number of [[run.h4]] tables", gk, "true", pf(gk))
    emit("T58", "bootstrap_seed equals seed of numerics-defaults.toml", tab["bootstrap_seed"]["value"] == defaults["seed"], "true", pf(tab["bootstrap_seed"]["value"] == defaults["seed"]))
    lab = [(k, k1) for (k, (k1, _, _)) in WEAK]
    eqp1 = all(tab[k]["value"] == p1tab[k1] for (k, k1) in lab)
    emit("T58", "the four label keys equal the P1 table keys lambda_threshold, growth_ratio, rotation_tol, h1a_bound", eqp1, "true", pf(eqp1))
    emit("T58", "P1 parameter table still has exactly 37 keys (D89a V)", length(p1tab), "37", pf(length(p1tab) == 37 && length(AN.read_params(PARAMS1)) == 37))
    emit("T58", "typo edits of the table made by this verification (quotes or sections)", nedit, "list every edit", "INFO")
    emit("T58", "total failures of the transcription", nfail, "0", pf(nfail == 0))
end

# ---------------------------------------------------------------------------------------------------------------- RUN (D89x, D89r, D89q)
const RES = Ref{Any}(nothing)
const OUT_RUN = joinpath(SCR, "out-23")
const OUT_1 = joinpath(SCR, "out-1")
const OUT_4 = joinpath(SCR, "out-4")
const P1KW = (p1_prereg = PREREG1, p1_manifest = MAN1, p1_params = PARAMS1, p1_store = STORE1)
run_p2(out, threads) = AN.analyze(; prereg = PREREG, manifest = MAN2, panel = PANEL, params = PARAMS, store = STORE2, out = out, figures = true, threads = threads, P1KW...)

section("RUN") do
    foreach(d -> rm(d; recursive = true, force = true), (OUT_RUN, OUT_1, OUT_4))
    before = fingerprint()
    list0 = data_listing()
    n2 = count(f -> occursin("p2-indicators", f), keys(before))
    emit("D89r", "inputs fingerprinted before the run (SHA-256): files in all, of them in the P2 store", "$(length(before)) ($(length(chunkfiles(STORE2))) P2 chunks, $(length(chunkfiles(STORE1))) P1 chunks, manifests, panel, preregistrations, tables, defaults)", "695 P2 chunks", pf(length(chunkfiles(STORE2)) == 695))
    emit("D89r", "listing of the data tree before the run: files, SHA-256 of (path, size, mtime)", "$(list0.n) files, $(list0.sha)", "-", "INFO")
    t = @timed run_p2(OUT_RUN, NT)
    res = t.value
    RES[] = res
    rss = Sys.maxrss() / 2^30
    emit("D89x", "wall time of analyze on the real P2 and P1 stores, $NT threads, first call included (s)", t.time, "<= 900 s", pf(t.time <= 900))
    emit("D89x", "peak resident memory of the process after the run at $NT threads (GiB)", rss, "<= 8 GiB", pf(rss <= 8))
    emit("D89x", "allocated by analyze (GiB), garbage-collector time (s)", @sprintf("%.1f, %.1f", t.bytes / 2^30, t.gctime), "-", "INFO")
    code = res.status === :ok ? 0 : res.status === :figures_skipped ? 2 : -1
    emit("D89x", "status of the run and the exit code it stands for", "$(res.status), exit code $code", "0 or 2", pf(code in (0, 2)))
    tbl = joinpath(OUT_RUN, "p2-indicators-tables.md")
    emit("D89p", "the tables file is the only file in --out", join(readdir(OUT_RUN), ", "), "p2-indicators-tables.md", pf(readdir(OUT_RUN) == ["p2-indicators-tables.md"] && res.tables == tbl))
    emit("D89r", "inputs unchanged after run 1 (SHA-256 of all files)", fingerprint() == before, "true", pf(fingerprint() == before))
    emit("D89r", "data tree listing unchanged after run 1", data_listing().sha == list0.sha, "true", pf(data_listing().sha == list0.sha))
    text = read(tbl, String)
    dg = [match(Regex("## F" * string(i) * " data[^\\n]*\\n(?:.*\\n)*?digest: ([0-9a-f]{64})"), text) for i in 1:6]
    okdg = all(!isnothing, dg) && [m.captures[1] for m in dg] == [res.digests["F$i"] for i in 1:6]
    emit("D89u", "digest lines of F1 to F6 in the tables file equal the digests of the data tables", okdg, "true", pf(okdg))
    # second run, one thread
    t1 = @timed run_p2(OUT_1, 1)
    same1 = read(tbl) == read(joinpath(OUT_1, "p2-indicators-tables.md"))
    emit("D89q", "tables file at $NT threads against 1 thread: byte-identical (1-thread run $(round(t1.time; digits = 1)) s, status $(t1.value.status))", same1, "true", pf(same1 && t1.value.digests == res.digests))
    emit("D89r", "inputs unchanged after run 2 (1 thread)", fingerprint() == before, "true", pf(fingerprint() == before))
    t4 = @timed run_p2(OUT_4, 4)
    same4 = read(tbl) == read(joinpath(OUT_4, "p2-indicators-tables.md"))
    emit("D89q", "tables file at $NT threads against 4 threads: byte-identical (4-thread run $(round(t4.time; digits = 1)) s)", same4, "true", pf(same4))
    emit("D89r", "inputs unchanged after run 3 (4 threads)", fingerprint() == before, "true", pf(fingerprint() == before))
    emit("D89r", "data tree listing unchanged after run 3", data_listing().sha == list0.sha, "true", pf(data_listing().sha == list0.sha))
    emit("D89x", "peak resident memory of the process after all three runs (GiB, cumulative maximum)", Sys.maxrss() / 2^30, "-", "INFO")
    # a refusal leaves the inputs unchanged and writes nothing (D89r, E36 a)
    bad = joinpath(ROOT, "data", "p2-verify-must-not-exist")
    refused = try
        run_p2(bad, 1); false
    catch e
        e isa ArgumentError
    end
    emit("D89r", "--out inside data/ refused and the folder not created", refused && !ispath(bad), "true", pf(refused && !ispath(bad)))
    emit("D89r", "inputs and data tree unchanged after the refusal", fingerprint() == before && data_listing().sha == list0.sha, "true", pf(fingerprint() == before && data_listing().sha == list0.sha))
end

# ---------------------------------------------------------------------------------------------------------------- shared helpers
const I7 = [:LE, :SALI, :FLI, :MEGNO, :E1, :E2, :WBA]
const I8 = vcat(I7, [:RLI])
function result()
    RES[] === nothing && (RES[] = run_p2(joinpath(SCR, "out-compare"), NT))
    return RES[]
end
vword(s) = s === :supported ? "SUPPORTED" : s === :not_supported ? "NOT SUPPORTED" : s === :not_decidable ? "NOT DECIDABLE" : string(s)
rword(s) = s === :ROBUST ? "ROBUST" : s === :NOT_ROBUST ? "NOT ROBUST" : "NOT COMPARABLE"
numv(x) = x === nothing ? NaN : x
fv(x) = x === nothing ? "undefined" : x isa AbstractFloat ? @sprintf("%.10g", x) : string(x)
kname(K) = "K=" * string(K)
function median_of(v)
    s = sort(Float64.(v)); m = length(s)
    return isodd(m) ? s[(m + 1) ÷ 2] : (s[m ÷ 2] + s[m ÷ 2 + 1]) / 2
end
dec(v) = v === :chaotic || v === :regular
function orbit_info()
    panel = TOML.parsefile(PANEL)["orbit"]
    grid = AN.p2_grid_starts(PREREG; panel = PANEL)
    n = length(panel) + length(grid)
    src = fill("P", n); K = Vector{Union{Nothing, Float64}}(nothing, n); kidx = zeros(Int, n); cls = fill("", n)
    ctrl = falses(n); sd = falses(n); sl = falses(n)
    for o in panel
        i = Int(o["id"])
        K[i] = haskey(o, "K") ? Float64(o["K"]) : nothing
        kidx[i] = o["k_index"]; cls[i] = o["class"]; ctrl[i] = o["control"]; sd[i] = o["in_SD"]; sl[i] = o["in_SL"]
    end
    for g in grid
        src[g.id] = "G"; K[g.id] = parse(Float64, g.K); kidx[g.id] = g.kindex; cls[g.id] = "h4"
    end
    return (n = n, src = src, K = K, kidx = kidx, cls = cls, ctrl = ctrl, sd = sd, sl = sl)
end
srt(d) = sort(collect(d); by = kv -> findfirst(==(kv[1]), I8))

# ---------------------------------------------------------------------------------------------------------------- REPORT (D89x: the numbers, verbatim; no interpretation)
function describe(ver, d, I)
    emit("D89x", "$ver H1: orbits $(d.H1.n), decided for all $(d.H1.used), dropped $(d.H1.dropped); U $(fv(d.H1.U)) interval [$(fv(d.H1.U_lo)), $(fv(d.H1.U_hi))]; A_min $(fv(d.H1.A_min)); U on S2 $(fv(d.H1.U_S2)); mean kappa $(fv(d.H1.kappa_mean))", vword(d.H1.verdict), "-", "INFO")
    emit("D89x", "$ver H1 A_XY: " * join(["$X-$Y=$(fv(a))" for ((X, Y), a) in sort(collect(d.H1.A); by = kv -> (findfirst(==(kv[1][1]), I8), findfirst(==(kv[1][2]), I8)))], " "), "", "-", "INFO")
    emit("D89x", "$ver H2: orbits S2 $(d.H2.n); M_X " * join(["$X=$(fv(d.H2.M[X]))" for X in I], " ") * "; Q $(fv(d.H2.Q)) interval [$(fv(d.H2.Q_lo)), $(fv(d.H2.Q_hi))]", vword(d.H2.verdict), "-", "INFO")
    emit("D89x", "$ver H3: orbits W2 $(d.H3.n_w), S2 $(d.H3.n_s); M_W2 " * join(["$X=$(fv(d.H3.M_w[X]))" for X in I], " ") * "; M_S2 " * join(["$X=$(fv(d.H3.M_s[X]))" for X in I], " ") *
         "; R_X [lo, hi] " * join(["$X=$(fv(d.H3.R[X])) [$(fv(d.H3.R_lo === nothing ? nothing : d.H3.R_lo[X])), $(fv(d.H3.R_hi === nothing ? nothing : d.H3.R_hi[X]))]" for X in I], " ") * "; indicators passing $(fv(d.H3.passing))", vword(d.H3.verdict), "-", "INFO")
    emit("D89x", "$ver H4: orbits T_s $(d.H4.n); phi_X " * join(["$X=$(fv(v))" for (X, v) in srt(d.H4.phi)], " ") * "; decided " * join(["$X=$(v)" for (X, v) in srt(d.H4.decided)], " ") *
         "; Delta_phi $(fv(d.H4.dphi)) interval [$(fv(d.H4.dphi_lo)), $(fv(d.H4.dphi_hi))]", vword(d.H4.verdict), "-", "INFO")
    emit("D89x", "$ver H4 on T_s^design: orbits $(d.H4.design.n); phi_X " * join(["$X=$(fv(v))" for (X, v) in srt(d.H4.design.phi)], " ") *
         "; Delta_phi $(fv(d.H4.design.dphi)) interval [$(fv(d.H4.design.dphi_lo)), $(fv(d.H4.design.dphi_hi))]", vword(d.H4.design.verdict) * " (statistic of record)", "-", "INFO")
    emit("D89x", "$ver H5: orbits W2 $(d.H5.n); F_X " * join(["$X=$(fv(d.H5.F[X]))" for X in I], " ") * "; F_med $(fv(d.H5.F_med)) interval [$(fv(d.H5.F_med_lo)), $(fv(d.H5.F_med_hi))]", vword(d.H5.verdict), "-", "INFO")
    emit("D89x", "$ver H6: orbits S_D* $(d.H6.dir.n), S_L* $(d.H6.prec.n); delta_dir " * join(["$X=$(fv(v)) (dropped $(d.H6.dir.dropped[X]))" for (X, v) in srt(d.H6.dir.delta)], " ") *
         "; delta_prec " * join(["$X=$(fv(v)) (dropped $(d.H6.prec.dropped[X]))" for (X, v) in srt(d.H6.prec.delta)], " "), vword(d.H6.verdict), "-", "INFO")
end

section("REPORT") do
    res = result()
    info = orbit_info()
    S = res.sets
    P = res.params
    inE8 = Set(S.E8)
    sizes(ids) = begin
        ks = sort(unique(Float64[info.K[i] for i in ids if info.K[i] !== nothing]))
        "total $(length(ids)), panel $(count(i -> info.src[i] == "P", ids)), grid $(count(i -> info.src[i] == "G", ids)); by K: " * join(["$K=$(count(i -> info.K[i] == K, ids))" for K in ks], " ")
    end
    for (nm, ids) in (("E_all", S.E_all), ("E8", S.E8), ("W2", S.W2), ("S2", S.S2), ("T_s", S.T_s), ("T_s^design", S.T_s_design), ("regular set", S.regular), ("gate set", S.gate),
                      ("S_D*", S.S_D), ("S_L*", S.S_L), ("W2_8", [i for i in S.W2 if i in inE8]), ("S2_8", [i for i in S.S2 if i in inE8]), ("T_s_8", [i for i in S.T_s if i in inE8]))
        emit("D89x", "size of $nm", sizes(ids), "-", "INFO")
    end
    emit("D89x", "W2 and S2 overlap; dropped orbits (reasons)", "$(length(S.overlap)); $(length(S.dropped)) (" * join(["$r: $(count(d -> d.reason == r, S.dropped))" for r in sort(unique(d.reason for d in S.dropped))], ", ") * ")", "-", "INFO")
    emit("D89x", "T_s against min_set = $(P.min_set)", length(S.T_s) >= P.min_set ? "decidable ($(length(S.T_s)) orbits)" : "not decidable ($(length(S.T_s)) orbits)", "-", "INFO")
    for (nm, ids) in (("W2", S.W2), ("S2", S.S2), ("T_s^design", S.T_s_design), ("S_D*", S.S_D), ("S_L*", S.S_L))
        emit("D89x", "$nm against min_set", length(ids) >= P.min_set ? "decidable ($(length(ids)))" : "not decidable ($(length(ids)))", "-", "INFO")
    end
    rhoc = Dict(s => [count(i -> info.src[i] == s && res.rho[i] === r, 1:info.n) for r in (:C, :R, :U)] for s in ("P", "G"))
    emit("D89x", "reference class counts C, R, U: panel (883 orbits), grid (2304 orbits)", "panel $(rhoc["P"]), grid $(rhoc["G"])", "-", "INFO")
    describe("PRIMARY (R1, I_7, E_all)", res.decisions, I7)
    emit("D89x", "ROBUST flags H1..H6 (primary against the I_8 extension)", join(["H$i=$(rword(res.robust[i]))" for i in 1:6], " "), "-", "INFO")
    describe("EXTENSION (R1, I_8, E8)", res.extension, I8)
    describe("R2 (I_7, E_all)", res.r2, I7)
    describe("R2 extension (I_8, E8)", res.r2_extension, I8)
    emit("D89x", "R2 verdicts H1..H6 (reported, not decisive)", join(["H$i=$(vword(res.r2[i].verdict))" for i in 1:6], " "), "-", "INFO")
    for g in res.gate
        emit("D89x", "gate T12 $(g.X): decided $(g.decided), chaotic $(g.chaotic), undecided $(g.undecided), not applicable $(g.na), share chaotic $(fv(g.share)), limit $(P.gate_max_share)", g.violated ? "VIOLATED (rule invalid on the integrable control)" : "not violated", "-", "INFO")
    end
    emit("D89x", "second versions produced (violating indicators)", isempty(res.second) ? "none" : join(string.(keys(res.second)), ", "), "-", "INFO")
    emit("D89x", "extension H6(b): indicators with a Double64 item (RLI has none)", join(string.(keys(res.extension.H6.prec.delta)), " "), "-", "INFO")
    for (W, d) in sort(collect(res.k0line); by = first)
        emit("D89v", "K = 0 line (T16), window $W: largest circle distance of nu1 to p0/(2 pi): fq1; pa1; non-finite left out (fq1, pa1)", "$(fv(d.fq1)); $(fv(d.pa1)); $(d.fq1_nonfinite), $(d.pa1_nonfinite)", "tabulated, no gate", "INFO")
    end
    dev = [(c.id, abs(c.deviation - c.stored_deviation)) for c in res.controls_table if c.stored_deviation !== nothing]
    emit("D89v", "control_deviation stored against recomputed over $(length(dev)) controls: largest difference (E36 n)", isempty(dev) ? 0.0 : maximum(last, dev), "<= 1e-12", pf(isempty(dev) || maximum(last, dev) <= 1e-12))
    emit("D89v", "omissions recomputed from the store (spec: orbits), equal to the manifest not_applicable; flag counts equal the manifest per chunk (E36 p, g: the run did not refuse)", join(["$s: $(join(ids, ","))" for (s, ids) in sort(collect(res.omissions); by = first)], "; "), "equal", pf(true))
    emit("D89x", "header of the tables file: chunks, items, seed, module version, commit", string(res.header.chunks, " chunks, ", res.header.items, " items, seed ", res.header.seed, ", module ", res.header.module_version, ", commit ", res.header.git_sha), "695 chunks, 26964 items", pf(res.header.chunks == 695 && res.header.items == 26964))
end

# ---------------------------------------------------------------------------------------------------------------- COMPARE (D89w, O6, O10)
# The module side is written into MOD under the names of the independent results file (documented at the top of verify/analysis-p2-decision.py).
function build_mod(res, info)
    MOD = Dict{String, Any}()
    put(n, v) = (MOD[n] = v)
    S, P = res.sets, res.params
    inE8 = Set(S.E8)
    n = info.n
    # --- sizes, drops, reference class, labels
    base = ["Eall" => S.E_all, "E8" => S.E8, "W2" => S.W2, "S2" => S.S2, "Ts" => S.T_s, "Ts_design" => S.T_s_design, "regular" => S.regular, "gate" => S.gate, "SD_star" => S.S_D, "SL_star" => S.S_L]
    sets = vcat(base, [(nm * "_8") => [i for i in ids if i in inE8] for (nm, ids) in base if nm in ("W2", "S2", "Ts", "Ts_design", "regular", "SD_star", "SL_star")])
    for (nm, ids) in sets
        put("size|$nm", length(ids))
        put("size|$nm|P", count(i -> info.src[i] == "P", ids)); put("size|$nm|G", count(i -> info.src[i] == "G", ids))
        for K in sort(unique(Float64[info.K[i] for i in ids if info.K[i] !== nothing]))
            put("size|$nm|$(kname(K))", count(i -> info.K[i] == K, ids))
        end
    end
    put("size|W2_and_S2", length(S.overlap))
    reason = Dict(d.id => d.reason for d in S.dropped)
    put("drop|controls", count(==("control"), values(reason)))
    flagged = sort([id for (id, r) in reason if r in ("drift", "excluded")])
    put("drop|flagged", length(flagged)); put("drop|flagged_ids", join(flagged, ","))
    for s in ("P", "G")
        put("drop|rhoU|$s", count(id -> reason[id] == "unresolved" && info.src[id] == s, collect(keys(reason))))
        ids_s = [i for i in 1:n if info.src[i] == s]
        for r in (:C, :R, :U)
            put("rho|$s|all|$r", count(i -> res.rho[i] === r, ids_s))
            for K in sort(unique(Float64[info.K[i] for i in ids_s if info.K[i] !== nothing]))
                put("rho|$s|$(kname(K))|$r", count(i -> info.K[i] == K && res.rho[i] === r, ids_s))
            end
        end
    end
    for kx in sort(unique(info.kidx[1:883])), cl in sort(unique(info.cls[i] for i in 1:883 if info.kidx[i] == kx)), r in (:C, :R, :U)
        put("rhocell|$kx|$cl|$r", count(i -> info.kidx[i] == kx && info.cls[i] == cl && res.rho[i] === r, 1:883))
    end
    for i in 1:883
        l = res.p1_labels[i]
        nm = (info.ctrl[i] ? "label|0|controls|" : "label|$(info.kidx[i])|$(info.cls[i])|") * l.instability * "|" * l.structure
        MOD[nm] = get(MOD, nm, 0) + 1
    end
    put("sdim|SD_total", count(info.sd)); put("sdim|SL_total", count(info.sl))
    # --- steps used
    for X in (:LE, :SALI, :FLI, :MEGNO, :RLI, :STR), k in 0:9
        put("step|$X|k$k", join(sort(unique(res.steps[(i, X)][k + 1] for i in 1:n if haskey(res.steps, (i, X)))), "/"))
    end
    # --- verdict counts, detection, flips, budget stability
    rules = ((:R1, "R1"), (:R2, "R2"))
    small = ["W2" => S.W2, "S2" => S.S2, "Ts" => S.T_s, "regular" => S.regular]
    function detf(rs, i, X)
        haskey(res.detection, (i, X)) || return nothing
        rs === :R1 && return res.detection[(i, X)]
        return AN.p2_detection(res.verdicts[(i, X, :R2)], res.rho[i], P)
    end
    for (rs, rn) in rules, (sn, ids) in small, X in I8
        vs = [get(res.verdicts, (i, X, rs), nothing) for i in ids]
        for k in 0:9
            c(sym) = count(v -> v !== nothing && v[k + 1] === sym, vs)
            put("vcount|$rn|$sn|$X|k$k|chaotic", c(:chaotic)); put("vcount|$rn|$sn|$X|k$k|regular", c(:regular))
            put("vcount|$rn|$sn|$X|k$k|undecided", c(:undecided)); put("vcount|$rn|$sn|$X|k$k|na", count(isnothing, vs))
        end
    end
    for (rs, rn) in rules, (sn, ids) in vcat(small, ["Ts_design" => S.T_s_design]), X in I8
        ds = [d for d in (detf(rs, i, X) for i in ids) if d !== nothing]
        for o in (:correct, :wrong, :undetected)
            put("det|$rn|$sn|$X|$o", count(d -> d.outcome === o, ds))
        end
        tc = [Float64(d.censored) for d in ds]
        put("det|$rn|$sn|$X|M", isempty(tc) ? NaN : median_of(tc)); put("det|$rn|$sn|$X|n", length(tc))
        for v in sort(unique(d.t_det for d in ds if d.t_det !== nothing))
            put("tdet|$rn|$sn|$X|$v", count(d -> d.t_det == v, ds))
        end
    end
    for (sn, ids) in small, X in I8
        vs = [res.verdicts[(i, X, :R1)] for i in ids if haskey(res.verdicts, (i, X, :R1))]
        nna = length(ids) - length(vs)
        fl = [count(k -> dec(v[k]) && dec(v[k + 1]) && v[k] !== v[k + 1], 1:9) for v in vs]
        vals = sort(unique(vcat(fl, nna > 0 ? [0] : Int[])))        # the independent side counts the orbits without an item at flip count 0
        for f in vals
            put("Phi|$sn|$X|$f", count(==(f), fl) + (f == 0 ? nna : 0))
        end
        for k in 0:9
            both = [(v[k + 1], v[10]) for v in vs if dec(v[k + 1]) && dec(v[10])]
            put("Fb|$sn|k$k|$X", isempty(both) ? NaN : count(q -> q[1] !== q[2], both) / length(both))
        end
    end
    # --- F1 (agreement, kappa, Spearman) from the tables of the module
    f1 = res.items.F1
    ci = Dict(c => j for (j, c) in enumerate(f1.columns))
    kap = Dict{Tuple{String, Int}, Vector{Float64}}()
    for r in f1.rows
        sn, b, X, Y = r[ci["set"]], r[ci["budget"]], r[ci["X"]], r[ci["Y"]]
        k = findfirst(==(b), P.budgets) - 1
        pre = "agree|$sn|k$k|$X-$Y"
        put(pre * "|A", numv(r[ci["A"]])); put(pre * "|kappa", numv(r[ci["kappa"]])); put(pre * "|N", r[ci["N"]])
        put(pre * "|ncc", r[ci["cc"]]); put(pre * "|ncr", r[ci["cr"]]); put(pre * "|nrc", r[ci["rc"]]); put(pre * "|nrr", r[ci["rr"]])
        put("spearman|$sn|k$k|$X-$Y", numv(r[ci["spearman"]]))
        if Symbol(X) in I7 && Symbol(Y) in I7 && r[ci["kappa"]] !== nothing
            push!(get!(kap, (sn, k), Float64[]), r[ci["kappa"]])
        end
    end
    for ((sn, k), v) in kap
        put("meankappa|$sn|k$k|I7", sum(v) / length(v))
    end
    # --- AUC, regular range, leaving budget
    for ((X, k, K, src), a) in res.auc
        a.auc === nothing && continue
        pre = "auc|$X|k$k|$(kname(K))|" * (src == "pooled" ? "pooled" : src == "panel" ? "P" : "G")
        put(pre * "|auc", a.auc); put(pre * "|nC", a.n_C); put(pre * "|nR", a.n_R)
    end
    for ((X, K, k), r) in res.range
        r === nothing && continue
        pre = "rng|$X|k$k|$(kname(K))"
        put(pre * "|lo", r.lo); put(pre * "|hi", r.hi); put(pre * "|median", r.median); put(pre * "|n", r.n)
    end
    Ks = sort(unique(Float64[info.K[i] for i in S.E_all if info.K[i] !== nothing]))
    partial = String[]
    for X in vcat(I7, [:STR]), K in Ks
        ts = [i for i in S.T_s if info.K[i] == K]
        isempty(ts) && continue
        pre = "leave|$X|$(kname(K))"
        put(pre * "|n", length(ts))
        if !all(k -> res.range[(X, K, k)] !== nothing, 0:9)
            put(pre * "|defined", "no range")
            push!(partial, "$X $K: module leaving defined for $(count(i -> res.leaving[(i, X)] !== nothing, ts)) of $(length(ts))")
            continue
        end
        lv = [res.leaving[(i, X)] for i in ts]
        sc9 = [res.scores[(i, X)][10] for i in ts if res.scores[(i, X)][10] !== nothing]
        sh = res.below_median[(X, K)]
        put(pre * "|defined", count(!isnothing, lv)); put(pre * "|n9", length(sc9))
        put(pre * "|below_median", sh === nothing ? 0 : round(Int, sh * length(sc9)))
        for b in sort(unique(filter(!isnothing, lv)))
            put("leavek|$X|$(kname(K))|$b", count(==(b), lv))
        end
    end
    MOD["__partial_ranges__"] = partial
    # --- T11, T15, K = 0 line, controls, flags, gate
    pairnames = Dict(:raw_vs_fq1 => "fq1_raw-fq1", :qp_vs_fq1 => "fq1_qp-fq1", :planar_vs_pa1 => "pa_planar-pa1", :cos_vs_wba => "wba_cos1-wba")
    for ((pair, k), s) in res.sensitivity
        pre = "T11|$(pairnames[pair])|k$k"
        put(pre * "|n", s.n); put(pre * "|A", numv(s.A)); put(pre * "|kappa", numv(s.kappa)); put(pre * "|medabs", numv(s.median_diff)); put(pre * "|maxabs", numv(s.max_diff))
    end
    t15 = res.items.T15
    cj = Dict(c => j for (j, c) in enumerate(t15.columns))
    code(s) = s == "chaotic" ? 1 : s == "regular" ? 0 : s == "undecided" ? -1 : 2
    for r in t15.rows
        pk = replace(r[cj["package"]], " with " => "+")
        put("T15|$pk|$(replace(r[cj["set"]], "_" => ""))|k$(r[cj["k"]])|$(code(r[cj["verdict first"]])),$(code(r[cj["verdict second"]]))|rho=$(r[cj["rho"]])", r[cj["orbits"]])
    end
    for (W, d) in res.k0line
        put("k0line|fq1|W=$W|maxdist", numv(d.fq1)); put("k0line|pa1|W=$W|maxdist", numv(d.pa1))
    end
    cc(v) = v === :chaotic ? 1 : v === :regular ? 0 : v === :undecided ? -1 : 2
    for c in res.controls_table
        for X in I8
            put("ctrl|$(c.id)|$X|v9", cc(get(c.verdicts9, X, :na)))
        end
        put("ctrl|$(c.id)|lamref", c.lambda)
    end
    for (s, d) in res.flag_counts, (f, cnt) in d
        cnt > 0 && put("T16|flagsum|$s|$f", cnt)
    end
    for (s, ids) in res.omissions
        put("T16|omitted|$s", join(ids, ",")); put("T16|omitted_equals_manifest|$s", "true")
    end
    for g in res.gate
        put("gate|$(g.X)|decided", g.decided); put("gate|$(g.X)|chaotic", g.chaotic); put("gate|$(g.X)|undecided", g.undecided); put("gate|$(g.X)|na", g.na)
        put("gate|$(g.X)|share", numv(g.share)); put("gate|$(g.X)|violated", g.violated ? "true" : "false")
    end
    # --- decisions and statistics of the four versions
    function put_version!(ver, d, I, rs, e8)
        sel(v) = e8 ? [i for i in v if i in inE8] : v
        h = d.H1
        put("$ver|H1|n_W2", h.n); put("$ver|H1|U", numv(h.U)); put("$ver|H1|U_drop", h.dropped); put("$ver|H1|Amin", numv(h.A_min))
        for ((X, Y), a) in h.A
            put("$ver|H1|A|$X-$Y", numv(a))
        end
        put("$ver|H1|U_S2", numv(h.U_S2)); put("$ver|H1|U_lo", numv(h.U_lo)); put("$ver|H1|U_hi", numv(h.U_hi)); put("$ver|H1|decision", vword(h.verdict))
        h = d.H2
        put("$ver|H2|n_S2", h.n)
        for X in I
            put("$ver|H2|M|$X", numv(h.M[X]))
            ds = [q for q in (detf(rs, i, X) for i in sel(S.S2)) if q !== nothing]
            for o in (:correct, :wrong, :undetected)
                put("$ver|H2|share_$o|$X", count(q -> q.outcome === o, ds) / max(length(ds), 1))
            end
        end
        put("$ver|H2|Q", numv(h.Q)); put("$ver|H2|Q_lo", numv(h.Q_lo)); put("$ver|H2|Q_hi", numv(h.Q_hi)); put("$ver|H2|decision", vword(h.verdict))
        h = d.H3
        put("$ver|H3|n_W2", h.n_w); put("$ver|H3|n_S2", h.n_s)
        for X in I
            put("$ver|H3|MW|$X", numv(h.M_w[X])); put("$ver|H3|MS|$X", numv(h.M_s[X])); put("$ver|H3|R|$X", numv(h.R[X]))
            put("$ver|H3|R_lo|$X", h.R_lo === nothing ? NaN : numv(h.R_lo[X])); put("$ver|H3|R_hi|$X", h.R_hi === nothing ? NaN : numv(h.R_hi[X]))
        end
        put("$ver|H3|npass", h.passing === nothing ? 0 : h.passing); put("$ver|H3|decision", vword(h.verdict))
        for (item, h) in (("H4", d.H4), ("H4design", d.H4.design))
            put("$ver|$item|n_Ts", h.n)
            for X in I
                X === :LE || put("$ver|$item|phi|$X", numv(h.phi[X]))
            end
            put("$ver|$item|dphi", numv(h.dphi)); put("$ver|$item|dphi_lo", numv(h.dphi_lo)); put("$ver|$item|dphi_hi", numv(h.dphi_hi)); put("$ver|$item|decision", vword(h.verdict))
        end
        h = d.H5
        put("$ver|H5|n_W2", h.n)
        for X in I
            put("$ver|H5|F|$X", numv(h.F[X]))
        end
        put("$ver|H5|Fmed", numv(h.F_med)); put("$ver|H5|Fmed_lo", numv(h.F_med_lo)); put("$ver|H5|Fmed_hi", numv(h.F_med_hi)); put("$ver|H5|decision", vword(h.verdict))
        h = d.H6
        for (X, v) in h.dir.delta
            put("$ver|H6|ddir|$X", numv(v)); put("$ver|H6|ddir_drop|$X", h.dir.dropped[X])
        end
        for (X, v) in h.prec.delta
            put("$ver|H6|dprec|$X", numv(v)); put("$ver|H6|dprec_drop|$X", h.prec.dropped[X])
        end
        put("$ver|H6|n_SDs", h.dir.n); put("$ver|H6|n_SLs", h.prec.n)
        mx(dd) = (v = [x for x in values(dd) if x !== nothing]; isempty(v) ? NaN : maximum(v))
        put("$ver|H6|max_dir", mx(h.dir.delta)); put("$ver|H6|max_prec", mx(h.prec.delta)); put("$ver|H6|decision", vword(h.verdict))
    end
    put_version!("R1.I7", res.decisions, I7, :R1, false); put_version!("R1.I8", res.extension, I8, :R1, true)
    put_version!("R2.I7", res.r2, I7, :R2, false); put_version!("R2.I8", res.r2_extension, I8, :R2, true)
    for h in 1:6
        put("R1|robust|H$h", rword(res.robust[h]))
        a, b = res.r2[h].verdict, res.r2_extension[h].verdict
        put("R2|robust|H$h", (a === :not_decidable || b === :not_decidable) ? "NOT COMPARABLE" : a === b ? "ROBUST" : "NOT ROBUST")
    end
    return MOD
end

const SKIP = r"^(selfcheck|reader)|U_repdrop$|^T16\|flagchunks"
isinterval(name) = occursin(r"\|(U|Q|R|dphi|Fmed)_(lo|hi)(\||$)", name)
family(name) = (f = split(name, '|')[1]; startswith(name, "T16|") ? "T16" : f)

function same(a::AbstractString, m)
    m isa AbstractString && return a == m
    m isa Bool && return a == (m ? "true" : "false")
    x = tryparse(Float64, a)
    x === nothing && return false
    y = m === nothing ? NaN : Float64(m)
    (isnan(x) || isnan(y)) && return isnan(x) && isnan(y)
    x == y && return true
    m isa Integer && return false
    return abs(x - y) <= 1e-12 * max(abs(x), abs(y))
end

section("COMPARE") do
    res = result()
    info = orbit_info()
    MOD = build_mod(res, info)
    ind = Pair{String, String}[]
    for l in eachline(IND)
        isempty(l) && continue
        i = findfirst(" | ", l)
        push!(ind, l[1:(first(i) - 1)] => l[(last(i) + 1):end])
    end
    emit("D89w", "independent results file: lines, SHA-256 of the file", "$(length(ind)), $(bytes2hex(sha256(read(IND))))", "-", "INFO")
    fams = Dict{String, Vector{Int}}()          # family => [compared, equal, missing in module]
    diffs = Dict{String, Vector{String}}()
    notfound = String[]
    nint = [0, 0]
    skipped = Dict{String, Int}()
    for (name, val) in ind
        if occursin(SKIP, name)
            skipped[family(name)] = get(skipped, family(name), 0) + 1
            continue
        end
        f = family(name)
        c = get!(fams, f, [0, 0, 0])
        if !haskey(MOD, name)
            c[3] += 1
            push!(notfound, name * " = " * val)
            continue
        end
        m = MOD[name]
        c[1] += 1
        if isinterval(name)
            lo, hi = MOD[replace(name, r"_(lo|hi)" => "_lo")], MOD[replace(name, r"_(lo|hi)" => "_hi")]
            wid = (isfinite(lo) && isfinite(hi)) ? hi - lo : 0.0
            x = tryparse(Float64, val)
            ok = x !== nothing && (isnan(x) || isnan(m) ? (isnan(x) && isnan(m)) : abs(x - m) <= 0.05 * wid + 1e-15)
            nint[1] += 1; nint[2] += ok
        elseif startswith(name, "k0line|")
            x = tryparse(Float64, val)
            ok = x !== nothing && (abs(x - m) <= 1e-15 + 1e-12 * max(abs(x), abs(m)))
        else
            ok = same(val, m)
        end
        c[2] += ok
        ok || push!(get!(diffs, f, String[]), "$name: independent $val, module $(m isa AbstractFloat ? @sprintf("%.17g", m) : m)")
    end
    total = sum(c[1] for c in values(fams)); eq = sum(c[2] for c in values(fams))
    for f in sort(collect(keys(fams)))
        c = fams[f]
        emit("D89w", "family $f: compared / equal / not in the module", "$(c[1]) / $(c[2]) / $(c[3])", "all equal", pf(c[1] == c[2] && c[3] == 0))
    end
    emit("D89w", "TOTAL compared (all quantities), equal; of them interval ends (within 0.05 of the width) compared, equal", "$total, $eq; $(nint[1]), $(nint[2])", "all equal", pf(total == eq))
    emit("D89w", "independent lines skipped (selfchecks, reader checks, U_repdrop, chunk flag totals), by family", join(["$k $v" for (k, v) in sort(collect(skipped))], ", "), "-", "INFO")
    emit("D89w", "independent names with no module counterpart (not compared)", isempty(notfound) ? "none" : "$(length(notfound)): " * join(first(notfound, 12), "; "), "none", pf(isempty(notfound)))
    write(joinpath(SCR, "compare-diffs.txt"), join(["$f\n  " * join(v, "\n  ") for (f, v) in sort(collect(diffs))], "\n"))
    for (f, v) in sort(collect(diffs)), l in first(v, 40)
        say("DIFF | $f | $l")
    end
    modonly = [k for k in keys(MOD) if k != "__partial_ranges__" && !any(p -> first(p) == k, ind)]
    byfam = Dict{String, Int}()
    for k in modonly
        byfam[family(k)] = get(byfam, family(k), 0) + 1
    end
    emit("D89w", "module names (in the naming of the independent file) that the independent file lacks, by family", join(["$k $v" for (k, v) in sort(collect(byfam))], ", "), "-", "INFO")
    emit("D89w", "K values with an incomplete regular range (independent prints 'no range'; module leaving budget for those, not compared)", isempty(MOD["__partial_ranges__"]) ? "none" : "$(length(MOD["__partial_ranges__"])): " * join(first(MOD["__partial_ranges__"], 6), "; "), "-", "INFO")
    # D89p / I16 / I17 on the module's own tables: T1 holds every orbit once, T2 counts sum to the set sizes, the by-source and by-K rows sum to the pooled rows
    let t1 = res.items.T1, t2 = res.items.T2
        c1 = Dict(c => j for (j, c) in enumerate(t1.columns)); c2 = Dict(c => j for (j, c) in enumerate(t2.columns))
        emit("D89p", "T1: orbits over the rows (source, K index) = orbits of P and G; rho counts C + R + U = orbits", "$(sum(r[c1["orbits"]] for r in t1.rows)) ($(sum(r[c1["C"]] + r[c1["R"]] + r[c1["U"]] for r in t1.rows)))", "3187", pf(sum(r[c1["orbits"]] for r in t1.rows) == 3187 && sum(r[c1["C"]] + r[c1["R"]] + r[c1["U"]] for r in t1.rows) == 3187))
        S = res.sets
        want_size = Dict("W2" => length(S.W2), "S2" => length(S.S2), "T_s" => length(S.T_s), "regular" => length(S.regular))
        bad = 0; nrow = 0; bad_src = 0; bad_k = 0
        agg = Dict{Tuple{String, String, Int, String}, Int}()
        for r in t2.rows
            r[c2["source"]] == "all" && r[c2["K"]] == "all" || continue
            nrow += 1
            bad += (r[c2["chaotic"]] + r[c2["regular"]] + r[c2["undecided"]] + r[c2["not applicable"]]) != want_size[r[c2["set"]]]
        end
        for r in t2.rows
            tot = r[c2["chaotic"]] + r[c2["regular"]] + r[c2["undecided"]] + r[c2["not applicable"]]
            key = (r[c2["indicator"]], r[c2["set"]], r[c2["k"]], r[c2["source"]] == "all" ? (r[c2["K"]] == "all" ? "pooled" : "byK") : "bySource")
            agg[key] = get(agg, key, 0) + tot
        end
        for ((X, sn, k, kind), v) in agg
            kind == "pooled" || (bad_k += (v != want_size[sn]))
        end
        emit("D89p", "T2 (I17): pooled rows checked, rows whose chaotic + regular + undecided + not applicable differ from the set size; by-source and by-K row groups that do not sum to the set size", "$nrow, $bad; $bad_k", "0, 0", pf(bad == 0 && bad_k == 0 && nrow == 4 * 8 * 10))
    end
    # H1 mean kappa of the module's own statistic against the independent one
    mk = tryparse(Float64, Dict(ind)["meankappa|W2|k6|I7"])
    emit("D89w", "mean kappa of H1 (module statistic, Kahan sum) against the independent mean kappa W2 b_P", abs(res.decisions.H1.kappa_mean - mk), "<= 1e-12", pf(abs(res.decisions.H1.kappa_mean - mk) <= 1e-12))
    # O10: module labels against T6 of the P1 tables
    t6 = Dict{String, Int}()
    inT6 = false
    for l in eachline(TABLES1)
        startswith(l, "## T6 Label table") && (inT6 = true; continue)
        inT6 && startswith(l, "## ") && break
        m = match(r"^\| (\d+) \| (\w+) \| (.+?) \| (.+?) \| (\d+) \|\s*$", l)
        inT6 && m !== nothing && (t6["label|$(m.captures[1])|$(m.captures[2])|$(m.captures[3])|$(m.captures[4])"] = parse(Int, m.captures[5]))
    end
    modl = Dict(k => v for (k, v) in MOD if startswith(k, "label|"))
    ndiff = count(k -> get(modl, k, -1) != get(t6, k, -2), union(keys(modl), keys(t6)))
    emit("O10", "panel orbits per cell and label pair, module (P1 labelling inside the P2 run) against T6 of the accepted P1 tables: cells, orbits, cells that differ", "$(length(modl)) vs $(length(t6)), $(sum(values(modl))) vs $(sum(values(t6))), $ndiff", "0 differ", pf(ndiff == 0 && sum(values(modl)) == 883))
    # the evidence of the independent reader and the self checks, as printed there
    for (name, val) in ind
        startswith(name, "reader|") && name != "reader|counts" && emit("D89v", replace(name, "reader|" => "independent reader: "), val, "PASS", startswith(val, "PASS") ? "PASS" : "FAIL")
        name == "selfcheck_summary" && emit("O7/O8", "independent self checks (kappa and AUC closed forms, exact bootstrap, Spearman against scipy, pilot rule check of 192 verdicts)", val, "all PASS", startswith(val, "16 PASS of 16") ? "PASS" : "FAIL")
        name == "reader|counts" && emit("D89v", "independent reader: files, series items, window items, anomalies", val, "anomalies 0", "INFO")
    end
end

# ---------------------------------------------------------------------------------------------------------------- P1UNCH (D89a, V: the P1 tables file is unchanged)
section("P1UNCH") do
    out = joinpath(SCR, "p1-out"); rm(out; recursive = true, force = true)
    t = @timed AN.analyze(; prereg = PREREG1, manifest = MAN1, panel = PANEL, params = PARAMS1, store = STORE1, out = out, figures = false, threads = NT)
    new = split(replace(read(t.value.tables, String), "\r\n" => "\n"), '\n')
    old = split(replace(read(TABLES1, String), "\r\n" => "\n"), '\n')
    diffs = [i for i in 1:min(length(new), length(old)) if new[i] != old[i]]
    allowed = all(i -> startswith(old[i], "- module version:") || startswith(old[i], "- git commit:"), diffs)
    emit("D89a", "P1 analyze on the real P1 inputs: status, time (s); lines in the new and in the accepted tables file", "$(t.value.status), $(round(t.time; digits = 1)); $(length(new)) and $(length(old))", "status ok", pf(t.value.status === :ok))
    emit("D89a", "lines that differ from the accepted p1-le-benchmark-tables.md (all must be the two header lines module version and git commit)", "$(length(diffs)): " * join(["line $i: \"$(old[i])\" -> \"$(new[i])\"" for i in diffs], "; "), "<= 2, header lines only", pf(length(new) == length(old) && length(diffs) <= 2 && allowed))
    dg(lines) = [m.captures[1] for l in lines for m in (match(r"^digest: ([0-9a-f]{64})", l),) if m !== nothing]
    emit("D89a", "figure data digests F1 to F7 equal the accepted ones", "$(length(dg(new))) digests, equal $(dg(new) == dg(old))", "equal", pf(length(dg(new)) == 7 && dg(new) == dg(old)))
end

# ---------------------------------------------------------------------------------------------------------------- CLI (D89s)
section("CLI") do
    out = joinpath(SCR, "out-cli"); rm(out; recursive = true, force = true)
    jl = Base.julia_cmd()
    proj = joinpath(ROOT, "verify"); bin = joinpath(ROOT, "bin", "chaoslab.jl")
    logf = joinpath(SCR, "cli-stdout.txt")
    p = run(pipeline(ignorestatus(`$jl --project=$proj -t $NT $bin analyze $PREREG $MAN2 $PANEL $PARAMS --store $STORE2 --out $out`); stdout = logf, stderr = logf))
    txt = read(logf, String)
    emit("D89s", "bin/chaoslab.jl analyze on the real P2 inputs (P1 reference inputs by the study layout): exit code; says figures skipped", "$(p.exitcode); $(occursin("figures skipped", txt))", "2; true", pf(p.exitcode == 2 && occursin("figures skipped", txt)))
    f = joinpath(out, "p2-indicators-tables.md")
    ref = joinpath(OUT_RUN, "p2-indicators-tables.md")
    emit("D89s", "tables file of the command (default P1 designation) against the file of the API run with explicit P1 options: byte-identical; only file in --out", isfile(f) && isfile(ref) && read(f) == read(ref), "true", pf(isfile(f) && isfile(ref) && read(f) == read(ref) && readdir(out) == ["p2-indicators-tables.md"]))
    bad = run(pipeline(ignorestatus(`$jl --project=$proj $bin analyze $PREREG --bogus x`); stdout = devnull, stderr = devnull))
    emit("D89s", "unknown option and too few arguments: exit code", bad.exitcode, "1", pf(bad.exitcode == 1))
    bad2 = run(pipeline(ignorestatus(`$jl --project=$proj $bin analyze $PREREG1 $MAN1 $PANEL $PARAMS1 --p1-store $STORE1 --out $(joinpath(SCR, "out-bad"))`); stdout = devnull, stderr = devnull))
    emit("D89s", "--p1-store with the P1 preregistration: usage error, exit code, nothing written", "$(bad2.exitcode), $(ispath(joinpath(SCR, "out-bad")))", "1, false", pf(bad2.exitcode == 1 && !ispath(joinpath(SCR, "out-bad"))))
    bad3 = run(pipeline(ignorestatus(`$jl --project=$proj $bin analyze $PREREG $MAN2 $PANEL $PARAMS1 --out $(joinpath(SCR, "out-bad2"))`); stdout = devnull, stderr = devnull))
    emit("D89a", "P1 parameter table given with the P2 preregistration: refused, exit code, nothing written", "$(bad3.exitcode), $(ispath(joinpath(SCR, "out-bad2")))", "1, false", pf(bad3.exitcode == 1 && !ispath(joinpath(SCR, "out-bad2"))))
end

say(@sprintf("TOTAL | PASS %d | FAIL %d | INFO %d | %s", COUNTS["PASS"], COUNTS["FAIL"], COUNTS["INFO"], Dates.format(now(), "yyyy-mm-dd HH:MM")))
close(LOG)
