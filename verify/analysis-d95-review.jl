# Stage 04 review fixes of Amendment D95 (test side of study p3-abstention): section REV, one probe per fixed finding (C1, H2, H3, H4, M5 to M10, L11 to L14, L16, L17, V1).
# Run from the repo root (output form: `REV <id> | <what> | PASS` or `FAIL`):
#   julia --project=verify -t 4 verify/analysis-d95-review.jl [<scratch root>] [<campaign folder name, default campaign>]
# Reads the toy campaign <root>/<campaign> (INPUTS.md there; read only); writes only <root>/review/. Nothing in data/, src/, test/ is written.
using Chaoslab, SHA, TOML
const AN = Chaoslab.Analysis
const XP = Chaoslab.Experiments

const ROOT = length(ARGS) >= 1 ? ARGS[1] : "C:/Users/tesse/AppData/Local/Temp/claude/C--Users-tesse-Documents-Vibe-Coding-Chaos/5df15ee1-47df-4875-85bd-1bd4d7c665a6/scratchpad/d95v"
const CAMP = joinpath(ROOT, length(ARGS) >= 2 ? ARGS[2] : "campaign")
const REV = joinpath(ROOT, "review")
rm(REV; recursive = true, force = true); mkpath(REV)
P(x) = joinpath(CAMP, x)
base = (prereg = P("inputs/p3-abstention-preregistration.md"), params = P("inputs/p3-abstention-calibration-parameters.toml"), test_params = P("inputs/p3-abstention-test-parameters.toml"),
        panel = P("inputs/p3-abstention-test-panel.toml"), p2_manifest = P("p2/p2-manifest.json"), p2_prereg = P("p2/p2-preregistration.md"), p2_params = P("p2/p2-parameters.toml"),
        calibration_record = P("inputs/p3-abstention-calibration-record.json"))
manifest, p1p, p1q = P("inputs/p3-abstention-test-manifest.json"), P("p1/p1-preregistration.md"), P("p1/p1-parameters.toml")
LOG3 = read(P("round-log.txt"), String)
ents = [String(e) for e in split(LOG3, r"(?<=\n)(?=round=)") if !isempty(e)]

nresult = Ref(0); nfail = Ref(0)
function emit(id, what, ok)
    nresult[] += 1; ok || (nfail[] += 1)
    println("REV $id | $what | ", ok ? "PASS" : "FAIL"); flush(stdout)
end
# the message of the ArgumentError that f throws, "" when it returns, "OTHER: ..." for any other error
function refusal(f)
    try
        f(); return ""
    catch e
        return e isa ArgumentError ? e.msg : "OTHER: " * first(sprint(showerror, e), 160)
    end
end
has(msg, parts...) = !isempty(msg) && !startswith(msg, "OTHER") && all(p -> occursin(p, msg), parts)
probe(f, id, what) = try; emit(id, what, f() === true); catch e; emit(id, what * " (probe error: " * first(sprint(showerror, e), 120) * ")", false); end

cr(store, log; kw...) = AN.count_round(; base..., store, round_log = log, kw...)
an(store, log, out; kw...) = AN.analyze(; base..., store, round_log = log, manifest, p1_prereg = p1p, p1_params = p1q, out, kw...)
ctx(store) = AN.p3t_context(; base..., store, on_read = nothing, command = "count-round")
chunkfiles(dir) = sort([f for f in readdir(dir) if occursin(r"^chunk-\d{5}\.h5$", f)])
function cpstore(dst, files)       # a store with the given chunk files of the campaign store
    rm(dst; recursive = true, force = true); mkpath(dst)
    for f in files; cp(joinpath(P("store"), f), joinpath(dst, f)); end
    return dst
end
function refround(f)               # the round of a ref chunk, 0 for another arm
    a = XP.read_chunk_attrs(joinpath(P("store"), f)); m = match(r"cell=k\d+:r(\d+)", String(a["chunk_key"]))
    return "ref" in String.(a["arms"]) && m !== nothing ? parse(Int, m.captures[1]) : 0
end
refs = Dict(f => refround(f) for f in chunkfiles(P("store")))
stage(r, name) = cpstore(joinpath(REV, name), [f for (f, k) in sort!(collect(refs)) if 1 <= k <= r])
rewrite(path, edit) = (ch = XP.read_chunk(path); edit(ch.attrs, ch.orbits); rm(path); XP.write_chunk(path, ch.attrs, ch.orbits))
setroot(path, name, v) = XP.HDF5.h5open(h -> (XP.HDF5.delete_attribute(h, name); XP.HDF5.write_attribute(h, name, v)), path, "r+")
function dupchunk(store, src, dst, cell)   # a copy of chunk file src as chunk file dst with another cell in the chunk key
    cp(joinpath(store, src), joinpath(store, dst))
    k = String(XP.read_chunk_attrs(joinpath(store, dst))["chunk_key"])
    setroot(joinpath(store, dst), "chunk_key", replace(k, r"cell=\S+" => "cell=" * cell))
end
L1ref = first(f for (f, k) in sort!(collect(refs)) if k == 1 && String(XP.read_chunk_attrs(P("store/" * f))["method"]) == "L1")

# ---------------------------------------------------------------- C1
probe("C1", "record in the shape of p3_record (path j = 0, 1, 2 with a null power at j = 0; power_expected a row): the planned powers are read at the target point and in the row; 1e-13 apart and null refuse") do
    p = AN.p3t_read_params(base.test_params; prereg = base.prereg, record = base.calibration_record, params = base.params)
    runs = 1000
    hbrow(j, pw) = Dict{String, Any}("j" => j, "rho" => "1/1", "allowed" => pw !== nothing, "successes" => pw === nothing ? nothing : round(Int, pw * runs), "runs" => runs, "power" => pw, "se" => pw === nothing ? nothing : 0.05)
    function mk(; hb = p.power_hb_planned, hc = p.power_hc_planned, hbs = round(Int, p.power_hb_planned * runs), hcs = round(Int, p.power_hc_planned * runs), nullhc = false)
        rec = deepcopy(XP.read_json(base.calibration_record))
        rec["pp_hb"]["targets"]["j"] = 2
        rec["pp_hb"]["path"] = Any[hbrow(0, nothing), hbrow(1, 0.1), hbrow(2, hb)]
        rec["pp_hb"]["path"][3]["successes"] = hbs
        rec["pp_hc"]["power_expected"] = nullhc ? nothing : Dict{String, Any}("nf" => p.target_f, "successes" => hcs, "runs" => runs, "power" => hc, "se" => 0.01)
        return rec
    end
    chk(rec) = refusal(() -> AN._p3t_check_record_values("test-parameters.toml", p, rec))
    ok = chk(mk())
    near = chk(mk(; hb = p.power_hb_planned * (1 + 1e-13), hbs = round(Int, p.power_hb_planned * runs) + 1))
    nul = chk(mk(; nullhc = true))
    hcoff = chk(mk(; hc = p.power_hc_planned * (1 + 1e-13), hcs = round(Int, p.power_hc_planned * runs) + 1))
    println("    C1 detail: shape: \"", first(ok, 90), "\" | near: ", first(near, 90), " | null: ", first(nul, 90), " | hc 1e-13: ", first(hcoff, 90))
    return ok == "" && has(near, "power_hb_planned") && has(nul, "power_hc_planned") && has(hcoff, "power_hc_planned")
end

# ---------------------------------------------------------------- H2
probe("H2", "a second ref item with the same content (another cell) refuses naming the orbit and the item (before: the second silently overwrote the first)") do
    s = stage(1, "h2-store"); dupchunk(s, L1ref, "chunk-00999.h5", "k9:r1")
    m = refusal(() -> cr(s, joinpath(REV, "h2.log")))
    println("    H2 detail: ", first(m, 140))
    return has(m, "stored twice", "L1") && occursin(r"orbit \d+", m) && !isfile(joinpath(REV, "h2.log"))
end

# ---------------------------------------------------------------- H3
probe("H3", "the drift flag on an ind-arm item (fq1) drops the orbit with reason drift (analyze, drops of the results file)") do
    dropsof(f) = Dict(Int(d["id"]) => d["reason"] for d in XP.read_json(f)["sets"]["drops"])
    d0 = dropsof(joinpath(CAMP, "module-out", "p3-abstention-test-results.json"))
    fq = last(f for f in chunkfiles(P("store")) if String(XP.read_chunk_attrs(P("store/" * f))["method"]) == "fq1")
    s = cpstore(joinpath(REV, "h3-store"), chunkfiles(P("store")))
    etest = Set(Int.(XP.read_json(joinpath(CAMP, "module-out", "p3-abstention-test-results.json"))["sets"]["ids"]["E_test"]))
    id = first(o.id for o in XP.read_chunk(joinpath(s, fq)).orbits if o.id in etest)
    rewrite(joinpath(s, fq), (a, os) -> push!(first(o for o in os if o.id == id).attrs["flags"], "drift"))
    out = joinpath(REV, "h3-out"); log = joinpath(REV, "h3.log"); write(log, LOG3)
    m = refusal(() -> an(s, log, out))
    d1 = m == "" ? dropsof(joinpath(out, "p3-abstention-test-results.json")) : Dict{Int, String}()
    println("    H3 detail: orbit $id; refusal \"", first(m, 100), "\"; reason now: ", get(d1, id, "not dropped"))
    return m == "" && get(d1, id, "") == "drift"
end

# ---------------------------------------------------------------- H4, M5
c1 = ctx(stage(1, "r1-store")); st1 = AN._p3t_extend_state(nothing, c1, 1, 1, nothing)
crow = first(r for r in st1.rows if r.rho === :C && !r.excluded && !r.drift)
probe("H4", "a non-finite lambda_ns of an orbit that needs it refuses naming the orbit (orbit $(crow.id))") do
    it1 = st1.items[(crow.id, "L1")]; i = AN._p3_step(it1, c1.q.n_s)
    s = stage(1, "h4-store"); path = joinpath(s, "chunk-" * lpad(it1.chunk, 5, '0') * ".h5")
    rewrite(path, (a, os) -> (first(o for o in os if o.id == crow.id).lyap["L1"][i, 1] = NaN))
    m = refusal(() -> cr(s, joinpath(REV, "h4.log")))
    println("    H4 detail: ", first(m, 140))
    return has(m, "orbit $(crow.id)") && !isfile(joinpath(REV, "h4.log"))
end
probe("M5", "lambda_ref and lambda_ns are the exact binary64 values (Rational{BigInt} of the Float64), not the shortest decimal text") do
    n, odd = 0, 0
    for r in st1.rows
        (r.excluded || r.lref === nothing) && continue
        it1 = st1.items[(r.id, "L1")]; a, b = it1.lam[AN._p3_step(it1, c1.q.N_1)], it1.lam[AN._p3_step(it1, c1.q.n_s)]
        n += 1; odd += AN.p3_rat(a) != Rational{BigInt}(a) || AN.p3_rat(b) != Rational{BigInt}(b)
        (r.lref == Rational{BigInt}(a) && r.lns == Rational{BigInt}(b)) || return false
    end
    println("    M5 detail: $n rows checked, $odd with a lambda whose two readings differ")
    return n > 0 && odd > 0
end

# ---------------------------------------------------------------- M6
probe("M6", "append: the new log is old + entry, no temporary file stays; a failure at the last step leaves the log bytes unchanged and no temporary file") do
    d = joinpath(REV, "m6"); mkpath(d); log = joinpath(d, "round-log.txt")
    r1 = cr(stage(1, "m6-r1"), log)
    ok1 = read(log, String) == ents[1] && r1.action === :appended && readdir(d) == ["round-log.txt"]
    chmod(log, 0o444)      # the last step (rename over the log) is refused on a read-only log; an in-place append is refused too
    m = refusal(() -> cr(stage(2, "m6-r2"), log))
    chmod(log, 0o644)
    ok2 = m != "" && read(log, String) == ents[1] && readdir(d) == ["round-log.txt"]
    println("    M6 detail: append ok=$ok1; failure injected: \"", first(m, 80), "\"; log unchanged and folder clean=$ok2")
    return ok1 && ok2
end

# ---------------------------------------------------------------- M7
probe("M7", "a round log that is not valid UTF-8 refuses naming the log (count_round on an open log and analyze on the closed log)") do
    bad(t) = vcat(Vector{UInt8}(codeunits(t)), UInt8[0xff, 0x0a])
    l1, l3 = joinpath(REV, "m7-open.log"), joinpath(REV, "m7-closed.log")
    write(l1, bad(ents[1])); write(l3, bad(LOG3))
    m1 = refusal(() -> cr(stage(1, "m7-r1"), l1))
    m3 = refusal(() -> an(P("store"), l3, joinpath(REV, "m7-out")))
    println("    M7 detail: ", first(m1, 80), " | ", first(m3, 80))
    return has(m1, "m7-open.log", "UTF-8") && has(m3, "m7-closed.log", "UTF-8") && !isdir(joinpath(REV, "m7-out"))
end

# ---------------------------------------------------------------- M8
probe("M8", "D_B is the set of orbits with an L5 Double64 item of the arm chk: an L5 Double64 item of another arm is not counted") do
    c = ctx(P("store")); inD = [Int(o["id"]) for o in c.panel_orbits if o["in_D"] === true]
    N, d1 = c.p2.N_1, c.p2.directions[1]
    key(id) = AN.p2_key(id, "L5", "Double64", N, 1, d1)
    idx = Dict{Tuple, Any}(key(inD[1]) => (arms = ["d64"],), key(inD[2]) => (arms = ["chk"],))
    r = AN.p3t_dB(c, idx)
    println("    M8 detail: d64-arm orbit $(inD[1]), chk-arm orbit $(inD[2]); D_B = $r")
    return r == [inD[2]]
end

# ---------------------------------------------------------------- M9
probe("M9", "ladder test in exact rationals: a pair that binary64 subtraction rounds to the bound is decided by the exact difference") do
    tol = 0.2; found = nothing
    for _ in 1:200000
        a, b = 0.3 + (rand() - 0.5) * 1e-15, 0.1 + (rand() - 0.5) * 1e-15
        (abs(a - b) > tol) != (abs(Rational{BigInt}(a) - Rational{BigInt}(b)) > AN.p3_rat(tol)) && (found = (a, b); break)
    end
    found === nothing && return false
    a, b = found
    want = abs(Rational{BigInt}(a) - Rational{BigInt}(b)) > AN.p3_rat(tol)
    println("    M9 detail: pair ($a, $b): binary64 says $(abs(a - b) > tol), exact says $want")
    return AN._p3t_ladder((a, b, b), tol) == want
end

# ---------------------------------------------------------------- M10
probe("M10", "_is_p3t reads the first positional argument, not the value of an option") do
    p3, p2 = base.prereg, base.p2_prereg
    a1 = AN._is_p3t(["analyze", "--out", "outdir", p3, "m.json", "panel.toml", "q.toml"])
    a2 = AN._is_p3t(["analyze", "--store", p3, p2, "m.json", "panel.toml", "q.toml"])
    a3 = AN._is_p3t(["analyze", p3, "q.toml", "t.toml", "--threads", "4"])
    println("    M10 detail: option value first: $a1 (want true); p3 as an option value: $a2 (want false); plain: $a3 (want true)")
    return a1 === true && a2 === false && a3 === true
end

# ---------------------------------------------------------------- L11
probe("L11", "analyze for the p3 slug refuses the keywords figures and seed") do
    m1 = refusal(() -> an(P("store"), P("round-log.txt"), joinpath(REV, "l11-out-a"); figures = true))
    m2 = refusal(() -> an(P("store"), P("round-log.txt"), joinpath(REV, "l11-out-b"); seed = 5))
    println("    L11 detail: ", first(m1, 90), " | ", first(m2, 90))
    return has(m1, "figures") && has(m2, "seed") && !isdir(joinpath(REV, "l11-out-a")) && !isdir(joinpath(REV, "l11-out-b"))
end

# ---------------------------------------------------------------- L12
probe("L12", "round-down sensitivity: capped list of the round-down budgets, LE round-down budget limited to N_c") do
    mk(rd, none) = (X = Dict(:SALI => (b = nothing, capped = true, round_down = rd, round_down_none = none), :LE => (b = 640, capped = false, round_down = 640, round_down_none = false)), le_read = 100)
    qt = (voter_bits = ["SALI"], N_c = 100)
    up, ok_down, none_down = AN._p3t_capped(qt, mk(64, false), false), AN._p3t_capped(qt, mk(64, false), true), AN._p3t_capped(qt, mk(nothing, true), true)
    le = (AN.p3t_budget_fn(mk(64, false), true, 100)(:LE), AN.p3t_budget_fn(mk(64, false), true, 1000)(:LE), AN.p3t_budget_fn(mk(64, false), false, 100)(:LE))
    println("    L12 detail: capped up=$up, round-down available=$ok_down, none=$none_down; LE round-down budget $le")
    return up == [:SALI] && isempty(ok_down) && none_down == [:SALI] && le == (100, 640, 100)
end

# ---------------------------------------------------------------- L13
probe("L13", "a chunk file without the root attribute tau refuses naming it; a panel orbit without x0_dec or with a round outside 1 to r_max refuses naming the orbit and the key") do
    s = stage(1, "l13-store"); rewrite(joinpath(s, L1ref), (a, os) -> delete!(a, "tau"))
    m1 = refusal(() -> cr(s, joinpath(REV, "l13a.log")))
    txt = read(base.panel, String)
    pan(name, f) = (p = joinpath(REV, name); write(p, f(txt)); p)
    m2 = refusal(() -> AN.count_round(; base..., panel = pan("l13-x0.toml", t -> replace(t, r"x0_dec = \"[^\"]*\"\n" => ""; count = 1)), store = stage(1, "l13-s2"), round_log = joinpath(REV, "l13b.log")))
    m3 = refusal(() -> AN.count_round(; base..., panel = pan("l13-round.toml", t -> replace(t, "round = 1\n" => "round = 7\n"; count = 1)), store = stage(1, "l13-s3"), round_log = joinpath(REV, "l13c.log")))
    println("    L13 detail: ", first(m1, 80), " | ", first(m2, 80), " | ", first(m3, 80))
    return has(m1, "tau") && has(m2, "x0_dec", "orbit") && has(m3, "round", "orbit") && !isfile(joinpath(REV, "l13a.log"))
end

# ---------------------------------------------------------------- L14
probe("L14", "closed log: analyze refuses a ref chunk of a round above the last entry (count_round on a closed log stays tolerant, E62 c (i))") do
    s = cpstore(joinpath(REV, "l14-store"), chunkfiles(P("store")))
    dupchunk(s, L1ref, "chunk-00999.h5", "k9:r4")
    m = refusal(() -> an(s, P("round-log.txt"), joinpath(REV, "l14-out")))
    println("    L14 detail: ", first(m, 140))
    return has(m, "round 4") && !isdir(joinpath(REV, "l14-out"))
end

# ---------------------------------------------------------------- L16
probe("L16", "a ref chunk whose cell is not k<K index>:r<round> refuses naming the file") do
    s = stage(1, "l16-store"); dupchunk(s, L1ref, "chunk-00998.h5", "k9")
    m = refusal(() -> cr(s, joinpath(REV, "l16.log")))
    println("    L16 detail: ", first(m, 140))
    return has(m, "chunk-00998.h5")
end

# ---------------------------------------------------------------- L17
probe("L17", "T6 prints a McNemar p-value of a large exact rational with 12 significant digits (no fraction)") do
    x = (big(2)^200 + 1) // big(2)^201
    s = AN._c(AN._p3t_pshow(x))
    println("    L17 detail: ", s)
    return !occursin('/', s) && length(s) <= 20 && startswith(s, "0.5")
end

# ---------------------------------------------------------------- V1
probe("V1", "count_round with an altered entry 3 (decision=closf) and closure chunks in the store names entry 3 (before: the closure-chunk guard answered first)") do
    log = joinpath(REV, "v1.log"); write(log, replace(LOG3, "decision=close" => "decision=closf"; count = 1))
    m = refusal(() -> cr(P("store"), log))
    println("    V1 detail: ", first(m, 140))
    return has(m, "round log entry 3")
end

println("REV summary | $(nresult[]) probes, $(nresult[] - nfail[]) PASS, $(nfail[]) FAIL | ", nfail[] == 0 ? "PASS" : "FAIL")
