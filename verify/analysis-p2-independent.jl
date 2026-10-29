# Independent reader and store-identity checker for study p2-indicators (oracle O6 reader, O9, O10 of the analysis spec, Amendment D89).
# Uses HDF5, SHA, TOML and the standard library only: no Chaoslab code is loaded and no src/Analysis file was read.
# Written from analysis-spec.md (E23 to E26, E29, E36, O9, O10), schema/run-store.md and the preregistration.
#
# Usage (from the repository root):
#   julia --project=verify -t 8 verify/analysis-p2-independent.jl [export file]
# Inputs (fixed paths, read-only): data/p2-indicators, data/p1-le-benchmark, the P2 manifest, the panel file, the P2 parameter table,
# the P1 parameter table and the P2 preregistration (for the [[run.h4]] tables).
# Output: a JSON-lines export (default C:\...\scratchpad\p2ind\p2-export.jsonl, argument 1), read by verify/analysis-p2-decision.py.
#
# Export format. One JSON object per line; NaN, Infinity and -Infinity are written as such (Python json reads them).
#   {"k":"orbit", "id", "src":"P"|"G", "cls" (panel class, "h4" for the grid), "control":bool, "ckind", "K":float|null, "kidx":int,
#    "i", "j" (grid only), "in_SD", "in_SL", "ic", "params", "flags":[orbit flag names, union over all items],
#    "drift_recomputed":bool, "inst", "struct" (P1 labels, panel only, final after all P1 tests), "forced":[...],
#    "exponent" (controls), "cdev", "cscope" (stored control_deviation and control_scope of the L1 item at N_1, controls only)}
#   {"k":"s", "o", "spec", "prec", "hor", "dir", "t":[recorded steps], "c":{column name: values}}      series item
#       spec L1: c = {"lam1"}; L5 (grid, N_1): c = {"lam1","xlift"} (xlift = wrapped x + 2 pi winding, binary64);
#       tg: c = {"FLI.fli","MEGNO.Ybar20_1","SALI.sali","SALI.L"}; pr: {"RLI.RLI"}; st: {"stretching.D"}
#   {"k":"w", "o", "spec", "prec", "hor", "rows":[labels], "c":{column: one value per row}, "flags":[int per row], "bits":[names],
#    "sc":{scalar name: value}}                                                                     window item (fq1, fq1_raw, fq1_qp, pa1, pa_planar, bk1)
#       fq1 and pa1 columns D, nu1; bk1 column dig; (all window specs keep every column named D, nu1, dig that exist)
#   {"k":"chunkflags", "ordinal", "spec", "prec", "recomputed":{flag: orbit count}, "manifest":{flag: count}}
#   {"k":"na", "spec", "recomputed":[orbit ids omitted], "manifest":[...]}
#   {"k":"meta", ...}  counts and the O9 results.
# Direction text of an item is the stored attribute `direction` ("" if absent). Items are found by content, never by chunk number.
using HDF5, SHA, TOML

const T0 = time()
const ROOT = normpath(joinpath(@__DIR__, ".."))
const OUTF = length(ARGS) >= 1 ? ARGS[1] : "C:/Users/tesse/AppData/Local/Temp/claude/C--Users-tesse-Documents-Vibe-Coding-Chaos/5df15ee1-47df-4875-85bd-1bd4d7c665a6/scratchpad/p2ind/p2-export.jsonl"
const P2STORE = joinpath(ROOT, "data/p2-indicators")
const P1STORE = joinpath(ROOT, "data/p1-le-benchmark")
const MANIF = joinpath(ROOT, "workspaces/study/stages/03-run/output/p2-indicators-manifest.json")
const PANELF = joinpath(ROOT, "workspaces/study/stages/02-panel/output/p1-le-benchmark-panel.toml")
const PARF = joinpath(ROOT, "workspaces/study/stages/04-analyze/references/p2-indicators-analysis-parameters.toml")
const P1PARF = joinpath(ROOT, "workspaces/study/stages/04-analyze/references/p1-le-benchmark-analysis-parameters.toml")
const PREREG = joinpath(ROOT, "workspaces/study/stages/01-design/output/p2-indicators-preregistration.md")
const P1TABLES = joinpath(ROOT, "workspaces/study/stages/04-analyze/output/p1-le-benchmark-tables.md")
const P1RESULTS = joinpath(ROOT, "workspaces/module-build/stages/04-verify/output/analysis-independent-results.txt")

say(msg...) = println("[", round(time() - T0; digits = 1), " s] ", msg...)
const ANOM = String[]
anom(msg) = (push!(ANOM, msg); say("ANOMALY: ", msg))
const CHECKS = Tuple{String,Bool,String}[]
check(name, ok, detail = "") = (push!(CHECKS, (name, ok, detail)); say(ok ? "PASS " : "FAIL ", name, " ", detail))

# ------------------------------------------------------------------ parameters (the table is the only source of study numbers)
const PAR = TOML.parsefile(PARF)
pv(k) = PAR[k]["value"]
const N1, NC, NDIR = Int(pv("N_1")), Int(pv("N_c")), Int(pv("N_dir"))
const DIRS = String.(pv("directions"))
const WINDOWS = Int.(pv("windows"))
const BUDGETS = Int.(pv("budgets"))
const GRIDK = String.(pv("grid_k_values"))
const NS = Int(pv("n_s"))
const DRIFT = Float64(pv("drift_bound"))
const P1PAR = TOML.parsefile(P1PARF)

const PANEL = TOML.parsefile(PANELF)
const ORB = sort(PANEL["orbit"]; by = o -> o["id"])
const PORB = Dict(o["id"] => o for o in ORB)
const MAXPANEL = maximum(o["id"] for o in ORB)

# ------------------------------------------------------------------ json writer
jnum(x::Real) = x isa Integer ? string(x) : (isnan(x) ? "NaN" : (isinf(x) ? (x > 0 ? "Infinity" : "-Infinity") : repr(Float64(x))))
jstr(s::AbstractString) = "\"" * replace(String(s), "\\" => "\\\\", "\"" => "\\\"") * "\""
jvec(v) = "[" * join((x isa AbstractString ? jstr(x) : jnum(x) for x in v), ",") * "]"
jdict(d::AbstractDict) = "{" * join((jstr(k) * ":" * (v isa AbstractVector ? jvec(v) : (v isa AbstractString ? jstr(v) : (v isa Bool ? string(v) : (v === nothing ? "null" : jnum(v))))) for (k, v) in sort(collect(d); by = first)), ",") * "}"

# ================================================================== PART 1: P1 labels of the panel orbits (copied from verify/analysis-independent.jl, trimmed)
const POS, NONE, UNR = "positive growth supported", "none detected on the tested horizon", "unresolved"
struct Item1
    t::Vector{Int64}; lam::Vector{Float64}; st::Matrix{Float64}; flags::Vector{String}
end
tdist(a, b) = (d = a - b; d - 2pi * round(d / (2pi)))
function idx1(it, n::Integer)
    i = searchsortedfirst(it.t, n)
    return (i <= length(it.t) && it.t[i] == n) ? i : 0
end
lamat1(it, n) = (it === nothing ? NaN : (i = idx1(it, n); i == 0 ? NaN : it.lam[i]))

function labelof(P, it, N::Int, o; hfix::Int = 0)
    it === nothing && return (UNR, UNR)
    iN = idx1(it, N); i10 = idx1(it, N ÷ 10)
    (iN == 0 || i10 == 0) && return (UNR, UNR)
    lam = it.lam[iN]; lam10 = it.lam[i10]
    inst = UNR
    if isfinite(lam) && isfinite(lam10) && lam10 > 0
        r = 10 * lam / lam10
        if lam >= P["lambda_threshold"] && r >= P["growth_ratio"]
            inst = POS
        elseif lam < P["lambda_threshold"] && r < P["growth_ratio"]
            inst = NONE
        end
    end
    x0, p0 = Float64(o["x0"]), Float64(o["p0"])
    stc = UNR
    if get(o, "control", false) && get(o, "control_kind", "none") in ("fixed", "newton")
        iq = idx1(it, Int(o["period"]))
        if iq != 0
            res = hypot(tdist(it.st[iq, 1], x0), tdist(it.st[iq, 2], p0))
            res <= P["return_tol"] && (stc = "periodic")
        end
    end
    if stc == UNR
        if inst == POS
            stc = "irregular"
        elseif inst == NONE
            best = 0; bd = Inf
            hfix > 0 && (best = idx1(it, hfix))
            hfix == 0 && for j in 1:(iN - 1)
                d = abs(it.t[j] - N / 2)
                if d < bd; bd = d; best = j; end
            end
            if best > 0
                h = Int(it.t[best])
                rho1 = (it.st[best, 1] - x0) / (2pi * h)
                rho2 = (it.st[iN, 1] - it.st[best, 1]) / (2pi * (N - h))
                abs(rho1 - rho2) <= P["rotation_tol"] && (stc = "quasiperiodic")
            end
        end
    end
    return (inst, stc)
end

function p1_labels()
    P = P1PAR
    n1, nlong, nb = Int(P["N_1"]), Int(P["N_long"]), Int(P["N_B"])
    VM = ["L1", "L3", "L4", "L4b", "L5", "L5b", "L6", "L9"]
    BM = filter(!=("L9"), VM)
    ITEMS = Dict{Tuple{Int,String,String,Int,String,String},Item1}()
    ORBFLAGS = Dict{Int,Set{String}}()
    files = sort(filter(fn -> occursin(r"^chunk-\d{5}\.h5$", fn), readdir(P1STORE)))
    say("P1 labels: reading ", length(files), " P1 chunk files")
    for fn in files
        h5open(joinpath(P1STORE, fn), "r") do f
            a = Dict(nm => HDF5.read_attribute(f, nm) for nm in keys(HDF5.attributes(f)))
            m = a["method"]
            prec = startswith(a["precision"], "BigFloat") ? "BigFloat" : a["precision"]
            hor = Int(a["horizon"]); tau = Int(a["tau"]); dir = a["direction"]; mode = a["coordinate_mode"]
            ok = a["status"] == "ok"
            keep = ok && tau == 1 && !(mode == "wrapped" && m != "L2") && m != "L2"
            for on in keys(f["orbits"])
                g = f["orbits"][on]; id = parse(Int, on)
                fl = String.(HDF5.read_attribute(g, "flags"))
                union!(get!(ORBFLAGS, id, Set{String}()), fl)
                keep || continue
                t = read(g["t"]); L = read(g["lyap"][m])
                lam = ndims(L) == 1 ? collect(L) : collect(L[:, 1])
                st = m == "L5" ? read(g["state"]) : zeros(0, 2)
                key = (id, m, prec, hor, dir, mode)
                haskey(ITEMS, key) || (ITEMS[key] = Item1(t, lam, st, fl))
            end
        end
    end
    getit(id, m, prec, hor; dir = "pi/4", mode = "lifted") = get(ITEMS, (id, m, prec, hor, dir, mode), nothing)
    out = Dict{Int,NamedTuple}()
    for o in ORB
        id = o["id"]
        prim = Dict(m => getit(id, m, "Float64", n1) for m in VM)
        l5long = getit(id, "L5", "Float64", nlong)
        Nlab = l5long !== nothing ? nlong : n1
        litem = l5long !== nothing ? l5long : prim["L5"]
        rinst, rstruct = labelof(P, litem, Nlab, o)
        fl = sort(collect(get(ORBFLAGS, id, Set{String}())))
        forced = String[]
        if get(o, "in_SL", false)
            f64 = prim["L5"]; d64 = getit(id, "L5", "Double64", n1); bf = getit(id, "L5", "BigFloat", nb)
            if f64 !== nothing && d64 !== nothing && bf !== nothing
                lf, ld, lb = lamat1(f64, nb), lamat1(d64, nb), lamat1(bf, nb)
                tc = intersect(f64.t, d64.t, bf.t); tc = tc[(tc .> 0) .& (tc .< nb)]
                hc = tc[argmin(abs.(tc .- nb / 2))]
                pf, pd, pb = labelof(P, f64, nb, o; hfix = Int(hc)), labelof(P, d64, nb, o; hfix = Int(hc)), labelof(P, bf, nb, o; hfix = Int(hc))
                (abs(lf - ld) > P["ladder_tol"] || abs(ld - lb) > P["ladder_tol"] || !(pf == pd == pb)) && push!(forced, "ladder")
            else
                anom("P1 orbit $id (S_L): ladder item missing")
            end
        end
        vals_B = [lamat1(prim[m], n1) for m in BM]
        dqb = maximum(vals_B) - minimum(vals_B)
        d9 = abs(lamat1(prim["L9"], n1) - lamat1(prim["L5"], n1))
        (dqb > P["h1a_bound"] || d9 > P["h1b_bound"]) && push!(forced, "estimator")
        if get(o, "in_SD", false)
            labs = Set{Tuple{String,String}}(); nd = 0
            for ((kid, km, kp, kh, kd, kmode), it) in ITEMS
                (kid == id && km == "L5" && kp == "Float64" && kh == n1 && kmode == "lifted") || continue
                nd += 1; push!(labs, labelof(P, it, n1, o))
            end
            (nd >= 2 && length(labs) > 1) && push!(forced, "direction")
        end
        prim["L5"] !== nothing && ("diagnostics exceeded" in prim["L5"].flags) && push!(forced, "diagnostics")
        "excluded" in fl && push!(forced, "excluded")
        order = ["ladder", "estimator", "direction", "diagnostics", "excluded"]
        forced = filter(x -> x in forced, order)
        inst, stc = isempty(forced) ? (rinst, rstruct) : (UNR, UNR)
        out[id] = (inst = inst, stc = stc, forced = forced, nlab = Nlab)
    end
    return out
end

# ================================================================== PART 2: the P2 store
struct SItem
    t::Vector{Int64}; c::Dict{String,Vector{Float64}}
end
struct WItem
    rows::Vector{String}; c::Dict{String,Vector{Float64}}; flags::Vector{Int64}; bits::Vector{String}; sc::Dict{String,Float64}
end
mtxt = read(MANIF, String)
# tiny JSON parser for the chunk objects (strings, numbers, bool, null, arrays, objects)
function jparse(s::String)
    pos = Ref(1)
    ws() = (while pos[] <= lastindex(s) && isspace(s[pos[]]); pos[] = nextind(s, pos[]); end)
    function val()
        ws(); c = s[pos[]]
        if c == '{'
            d = Dict{String,Any}(); pos[] += 1; ws()
            if s[pos[]] == '}'; pos[] += 1; return d; end
            while true
                ws(); k = str(); ws(); pos[] += 1  # colon
                d[k] = val(); ws()
                if s[pos[]] == ','; pos[] += 1; else pos[] += 1; break; end
            end
            return d
        elseif c == '['
            a = Any[]; pos[] += 1; ws()
            if s[pos[]] == ']'; pos[] += 1; return a; end
            while true
                push!(a, val()); ws()
                if s[pos[]] == ','; pos[] += 1; else pos[] += 1; break; end
            end
            return a
        elseif c == '"'; return str()
        elseif startswith(SubString(s, pos[]), "true"); pos[] += 4; return true
        elseif startswith(SubString(s, pos[]), "false"); pos[] += 5; return false
        elseif startswith(SubString(s, pos[]), "null"); pos[] += 4; return nothing
        else
            m = match(r"^-?[0-9.eE+-]+", SubString(s, pos[]))
            pos[] += length(m.match)
            x = parse(Float64, m.match)
            return (occursin(r"^-?\d+$", m.match)) ? Int(x) : x
        end
    end
    function str()
        pos[] += 1; buf = IOBuffer()
        while s[pos[]] != '"'
            if s[pos[]] == '\\'
                pos[] += 1; c = s[pos[]]
                write(buf, c == 'n' ? '\n' : c == 't' ? '\t' : c)
            else
                write(buf, s[pos[]])
            end
            pos[] = nextind(s, pos[])
        end
        pos[] += 1
        return String(take!(buf))
    end
    return val()
end
const MTOP = jparse(mtxt)  # whole manifest (800 kB), fine
const MCH = Dict{Int,Dict{String,Any}}(Int(d["ordinal"]) => d for d in MTOP["chunks"])
say("manifest: study ", MTOP["study"], ", chunks ", length(MCH), ", status ", MTOP["status"])

function chunk_attrs(f)
    Dict(nm => HDF5.read_attribute(f, nm) for nm in keys(HDF5.attributes(f)))
end

const SERIES_COLS = Dict("L1" => ["lam1"], "tg" => ["FLI.fli", "MEGNO.Ybar20_1", "SALI.sali", "SALI.L"], "pr" => ["RLI.RLI"], "st" => ["stretching.D"])
const SI = Dict{Tuple{Int,String,String,Int,String},SItem}()          # (orbit, spec, prec, horizon, direction)
const WI = Dict{Tuple{Int,String,String,Int},WItem}()
const ORBINFO = Dict{Int,Dict{String,Any}}()
const FULLTG = Dict{Tuple{Int,Int},Matrix{Float64}}()                  # S_D prefix check: (orbit, horizon) -> full tangent chain, direction pi/4 Float64
const FULLL1 = Dict{Tuple{Int,Int},Vector{Float64}}()
const DRIFTMAX = Dict{Int,Float64}()
const CFLAGS = Any[]
const PRESENT = Dict{Tuple{String,String},Set{Int}}()                  # (spec, prec) -> orbit ids with an item (for the omission check)
const EXPECT_CHUNK_ATTR = Dict{String,Any}()

function read_p2()
    files = sort(filter(fn -> occursin(r"^chunk-\d{5}\.h5$", fn), readdir(P2STORE)))
    say("reading ", length(files), " P2 chunk files")
    sd_ids = Set(o["id"] for o in ORB if get(o, "in_SD", false))
    nfile = 0
    prereg = MTOP["prereg"]["header"]
    for fn in files
        nfile += 1
        ordn = parse(Int, fn[7:11])
        h5open(joinpath(P2STORE, fn), "r") do f
            a = chunk_attrs(f)
            me = get(MCH, ordn, nothing)
            me === nothing && anom("$fn not in manifest")
            a["study"] == "p2-indicators" || anom("$fn study")
            get(a, "schema_version", 0) == 2 || anom("$fn schema_version")
            Bool(a["complete"]) || anom("$fn incomplete")
            Bool(a["pilot"]) && anom("$fn pilot")
            a["status"] == "ok" || anom("$fn status $(a["status"])")
            a["prereg_sha256"] == prereg || anom("$fn prereg_sha256")
            a["plan_sha256"] == MTOP["plan_sha256"] || anom("$fn plan_sha256")
            a["panel_sha256"] == MTOP["panel"]["sha256_start"] || anom("$fn panel_sha256")
            me !== nothing && a["chunk_hash"] != me["hash"] && anom("$fn chunk_hash differs from manifest")
            me !== nothing && me["status"] != "ok" && anom("$fn manifest status $(me["status"])")
            m = a["method"]; prec = a["precision"]
            hor = Int(a["horizon"]); dir = get(a, "direction", ""); mode = a["coordinate_mode"]
            mode == "wrapped" || anom("$fn mode $mode")
            kind = get(a, "kind", "")
            Int(a["tau"]) == 1 || anom("$fn tau")
            orbs = keys(f["orbits"])
            me !== nothing && length(orbs) != me["orbits"] && anom("$fn orbit count $(length(orbs)) vs manifest $(me["orbits"])")
            fcount = Dict{String,Int}()
            bitnames = String[]
            for on in orbs
                g = f["orbits"][on]; id = parse(Int, on)
                fl = String.(HDF5.read_attribute(g, "flags"))
                info = get!(ORBINFO, id, Dict{String,Any}("flags" => Set{String}(), "params" => String(HDF5.read_attribute(g, "params")),
                                                          "ic" => String(HDF5.read_attribute(g, "ic")), "dc" => String(HDF5.read_attribute(g, "design_class"))))
                union!(info["flags"], fl)
                info["params"] == String(HDF5.read_attribute(g, "params")) || anom("orbit $id params differ between chunks")
                info["ic"] == String(HDF5.read_attribute(g, "ic")) || anom("orbit $id ic differs between chunks")
                push!(get!(PRESENT, (m, prec), Set{Int}()), id)
                t = read(g["t"])
                if kind == ""                                      # catalog method: L1, L5
                    L = read(g["lyap"][m]); lam = ndims(L) == 1 ? collect(L) : collect(L[:, 1])
                    if m == "L1"
                        key = (id, "L1", prec, hor, dir)
                        haskey(SI, key) && anom("duplicate item $key")
                        SI[key] = SItem(collect(t), Dict("lam1" => lam))
                        if prec == "Float64" && hor == N1 && dir == DIRS[1]
                            S = read(g["state"]); W = read(g["winding"]); DRIFTMAX[id] = maximum(abs.(S[:, 2] .+ 2pi .* Float64.(W[:, 2])))   # lifted p = wrapped p + 2 pi winding
                            if haskey(HDF5.attributes(g), "control_deviation")
                                info["cdev"] = Float64(HDF5.read_attribute(g, "control_deviation"))
                                info["cscope"] = String(HDF5.read_attribute(g, "control_scope"))
                            end
                        end
                        id in sd_ids && prec == "Float64" && dir == DIRS[1] && (FULLL1[(id, hor)] = lam)
                    elseif m == "L5"
                        S = read(g["state"]); W = read(g["winding"])
                        key = (id, "L5", prec, hor, dir)
                        SI[key] = SItem(collect(t), Dict("lam1" => lam, "xlift" => S[:, 1] .+ 2pi .* Float64.(W[:, 1])))
                    else
                        anom("$fn unexpected method $m")
                    end
                elseif kind in ("tangent", "pair", "stretching")
                    ds = g["indicators"][m]["series"]
                    cols = String.(HDF5.read_attribute(ds, "columns")); M = read(ds)
                    if m == "st" && dir != DIRS[1]
                        continue
                    end
                    c = Dict{String,Vector{Float64}}()
                    for cn in SERIES_COLS[m]
                        j = findfirst(==(cn), cols); j === nothing && error("$fn column $cn")
                        c[cn] = collect(M[:, j])
                    end
                    SI[(id, m, prec, hor, dir)] = SItem(collect(t), c)
                    if m == "tg" && id in sd_ids && prec == "Float64" && dir == DIRS[1]
                        FULLTG[(id, hor)] = Matrix{Float64}(M)
                    end
                elseif kind in ("frequency", "phase_advance", "birkhoff")
                    grp = g["indicators"][m]
                    ds = grp["final"]; cols = String.(HDF5.read_attribute(ds, "columns")); rows = String.(HDF5.read_attribute(ds, "rows"))
                    M = read(ds)
                    c = Dict{String,Vector{Float64}}()
                    for cn in ("D", "nu1", "dig")
                        j = findfirst(==(cn), cols); j === nothing || (c[cn] = collect(M[:, j]))
                    end
                    fds = grp["flags"]; fw = collect(Int64, read(fds)); bitnames = String.(HDF5.read_attribute(fds, "bits"))
                    sc = Dict{String,Float64}()
                    if haskey(grp, "scalars")
                        sds = grp["scalars"]; scn = String.(HDF5.read_attribute(sds, "columns")); sv = read(sds)
                        for (j, nme) in enumerate(scn); sc[nme] = sv[j]; end
                    end
                    key = (id, m, prec, hor)
                    haskey(WI, key) && anom("duplicate window item $key")
                    WI[key] = WItem(rows, c, fw, bitnames, sc)
                    for (j, b) in enumerate(bitnames)
                        any(w -> (w >> (j - 1)) & 1 == 1, fw) && (fcount[b] = get(fcount, b, 0) + 1)
                    end
                else
                    anom("$fn unknown kind $kind")
                end
            end
            if kind in ("frequency", "phase_advance", "birkhoff")
                man = me === nothing ? Dict{String,Any}() : get(me, "indicator_flags", Dict{String,Any}())
                push!(CFLAGS, (ordn, m, prec, fcount, Dict{String,Int}(k => Int(v) for (k, v) in man)))
            end
        end
        nfile % 100 == 0 && say("  files read ", nfile)
    end
    return nfile
end

# ------------------------------------------------------------------ main
labels = p1_labels()
say("P1 labels done: ", length(labels))
nfile = read_p2()
say("P2 items: series ", length(SI), ", windows ", length(WI), " (total ", length(SI) + length(WI), "), orbits ", length(ORBINFO))

# ---- O10: panel label counts per cell against T6 of the P1 tables
let
    cnt = Dict{Tuple{Int,String,String,String},Int}()
    for o in ORB
        l = labels[o["id"]]
        key = get(o, "control", false) ? (0, "controls", l.inst, l.stc) : (Int(o["k_index"]), String(o["class"]), l.inst, l.stc)
        cnt[key] = get(cnt, key, 0) + 1
    end
    t6 = Dict{Tuple{Int,String,String,String},Int}()
    txt = read(P1TABLES, String)
    i0 = findfirst("## T6 Label table", txt)
    sec = txt[first(i0):end]
    stop = findnext("\n## ", sec, 5); sec = stop === nothing ? sec : sec[1:first(stop)]
    for ln in split(sec, '\n')
        c = strip.(split(ln, '|'))
        length(c) >= 7 && occursin(r"^\d+$", c[2]) && (t6[(parse(Int, c[2]), String(c[3]), String(c[4]), String(c[5]))] = parse(Int, c[6]))
    end
    check("O10 panel label counts per (K index, class, instability, structure) equal T6 of the P1 tables", cnt == t6, "cells $(length(cnt)) vs T6 $(length(t6)), orbits $(sum(values(cnt))) vs $(sum(values(t6)))")
    if isfile(P1RESULTS)
        old = Dict{Int,Tuple{String,String}}()
        for ln in eachline(P1RESULTS)
            startswith(ln, "ORBIT|") || continue
            c = split(ln, '|'); old[parse(Int, c[2])] = (String(c[15]), String(c[16]))
        end
        diffs = [id for id in keys(labels) if haskey(old, id) && old[id] != (labels[id].inst, labels[id].stc)]
        check("O10 final labels equal those of the D60 independent recomputation (analysis-independent-results.txt)", isempty(diffs) && length(old) == length(labels), "listed $(length(old)), differing $(length(diffs))")
    end
end

# ---- E24 identity
let
    ktext = Dict(i => GRIDK[i] for i in eachindex(GRIDK))
    # parse the [[run.h4]] tables of the preregistration
    txt = read(PREREG, String)
    i0 = findfirst("\n[[run.h4]]\n", txt)
    body = txt[(first(i0) + 1):end]
    j = findfirst("\n```", body)
    tabs = TOML.parse(body[1:first(j)])["run"]["h4"]
    check("E24 [[run.h4]] table count equals grid_orbit_count", length(tabs) == Int(pv("grid_orbit_count")), "$(length(tabs))")
    nbad = 0
    GRIDID = Dict{Int,Tuple{Int,Int,Int}}()
    for (n, tb) in enumerate(tabs)
        id = MAXPANEL + n
        m = match(r"^grid-k(\d+)-i(\d+)-j(\d+)$", tb["label"]); m === nothing && (nbad += 1; continue)
        ki, ii, jj = parse(Int, m[1]), parse(Int, m[2]), parse(Int, m[3])
        kk = findfirst(==(tb["K"]), GRIDK)
        info = get(ORBINFO, id, nothing)
        if kk != ki || info === nothing || parse(Float64, info["params"]) != parse(Float64, tb["K"]) || info["ic"] != tb["x0"] * "," * tb["p0"]
            nbad += 1; anom("E24 grid orbit $id mismatch"); continue
        end
        GRIDID[id] = (ki, ii, jj)
    end
    check("E24 all grid orbits: stored params and ic equal the table K, x0, p0; K index and label recovered", nbad == 0 && length(GRIDID) == 2304, "bad $nbad")
    nbp = 0
    for o in ORB
        info = get(ORBINFO, o["id"], nothing)
        if info === nothing || info["ic"] != o["x0_dec"] * "," * o["p0_dec"]; nbp += 1; continue; end
        if haskey(o, "K") && parse(Float64, info["params"]) != Float64(o["K"]); nbp += 1; end
    end
    check("E24 all panel orbits: stored ic equals panel decimal texts, params equals panel K", nbp == 0, "bad $nbp")
    global GRIDIDS = GRIDID
end

# ---- O9 (a) prefix consistency for S_D: horizon 1e6 items equal horizon 1e7 items at every common recorded step, bit for bit
let
    nvals = 0; ndiff = 0; norb = 0
    for o in ORB
        get(o, "in_SD", false) || continue
        id = o["id"]; norb += 1
        a = SI[(id, "L1", "Float64", N1, DIRS[1])]; b = SI[(id, "L1", "Float64", NDIR, DIRS[1])]
        ia = Dict(t => i for (i, t) in enumerate(a.t))
        for (j, t) in enumerate(b.t)
            haskey(ia, t) || continue
            nvals += 1
            reinterpret(UInt64, a.c["lam1"][ia[t]]) == reinterpret(UInt64, b.c["lam1"][j]) || (ndiff += 1)
        end
        A = FULLTG[(id, N1)]; B = FULLTG[(id, NDIR)]
        ta = SI[(id, "tg", "Float64", N1, DIRS[1])].t; tb = SI[(id, "tg", "Float64", NDIR, DIRS[1])].t
        ia = Dict(t => i for (i, t) in enumerate(ta))
        for (j, t) in enumerate(tb)
            haskey(ia, t) || continue
            for c in 1:size(A, 2)
                nvals += 1
                reinterpret(UInt64, A[ia[t], c]) == reinterpret(UInt64, B[j, c]) || (ndiff += 1)
            end
        end
    end
    check("O9a prefix consistency (L1 and tangent chain, horizon 1e6 against 1e7) for S_D orbits", ndiff == 0 && norb == 172, "orbits $norb, values compared $nvals, differences $ndiff")
end

# ---- O9 (b) flag counts per chunk against the manifest
let
    nbad = 0; nchunks = 0
    for (ordn, m, prec, fc, man) in CFLAGS
        nchunks += 1
        a = Dict(k => v for (k, v) in fc if v > 0); b = Dict(k => v for (k, v) in man if v > 0)
        if a != b; nbad += 1; anom("flag counts chunk $ordn $m: recomputed $a manifest $b"); end
    end
    nlist = count(c -> !isempty(c[5]), CFLAGS)
    check("O9b flag counts recomputed from the stored flag words equal the manifest indicator_flags per chunk", nbad == 0, "window chunks $nchunks, with listed flags $nlist, differing $nbad")
end

# ---- O9 (c) omissions (not_applicable) and item presence
let
    manifest_na = Dict{String,Vector{Int}}()
    for (ordn, me) in MCH
        if haskey(me, "not_applicable") && me["not_applicable"] isa Dict
            for (sp, ids) in me["not_applicable"]; append!(get!(manifest_na, sp, Int[]), Int.(ids)); end
        end
    end
    allp = Set(1:3187); panelnoy5 = Set(1:882)
    ok = true; det = String[]
    for (sp, expectset) in (("fq1", [allp, "Float64"]), ("bk1", [allp, "Float64"]))
        pres = get(PRESENT, (sp, "Float64"), Set{Int}())
        omitted = sort(collect(setdiff(allp, pres)))
        man = sort(unique(get(manifest_na, sp, Int[])))
        push!(det, "$sp omitted $(omitted) manifest $(man)")
        omitted == man || (ok = false)
    end
    for sp in ("pa1", "pr", "L1", "tg", "st")
        pres = get(PRESENT, (sp, "Float64"), Set{Int}())
        # pr only exists for S_L and G
        want = sp == "pr" ? Set(vcat([o["id"] for o in ORB if get(o, "in_SL", false)], collect(884:3187))) : allp
        pres == want || (ok = false; push!(det, "$sp presence differs: missing $(length(setdiff(want, pres))) extra $(length(setdiff(pres, want)))"))
    end
    check("O9c omissions recomputed from the store equal the manifest not_applicable; E25 item presence for orbit sets", ok, join(det, "; "))
    global MANIFEST_NA = manifest_na
end

# ---- O9 (d) control deviation, drift flag
let
    nbad = 0; ncontrols = 0; maxdiff = 0.0
    for o in ORB
        get(o, "control", false) || continue
        id = o["id"]; info = ORBINFO[id]
        haskey(info, "cdev") || continue
        ncontrols += 1
        it = SI[(id, "L1", "Float64", N1, DIRS[1])]
        i = findfirst(==(N1), it.t)
        dev = abs(it.c["lam1"][i] - Float64(o["exponent"]))
        d = abs(dev - info["cdev"]); maxdiff = max(maxdiff, d)
        d > 1e-12 && (nbad += 1; anom("control $id: stored control_deviation $(info["cdev"]) recomputed $dev"))
    end
    check("O9d stored control_deviation equals its recomputation to 1e-12", nbad == 0 && ncontrols > 0, "controls with the attribute $ncontrols, largest difference $maxdiff")
    ndr = 0
    for o in ORB
        id = o["id"]
        rec = DRIFTMAX[id] > DRIFT; st = "drift" in ORBINFO[id]["flags"]
        rec == st || (ndr += 1; anom("drift flag of orbit $id: stored $st recomputed $rec"))
    end
    for id in 884:3187
        (DRIFTMAX[id] > DRIFT) == ("drift" in ORBINFO[id]["flags"]) || (ndr += 1)
    end
    check("E36k drift flags equal the recomputation from the stored points (max |p| above drift_bound)", ndr == 0, "flagged orbits " * string(sort([id for id in keys(ORBINFO) if "drift" in ORBINFO[id]["flags"]])))
end

# ---- non-finite counts in used series (E36 j, informational)
let
    nonfin = Dict{String,Int}()
    for ((id, sp, prec, hor, dir), it) in SI
        for (cn, v) in it.c
            n = count(!isfinite, v); n > 0 && (nonfin["$sp/$prec/$cn"] = get(nonfin, "$sp/$prec/$cn", 0) + n)
        end
    end
    say("non-finite values in exported series: ", nonfin)
    nw = Dict{String,Int}()
    for ((id, sp, prec, hor), it) in WI
        for (cn, v) in it.c
            for (r, x) in enumerate(v)
                isfinite(x) || (nw["$sp/$prec/$cn"] = get(nw, "$sp/$prec/$cn", 0) + 1)
            end
        end
        if haskey(it.sc, "complete") && it.sc["complete"] != 1.0; anom("window item $id $sp $prec complete != 1"); end
    end
    say("non-finite values in window tables (with or without flag): ", nw)
end

# ------------------------------------------------------------------ export
say("writing export ", OUTF)
open(OUTF, "w") do io
    for o in ORB
        id = o["id"]; info = ORBINFO[id]; l = labels[id]
        d = Dict{String,Any}("k" => "orbit", "id" => id, "src" => "P", "cls" => String(o["class"]), "control" => Bool(get(o, "control", false)),
                             "ckind" => String(get(o, "control_kind", "none")), "kidx" => Int(o["k_index"]), "in_SD" => Bool(get(o, "in_SD", false)),
                             "in_SL" => Bool(get(o, "in_SL", false)), "ic" => info["ic"], "params" => info["params"], "flags" => sort(collect(info["flags"])),
                             "drift_recomputed" => DRIFTMAX[id] > DRIFT, "inst" => l.inst, "struct" => l.stc, "forced" => l.forced, "nlab" => l.nlab,
                             "K" => haskey(o, "K") ? Float64(o["K"]) : nothing)
        haskey(o, "exponent") && (d["exponent"] = Float64(o["exponent"]))
        haskey(info, "cdev") && (d["cdev"] = info["cdev"]; d["cscope"] = info["cscope"])
        println(io, jdict(d))
    end
    for id in (MAXPANEL + 1):3187
        info = ORBINFO[id]; (ki, ii, jj) = GRIDIDS[id]
        d = Dict{String,Any}("k" => "orbit", "id" => id, "src" => "G", "cls" => "h4", "control" => false, "ckind" => "none", "kidx" => ki, "i" => ii, "j" => jj,
                             "in_SD" => false, "in_SL" => false, "ic" => info["ic"], "params" => info["params"], "flags" => sort(collect(info["flags"])),
                             "drift_recomputed" => DRIFTMAX[id] > DRIFT, "K" => parse(Float64, info["params"]))
        println(io, jdict(d))
    end
    for key in sort(collect(keys(SI)))
        (id, sp, prec, hor, dir) = key; it = SI[key]
        println(io, "{\"k\":\"s\",\"o\":", id, ",\"spec\":", jstr(sp), ",\"prec\":", jstr(prec), ",\"hor\":", hor, ",\"dir\":", jstr(dir),
                ",\"t\":", jvec(it.t), ",\"c\":", jdict(it.c), "}")
    end
    for key in sort(collect(keys(WI)))
        (id, sp, prec, hor) = key; it = WI[key]
        println(io, "{\"k\":\"w\",\"o\":", id, ",\"spec\":", jstr(sp), ",\"prec\":", jstr(prec), ",\"hor\":", hor, ",\"rows\":", jvec(it.rows),
                ",\"c\":", jdict(it.c), ",\"flags\":", jvec(it.flags), ",\"bits\":", jvec(it.bits), ",\"sc\":", jdict(it.sc), "}")
    end
    for (ordn, m, prec, fc, man) in sort(CFLAGS; by = first)
        println(io, "{\"k\":\"chunkflags\",\"ordinal\":", ordn, ",\"spec\":", jstr(m), ",\"prec\":", jstr(prec), ",\"recomputed\":", jdict(fc), ",\"manifest\":", jdict(man), "}")
    end
    for sp in sort(collect(keys(MANIFEST_NA)))
        om = sort(collect(setdiff(Set(1:3187), get(PRESENT, (sp, "Float64"), Set{Int}()))))
        println(io, "{\"k\":\"na\",\"spec\":", jstr(sp), ",\"recomputed\":", jvec(om), ",\"manifest\":", jvec(sort(unique(MANIFEST_NA[sp]))), "}")
    end
    for (nm, ok, det) in CHECKS
        println(io, "{\"k\":\"meta\",\"check\":", jstr(nm), ",\"ok\":", ok, ",\"detail\":", jstr(det), "}")
    end
    println(io, "{\"k\":\"meta\",\"files\":", nfile, ",\"series_items\":", length(SI), ",\"window_items\":", length(WI), ",\"anomalies\":", length(ANOM), "}")
end
say("export written; checks ", count(c -> c[2], CHECKS), " pass of ", length(CHECKS), "; anomalies ", length(ANOM))
for a in ANOM[1:min(end, 40)]; println("ANOM: ", a); end
