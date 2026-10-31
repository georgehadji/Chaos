# Independent reader and store-identity checker for study p3-abstention (oracle O15 reader and O19 of the analysis spec, Amendment D93).
# Uses HDF5, SHA, TOML and the standard library only: no Chaoslab code is loaded and no src/Analysis file was read (nor test/analysis/an_p3_ref.jl).
# Written from analysis-spec.md (Amendment D93 E47 to E51, O15, O19; D89 E23 to E26, E29), schema/run-store.md and the preregistrations.
#
# Usage (from the repository root):
#   julia --project=verify -t 4 verify/analysis-p3-export.jl <output folder>
# Inputs (fixed paths, read only): data/p1-le-benchmark and data/p2-indicators with their manifests, the panel file, the P2 and P1 parameter tables,
# the P2 preregistration (for the [[run.h4]] tables). Nothing is written under data/ or in the repository: only the output folder.
# Outputs (folder of argument 1): p3-export.jsonl (one JSON object per line, NaN and +-Infinity written as such) and p3-o19.json (the O19 checks).
#
# Export format.
#   {"k":"meta", ...}
#   {"k":"orbit", "id", "src":"P"|"G", "cls", "control", "K":float|null, "kidx", "in_SD", "in_SL",
#    "ic" (stored text of the relabelling item), "relab": {"lam5_N1","lam5_N1_10","x_0","x_N1_2","x_N1","lam1_N1_rel"} (null where an item or a step is missing),
#    "relab_missing": [items missing], "diag5": bool, "diag1": bool (flag "diagnostics exceeded" on the L5 / L1 item read for the relabelling),
#    "lam1_N1", "lam1_ns" (p2 wrapped L1 item at N_1 and n_s; null when missing), "fl_p2": [orbit flags, union over the p2 Float64 items read by E49 and E51],
#    "fl_p2_all": [union over every p2 item], "fl_p1": [union over every P1 item of the orbit],
#    "old": [instability, structure, forced] (panel non-controls only: the base P1 labelling of D89 E23 a, final after all tests)}
#   {"k":"s", "o", "spec", "prec", "hor", "dir", "t":[recorded steps <= 200000], "c":{column: values}}   series item (L1: lam1; tg: four columns)
#   {"k":"w", "o", "spec", "prec", "hor", "rows":[labels], "c":{D | dig: one value per row}, "flags":[word per row], "bits":[names]}  window item (fq1, pa1, bk1)
using HDF5, SHA, TOML

const T0 = time()
const ROOT = normpath(joinpath(@__DIR__, ".."))
length(ARGS) >= 1 || error("usage: julia --project=verify verify/analysis-p3-export.jl <output folder>")
const OUTD = ARGS[1]
occursin(r"/data(/|$)", replace(abspath(OUTD), "\\" => "/")) && error("output folder inside data/")
mkpath(OUTD)
const P2STORE = joinpath(ROOT, "data/p2-indicators")
const P1STORE = joinpath(ROOT, "data/p1-le-benchmark")
const MANIF2 = joinpath(ROOT, "workspaces/study/stages/03-run/output/p2-indicators-manifest.json")
const MANIF1 = joinpath(ROOT, "workspaces/study/stages/03-run/output/p1-le-benchmark-manifest.json")
const PANELF = joinpath(ROOT, "workspaces/study/stages/02-panel/output/p1-le-benchmark-panel.toml")
const PARF = joinpath(ROOT, "workspaces/study/stages/04-analyze/references/p2-indicators-analysis-parameters.toml")
const P1PARF = joinpath(ROOT, "workspaces/study/stages/04-analyze/references/p1-le-benchmark-analysis-parameters.toml")
const PREREG2 = joinpath(ROOT, "workspaces/study/stages/01-design/output/p2-indicators-preregistration.md")

say(msg...) = println("[", round(time() - T0; digits = 1), " s] ", msg...)
const ANOM = String[]
anom(msg) = (push!(ANOM, msg); say("ANOMALY: ", msg))
const CHECKS = Tuple{String,Bool,String}[]
check(name, ok, detail = "") = (push!(CHECKS, (name, ok, detail)); say(ok ? "PASS " : "FAIL ", name, " ", detail))

const PAR = TOML.parsefile(PARF)
pv(k) = PAR[k]["value"]
const N1, NC = Int(pv("N_1")), Int(pv("N_c"))
const NS = Int(pv("n_s"))
const DIR1 = String(pv("directions")[1])
const GRIDK = String.(pv("grid_k_values"))
const P1PAR = TOML.parsefile(P1PARF)
const PANEL = TOML.parsefile(PANELF)
const ORB = sort(PANEL["orbit"]; by = o -> o["id"])
const PORB = Dict(o["id"] => o for o in ORB)
const MAXPANEL = maximum(o["id"] for o in ORB)
const TRUNC = 200_000   # series rows beyond this step are not exported (no verdict budget reads beyond step 112202)

# ------------------------------------------------------------------ json writer
jnum(x::Real) = x isa Integer ? string(x) : (isnan(x) ? "NaN" : (isinf(x) ? (x > 0 ? "Infinity" : "-Infinity") : repr(Float64(x))))
jstr(s::AbstractString) = "\"" * replace(String(s), "\\" => "\\\\", "\"" => "\\\"") * "\""
function jval(v)
    v === nothing && return "null"
    v isa Bool && return string(v)
    v isa AbstractString && return jstr(v)
    v isa Real && return jnum(v)
    v isa AbstractVector && return "[" * join((jval(x) for x in v), ",") * "]"
    v isa AbstractDict && return "{" * join((jstr(string(k)) * ":" * jval(x) for (k, x) in sort(collect(v); by = first)), ",") * "}"
    error("jval: $(typeof(v))")
end

# ------------------------------------------------------------------ manifests (small JSON reader)
function jparse(s::String)
    pos = Ref(1)
    ws() = (while pos[] <= lastindex(s) && isspace(s[pos[]]); pos[] = nextind(s, pos[]); end)
    function val()
        ws(); c = s[pos[]]
        if c == '{'
            d = Dict{String,Any}(); pos[] += 1; ws()
            if s[pos[]] == '}'; pos[] += 1; return d; end
            while true
                ws(); k = str(); ws(); pos[] += 1
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
const M2 = jparse(read(MANIF2, String))
const M1 = jparse(read(MANIF1, String))
say("manifests: p2 study ", M2["study"], " chunks ", length(M2["chunks"]), " status ", M2["status"], "; P1 study ", M1["study"], " chunks ", length(M1["chunks"]), " status ", M1["status"])
const HASH_P2M = bytes2hex(sha256(read(MANIF2)))
const HASH_P1M = bytes2hex(sha256(read(MANIF1)))
const HASH_PANEL = bytes2hex(sha256(read(PANELF)))
check("E48 p2 manifest study and status", M2["study"] == "p2-indicators" && M2["status"] == "complete", "")
check("E48 P1 manifest study and status", M1["study"] == "p1-le-benchmark" && M1["status"] == "complete", "")
check("E48 every chunk of both manifests ok", all(c -> c["status"] == "ok", M2["chunks"]) && all(c -> c["status"] == "ok", M1["chunks"]), "p2 $(length(M2["chunks"])) P1 $(length(M1["chunks"]))")
check("E48 panel hash in the P1 manifest equals the p2 manifest's and the panel file's", M1["panel"]["sha256_start"] == M2["panel"]["sha256_start"] == HASH_PANEL, HASH_PANEL)

chunk_attrs(f) = Dict(nm => HDF5.read_attribute(f, nm) for nm in keys(HDF5.attributes(f)))

# ================================================================== P1 store
struct Item1
    t::Vector{Int64}; lam::Vector{Float64}; st::Matrix{Float64}; flags::Vector{String}; ic::String
end
tdist(a, b) = (d = a - b; d - 2pi * round(d / (2pi)))
function idx1(it, n::Integer)
    i = searchsortedfirst(it.t, n)
    return (i <= length(it.t) && it.t[i] == n) ? i : 0
end
lamat1(it, n) = (it === nothing ? NaN : (i = idx1(it, n); i == 0 ? NaN : it.lam[i]))

const POS, NONE, UNR = "positive growth supported", "none detected on the tested horizon", "unresolved"
# the base labelling of one panel orbit (D89 E23 a: E14 with its four tests), as the existing independent recomputation of this folder does it
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

const P1ITEMS = Dict{Tuple{Int,String,String,Int,String,String},Item1}()
const P1DUP = Dict{Tuple{Int,String,String,Int,String,String},Int}()
const P1ORBFL = Dict{Int,Set{String}}()
const P1DIAG = Ref(0)      # L1 or L5 items of the P1 store (any variant) with the flag "diagnostics exceeded"
function read_p1()
    files = sort(filter(fn -> occursin(r"^chunk-\d{5}\.h5$", fn), readdir(P1STORE)))
    say("reading ", length(files), " P1 chunk files")
    for fn in files
        h5open(joinpath(P1STORE, fn), "r") do f
            a = chunk_attrs(f)
            a["study"] == "p1-le-benchmark" || anom("$fn study")
            Bool(a["complete"]) || anom("$fn incomplete")
            m = a["method"]
            prec = startswith(a["precision"], "BigFloat") ? "BigFloat" : a["precision"]
            hor = Int(a["horizon"]); tau = Int(a["tau"]); dir = get(a, "direction", ""); mode = a["coordinate_mode"]
            ok = a["status"] == "ok"
            for on in keys(f["orbits"])
                g = f["orbits"][on]; id = parse(Int, on)
                fl = String.(HDF5.read_attribute(g, "flags"))
                union!(get!(P1ORBFL, id, Set{String}()), fl)
                (m == "L1" || m == "L5") && ("diagnostics exceeded" in fl) && (P1DIAG[] += 1)
                ((m in ("L1", "L5", "L3", "L4", "L4b", "L5b", "L6", "L9")) && ok && tau == 1 && mode == "lifted") || continue
                t = read(g["t"]); L = read(g["lyap"][m])
                lam = ndims(L) == 1 ? collect(L) : collect(L[:, 1])
                st = (m == "L5") ? read(g["state"]) : zeros(0, 2)
                key = (id, m, prec, hor, dir, mode)
                if haskey(P1ITEMS, key)
                    P1DUP[key] = get(P1DUP, key, 1) + 1
                else
                    P1ITEMS[key] = Item1(collect(t), lam, st, fl, String(HDF5.read_attribute(g, "ic")))
                end
            end
        end
    end
end

# ================================================================== p2 store
const SER = Dict{Tuple{Int,String,String,Int,String},Any}()       # (orbit, spec, prec, horizon, dir) -> (t, columns, flags, extra)
const WIN = Dict{Tuple{Int,String,String,Int},Any}()               # (orbit, spec, prec, horizon) -> (rows, cols, flagwords, bits, flags)
const P2INFO = Dict{Int,Dict{String,Any}}()
const P2L5 = Dict{Int,Any}()                                       # grid L5 F64 N_1 wrapped dir pi/4 -> (t, lam, xlift, flags, ic)
const P2DIAG = Ref(0)
const P2DUP = Ref(0)
const SERIES_COLS = Dict("L1" => ["lam1"], "tg" => ["FLI.fli", "MEGNO.Ybar20_1", "SALI.sali", "SALI.L"])
function read_p2()
    files = sort(filter(fn -> occursin(r"^chunk-\d{5}\.h5$", fn), readdir(P2STORE)))
    say("reading ", length(files), " p2 chunk files")
    for fn in files
        h5open(joinpath(P2STORE, fn), "r") do f
            a = chunk_attrs(f)
            a["study"] == "p2-indicators" || anom("$fn study")
            Bool(a["complete"]) || anom("$fn incomplete")
            a["status"] == "ok" || anom("$fn status")
            m = a["method"]; prec = a["precision"]
            hor = Int(a["horizon"]); dir = get(a, "direction", ""); mode = a["coordinate_mode"]
            kind = get(a, "kind", "")
            for on in keys(f["orbits"])
                g = f["orbits"][on]; id = parse(Int, on)
                fl = String.(HDF5.read_attribute(g, "flags"))
                info = get!(P2INFO, id, Dict{String,Any}("fl_all" => Set{String}(), "fl_read" => Set{String}(), "params" => String(HDF5.read_attribute(g, "params")),
                                                         "ic" => String(HDF5.read_attribute(g, "ic")), "dc" => String(HDF5.read_attribute(g, "design_class"))))
                union!(info["fl_all"], fl)
                (m == "L1" || m == "L5") && ("diagnostics exceeded" in fl) && (P2DIAG[] += 1)
                t = read(g["t"])
                tr = findlast(<=(TRUNC), t); tr === nothing && (tr = 0)
                if kind == ""
                    (m == "L1" || m == "L5") || continue
                    L = read(g["lyap"][m]); lam = ndims(L) == 1 ? collect(L) : collect(L[:, 1])
                    if m == "L1" && dir == DIR1 && (prec == "Float64" ? hor == N1 : true)
                        key = (id, "L1", prec, hor, dir)
                        haskey(SER, key) && (P2DUP[] += 1)
                        ex = Dict{String,Any}()
                        if hor == N1
                            iN = findfirst(==(N1), t); is = findfirst(==(NS), t)
                            ex["lam_N1"] = iN === nothing ? nothing : lam[iN]
                            ex["lam_ns"] = is === nothing ? nothing : lam[is]
                        end
                        SER[key] = (collect(t[1:tr]), Dict("lam1" => lam[1:tr]), fl, ex)
                        prec == "Float64" && hor == N1 && union!(info["fl_read"], fl)
                    elseif m == "L5" && prec == "Float64" && hor == N1 && dir == DIR1 && mode == "wrapped"
                        S = read(g["state"]); W = read(g["winding"])
                        haskey(P2L5, id) && (P2DUP[] += 1)
                        P2L5[id] = (collect(t), lam, S[:, 1] .+ 2pi .* Float64.(W[:, 1]), fl, String(HDF5.read_attribute(g, "ic")))
                        union!(info["fl_read"], fl)
                    end
                elseif kind == "tangent"
                    (m == "tg" && dir == DIR1) || continue
                    ds = g["indicators"][m]["series"]
                    cols = String.(HDF5.read_attribute(ds, "columns")); M = read(ds)
                    c = Dict{String,Vector{Float64}}()
                    for cn in SERIES_COLS[m]
                        j = findfirst(==(cn), cols); j === nothing && error("$fn column $cn")
                        c[cn] = collect(M[1:tr, j])
                    end
                    key = (id, "tg", prec, hor, dir)
                    haskey(SER, key) && (P2DUP[] += 1)
                    SER[key] = (collect(t[1:tr]), c, fl, Dict{String,Any}())
                    prec == "Float64" && hor == N1 && union!(info["fl_read"], fl)
                elseif kind in ("frequency", "phase_advance", "birkhoff")
                    (m == "fq1" || m == "pa1" || m == "bk1") || continue
                    grp = g["indicators"][m]
                    ds = grp["final"]; cols = String.(HDF5.read_attribute(ds, "columns")); rows = String.(HDF5.read_attribute(ds, "rows"))
                    M = read(ds)
                    c = Dict{String,Vector{Float64}}()
                    for cn in ("D", "dig")
                        j = findfirst(==(cn), cols); j === nothing || (c[cn] = collect(M[:, j]))
                    end
                    fds = grp["flags"]; fw = collect(Int64, read(fds)); bitnames = String.(HDF5.read_attribute(fds, "bits"))
                    key = (id, m, prec, hor)
                    haskey(WIN, key) && (P2DUP[] += 1)
                    WIN[key] = (rows, c, fw, bitnames, fl)
                    prec == "Float64" && union!(info["fl_read"], fl)
                end
            end
        end
    end
end

read_p1()
say("P1 items kept: ", length(P1ITEMS), "; duplicate keys: ", length(P1DUP))
read_p2()
say("p2 series items ", length(SER), ", window items ", length(WIN), ", orbits ", length(P2INFO), ", duplicate keys ", P2DUP[])

# ------------------------------------------------------------------ grid identity from the p2 preregistration ([[run.h4]] tables)
txt = read(PREREG2, String)
i0 = findfirst("\n[[run.h4]]\n", txt)
body = txt[(first(i0) + 1):end]
j0 = findfirst("\n```", body)
const TABS = TOML.parse(body[1:first(j0)])["run"]["h4"]
check("E24 [[run.h4]] table count", length(TABS) == Int(pv("grid_orbit_count")), string(length(TABS)))
const GRIDID = Dict{Int,Tuple{Int,Int,Int,String}}()
let nbad = 0
    for (n, tb) in enumerate(TABS)
        id = MAXPANEL + n
        m = match(r"^grid-k(\d+)-i(\d+)-j(\d+)$", tb["label"]); m === nothing && (nbad += 1; continue)
        ki, ii, jj = parse(Int, m[1]), parse(Int, m[2]), parse(Int, m[3])
        kk = findfirst(==(tb["K"]), GRIDK)
        info = get(P2INFO, id, nothing)
        if kk != ki || info === nothing || parse(Float64, info["params"]) != parse(Float64, tb["K"]) || info["ic"] != tb["x0"] * "," * tb["p0"]
            nbad += 1; continue
        end
        GRIDID[id] = (ki, ii, jj, tb["K"])
    end
    check("E24 grid identity (params, ic, K index)", nbad == 0 && length(GRIDID) == 2304, "bad $nbad")
    nbp = 0
    for o in ORB
        info = get(P2INFO, o["id"], nothing)
        if info === nothing || info["ic"] != o["x0_dec"] * "," * o["p0_dec"]; nbp += 1; continue; end
        if haskey(o, "K") && parse(Float64, info["params"]) != Float64(o["K"]); nbp += 1; end
    end
    check("E24 panel identity (ic, params)", nbp == 0, "bad $nbp")
end

# ================================================================== orbit records and the relabelling inputs (E49)
p1get(id, m, hor) = get(P1ITEMS, (id, m, "Float64", hor, DIR1, "lifted"), nothing)

o19 = Dict{String,Any}()
panel_missing = Int[]; grid_missing = Int[]
open(joinpath(OUTD, "p3-export.jsonl"), "w") do io
    println(io, jval(Dict("k" => "meta", "N_1" => N1, "N_c" => NC, "n_s" => NS, "dir" => DIR1, "trunc" => TRUNC, "p2_manifest_sha256" => HASH_P2M,
                          "p1_manifest_sha256" => HASH_P1M, "panel_sha256" => HASH_PANEL)))
    for id in 1:(MAXPANEL + length(TABS))
        isP = id <= MAXPANEL
        o = isP ? PORB[id] : nothing
        info = get(P2INFO, id, nothing)
        d = Dict{String,Any}("k" => "orbit", "id" => id, "src" => isP ? "P" : "G")
        relab = Dict{String,Any}("lam5_N1" => nothing, "lam5_N1_10" => nothing, "x_0" => nothing, "x_N1_2" => nothing, "x_N1" => nothing, "lam1_N1_rel" => nothing)
        miss = String[]
        diag5 = false; diag1 = false; ic = nothing
        if isP
            d["cls"] = o["class"]; d["control"] = o["control"]; d["K"] = haskey(o, "K") ? Float64(o["K"]) : nothing
            d["kidx"] = Int(o["k_index"]); d["in_SD"] = get(o, "in_SD", false); d["in_SL"] = get(o, "in_SL", false)
            if !o["control"]
                l5 = p1get(id, "L5", N1); l1 = p1get(id, "L1", N1)
                if l5 === nothing
                    push!(miss, "L5")
                else
                    ic = l5.ic; diag5 = "diagnostics exceeded" in l5.flags
                    i10 = idx1(l5, N1 ÷ 10); ih = idx1(l5, N1 ÷ 2); iN = idx1(l5, N1)
                    i10 == 0 ? push!(miss, "L5 step N_1/10") : (relab["lam5_N1_10"] = l5.lam[i10])
                    ih == 0 ? push!(miss, "L5 step N_1/2") : (relab["x_N1_2"] = l5.st[ih, 1])
                    iN == 0 ? push!(miss, "L5 step N_1") : (relab["lam5_N1"] = l5.lam[iN]; relab["x_N1"] = l5.st[iN, 1])
                    relab["x_0"] = parse(Float64, split(l5.ic, ",")[1])
                end
                if l1 === nothing
                    push!(miss, "L1")
                else
                    diag1 = "diagnostics exceeded" in l1.flags
                    iN1 = idx1(l1, N1)
                    iN1 == 0 ? push!(miss, "L1 step N_1") : (relab["lam1_N1_rel"] = l1.lam[iN1])
                end
                isempty(miss) || push!(panel_missing, id)
            end
        else
            ki, ii, jj, ktxt = GRIDID[id]
            d["cls"] = "h4"; d["control"] = false; d["K"] = parse(Float64, ktxt); d["kidx"] = ki; d["in_SD"] = false; d["in_SL"] = false
            l5 = get(P2L5, id, nothing)
            l1 = get(SER, (id, "L1", "Float64", N1, DIR1), nothing)
            if l5 === nothing
                push!(miss, "L5")
            else
                t, lam, xl, fl, icx = l5
                ic = icx; diag5 = "diagnostics exceeded" in fl
                i10 = findfirst(==(N1 ÷ 10), t); ih = findfirst(==(N1 ÷ 2), t); iN = findfirst(==(N1), t)
                i10 === nothing ? push!(miss, "L5 step N_1/10") : (relab["lam5_N1_10"] = lam[i10])
                ih === nothing ? push!(miss, "L5 step N_1/2") : (relab["x_N1_2"] = xl[ih])
                iN === nothing ? push!(miss, "L5 step N_1") : (relab["lam5_N1"] = lam[iN]; relab["x_N1"] = xl[iN])
                relab["x_0"] = parse(Float64, split(icx, ",")[1])
            end
            if l1 === nothing
                push!(miss, "L1")
            else
                diag1 = "diagnostics exceeded" in l1[3]
                v = get(l1[4], "lam_N1", nothing)
                v === nothing ? push!(miss, "L1 step N_1") : (relab["lam1_N1_rel"] = v)
            end
            isempty(miss) || push!(grid_missing, id)
        end
        d["ic"] = ic; d["relab"] = relab; d["relab_missing"] = miss; d["diag5"] = diag5; d["diag1"] = diag1
        l1p2 = info === nothing ? nothing : get(SER, (id, "L1", "Float64", N1, DIR1), nothing)
        d["lam1_N1"] = l1p2 === nothing ? nothing : get(l1p2[4], "lam_N1", nothing)
        d["lam1_ns"] = l1p2 === nothing ? nothing : get(l1p2[4], "lam_ns", nothing)
        d["fl_p2"] = info === nothing ? String[] : sort(collect(info["fl_read"]))
        d["fl_p2_all"] = info === nothing ? String[] : sort(collect(info["fl_all"]))
        d["fl_p1"] = isP ? sort(collect(get(P1ORBFL, id, Set{String}()))) : String[]
        if isP && !o["control"]
            # the base P1 labelling (final label after all tests), as the existing recomputation of this folder does it
            P = P1PAR
            n1, nlong = Int(P["N_1"]), Int(P["N_long"])
            l5l = get(P1ITEMS, (id, "L5", "Float64", nlong, DIR1, "lifted"), nothing)
            prim5 = p1get(id, "L5", n1)
            Nlab = l5l !== nothing ? nlong : n1
            litem = l5l !== nothing ? l5l : prim5
            rinst, rstruct = labelof(P, litem, Nlab, o)
            forced = String[]
            if get(o, "in_SL", false)
                nb = Int(P["N_B"])
                d64 = get(P1ITEMS, (id, "L5", "Double64", n1, DIR1, "lifted"), nothing); bf = get(P1ITEMS, (id, "L5", "BigFloat", nb, DIR1, "lifted"), nothing)
                if prim5 !== nothing && d64 !== nothing && bf !== nothing
                    lf, ld, lb = lamat1(prim5, nb), lamat1(d64, nb), lamat1(bf, nb)
                    tc = intersect(prim5.t, d64.t, bf.t); tc = tc[(tc .> 0) .& (tc .< nb)]
                    hc = tc[argmin(abs.(tc .- nb / 2))]
                    pf, pd, pb = labelof(P, prim5, nb, o; hfix = Int(hc)), labelof(P, d64, nb, o; hfix = Int(hc)), labelof(P, bf, nb, o; hfix = Int(hc))
                    (abs(lf - ld) > P["ladder_tol"] || abs(ld - lb) > P["ladder_tol"] || !(pf == pd == pb)) && push!(forced, "ladder")
                else
                    anom("P1 orbit $id (S_L): ladder item missing")
                end
            end
            VM = ["L1", "L3", "L4", "L4b", "L5", "L5b", "L6"]
            vals_B = [lamat1(p1get(id, m, n1), n1) for m in VM]
            dqb = maximum(vals_B) - minimum(vals_B)
            d9 = abs(lamat1(p1get(id, "L9", n1), n1) - lamat1(prim5, n1))
            (dqb > P["h1a_bound"] || d9 > P["h1b_bound"]) && push!(forced, "estimator")
            if get(o, "in_SD", false)
                labs = Set{Tuple{String,String}}(); nd = 0
                for ((kid, km, kp, kh, kd, kmode), it) in P1ITEMS
                    (kid == id && km == "L5" && kp == "Float64" && kh == n1 && kmode == "lifted") || continue
                    nd += 1; push!(labs, labelof(P, it, n1, o))
                end
                (nd >= 2 && length(labs) > 1) && push!(forced, "direction")
            end
            prim5 !== nothing && ("diagnostics exceeded" in prim5.flags) && push!(forced, "diagnostics")
            "excluded" in get(P1ORBFL, id, Set{String}()) && push!(forced, "excluded")
            order = ["ladder", "estimator", "direction", "diagnostics", "excluded"]
            forced = filter(x -> x in forced, order)
            inst, stc = isempty(forced) ? (rinst, rstruct) : (UNR, UNR)
            d["old"] = [inst, stc, forced]
        end
        println(io, jval(d))
    end
    ns = 0
    for (key, (t, c, fl, extra)) in sort(collect(SER); by = x -> x[1])
        id, spec, prec, hor, dir = key
        println(io, jval(Dict("k" => "s", "o" => id, "spec" => spec, "prec" => prec, "hor" => hor, "dir" => dir, "t" => t, "c" => c)))
        ns += 1
    end
    nw = 0
    for (key, (rows, c, fw, bits, fl)) in sort(collect(WIN); by = x -> x[1])
        id, spec, prec, hor = key
        println(io, jval(Dict("k" => "w", "o" => id, "spec" => spec, "prec" => prec, "hor" => hor, "rows" => rows, "c" => c, "flags" => fw, "bits" => bits)))
        nw += 1
    end
    say("export: series ", ns, ", windows ", nw)
end

# ================================================================== O19 checks
let
    nonctl = [o["id"] for o in ORB if !o["control"]]
    check("O19a P1 store holds the L5 and L1 items (Float64, N_1, lifted, direction pi/4) and the steps N_1/10, N_1/2, N_1 for every non-control panel orbit",
          isempty(panel_missing) && length(nonctl) == 870, "non-control panel orbits $(length(nonctl)), with a missing item or step: $(length(panel_missing))")
    ok1 = all(id -> haskey(SER, (id, "L1", "Float64", N1, DIR1)), 1:3187)
    ok5 = all(id -> haskey(P2L5, id), 884:3187)
    check("O19b p2 store holds the L1 item (Float64, N_1, wrapped) for all 883 panel and 2304 grid orbits", ok1, "")
    check("O19b p2 store holds the L5 item (Float64, N_1, wrapped) for the 2304 grid orbits", ok5, "")
    check("O19b p2 store holds no such L5 item for the panel orbits (as the coordinator reported)", !any(id -> haskey(P2L5, id), 1:MAXPANEL), "")
    check("O19c no L1 or L5 item of the P1 store carries \"diagnostics exceeded\"", P1DIAG[] == 0, "items with the flag: $(P1DIAG[])")
    check("O19c no L1 or L5 item of the p2 store carries \"diagnostics exceeded\"", P2DIAG[] == 0, "items with the flag: $(P2DIAG[])")
    check("no duplicate item keys among the items read (P1 lifted items kept; p2 series, windows)", isempty(P1DUP) && P2DUP[] == 0, "P1 $(length(P1DUP)) p2 $(P2DUP[])")
    o19["P1_duplicate_keys"] = [string(k) * " x" * string(v) for (k, v) in P1DUP]
    o19["panel_missing_ids"] = panel_missing
    o19["grid_missing_ids"] = grid_missing
end
o19["checks"] = [Dict("name" => n, "ok" => ok, "detail" => dt) for (n, ok, dt) in CHECKS]
o19["anomalies"] = ANOM
open(joinpath(OUTD, "p3-o19.json"), "w") do io
    println(io, jval(o19))
end
say("done; checks ", count(c -> c[2], CHECKS), " pass of ", length(CHECKS), "; anomalies ", length(ANOM))
