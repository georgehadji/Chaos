# Independent recomputation of the analysis of study p1-le-benchmark (oracle O3 of the analysis spec).
# Uses HDF5, SHA, TOML, Random, Dates and Base only: no Chaoslab code is loaded and no module source was read.
# Usage: julia --project=verify verify/analysis-independent.jl <store dir> <manifest json> <panel toml> <parameter toml> <out file>
# Written from analysis-spec.md, run-store.md and the preregistration; the readings a-j of the task owner are followed.
using HDF5, SHA, TOML, Random, Dates

const T0 = time()
length(ARGS) == 5 || error("usage: analysis-independent.jl <store dir> <manifest json> <panel toml> <parameter toml> <out file>")
const STORE, MANIFEST, PANELF, PARF, OUTF = ARGS
say(msg...) = println("[", round(time() - T0; digits = 1), " s] ", msg...)
const ANOM = String[]
anom(msg) = (push!(ANOM, msg); say("ANOMALY: ", msg))

# ---------------------------------------------------------------- parameters, panel, manifest
const PAR = TOML.parsefile(PARF)
const N1, NLONG, NB, NH2 = Int(PAR["N_1"]), Int(PAR["N_long"]), Int(PAR["N_B"]), Int(PAR["n_H2"])
const MINSET = Int(PAR["min_set"])
const VM = ["L1", "L3", "L4", "L4b", "L5", "L5b", "L6", "L9"]
const BM = filter(!=("L9"), VM)
const POS, NONE, UNR = "positive growth supported", "none detected on the tested horizon", "unresolved"

const PANEL = TOML.parsefile(PANELF)
const ORB = sort(PANEL["orbit"]; by = o -> o["id"])
const MAXPANEL = maximum(o["id"] for o in ORB)
const PANELIDS = Set(o["id"] for o in ORB)
const SEED0 = parse(UInt64, string(ORB[1]["run_seed"]))
say("panel orbits ", length(ORB), ", max id ", MAXPANEL, ", seed ", SEED0)

mtxt = read(MANIFEST, String)
const STUDY = let m = match(r"\"study\"\s*:\s*\"([^\"]+)\"", mtxt)
    m === nothing ? error("no study in manifest") : String(m.captures[1])
end
nman = length(collect(eachmatch(r"\"ordinal\"\s*:\s*\d+", mtxt)))
say("study ", STUDY, ", manifest lists ", nman, " chunk ordinals")

# ---------------------------------------------------------------- store reader
struct Item
    t::Vector{Int64}
    lam::Vector{Float64}      # first exponent, column 1 of lyap/<method>
    st::Matrix{Float64}       # (n, 2) recorded point, only for method L5 (else 0 x 2)
    flags::Vector{String}
    ctrldev::Float64          # stored control_deviation or NaN
end
const ITEMS = Dict{Tuple{Int,String,String,Int,String,String},Item}()
const ORBFLAGS = Dict{Int,Set{String}}()
const ORBPARAMS = Dict{Int,String}()
const EPS0 = Float64(PAR["h1c_eps"])
nfile = 0; nkept = 0; nskipped_nonok = 0; nblank = 0

files = sort(filter(fn -> occursin(r"^chunk-\d{5}\.h5$", fn), readdir(STORE)))
say("reading ", length(files), " chunk files")
for fn in files
    global nfile += 1
    h5open(joinpath(STORE, fn), "r") do f
        a = Dict(nm => HDF5.read_attribute(f, nm) for nm in keys(HDF5.attributes(f)))
        a["study"] == STUDY || anom("$fn: study $(a["study"])")
        get(a, "schema_version", 0) == 2 || anom("$fn: schema_version")
        Bool(a["complete"]) || anom("$fn: not complete")
        Bool(a["pilot"]) && anom("$fn: pilot file")
        m = a["method"]
        prec = startswith(a["precision"], "BigFloat") ? "BigFloat" : a["precision"]
        hor = Int(a["horizon"]); tau = Int(a["tau"]); dir = a["direction"]; mode = a["coordinate_mode"]
        pert = m == "L2" ? parse(Float64, string(a["perturbation"])) : NaN
        ok = a["status"] == "ok"
        ok || (global nskipped_nonok += 1)
        keep = ok && tau == 1 && !(mode == "wrapped" && m != "L2") && (m != "L2" || pert == EPS0)
        for on in keys(f["orbits"])
            g = f["orbits"][on]; id = parse(Int, on)
            fl = String.(HDF5.read_attribute(g, "flags"))
            union!(get!(ORBFLAGS, id, Set{String}()), fl)
            haskey(ORBPARAMS, id) || (ORBPARAMS[id] = string(HDF5.read_attribute(g, "params")))
            keep || continue
            t = read(g["t"]); L = read(g["lyap"][m])
            lam = ndims(L) == 1 ? collect(L) : collect(L[:, 1])
            st = m == "L5" ? read(g["state"]) : zeros(0, 2)
            length(t) == length(lam) || anom("$fn orbit $id: length of t and lyap differ")
            if any(==(0), t) || any(!isfinite, lam)
                global nblank += 1
                anom("$fn orbit $id $m: blank or non-finite rows (t zeros $(count(==(0), t)), non-finite lam $(count(!isfinite, lam)))")
            end
            cd = haskey(HDF5.attributes(g), "control_deviation") ? Float64(HDF5.read_attribute(g, "control_deviation")) : NaN
            key = (id, m, prec, hor, dir, mode)
            it = Item(t, lam, st, fl, cd)
            if haskey(ITEMS, key)
                old = ITEMS[key]
                (old.t == t && isequal(old.lam, lam)) || anom("duplicate item $key differs")
            else
                ITEMS[key] = it
                global nkept += 1
            end
        end
    end
    nfile % 200 == 0 && say("  files read ", nfile)
end
say("items kept ", nkept, ", non-ok files ", nskipped_nonok, ", files ", nfile, " (manifest ordinals ", nman, ")")
nfile == nman || anom("store has $nfile chunk files, manifest lists $nman ordinals")

# ---------------------------------------------------------------- helpers
function idx(it::Item, n::Integer)
    i = searchsortedfirst(it.t, n)
    return (i <= length(it.t) && it.t[i] == n) ? i : 0
end
lamat(it, n) = (it === nothing ? NaN : (i = idx(it, n); i == 0 ? NaN : it.lam[i]))
getit(id, m, prec, hor; dir = "pi/4", mode = (m == "L2" ? "wrapped" : "lifted")) = get(ITEMS, (id, m, prec, hor, dir, mode), nothing)

tdist(a, b) = (d = a - b; d - 2pi * round(d / (2pi)))
function typeq7(sorted::Vector{Float64}, p)
    n = length(sorted); h = (n - 1) * p + 1; lo = floor(Int, h)
    lo >= n && return sorted[n]
    return sorted[lo] + (h - lo) * (sorted[lo + 1] - sorted[lo])
end
function medsorted!(v::Vector{Float64})
    sort!(v); m = length(v)
    return isodd(m) ? v[(m + 1) >> 1] : (v[m >> 1] + v[(m >> 1) + 1]) / 2
end
median_of(v) = medsorted!(collect(Float64, v))
function seedfor(label::String)
    h = sha256(codeunits(label))
    u = zero(UInt64)
    for i in 1:8; u = (u << 8) | UInt64(h[i]); end
    return u + SEED0     # wrapping UInt64 addition
end
function bootmedians(v::Vector{Float64}, label::String, B::Int)
    rng = Xoshiro(seedfor(label)); m = length(v)
    out = Vector{Float64}(undef, B); buf = Vector{Float64}(undef, m)
    for b in 1:B
        for i in 1:m; buf[i] = v[rand(rng, 1:m)]; end
        out[b] = medsorted!(buf)
    end
    return out
end
function interval(vals::Vector{Float64})
    lvl = Float64(PAR["bootstrap_level"]); s = sort(vals)
    return typeq7(s, (1 - lvl) / 2), typeq7(s, 1 - (1 - lvl) / 2)
end

# ---------------------------------------------------------------- labels (E14, readings e)
# label at horizon N from an L5 item; returns (inst, struct)
function labelof(it, N::Int, o; hfix::Int = 0)
    it === nothing && return (UNR, UNR)
    iN = idx(it, N); i10 = idx(it, N ÷ 10)
    (iN == 0 || i10 == 0) && return (UNR, UNR)
    lam = it.lam[iN]; lam10 = it.lam[i10]
    inst = UNR
    if isfinite(lam) && isfinite(lam10) && lam10 > 0
        r = 10 * lam / lam10
        if lam >= PAR["lambda_threshold"] && r >= PAR["growth_ratio"]
            inst = POS
        elseif lam < PAR["lambda_threshold"] && r < PAR["growth_ratio"]
            inst = NONE
        end
    end
    x0, p0 = Float64(o["x0"]), Float64(o["p0"])
    stc = UNR
    if get(o, "control", false) && get(o, "control_kind", "none") in ("fixed", "newton")
        iq = idx(it, Int(o["period"]))
        if iq != 0
            res = hypot(tdist(it.st[iq, 1], x0), tdist(it.st[iq, 2], p0))
            res <= PAR["return_tol"] && (stc = "periodic")
        end
    end
    if stc == UNR
        if inst == POS
            stc = "irregular"
        elseif inst == NONE
            # recorded step nearest N/2 among steps below N
            best = 0; bd = Inf
            hfix > 0 && (best = idx(it, hfix))
            hfix == 0 && for j in 1:(iN - 1)
                d = abs(it.t[j] - N / 2)
                if d < bd; bd = d; best = j; end
            end
            if best > 0
                h = Int(it.t[best])
                rho1 = (it.st[best, 1] - x0) / (2pi * h)
                rho2 = (it.st[iN, 1] - it.st[best, 1]) / (2pi * (N - h))
                abs(rho1 - rho2) <= PAR["rotation_tol"] && (stc = "quasiperiodic")
            end
        end
    end
    return (inst, stc)
end

# ---------------------------------------------------------------- per-orbit quantities
struct Rec
    id::Int; cls::String; K::Float64; N::Int
    lamref::Float64; dqb::Float64; d9::Float64; d2::Float64; R::Float64; nconv::Float64; ndef::Bool
    rinst::String; rstruct::String; inst::String; stc::String
    forced::Vector{String}; ladder::Union{Nothing,Bool}; flags::Vector{String}
end
function nconv_method(it, ref)
    (isfinite(ref) && ref > 0) || return (Float64(N1), false)
    band = PAR["h3_band"] * ref
    n = length(it.t); iN = idx(it, N1)
    out = 0
    for i in 1:n
        (it.t[i] > N1) && break
        abs(it.lam[i] - ref) > band && (out = i)
    end
    (iN == 0 || out >= iN) && return (Float64(N1), false)   # outside the band at N_1 (or N_1 not recorded)
    return (Float64(it.t[out + 1]), true)
end

recs = Rec[]
nbad_main = 0
for o in ORB
    id = o["id"]; cls = String(o["class"]); K = Float64(get(o, "K", NaN))
    prim = Dict(m => getit(id, m, "Float64", N1) for m in VM)
    l2 = getit(id, "L2", "Float64", N1)
    if any(isnothing, values(prim)) || l2 === nothing
        anom("orbit $id: missing main item(s) " * join([m for m in VM if prim[m] === nothing], ",") * (l2 === nothing ? ",L2" : ""))
        global nbad_main += 1
    end
    l5long = getit(id, "L5", "Float64", NLONG)
    Nlab = l5long !== nothing ? NLONG : N1
    litem = l5long !== nothing ? l5long : prim["L5"]
    lamref = lamat(litem, Nlab)
    lN(m) = lamat(prim[m], N1)
    vals_B = [lN(m) for m in BM]
    dqb = maximum(vals_B) - minimum(vals_B)
    d9 = abs(lN("L9") - lN("L5"))
    d2 = abs(lamat(l2, N1) - lN("L5"))
    vH2 = [lamat(prim[m], NH2) for m in VM]
    R = (maximum(vH2) - minimum(vH2)) / lamref
    # n_conv
    ndef = true; nc = 0.0
    for m in VM
        (n, d) = prim[m] === nothing ? (Float64(N1), false) : nconv_method(prim[m], lamref)
        nc = max(nc, n); ndef &= d
    end
    # labels
    rinst, rstruct = labelof(litem, Nlab, o)
    forced = String[]
    ladder = nothing
    fl = sort(collect(get(ORBFLAGS, id, Set{String}())))
    # (a) ladder, S_L only
    if get(o, "in_SL", false)
        f64 = prim["L5"]; d64 = getit(id, "L5", "Double64", N1); bf = getit(id, "L5", "BigFloat", NB)
        if f64 === nothing || d64 === nothing || bf === nothing
            anom("orbit $id (S_L): ladder item missing (F64 $(f64 !== nothing), D64 $(d64 !== nothing), BF $(bf !== nothing))")
        else
            lf, ld, lb = lamat(f64, NB), lamat(d64, NB), lamat(bf, NB)
            tc = intersect(f64.t, d64.t, bf.t); tc = tc[(tc .> 0) .& (tc .< NB)]
            hc = tc[argmin(abs.(tc .- NB / 2))]
            pf, pd, pb = labelof(f64, NB, o; hfix = Int(hc)), labelof(d64, NB, o; hfix = Int(hc)), labelof(bf, NB, o; hfix = Int(hc))
            ladder = (pf == pd == pb)
            (abs(lf - ld) > PAR["ladder_tol"] || abs(ld - lb) > PAR["ladder_tol"] || !ladder) && push!(forced, "ladder")
        end
    end
    # (b) estimators
    (dqb > PAR["h1a_bound"] || d9 > PAR["h1b_bound"]) && push!(forced, "estimator")
    # (c) direction, S_D only
    if get(o, "in_SD", false)
        labs = Set{Tuple{String,String}}()
        nd = 0
        for ((kid, km, kp, kh, kd, kmode), it) in ITEMS
            (kid == id && km == "L5" && kp == "Float64" && kh == N1 && kmode == "lifted") || continue
            nd += 1; push!(labs, labelof(it, N1, o))
        end
        (nd >= 2 && length(labs) > 1) && push!(forced, "direction")
        nd >= 2 || anom("orbit $id (S_D): only $nd direction item(s)")
    end
    # (d) diagnostics on the primary L5 item at N_1
    prim["L5"] !== nothing && ("diagnostics exceeded" in prim["L5"].flags) && push!(forced, "diagnostics")
    # (e) excluded
    "excluded" in fl && push!(forced, "excluded")
    order = ["ladder", "estimator", "direction", "diagnostics", "excluded"]
    forced = filter(x -> x in forced, order)
    inst, stc = isempty(forced) ? (rinst, rstruct) : (UNR, UNR)
    push!(recs, Rec(id, cls, K, Nlab, lamref, dqb, d9, d2, R, nc, ndef, rinst, rstruct, inst, stc, forced, ladder, fl))
end
say("per-orbit quantities done for ", length(recs), " orbits (", nbad_main, " with missing main items)")
REC = Dict(r.id => r for r in recs)
PORB = Dict(o["id"] => o for o in ORB)

# ---------------------------------------------------------------- cross-checks of stored facts (anomalies only)
for r in recs
    mp = 0.0
    for h in (N1, NLONG)
        it = getit(r.id, "L5", "Float64", h)
        it === nothing || (mp = max(mp, maximum(abs, it.st[:, 2])))
    end
    (mp > PAR["drift_bound"]) == ("drift" in r.flags) || anom("orbit $(r.id): drift flag $("drift" in r.flags) vs max|p| $mp over the recorded points of the lifted L5 items (N_1, N_long)")
    o = PORB[r.id]
    if get(o, "control", false) && haskey(o, "exponent") && o["exponent"] > 0
        for m in VM
            it = getit(r.id, m, "Float64", N1)
            (it === nothing || !isfinite(it.ctrldev)) && continue
            dev = abs(lamat(it, N1) - o["exponent"])
            abs(dev - it.ctrldev) > 1e-12 && anom("control $(r.id) $m: stored control_deviation $(it.ctrldev) vs recomputed $dev")
        end
    end
end

# ---------------------------------------------------------------- sets (readings d)
function isflagged(r)
    return ("excluded" in r.flags) ? "excluded" : (("drift" in r.flags) ? "drift" : "")
end
Kc = Float64(PAR["K_c"])
inE(r) = (r.cls == "C" && r.K > 0) || (r.cls == "control" && get(PORB[r.id], "exponent", 0.0) > 0)
inW(r) = r.inst == POS && ((r.cls in ("S", "X") && r.K > 0) || (r.cls == "C" && 0 < r.K <= Kc))
inS(r) = r.cls == "C" && (abs(r.K - 2pi) <= 1e-9 || abs(r.K - 10) <= 1e-9)
dropped = Tuple{Int,String}[]
Eset = Int[]; Wset = Int[]; Sset = Int[]
for r in recs
    mem = (inE(r), inW(r), inS(r))
    why = isflagged(r)
    if any(mem) && why != ""
        push!(dropped, (r.id, why))
    elseif why == ""
        mem[1] && push!(Eset, r.id); mem[2] && push!(Wset, r.id); mem[3] && push!(Sset, r.id)
    end
end
sort!(Eset); sort!(Wset); sort!(Sset)

# ---------------------------------------------------------------- hypotheses
out = String[]
idl(v) = isempty(v) ? "-" : join(v, ",")
for r in recs
    push!(out, join(["ORBIT", r.id, r.cls, repr(r.K), r.N, repr(r.lamref), repr(r.dqb), repr(r.d9), repr(r.d2), repr(r.R), repr(r.nconv),
                     r.ndef, r.rinst, r.rstruct, r.inst, r.stc, idl(r.forced), r.ladder === nothing ? "-" : string(r.ladder), idl(r.flags)], "|"))
end
push!(out, "SET|E_ch|" * idl(Eset)); push!(out, "SET|W|" * idl(Wset)); push!(out, "SET|S|" * idl(Sset))
push!(out, "SET|dropped|" * idl(["$(a):$(b)" for (a, b) in sort(dropped)]))
push!(out, "SET|S_nonpositive|" * string(count(i -> !(REC[i].lamref > 0), Sset)))
push!(out, "SET|S_notpositive_label|" * string(count(i -> REC[i].inst != POS, Sset)))

# H1
nE = length(Eset)
pa = count(i -> REC[i].dqb <= PAR["h1a_bound"], Eset)
pb = count(i -> REC[i].d9 <= PAR["h1b_bound"], Eset)
pc = count(i -> REC[i].d2 <= PAR["h1c_bound"], Eset)
v1 = nE < MINSET ? "not_decidable" :
     (pa / nE >= PAR["h1a_fraction"] && pb / nE >= PAR["h1b_fraction"] && pc / nE >= PAR["h1c_fraction"]) ? "supported" : "falsified"
push!(out, "H1|$nE|$pa|$pb|$pc|$v1")

# H2
B = Int(PAR["bootstrap_resamples"])
RW = [REC[i].R for i in Wset]
any(!isfinite, RW) && anom("non-finite R(n_H2) in W")
RW = filter(isfinite, RW)
if length(RW) >= 1
    medW = median_of(RW); loW, hiW = interval(bootmedians(RW, "H2/W", B))
else
    medW = loW = hiW = NaN
end
v2 = length(RW) < MINSET ? "not_decidable" : (medW >= PAR["h2_median"] && loW >= PAR["h2_interval"]) ? "supported" : "not_supported"
push!(out, "H2|$(length(RW))|$(repr(medW))|$(repr(loW))|$(repr(hiW))|$v2")

# H3
CW = [REC[i].nconv for i in Wset]; CS = [REC[i].nconv for i in Sset]
if length(CW) >= 1 && length(CS) >= 1
    mW, mS = median_of(CW), median_of(CS); ratio = mW / mS
    bw = bootmedians(CW, "H3/W", B); bs = bootmedians(CS, "H3/S", B)
    loR, hiR = interval(bw ./ bs)
else
    mW = mS = ratio = loR = hiR = NaN
end
v3 = (length(CW) < MINSET || length(CS) < MINSET) ? "not_decidable" :
     (ratio >= PAR["h3_ratio"] && loR > PAR["h3_interval"]) ? "supported" : "not_supported"
push!(out, "H3|$(length(CW))|$(length(CS))|$(repr(mW))|$(repr(mS))|$(repr(ratio))|$(repr(loR))|$(repr(hiR))|$v3")

# H4
h4ids = sort(unique(id for (id, _, _, _, _, _) in keys(ITEMS) if id > MAXPANEL))
Th, Tc = Int(PAR["h4_T"]), Int(PAR["h4_T_cross"])
h4 = Dict{String,Any}()
for id in h4ids
    K = parse(Float64, ORBPARAMS[id])
    tag = abs(K - pi / 2) <= 1e-9 ? "halfpi" : abs(K - 2pi) <= 1e-9 ? "2pi" : "other"
    f1 = getit(id, "L1", "Float64", Th); f5 = getit(id, "L5", "Float64", Th); d5 = getit(id, "L5", "Double64", Tc)
    h4[tag] = (L1 = lamat(f1, Th), L5 = lamat(f5, Th), cf = lamat(f5, Tc), cd = lamat(d5, Tc))
end
say("H4 starts: ", h4ids, " tags ", collect(keys(h4)))
if haskey(h4, "halfpi") && haskey(h4, "2pi")
    a, b = h4["halfpi"], h4["2pi"]
    okh = abs(a.L1 - PAR["h4_ref_halfpi"]) <= PAR["h4_tol_halfpi"] && abs(a.L5 - PAR["h4_ref_halfpi"]) <= PAR["h4_tol_halfpi"]
    ok2 = abs(b.L1 - PAR["h4_ref_2pi"]) <= PAR["h4_tol_2pi"] && abs(b.L5 - PAR["h4_ref_2pi"]) <= PAR["h4_tol_2pi"]
    v4 = (okh && ok2) ? "supported" : "did_not_reproduce"
    push!(out, "H4|$(repr(a.L1))|$(repr(a.L5))|$(repr(b.L1))|$(repr(b.L5))|$(repr(a.cf))|$(repr(a.cd))|$(repr(b.cf))|$(repr(b.cd))|$v4")
else
    anom("H4 starts not found for both K")
    push!(out, "H4|NaN|NaN|NaN|NaN|NaN|NaN|NaN|NaN|not_decidable")
end

# H5
live = [r for r in recs if !("excluded" in r.flags)]
sl = [r for r in live if get(PORB[r.id], "in_SL", false)]
h5a_n = length(sl); h5a_k = count(r -> r.ladder === true, sl)
sl_missing = count(r -> r.ladder === nothing, sl)
sl_missing > 0 && anom("H5a: $sl_missing S_L orbits without ladder data count as not agreeing")
h5b_n = length(live); h5b_k = count(r -> r.inst == UNR || r.stc == UNR, live)
nonec = [r for r in live if r.inst == NONE]
h5c_n = length(nonec); h5c_k = count(r -> r.stc == UNR, nonec)
dec = (h5a_n >= MINSET, h5b_n >= MINSET, h5c_n >= MINSET)
pass = (h5a_k / h5a_n >= PAR["h5a_min"], h5b_k / h5b_n <= PAR["h5b_max"], h5c_n > 0 && h5c_k / h5c_n <= PAR["h5c_max"])
v5 = all(dec) && all(pass) ? "supported" : any(dec[i] && !pass[i] for i in 1:3) ? "not_supported" : "not_decidable"
push!(out, "H5|$h5a_n|$h5a_k|$h5b_n|$h5b_k|$h5c_n|$h5c_k|$v5")

# cells (all panel orbits)
cells = Dict{Tuple{Int,String,String,String},Int}()
for r in recs
    o = PORB[r.id]
    key = r.cls == "control" ? (0, "controls", r.inst, r.stc) : (Int(o["k_index"]), r.cls, r.inst, r.stc)
    cells[key] = get(cells, key, 0) + 1
end
for key in sort(collect(keys(cells)))
    push!(out, "CELL|$(key[1])|$(key[2])|$(key[3])|$(key[4])|$(cells[key])")
end

# K = 0 closed form
a0 = pi / 4
for r in recs
    (r.cls == "T" && r.K == 0.0) || continue
    dev = 0.0; miss = false
    for m in VM
        it = getit(r.id, m, "Float64", N1)
        if it === nothing; miss = true; continue; end
        for (n, l) in zip(it.t, it.lam)
            cf = log(sqrt((cos(a0) + n * sin(a0))^2 + sin(a0)^2)) / n
            dev = max(dev, abs(l - cf))
        end
    end
    push!(out, "K0|$(r.id)|$(miss ? "NaN" : repr(dev))")
    # extra, not in the required format: deviation without L9, and L9 against its own closed form ln(sigma_max)/n
    dv = 0.0; d9 = 0.0
    for m in VM
        it = getit(r.id, m, "Float64", N1); it === nothing && continue
        for (n, l) in zip(it.t, it.lam)
            if m == "L9"
                d9 = max(d9, abs(l - log((n + sqrt(n^2 + 4.0)) / 2) / n))
            else
                dv = max(dv, abs(l - log(sqrt((cos(a0) + n * sin(a0))^2 + sin(a0)^2)) / n))
            end
        end
    end
    push!(out, "K0X|$(r.id)|$(repr(dv))|$(repr(d9))")
end

open(OUTF, "w") do io
    for l in out; println(io, l); end
end
say("wrote ", length(out), " lines to ", OUTF, " (", count(startswith("ORBIT|"), out), " ORBIT lines)")
say("anomalies: ", length(ANOM))
println("TOTAL WALL TIME ", round(time() - T0; digits = 1), " s")
