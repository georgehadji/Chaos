# Study p3-abstention (Amendment D95): the data of `analyze` (E61 j, k; E64; E65): the items of the four arms read once into an index, the identity check, the labels and drops, the sets, the
# strata, the verdicts and the outputs of the rules R_a and R_ab on the test orbits. Nothing is recomputed from orbit series.

p3t_theta(t) = (V = Symbol.(t.theta_voters), qC = t.theta_q_c, qR = t.theta_q_r, w = t.theta_w, v = t.theta_v)

# All items of the ok chunks: the ref items already read for the log (`st`) and the datasets of the closure arms; the index by content (a repeated content refuses).
function p3t_load_items(c, st, on_read)
    items = Any[it for it in values(st.items)]
    for ch in c.chunks
        (ch.ref || ch.status != "ok") && continue
        on_read === nothing || on_read(:datasets, ch.path)
        append!(items, _p2_chunk_items(read_chunk(ch.path), ch.ordinal))
    end
    return items, p2_index(items)
end

# E61 j: the stored start text and K of every item against the panel file (D89 E24).
function p3t_check_identity(items, panel_orbits)
    byid = Dict(Int(o["id"]) => o for o in panel_orbits)
    for it in items
        o = get(byid, it.orbit, nothing)
        o === nothing && continue
        a = it.attrs
        lbl = "chunk $(_chunk_name(it.chunk)), orbit $(it.orbit), item $(it.id)"
        haskey(a, "ic") && a["ic"] != "$(o["x0_dec"]),$(o["p0_dec"])" && _fail("$lbl: stored ic \"$(a["ic"])\" differs from the start of the panel file \"$(o["x0_dec"]),$(o["p0_dec"])\"")
        haskey(a, "params") && tryparse(Float64, String(a["params"])) != Float64(o["K"]) && _fail("$lbl: stored params \"$(a["params"])\" differs from the K $(o["K"]) of the panel file")
    end
end

# The verdict items of E65 that an orbit lacks (Float64, arm ind) and the cells of the ind chunks that failed twice.
function _p3t_missing_items(idx, p2, id)
    d1, F = p2.directions[1], "Float64"
    need = (("tg", p2.N_1, d1), ("fq1", p2.N_c, nothing), ("pa1", p2.N_c + 1, nothing), ("bk1", p2.N_c, nothing))
    return String[s for (s, h, d) in need if !haskey(idx, p2_key(id, s, F, h, 1, d))]
end

_p3t_failed_ind_cells(chunks) = Set(parse(Int, c.cell[2:end]) for c in chunks if "ind" in c.arms && c.status == "failed" && occursin(r"^k\d+$", c.cell))

# The orbits with the drift flag on any item read (E64 b); the counts of the round log use the flag of L1 and L5 only (E62 b), the difference is reported (I46).
_p3t_drifted(idx) = Set(k[1] for (k, it) in idx if "drift" in it.flags)

# The reason an orbit leaves E_test (E64 b), or "" to keep it. An item that E65 needs and an orbit lacks, outside a chunk that failed twice, refuses (E61 k).
function _p3t_reason(r, idx, p2, failed_ind, drifted)
    r.excluded && return "excluded"
    (r.drift || r.id in drifted) && return "drift"
    r.rho === :U && return "unresolved"
    miss = _p3t_missing_items(idx, p2, r.id)
    isempty(miss) && return ""
    r.kindex in failed_ind && return "no verdict data"
    _fail("item $(join(miss, ", ")) is missing from the store for orbit $(r.id), which the manifest does not list as excluded (E61 k)")
end

# The verdict tables of the orbits `ids` (E65): table R1 on the Float64 items, ten budgets per indicator of I_7.
p3t_verdict_tables(ctx, ids) = p3_verdicts(ctx, ids)

"""
    p3t_orbits(c, st, idx, theta) -> NamedTuple

The orbits of the test run (E64, E65): `rows` (the orbits of E_test: `id, rho, round, kindex, lref, lns, Ts, W2, S2, V, out, in_D`), `drops` (`id, reason, round, k_index`), `gate` (the K = 0
orbits of round 1 with their output, `nothing` without verdict items), the set sizes and the Double64 verdicts of D.
"""
function p3t_orbits(c, st, idx, theta)
    q, p2 = c.q, c.p2
    ctx = (p = p2, idx = idx)
    failed_ind, drifted = _p3t_failed_ind_cells(c.chunks), _p3t_drifted(idx)
    panel = Dict(Int(o["id"]) => o for o in c.panel_orbits)
    rows, drops, kept = NamedTuple[], NamedTuple[], NamedTuple[]
    for r in sort(st.rows; by = x -> x.id)
        why = _p3t_reason(r, idx, p2, failed_ind, drifted)
        isempty(why) ? push!(kept, r) : push!(drops, (id = r.id, reason = why, round = r.round, k_index = r.kindex))
    end
    V = p3t_verdict_tables(ctx, [r.id for r in kept])
    kP = _p3_kP(q)
    fs, lw, ls = p3_rat(q.f_s), p3_rat(q.lam_w2), p3_rat(q.lam_s2)
    for r in kept
        Vr = Dict{Symbol, Vector{Symbol}}(X => V[(r.id, X)] for X in P2_I7)
        C = r.rho === :C
        Ts = C && r.lns <= fs * r.lref
        push!(rows, (id = r.id, rho = r.rho, round = r.round, kindex = r.kindex, lref = r.lref, lns = r.lns, Ts = Ts, W2 = C && (Ts || r.lref <= lw), S2 = C && r.lref >= ls,
                     V = Vr, out = p3_rule_a(Vr, kP, theta), in_D = panel[r.id]["in_D"] === true))
    end
    return (rows = rows, drops = drops, gate = _p3t_gate(c, idx, ctx, theta, kP), vd = _p3t_double(rows, ctx, theta))
end

# The gate set (E64 c): the K = 0 orbits of round 1; their number must equal `gate_size`. An orbit without verdict items has the output `nothing` (HA is then not decidable, P100).
function _p3t_gate(c, idx, ctx, theta, kP)
    ids = sort!(Int[Int(o["id"]) for o in c.panel_orbits if Float64(o["K"]) == 0 && Int(o["round"]) == 1])
    length(ids) == c.q.gate_size || _fail("the gate set has $(length(ids)) orbits, the parameter table says gate_size = $(c.q.gate_size)")
    have = [id for id in ids if isempty(_p3t_missing_items(idx, c.p2, id))]
    V = p3t_verdict_tables(ctx, have)
    out = Union{Symbol, Nothing}[id in have ? p3_rule_a(Dict{Symbol, Vector{Symbol}}(X => V[(id, X)] for X in P2_I7), kP, theta) : nothing for id in ids]
    return (ids = ids, out = out, missing_ids = [id for id in ids if !(id in have)])
end

# The Double64 verdicts of the voters of V* for the orbits of E_test with the flag in_D that carry the four Double64 items (E64 c): the set D.
function _p3t_double(rows, ctx, theta)
    out = Dict{Int, Dict{Symbol, Vector{Symbol}}}()
    for r in rows
        r.in_D || continue
        vals = [p2_values(ctx, r.id, X; prec = "Double64") for X in theta.V]
        any(isnothing, vals) && continue
        out[r.id] = Dict{Symbol, Vector{Symbol}}(X => Symbol[p2_verdict(X, :R1, v, ctx.p) for v in vs] for (X, vs) in zip(theta.V, vals))
    end
    return out
end

"The orbits of D_B (E64 c): the orbits that carry an item of the arm chk (the L5 item in Double64, in a chunk of that arm); their flag in_D must be set."
function p3t_dB(c, idx)
    p2 = c.p2
    ids = sort!(unique(Int[k[1] for (k, it) in idx if k[2] == "L5" && k[3] == "Double64" && "chk" in it.arms]))
    flag = Dict(Int(o["id"]) => o["in_D"] === true for o in c.panel_orbits)
    for id in ids
        get(flag, id, false) || _fail("orbit $id carries a check item (arm chk) but the flag in_D of the panel file is not set (E64 c)")
    end
    return ids
end

# R_ab (E52) on D: the output of every orbit of D, the counts of the nine pairs and the overrides (T9).
function p3t_rab(rows, vd, theta, kP)
    pairs_ = Dict{String, Int}()
    over, ids, outs = 0, Int[], Dict{Int, Symbol}()
    for r in rows
        haskey(vd, r.id) || continue
        ab = p3_rule_ab(r.V, vd[r.id], kP, theta)
        key = "$(r.out)|$ab"
        pairs_[key] = get(pairs_, key, 0) + 1
        r.out !== :unresolved && ab === :unresolved && (over += 1)
        push!(ids, r.id); outs[r.id] = ab
    end
    return (pairs = pairs_, overrides = over, orbits = ids, out = outs)
end
