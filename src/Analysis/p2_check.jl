# Study p2-indicators (Amendment D89): the input checks of E36, before any output. Every failure is an ArgumentError that names the cause.
# The checks of the file level (hashes, manifest, chunk attributes) come first; those of the stored items (E36 g to k, n, p) follow the reading of the store.

const P2_SLUG = "p2-indicators"
const P1_SLUG = "p1-le-benchmark"

function _p2_check_prereg(prereg)
    isfile(prereg) || _fail("preregistration $prereg does not exist")
    try
        return check_prereg(prereg).header_sha256
    catch e
        e isa ArgumentError ? _fail("preregistration check failed: " * e.msg) : rethrow()
    end
end

# The hash rows of the manifest (header, start, end) all equal the hash of the file (E36 b).
function _p2_check_manifest(man, ph, p)
    _check_manifest(man)
    man["study"] == P2_SLUG || _fail("manifest study \"$(man["study"])\" differs from the preregistration slug \"$P2_SLUG\"")
    pr = get(man, "prereg", Dict())
    for k in ("header", "sha256_start", "sha256_end")
        get(pr, k, "") == ph || _fail("manifest prereg $k hash $(get(pr, k, "absent")) differs from the preregistration header hash $ph")
    end
    get(get(man, "git", Dict()), "clean", false) === true || _fail("the manifest says the git tree of the run was not clean (git.clean is not true)")
    ss = get(man, "sessions", Any[])
    all(s -> get(s, "clean", false) === true, ss) || _fail("a session of the manifest was not clean")
    isempty(get(man, "control_stops", Any[])) || _fail("the manifest lists control stops: $(length(man["control_stops"]))")
    pc = get(man, "plan_counts", Dict())
    length(man["chunks"]) == get(pc, "chunks", -1) || _fail("the manifest lists $(length(man["chunks"])) chunks but plan_counts.chunks is $(get(pc, "chunks", "absent"))")
    get(man, "plan_sha256", "") == p.plan_sha256 || _fail("plan hash $(get(man, "plan_sha256", "absent")) of the manifest differs from plan_sha256 of the parameter table")
    return pc
end

# The panel file against the manifest (reused from P1) and the table; the build hash of the manifest against the table (E36 e).
function _p2_check_panel(panel, man, p)
    h = _check_panel(panel, man)
    h == p.panel_sha256 || _fail("panel file hash $h differs from panel_sha256 of the parameter table")
    get(get(man, "panel", Dict()), "build_hash", "") == p.panel_build_sha256 || _fail("panel build hash of the manifest differs from panel_build_sha256 of the parameter table")
    return h
end

# The chunk attributes beyond the P1 checks: plan hash, chunk status, study (E36 d, f).
function _p2_check_chunk_attrs(man, store, p)
    for c in man["chunks"]
        a = Experiments.read_chunk_attrs(chunk_path(store, c["ordinal"]))
        name = _chunk_name(c["ordinal"])
        get(a, "plan_sha256", p.plan_sha256) == p.plan_sha256 || _fail("plan_sha256 of $name differs from the plan hash of the parameter table")
        get(a, "status", "ok") == "ok" || _fail("chunk $name has status \"$(a["status"])\", not ok")
        get(a, "study", P2_SLUG) == P2_SLUG || _fail("chunk $name belongs to study \"$(a["study"])\"")
    end
end

"""
    p2_check_inputs(prereg, manifest, panel, params, store, p1_params) -> (p, man, hashes, raw)

The checks of E36 a to f and l that need the files only: preregistration hash, parameter table, manifest, panel file, chunk files. Nothing is written.
"""
function p2_check_inputs(prereg, manifest, panel, params_file, store, p1_params)
    ph = _p2_check_prereg(prereg)
    p = p2_read_params(params_file; prereg = prereg, p1_params = p1_params)
    isfile(manifest) || _fail("manifest $manifest does not exist")
    man = read_json(manifest)
    _p2_check_manifest(man, ph, p)
    pnl = _p2_check_panel(panel, man, p)
    _check_store(man, store, ph, pnl)
    _p2_check_chunk_attrs(man, store, p)
    files = [prereg, manifest, panel, params_file]
    hashes = (prereg_header_sha256 = ph, panel_sha256 = pnl, manifest_sha256 = bytes2hex(sha256(read(manifest))), params_sha256 = bytes2hex(sha256(read(params_file))),
              plan_sha256 = p.plan_sha256)
    return (p = p, man = man, hashes = hashes, raw = Dict(f => bytes2hex(sha256(read(f))) for f in files))
end

# ---------------------------------------------------------------- the stored items

_lbl(it) = "$(_chunk_name(it.chunk)), orbit $(it.orbit), spec $(it.id)"

# Recorded steps (E26): the count matches the attribute, the last step is the horizon, and every orbit of a chunk has the grid of the first (E36 i).
function _chk_steps(it, first_t)
    n = length(it.t)
    Int(get(it.attrs, "steps_recorded", n)) == n || _fail("$(_lbl(it)): $n recorded steps but steps_recorded is $(it.attrs["steps_recorded"]); a recorded step is missing")
    it.t == first_t || _fail("$(_lbl(it)): the recorded step grid differs from the grid of the first orbit of the chunk; a recorded step is missing")
    (n > 0 && issorted(it.t) && it.t[end] == it.horizon) || _fail("$(_lbl(it)): the last recorded step is not the horizon $(it.horizon)")
end

# Values of the stored point and the drift flag against the recomputation from the points (E36 j, k).
function _chk_points(it, p)
    all(isfinite, it.state) || _fail("non-finite value in the stored points of $(_lbl(it))")
    it.lam === nothing || all(isfinite, it.lam) || _fail("non-finite value in the exponents of $(_lbl(it))")
    pp = it.winding === nothing ? view(it.state, :, 2) : it.state[:, 2] .+ 2pi .* it.winding[:, 2]
    drifted = any(x -> abs(x) > p.drift_bound, pp)
    drifted == ("drift" in it.flags) || _fail("drift flag of $(_lbl(it)) differs from the recomputation from the stored points")
end

const _P2_SECOND_BIT = Dict("frequency" => "degenerate", "phase_advance" => "undefined window", "birkhoff" => "zero")

# A table of a window kind: finite unless the flag of the row says undefined (table E31 of D86); `complete` of a phase advance or Birkhoff item is 1.
function _chk_window(it)
    sec = _P2_SECOND_BIT[it.kind]
    (length(it.bits) >= 2 && it.bits[2] == sec) || _fail("$(_lbl(it)): the second flag name is not \"$sec\"; the bit values cannot be read")
    for r in axes(it.final, 1)
        w = it.flagwords[r]
        any(isinf, view(it.final, r, :)) && _fail("non-finite (infinite) value in row $(it.rows[r]) of $(_lbl(it))")
        (w & 2) != 0 && continue
        for c in axes(it.final, 2)
            isnan(it.final[r, c]) || continue
            it.kind == "birkhoff" && it.fcolumns[c] == "reldig" && (w & 4) != 0 && continue
            _fail("non-finite value in row $(it.rows[r]), column $(it.fcolumns[c]) of $(_lbl(it)) without the flag $sec")
        end
        it.kind == "frequency" && (w & 4) != 0 && _fail("row $(it.rows[r]) of $(_lbl(it)) has the flag failed")
    end
    j = findfirst(==("complete"), it.scolumns)
    j === nothing || it.scalars[j] == 1 || _fail("complete of $(_lbl(it)) is not 1")
end

function _chk_item(it, first_t, p)
    _chk_steps(it, first_t)
    _chk_points(it, p)
    it.series === nothing || all(isfinite, it.series) || _fail("non-finite value in the series of $(_lbl(it))")
    it.final === nothing || _chk_window(it)
end

# Identity of the orbits (E24): the stored texts against the panel file and the grid tables.
function _chk_identity(it, panel_by_id, grid_by_id)
    a = it.attrs
    if haskey(grid_by_id, it.orbit)
        g = grid_by_id[it.orbit]
        tryparse(Float64, String(a["params"])) == parse(Float64, g.K) || _fail("$(_lbl(it)): stored params \"$(a["params"])\" differs from the K \"$(g.K)\" of its h4 table")
        a["ic"] == "$(g.x0),$(g.p0)" || _fail("$(_lbl(it)): stored ic \"$(a["ic"])\" differs from \"$(g.x0),$(g.p0)\" of its h4 table")
    elseif haskey(panel_by_id, it.orbit)
        o = panel_by_id[it.orbit]
        a["ic"] == "$(o["x0_dec"]),$(o["p0_dec"])" || _fail("$(_lbl(it)): stored ic \"$(a["ic"])\" differs from the start of the panel file \"$(o["x0_dec"]),$(o["p0_dec"])\"")
        haskey(o, "K") && (tryparse(Float64, String(a["params"])) == Float64(o["K"]) ||
                           _fail("$(_lbl(it)): stored params \"$(a["params"])\" differs from the K $(o["K"]) of the panel file"))
    else
        _fail("$(_lbl(it)): the orbit is neither a panel orbit nor a grid orbit")
    end
end

# The flags of the manifest per chunk (the orbits with the flag in at least one row) against the stored words (E36 p).
function _chk_flag_counts(byc, man)
    for c in man["chunks"]
        its = get(byc, Int(c["ordinal"]), NamedTuple[])
        counts = Dict{String, Int}()
        for it in its, (j, nm) in enumerate(it.bits)
            it.flagwords === nothing && continue
            any(w -> (w & (Int64(1) << (j - 1))) != 0, it.flagwords) && (counts[nm] = get(counts, nm, 0) + 1)
        end
        want = Dict{String, Int}(k => Int(v) for (k, v) in get(c, "indicator_flags", Dict()) if v != 0)
        counts == want || _fail("flag counts of $(_chunk_name(c["ordinal"])) recomputed from the stored words $(sort(collect(counts))) differ from indicator_flags of the manifest $(sort(collect(want)))")
    end
end

"""
    p2_check_items(items, c, panel_orbits, grid) -> Dict

The item-level checks of E36 h to k: steps, finite values, drift flags, identity of the orbits. Returns the items grouped by chunk (the flag counts of E36 p follow the coverage check).
"""
function p2_check_items(items, c, panel_orbits, grid)
    byc = Dict{Int, Vector{NamedTuple}}()
    for it in items
        push!(get!(byc, it.chunk, NamedTuple[]), it)
    end
    pb = Dict(Int(o["id"]) => o for o in panel_orbits)
    gb = Dict(g.id => g for g in grid)
    for its in values(byc), it in its
        _chk_item(it, its[1].t, c.p)
        _chk_identity(it, pb, gb)
    end
    return byc
end
