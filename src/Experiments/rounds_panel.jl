# Panel builder with the random screening source and the draw in rounds (Amendment D94: E32 to E37, R23, R24).
# Orbit records are TOML-ready dictionaries, as in cells.jl; a round cell is (round, K index, class).

const ROUND_STREAM_RULE = "first 8 bytes, big-endian, of SHA-256 over the UTF-8 text \"chaoslab-screen-1|<seed>|<K index>\" (screening points), " *
    "\"chaoslab-panel-1|<seed>|<K index>|<class>|<round>\" (round cells), \"chaoslab-d-1|<seed>|<round>|<K index>|<class>\" (subset D) and " *
    "\"chaoslab-db-1|<seed>|<rounds run>|<K index>\" (subset D_B, drawn by the run plan); StableRNGs.StableRNG (package version in stablerngs_version)"

# The rounds in which the cells of K index i are drawn: the gate rounds at K = 0, every round otherwise (E35).
_active(inp, i) = inp.K[i] == 0 ? inp.gate_rounds : collect(1:inp.r_max)

# Fields every orbit of a round cell carries on top of the base ones; `in_D` follows E37.
function _stamp_round!(inp, orbs, i, c, r)
    f = get(get(inp.d_fraction, i, Dict{String, Rational{BigInt}}()), c, zero(Rational{BigInt}))
    sel = Set(d_indices(inp.seed, r, i, c, length(orbs), f))
    for d in orbs
        d["round"], d["cell_seed"] = r, string(round_stream_seed(inp.seed, i, c, r))
        d["in_D"] = d["draw_index"] in sel
    end
    return orbs
end

function _record(r, i, c, q, pool, orbs, status, reason)
    return (round = r, k_index = i, class = c, q = q, pool = pool, orbits = orbs, status = status, reason = reason)
end

# A pool class: the orbits start at the screening points drawn in `cell` (point id kept).
function _pool_record(inp, i, r, cell, pts)
    orbs = [_with_start!(_drawn_base(inp, i, cell.class, t), pts[j][1], pts[j][2]) for (t, j) in enumerate(cell.drawn)]
    for (d, j) in zip(orbs, cell.drawn)
        d["point_id"] = j
    end
    return _record(r, i, cell.class, cell.q, cell.pool, _stamp_round!(inp, orbs, i, cell.class, r), cell.status, cell.reason)
end

function _free_record(inp, i, r, c, q, ps, pids)
    rng = StableRNG(round_stream_seed(inp.seed, i, c, r))
    rec = q == 0 ? (qualifying = 0, orbits = Dict{String, Any}[]) :
          c == "R" ? _draw_R(inp, i, ps, pids; rng = rng, n = q) : _draw_X(inp, i; rng = rng, n = q)
    status, reason = _round_status(c, length(rec.orbits), q)
    return _record(r, i, c, q, rec.qualifying, _stamp_round!(inp, rec.orbits, i, c, r), status, reason)
end

# One round cell of K index i (E34, E35): undefined at K = 0, a pool class from its draw, or an R or X draw.
function _round_cell(inp, r, i, c, pools, pts, ps, pids)
    zero_k = inp.K[i] == 0
    if zero_k && c in _ZERO_UNDEFINED
        return _record(r, i, c, 0, 0, Dict{String, Any}[], "undefined", "class undefined at K = 0")
    elseif c in ("T", "C", "S")
        return _pool_record(inp, i, r, only(x for x in pools if x.round == r && x.class == c), pts)
    end
    q = (zero_k && !(r in inp.gate_rounds)) ? 0 : get(inp.quota[i], c, 0)
    return _free_record(inp, i, r, c, q, ps, pids)
end

# All round cells in the order of E36: round, K index, class (T, R, C, S, X).
function _round_records(inp, psets, pids, bits, pts)
    pools = [draw_pool_rounds(inp.seed, i, bits[i]; quota = inp.quota[i], r_max = inp.r_max, active = _active(inp, i)) for i in eachindex(inp.K)]
    return [_round_cell(inp, r, i, c, pools[i], pts[i], psets[i], pids[i]) for r in 1:inp.r_max for i in eachindex(inp.K) for c in CLASSES]
end

# R24 and I29: a point id is drawn once per K and carries the bit of its class; no two orbits share a start.
function _check_round_draws(recs, bits, orbs)
    for i in eachindex(bits)
        ids = [d["point_id"] for rc in recs if rc.k_index == i for d in rc.orbits if haskey(d, "point_id")]
        allunique(ids) || error("build_panel: a screening point of K index $i is drawn twice")
    end
    for rc in recs, d in rc.orbits
        haskey(d, "point_id") && bits[rc.k_index][d["point_id"]] & _POOL_BIT[rc.class] == 0 && error("build_panel: a drawn point lacks the bit of its class")
    end
    seen = Dict{Tuple{Float64, Float64}, Int}()
    for d in orbs
        d["control"] && continue
        k = (d["x0"], d["p0"])
        haskey(seen, k) && throw(ArgumentError("orbits $(seen[k]) and $(d["id"]) have the same start $k (I29)"))
        seen[k] = d["id"]
    end
end

# Maximal runs `[first, last]` of consecutive point ids that carry `bit`.
function _runs(bits, bit)
    runs = Vector{Int}[]
    for j in eachindex(bits)
        bits[j] & bit == 0 && continue
        !isempty(runs) && runs[end][2] == j - 1 ? (runs[end][2] = j) : push!(runs, [j, j])
    end
    return runs
end

# The pool sidecar (R23 item 3): per K the number of points, the SHA-256 of the class bit vector and the bits as runs of point ids. No screening value.
function _pool_sidecar(bits)
    k = Dict{String, Any}(string(i) => Dict{String, Any}("points" => length(b), "bits_sha256" => bytes2hex(sha256(b)),
                                                         "T" => _runs(b, _BIT_T), "C" => _runs(b, _BIT_C), "S" => _runs(b, _BIT_S))
                          for (i, b) in enumerate(bits))
    return Dict{String, Any}("k" => k)
end

_tab(t, f) = Dict{String, Any}(string(i) => Dict{String, Any}(c => f(v) for (c, v) in row) for (i, row) in t)

# Manifest keys of the amendment: source, sizes, rounds, the stream rule, the pool counts and one record per round cell.
function _round_manifest(inp, recs, bits)
    pools = [Dict{String, Any}("k_index" => i, "class" => c, "count" => count(x -> x & _POOL_BIT[c] != 0, b))
             for (i, b) in enumerate(bits) for c in ("T", "C", "S")]
    cells = [Dict{String, Any}("round" => rc.round, "k_index" => rc.k_index, "class" => rc.class, "pool" => rc.pool, "q" => rc.q,
                               "drawn" => length(rc.orbits), "status" => rc.status, "reason" => rc.reason) for rc in recs]
    return Dict{String, Any}("screen_source" => inp.screen_source, "screen_points" => inp.screen_points, "r_max" => inp.r_max,
        "gate_rounds" => inp.gate_rounds, "db_size" => inp.db_size, "quota" => _tab(inp.quota, identity),
        "d_fraction" => _tab(inp.d_fraction, f -> "$(numerator(f))/$(denominator(f))"), "stream_rule" => ROUND_STREAM_RULE,
        "pools" => pools, "round_cells" => cells)
end

# The cell records of the base manifest, one per (K index, class), summed over the rounds.
function _round_totals(inp, recs, bits)
    rows = Dict{String, Any}[]
    for i in eachindex(inp.K), c in CLASSES
        rs = [rc for rc in recs if rc.k_index == i && rc.class == c]
        n, q = sum(length(rc.orbits) for rc in rs), sum(rc.q for rc in rs)
        undefined = all(rc.status == "undefined" for rc in rs)
        status, reason = undefined ? ("undefined", "class undefined at K = 0") : _round_status(c, n, q)
        qual = c in ("T", "C", "S") ? count(x -> x & _POOL_BIT[c] != 0, bits[i]) : 0
        push!(rows, Dict{String, Any}("k_index" => i, "class" => c, "qualifying" => qual, "drawn" => n, "status" => status, "reason" => reason))
    end
    return rows
end

# Screening of every K on the random points; the points and the screening values are dropped as soon as the bits are set (I34).
function _screen_random(inp, psets, chunk, nt)
    pts, bits = Vector{SVector{2, Float64}}[], Vector{UInt8}[]
    for i in eachindex(inp.K)
        p = random_screen_points(inp.seed, i, inp.screen_points)
        push!(pts, p)
        push!(bits, _screen_flags(inp.K[i], inp, psets[i], p, chunk, nt))
    end
    return pts, bits
end

# `_assemble` of the random source: everything computed at the working precision, nothing written.
function _assemble_rounds(path, pr, nt, chunk)
    inp = pr.inputs
    psets = [_search_k(inp, i, nt) for i in eachindex(inp.K)]
    rows, pids = _periodic_rows(psets)
    pts, bits = _screen_random(inp, psets, chunk, nt)
    recs = _round_records(inp, psets, pids, bits, pts)
    controls = _controls(inp, psets)
    orbs, sub = _number_orbits!(inp, recs, controls)
    _check_round_draws(recs, bits, orbs)
    man = _manifest(path, pr, psets, rows, _round_totals(inp, recs, bits), orbs, sub, filter(d -> d["control"], orbs), _root())
    merge!(man, _round_manifest(inp, recs, bits))
    return Dict{String, Any}("manifest" => man, "orbit" => orbs, "pools_sidecar" => _pool_sidecar(bits))
end
