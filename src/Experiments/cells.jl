# Panel cells: screening of the grid and the draw of the starts of one K (spec Algorithm 4 to 6, E12 to E15).
# Orbit records are TOML-ready dictionaries; ids, per-orbit seeds and subset flags are added once the order is known (E20).

const CLASSES = ("T", "R", "C", "S", "X")            # canonical order of the classes inside a K (E20)
const DRAW_ORDER = ("S", "C", "T", "R", "X")         # S before C because of the S-to-C exclusion (E14)
const _BIT_T, _BIT_C, _BIT_S = 0x01, 0x02, 0x04

# Cells whose screening run is held at once; at most the free host memory divided by four (spec Compute, chunking).
function _default_chunk(inp, n)
    per_cell = screen_bytes(inp.N_scr, inp.w) + 64
    return max(1, min(n, Int(Sys.free_memory() ÷ 4 ÷ per_cell)))
end

function _cell_flags(inp, ps, buf, z)
    lam, f = buf.lambda, 0x00
    if lam <= inp.theta_T && buf.S_tenth > 0 &&
       pred_T(lam, buf.S_tenth, buf.S_N, periodic_distance(z, ps); theta_T = inp.theta_T, g_T = inp.g_T, r_T = inp.r_T)
        f |= _BIT_T
    end
    if lam >= inp.theta_C && pred_C(lam, periodic_distance(z, ps; elliptic_only = true); theta_C = inp.theta_C, r_C = inp.r_C)
        f |= _BIT_C
    end
    pred_S(lam, buf.blocks; theta_S = inp.theta_S, w = inp.w, f_S = inp.f_S, L_S = inp.L_S) && (f |= _BIT_S)
    return f
end

"Class bits of every grid cell of one K (T, C, S); the screening values are dropped as soon as the bits are set (I4)."
function _screen_flags(K::Float64, inp, ps, centres, chunk, nt)
    n = length(centres)
    flags = zeros(UInt8, n)
    ch = chunk === nothing ? _default_chunk(inp, n) : max(1, Int(chunk))
    for lo in 1:ch:n
        m = min(n, lo + ch - 1) - lo + 1
        _par(() -> screen_buffer(inp.N_scr, inp.w), m, nt) do buf, k
            i = lo + k - 1
            screen_point!(buf, K, centres[i][1], centres[i][2], inp.N_scr, inp.w, inp.alpha, inp.tau)
            flags[i] = _cell_flags(inp, ps, buf, centres[i])
        end
    end
    return flags
end

_ktext(v) = v isa AbstractString ? String(v) : string(v)

# Fields every drawn orbit carries; `start` adds x0, p0 and their decimal texts.
function _drawn_base(inp, i, cls, t)
    seed = inp.seed
    return Dict{String, Any}("class" => cls, "control" => false, "control_kind" => "none", "system" => "standard", "k_index" => i,
                             "K_text" => _ktext(inp.raw["K"][i]), "K" => inp.K[i], "run_seed" => string(seed),
                             "cell_seed" => string(cell_stream_seed(seed, i, cls)), "draw_index" => t,
                             "in_SD" => false, "in_SL" => false)
end

function _with_start!(d, xb, pb)
    d["x0"], d["p0"] = Float64(xb), Float64(pb)
    d["x0_dec"], d["p0_dec"] = string(xb), string(pb)
    return d
end

_with_start!(d, x::Float64, p::Float64) = _with_start!(d, BigFloat(x), BigFloat(p))

# Draw n_c distinct grid cells among the candidate ids (ascending) for T, S or C; the start is the grid centre itself.
function _draw_screened(inp, i, cls, cands, centres)
    rng = cell_rng(inp.seed, i, cls)
    ids = draw_candidates(rng, cands, inp.n_c)
    orbs = [_with_start!(_drawn_base(inp, i, cls, t), centres[g + 1][1], centres[g + 1][2]) for (t, g) in enumerate(ids)]
    return (class = cls, qualifying = length(cands), orbits = orbs, grid_ids = ids)
end

# Elliptic periodic points of one K in table order, with (periodic id, point index) of each.
function _elliptic_points(ps, pids)
    pts = SVector{2, BigFloat}[]
    refs = Tuple{Int, Int}[]
    for (o, id) in zip(ps.orbits, pids)
        o.class === :elliptic || continue
        for (k, pt) in enumerate(o.points)
            push!(pts, pt)
            push!(refs, (id, k))
        end
    end
    return pts, refs
end

function _draw_R(inp, i, ps, pids)
    rng = cell_rng(inp.seed, i, "R")
    orbs = Dict{String, Any}[]
    if inp.K[i] == 0
        for t in 1:inp.n_c
            r = place_R0(rng, inp.q_max)
            d = _with_start!(_drawn_base(inp, i, "R", t), r.x, r.p)
            d["a"], d["b"] = r.a, r.b
            push!(orbs, d)
        end
        return (class = "R", qualifying = 0, orbits = orbs, grid_ids = Int[])
    end
    pts, refs = _elliptic_points(ps, pids)
    for t in 1:(isempty(pts) ? 0 : inp.n_c)
        r = place_R(rng, pts, inp.delta_R)
        d = _with_start!(_drawn_base(inp, i, "R", t), r.xb, r.pb)
        d["source_orbit"], d["source_point"], d["phi"] = refs[r.index][1], refs[r.index][2], r.phi
        push!(orbs, d)
    end
    return (class = "R", qualifying = length(pts), orbits = orbs, grid_ids = Int[])
end

function _draw_X(inp, i)
    rng = cell_rng(inp.seed, i, "X")
    orbs = map(1:inp.n_c) do t
        r = place_X(rng, inp.log10_r[1], inp.log10_r[2])
        d = _with_start!(_drawn_base(inp, i, "X", t), r.x, r.p)
        d["r"], d["phi"] = r.r, r.phi
        d
    end
    return (class = "X", qualifying = 0, orbits = orbs, grid_ids = Int[])
end

"Status and reason of a cell with `n` orbits drawn (spec Algorithm 6)."
function _status(cls, n, n_c)
    n == n_c && return "populated", ""
    n > 0 && return "short", "fewer qualifying grid cells than n_c"
    return "empty", cls == "R" ? "no elliptic periodic point" : "no qualifying grid cell"
end

_undefined(cls) = (class = cls, qualifying = 0, orbits = Dict{String, Any}[], grid_ids = Int[], undefined = true)

"""
Draw all cells of K index `i`; returns the cell records in canonical class order. `flags` are the screening bits of the grid.
"""
function _draw_k(inp, i, ps, pids, flags, centres)
    ids(bit) = [g - 1 for g in eachindex(flags) if flags[g] & bit != 0]
    zero_k = inp.K[i] == 0
    S = zero_k ? _undefined("S") : _draw_screened(inp, i, "S", ids(_BIT_S), centres)
    cC = zero_k ? Int[] : setdiff(ids(_BIT_C), S.grid_ids)
    C = zero_k ? _undefined("C") : _draw_screened(inp, i, "C", cC, centres)
    T = _draw_screened(inp, i, "T", ids(_BIT_T), centres)
    R = _draw_R(inp, i, ps, pids)
    X = zero_k ? _undefined("X") : _draw_X(inp, i)
    byclass = Dict("T" => T, "R" => R, "C" => C, "S" => S, "X" => X)
    return [byclass[c] for c in CLASSES]
end
