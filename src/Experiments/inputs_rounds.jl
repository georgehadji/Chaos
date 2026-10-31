# Panel inputs of Amendment D94 (R22): screening source, rounds, per-class sizes, fractions of D, size of D_B, Newton lattice, controls.

const _CLASS5 = ("T", "R", "C", "S", "X")
const _ZERO_UNDEFINED = ("C", "S", "X")       # classes that do not exist at K = 0

function _check_source(d)
    src = get(d, "screen_source", "lattice")
    (src isa AbstractString && src in ("lattice", "random")) || _bad(("screen_source",), "must be \"lattice\" or \"random\", got $(repr(src))")
    if src == "random"
        haskey(d, "screen_points") || _bad(("screen_points",), "is required with screen_source = \"random\"")
        return String(src), _int("screen_points", d["screen_points"])
    end
    haskey(d, "screen_points") && _bad(("screen_points",), "is refused with screen_source = \"lattice\"")
    return String(src), 0
end

function _check_rounds_keys(d)
    r_max = haskey(d, "r_max") ? _int("r_max", d["r_max"]) : 1
    for k in ("d_fraction", "db_size")
        haskey(d, k) && !haskey(d, "r_max") && _bad((k, "r_max"), "`$k` needs `r_max` (rounds are asked for explicitly)")
    end
    g = get(d, "gate_rounds", [1])
    (g isa AbstractVector && !isempty(g)) || _bad(("gate_rounds",), "must be a non-empty list of integers, got $(repr(g))")
    gate = [_int("gate_rounds", x) for x in g]
    allunique(gate) || _bad(("gate_rounds",), "has a repeated round")
    all(<=(r_max), gate) || _bad(("gate_rounds", "r_max"), "has a round above r_max = $r_max")
    db = haskey(d, "db_size") ? _int("db_size", d["db_size"]; min = 0) : 0
    ctl = get(d, "controls", "standard")
    (ctl isa AbstractString && ctl in ("none", "standard")) || _bad(("controls",), "must be \"none\" or \"standard\", got $(repr(ctl))")
    return (r_max = r_max, gate_rounds = gate, db_size = db, controls = String(ctl))
end

function _kindex(key, s, nK)
    i = s isa AbstractString ? tryparse(Int, s) : nothing
    (i !== nothing && 1 <= i <= nK) || _bad((key,), "K index $(repr(s)) is not an integer in 1 to $nK")
    return i
end

function _table_rows(key, t, nK)
    t isa AbstractDict || _bad((key,), "must be a table by K index and class")
    rows = Dict{Int, Dict{String, Any}}()
    for (ks, row) in t
        i = _kindex(key, ks, nK)
        row isa AbstractDict || _bad((key,), "the entry of K index $i must be a table by class")
        for c in keys(row)
            c in _CLASS5 || _bad((key,), "unknown class $(repr(c)) at K index $i")
        end
        rows[i] = Dict{String, Any}(String(c) => v for (c, v) in row)
    end
    return rows
end

# The size table: by K index and class; a class that exists at the K is required, an undefined class may be absent or 0.
function _check_quota(d, K, n_c)
    nK = length(K)
    exists(i, c) = !(K[i] == 0 && c in _ZERO_UNDEFINED)
    haskey(d, "quota") || return Dict{Int, Dict{String, Int}}(i => Dict(c => n_c for c in _CLASS5 if exists(i, c)) for i in 1:nK)
    rows = _table_rows("quota", d["quota"], nK)
    out = Dict{Int, Dict{String, Int}}()
    for i in 1:nK
        haskey(rows, i) || _bad(("quota",), "misses K index $i")
        out[i] = Dict{String, Int}()
        for c in _CLASS5
            if !haskey(rows[i], c)
                exists(i, c) && _bad(("quota",), "misses class $c at K index $i")
                continue
            end
            n = rows[i][c]
            (n isa Integer && !(n isa Bool) && n >= 0) || _bad(("quota",), "size of K index $i class $c must be an integer >= 0, got $(repr(n))")
            (n > 0 && !exists(i, c)) && _bad(("quota",), "class $c is undefined at K = 0 (K index $i) and cannot have size $n")
            out[i][c] = Int(n)
        end
    end
    return out
end

# A rational from an integer, a decimal text or float (its shortest decimal text), or `a/b`: exact, never through a binary64 product (I36).
function _ratio(where, v)
    v isa Bool && _bad(("d_fraction",), "$where must be a rational, not a bool")
    v isa Integer && return Rational{BigInt}(v)
    s = v isa AbstractFloat ? (isfinite(v) ? string(v) : "?") : v isa AbstractString ? String(strip(v)) : "?"
    m = match(r"^(-?\d+)\s*/\s*(\d+)$", s)
    if m !== nothing
        b = parse(BigInt, m.captures[2])
        b > 0 || _bad(("d_fraction",), "$where has a zero denominator")
        return Rational{BigInt}(parse(BigInt, m.captures[1]), b)
    end
    m = match(r"^(-?)(\d+)(?:\.(\d+))?(?:[eE]([+-]?\d+))?$", s)
    m === nothing && _bad(("d_fraction",), "$where must be a rational (decimal or a/b), got $(repr(v))")
    frac = m.captures[3] === nothing ? "" : string(m.captures[3])
    ex = m.captures[4] === nothing ? 0 : parse(Int, m.captures[4])
    r = Rational{BigInt}(parse(BigInt, string(m.captures[2]) * frac), BigInt(10)^length(frac)) *
        Rational{BigInt}(big(10)^max(ex, 0), big(10)^max(-ex, 0))
    return m.captures[1] == "-" ? -r : r
end

function _check_fractions(d, nK)
    haskey(d, "d_fraction") || return Dict{Int, Dict{String, Rational{BigInt}}}()
    rows = _table_rows("d_fraction", d["d_fraction"], nK)
    out = Dict{Int, Dict{String, Rational{BigInt}}}()
    for (i, row) in rows
        out[i] = Dict{String, Rational{BigInt}}()
        for (c, v) in row
            f = _ratio("fraction of K index $i class $c", v)
            (0 <= f <= 1) || _bad(("d_fraction",), "fraction of K index $i class $c must lie in [0, 1], got $(repr(v))")
            out[i][c] = f
        end
    end
    return out
end

# Amendment D94 keys of the panel input (R22); the base reader above handles every other key.
function _check_rounds(d, K, n_c, grid_side)
    src, pts = _check_source(d)
    rk = _check_rounds_keys(d)
    ngs = haskey(d, "newton_grid_side") ? _int("newton_grid_side", d["newton_grid_side"]) : grid_side
    src == "lattice" && haskey(d, "quota") && _bad(("quota", "screen_source"), "`quota` needs screen_source = \"random\"")
    return (screen_source = src, screen_points = pts, rk..., newton_grid_side = ngs, quota = _check_quota(d, K, n_c),
            d_fraction = _check_fractions(d, length(K)))
end
