# Study p2-indicators (Amendment D89): the items each orbit must have (E25, E36 g), the value records read at the recorded steps of E26, and the orbit table
# with the reference class and lambda_ref (E23).

# The rows of E25 as requirements: (spec id, precision, horizon, direction, orbits that must have the item).
function p2_requirements(p, panel_orbits, grid)
    P = Int[Int(o["id"]) for o in panel_orbits]
    G = Int[g.id for g in grid]
    sd = Int[Int(o["id"]) for o in panel_orbits if o["in_SD"] === true]
    sl = Int[Int(o["id"]) for o in panel_orbits if o["in_SL"] === true]
    d1, F = p.directions[1], "Float64"
    req(id, prec, hor, dir, orbs) = (id = id, precision = prec, horizon = hor, direction = dir, orbits = orbs)
    out = [req("L1", F, p.N_1, d1, vcat(P, G)), req("L5", F, p.N_1, d1, G), req("tg", F, p.N_1, d1, vcat(P, G)), req("pr", F, p.N_c, nothing, vcat(sl, G)),
           req("st", F, p.N_c, d1, vcat(P, G)), req("fq1", F, p.N_c, nothing, vcat(P, G)), req("pa1", F, p.N_c + 1, nothing, vcat(P, G)), req("bk1", F, p.N_c, nothing, vcat(P, G))]
    for d in p.directions, id in ("L1", "tg")
        push!(out, req(id, F, p.N_dir, d, sd))
    end
    append!(out, [req("L1", "Double64", p.N_dir, d1, sl), req("tg", "Double64", p.N_dir, d1, sl), req("fq1", "Double64", p.N_c, nothing, sl),
                  req("bk1", "Double64", p.N_c, nothing, sl), req("pa1", "Double64", p.N_c + 1, nothing, sl), req("fq1_raw", F, p.N_c, nothing, sl),
                  req("fq1_qp", F, p.N_c, nothing, sl), req("pa_planar", F, p.N_c + 1, nothing, sl)])
    return out
end

# The orbits the manifest lists as not applicable for a spec (all chunks).
function _na_lists(man)
    na = Dict{String, Set{Int}}()
    for c in man["chunks"], (s, ids) in get(c, "not_applicable", Dict())
        union!(get!(na, String(s), Set{Int}()), Int.(ids))
    end
    return na
end

"""
    p2_check_coverage(idx, man, reqs) -> Dict{String, Vector{Int}}

An item of E25 missing for an orbit that the manifest does not list in `not_applicable`, or listed there and present, refuses (E36 g). Returns the omissions per spec.
"""
function p2_check_coverage(idx, man, reqs)
    na = _na_lists(man)
    omitted = Dict{String, Vector{Int}}()
    for r in reqs, o in r.orbits
        has = haskey(idx, p2_key(o, r.id, r.precision, r.horizon, 1, r.direction))
        listed = o in get(na, r.id, Set{Int}())
        what = "orbit $o, spec $(r.id), $(r.precision), horizon $(r.horizon), direction $(r.direction)"
        !has && !listed && _fail("item missing from the store: $what; the manifest does not list it in not_applicable")
        has && listed && _fail("the manifest lists $what in not_applicable but the store has its item")
        has || push!(get!(omitted, r.id, Int[]), o)
    end
    return Dict(k => sort!(unique(v)) for (k, v) in omitted)
end

# ---------------------------------------------------------------- value records (E26, E27)

_P2_SPEC = Dict(:LE => "L1", :SALI => "tg", :FLI => "tg", :MEGNO => "tg", :RLI => "pr", :E1 => "fq1", :E2 => "pa1", :WBA => "bk1", :STR => "st")
_P2_SERIES = (:LE, :SALI, :FLI, :MEGNO, :STR)

# The item of indicator X for an orbit: precision, horizon and direction default to the primary item (E25).
function _p2_find(ctx, id, X; prec = "Float64", horizon = nothing, dir = nothing)
    p = ctx.p
    hor = horizon !== nothing ? horizon : X in (:LE, :SALI, :FLI, :MEGNO) ? p.N_1 : X === :E2 ? p.N_c + 1 : p.N_c
    d = X in _P2_SERIES ? (dir === nothing ? p.directions[1] : dir) : nothing
    return get(ctx.idx, p2_key(id, _P2_SPEC[X], prec, hor, 1, d), nothing)
end

function _p2_col(it, names, what)
    j = findfirst(==(what), names)
    j === nothing && _fail("$(_lbl(it)) has no column $what")
    return j
end

_p2_row(it, n::Integer, what) = (i = searchsortedfirst(it.t, n); (i <= length(it.t) && it.t[i] == n) ? i : _fail("recorded step $n is missing in $(_lbl(it)) ($what)"))

function _vals_le(it, p)
    out = NamedTuple[]
    for b in p.budgets
        i = _p2_row(it, p2_budget_step(it.t, b), "budget $b")
        n = it.t[i]
        ip = _p2_row(it, p2_nearest_step(it.t, n / 10), "n-prime of $n")
        push!(out, (n = n, lambda = it.lam[i], r = p2_ratio(p2_growth(:LE, n, it.lam[i]), p2_growth(:LE, it.t[ip], it.lam[ip]))))
    end
    return out
end

_P2_TGCOL = Dict(:SALI => "SALI.L", :FLI => "FLI.fli", :MEGNO => "MEGNO.Ybar20_1")

function _vals_tg(it, p, X)
    j = _p2_col(it, it.columns, _P2_TGCOL[X])
    js = X === :SALI ? _p2_col(it, it.columns, "SALI.sali") : j
    out = NamedTuple[]
    for b in p.budgets
        i = _p2_row(it, p2_budget_step(it.t, b), "budget $b")
        n = it.t[i]
        ip = _p2_row(it, p2_nearest_step(it.t, n / 10), "n-prime of $n")
        r = p2_ratio(p2_growth(X, n, it.series[i, j]), p2_growth(X, it.t[ip], it.series[ip, j]))
        push!(out, X === :SALI ? (n = n, r = r, sali = it.series[i, js]) : X === :FLI ? (n = n, r = r, fli = it.series[i, j]) : (n = n, r = r, ybar = it.series[i, j]))
    end
    return out
end

function _vals_col(it, p, name, field)
    j = _p2_col(it, it.columns, name)
    out = NamedTuple[]
    for b in p.budgets
        i = _p2_row(it, p2_budget_step(it.t, b), "budget $b")
        push!(out, NamedTuple{(:n, field)}((it.t[i], it.series[i, j])))
    end
    return out
end

# A window row `W=<W>` of a table (E26): the window is b / 2 and must be an element of the windows.
function _vals_window(it, p, X)
    j = _p2_col(it, it.fcolumns, X === :WBA ? "dig" : "D")
    out = NamedTuple[]
    for b in p.budgets
        p2_window_row(p.windows, b)
        W = b ÷ 2
        if X === :WBA
            rs = [r for r in eachindex(it.rows) if startswith(it.rows[r], "W=$W ")]
            length(rs) == 4 || _fail("window W=$W of $(_lbl(it)) has $(length(rs)) rows, not four")
            push!(out, (W = W, dig = Float64[it.final[r, j] for r in rs], flags = Int64[it.flagwords[r] for r in rs]))
        else
            r = findfirst(==("W=$W"), it.rows)
            r === nothing && _fail("window row W=$W is missing in $(_lbl(it))")
            push!(out, (W = W, D = it.final[r, j], flags = it.flagwords[r]))
        end
    end
    return out
end

"""
    p2_values(ctx, id, X; prec, horizon, dir) -> Union{Nothing, Vector}

The value records of indicator `X` (an indicator of I_8 or `:STR`) for orbit `id` at the ten budgets, read at the recorded steps of E26; `nothing` if the orbit has no item.
"""
function p2_values(ctx, id, X::Symbol; kw...)
    it = _p2_find(ctx, id, X; kw...)
    it === nothing && return nothing
    p = ctx.p
    X === :LE && return _vals_le(it, p)
    X in (:SALI, :FLI, :MEGNO) && return _vals_tg(it, p, X)
    X === :RLI && return _vals_col(it, p, "RLI.RLI", :rli)
    X === :STR && return _vals_col(it, p, "stretching.D", :D)
    return _vals_window(it, p, X)
end

p2_step_of(v) = hasproperty(v, :n) ? v.n : v.W

# ---------------------------------------------------------------- the orbit table

# lambda_1 at a recorded step of the L1 item of an orbit.
_lam_at(ctx, id, n) = (it = ctx.idx[p2_key(id, "L1", "Float64", ctx.p.N_1, 1, ctx.p.directions[1])]; it.lam[_p2_row(it, n, "lambda_1")])

# The grid orbit label by the reduced rule of E23 b, with the lifted first coordinate (P19).
function _grid_label(ctx, g)
    p = ctx.p
    it = ctx.idx[p2_key(g.id, "L5", "Float64", p.N_1, 1, p.directions[1])]
    h = p.N_1 ÷ 2
    xl(n) = (i = _p2_row(it, n, "grid label"); it.state[i, 1] + 2 * pi * it.winding[i, 1])
    x0 = parse(Float64, first(split(String(it.attrs["ic"]), ',')))
    r1, r2 = p2_rotation_numbers(x0, xl(h), xl(p.N_1), h, p.N_1)
    lam5 = it.lam[_p2_row(it, p.N_1, "grid label")]
    return p2_grid_label((lam5 = lam5, lam5_tenth = it.lam[_p2_row(it, p.N_1 ÷ 10, "grid label")], lam1 = _lam_at(ctx, g.id, p.N_1), rho1 = r1, rho2 = r2), p)
end

"""
    p2_orbit_table(ctx, labels) -> (orbits, grid_rows)

One record per orbit of P and G for `p2_sets` (E23, E24, E29): the reference class, lambda_ref = lambda_1(N_1) and lambda_ns = lambda_1(n_s) of the P2 L1 item (P13).
`labels` maps a panel orbit id to the label of the P1 analysis. The flags of an orbit are the union over its items.
"""
function p2_orbit_table(ctx, labels)
    p = ctx.p
    fl = Dict{Int, Set{String}}()
    for it in ctx.items
        union!(get!(fl, it.orbit, Set{String}()), it.flags)
    end
    flags(id) = sort!(collect(get(fl, id, Set{String}())))
    out = NamedTuple[]
    for o in ctx.panel
        id = Int(o["id"])
        lab = labels[id]
        push!(out, (id = id, source = :panel, kindex = Int(o["k_index"]), K = haskey(o, "K") ? Float64(o["K"]) : nothing, class = String(o["class"]), control = o["control"] === true,
                    flags = flags(id), rho = p2_rho(lab.instability, lab.structure), lambda_ref = _lam_at(ctx, id, p.N_1), lambda_ns = _lam_at(ctx, id, p.n_s),
                    in_SD = o["in_SD"] === true, in_SL = o["in_SL"] === true, instability = lab.instability, structure = lab.structure))
    end
    rows = NamedTuple[]
    for g in ctx.grid
        lab = _grid_label(ctx, g)
        push!(out, (id = g.id, source = :grid, kindex = g.kindex, K = parse(Float64, g.K), class = "h4", control = false, flags = flags(g.id), rho = lab.rho,
                    lambda_ref = _lam_at(ctx, g.id, p.N_1), lambda_ns = _lam_at(ctx, g.id, p.n_s), in_SD = false, in_SL = false,
                    instability = lab.instability, structure = lab.structure))
        push!(rows, (id = g.id, label = g.label, kindex = g.kindex, i = g.i, j = g.j, K = parse(Float64, g.K), x0 = parse(Float64, g.x0), p0 = parse(Float64, g.p0),
                     instability = lab.instability, structure = lab.structure, rho = lab.rho, ref_weak = true))
    end
    return out, rows
end
