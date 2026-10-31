# Study p3-abstention (Amendment D93): the relabelling C0 on the stores (E49): the inputs of the reduced rule of every orbit of the panel (P1 items) and of the grid (p2 items), the
# data that are missing (counted, not refused), lambda_ref and lambda_ns (E49 d) and the change counts (E49 f).

const P3_DIAG = "diagnostics exceeded"

# The row of recorded step `n` in an item, or `nothing`.
_p3_step(it, n) = (i = searchsortedfirst(it.t, n); (i <= length(it.t) && it.t[i] == n) ? i : nothing)

# The L5 and L1 inputs of C0 from the two items and the lifted first coordinate `xl(it, row)`; `missing` lists what is lacking (E49 e), `data` is `nothing` then.
function _p3_c0_inputs(it5, it1, N, xl)
    miss = String[]
    r5 = Union{Int, Nothing}[nothing, nothing, nothing]
    if it5 === nothing
        push!(miss, "L5")
    else
        r5 = Union{Int, Nothing}[_p3_step(it5, n) for n in (N ÷ 10, N ÷ 2, N)]
        for (r, name) in zip(r5, ("L5 step N_1/10", "L5 step N_1/2", "L5 step N_1"))
            r === nothing && push!(miss, name)
        end
    end
    r1 = it1 === nothing ? nothing : _p3_step(it1, N)
    it1 === nothing ? push!(miss, "L1") : r1 === nothing && push!(miss, "L1 step N_1")
    isempty(miss) || return (data = nothing, missing = miss)
    x0 = parse(Float64, first(split(String(it5.attrs["ic"]), ',')))
    h = N ÷ 2
    rho1, rho2 = p2_rotation_numbers(x0, xl(it5, r5[2]), xl(it5, r5[3]), h, N)
    d = (lam5 = it5.lam[r5[3]], lam5_tenth = it5.lam[r5[1]], lam1 = it1.lam[r1], rho1 = rho1, rho2 = rho2, diag5 = P3_DIAG in it5.flags, diag1 = P3_DIAG in it1.flags)
    return (data = d, missing = String[])
end

_p3_lifted_p1(it, r) = it.state[r, 1]
_p3_lifted_p2(it, r) = it.state[r, 1] + 2 * pi * it.winding[r, 1]

# One relabelled orbit: its row, or the list of what is missing.
function _p3_relabel_orbit(id, source, kindex, it5, it1, N, xl, p2)
    r = _p3_c0_inputs(it5, it1, N, xl)
    r.data === nothing && return (row = nothing, missing = [(id = id, source = source, kindex = kindex, item = m) for m in r.missing], d = nothing)
    lab = p3_relabel(r.data, p2)
    row = (id = id, source = source, kindex = kindex, instability = lab.instability, structure = lab.structure, rho = lab.rho, rho1 = r.data.rho1, rho2 = r.data.rho2,
           lam5 = r.data.lam5, lam1 = r.data.lam1)
    return (row = row, missing = NamedTuple[], d = r.data)
end

"""
    p3_relabel_orbits(panel_orbits, grid, idx1, idx2, p2, old) -> (rows, missing, changes)

C0 on the stores (E49): `rows` (`id, source, kindex, instability, structure, rho, rho1, rho2, lam5, lam1`) of the orbits that have their data, `missing` rows
`(id, source, kindex, item)`, and the change counts of E49 f: `changes.panel` and `changes.grid` map `(K index, old, new)` to a count (panel old = the base label of the P1
inputs, grid old = the D89 rule without the diagnostics condition), with `panel_changed` and `grid_changed`.
"""
function p3_relabel_orbits(panel_orbits, grid, idx1, idx2, p2, old)
    N, d1 = p2.N_1, p2.directions[1]
    rows, miss = NamedTuple[], NamedTuple[]
    panel, gridc = Dict{Tuple{Int, Symbol, Symbol}, Int}(), Dict{Tuple{Int, Symbol, Symbol}, Int}()
    add(d, k) = d[k] = get(d, k, 0) + 1
    for o in panel_orbits
        o["control"] === true && continue
        id, k = Int(o["id"]), Int(o["k_index"])
        r = _p3_relabel_orbit(id, :panel, k, get(idx1, item_key(id, "L5", "Float64", N, 1, nothing, d1, "lifted"), nothing),
                              get(idx1, item_key(id, "L1", "Float64", N, 1, nothing, d1, "lifted"), nothing), N, _p3_lifted_p1, p2)
        r.row === nothing ? append!(miss, r.missing) : (push!(rows, r.row); add(panel, (k, get(old, id, :missing), r.row.rho)))
    end
    for g in grid
        r = _p3_relabel_orbit(g.id, :grid, g.kindex, get(idx2, p2_key(g.id, "L5", "Float64", N, 1, d1), nothing), get(idx2, p2_key(g.id, "L1", "Float64", N, 1, d1), nothing), N, _p3_lifted_p2, p2)
        r.row === nothing ? append!(miss, r.missing) : (push!(rows, r.row); add(gridc, (g.kindex, p2_grid_label(r.d, p2).rho, r.row.rho)))
    end
    return (rows = rows, missing = miss, changes = (panel = panel, grid = gridc, panel_changed = sum(c for ((_, o, n), c) in panel if o !== n; init = 0),
                                                    grid_changed = sum(c for ((_, o, n), c) in gridc if o !== n; init = 0)))
end

"""
    p3_orbit_rows(panel_orbits, grid, ctx, rel, flags) -> Vector{NamedTuple}

The orbit rows of `p3_sets` for every orbit of P and G (E50): reference class `rho` from the relabelling (`:missing` where its data are missing), the stored flags, and `lambda_ref`
and `lambda_ns` from the p2 wrapped L1 item at N_1 and at n_s (E49 d). A control has neither (it is dropped first).
"""
function p3_orbit_rows(panel_orbits, grid, ctx, rel, flags)
    p = ctx.p
    rho = Dict(r.id => r.rho for r in rel.rows)
    fl(id) = sort!(collect(get(flags, id, Set{String}())))
    out = NamedTuple[]
    has(id) = haskey(ctx.idx, p2_key(id, "L1", "Float64", p.N_1, 1, p.directions[1]))
    lam(id) = has(id) ? (_lam_at(ctx, id, p.N_1), _lam_at(ctx, id, p.n_s)) : (NaN, NaN)
    for o in panel_orbits
        id, ctl = Int(o["id"]), o["control"] === true
        a, b = ctl ? (NaN, NaN) : lam(id)
        push!(out, (id = id, source = :panel, kindex = Int(o["k_index"]), class = String(o["class"]), control = ctl, flags = fl(id), rho = get(rho, id, :missing), lambda_ref = a, lambda_ns = b))
    end
    for g in grid
        a, b = lam(g.id)
        push!(out, (id = g.id, source = :grid, kindex = g.kindex, class = "h4", control = false, flags = fl(g.id), rho = get(rho, g.id, :missing), lambda_ref = a, lambda_ns = b))
    end
    return out
end
