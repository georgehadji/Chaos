# Study p3-abstention (Amendment D93): the verdicts of I_7 at the ten budgets (E51, table R1 of D89 unchanged), the orbit rows of the calibration and the rows of R_ab (P63).

"""
    p3_verdicts(ctx, ids) -> Dict{(id, X), Vector{Symbol}}

The verdicts of the seven indicators of I_7 at the ten budgets for the orbits `ids` (E51): table R1 of D89 on the value records of the p2 items at the recorded steps. An item that
the manifest lists as not applicable gives undecided verdicts. `ctx = (p, idx)` with the p2 table and the item index.
"""
function p3_verdicts(ctx, ids)
    p = ctx.p
    out = Dict{Tuple{Int, Symbol}, Vector{Symbol}}()
    for id in ids, X in P2_I7
        vals = p2_values(ctx, id, X)
        out[(id, X)] = vals === nothing ? fill(:undecided, length(p.budgets)) : Symbol[p2_verdict(X, :R1, v, p) for v in vals]
    end
    return out
end

"""
    p3_rows(orbits, S, verdicts, p2) -> Vector{NamedTuple}

One row per non-control orbit of P and G (`orbits` as for `p3_sets`, `S` its result): `id, source, kindex, in_cal, rho, S2, Ts, W2, gate, k0, V` where `V` holds the ten verdicts
of the seven indicators for the orbits of Cal and of the gate set and is `nothing` for the others; `k0` marks the K = 0 orbits (the panel orbits of the gate K index).
"""
function p3_rows(orbits, S, verdicts, p2)
    cal, s2, ts, w2, gate = Set(S.cal), Set(S.S2), Set(S.T_s), Set(S.W2), Set(S.gate)
    return [(id = o.id, source = o.source, kindex = o.kindex, in_cal = o.id in cal, rho = o.rho, S2 = o.id in s2, Ts = o.id in ts, W2 = o.id in w2, gate = o.id in gate,
             k0 = o.source === :panel && o.kindex == p2.gate_k_index,
             V = (o.id in cal || o.id in gate) ? Dict{Symbol, Vector{Symbol}}(X => verdicts[(o.id, X)] for X in P2_I7) : nothing) for o in orbits if !o.control]
end

"""
    p3_rab(ctx, rows, theta, kP) -> Vector{NamedTuple}

R_ab (E52, P63) on the orbits of Cal that carry the Double64 items of every indicator of `theta.V`: `(id, output)` at budget index `kP`. The Double64 verdicts are table R1 on the
Double64 items (horizon N_dir for the tangent-map indicators).
"""
function p3_rab(ctx, rows, theta, kP)
    p = ctx.p
    horizon(X) = X in (:SALI, :FLI, :MEGNO) ? (horizon = p.N_dir,) : (;)
    out = NamedTuple[]
    for o in rows
        o.in_cal || continue
        vals = [p2_values(ctx, o.id, X; prec = "Double64", horizon(X)...) for X in theta.V]
        any(isnothing, vals) && continue
        Vd = Dict{Symbol, Vector{Symbol}}(X => Symbol[p2_verdict(X, :R1, v, p) for v in vs] for (X, vs) in zip(theta.V, vals))
        push!(out, (id = o.id, output = p3_rule_ab(o.V, Vd, kP, theta)))
    end
    return out
end
