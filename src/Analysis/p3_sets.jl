# Study p3-abstention (Amendment D93): the relabelling rule C0 (E49 c) and the calibration orbits and sets (E50). Pure functions of orbit records.

"""
    p3_relabel(d, p2) -> (instability, structure, rho)

The reduced rule of the test (E49 c): the rule of the grid orbits of D89 (E23 b, `p2_grid_label`, the four label keys of the p2 table) plus the condition of the
preregistration that both axes are unresolved if a diagnostic fails: `d.diag5` (the L5 item read carries "diagnostics exceeded") or `d.diag1` (the L1 item).
`d = (lam5, lam5_tenth, lam1, rho1, rho2, diag5, diag1)`; `rho` is `:C`, `:R` or `:U`.
"""
function p3_relabel(d, p2)
    (d.diag5 || d.diag1) && return (instability = UNRESOLVED, structure = UNRESOLVED, rho = :U)
    return p2_grid_label(d, p2)
end

# Reason an orbit leaves Cal (E50): control, drift, excluded, unresolved, relabelling data missing; "" keeps it.
function _p3_drop_reason(o)
    o.control && return "control"
    "drift" in o.flags && return "drift"
    "excluded" in o.flags && return "excluded"
    o.rho === :U && return "unresolved"
    o.rho === :missing && return "missing"
    return ""
end

"""
    p3_sets(orbits, q, p2) -> NamedTuple

Cal and its sets (E50) from orbit rows `(id, source, kindex, class, control, flags, rho, lambda_ref, lambda_ns)` (`rho` in `:C :R :U :missing`): `cal`, `S2`, `T_s`,
`W2`, `R_W`, `regular`, `gate`, `overlap` (S2 and W2) as sorted id vectors and `dropped` rows `(id, reason)`. The sets of Cal with `rho = :C` use `lam_s2`, `f_s`,
`lam_w2` of the p3 table; the gate set (the panel orbits of K index `gate_k_index` in the classes `gate_classes` of the p2 table) is not intersected with Cal. An
`ArgumentError` naming the gate set if its size is not `gate_size`.
"""
function p3_sets(orbits, q, p2)
    dropped = NamedTuple[(id = o.id, reason = r) for o in orbits for r in (_p3_drop_reason(o),) if !isempty(r)]
    keep = [o for o in orbits if isempty(_p3_drop_reason(o))]
    ids(f, os = keep) = sort!(Int[o.id for o in os if f(o)])
    chaotic = [o for o in keep if o.rho === :C]
    Ts = ids(o -> o.lambda_ns <= q.f_s * o.lambda_ref, chaotic)
    W2 = sort!(unique(vcat(ids(o -> o.lambda_ref <= q.lam_w2, chaotic), Ts)))
    S2 = ids(o -> o.lambda_ref >= q.lam_s2, chaotic)
    gate = ids(o -> o.source === :panel && !o.control && o.kindex == p2.gate_k_index && o.class in p2.gate_classes, orbits)
    length(gate) == q.gate_size || throw(ArgumentError("the gate set has $(length(gate)) orbits, the parameter table says gate_size = $(q.gate_size)"))
    return (cal = ids(_ -> true), S2 = S2, T_s = Ts, W2 = W2, R_W = sort!(setdiff(W2, Ts)), regular = ids(o -> o.rho === :R), gate = gate,
            overlap = sort!(intersect(S2, W2)), dropped = sort!(dropped; by = d -> d.id))
end
