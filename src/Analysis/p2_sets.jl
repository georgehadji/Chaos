# Study p2-indicators (Amendment D89): the sets of E29 and the reference class rho of E23 (rule of the panel labels, reduced rule of the grid).

"True iff a set of `n` orbits is large enough for a statistic (E29): `n >= min_set`."
p2_decidable(n::Integer, p) = n >= p.min_set

# Reason an orbit leaves E_all (E29, I16): control, then drift or excluded, then unresolved; "" keeps it.
function _drop_reason(o)
    o.control && return "control"
    "drift" in o.flags && return "drift"
    "excluded" in o.flags && return "excluded"
    return o.rho === :U ? "unresolved" : ""
end

"""
    p2_sets(orbits, p) -> NamedTuple

The sets of E29 from orbit records `(id, source, kindex, K, class, control, flags, rho, lambda_ref, lambda_ns, in_SD, in_SL)`: `E_all`, `E8`, `W2`, `S2`, `T_s`,
`T_s_design`, `regular`, `gate`, `S_D`, `S_L` (S_D* and S_L*), `overlap` of S2 and W2 and `dropped` (`(id, reason)`, each orbit once), all sorted by id.
An `ArgumentError` naming the gate set if its size differs from `gate_size`.
"""
function p2_sets(orbits, p)
    dropped = NamedTuple[(id = o.id, reason = r) for o in orbits for r in (_drop_reason(o),) if !isempty(r)]
    keep = [o for o in orbits if isempty(_drop_reason(o))]
    ids(f, os = keep) = sort!(Int[o.id for o in os if f(o)])
    chaotic = [o for o in keep if o.rho === :C]
    Ts = ids(o -> o.lambda_ns <= p.f_s * o.lambda_ref, chaotic)
    S2 = ids(o -> o.lambda_ref >= p.s2_lambda_min, chaotic)
    W2 = sort!(unique(vcat(ids(o -> o.lambda_ref <= p.w2_lambda_max, chaotic), Ts)))
    gate = ids(o -> o.source === :panel && !o.control && o.kindex == p.gate_k_index && o.class in p.gate_classes, orbits)
    length(gate) == p.gate_size || throw(ArgumentError("the gate set has $(length(gate)) orbits, the parameter table says gate_size = $(p.gate_size)"))
    return (E_all = ids(_ -> true), E8 = ids(o -> o.in_SL || o.source === :grid), W2 = W2, S2 = S2, T_s = Ts,
            T_s_design = ids(o -> o.source === :panel && o.class == "S" && o.rho === :C), regular = ids(o -> o.rho === :R), gate = gate,
            S_D = ids(o -> o.in_SD), S_L = ids(o -> o.in_SL), overlap = sort!(intersect(S2, W2)), dropped = sort!(dropped; by = d -> d.id))
end

"The reference class of E23 from the two label axes: C iff (positive growth supported, irregular), R iff (none detected on the tested horizon, quasiperiodic), else U."
function p2_rho(instability::AbstractString, structure::AbstractString)
    instability == POSITIVE && structure == "irregular" && return :C
    instability == NONE_DETECTED && structure == "quasiperiodic" && return :R
    return :U
end

"""
    p2_grid_label((lam5, lam5_tenth, lam1, rho1, rho2), p) -> (instability, structure, rho)

The reduced label rule of the grid orbits (E23 b): the instability and structure axes of the base labels with the parameters of the P2 table, and the
disagreement test of L1 against L5 (unresolved on both axes above `label_estimator_bound`). A grid orbit is never periodic.
"""
function p2_grid_label(d, p)
    q = (lambda_threshold = p.label_lambda_min, growth_ratio = p.label_growth_ratio, rotation_tol = p.label_rotation_tol, return_tol = 0.0)
    b = base_label((lam_N = d.lam5, lam_N10 = d.lam5_tenth, rho1 = d.rho1, rho2 = d.rho2, periodic = false, residual = nothing), q)
    bad = abs(d.lam1 - d.lam5) > p.label_estimator_bound
    inst, str = bad ? (UNRESOLVED, UNRESOLVED) : (b.instability, b.structure)
    return (instability = inst, structure = str, rho = p2_rho(inst, str))
end

"Rotation numbers of the two halves (E23 b): `((x_h - x_0) / (2 pi h), (x_N - x_h) / (2 pi (N - h)))` from the lifted first coordinate."
p2_rotation_numbers(x0, xh, xN, h, N) = ((xh - x0) / (2 * pi * h), (xN - xh) / (2 * pi * (N - h)))

"The lifted coordinate of a wrapped value and its winding count: `x + 2 pi winding` in binary64 (P19)."
p2_lift(x, winding) = x .+ 2 * pi .* winding
