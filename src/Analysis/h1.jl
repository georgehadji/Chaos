# H1 (E4 to E6, E8) and the exact-exponent controls (E9).

"The variational set V of the preregistration and its subset B (V without the SVD method L9)."
const V_METHODS = ("L1", "L3", "L4", "L4b", "L5", "L5b", "L6", "L9")
const B_METHODS = ("L1", "L3", "L4", "L4b", "L5", "L5b", "L6")

"""
    h1_orbit(lam, p) -> NamedTuple

`lam` maps method ids (V and L2) to lambda_1 at N_1. Delta_QB over B, |L9 - L5| and |L2 - L5| with their pass flags (bounds inclusive).
"""
function h1_orbit(lam::AbstractDict, p)
    qb = [lam[m] for m in B_METHODS]
    delta_qb = maximum(qb) - minimum(qb)
    diff_l9, diff_l2 = abs(lam["L9"] - lam["L5"]), abs(lam["L2"] - lam["L5"])
    return (delta_qb = delta_qb, diff_l9 = diff_l9, diff_l2 = diff_l2,
            pass_a = delta_qb <= p.h1a_bound, pass_b = diff_l9 <= p.h1b_bound, pass_c = diff_l2 <= p.h1c_bound)
end

_fraction(k, n) = n == 0 ? NaN : k / n

"""
    decide_h1(orbits, p) -> (verdict, n, fractions, passes)

`orbits` carry `pass_a`, `pass_b`, `pass_c`. Supported if all three fractions reach their limits, falsified if one does not, not decidable
(with the count) when fewer than `p.min_set` orbits are eligible.
"""
function decide_h1(orbits, p)
    n = length(orbits)
    passes = (a = count(o -> o.pass_a, orbits), b = count(o -> o.pass_b, orbits), c = count(o -> o.pass_c, orbits))
    fr = (a = _fraction(passes.a, n), b = _fraction(passes.b, n), c = _fraction(passes.c, n))
    verdict = n < p.min_set ? :not_decidable : (fr.a >= p.h1a_fraction && fr.b >= p.h1b_fraction && fr.c >= p.h1c_fraction ? :supported : :falsified)
    return (verdict = verdict, n = n, fractions = fr, passes = passes)
end

"Absolute deviations |lambda_m - exact| for every method of `lam`."
control_deviations(lam::AbstractDict, exact::Real) = Dict{String, Float64}(String(m) => abs(v - exact) for (m, v) in lam)

_control_bound(m, p) = m == "L9" ? p.control_bound_L9 : p.control_bound

"""
    check_controls(rows, p) -> (ok, defects, info, mismatches)

Rows carry `id`, `exponent`, `scope` ("bound" or "info"), `lam` (method => lambda_1 at N_1) and `stored` (the `control_deviation` of each item or
`nothing`). A "bound" control with positive exponent whose deviation exceeds the bound of its method is a defect; "info" controls are listed;
a stored deviation that differs from the recomputed one by more than 1e-12 is a mismatch.
"""
function check_controls(rows, p)
    defects, info, mismatches = NamedTuple[], NamedTuple[], Int[]
    for r in rows
        dev = control_deviations(Dict(m => r.lam[m] for m in V_METHODS if haskey(r.lam, m)), r.exponent)
        if r.scope == "bound" && r.exponent > 0
            for m in V_METHODS
                haskey(dev, m) && dev[m] > _control_bound(m, p) && push!(defects, (id = r.id, method = m, deviation = dev[m], bound = _control_bound(m, p)))
            end
        elseif r.scope == "info"
            push!(info, (id = r.id, deviations = dev))
        end
        r.stored === nothing || any(m -> haskey(r.stored, m) && abs(r.stored[m] - dev[m]) > 1e-12, keys(dev)) && push!(mismatches, r.id)
    end
    return (ok = isempty(defects), defects = defects, info = info, mismatches = mismatches)
end
