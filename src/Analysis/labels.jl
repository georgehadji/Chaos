# Labels of D36 and D37 (E14) and H5 (E15).

const NONE_DETECTED = "none detected on the tested horizon"
const UNRESOLVED = "unresolved"

"r(N) = N lambda(N) / [(N/10) lambda(N/10)] = 10 lambda(N) / lambda(N/10); `nothing` when lambda(N/10) is not positive."
growth_ratio(lam_N::Real, lam_N10::Real) = lam_N10 > 0 ? 10 * lam_N / lam_N10 : nothing

function _instability(d, p)
    r = growth_ratio(d.lam_N, d.lam_N10)
    r === nothing && return UNRESOLVED
    d.lam_N >= p.lambda_threshold && r >= p.growth_ratio && return POSITIVE
    d.lam_N < p.lambda_threshold && r < p.growth_ratio && return NONE_DETECTED
    return UNRESOLVED
end

"""
    base_label(d, p) -> (instability, structure)

The two axes before any disagreement test. `d` carries `lam_N`, `lam_N10` (L5), `rho1`, `rho2` (rotation numbers of the two halves),
`periodic` (a periodic point of the panel) and `residual` (its return residual).
"""
function base_label(d, p)
    inst = _instability(d, p)
    str = if d.periodic && d.residual !== nothing && d.residual <= p.return_tol
        "periodic"
    elseif inst == NONE_DETECTED && d.rho1 !== nothing && d.rho2 !== nothing && abs(d.rho1 - d.rho2) <= p.rotation_tol
        "quasiperiodic"
    elseif inst == POSITIVE
        "irregular"
    else
        UNRESOLVED
    end
    return (instability = inst, structure = str)
end

function _ladder_forces(l, p)
    return abs(l.lam_f64 - l.lam_d64) > p.ladder_tol || abs(l.lam_d64 - l.lam_bf) > p.ladder_tol || !(l.labels.f64 == l.labels.d64 == l.labels.bf)
end

"""
    orbit_label(d, p) -> (instability, structure, forced, raw, ladder_same)

`base_label` plus the tests that force both axes to "unresolved": ladder (`d.ladder`, S_L only), estimators (`estimator_fail`), direction
(`d.direction`, the five labels, S_D only), diagnostics (`diag_exceeded`) and exclusion. `forced` lists the tests that fired, `ladder_same` is
whether the per-precision labels agree (before the ladder test acts; `nothing` outside S_L).
"""
function orbit_label(d, p)
    raw = base_label(d, p)
    ladder_same = d.ladder === nothing ? nothing : (d.ladder.labels.f64 == d.ladder.labels.d64 == d.ladder.labels.bf)
    forced = Symbol[]
    d.ladder !== nothing && _ladder_forces(d.ladder, p) && push!(forced, :ladder)
    d.estimator_fail && push!(forced, :estimator)
    d.direction !== nothing && !all(==(first(d.direction)), d.direction) && push!(forced, :direction)
    d.diag_exceeded && push!(forced, :diagnostics)
    d.excluded && push!(forced, :excluded)
    lab = isempty(forced) ? raw : (instability = UNRESOLVED, structure = UNRESOLVED)
    return (instability = lab.instability, structure = lab.structure, forced = forced, raw = raw, ladder_same = ladder_same)
end

function _h5_sub(k, n, limit, atleast, p)
    f = _fraction(k, n)
    n < p.min_set && return (fraction = f, n = n, pass = nothing, verdict = :not_decidable)
    ok = atleast ? f >= limit : f <= limit
    return (fraction = f, n = n, pass = ok, verdict = ok ? :pass : :fail)
end

"""
    decide_h5(orbits, p) -> (verdict, h5a, h5b, h5c)

Orbits carry `in_SL`, `ladder_same`, the final `instability` and `structure` and `excluded` (left out of every fraction). H5a: ladder
stability on S_L; H5b: share of orbits with an unresolved axis; H5c: share of unresolved structure among "none detected". Supported only if
all three are decidable and pass.
"""
function decide_h5(orbits, p)
    os = [o for o in orbits if !o.excluded]
    sl = [o for o in os if o.in_SL]
    none = [o for o in os if o.instability == NONE_DETECTED]
    a = _h5_sub(count(o -> o.ladder_same === true, sl), length(sl), p.h5a_min, true, p)
    b = _h5_sub(count(o -> o.instability == UNRESOLVED || o.structure == UNRESOLVED, os), length(os), p.h5b_max, false, p)
    c = _h5_sub(count(o -> o.structure == UNRESOLVED, none), length(none), p.h5c_max, false, p)
    subs = (a, b, c)
    verdict = any(s -> s.pass === false, subs) ? :not_supported : (any(s -> s.pass === nothing, subs) ? :not_decidable : :supported)
    return (verdict = verdict, h5a = a, h5b = b, h5c = c)
end
