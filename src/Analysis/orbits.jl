# Per-orbit quantities (E3 to E6, E10, E11, E14): everything the sets and the decisions need from the items of one panel orbit.

_missing(id, what) = throw(ArgumentError("item missing from the store: orbit $id, $what"))

# Item by property; a required item that is absent refuses the run.
function _need(idx, id, m, precision, horizon; tau = 1, eps = nothing, dir = "pi/4", mode = "lifted")
    it = get(idx, item_key(id, m, precision, horizon, tau, eps, dir, mode), nothing)
    it === nothing && _missing(id, "method $m, $precision, horizon $horizon, tau $tau, eps $eps, direction $dir, $mode")
    return it
end
_has(idx, id, m, precision, horizon; kw...) = get(idx, item_key(id, m, precision, horizon, get(kw, :tau, 1), get(kw, :eps, nothing), get(kw, :dir, "pi/4"), get(kw, :mode, "lifted")), nothing) !== nothing

# Primary item (E1): V methods lifted, L2 wrapped with the primary perturbation size.
_primary(idx, id, m, horizon, p) = m == "L2" ? _need(idx, id, "L2", "Float64", horizon; eps = p.h1c_eps, mode = "wrapped") : _need(idx, id, m, "Float64", horizon)

function _step(it, n::Integer)
    i = searchsortedfirst(it.t, n)
    (i <= length(it.t) && it.t[i] == n) || throw(ArgumentError("recorded step $n is missing in orbit $(it.orbit), method $(it.method), horizon $(it.horizon)"))
    return i
end
_val(it, n) = it.lam[_step(it, n)]

# Torus distance of a coordinate difference to the nearest multiple of 2 pi.
_delta(s) = abs(s - 2 * pi * round(s / (2 * pi)))

# Rotation numbers of the two halves (E14, P5) from the lifted x; the half point is the recorded step nearest N/2 (N/2 itself where recorded).
function _rotation(it, N::Integer, x0, h::Union{Nothing, Integer} = nothing)
    x = view(it.state, :, 1)
    h === nothing && (h = it.t[argmin(abs.(it.t .- N / 2))])
    return ((x[_step(it, h)] - x0) / (2 * pi * h), (x[_step(it, N)] - x[_step(it, h)]) / (2 * pi * (N - h)))
end

_is_periodic(o) = o["control"] && o["control_kind"] in ("fixed", "newton")

# Return residual of an exact periodic point: torus distance of the stored point at the known period to the start.
function _return_residual(it, o)
    i = _step(it, Int(o["period"]))
    return sqrt(_delta(it.state[i, 1] - o["x0"])^2 + _delta(it.state[i, 2] - o["p0"])^2)
end

# The L5 part of a label input at horizon N from one item.
function _label_core(it, N::Integer, o, h::Union{Nothing, Integer} = nothing)
    rho1, rho2 = _rotation(it, N, o["x0"], h)
    per = _is_periodic(o)
    return (lam_N = _val(it, N), lam_N10 = _val(it, N ÷ 10), rho1 = rho1, rho2 = rho2, periodic = per, residual = per ? _return_residual(it, o) : nothing)
end

_pair(lab) = (lab.instability, lab.structure)

# Ladder data of an S_L orbit at N_B (E14 a): the labels and exponents of L5 at Float64, Double64 and BigFloat.
function _ladder(idx, o, p, bf::AbstractString)
    id, NB = Int(o["id"]), p.N_B
    its = (_need(idx, id, "L5", "Float64", p.N_1), _need(idx, id, "L5", "Double64", p.N_1), _need(idx, id, "L5", bf, NB))
    common = intersect(its[1].t, its[2].t, its[3].t)          # one half point for the three precisions: a step all three items recorded
    h = common[argmin(abs.(common .- NB / 2))]
    cores = [_label_core(it, NB, o, h) for it in its]
    labs = [_pair(base_label(c, p)) for c in cores]
    return (lam_f64 = cores[1].lam_N, lam_d64 = cores[2].lam_N, lam_bf = cores[3].lam_N, labels = (f64 = labs[1], d64 = labs[2], bf = labs[3]))
end

# Labels of an S_D orbit at N_1 for each initial direction (E14 c); `nothing` when the store has fewer than two directions.
function _direction_labels(idx, o, p, dirs)
    ds = get(dirs, Int(o["id"]), String[])
    length(ds) < 2 && throw(ArgumentError("orbit $(o["id"]) is in S_D but the store has fewer than two initial directions for it"))
    return [_pair(base_label(_label_core(_need(idx, Int(o["id"]), "L5", "Float64", p.N_1; dir = d), p.N_1, o), p)) for d in ds]
end

"""
    orbit_record(idx, o, p, flags, dirs, bf) -> NamedTuple

All per-orbit quantities of panel orbit `o` (a panel-file table): the H1 metrics, lambda_ref, R(n), the convergence time, the label with its
tests and the flags. `flags` maps orbit ids to the union of the flags of their items, `dirs` to the directions of their L5 items.
"""
function orbit_record(idx, o, p, flags, dirs, bf)
    id, sd, sl = Int(o["id"]), Bool(o["in_SD"]), Bool(o["in_SL"])
    lam1 = Dict{String, Float64}(m => _val(_primary(idx, id, m, p.N_1, p), p.N_1) for m in vcat(collect(V_METHODS), "L2"))
    h1 = h1_orbit(lam1, p)
    N = sd ? p.N_long : p.N_1          # an S_D orbit without its long item refuses the run (E3, I4)
    l5 = _primary(idx, id, "L5", N, p)
    main = [_primary(idx, id, m, p.N_1, p) for m in V_METHODS]
    series = Dict{String, Vector{Float64}}(m => it.lam for (m, it) in zip(V_METHODS, main))
    lam_ref = _val(l5, N)
    steps = main[1].t
    R = relative_spread(series, lam_ref)
    conv = orbit_convergence(steps, series, lam_ref, p.h3_band)
    fl = get(flags, id, Set{String}())
    l5main = _primary(idx, id, "L5", p.N_1, p)
    d = (; _label_core(l5, N, o)..., estimator_fail = !(h1.pass_a && h1.pass_b), diag_exceeded = "diagnostics exceeded" in l5main.flags || "diagnostics exceeded" in l5.flags,
         excluded = "excluded" in fl, ladder = sl ? _ladder(idx, o, p, bf === nothing ? throw(ArgumentError("orbit $id is in S_L but the store has no BigFloat item")) : bf) : nothing, direction = sd ? _direction_labels(idx, o, p, dirs) : nothing)
    label = orbit_label(d, p)
    applicable = Symbol[s for (s, on) in ((:ladder, d.ladder !== nothing), (:direction, d.direction !== nothing), (:estimator, true), (:diagnostics, true)) if on]
    return (id = id, horizon = N, lam1 = lam1, h1 = h1, lam_ref = lam_ref, steps = steps, R = R, R_nH2 = R[_step(l5main, p.n_H2)], conv = conv,
            flags = sort!(collect(fl)), label = label, applicable = applicable, series = series, in_SD = sd, in_SL = sl)
end
