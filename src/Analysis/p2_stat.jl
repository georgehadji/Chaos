# Study p2-indicators (Amendment D89): the statistics of H1 to H6 and their decision rules (E32, E33). A record is (id, V, T): V maps an indicator to its ten verdicts,
# T maps it to the censored detection time; an indicator without an item is absent from both.

_kp(p) = something(findfirst(==(p.b_primary), p.budgets), 0)

_med_or_nothing(v) = isempty(v) ? nothing : p2_median(v)

# H1: U and the agreement statistics at the primary budget.
function _stat_h1(recs, p, I, kp)
    V = Dict(X => Symbol[haskey(r.V, X) ? r.V[X][kp] : :na for r in recs] for X in I)
    u = p2_unanimity(V, I)
    return (n = length(recs), used = u.used, dropped = u.dropped, U = u.U, A_min = u.A_min, A = u.A, kappa_mean = u.kappa_mean)
end

# H2: M_X over the set and Q = max M_X / min M_X.
function _medians(recs, I)
    return Dict{Symbol, Union{Nothing, Float64}}(X => _med_or_nothing([r.T[X] for r in recs if haskey(r.T, X)]) for X in I)
end

_ratio_or_nothing(a, b) = (a === nothing || b === nothing) ? nothing : a / b

function _spread(M)
    v = Float64[x for x in values(M) if x !== nothing]
    return (isempty(v) || length(v) != length(M)) ? nothing : maximum(v) / minimum(v)
end

# max - min of the defined values; `nothing` unless every entry is defined.
function _range_of(phi)
    v = Float64[x for x in values(phi) if x !== nothing]
    return (isempty(v) || length(v) != length(phi)) ? nothing : maximum(v) - minimum(v)
end

# H4: phi_X, the share regular among the orbits decided at the last budget, for the indicators without LE.
function _stat_h4(recs, p, I)
    nb = length(p.budgets)
    J = [X for X in I if X !== :LE]
    ok(r, X) = haskey(r.V, X) && _decided(r.V[X][nb])
    dec = Dict{Symbol, Int}(X => count(r -> ok(r, X), recs) for X in J)
    phi = Dict{Symbol, Union{Nothing, Float64}}(X => dec[X] == 0 ? nothing : count(r -> ok(r, X) && r.V[X][nb] === :regular, recs) / dec[X] for X in J)
    return (n = length(recs), phi = phi, decided = dec, dphi = _range_of(phi))
end

# H5: F_X, the share of verdicts that differ between b_P and the last budget among the orbits decided at both; F_med, the median of the F_X.
function _stat_h5(recs, p, I, kp)
    nb = length(p.budgets)
    dec = Dict{Symbol, Int}()
    F = Dict{Symbol, Union{Nothing, Float64}}()
    for X in I
        sel = [r for r in recs if haskey(r.V, X) && _decided(r.V[X][kp]) && _decided(r.V[X][nb])]
        dec[X] = length(sel)
        F[X] = isempty(sel) ? nothing : count(r -> r.V[X][kp] !== r.V[X][nb], sel) / length(sel)
    end
    v = [F[X] for X in I]
    return (n = length(recs), F = F, decided = dec, F_med = (isempty(v) || any(isnothing, v)) ? nothing : p2_median(Float64.(v)))
end

"""
    p2_stat(h, sets, p; I = I_7, kp) -> NamedTuple

Statistics of hypothesis `h` (`:H1` to `:H5`) on the tuple of record sets (E32). `nothing` marks an undefined quantity.
H1 `(n, used, dropped, U, A_min, A, kappa_mean)`, H2 `(n, M, Q)`, H3 (sets W, S) `(n_w, n_s, M_w, M_s, R)`, H4 `(n, phi, decided, dphi)`, H5 `(n, F, decided, F_med)`.
"""
function p2_stat(h::Symbol, sets::Tuple, p; I = P2_I7, kp = _kp(p))
    h === :H1 && return _stat_h1(sets[1], p, I, kp)
    h === :H4 && return _stat_h4(sets[1], p, I)
    h === :H5 && return _stat_h5(sets[1], p, I, kp)
    if h === :H2
        M = _medians(sets[1], I)
        return (n = length(sets[1]), M = M, Q = _spread(M))
    elseif h === :H3
        w, s = sets
        Mw, Ms = _medians(w, I), _medians(s, I)
        R = Dict{Symbol, Union{Nothing, Float64}}(X => _ratio_or_nothing(Mw[X], Ms[X]) for X in I)
        return (n_w = length(w), n_s = length(s), M_w = Mw, M_s = Ms, R = R)
    end
    throw(ArgumentError("unknown hypothesis $h"))
end

# The fraction of the decided rows whose verdicts differ, and the number of rows left out (P23).
function _delta(rows, differs, decided)
    dec = [r for r in rows if decided(r)]
    return (isempty(dec) ? nothing : count(differs, dec) / length(dec), length(rows) - length(dec))
end

_ordered(ks) = sort!(collect(ks); by = X -> something(findfirst(==(X), P2_I8), 99))

"""
    p2_h6_stats(dir, prec, p) -> (dir, prec)

H6 statistics (E32): `delta_dir` from rows `(id, v)` with the five verdicts at the last budget per indicator, `delta_prec` from rows `(id, f64, d64)`.
Each part is `(n, delta, dropped)`, dictionaries by indicator; an orbit undecided in some direction or precision is dropped and counted.
"""
function p2_h6_stats(dir, prec, p)
    dd, dr = Dict{Symbol, Any}(), Dict{Symbol, Int}()
    for X in _ordered(union(Set{Symbol}(), (keys(r.v) for r in dir)...))
        dd[X], dr[X] = _delta(dir, r -> !all(==(r.v[X][1]), r.v[X]), r -> all(_decided, r.v[X]))
    end
    pd, pr = Dict{Symbol, Any}(), Dict{Symbol, Int}()
    for X in _ordered(union(Set{Symbol}(), (keys(r.f64) for r in prec)...))
        pd[X], pr[X] = _delta(prec, r -> r.f64[X] !== r.d64[X], r -> _decided(r.f64[X]) && _decided(r.d64[X]))
    end
    return (dir = (n = length(dir), delta = dd, dropped = dr), prec = (n = length(prec), delta = pd, dropped = pr))
end

"The flip count Phi_X: changes between consecutive budgets where both verdicts are decided (E32, H5)."
p2_flips(V::AbstractVector{Symbol}) = count(k -> _decided(V[k]) && _decided(V[k + 1]) && V[k] !== V[k + 1], 1:(length(V) - 1))

_sup(c) = c ? :supported : :not_supported
_vals(x) = x isa AbstractDict ? collect(values(x)) : collect(x)

function _decide_h3(s, p)
    (s.n_w < p.min_set || s.n_s < p.min_set || s.R_lo === nothing) && return :not_decidable
    ks = s.R isa AbstractDict ? collect(keys(s.R)) : collect(eachindex(s.R))
    (any(k -> s.R[k] === nothing || s.R_lo[k] === nothing, ks)) && return :not_decidable
    return _sup(count(k -> s.R[k] >= p.h3_r_min && s.R_lo[k] > p.h3_r_lower_gt, ks) >= p.h3_count_min)
end

function _decide_h6(s, p)
    da, db = s.dir.n >= p.min_set, s.prec.n >= p.min_set
    reach(part) = (v = filter(!isnothing, _vals(part.delta)); !isempty(v) && maximum(v) >= p.h6_flip_min)
    ((da && reach(s.dir)) || (db && reach(s.prec))) && return :supported
    return (da && db) ? :not_supported : :not_decidable
end

"""
    p2_decide(h, s, p) -> NamedTuple

The decision of hypothesis `h` (E32) from its statistics and the lower end of its interval (`U_lo`, `Q_lo`, `R_lo`, `dphi_lo`, `F_med_lo`; H6 `(dir, prec)`):
the inputs echoed plus `verdict` in `:supported`, `:not_supported`, `:not_decidable`. A set below `min_set` or an undefined quantity is not decidable.
"""
function p2_decide(h::Symbol, s::NamedTuple, p)
    v = if h === :H1
        (s.n < p.min_set || s.U === nothing || s.U_lo === nothing || s.A_min === nothing) ? :not_decidable :
            _sup(s.U >= p.h1_u_min && s.U_lo >= p.h1_u_lower_min && s.A_min <= p.h1_amin_max)
    elseif h === :H2
        (s.n < p.min_set || s.Q === nothing || s.Q_lo === nothing) ? :not_decidable : _sup(s.Q >= p.h2_q_min && s.Q_lo > p.h2_q_lower_gt)
    elseif h === :H3
        _decide_h3(s, p)
    elseif h === :H4
        (s.n < p.min_set || s.dphi === nothing || s.dphi_lo === nothing) ? :not_decidable : _sup(s.dphi >= p.h4_dphi_min && s.dphi_lo >= p.h4_dphi_lower_min)
    elseif h === :H5
        (s.n < p.min_set || s.F_med === nothing || s.F_med_lo === nothing) ? :not_decidable : _sup(s.F_med >= p.h5_fmed_min && s.F_med_lo >= p.h5_fmed_lower_min)
    elseif h === :H6
        _decide_h6(s, p)
    else
        throw(ArgumentError("unknown hypothesis $h"))
    end
    return merge(s, (verdict = v,))
end

"ROBUST iff the primary and the extension verdicts are both decidable and equal, NOT ROBUST iff both are decidable and differ, NOT COMPARABLE otherwise (E33)."
function p2_robust(primary::Symbol, extension::Symbol)
    (primary === :not_decidable || extension === :not_decidable) && return :NOT_COMPARABLE
    return primary === extension ? :ROBUST : :NOT_ROBUST
end
