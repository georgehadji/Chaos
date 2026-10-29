# Study p2-indicators (Amendment D89): the items of the hypotheses T3 to T8, the R2 version T10, the gate T12 and the extension T13 (E32, E33, E37).
# A row is [source, K index, statistic, value, lower, upper]; the pooled rows come first, then the rows by source and by K index with counts and point values.

_vstr(v) = replace(string(v), "_" => " ")
_row(name, v, lo = nothing, hi = nothing) = Any["all", "all", name, v, lo, hi]
_byx(d, name, I) = [_row("$name $X", d[X]) for X in I if haskey(d, X)]

function _flat_h1(d, I)
    A = [_row("A $X-$Y", d.A[(X, Y)]) for (X, Y) in sort!(collect(keys(d.A)); by = pr -> (findfirst(==(pr[1]), P2_I8), findfirst(==(pr[2]), P2_I8)))]
    return vcat([_row("verdict", _vstr(d.verdict)), _row("orbits", d.n), _row("orbits decided for every indicator", d.used), _row("orbits dropped", d.dropped),
                 _row("U", d.U, d.U_lo, d.U_hi), _row("A_min", d.A_min), _row("mean kappa", d.kappa_mean), _row("U on S2 (no rule)", d.U_S2)], A)
end

_flat_h2(d, I) = vcat([_row("verdict", _vstr(d.verdict)), _row("orbits", d.n)], _byx(d.M, "M", I), [_row("Q", d.Q, d.Q_lo, d.Q_hi)])

function _flat_h3(d, I)
    R = [_row("R $X", d.R[X], d.R_lo === nothing ? nothing : d.R_lo[X], d.R_hi === nothing ? nothing : d.R_hi[X]) for X in I]
    return vcat([_row("verdict", _vstr(d.verdict)), _row("orbits W2", d.n_w), _row("orbits S2", d.n_s)], _byx(d.M_w, "M_W2", I), _byx(d.M_s, "M_S2", I), R,
                [_row("indicators passing", d.passing)])
end

function _flat_h4(d, I; prefix = "")
    J = [X for X in I if X !== :LE]
    own = vcat([_row(prefix * "orbits", d.n)], _byx(d.phi, prefix * "phi", J), _byx(d.decided, prefix * "decided", J), [_row(prefix * "Delta_phi", d.dphi, d.dphi_lo, d.dphi_hi)])
    isempty(prefix) || return own
    return vcat([_row("verdict", _vstr(d.verdict))], own, [_row("design verdict", _vstr(d.design.verdict))], _flat_h4(d.design, I; prefix = "design "))
end

_flat_h5(d, I) = vcat([_row("verdict", _vstr(d.verdict)), _row("orbits", d.n)], _byx(d.F, "F", I), _byx(d.decided, "decided", I), [_row("F_med", d.F_med, d.F_med_lo, d.F_med_hi)])

function _flat_h6(d, I)
    rows = Any[_row("verdict", _vstr(d.verdict)), _row("orbits S_D*", d.dir.n)]
    append!(rows, _byx(d.dir.delta, "delta_dir", P2_I8), _byx(d.dir.dropped, "dropped dir", P2_I8), [_row("orbits S_L*", d.prec.n)])
    append!(rows, _byx(d.prec.delta, "delta_prec", P2_I8), _byx(d.prec.dropped, "dropped prec", P2_I8))
    return rows
end

const _FLAT = Dict(:H1 => _flat_h1, :H2 => _flat_h2, :H3 => _flat_h3, :H4 => _flat_h4, :H5 => _flat_h5, :H6 => _flat_h6)

# Point values by source and by K index (E37): the statistic of the hypothesis on the part of its sets that lies in the group.
function _break(h, orbits, S, vt, p, I, rule, suffix)
    orb = Dict(o.id => o for o in orbits)
    sel(v) = isempty(suffix) ? v : sort!(intersect(v, S.E8))
    rows = Vector{Any}[]
    for (src, kk, pred) in _groups(orbits)[2:end]
        part(v) = p2_records(vt, [id for id in sel(v) if pred(orb[id])], I, rule)
        put(name, v) = push!(rows, Any[src, kk, name, v, nothing, nothing])
        if h === :H1
            s = p2_stat(:H1, (part(S.W2),), p; I = I); put("orbits", s.n); put("U", s.U); put("A_min", s.A_min)
        elseif h === :H2
            s = p2_stat(:H2, (part(S.S2),), p; I = I); put("orbits", s.n); put("Q", s.Q)
        elseif h === :H3
            s = p2_stat(:H3, (part(S.W2), part(S.S2)), p; I = I); put("orbits W2", s.n_w); put("orbits S2", s.n_s); foreach(X -> put("R $X", s.R[X]), I)
        elseif h === :H4
            s = p2_stat(:H4, (part(S.T_s),), p; I = I); put("orbits", s.n); put("Delta_phi", s.dphi)
        elseif h === :H5
            s = p2_stat(:H5, (part(S.W2),), p; I = I); put("orbits", s.n); put("F_med", s.F_med)
        end
    end
    return rows
end

# H2: the shares of the three outcomes of E28 per indicator on S2. H5: F_X(b, A) at every budget for W2, S2, T_s and the regular set, and the distribution of the flip count Phi_X on W2 (M5, E32).
function _extra(h, S, vt, p, I, rule)
    rows = Vector{Any}[]
    if h === :H2
        for X in I
            d = [vt.detection[(id, X)] for id in S.S2 if haskey(vt.detection, (id, X))]
            for o in (:correct, :wrong, :undetected)
                push!(rows, _row("share $o $X", isempty(d) ? nothing : count(x -> x.outcome === o, d) / length(d)))
            end
        end
    elseif h === :H5
        nb = length(p.budgets)
        for (name, ids) in (("W2", S.W2), ("S2", S.S2), ("T_s", S.T_s), ("regular", S.regular)), k in 0:(nb - 1)
            F = p2_stat(:H5, (p2_records(vt, ids, I, rule),), p; I = I, kp = k + 1).F
            foreach(X -> push!(rows, _row("F $X $name b_$k", F[X])), I)
        end
        recs = p2_records(vt, S.W2, I, rule)
        for X in I, f in 0:(nb - 1)
            push!(rows, _row("Phi $X = $f on W2", count(r -> haskey(r.V, X) && p2_flips(r.V[X]) == f, recs)))
        end
    end
    return rows
end

"""
    p2_hypothesis_items(ctx, orbits, S, vt, dec, I; rule, suffix) -> NamedTuple

The items T3 to T8 (H1 to H6) of a set of decisions `dec`: the decision rows, then the rows by source and by K index.
"""
function p2_hypothesis_items(p, orbits, S, vt, dec, I; rule = :R1, suffix = "")
    cols = ["source", "K index", "statistic", "value", "lower", "upper"]
    tab(h) = _tbl(cols, vcat(_FLAT[h](dec[h], I), h === :H6 ? Vector{Any}[] : _break(h, orbits, S, vt, p, I, rule, suffix), _extra(h, S, vt, p, I, rule)))
    return (T3 = tab(:H1), T4 = tab(:H2), T5 = tab(:H3), T6 = tab(:H4), T7 = tab(:H5), T8 = tab(:H6))
end

# The statistics and decisions of H1 to H6 stacked in one table (T10 under R2, T13 with the ROBUST flags).
function p2_stacked(dec, I; robust = nothing)
    rows = Vector{Any}[]
    for h in (:H1, :H2, :H3, :H4, :H5, :H6)
        for r in _FLAT[h](dec[h], I)
            push!(rows, Any[string(h), r...])
        end
        robust === nothing || push!(rows, Any[string(h), "all", "all", "robust", _vstr(robust[h]), nothing, nothing])
    end
    return _tbl(["hypothesis", "source", "K index", "statistic", "value", "lower", "upper"], rows)
end

function p2_t12(gate)
    return _tbl(["indicator", "decided", "chaotic", "undecided", "not applicable", "share chaotic", "violated"],
                [[String(g.X), g.decided, g.chaotic, g.undecided, g.na, g.share, g.violated] for g in gate])
end
