# Study p2-indicators (Amendment D89): verdicts and scores of every orbit (E26 to E31), the decisions H1 to H6 in this order, the I_8 extension with the ROBUST flags,
# the R2 version and the rule-validity gate with its second versions (E32, E33).

const P2_ALL = vcat(P2_I8, :STR)

# E33 says the R2 version recomputes every statistic with the verdicts of R2. The stage 02 tests pin a narrower reading: under R2 only the verdict series change; the detection times
# (M_X, Q, R_X) and the verdicts of H6 stay those of R1 (an_p2_refstore.jl, `p2_ref_recs` and `p2_ref_h6`). `true` gives the reading of the spec; see the implementation notes.
const P2_R2_REDETECT = false

"""
    p2_verdict_tables(ctx, orbits, sets) -> NamedTuple

For every orbit and indicator that has an item: the verdicts of table R1 and R2 (`verdicts[(id, X, rule)]`), the chaos scores (`scores[(id, X)]`, also `:STR`), the
steps or windows used (`steps[(id, X)]`) and, for the orbits of E_all, the detection of E28 under both rules (`detection[(id, X)]` is the R1 one, `detection2` the R2 one).
"""
function p2_verdict_tables(ctx, orbits, sets)
    p = ctx.p
    verdicts = Dict{Tuple{Int, Symbol, Symbol}, Vector{Symbol}}()
    scores = Dict{Tuple{Int, Symbol}, Vector{Any}}()
    steps = Dict{Tuple{Int, Symbol}, Vector{Int}}()
    for o in orbits, X in P2_ALL
        vals = p2_values(ctx, o.id, X)
        vals === nothing && continue
        steps[(o.id, X)] = Int[p2_step_of(v) for v in vals]
        scores[(o.id, X)] = Any[p2_score(X, v) for v in vals]
        X === :STR && continue
        for rule in (:R1, :R2)
            verdicts[(o.id, X, rule)] = Symbol[p2_verdict(X, rule, v, p) for v in vals]
        end
    end
    rho = Dict(o.id => o.rho for o in orbits)
    det1, det2 = Dict{Tuple{Int, Symbol}, Any}(), Dict{Tuple{Int, Symbol}, Any}()
    for id in sets.E_all, X in P2_I8
        haskey(verdicts, (id, X, :R1)) || continue
        det1[(id, X)] = p2_detection(verdicts[(id, X, :R1)], rho[id], p)
        det2[(id, X)] = p2_detection(verdicts[(id, X, :R2)], rho[id], p)
    end
    return (verdicts = verdicts, scores = scores, steps = steps, detection = det1, detection2 = det2)
end

# Records of E32 for a list of orbit ids: the verdicts of the rule and the censored detection times of that rule.
function p2_records(vt, ids, I, rule)
    det = (rule === :R1 || !P2_R2_REDETECT) ? vt.detection : vt.detection2
    return [(id = id, V = Dict{Symbol, Vector{Symbol}}(X => vt.verdicts[(id, X, rule)] for X in I if haskey(vt.verdicts, (id, X, rule))),
             T = Dict{Symbol, Int}(X => det[(id, X)].censored for X in I if haskey(det, (id, X)))) for id in ids]
end

# Rows of H6 (E32): the verdicts at the last budget in each direction (L1 and tangent items of horizon N_dir) and in the two precisions (Float64 primary against Double64).
function p2_h6_rows(ctx, vt, ids_d, ids_l, I, rule)
    p = ctx.p
    rule = P2_R2_REDETECT ? rule : :R1
    v9(id, X; kw...) = (vals = p2_values(ctx, id, X; kw...); p2_verdict(X, rule, vals[end], p))
    Idir = [X for X in I if X in (:LE, :SALI, :FLI, :MEGNO)]
    Iprec = [X for X in I if X !== :RLI]
    hz(X) = X in (:LE, :SALI, :FLI, :MEGNO) ? (horizon = p.N_dir,) : (;)
    dirs = [(id = id, v = Dict(X => Symbol[v9(id, X; horizon = p.N_dir, dir = d) for d in p.directions] for X in Idir)) for id in ids_d]
    prec = [(id = id, f64 = Dict(X => vt.verdicts[(id, X, rule)][end] for X in Iprec), d64 = Dict(X => v9(id, X; prec = "Double64", hz(X)...) for X in Iprec)) for id in ids_l]
    return dirs, prec
end

_lower(iv, k) = iv === nothing ? nothing : iv[k].lo
_upper(iv, k) = iv === nothing ? nothing : iv[k].hi

# The six decisions on one set rule, indicator set and verdict rule; the bootstraps run over the labels (jobs), the thread count does not change a number.
function p2_decisions(ctx, vt, S, I; srule = "E_all", rule = :R1, var = "all", suffix = "", threads = 1)
    p = ctx.p
    sel(v) = isempty(suffix) ? v : sort!(intersect(v, S.E8))
    rW, rS, rT = (p2_records(vt, sel(x), I, rule) for x in (S.W2, S.S2, S.T_s))
    rD = p2_records(vt, S.T_s_design, I, rule)
    lab(item, set) = p2_label(item, set * suffix, srule, rule, var)
    jobs = [(:H1, (rW,), lab(:H1, "W2")), (:H2, (rS,), lab(:H2, "S2")), (:H3, (rW, rS), (lab(:H3, "W2"), lab(:H3, "S2"))), (:H4, (rT,), lab(:H4, "T_s")),
            (:H4, (rD,), lab(:H4design, "T_s_design")), (:H5, (rW,), lab(:H5, "W2"))]
    b = _pmap(j -> p2_bootstrap(j[1], j[2], j[3], p; I = I), jobs, threads)
    h1 = _dec_h1(p, I, rW, rS, b[1].interval)
    h2 = _dec_h2(p, I, rS, b[2].interval)
    h3 = _dec_h3(p, I, rW, rS, b[3].interval)
    h4 = merge(_dec_h4(p, I, rT, b[4].interval), (design = _dec_h4(p, I, rD, b[5].interval),))
    h5 = _dec_h5(p, I, rW, b[6].interval)
    dirs, prec = p2_h6_rows(ctx, vt, sel(S.S_D), sel(S.S_L), I, rule)
    return (H1 = h1, H2 = h2, H3 = h3, H4 = h4, H5 = h5, H6 = p2_decide(:H6, p2_h6_stats(dirs, prec, p), p))
end

function _dec_h1(p, I, rW, rS, iv)
    s = p2_stat(:H1, (rW,), p; I = I)
    return p2_decide(:H1, merge(s, (U_lo = _lower(iv, :U), U_hi = _upper(iv, :U), U_S2 = p2_stat(:H1, (rS,), p; I = I).U)), p)
end

function _dec_h2(p, I, rS, iv)
    return p2_decide(:H2, merge(p2_stat(:H2, (rS,), p; I = I), (Q_lo = _lower(iv, :Q), Q_hi = _upper(iv, :Q))), p)
end

function _dec_h3(p, I, rW, rS, iv)
    s = p2_stat(:H3, (rW, rS), p; I = I)
    lo = iv === nothing ? nothing : Dict{Symbol, Any}(X => iv[Symbol(:R_, X)].lo for X in I)
    hi = iv === nothing ? nothing : Dict{Symbol, Any}(X => iv[Symbol(:R_, X)].hi for X in I)
    d = p2_decide(:H3, merge(s, (R_lo = lo, R_hi = hi)), p)
    passing = lo === nothing ? nothing : count(X -> s.R[X] !== nothing && s.R[X] >= p.h3_r_min && lo[X] > p.h3_r_lower_gt, I)
    return merge(d, (passing = passing,))
end

function _dec_h4(p, I, rT, iv)
    s = p2_stat(:H4, (rT,), p; I = I)
    return p2_decide(:H4, merge(s, (dphi_lo = s.dphi === nothing ? nothing : _lower(iv, :dphi), dphi_hi = s.dphi === nothing ? nothing : _upper(iv, :dphi))), p)
end

function _dec_h5(p, I, rW, iv)
    s = p2_stat(:H5, (rW,), p; I = I)
    return p2_decide(:H5, merge(s, (F_med_lo = s.F_med === nothing ? nothing : _lower(iv, :F_med), F_med_hi = s.F_med === nothing ? nothing : _upper(iv, :F_med))), p)
end

"""
    p2_gate(vt, S, p) -> Vector{NamedTuple}

The rule-validity gate (E33, T12): per indicator of I_8 the orbits of the gate set decided and chaotic at the last budget, undecided and without the item; the gate is violated
if the share chaotic among the decided exceeds `gate_max_share`.
"""
function p2_gate(vt, S, p)
    rows = NamedTuple[]
    for X in P2_I8
        vs = [get(vt.verdicts, (id, X, :R1), nothing) for id in S.gate]
        v9 = [v[end] for v in vs if v !== nothing]
        dec, ch = count(_decided, v9), count(==(:chaotic), v9)
        share = dec == 0 ? nothing : ch / dec
        push!(rows, (X = X, decided = dec, chaotic = ch, undecided = count(==(:undecided), v9), na = count(isnothing, vs), share = share,
                     violated = share !== nothing && share > p.gate_max_share))
    end
    return rows
end

"""
    p2_hypotheses(ctx, vt, S, threads) -> NamedTuple

Primary decisions (I_7 on E_all), the extension (I_8 on E8) and the ROBUST flags, the R2 versions, the gate and the second versions without each violating indicator (E32, E33).
"""
function p2_hypotheses(ctx, vt, S, threads)
    p = ctx.p
    dec = p2_decisions(ctx, vt, S, P2_I7; threads = threads)
    ext = p2_decisions(ctx, vt, S, P2_I8; srule = "E8", suffix = "_8", threads = threads)
    robust = NamedTuple{(:H1, :H2, :H3, :H4, :H5, :H6)}(Tuple(p2_robust(dec[h].verdict, ext[h].verdict) for h in (:H1, :H2, :H3, :H4, :H5, :H6)))
    r2 = p2_decisions(ctx, vt, S, P2_I7; rule = :R2, threads = threads)
    r2ext = p2_decisions(ctx, vt, S, P2_I8; srule = "E8", suffix = "_8", rule = :R2, threads = threads)
    gate = p2_gate(vt, S, p)
    second = Dict{Symbol, Any}()
    for g in gate
        g.violated || continue
        second[g.X] = (decisions = p2_decisions(ctx, vt, S, [X for X in P2_I7 if X !== g.X]; var = "minus-$(g.X)", threads = threads),
                       extension = p2_decisions(ctx, vt, S, [X for X in P2_I8 if X !== g.X]; srule = "E8", suffix = "_8", var = "minus-$(g.X)", threads = threads))
    end
    return (decisions = dec, extension = ext, robust = robust, r2 = r2, r2_extension = r2ext, gate = gate, second = second)
end
