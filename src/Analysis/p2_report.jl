# Study p2-indicators (Amendment D89): the quantities reported without a decision rule (E31, E35): AUC, regular range, leaving budget, sensitivities, packages, flags,
# controls, the K = 0 line and the cost. Every number is read from the verdict tables or from stored values; nothing is recomputed from an orbit (I18).

_kindexes(orb, ids) = sort!(unique(orb[id].kindex for id in ids))

"AUC(X, k, K) with the class counts, pooled over the panel and the grid and by source (E31); the key is `(X, k, K index, source)`, k from 0."
function p2_auc_table(vt, orbits, S, p)
    orb = Dict(o.id => o for o in orbits)
    out = Dict{Tuple{Symbol, Int, Int, String}, NamedTuple}()
    for X in P2_ALL, k in 0:(length(p.budgets) - 1), kidx in _kindexes(orb, S.E_all), src in ("pooled", "panel", "grid")
        cs, rs = Float64[], Float64[]
        for id in S.E_all
            o = orb[id]
            (o.kindex == kidx && (src == "pooled" || (src == "panel") == (o.source === :panel))) || continue
            sc = get(vt.scores, (id, X), nothing)
            (sc === nothing || sc[k + 1] === nothing) && continue
            push!(o.rho === :C ? cs : rs, sc[k + 1])
        end
        out[(X, k, kidx, src)] = (n_C = length(cs), n_R = length(rs), auc = p2_auc(cs, rs; min_set = p.min_set))
    end
    return out
end

"""
    p2_range_tables(vt, orbits, S, p) -> (range, leaving, below_median)

Regular range of the scores of the regular set at each K index and budget, the leaving budget of each sticky orbit and the share of the sticky orbits below the regular median
at the last budget (M10, P25), for the indicators of I_7 and STR.
"""
function p2_range_tables(vt, orbits, S, p)
    orb = Dict(o.id => o for o in orbits)
    nb, KI = length(p.budgets), _kindexes(orb, S.E_all)
    Xs = vcat(P2_I7, :STR)
    range = Dict{Tuple{Symbol, Int, Int}, Any}()
    for X in Xs, kidx in KI, k in 0:(nb - 1)
        sc = Float64[vt.scores[(id, X)][k + 1] for id in S.regular if orb[id].kindex == kidx && haskey(vt.scores, (id, X)) && vt.scores[(id, X)][k + 1] !== nothing]
        range[(X, kidx, k)] = p2_regular_range(sc, p)
    end
    leaving = Dict{Tuple{Int, Symbol}, Any}()
    for id in S.T_s, X in Xs
        haskey(vt.scores, (id, X)) || continue
        up = [range[(X, orb[id].kindex, k)] === nothing ? nothing : range[(X, orb[id].kindex, k)].hi for k in 0:(nb - 1)]
        leaving[(id, X)] = any(isnothing, up) ? nothing : p2_leaving_budget(vt.scores[(id, X)], up, p)
    end
    below = Dict{Tuple{Symbol, Int}, Any}()
    for X in Xs, kidx in KI
        sc = Float64[vt.scores[(id, X)][nb] for id in S.T_s if orb[id].kindex == kidx && haskey(vt.scores, (id, X)) && vt.scores[(id, X)][nb] !== nothing]
        rg = range[(X, kidx, nb - 1)]
        below[(X, kidx)] = (rg === nothing || isempty(sc)) ? nothing : p2_share_below(sc, rg.median)
    end
    return range, leaving, below
end

# One pair of the sensitivities (T11): the verdicts of two versions of an indicator at one budget on the orbits of S_L*, the agreement and the differences of the score.
function _sens_pair(ctx, ids, a, b, k)
    p = ctx.p
    W = p.budgets[k + 1] ÷ 2
    va, vb, dd = Symbol[], Symbol[], Float64[]
    for id in ids
        if a == "cos"
            it = ctx.idx[p2_key(id, "bk1", "Float64", p.N_c, 1, nothing)]
            j = _p2_col(it, it.fcolumns, "dig")
            r4 = [findfirst(==("W=$W $nm"), it.rows) for nm in P2_BK_ROWS]
            any(isnothing, r4) && _fail("window W=$W of $(_lbl(it)) has no row for an observable")
            d1, dmin = it.final[r4[1], j], minimum(it.final[r, j] for r in r4)
            push!(va, p2_verdict(:WBA, :R1, (dig = [d1], flags = [it.flagwords[r4[1]]]), p))
            push!(vb, p2_verdict(:WBA, :R1, (dig = Float64[it.final[r, j] for r in r4], flags = Int64[it.flagwords[r] for r in r4]), p))
            isfinite(d1) && isfinite(dmin) && push!(dd, abs(d1 - dmin))
        else
            X = a == "pa_planar" ? :E2 : :E1
            hor = X === :E2 ? p.N_c + 1 : p.N_c
            ia, ib = ctx.idx[p2_key(id, a, "Float64", hor, 1, nothing)], ctx.idx[p2_key(id, b, "Float64", hor, 1, nothing)]
            r = findfirst(==("W=$W"), ia.rows)
            r === nothing && _fail("window row W=$W is missing in $(_lbl(ia))")
            j = _p2_col(ia, ia.fcolumns, "D")
            push!(va, p2_verdict(X, :R1, (D = ia.final[r, j], flags = ia.flagwords[r]), p))
            push!(vb, p2_verdict(X, :R1, (D = ib.final[r, j], flags = ib.flagwords[r]), p))
            isfinite(ia.final[r, j]) && isfinite(ib.final[r, j]) && push!(dd, abs(ia.final[r, j] - ib.final[r, j]))
        end
    end
    t = p2_contingency(va, vb)
    return (n = t.n, A = t.A, kappa = t.kappa, median_diff = isempty(dd) ? nothing : p2_median(dd), max_diff = isempty(dd) ? nothing : maximum(dd))
end

const P2_BK_ROWS = ["cos(1 x[1])", "sin(1 x[1])", "cos(1 x[2])", "sin(1 x[2])"]

"T11 (E35): `(pair, k)` => `(n, A, kappa, median_diff, max_diff)` on S_L* for the four pairs of versions."
function p2_sensitivity(ctx, S)
    out = Dict{Tuple{Symbol, Int}, NamedTuple}()
    for (pair, a, b) in ((:raw_vs_fq1, "fq1_raw", "fq1"), (:qp_vs_fq1, "fq1_qp", "fq1"), (:planar_vs_pa1, "pa_planar", "pa1"), (:cos_vs_wba, "cos", "bk1")), k in 0:(length(ctx.p.budgets) - 1)
        out[(pair, k)] = _sens_pair(ctx, S.S_L, a, b, k)
    end
    return out
end

"T15 (E35): the counts of the joint verdicts of FLI with RLI (on E8) and MEGNO with SALI at b_P and the last budget on W2, S2, T_s and the regular set; and the same by reference class."
function p2_packages(vt, orbits, S, p)
    rho = Dict(o.id => o.rho for o in orbits)
    sets = (("W2", S.W2), ("S2", S.S2), ("T_s", S.T_s), ("regular", S.regular))
    out = Dict{Tuple{String, Tuple{Symbol, Symbol}, Int}, Dict{Tuple{Symbol, Symbol}, Int}}()
    byrho = Dict{Tuple{String, Tuple{Symbol, Symbol}, Int}, Dict{Tuple{Symbol, Symbol, Symbol}, Int}}()
    for (name, ids) in sets, pair in ((:FLI, :RLI), (:MEGNO, :SALI)), k in (_kp(p) - 1, length(p.budgets) - 1)
        cnt, cr = Dict{Tuple{Symbol, Symbol}, Int}(), Dict{Tuple{Symbol, Symbol, Symbol}, Int}()
        for id in (pair == (:FLI, :RLI) ? [i for i in ids if i in S.E8] : ids)
            (haskey(vt.verdicts, (id, pair[1], :R1)) && haskey(vt.verdicts, (id, pair[2], :R1))) || continue
            key = (vt.verdicts[(id, pair[1], :R1)][k + 1], vt.verdicts[(id, pair[2], :R1)][k + 1])
            cnt[key] = get(cnt, key, 0) + 1
            cr[(key..., rho[id])] = get(cr, (key..., rho[id]), 0) + 1
        end
        out[(name, pair, k)], byrho[(name, pair, k)] = cnt, cr
    end
    return out, byrho
end

"The controls (E35, T16): lambda_1(N_1) against the analytic exponent with the stored deviation and scope. The recomputed deviation must equal the stored one to 1e-12 (E36 n)."
function p2_controls(ctx, vt = nothing)
    p = ctx.p
    rows = NamedTuple[]
    for o in ctx.panel
        o["control"] === true || continue
        id = Int(o["id"])
        it = ctx.idx[p2_key(id, "L1", "Float64", p.N_1, 1, p.directions[1])]
        lam = it.lam[_p2_row(it, p.N_1, "control")]
        dev, stored = abs(lam - Float64(o["exponent"])), get(it.attrs, "control_deviation", nothing)
        stored === nothing || abs(dev - Float64(stored)) <= 1e-12 || _fail("control_deviation $stored of control $id differs from the deviation $dev recomputed from the stored series")
        v9 = vt === nothing ? Dict{Symbol, Symbol}() : Dict{Symbol, Symbol}(X => vt.verdicts[(id, X, :R1)][end] for X in P2_I8 if haskey(vt.verdicts, (id, X, :R1)))
        push!(rows, (id = id, scope = haskey(it.attrs, "control_scope") ? String(it.attrs["control_scope"]) : nothing, exponent = Float64(o["exponent"]), lambda = lam,
                     deviation = dev, stored_deviation = stored === nothing ? nothing : Float64(stored), verdicts9 = v9))
    end
    return rows
end

"The K = 0 line (E35, T16): per window, the largest circle distance of the stored frequency nu1 of fq1 and of pa1 to p0 / (2 pi) modulo 1 over the gate orbits."
function p2_k0line(ctx, S)
    p = ctx.p
    cd(a, b) = abs(mod(a - b + 0.5, 1.0) - 0.5)
    p0 = Dict(Int(o["id"]) => Float64(o["p0"]) for o in ctx.panel)
    out = Dict{Int, NamedTuple}()
    for W in p.windows
        dist(spec, hor) = maximum(cd(it.final[findfirst(==("W=$W"), it.rows), _p2_col(it, it.fcolumns, "nu1")], mod(p0[id] / (2 * pi), 1.0))
                                  for id in S.gate for it in (ctx.idx[p2_key(id, spec, "Float64", hor, 1, nothing)],))
        out[W] = (fq1 = dist("fq1", p.N_c), pa1 = dist("pa1", p.N_c + 1))
    end
    return out
end

"The cost lines of the manifest (T17): per spec the measured and budgeted ns per orbit-step, the bytes per orbit and the transform seconds. Nothing is measured here."
function p2_cost(man)
    out = Dict{String, NamedTuple}()
    for (id, e) in get(man, "indicators", Dict())
        secs = 0.0
        for c in man["chunks"]
            m = match(r"^method=(\S+)", String(get(c, "key", "")))
            m !== nothing && m.captures[1] == id && (secs += Float64(get(c, "transform_seconds", 0.0)))
        end
        line = measured = nothing
        for a in sort!(collect(keys(get(man, "arms", Dict()))))
            s = get(get(man["arms"][a], "specs", Dict()), id, nothing)
            s === nothing && continue
            line, measured = s["budget_line_ns"], s["measured_ns_per_orbit_step"]
            break
        end
        out[String(id)] = (measured_ns = measured, budget_line_ns = line, bytes_per_orbit = get(get(e, "bytes_per_orbit", Dict()), "total", nothing), transform_seconds = secs)
    end
    return out
end

"Flag counts per spec (T16): the sums over the chunks of the orbits with each flag, recomputed from the stored words (they equal the manifest's, E36 p)."
function p2_flag_counts(byc, man)
    out = Dict{String, Dict{String, Int}}()
    for c in man["chunks"], it in get(byc, Int(c["ordinal"]), NamedTuple[])
        it.kind === nothing && continue
        d = get!(out, it.id, Dict{String, Int}())
        it.flagwords === nothing && continue
        for (j, nm) in enumerate(it.bits)
            any(w -> (w & (Int64(1) << (j - 1))) != 0, it.flagwords) && (d[nm] = get(d, nm, 0) + 1)
        end
    end
    return out
end
