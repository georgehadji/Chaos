# Amendment D89, stage 02: the reference computation on the synthetic study. It reads the chunk files with the real store reader (Experiments.read_chunk),
# not with the module, applies the rules of an_p2_ref.jl and the designed ranks, and returns the expected value of every table the pipeline must produce.

# orbit table for the set rules (E29): one NamedTuple per orbit of P and G
function p2_ref_orbits(st)
    d = st.design
    out = NamedTuple[]
    for id in 1:44
        if id <= 20
            o = AN_ORBITS[id]
            push!(out, (id = id, source = :panel, kindex = Int(o["k_index"]), K = get(o, "K", nothing), class = o["class"], control = o["control"] === true,
                        flags = id in d.drift ? ["drift"] : String[], rho = d.rho[id], lambda_ref = d.lref[id], lambda_ns = d.lns[id], in_SD = o["in_SD"] === true, in_SL = o["in_SL"] === true))
        else
            g = st.h4[id - 20]
            push!(out, (id = id, source = :grid, kindex = g.kindex, K = parse(Float64, g.K), class = "h4", control = false, flags = String[], rho = d.rho[id],
                        lambda_ref = d.lref[id], lambda_ns = d.lns[id], in_SD = false, in_SL = false))
        end
    end
    return out
end

function p2ref_sets(orbits, p)
    ids(f) = sort!([o.id for o in orbits if f(o)])
    dropped = NamedTuple[]
    keep = NamedTuple[]
    for o in orbits
        if o.control
            push!(dropped, (id = o.id, reason = "control"))
        elseif "drift" in o.flags
            push!(dropped, (id = o.id, reason = "drift"))
        elseif "excluded" in o.flags
            push!(dropped, (id = o.id, reason = "excluded"))
        elseif o.rho === :U
            push!(dropped, (id = o.id, reason = "unresolved"))
        else
            push!(keep, o)
        end
    end
    sort!(dropped; by = x -> x.id)
    E = Set(o.id for o in keep)
    inE(f) = sort!([o.id for o in keep if f(o)])
    C = filter(o -> o.rho === :C, keep)
    S2 = sort!([o.id for o in C if o.lambda_ref >= p.s2_lambda_min])
    Ts = sort!([o.id for o in C if o.lambda_ns <= p.f_s * o.lambda_ref])
    W2 = sort!(unique(vcat([o.id for o in C if o.lambda_ref <= p.w2_lambda_max], Ts)))
    gate = ids(o -> o.source === :panel && !o.control && o.kindex == p.gate_k_index && o.class in p.gate_classes)
    return (E_all = inE(_ -> true), E8 = inE(o -> o.in_SL || o.source === :grid), W2 = W2, S2 = S2, T_s = Ts,
            T_s_design = inE(o -> o.source === :panel && o.class == "S" && o.rho === :C), regular = inE(o -> o.rho === :R), gate = gate,
            S_D = inE(o -> o.in_SD), S_L = inE(o -> o.in_SL), overlap = sort!(intersect(S2, W2)), dropped = dropped)
end

# ---------------------------------------------------------------- values read from the chunk files at the recorded steps of E26
function p2_ref_items(st)
    items = Dict{Tuple, Any}()
    for (i, sp) in enumerate(st.specs), o in AX.read_chunk(p2_chunk(st, i)).orbits
        items[(o.id, sp.id, sp.precision, sp.horizon, sp.dindex)] = o
    end
    return items
end

# the value record `v` (see the test plan) of indicator X for orbit id at budget index k, from the item named by `src`; nothing if the orbit has no item
function p2ref_value(items, p, id, X::Symbol, k; prec = "Float64", horizon = nothing, dindex = 1)
    b = p.budgets[k + 1]
    colof(o, nm) = findfirst(==(nm), o.indicators[first(keys(o.indicators))].columns)
    if X === :LE
        o = get(items, (id, "L1", prec, something(horizon, p.N_1), dindex), nothing)
        o === nothing && return nothing
        t, lam = o.t, o.lyap["L1"][:, 1]
        n = p2ref_n(t, b); np = p2ref_nearest(t, n / 10)
        i, ip = searchsortedfirst(t, n), searchsortedfirst(t, np)
        return (n = n, lambda = lam[i], r = p2ref_ratio(p2ref_g(:LE, n, lam[i]), p2ref_g(:LE, np, lam[ip])), np = np)
    elseif X in (:SALI, :FLI, :MEGNO)
        o = get(items, (id, "tg", prec, something(horizon, p.N_1), dindex), nothing)
        o === nothing && return nothing
        S = o.indicators["tg"].series
        t = o.t
        n = p2ref_n(t, b); np = p2ref_nearest(t, n / 10)
        i, ip = searchsortedfirst(t, n), searchsortedfirst(t, np)
        c = X === :SALI ? "SALI.L" : X === :FLI ? "FLI.fli" : "MEGNO.Ybar20_1"
        j = findfirst(==(c), o.indicators["tg"].columns)
        r = p2ref_ratio(p2ref_g(X, n, S[i, j]), p2ref_g(X, np, S[ip, j]))
        extra = X === :SALI ? (sali = S[i, findfirst(==("SALI.sali"), o.indicators["tg"].columns)],) : X === :FLI ? (fli = S[i, j], n = n) : (ybar = S[i, j],)
        return merge((n = n, r = r, np = np), extra)
    elseif X === :RLI
        o = get(items, (id, "pr", prec, p.N_c, dindex), nothing)
        o === nothing && return nothing
        t = o.t; n = p2ref_n(t, b)
        return (n = n, rli = o.indicators["pr"].series[searchsortedfirst(t, n), findfirst(==("RLI.RLI"), o.indicators["pr"].columns)])
    elseif X in (:E1, :E2)
        o = get(items, (id, X === :E1 ? "fq1" : "pa1", prec, X === :E1 ? p.N_c : p.N_c + 1, 1), nothing)
        o === nothing && return nothing
        tb = o.indicators[X === :E1 ? "fq1" : "pa1"]
        r = findfirst(==("W=$(b ÷ 2)"), tb.rows)
        return (D = tb.final[r, findfirst(==("D"), tb.fcolumns)], flags = tb.flags[r], W = b ÷ 2)
    elseif X === :WBA
        o = get(items, (id, "bk1", prec, p.N_c, 1), nothing)
        o === nothing && return nothing
        tb = o.indicators["bk1"]
        rs = [findfirst(==("W=$(b ÷ 2) $nm"), tb.rows) for nm in P2_BKROWS]
        j = findfirst(==("dig"), tb.fcolumns)
        return (dig = [tb.final[r, j] for r in rs], flags = [tb.flags[r] for r in rs], W = b ÷ 2)
    end
end

# STR: the value D at n(b) of the stretching item
function p2ref_str(items, p, id, k)
    o = get(items, (id, "st", "Float64", p.N_c, 1), nothing)
    o === nothing && return nothing
    t = o.t
    n = p2ref_n(t, p.budgets[k + 1])
    return (D = o.indicators["st"].series[searchsortedfirst(t, n), 1], n = n)
end

# the whole expected table: verdicts R1 and R2, chaos scores, steps, detection, for every orbit and indicator that has an item
function p2_ref_all(st)
    p, d = st.p, st.design
    items = p2_ref_items(st)
    V1, V2, S, N, Det, Det2 = Dict{Tuple{Int, Symbol}, Vector{Symbol}}(), Dict{Tuple{Int, Symbol}, Vector{Symbol}}(), Dict{Tuple{Int, Symbol}, Vector{Any}}(),
                              Dict{Tuple{Int, Symbol}, Vector{Int}}(), Dict{Tuple{Int, Symbol}, Any}(), Dict{Tuple{Int, Symbol}, Any}()
    for id in 1:44, X in P2_X8
        vals = [p2ref_value(items, p, id, X, k) for k in 0:9]
        any(isnothing, vals) && continue
        V1[(id, X)] = [p2ref_verdict(X, :R1, v, p) for v in vals]
        V2[(id, X)] = [p2ref_verdict(X, :R2, v, p) for v in vals]
        S[(id, X)] = Any[p2ref_score(X, v) for v in vals]
        N[(id, X)] = [haskey(v, :n) ? v.n : v.W for v in vals]
        rho = d.rho[id]
        rho in (:C, :R) && (Det[(id, X)] = p2ref_detection(V1[(id, X)], rho, p); Det2[(id, X)] = p2ref_detection(V2[(id, X)], rho, p))     # detection under R1 and under R2 (E33)
    end
    for id in 1:44
        sv = [p2ref_str(items, p, id, k) for k in 0:9]
        any(isnothing, sv) && continue
        S[(id, :STR)] = Any[p2ref_score(:STR, v) for v in sv]
        N[(id, :STR)] = [v.n for v in sv]
    end
    return (items = items, V1 = V1, V2 = V2, S = S, N = N, Det = Det, Det2 = Det2)
end

# orbit records of the set statistics for a list of ids (rule :R1 or :R2)
function p2_ref_recs(ref, ids, X_list; rule = :R1)
    V = rule === :R1 ? ref.V1 : ref.V2
    return [(id = id, V = Dict(X => V[(id, X)] for X in X_list if haskey(V, (id, X))),
             T = Dict(X => (rule === :R1 ? ref.Det : ref.Det2)[(id, X)].censored for X in X_list if haskey(ref.Det, (id, X)))) for id in ids]
end

# ---------------------------------------------------------------- expected decisions of the pipeline (E32, E33), computed from the designed study by the reference rules
p2_lab(item, set, srule, rule, var) = "p2|$item|$set|$srule|$rule|$var"

function p2_ref_h6(ref, p, ids_d, ids_l, I_dir, I_prec; rule = :R1)
    items = ref.items
    dirs = [(id = id, v = Dict(X => [(v = p2ref_value(items, p, id, X, 9; horizon = p.N_dir, dindex = d); p2ref_verdict(X, rule, v, p)) for d in 1:5] for X in I_dir)) for id in ids_d]
    prec = [(id = id, f64 = Dict(X => (rule === :R1 ? ref.V1 : ref.V2)[(id, X)][10] for X in I_prec),
             d64 = Dict(X => p2ref_verdict(X, rule, p2ref_value(items, p, id, X, 9; prec = "Double64", horizon = p.N_dir), p) for X in I_prec)) for id in ids_l]
    function frac(rows, getA, getB, X)
        dec = [r for r in rows if getA(r, X) in (c_, r_) && getB(r, X) in (c_, r_)]
        return (length(dec) == 0 ? nothing : count(r -> getA(r, X) !== getB(r, X), dec) / length(dec), length(rows) - length(dec))
    end
    dd = Dict{Symbol, Any}(); dr = Dict{Symbol, Int}()
    for X in I_dir
        dec = [r for r in dirs if all(v -> v in (c_, r_), r.v[X])]
        dd[X] = isempty(dec) ? nothing : count(r -> !all(==(r.v[X][1]), r.v[X]), dec) / length(dec)
        dr[X] = length(dirs) - length(dec)
    end
    pd = Dict{Symbol, Any}(); pr = Dict{Symbol, Int}()
    for X in I_prec
        f, dropn = frac(prec, (r, X) -> r.f64[X], (r, X) -> r.d64[X], X)
        pd[X] = f; pr[X] = dropn
    end
    return (dir = (n = length(dirs), delta = dd, dropped = dr), prec = (n = length(prec), delta = pd, dropped = pr), rows = (dirs, prec))
end

# the six decisions for one set rule: sets named by `names` (W2, S2, T_s, T_s_design, and the direction and ladder subsets), indicators I, verdict rule `rule`
function p2_ref_hyp(st, ref, S, I; srule = "E_all", rule = :R1, var = "all", suffix = "", ids = nothing, h6 = true)
    p = st.p
    rec(ids_) = p2_ref_recs(ref, ids_, I; rule = rule)
    sel(v) = suffix == "" ? v : sort!(intersect(v, S.E8))
    W2, S2, Ts, Td = sel(S.W2), sel(S.S2), sel(S.T_s), S.T_s_design
    rW, rS, rT, rD = rec(W2), rec(S2), rec(Ts), rec(Td)
    nm(x) = x * suffix
    out = Dict{Symbol, Any}()
    # H1
    s1 = p2ref_stat(:H1, (rW,), p; I = I)
    b1 = p2ref_bootstrap(:H1, (rW,), p2_lab("H1", nm("W2"), srule, rule, var), p; I = I)
    lo1 = b1.interval === nothing ? nothing : b1.interval[:U].lo
    out[:H1] = merge(s1, (U_lo = lo1, U_hi = b1.interval === nothing ? nothing : b1.interval[:U].hi,
                          verdict = p2ref_decide(:H1, (n = s1.n, U = s1.U, U_lo = lo1, A_min = s1.A_min), p),
                          U_S2 = p2ref_stat(:H1, (rS,), p; I = I).U))
    # H2
    s2 = p2ref_stat(:H2, (rS,), p; I = I)
    b2 = p2ref_bootstrap(:H2, (rS,), p2_lab("H2", nm("S2"), srule, rule, var), p; I = I)
    lo2 = b2.interval === nothing ? nothing : b2.interval[:Q].lo
    out[:H2] = merge(s2, (Q_lo = lo2, Q_hi = b2.interval === nothing ? nothing : b2.interval[:Q].hi, verdict = p2ref_decide(:H2, (n = s2.n, Q = s2.Q, Q_lo = lo2), p)))
    # H3
    s3 = p2ref_stat(:H3, (rW, rS), p; I = I)
    b3 = p2ref_bootstrap(:H3, (rW, rS), (p2_lab("H3", nm("W2"), srule, rule, var), p2_lab("H3", nm("S2"), srule, rule, var)), p; I = I)
    Rlo = b3.interval === nothing ? nothing : Dict(X => b3.interval[Symbol(:R_, X)].lo for X in I)
    Rv = [s3.R[X] for X in I]
    out[:H3] = merge(s3, (R_lo = Rlo, verdict = Rlo === nothing ? :not_decidable : p2ref_decide(:H3, (n_w = s3.n_w, n_s = s3.n_s, R = Rv, R_lo = [Rlo[X] for X in I]), p),
                          passing = Rlo === nothing ? nothing : count(X -> s3.R[X] >= p.h3_r_min && Rlo[X] > p.h3_r_lower_gt, I)))
    # H4 (T_s and T_s_design)
    function h4(rr, setname)
        s = p2ref_stat(:H4, (rr,), p; I = I)
        b = p2ref_bootstrap(:H4, (rr,), p2_lab(setname == "T_s_design" ? "H4design" : "H4", nm(setname), srule, rule, var), p; I = I)
        lo = (b.interval === nothing || s.dphi === nothing) ? nothing : b.interval[:dphi].lo
        return merge(s, (dphi_lo = lo, verdict = p2ref_decide(:H4, (n = s.n, dphi = s.dphi, dphi_lo = lo), p)))
    end
    out[:H4] = h4(rT, "T_s")
    out[:H4design] = h4(rD, "T_s_design")
    # H5
    s5 = p2ref_stat(:H5, (rW,), p; I = I)
    b5 = p2ref_bootstrap(:H5, (rW,), p2_lab("H5", nm("W2"), srule, rule, var), p; I = I)
    lo5 = (b5.interval === nothing || s5.F_med === nothing) ? nothing : b5.interval[:F_med].lo
    out[:H5] = merge(s5, (F_med_lo = lo5, verdict = p2ref_decide(:H5, (n = s5.n, F_med = s5.F_med, F_med_lo = lo5), p)))
    # H6: the direction subset (four indicators of the tangent chain and LE) and the ladder subset (indicators with a Double64 item)
    if h6
        sd, sl = sel(S.S_D), sel(S.S_L)
        I_dir = [X for X in I if X in (:LE, :SALI, :FLI, :MEGNO)]
        I_prec = [X for X in I if X !== :RLI]
        hh = p2_ref_h6(ref, p, sd, sl, I_dir, I_prec; rule = rule)
        out[:H6] = merge(hh, (verdict = p2ref_decide(:H6, hh, p),))
    end
    return NamedTuple(out)
end

# ---------------------------------------------------------------- reported items (E31, E33, E35), computed from the designed study
# AUC per (indicator, budget index, K index, source): the orbits of E_all with rho C and R at that K that have an available score
function p2_ref_auc(st, ref, S)
    p = st.p
    orbs = Dict(o.id => o for o in p2_ref_orbits(st))
    out = Dict{Tuple{Symbol, Int, Int, String}, Any}()
    inE = Set(S.E_all)
    for X in vcat(P2_X8, :STR), k in 0:9, kidx in 1:2, src in ("pooled", "panel", "grid")
        cs, rs = Float64[], Float64[]
        for id in S.E_all
            o = orbs[id]
            o.kindex == kidx || continue
            (src == "pooled" || (src == "panel") == (o.source === :panel)) || continue
            sc = get(ref.S, (id, X), nothing)
            sc === nothing && continue
            v = sc[k + 1]
            v === nothing && continue
            o.rho === :C ? push!(cs, v) : push!(rs, v)
        end
        out[(X, k, kidx, src)] = (n_C = length(cs), n_R = length(rs), auc = p2ref_auc(cs, rs; min_set = p.min_set))
    end
    return out
end

# the sensitivities of T11 on S_L*: for each pair and budget, the number decided in both, A_XY, kappa, the median and the largest absolute difference of the score D (dig)
function p2_ref_sens(st, ref, S)
    p, it = st.p, ref.items
    out = Dict{Tuple{Symbol, Int}, Any}()
    col(o, id, name) = (tb = o.indicators[id]; findfirst(==(name), tb.fcolumns))
    for (pair, (a, b)) in ((:raw_vs_fq1, ("fq1_raw", "fq1")), (:qp_vs_fq1, ("fq1_qp", "fq1")), (:planar_vs_pa1, ("pa_planar", "pa1")), (:cos_vs_wba, ("bk1", "bk1")))
        for k in 0:9
            W = p.budgets[k + 1] ÷ 2
            va, vb, dd = Symbol[], Symbol[], Float64[]
            for id in S.S_L
                hz(s) = s in ("pa1", "pa_planar") ? p.N_c + 1 : p.N_c
                if pair === :cos_vs_wba
                    o = it[(id, "bk1", "Float64", p.N_c, 1)]; tb = o.indicators["bk1"]
                    rows = [findfirst(==("W=$W $nm"), tb.rows) for nm in P2_BKROWS]; j = findfirst(==("dig"), tb.fcolumns)
                    d1, f1 = tb.final[rows[1], j], tb.flags[rows[1]]
                    dmin = minimum(tb.final[r, j] for r in rows)
                    push!(va, p2ref_verdict(:WBA, :R1, (dig = [d1], flags = [f1]), p)); push!(vb, p2ref_verdict(:WBA, :R1, (dig = [tb.final[r, j] for r in rows], flags = [tb.flags[r] for r in rows]), p))
                    isfinite(d1) && isfinite(dmin) && push!(dd, abs(d1 - dmin))
                else
                    X = a == "pa_planar" ? :E2 : :E1
                    oa, ob = it[(id, a, "Float64", hz(a), 1)], it[(id, b, "Float64", hz(b), 1)]
                    ta, tb = oa.indicators[a], ob.indicators[b]
                    r = findfirst(==("W=$W"), ta.rows); j = findfirst(==("D"), ta.fcolumns)
                    push!(va, p2ref_verdict(X, :R1, (D = ta.final[r, j], flags = ta.flags[r]), p)); push!(vb, p2ref_verdict(X, :R1, (D = tb.final[r, j], flags = tb.flags[r]), p))
                    isfinite(ta.final[r, j]) && isfinite(tb.final[r, j]) && push!(dd, abs(ta.final[r, j] - tb.final[r, j]))
                end
            end
            t = p2ref_table(va, vb)
            out[(pair, k)] = (n = t.n, A = t.A, kappa = t.kappa, median_diff = isempty(dd) ? nothing : p2ref_median(dd), max_diff = isempty(dd) ? nothing : maximum(dd))
        end
    end
    return out
end

# packages (T15): the joint verdicts of FLI with RLI and of MEGNO with SALI at b_P and N_c on each set
function p2_ref_packages(st, ref, S)
    out = Dict{Tuple{String, Tuple{Symbol, Symbol}, Int}, Dict{Tuple{Symbol, Symbol}, Int}}()
    sets = Dict("W2" => S.W2, "S2" => S.S2, "T_s" => S.T_s, "regular" => S.regular)
    for (name, ids) in sets, pair in ((:FLI, :RLI), (:MEGNO, :SALI)), k in (P2_KP - 1, 9)
        ii = pair == (:FLI, :RLI) ? [i for i in ids if i in S.E8] : ids
        cnt = Dict{Tuple{Symbol, Symbol}, Int}()
        for id in ii
            haskey(ref.V1, (id, pair[1])) && haskey(ref.V1, (id, pair[2])) || continue
            key = (ref.V1[(id, pair[1])][k + 1], ref.V1[(id, pair[2])][k + 1])
            cnt[key] = get(cnt, key, 0) + 1
        end
        out[(name, pair, k)] = cnt
    end
    return out
end

# the rule-validity gate (T12): per indicator the counts on the gate set at k = 9
function p2_ref_gate(st, ref, S)
    p = st.p
    rows = NamedTuple[]
    for X in P2_X8
        vs = [get(ref.V1, (id, X), nothing) for id in S.gate]
        na = count(isnothing, vs)
        v9 = [v[10] for v in vs if v !== nothing]
        dec = count(in((c_, r_)), v9); ch = count(==(c_), v9); un = count(==(u_), v9)
        share = dec == 0 ? nothing : ch / dec
        push!(rows, (X = X, decided = dec, chaotic = ch, undecided = un, na = na, share = share, violated = share !== nothing && share > p.gate_max_share))
    end
    return rows
end

# the K = 0 line (T16): the largest circle distance of the stored frequency nu1 to p0 / (2 pi) modulo 1 over the gate orbits, per window
function p2_ref_k0(st, ref, S)
    p = st.p
    cd(a, b) = abs(mod(a - b + 0.5, 1.0) - 0.5)
    out = Dict{Int, NamedTuple}()
    for (w, W) in enumerate(p.windows)
        f = [cd(ref.items[(id, "fq1", "Float64", p.N_c, 1)].indicators["fq1"].final[w, 1], p2_nu0(AN_ORBITS[id]["p0"])) for id in S.gate]
        g = [cd(ref.items[(id, "pa1", "Float64", p.N_c + 1, 1)].indicators["pa1"].final[w, 1], p2_nu0(AN_ORBITS[id]["p0"])) for id in S.gate]
        out[W] = (fq1 = maximum(f), pa1 = maximum(g))
    end
    return out
end
