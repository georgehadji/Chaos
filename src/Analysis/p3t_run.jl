# Study p3-abstention (Amendment D95): HB at the pipeline level (E68, E69 a, b), HC rows (E70) and the power at the achieved sizes (E69 c, E70). Exact rationals; the one binary64 part is the
# power recomputation, which uses the simulations of D93 with the streams (pp-achieved-hb, run) and (pp-achieved-hc, run).

# The paired indicators of HB (E68) for the budget function `bx` (indicator -> budget or `nothing`): `cmp[(set, X)] = (rule, x)`, the coverage and the sizes of T_s and W2.
function p3t_hb_tables(rows, qt, bx; le::Bool = true)
    sets = (Ts = sort([r for r in rows if r.Ts]; by = r -> r.id), W2 = sort([r for r in rows if r.W2]; by = r -> r.id))
    voters = Symbol.(qt.voter_bits)
    cmp = Dict{Tuple{Symbol, Symbol}, Any}()
    for (A, Xs) in ((:Ts, voters), (:W2, le ? vcat(:LE, voters) : voters))
        for X in Xs
            b = bx(X)
            kx = b === nothing ? nothing : findfirst(==(b), qt.budgets)
            cmp[(A, X)] = (rule = Bool[r.out === :regular_supported for r in sets[A]], x = Bool[kx !== nothing && r.V[X][kx] === :regular for r in sets[A]])
        end
    end
    share(s) = isempty(s) ? 0 // 1 : P3T_R(count(r -> r.out !== :unresolved, s), length(s))
    return cmp, Dict{Symbol, Any}(:Ts => share(sets.Ts), :W2 => share(sets.W2)), (Ts = length(sets.Ts), W2 = length(sets.W2))
end

# The budget function of the matched budgets of `M` (E69 a): the ceiling rule, or the round-down budgets of the sensitivity (i).
p3t_budget_fn(M, down::Bool) = X -> down ? M.X[X].round_down : (X === :LE ? M.le_read : M.X[X].b)

"HB (E68) for a matching `M`: the decision with the comparators LE included, and the same without LE (reported, no effect)."
function p3t_hb_run(rows, qt, M, theta; down::Bool = false)
    bx = p3t_budget_fn(M, down)
    cmp, cov, sizes = p3t_hb_tables(rows, qt, bx)
    capped = Symbol[X for X in Symbol.(qt.voter_bits) if M.X[X].capped]
    hb = p3t_hb(cmp, cov, sizes, qt; capped)
    cmp12, cov12, _ = p3t_hb_tables(rows, qt, bx; le = false)
    return (hb = hb, cov = cov, sizes = sizes, without_le = p3t_hb(cmp12, cov12, sizes, qt; capped), bx = bx)
end

# The measured costs of the manifest (thread-ns per orbit-step per spec, P114) as the cost keys they replace; the specs without a value are listed.
function p3t_measured_costs(man)
    vals, miss = Dict{Symbol, Float64}(), String[]
    for (key, spec) in ((:cost_l1, "L1"), (:cost_tg, "tg"), (:cost_pa1, "pa1"), (:cost_bk1, "bk1"), (:cost_fq1, "fq1"))
        v = nothing
        for a in ("ref", "ind")
            s = get(get(get(get(man, "arms", Dict()), a, Dict()), "specs", Dict()), spec, nothing)
            s !== nothing && haskey(s, "measured_ns_per_orbit_step") && (v = Float64(s["measured_ns_per_orbit_step"]))
        end
        v === nothing ? push!(miss, spec) : (vals[key] = v)
    end
    return vals, miss
end

_p3t_sens(decidable, verdict, failing, reason) = (decidable = decidable, verdict = verdict, failing = failing, reason = reason)

# One sensitivity (E69 b): the HB comparisons recomputed with the costs `costs` (a NamedTuple or `nothing`) or the round-down budgets; not decidable if a comparator has no budget.
function _p3t_variant(rows, qt, theta, costs, down)
    M = p3_matched(costs === nothing ? qt : merge(qt, costs), theta.V)
    bad = Symbol[X for X in vcat(:LE, Symbol.(qt.voter_bits)) if p3t_budget_fn(M, down)(X) === nothing]
    isempty(bad) || return _p3t_sens(false, "not decidable", String[], down ? "no round-down budget for $(join(bad, ", "))" : "the matching caps the comparators $(join(bad, ", "))")
    r = p3t_hb_run(rows, qt, M, theta; down)
    return _p3t_sens(true, r.hb.decision, r.hb.failing, "")
end

"The sensitivities (i) to (iv) of E69 b: round-down budgets, pilot costs, costs measured in the test run, without LE. Each is reported only."
function p3t_sensitivities(rows, qt, theta, man, base)
    pilot = (cost_l1 = qt.pilot_cost_l1, cost_tg = qt.pilot_cost_tg, cost_pa1 = qt.pilot_cost_pa1, cost_bk1 = qt.pilot_cost_bk1, cost_fq1 = qt.pilot_cost_fq1)
    vals, miss = p3t_measured_costs(man)
    measured = isempty(miss) ? _p3t_variant(rows, qt, theta, NamedTuple(vals), false) :
               _p3t_sens(false, "not decidable", String[], "the manifest has no measured cost per orbit-step for the spec $(join(miss, ", ")) (P114)")
    return (round_down = _p3t_variant(rows, qt, theta, nothing, true), pilot = _p3t_variant(rows, qt, theta, pilot, false), measured = measured,
            without_le = _p3t_sens(true, base.without_le.decision, base.without_le.failing, ""))
end

# ---------------------------------------------------------------- HC rows (E70)

"The rows of HC from the orbits of D: the Float64 votes of V* at b_P, the share of voters whose Float64 and Double64 verdicts are both decided and differ (P104 a)."
function p3t_hc_rows(rows, vd, theta, kP)
    out = NamedTuple[]
    for r in sort(rows; by = x -> x.id)
        haskey(vd, r.id) || continue
        votes = [r.V[X][kP + 1] for X in theta.V]
        diff_ = count(X -> (a = r.V[X][kP + 1]; b = vd[r.id][X][kP + 1]; a !== :undecided && b !== :undecided && a != b), theta.V)
        push!(out, (id = r.id, rho = r.rho, nC = count(==(:chaotic), votes), nR = count(==(:regular), votes), sP = P3T_R(diff_, length(theta.V))))
    end
    return out
end

# ---------------------------------------------------------------- the power at the achieved sizes (E69 c, E70)

# The model of the record (the rates written by D93, 12 significant digits) in the form of the simulation of E56 c.
function p3t_record_model(rec, qt)
    m = rec["pp_hb"]["model"]
    rate(d) = (p = Float64(d["p"]), d = Float64(d["d"]))
    cs(d) = (T = (p = Float64(d["T"]["p"]), s = Float64(d["T"]["s"])), RW = (p = Float64(d["RW"]["p"]), s = Float64(d["RW"]["s"])))
    return (rule = Dict(:T => rate(m["rule"]["T"]), :RW => rate(m["rule"]["RW"])), comp = Dict(Symbol(X) => cs(v) for (X, v) in m["comp"]), ts = Symbol.(qt.voter_bits), w2 = Symbol.(qt.indicators_i7))
end

"The power of the full HB decision at the final sizes `nT`, `nRW` (E69 c): `pp_runs` runs with the streams (pp-achieved-hb, run index); counts are integers, so the threads change nothing."
function p3t_power_hb(model, nT, nRW, qt, threads)
    S = (_p3_stratum(model, :T), _p3_stratum(model, :RW))
    D = _p3_decision_data(model, qt, nT + nRW)
    blocks = _p3_blocks(qt.pp_runs, threads)
    ok = _pmap(b -> _p3t_hb_block(S, D, qt, nT, nRW, b), blocks, threads)
    return merge((successes = sum(ok), runs = qt.pp_runs), _p3_power(sum(ok), qt.pp_runs), (planned = qt.power_hb_planned,))
end

function _p3t_hb_block(S, D, qt, nT, nRW, runs)
    cnt, dec, ok = zeros(Int, 2, length(P3_COMPS), 4), zeros(Int, 2), 0
    for i in runs
        _p3_orbits!(cnt, dec, p3_stream(qt.test_seed, "pp-achieved-hb", i), S, nT, nRW)
        ok += _p3_decide(cnt, dec, D, nT, nRW, nothing)
    end
    return ok
end

"The power of the HC decision at the achieved failures and non-failures (E70, P103): the simulation of D93 E57 with the AUCs of the record and the streams (pp-achieved-hc, run index)."
function p3t_power_hc(rec, nf, nn, qt, threads)
    hc = rec["pp_hc"]
    dI, dP = p3_binormal_shift(Float64(hc["auc_i_used"])), p3_binormal_shift(Float64(hc["auc_p"]))
    blocks = _p3_blocks(qt.pp_runs, threads)
    dc = p3_rat(qt.delta_c)
    ok = _pmap(b -> _p3t_hc_block(qt, dI, dP, nf, nn, dc, b), blocks, threads)
    return merge((successes = sum(ok), runs = qt.pp_runs), _p3_power(sum(ok), qt.pp_runs), (planned = qt.power_hc_planned,))
end

function _p3t_hc_block(qt, dI, dP, nf, nn, dc, runs)
    ok = 0
    for i in runs
        cI, cP, d = _p3_hc_run(p3_stream(qt.test_seed, "pp-achieved-hc", i), dI, dP, nf, nn)
        AI, AP = cI / d, cP / d
        SE = sqrt(p3_hanley_mcneil(AI, nf, nn) + p3_hanley_mcneil(AP, nf, nn))
        ((cI - cP) * denominator(dc) >= numerator(dc) * d && (AI - AP) - qt.z_one_sided * SE > 0) && (ok += 1)
    end
    return ok
end
