# Study p3-abstention (Amendment D93): the content of the calibration record of E58 (the items; the text is `p3_json`).

const P3_READINGS = [
    "P43: the command calibrate exists for this study only; analyze refuses it until the test-side amendment",
    "P44: the parameter table and the output are in stage 01 of the study; the table has no machine key",
    "P45: the L1 item of a panel orbit in the L1 against L5 test is the P1 lifted item at N_1",
    "P46: the diagnostics condition reads the flag diagnostics exceeded on the L5 and L1 items read; the P1 base labelling is run in memory for the change counts",
    "P47: the verdict rules of table R1 and the four label numbers come from the p2 table, validated against the p2 preregistration",
    "P48: C3 to C6 use all of Cal, K = 0 orbits included; e and the unresolved share on Cal with K > 0 are reported",
    "P49: a set below min_set gives the outcome sets too small; an element with no decided orbit is not feasible",
    "P50: every count-based comparison, J, Theta_tol, the ceiling rule and the order of C6 are exact rationals of the decimal texts",
    "P51: the optimism resamples the Cal orbits with replacement; an empty F skips a resample; drawn gate members count with multiplicity",
    "P52: R_W rates by counts by difference; shared_default when the rule has no false-regular orbit; LE is modelled on T_s too",
    "P53: the decided probability is clamped to the unit interval; the comparator probability is 0 if the rule rate is 1",
    "P54: the type-I check sets the rule rate to the rate of one comparator at a time",
    "P55: round half away from zero; n_W of the path is n_T plus n_RW; the target is ceil5 of the product",
    "P56: the n_f search does not stop at the expected deliverable number",
    "P57: decided means chaotic-supported or regular-supported",
    "P58: one named stream per purpose and unit; xoshiro256++ seeded by SplitMix64; a pair unit is written a,b; an integer in 1 to m by modulo with rejection; false-regular iff u < p",
    "P59: the program identity is the commit, whether the tree was clean, and the SHA-256 of the parameter table and of this record",
    "P62: the hash row of the preregistration is a self-consistency check; the recomputed hash is recorded",
    "P63: R_ab is evaluated on the Cal orbits that carry Double64 items and reported",
    "P64: HC uses the decided Float64 votes of V* at b_P without the stability filter; Hanley and McNeil variances at the empirical AUCs; a tie is a failure",
    "P65, P66: no limit on skipped resamples; a set below min_set inside a resample does not skip it",
    "P67: the n_f search cap is the smallest listed value at least the worst-case deliverable failures",
    "P68: a gate orbit with a missing verdict item refuses",
    "A5: with no failure or no non-failure on Cal the AUC is undefined and PP for HC is not run",
    "An item listed in not_applicable that a Cal orbit lacks gives undecided verdicts"]

_p3_counts(ids, rows, f) = (pooled = length(ids), panel = count(i -> f(i) === :panel, ids), grid = count(i -> f(i) === :grid, ids))

# The sets of E50 by source and the counts by K index of the orbits of Cal.
function _p3_sets_block(res)
    src = Dict(r.id => r.source for r in res.rows)
    S = res.sets
    blk(ids) = _p3_counts(ids, res.rows, i -> src[i])
    byk = [(source = s, k_index = k, count = count(r -> r.in_cal && r.source === s && r.kindex == k, res.rows)) for s in (:panel, :grid) for k in sort(unique(r.kindex for r in res.rows if r.source === s))]
    return (cal = blk(S.cal), S2 = blk(S.S2), T_s = blk(S.T_s), W2 = blk(S.W2), R_W = blk(S.R_W), regular = blk(S.regular), gate = blk(S.gate), overlap = length(S.overlap), cal_by_k = byk)
end

function _p3_relabel_block(res)
    ch = res.relabel.changes
    table(d) = [(k_index = k, old = o, new = n, count = c) for ((k, o, n), c) in sort(collect(d); by = x -> (x[1][1], string(x[1][2]), string(x[1][3])))]
    by_reason = Dict{String, Int}()
    for d in res.dropped; by_reason[d.reason] = get(by_reason, d.reason, 0) + 1; end
    rho = [(source = s, k_index = k, rho = r, count = count(x -> x.source === s && x.kindex == k && x.rho === r, res.relabel.rows)) for s in (:panel, :grid)
           for k in sort(unique(x.kindex for x in res.relabel.rows if x.source === s)) for r in (:C, :R, :U)]
    axes = [(instability = i, structure = s, count = count(x -> x.instability == i && x.structure == s, res.relabel.rows)) for (i, s) in unique((x.instability, x.structure) for x in res.relabel.rows)]
    return (rho_new = filter(x -> x.count > 0, rho), axes = axes, changes_panel = table(ch.panel), changes_grid = table(ch.grid), panel_changed = ch.panel_changed, grid_changed = ch.grid_changed,
            missing_data = [(id = m.id, source = m.source, k_index = m.kindex, item = m.item) for m in res.relabel.missing],
            dropped = (total = length(res.dropped), by_reason = by_reason, rows = [(id = d.id, source = d.source, k_index = d.kindex, reason = d.reason) for d in res.dropped]))
end

_p3_theta(t) = (index = t.index, V = t.V, qC = t.qC, qR = t.qR, w = t.w, v = t.v)

function _p3_selection_block(res, q)
    sel = res.selection
    base = (feasible = length(sel.feasible), tol = length(sel.tol), index = sel.star, star = nothing, J_min = sel.J_min, theta = nothing)
    sel.star === nothing && return base
    t, m = sel.thetas[sel.star], res.matched
    star = merge(_p3_theta(t), (cost = sel.stats[sel.star].cost, matched = Dict(X => (b = i.b, capped = i.capped, round_down = i.round_down, round_down_none = i.round_down_none) for (X, i) in m.X),
                               le_read = m.le_read, le_cost_unmatched = m.le_unmatched, d64_costs = p3_d64_costs(q, t.V)))
    return merge(base, (star = star, theta = (V = t.V, qC = t.qC, qR = t.qR, w = t.w, v = t.v)))
end

function _p3_hb_block(res)
    res.hb === nothing && return nothing
    h = res.hb
    return (rates = h.rates, model = h.model, path = [(j = x.j, rho = x.rho, nT = x.nT, nRW = x.nRW, nW = x.nW, allowed = x.allowed, successes = x.successes, runs = x.runs, power = x.power, se = x.se) for x in h.path],
            targets = h.targets.cap ? merge(h.targets, (rho_star = nothing,)) : h.targets, type1 = h.type1)    # PP step 5: in the cap case no rho* exists
end

function _p3_hc_block(res)
    res.hc === nothing && return nothing
    hc = res.hc
    return hc.undefined ? hc : merge(hc, (table = [(nf = r.nf, successes = r.successes, runs = r.runs, power = r.power, se = r.se) for r in hc.table],))
end

"""
    p3_record(res, q, fp, program) -> NamedTuple

The record of E58 as an ordered NamedTuple (the keys of the top level in the order of the file): `program, fingerprints, parameters, outcome, relabelling, sets, theta_space,
selection, c7, optimism, pp_hb, pp_hc, r_ab, readings`.
"""
function p3_record(res, q, fp, program)
    outcome = res.outcome === :selected ? "selected" : res.outcome === :no_feasible_rule ? "no feasible rule" : "sets too small"
    sel = res.selection
    return (program = program, fingerprints = fp, parameters = NamedTuple{P3_PARAM_KEYS}(Tuple(getfield(q, k) for k in P3_PARAM_KEYS)),
            outcome = (selection = outcome, optimism_skipped = res.optimism === nothing ? nothing : res.optimism.skipped, pp_hb = res.hb === nothing ? "not run" : string(res.hb.targets.outcome),
                       pp_hc = res.hc === nothing ? "not run" : res.hc.undefined ? "AUC not defined" : res.hc.underpowered ? "underpowered" : "target reached"),
            relabelling = _p3_relabel_block(res), sets = _p3_sets_block(res),
            theta_space = (size = length(sel.thetas), alone = collect(sel.alone), first_fail = [count(==(j), sel.first_fail) for j in 0:6], cost_cap = sel.cost_cap, sizes = sel.sizes), selection = _p3_selection_block(res, q),
            c7 = res.c7, optimism = res.optimism, pp_hb = _p3_hb_block(res), pp_hc = _p3_hc_block(res), r_ab = res.rab, readings = P3_READINGS)
end
