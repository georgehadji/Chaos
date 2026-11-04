# Study p3-abstention (Amendment D95): the machine file `p3-abstention-test-results.json` (E73). UTF-8, LF, fixed key order, integers, rationals as "a/b" in lowest terms, reals with 12 significant
# digits, no clock and no absolute path. Written with the JSON writer of D93.

const P3T_READINGS = [
    "P87, P88: the commands count-round and analyze; the test-side table, the round log and the tables are in the stage folders of the study with the slug prefix, outside data/",
    "P89: the test-side table has 18 keys with quote checks against the text of Amendment A1 and cross checks against the calibration record (E60)",
    "P91: the counts of the stopping rule are over the orbits with K > 0 of the rounds so far; drift orbits are removed and counted, unreadable orbits are counted as excluded",
    "P92, P93: one entry per call; both commands verify the whole log by recomputing every count from the ref datasets",
    "P94: a decision close at r = r_max regardless of the counts",
    "P95: five strata in the order T_s, W2 minus T_s, S2, the other orbits with rho = C, rho = R; an orbit is in the first it qualifies for",
    "P96, P97: exact type 7 percentiles at p = 1 - alpha_ci and alpha_ci; more than 5 percent dropped replicates give no bound",
    "P98: one stream per replicate, the strata in a fixed order inside the replicate",
    "P99, P100: a bound that is undefined, or a gate orbit without a verdict, makes HA not decidable",
    "P101: orbits with a ref label and no verdict data are dropped from E_test and every set, counted, and the difference to the log counts is reported",
    "P102: |A| of the final set is the denominator of FR_rule, FR_X and the coverage; an undecided comparator verdict is not false-regular",
    "P103: the power at the achieved sizes is recomputed from the rates of the calibration record with new streams and is never used for a decision",
    "P104: the denominator of s_P is |V*|; a D orbit without Double64 verdicts of V* leaves the set; the bootstrap does not stratify by design cell",
    "P105: the REVERSED outcome of HC uses the upper end at p = 1 - alpha_ci",
    "P106: not decidable first (floors, undefined bound, capped comparator, missing gate), then the decision; the flags attach to any decided outcome",
    "P108: the sensitivities are round-down budgets, pilot costs, costs measured in the test run, and without LE; each is reported only",
    "P109: the ladder bound of T3 is read from the P1 parameter table; T3 reads the L5 value at the recorded step N_B (the last recorded step not above it)",
    "P110: the data of F2 and F3 are written with digests; no figure is drawn",
    "P112, P113: T1 transcribes the by-source rows of the record; the round of an orbit is the panel attribute round",
    "P114: the measured costs (E78): tg, fq1, bk1 from the specs of the arm ind, pa1 from the specs of the arm ind_pa1, L1 as 1e9 times the thread seconds over the steps of the ok L1 chunk entries of the arm ref; the manifest status is complete, or complete with failed chunks when a chunk entry failed (E76); a spec without a value makes that sensitivity not decidable",
    "P115: the six arms ref, ind, ind_pa1, d64, d64_pa1, chk are read from the manifest; items are found by content and attributed to their arm by their chunk; tg is read at the horizon N_c of the arms ind and d64 (E75, E81); the coordinate mode of an item is read from the item"]

const P3T_DEVIATIONS = [
    "P120: the column labels agree of table T3 is omitted, a deviation of the analysis: the items of the arm chk run at their own horizon N_B with the label BigFloat(256) and have no steps at the scale of N_1, so the labels cannot be formed in three precisions",
    "P95: the strata of HA are not disjoint by definition (an orbit can be in T_s and S2); the first stratum an orbit qualifies for holds it",
    "P100: a gate orbit without verdict items makes HA not decidable and is listed",
    "P104: D orbits without Double64 verdicts of V* leave the HC set"]

_p3t_pq(x) = x === nothing ? nothing : x

# The comparisons of HB as rows for the machine file and T6, in the order of the sets: the voters of I_6 on T_s, then I_7 on W2.
function p3t_cmp_rows(hb, M, qt, bx)
    out = NamedTuple[]
    for (A, Xs, name) in ((:Ts, Symbol.(qt.voter_bits), "T_s"), (:W2, Symbol.(qt.indicators_i7), "W2")), X in Xs
        haskey(hb.rows, (A, X)) || continue
        r = hb.rows[(A, X)]
        b = bx(X)
        cost = b === nothing ? nothing : b * _p3_cost(qt, _p3_spec(qt, X))
        push!(out, (set = name, indicator = string(X), b = b, cost = cost, capped = X === :LE ? false : M.X[X].capped, cost_unmatched = X === :LE ? M.le_unmatched : false,
                    fr_x = r.fr_x, fr_rule = r.fr_rule, n11 = r.n11, n10 = r.n10, n01 = r.n01, n00 = r.n00, p = r.p, margin = r.margin, pass = r.pass))
    end
    return out
end

_p3t_fail_names(f) = String[(k[1] === :Ts ? "T_s" : "W2") * ":" * string(k[2]) for k in f]

function _p3t_bounds_json(b)
    return NamedTuple{keys(b)}(Tuple((value = v.value, dropped = v.dropped) for v in values(b)))
end

function _p3t_ha_json(A)
    h = A.ha
    a = h.point.a
    return (e = h.point.e, a = (E_test = a[:E], W2 = a[:W2], T_s = a[:Ts], S2 = a[:S2]), dw = h.point.dw, dt = h.point.dt, g = h.point.g, bounds = _p3t_bounds_json(h.bounds), strata = h.strata,
            gate_missing = A.O.gate.missing_ids, verdict = h.decision, failed = h.failed)
end

function _p3t_hb_json(A)
    r, hb = A.hbr, A.hbr.hb
    cmp = p3t_cmp_rows(hb, A.M, A.qt, r.bx)
    wo = p3t_cmp_rows(r.without_le, A.M, A.qt, r.bx)
    return (verdict = hb.decision, comparisons = cmp, failing = _p3t_fail_names(hb.failing), coverage = (T_s = r.cov[:Ts], W2 = r.cov[:W2]), coverage_ok = hb.coverage_ok,
            coverage_price = (T_s = 1 - r.cov[:Ts], W2 = 1 - r.cov[:W2]), underpowered = hb.underpowered, cost_unmatched_le = true, defects = string.(hb.defects),
            sizes = (T_s = r.sizes.Ts, W2 = r.sizes.W2), without_le = (verdict = r.without_le.decision, comparisons = wo, failing = _p3t_fail_names(r.without_le.failing)))
end

function _p3t_hc_json(A)
    h = A.hc
    return (verdict = h.decision, auc_i = h.auc_i, auc_p = h.auc_p, delta = h.delta, lower = h.lower, upper = h.upper, failures = h.failures, nonfailures = h.nonfailures,
            failure_ids = Int[r.id for r in h.rows if r.f], dropped = h.dropped, underpowered = h.underpowered, set_size = length(h.rows), target_failures = A.c.t.target_f)
end

_p3t_sens_json(s) = (decidable = s.decidable, verdict = s.verdict, failing = _p3t_fail_names(s.failing), reason = s.reason)

function _p3t_sets_block(A)
    O, E = A.O, A.O.rows
    ids(f) = sort!(Int[r.id for r in E if f(r)])
    counts(f) = count(f, E)
    sets = (E_test = ids(_ -> true), W2 = ids(r -> r.W2), T_s = ids(r -> r.Ts), S2 = ids(r -> r.S2), D = sort!(collect(keys(O.vd))), D_B = A.dB, gate = O.gate.ids)
    cnt = (E_test = length(E), W2 = counts(r -> r.W2), T_s = counts(r -> r.Ts), S2 = counts(r -> r.S2), D = length(sets.D), D_B = length(A.dB), gate = length(O.gate.ids), overlap = counts(r -> r.Ts && r.S2))
    kr = sort!(unique((r.kindex, r.round) for r in E))
    byk = [(k_index = k, round = rd, E_test = count(r -> r.kindex == k && r.round == rd, E), W2 = count(r -> r.kindex == k && r.round == rd && r.W2, E),
            T_s = count(r -> r.kindex == k && r.round == rd && r.Ts, E), S2 = count(r -> r.kindex == k && r.round == rd && r.S2, E)) for (k, rd) in kr]
    last = A.entries[end].counts
    return (counts = cnt, ids = sets, strata = A.ha.strata, drops = O.drops, by_k_round = byk, log_difference = (t_s = last["t_s"] - cnt.T_s, w2 = last["w2"] - cnt.W2))
end

function _p3t_outputs(A)
    cnt(rs) = (chaotic_supported = count(r -> r.out === :chaotic_supported, rs), regular_supported = count(r -> r.out === :regular_supported, rs), unresolved = count(r -> r.out === :unresolved, rs))
    E = A.O.rows
    return (E_test = cnt(E), W2 = cnt([r for r in E if r.W2]), T_s = cnt([r for r in E if r.Ts]), S2 = cnt([r for r in E if r.S2]), D = cnt([r for r in E if haskey(A.O.vd, r.id)]))
end

function _p3t_fingerprints(A, files)
    c = A.c
    h(f) = f === nothing ? nothing : get(c.raw, f, nothing)
    dir = pkgdir(parentmodule(@__MODULE__))::String
    g = git_state(dir)
    return (preregistration = c.prereg_hash, p3_parameters = h(files.params), test_parameters = h(files.test_params), calibration_record = h(files.calibration_record), panel = c.panel_hash,
            manifest = h(files.manifest), round_log = h(files.round_log), p2_manifest = h(files.p2_manifest), p2_parameters = h(files.p2_params), p2_preregistration = h(files.p2_prereg),
            p1_parameters = h(files.p1_params), p1_preregistration = h(files.p1_prereg), git_commit = g.sha, git_clean = _p3_git_clean(dir),
            chunks = Dict(x.file => c.raw[x.path] for x in c.chunks))
end

"The machine file of E73 as an ordered record (the text is `p3_json`)."
function p3t_results(A, files)
    hbj = _p3t_hb_json(A)
    dec = (HA = (verdict = A.ha.decision, failed = A.ha.failed), HB = (verdict = A.hb.decision, failed = hbj.failing), HC = (verdict = A.hc.decision, failed = String[]))
    program = (name = "Chaoslab.Analysis.analyze (study p3-abstention, test side)", generator = P3_PROGRAM.generator, quantile = P3_PROGRAM.quantile, module_version = string(package_version()))
    flags = (underpowered_hb = A.hb.underpowered, underpowered_hc = A.hc.underpowered, cost_unmatched_le = true,
             not_decidable = [h for (h, v) in (("HA", A.ha.decision), ("HB", A.hb.decision), ("HC", A.hc.decision)) if v == "not decidable"])
    return (program = program, fingerprints = _p3t_fingerprints(A, files), parameters = NamedTuple{P3_PARAM_KEYS}(Tuple(getfield(A.c.q, k) for k in P3_PARAM_KEYS)),
            test_parameters = A.c.t, log = (entries = [(round = e.round, decision = e.decision, counts = e.counts) for e in A.entries],), sets = _p3t_sets_block(A), outputs = _p3t_outputs(A),
            decisions = dec, ha = _p3t_ha_json(A), hb = hbj, hc = _p3t_hc_json(A),
            sensitivities = (round_down = _p3t_sens_json(A.sens.round_down), pilot = _p3t_sens_json(A.sens.pilot), measured = _p3t_sens_json(A.sens.measured), without_le = _p3t_sens_json(A.sens.without_le)),
            power = (hb = A.power.hb, hc = A.power.hc), flags = flags, r_ab = (pairs = A.rab.pairs, overrides = A.rab.overrides, orbits = A.rab.orbits),
            readings = vcat(P3T_READINGS, P3T_DEVIATIONS), exploratory = String[])
end
