# Study p3-abstention (Amendment D95): the pipeline of `analyze` for the test side (Algorithm of the amendment, steps 1 to 8). Nothing is written before every check has passed (E61); the
# decisions are made from the numbers of the tables only, in the order HA, HB, HC; no step uses a result of a later step.

# Steps 1 to 3: the checks, the verification of the closed log, the manifest and the items of the four arms, read once.
function _p3t_inputs(a, on_read)
    p3t_check_not_data(a.out, "--out")
    extra = String[f for f in (a.manifest, a.round_log, a.p1_prereg, a.p1_params) if f !== nothing]
    c = p3t_context(; prereg = a.prereg, params = a.params, test_params = a.test_params, store = a.store, panel = a.panel, p2_manifest = a.p2_manifest, p2_prereg = a.p2_prereg, p2_params = a.p2_params,
                    calibration_record = a.calibration_record, on_read, command = "analyze", files = extra, p1_params = a.p1_params,
                    after_hash = () -> p3t_read_closed_log(a.round_log))
    text = p3t_read_closed_log(a.round_log)
    k = length(p3t_split_log(text))
    st = _p3t_extend_state(nothing, c, 1, k, on_read)
    v = p3t_verify_log(text, c, st)
    v.closed || _fail("the round log $(a.round_log) is not closed: its last entry says continue, analyze needs a log whose last entry says close")
    _p3t_no_round_above(c, k, "after the last entry (round $k) of the closed round log")
    isfile(a.manifest) || _fail("manifest $(a.manifest) does not exist")
    man = read_json(a.manifest)
    p3t_check_manifest(man, c, a.store)
    p3t_check_closure_complete(c)
    items, idx = p3t_load_items(c, st, on_read)
    p3t_check_identity(items, c.panel_orbits)
    return (c = c, st = st, v = v, man = man, items = items, idx = idx)
end

# Steps 4 to 6: sets, verdicts, outputs and the three hypotheses with their bootstrap.
function _p3t_hypotheses(I, threads)
    c = I.c
    qt = merge(c.q, c.t)
    theta = p3t_theta(c.t)
    kP = _p3_kP(qt)
    O = p3t_orbits(c, I.st, I.idx, theta)
    M = p3_matched(qt, theta.V)
    ha = p3t_ha([(id = r.id, rho = r.rho, out = r.out, W2 = r.W2, Ts = r.Ts, S2 = r.S2) for r in O.rows], O.gate.out, qt; threads)
    hbr = p3t_hb_run(O.rows, qt, M, theta)
    rab = p3t_rab(O.rows, O.vd, theta, kP)
    hc = p3t_hc(p3t_hc_rows(O.rows, O.vd, theta, kP), qt; threads)
    return (qt = qt, theta = theta, kP = kP, O = O, M = M, ha = ha, hbr = hbr, hb = hbr.hb, rab = rab, hc = hc)
end

# Step 7: the power at the achieved sizes (only when a target is missed), the sensitivities, the figure data.
function _p3t_reports(I, H, rec, threads)
    qt = H.qt
    sizes = H.hbr.sizes
    phb = H.hb.underpowered ? p3t_power_hb(p3t_record_model(rec, qt), sizes.Ts, sizes.W2 - sizes.Ts, qt, threads) : nothing
    ok = H.hc.failures > 0 && H.hc.nonfailures > 0
    phc = H.hc.underpowered && ok ? p3t_power_hc(rec, H.hc.failures, H.hc.nonfailures, qt, threads) : nothing
    sens = p3t_sensitivities(H.O.rows, qt, H.theta, I.man, H.hbr)
    return (power = (hb = phb, hc = phc), sens = sens, f2 = p3t_f2(H.O.rows, qt, H.theta, H.M), f3 = p3t_f3(H.hc))
end

"""
    analyze_p3t(; prereg, manifest, panel, params, test_params, store, round_log, p2_manifest, p2_prereg, p2_params, calibration_record, p1_prereg, p1_params, out, threads, on_read) -> NamedTuple

The test-side analysis of study `p3-abstention` (Amendment D95): checks every input (an `ArgumentError` names the cause; nothing is written and the output folder is not created), verifies the closed
round log, labels the test orbits, forms the sets, decides HA, HB and HC by the frozen rules, and writes `<out>/p3-abstention-test-results.json` and `<out>/p3-abstention-tables.md`. Called
through [`analyze`](@ref) for the preregistration of this study. Returns `(status, results, tables, f2, f3, decisions)`.
"""
function analyze_p3t(; prereg, manifest, panel, params, store, out, threads::Integer = Threads.nthreads(), p1_prereg = nothing, p1_params = nothing, test_params, round_log, p2_manifest, p2_prereg,
                     p2_params, calibration_record, on_read = nothing)
    a = (prereg = prereg, manifest = manifest, panel = panel, params = params, store = store, out = out, p1_prereg = p1_prereg, p1_params = p1_params, test_params = test_params,
         round_log = round_log, p2_manifest = p2_manifest, p2_prereg = p2_prereg, p2_params = p2_params, calibration_record = calibration_record)
    I = _p3t_inputs(a, on_read)
    H = _p3t_hypotheses(I, threads)
    rec = read_json(calibration_record)
    R = _p3t_reports(I, H, rec, threads)
    dB = p3t_dB(I.c, I.idx)
    decisions = (HA = H.ha.decision, HB = H.hb.decision, HC = H.hc.decision)
    A = (c = I.c, qt = H.qt, theta = H.theta, st = I.st, entries = I.v.entries, man = I.man, items = I.items, idx = I.idx, O = H.O, dB = dB, M = H.M, ha = H.ha, hbr = H.hbr, hb = H.hb, hc = H.hc,
         rab = H.rab, sens = R.sens, power = R.power, f2 = R.f2, f3 = R.f3, decisions = decisions, byid = Dict(r.id => r for r in H.O.rows))
    files = (params = params, test_params = test_params, calibration_record = calibration_record, manifest = manifest, round_log = round_log, p2_manifest = p2_manifest, p2_params = p2_params,
             p2_prereg = p2_prereg, p1_params = p1_params, p1_prereg = p1_prereg)
    text = p3_json(p3t_results(A, files))
    tables = p3t_tables(A, files, bytes2hex(sha256(codeunits(text))), rec)
    p3t_recheck(I.c.raw)
    mkpath(out)
    rp, tp = joinpath(out, P3_SLUG * "-test-results.json"), joinpath(out, P3_SLUG * "-tables.md")
    write(rp, text); write(tp, tables)
    return (status = :ok, results = rp, tables = tp, f2 = R.f2, f3 = R.f3, decisions = decisions)
end
