# Study p3-abstention (Amendment D95): the test-side table of E60 (18 keys: the five pilot costs and the values of Amendment A1) and its checks against the preregistration, the
# calibration record and the 71-key table. Same file form as D93 E47 (`value`, `section`, `quote`, `text`). The module holds none of the values (I41).

const P3T_PILOT_KEYS = (:pilot_cost_l1, :pilot_cost_tg, :pilot_cost_pa1, :pilot_cost_bk1, :pilot_cost_fq1)
const P3T_PARAM_KEYS = (P3T_PILOT_KEYS..., :calibration_record_sha256, :p3_table_sha256, :target_t_s, :target_w2, :target_f, :theta_voters, :theta_q_c, :theta_q_r, :theta_w, :theta_v,
                        :power_hb_planned, :power_hc_planned, :cap_flag)
const P3T_INT_KEYS = (:target_t_s, :target_w2, :target_f, :theta_q_c, :theta_q_r, :theta_w, :theta_v)
const P3T_HASH_KEYS = (:calibration_record_sha256, :p3_table_sha256)

_p3terr(path, k, msg) = throw(ArgumentError("test-side parameter file $path: key $k $msg"))

function _p3t_value(path, k::Symbol, v)
    if k in P3T_INT_KEYS
        _isint(v) || _p3terr(path, k, "needs an integer, got $(repr(v))")
        return Int(v)
    elseif k in P3T_HASH_KEYS
        (v isa AbstractString && occursin(r"^[0-9a-f]{64}$", v)) || _p3terr(path, k, "needs a text of 64 lowercase hex digits, got $(repr(v))")
        return String(v)
    elseif k === :theta_voters
        (v isa AbstractVector && all(x -> x isa AbstractString, v)) || _p3terr(path, k, "needs a list of texts, got $(repr(v))")
        return String[x for x in v]
    elseif k === :cap_flag
        v isa Bool || _p3terr(path, k, "needs true or false, got $(repr(v))")
        return v
    end
    (v isa Real && !(v isa Bool)) || _p3terr(path, k, "needs a number, got $(repr(v))")
    return Float64(v)
end

# Missing and unknown keys are named; a repeated key is named by the parser (`_p2_parse`).
function _p3t_keys(path, d)
    unknown = sort!([k for k in keys(d) if !(Symbol(k) in P3T_PARAM_KEYS)])
    missing_ = [string(k) for k in P3T_PARAM_KEYS if !haskey(d, string(k))]
    (isempty(unknown) && isempty(missing_)) && return
    throw(ArgumentError("test-side parameter file $path: " * join(vcat(isempty(unknown) ? String[] : ["unknown key $(join(unknown, ", "))"],
                                                                      isempty(missing_) ? String[] : ["missing key $(join(missing_, ", "))"]), "; ")))
end

function _p3t_entry(path, k::Symbol, e)
    e isa AbstractDict || _p3terr(path, k, "must be a table with value, section, quote and text")
    haskey(e, "value") || _p3terr(path, k, "has no value")
    for f in ("section", "quote", "text")
        (haskey(e, f) && e[f] isa AbstractString && !isempty(e[f])) || _p3terr(path, k, "has no $f (a text)")
    end
    return (value = _p3t_value(path, k, e["value"]), section = String(e["section"]), qtext = String(e["quote"]), text = String(e["text"]))
end

# The quote in its section of the preregistration, the text in its quote and (a number, a hash, a list, a flag) equal to the value.
function _p3t_quotes(path, p, entries, prereg_text)
    secs = _p2_sections(prereg_text)
    for k in P3T_PARAM_KEYS
        e = entries[k]
        haskey(secs, e.section) || _p3terr(path, k, "names the section \"$(e.section)\", which the preregistration does not have")
        occursin(replace(e.qtext, "\r\n" => "\n"), secs[e.section]) || _p3terr(path, k, "has a quote that does not occur in the section \"$(e.section)\"")
        occursin(e.text, e.qtext) || _p3terr(path, k, "has the text \"$(e.text)\", which does not occur in its quote")
        v = getfield(p, k)
        if v isa Real && !(v isa Bool)
            _p2_same_number(path, k, e.text, v)
        else
            want = v isa AbstractVector ? join(v, ", ") : string(v)
            e.text == want || _p3terr(path, k, "has the text \"$(e.text)\", which differs from the value $(repr(v))")
        end
    end
end

function _p3t_ranges(path, p)
    for k in P3T_PILOT_KEYS
        (isfinite(getfield(p, k)) && getfield(p, k) > 0) || _p3terr(path, k, "must be positive, got $(getfield(p, k))")
    end
    for k in (:target_t_s, :target_w2, :target_f, :theta_q_c, :theta_q_r)
        getfield(p, k) >= 1 || _p3terr(path, k, "must be at least 1, got $(getfield(p, k))")
    end
    for k in (:theta_w, :theta_v)
        getfield(p, k) >= 0 || _p3terr(path, k, "must be at least 0, got $(getfield(p, k))")
    end
    for k in (:power_hb_planned, :power_hc_planned)
        0 <= getfield(p, k) <= 1 || _p3terr(path, k, "must lie in the closed unit interval, got $(getfield(p, k))")
    end
    isempty(p.theta_voters) && _p3terr(path, :theta_voters, "must not be empty")
end

# The raw `value` fields of the 71-key table that the checks of E60 need (the full table is read and checked by `p3_read_params`).
function _p3t_table71(params)
    d = _p2_parse(params)
    val(k) = haskey(d, k) && d[k] isa AbstractDict && haskey(d[k], "value") ? d[k]["value"] : throw(ArgumentError("71-key table $params: key $k is missing"))
    return (voter_bits = String.(val("voter_bits")), min_voters = val("min_voters"), w_range = Int.(val("w_range")), v_range = Int.(val("v_range")), budgets = Int.(val("budgets")),
            b_primary = val("b_primary"), n_t_min = val("n_t_min"), n_w_min = val("n_w_min"), y_t = val("y_t"), y_w = val("y_w"), r_max = val("r_max"), test_seed = val("test_seed"))
end

# Checks (a) and (b) of E60: the hashes of the record and of the 71-key table, and the outcome of the record.
function _p3t_check_record_hashes(path, p, record, params, rec)
    bytes2hex(sha256(read(record))) == p.calibration_record_sha256 || _p3terr(path, :calibration_record_sha256, "is $(p.calibration_record_sha256), but the SHA-256 of $record differs")
    h71 = bytes2hex(sha256(read(params)))
    p.p3_table_sha256 == h71 || _p3terr(path, :p3_table_sha256, "is $(p.p3_table_sha256), but the SHA-256 of the 71-key table $params is $h71")
    get(get(rec, "fingerprints", Dict()), "p3_parameters", "") == h71 || _p3terr(path, :p3_table_sha256, "differs from the value fingerprints.p3_parameters written in the calibration record")
    out = get(get(rec, "outcome", Dict()), "selection", "absent")
    out == "selected" || throw(ArgumentError("calibration record $record: the outcome is \"$out\", not \"selected\": the study stopped at calibration"))
end

# A planned power equals the row of the record exactly: the rational successes / runs or the decimal text of its power (E60 c); a missing row or null counts refuse naming the key.
function _p3t_same_power(path, k, row, planned, what)
    (row isa AbstractDict && get(row, "successes", nothing) !== nothing && get(row, "runs", nothing) !== nothing) ||
        _p3terr(path, k, "cannot be checked: the calibration record has no $what with successes and runs")
    want = p3_rat(planned)
    ok = want == p3_rat(row["successes"]) // p3_rat(row["runs"]) || (get(row, "power", nothing) !== nothing && want == p3_rat(row["power"]))
    ok || _p3terr(path, k, "is $planned, but the calibration record has $(row["power"]) ($(row["successes"]) of $(row["runs"]))")
end

# Check (c): theta*, the targets, the cap flag and the planned powers of the record equal the keys (the HB power is that of the path point j of the targets, p3_record's `pp_hb`).
function _p3t_check_record_values(path, p, rec)
    th = rec["selection"]["theta"]
    same(k, a, b) = a == b || _p3terr(path, k, "is $(repr(b)), but the calibration record has $(repr(a))")
    same(:theta_voters, String.(th["V"]), p.theta_voters)
    same(:theta_q_c, th["qC"], p.theta_q_c); same(:theta_q_r, th["qR"], p.theta_q_r); same(:theta_w, th["w"], p.theta_w); same(:theta_v, th["v"], p.theta_v)
    (rec["pp_hb"] === nothing || rec["pp_hc"] === nothing) && _p3terr(path, :target_f, "cannot be checked: the calibration record has no block pp_hb or pp_hc")
    tg, hc = rec["pp_hb"]["targets"], rec["pp_hc"]
    same(:target_t_s, tg["n_T"], p.target_t_s); same(:target_w2, tg["n_W"], p.target_w2); same(:cap_flag, tg["cap"], p.cap_flag)
    hc["target"] === nothing && _p3terr(path, :target_f, "cannot be checked: pp_hc.target of the calibration record is null")
    same(:target_f, hc["target"], p.target_f)
    pt = findfirst(x -> x["j"] == tg["j"], rec["pp_hb"]["path"])
    _p3t_same_power(path, :power_hb_planned, pt === nothing ? nothing : rec["pp_hb"]["path"][pt], p.power_hb_planned, "path point j = $(tg["j"])")
    _p3t_same_power(path, :power_hc_planned, hc["power_expected"], p.power_hc_planned, "row power_expected")
end

# Checks (d), (e), (f): theta* in its ranges, the targets between the floors and the caps (exact), the seed and `r_max` of the record.
function _p3t_check_ranges71(path, p, t, rec)
    all(X -> X in t.voter_bits, p.theta_voters) && length(unique(p.theta_voters)) >= t.min_voters || _p3terr(path, :theta_voters, "must be a subset of I_6 ($(join(t.voter_bits, ", "))) of at least $(t.min_voters) elements")
    n = length(p.theta_voters)
    1 <= p.theta_q_c <= n || _p3terr(path, :theta_q_c, "must lie in 1 to $n")
    (n ÷ 2 + 1) <= p.theta_q_r <= n || _p3terr(path, :theta_q_r, "must lie in $(n ÷ 2 + 1) to $n")
    p.theta_w in t.w_range || _p3terr(path, :theta_w, "is not in w_range $(t.w_range)")
    p.theta_v in t.v_range || _p3terr(path, :theta_v, "is not in v_range $(t.v_range)")
    (t.b_primary % 2^p.theta_w == 0 && (t.b_primary ÷ 2^p.theta_w) in t.budgets) || _p3terr(path, :theta_w, "gives b_primary / 2^w = $(t.b_primary) / $(2^p.theta_w), not an element of B")
    rm = p3_rat(t.r_max)
    p.target_t_s >= t.n_t_min || _p3terr(path, :target_t_s, "is below n_t_min = $(t.n_t_min)")
    p.target_w2 >= t.n_w_min || _p3terr(path, :target_w2, "is below n_w_min = $(t.n_w_min)")
    p.target_t_s <= t.y_t * rm || _p3terr(path, :target_t_s, "is above y_t x r_max = $(t.y_t * rm)")
    p.target_w2 <= t.y_w * rm || _p3terr(path, :target_w2, "is above y_w x r_max = $(t.y_w * rm)")
    rec["parameters"]["test_seed"] == t.test_seed || throw(ArgumentError("test_seed of the calibration record is $(rec["parameters"]["test_seed"]), but the 71-key table has $(t.test_seed)"))
    p3_rat(rec["parameters"]["r_max"]) == rm || throw(ArgumentError("r_max of the calibration record is $(rec["parameters"]["r_max"]), but the 71-key table has $(t.r_max)"))
end

"""
    p3t_read_params(path; prereg, record, params) -> NamedTuple

Read the test-side table of study `p3-abstention` (E60, 18 keys: `P3T_PARAM_KEYS`) and apply its checks: types, the quote of every key in its section of `prereg` (the five pilot costs in
the cost section, the values of Amendment A1 in the heading of that amendment), the text against the quote and the value, ranges, and the checks (a) to (f) against the calibration `record`
and the 71-key table `params`. An `ArgumentError` names the key (or the outcome of the record).
"""
function p3t_read_params(path::AbstractString; prereg, record, params)
    isfile(path) || throw(ArgumentError("test-side parameter file $path does not exist"))
    d = _p2_parse(path)
    _p3t_keys(path, d)
    entries = Dict(k => _p3t_entry(path, k, d[string(k)]) for k in P3T_PARAM_KEYS)
    p = NamedTuple{P3T_PARAM_KEYS}(Tuple(entries[k].value for k in P3T_PARAM_KEYS))
    _p3t_quotes(path, p, entries, read(prereg, String))
    _p3t_ranges(path, p)
    isfile(record) || throw(ArgumentError("calibration record $record does not exist"))
    _p3t_check_record_hashes(path, p, record, params, read_json(record))
    rec = read_json(record)
    _p3t_check_record_values(path, p, rec)
    _p3t_check_ranges71(path, p, _p3t_table71(params), rec)
    return p
end
