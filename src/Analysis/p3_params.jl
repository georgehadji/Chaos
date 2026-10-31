# Study p3-abstention (Amendment D93): the parameter table of E47 (71 keys). One TOML table per key with `value`, `section`, `quote` and `text` (no machine key).
# The module owns no scientific threshold (I33): every number comes from this table, which transcribes the preregistration; the verdict rules come from the p2 table.

const P3_SLUG = "p3-abstention"

const P3_PARAM_KEYS = (:budgets, :b_primary, :N_c, :N_1, :n_s, :f_s, :lam_w2, :lam_s2, :min_set, :gate_size,
    :indicators_i7, :voter_bits, :w_range, :v_range, :min_voters,
    :cost_l1, :cost_tg, :cost_pa1, :cost_bk1, :cost_fq1, :ratio_d64_tg, :ratio_d64_fq1, :ratio_d64_pa1, :ratio_d64_bk1, :spec_of,
    :e_cal, :e_max, :a_cal, :a_max, :s_cal, :delta_a, :alpha_ci, :g_max, :c_min, :alpha_hb, :rho_hb, :n_t_min, :n_w_min, :p_hb, :p_hc, :delta_c, :delta_tol,
    :b_opt, :rule_rate_floor, :rule_rate_factor, :auc_shrink, :auc_p, :nonfail_per_fail, :pp_runs, :test_seed, :r_max, :rho_path_start, :rho_path_second, :y_t, :y_w,
    :nf_start, :nf_second, :ceil_step, :z_one_sided, :shared_default, :d_expected_total, :d_worst_total,
    :screening_points_per_axis, :d_frac_high, :d_frac_low, :d_b_size, :safety_float, :safety_dd, :wall_cap_hours, :bootstrap_resamples,
    :p2_manifest_sha256)

const P3_INT_KEYS = (:b_primary, :N_c, :N_1, :n_s, :min_set, :gate_size, :min_voters, :n_t_min, :n_w_min, :b_opt, :nonfail_per_fail, :pp_runs, :test_seed, :y_t, :y_w,
                     :nf_start, :nf_second, :ceil_step, :d_expected_total, :d_worst_total, :screening_points_per_axis, :d_b_size, :bootstrap_resamples, :g_max)
const P3_INTVEC_KEYS = (:budgets, :w_range, :v_range)
const P3_STRVEC_KEYS = (:indicators_i7, :voter_bits)
const P3_SPECS = ("l1", "tg", "pa1", "bk1", "fq1")
# rates and shares in the closed unit interval; the AUC of the comparator in the open one; the rest of the reals strictly positive
const P3_UNIT_KEYS = (:f_s, :e_cal, :e_max, :a_cal, :a_max, :s_cal, :delta_a, :alpha_ci, :c_min, :alpha_hb, :rho_hb, :p_hb, :p_hc, :delta_c, :delta_tol, :rule_rate_floor,
                      :shared_default, :auc_shrink, :d_frac_high, :d_frac_low)
const P3_POSITIVE_KEYS = (:lam_w2, :lam_s2, :cost_l1, :cost_tg, :cost_pa1, :cost_bk1, :cost_fq1, :ratio_d64_tg, :ratio_d64_fq1, :ratio_d64_pa1, :ratio_d64_bk1, :rule_rate_factor,
                          :r_max, :rho_path_start, :rho_path_second, :z_one_sided, :safety_float, :safety_dd, :wall_cap_hours)

_p3err(path, k, msg) = throw(ArgumentError("parameter file $path: key $k $msg"))

"""
    p3_rat(x) -> Rational{BigInt}

The exact rational of the decimal text of a real (`207.6` is `2076 // 10`, not the nearest binary fraction); an integer is exact (I32).
"""
p3_rat(x::Integer) = Rational{BigInt}(x)
function p3_rat(x::Real)
    s = string(Float64(x))
    m = match(r"^(-?)(\d+)(?:\.(\d+))?(?:e([+-]?\d+))?$", s)
    m === nothing && throw(ArgumentError("p3_rat: cannot read the number $s"))
    frac = m.captures[3] === nothing ? "" : String(m.captures[3])
    num = parse(BigInt, String(m.captures[2]) * frac) // big(10)^length(frac)
    m.captures[1] == "-" && (num = -num)
    e = m.captures[4] === nothing ? 0 : parse(Int, String(m.captures[4]))
    return Rational{BigInt}(e >= 0 ? num * big(10)^e : num // big(10)^(-e))
end

function _p3_value(path, k::Symbol, v)
    if k in P3_INT_KEYS
        _isint(v) || _p3err(path, k, "needs an integer, got $(repr(v))")
        return Int(v)
    elseif k in P3_INTVEC_KEYS
        (v isa AbstractVector && all(_isint, v)) || _p3err(path, k, "needs a list of integers, got $(repr(v))")
        return Int[x for x in v]
    elseif k in P3_STRVEC_KEYS
        (v isa AbstractVector && all(x -> x isa AbstractString, v)) || _p3err(path, k, "needs a list of texts, got $(repr(v))")
        return String[x for x in v]
    elseif k === :spec_of
        (v isa AbstractDict && all(kv -> kv[1] isa AbstractString && kv[2] isa AbstractString, v)) || _p3err(path, k, "needs a table indicator = spec of texts, got $(repr(v))")
        return Dict{String, String}(String(a) => String(b) for (a, b) in v)
    elseif k === :p2_manifest_sha256
        v isa AbstractString || _p3err(path, k, "needs a text, got $(repr(v))")
        return String(v)
    end
    (v isa Real && !(v isa Bool)) || _p3err(path, k, "needs a number, got $(repr(v))")
    return Float64(v)
end

function _p3_entry(path, k::Symbol, e)
    e isa AbstractDict || _p3err(path, k, "must be a table with value, section, quote and text")
    haskey(e, "value") || _p3err(path, k, "has no value")
    for f in ("section", "quote", "text")
        (haskey(e, f) && e[f] isa AbstractString && !isempty(e[f])) || _p3err(path, k, "has no $f (a text)")
    end
    return (value = _p3_value(path, k, e["value"]), section = String(e["section"]), qtext = String(e["quote"]), text = String(e["text"]))
end

# Missing, unknown and repeated keys are named together.
function _p3_keys(path, d)
    unknown = sort!([k for k in keys(d) if !(Symbol(k) in P3_PARAM_KEYS)])
    missing_ = [string(k) for k in P3_PARAM_KEYS if !haskey(d, string(k))]
    (isempty(unknown) && isempty(missing_)) && return
    parts = String[]
    isempty(unknown) || push!(parts, "unknown key $(join(unknown, ", "))")
    isempty(missing_) || push!(parts, "missing key $(join(missing_, ", "))")
    throw(ArgumentError("parameter file $path: " * join(parts, "; ")))
end

# Quote in its section, text in its quote, a numeric text equal to its value (E47 check i, as D89 A2).
function _p3_quotes(path, p, entries, text)
    secs = _p2_sections(text)
    for k in P3_PARAM_KEYS
        e = entries[k]
        haskey(secs, e.section) || _p3err(path, k, "names the section \"$(e.section)\", which the preregistration does not have")
        occursin(replace(e.qtext, "\r\n" => "\n"), secs[e.section]) || _p3err(path, k, "has a quote that does not occur in the section \"$(e.section)\"")
        occursin(e.text, e.qtext) || _p3err(path, k, "has the text \"$(e.text)\", which does not occur in its quote")
        _p2_same_number(path, k, e.text, getfield(p, k))
    end
end

function _p3_ranges(path, p)
    for k in P3_PARAM_KEYS
        v = getfield(p, k)
        v isa Float64 && !isfinite(v) && _p3err(path, k, "is not finite")
        k in P3_UNIT_KEYS && !(0 <= v <= 1) && _p3err(path, k, "must lie in the closed unit interval, got $v")
        k in P3_POSITIVE_KEYS && !(v > 0) && _p3err(path, k, "must be positive, got $v")
        k in P3_INT_KEYS && k !== :test_seed && k !== :g_max && v < 1 && _p3err(path, k, "must be at least 1, got $v")
        (k in P3_INTVEC_KEYS || k in P3_STRVEC_KEYS) && isempty(v) && _p3err(path, k, "must not be empty")
        k in (:w_range, :v_range) && any(<(0), v) && _p3err(path, k, "must hold integers of at least 0")
    end
    0 < p.auc_p < 1 || _p3err(path, :auc_p, "must lie in the open unit interval, got $(p.auc_p)")
    p.g_max >= 0 || _p3err(path, :g_max, "must be at least 0")
    p.rho_path_second > p.rho_path_start || _p3err(path, :rho_path_second, "must exceed rho_path_start (the step of the rho_path is positive)")
    p.nf_second > p.nf_start || _p3err(path, :nf_second, "must exceed nf_start (the step of the n_f search is positive)")
    p.y_w >= p.y_t || _p3err(path, :y_w, "must be at least y_t")
    p.d_worst_total >= p.d_expected_total || _p3err(path, :d_worst_total, "must be at least d_expected_total")
end

# B has ten elements, b_primary is one of them, N_c the last; the voters are the indicators without the first (E47 v).
function _p3_consistency(path, p)
    length(p.budgets) == 10 || _p3err(path, :budgets, "has $(length(p.budgets)) entries, B has ten")
    (all(>(0), p.budgets) && issorted(p.budgets; lt = <=)) || _p3err(path, :budgets, "must be positive and strictly increasing")
    p.b_primary in p.budgets || _p3err(path, :b_primary, "= $(p.b_primary) is not an element of budgets")
    p.N_c == p.budgets[end] || _p3err(path, :N_c, "is $(p.N_c), not the last budget $(p.budgets[end])")
    p.voter_bits == p.indicators_i7[2:end] || _p3err(path, :voter_bits, "must equal indicators_i7 without its first element, got $(repr(p.voter_bits))")
    p.min_voters <= length(p.voter_bits) || _p3err(path, :min_voters, "is above the number of voters")
    for X in p.indicators_i7
        s = get(p.spec_of, X, nothing)
        s in P3_SPECS || _p3err(path, :spec_of, "has no valid spec ($(join(P3_SPECS, ", "))) for the indicator $X")
    end
end

# The mirror keys equal the p2 keys of the same names (E47 ii); the p2 numbers are compared as numbers.
function _p3_mirrors(path, p, p2)
    for (k, k2) in ((:budgets, :budgets), (:b_primary, :b_primary), (:N_c, :N_c), (:N_1, :N_1), (:n_s, :n_s), (:f_s, :f_s), (:min_set, :min_set), (:gate_size, :gate_size),
                    (:lam_w2, :w2_lambda_max), (:lam_s2, :s2_lambda_min))
        getfield(p, k) == getfield(p2, k2) || _p3err(path, k, "is $(repr(getfield(p, k))), but $k2 of the p2 parameter table is $(repr(getfield(p2, k2))) (a p2 mirror)")
    end
end

# The p2 parameter table, validated against the p2 preregistration (E47 iv); also returned for the verdict rules of table R1 and the label numbers.
_p3_p2_table(p2_params, p2_prereg, p1_params, defaults) = p2_read_params(p2_params; prereg = p2_prereg, p1_params = p1_params, defaults = defaults)

function _p3_read(path, prereg, p2_params, p2_prereg, p2_manifest, p1_params, defaults)
    isfile(path) || throw(ArgumentError("parameter file $path does not exist"))
    d = _p2_parse(path)
    _p3_keys(path, d)
    entries = Dict(k => _p3_entry(path, k, d[string(k)]) for k in P3_PARAM_KEYS)
    q = NamedTuple{P3_PARAM_KEYS}(Tuple(entries[k].value for k in P3_PARAM_KEYS))
    _p3_quotes(path, q, entries, read(prereg, String))
    _p3_ranges(path, q)
    _p3_consistency(path, q)
    p2 = _p3_p2_table(p2_params, p2_prereg, p1_params, defaults)
    _p3_mirrors(path, q, p2)
    isfile(p2_manifest) || throw(ArgumentError("p2 manifest $p2_manifest does not exist"))
    h = bytes2hex(sha256(read(p2_manifest)))
    q.p2_manifest_sha256 == h || _p3err(path, :p2_manifest_sha256, "is $(q.p2_manifest_sha256), but the SHA-256 of the p2 manifest is $h")
    return (q = q, p2 = p2)
end

"""
    p3_read_params(path; prereg, p2_params, p2_prereg, p2_manifest, p1_params = nothing, defaults = nothing) -> NamedTuple

Read the parameter table of study `p3-abstention` (E47, 71 keys) and apply its checks: types, the quote of every key against its section of `prereg`, the text against
its quote and (a number) against the value, the value ranges, `voter_bits` against `indicators_i7`, the mirrors against the p2 table (which is itself checked against
`p2_prereg`, D89 E22) and `p2_manifest_sha256` against the p2 manifest file. An `ArgumentError` names the key.
"""
p3_read_params(path::AbstractString; prereg, p2_params, p2_prereg, p2_manifest, p1_params = nothing, defaults = nothing) =
    _p3_read(path, prereg, p2_params, p2_prereg, p2_manifest, p1_params, defaults).q
