# Study p2-indicators (Amendment D89): the parameter table of E22. One TOML table per key with `value`, `section` and either `machine = true` or `quote`.
# The module owns no scientific threshold (I14): every number comes from this table, which transcribes the preregistration.

const P2_PARAM_KEYS = (:budgets, :b_primary, :N_c, :N_1, :N_dir, :n_s, :f_s, :windows, :directions,
    :ratio_threshold, :le_floor, :rli_threshold, :e1_d_threshold, :e2_d_threshold, :wba_digits_threshold,
    :r2_le_factor, :r2_sali_threshold, :r2_fli_factor, :r2_megno_threshold, :persist_divisor, :censor_factor,
    :s2_lambda_min, :w2_lambda_max, :gate_k_index, :gate_classes, :gate_size, :drift_bound,
    :label_lambda_min, :label_growth_ratio, :label_rotation_tol, :label_estimator_bound,
    :h1_u_min, :h1_u_lower_min, :h1_amin_max, :h2_q_min, :h2_q_lower_gt, :h3_r_min, :h3_r_lower_gt, :h3_count_min,
    :h4_dphi_min, :h4_dphi_lower_min, :h5_fmed_min, :h5_fmed_lower_min, :h6_flip_min,
    :bootstrap_resamples, :bootstrap_level, :bootstrap_seed, :min_set, :range_lower_percentile, :range_upper_percentile, :gate_max_share,
    :grid_k_values, :grid_i_count, :grid_j_count, :grid_orbit_count, :panel_sha256, :panel_build_sha256, :plan_sha256)

const P2_INT_KEYS = (:b_primary, :N_c, :N_1, :N_dir, :n_s, :persist_divisor, :censor_factor, :gate_k_index, :gate_size, :h3_count_min, :bootstrap_resamples,
                     :bootstrap_seed, :min_set, :grid_i_count, :grid_j_count, :grid_orbit_count)
const P2_INTVEC_KEYS = (:budgets, :windows)
const P2_STRVEC_KEYS = (:directions, :gate_classes, :grid_k_values)
const P2_STR_KEYS = (:panel_sha256, :panel_build_sha256, :plan_sha256)
const P2_MACHINE_KEYS = (:N_1, :N_c, :N_dir, :windows, :directions, :drift_bound, :grid_k_values, :grid_i_count, :grid_j_count, :grid_orbit_count, :panel_build_sha256)

_p2err(path, k, msg) = throw(ArgumentError("parameter file $path: key $k $msg"))

function _p2_parse(path)
    try
        return TOML.parsefile(path)
    catch e
        text = read(path, String)
        heads = [m.captures[1] for m in eachmatch(r"(?m)^\[([^\[\]\n]+)\][ \t]*$", text)]
        rep = unique([h for h in heads if count(==(h), heads) > 1])
        isempty(rep) || throw(ArgumentError("parameter file $path: repeated key $(join(rep, ", "))"))
        throw(ArgumentError("parameter file $path is not valid TOML: $(first(sprint(showerror, e), 200))"))
    end
end

_isint(x) = x isa Integer && !(x isa Bool)

# The typed value of one key (E22): an `ArgumentError` naming the key for a wrong type.
function _p2_value(path, k::Symbol, v)
    if k in P2_INT_KEYS
        _isint(v) || _p2err(path, k, "needs an integer, got $(repr(v))")
        return Int(v)
    elseif k in P2_INTVEC_KEYS
        (v isa AbstractVector && all(_isint, v)) || _p2err(path, k, "needs a list of integers, got $(repr(v))")
        return Int[x for x in v]
    elseif k in P2_STRVEC_KEYS
        (v isa AbstractVector && all(x -> x isa AbstractString, v)) || _p2err(path, k, "needs a list of texts, got $(repr(v))")
        return String[x for x in v]
    elseif k in P2_STR_KEYS
        v isa AbstractString || _p2err(path, k, "needs a text, got $(repr(v))")
        return String(v)
    end
    (v isa Real && !(v isa Bool)) || _p2err(path, k, "needs a number, got $(repr(v))")
    return Float64(v)
end

# Missing and unknown keys are named together, so that the table of the other study is recognised.
function _p2_keys(path, d)
    unknown = sort!([k for k in keys(d) if !(Symbol(k) in P2_PARAM_KEYS)])
    missing_ = [string(k) for k in P2_PARAM_KEYS if !haskey(d, string(k))]
    (isempty(unknown) && isempty(missing_)) && return
    parts = String[]
    isempty(unknown) || push!(parts, "unknown key $(join(unknown, ", "))")
    isempty(missing_) || push!(parts, "missing key $(join(missing_, ", "))")
    throw(ArgumentError("parameter file $path: " * join(parts, "; ")))
end

function _p2_entry(path, k::Symbol, e)
    e isa AbstractDict || _p2err(path, k, "must be a table with value, section and quote or machine")
    haskey(e, "value") || _p2err(path, k, "has no value")
    (haskey(e, "section") && e["section"] isa AbstractString && !isempty(e["section"])) || _p2err(path, k, "has no section text")
    machine = get(e, "machine", false)
    machine isa Bool || _p2err(path, k, "has a machine field that is not true or false")
    (machine || (haskey(e, "quote") && e["quote"] isa AbstractString && !isempty(e["quote"]))) || _p2err(path, k, "has neither machine = true nor a quote")
    return (value = _p2_value(path, k, e["value"]), raw = e["value"], section = String(e["section"]), machine = machine, qtext = machine ? "" : String(e["quote"]))
end

# Ranges (E22): finite reals, counts at least 1, a level strictly between 0 and 1, lists not empty.
function _p2_ranges(path, p)
    for k in P2_PARAM_KEYS
        v = getfield(p, k)
        v isa Float64 && !isfinite(v) && _p2err(path, k, "is not finite")
        k in P2_INT_KEYS && k !== :bootstrap_seed && k !== :gate_k_index && v < 1 && _p2err(path, k, "must be at least 1, got $v")
        (k in P2_INTVEC_KEYS || k in P2_STRVEC_KEYS) && isempty(v) && _p2err(path, k, "must not be empty")
        k in P2_INTVEC_KEYS && any(<(1), v) && _p2err(path, k, "must hold integers of at least 1")
    end
    p.gate_k_index >= 0 || _p2err(path, :gate_k_index, "must be at least 0")
    0 < p.bootstrap_level < 1 || _p2err(path, :bootstrap_level, "must lie strictly between 0 and 1, got $(p.bootstrap_level)")
    0 <= p.range_lower_percentile < p.range_upper_percentile <= 100 || _p2err(path, :range_lower_percentile, "and range_upper_percentile must satisfy 0 <= lower < upper <= 100")
end

# budgets[k] is twice windows[k], N_c twice the largest window and the last budget, b_primary an element of budgets (E22).
function _p2_budgets(path, p)
    length(p.budgets) == length(p.windows) || _p2err(path, :budgets, "has $(length(p.budgets)) entries, windows has $(length(p.windows))")
    for k in eachindex(p.budgets)
        p.budgets[k] == 2 * p.windows[k] || _p2err(path, :budgets, "entry $k is $(p.budgets[k]), not twice windows[$k] = $(p.windows[k])")
    end
    p.N_c == 2 * maximum(p.windows) || _p2err(path, :N_c, "is $(p.N_c), not twice the largest window $(maximum(p.windows))")
    p.N_c == p.budgets[end] || _p2err(path, :N_c, "is $(p.N_c), not the last budget $(p.budgets[end])")
    p.b_primary in p.budgets || _p2err(path, :b_primary, "= $(p.b_primary) is not an element of budgets")
end

# Section texts of the preregistration by heading, outside the code fences: the lines up to the next heading of the same or a higher level.
function _p2_sections(text::AbstractString)
    lines = split(replace(text, "\r\n" => "\n"), '\n')
    heads = Tuple{Int, Int, String}[]
    fence = false
    for (i, l) in enumerate(lines)
        startswith(l, "```") && (fence = !fence)
        m = fence ? nothing : match(r"^(#+)[ \t]+(.*?)[ \t]*$", l)
        m === nothing || push!(heads, (i, length(m.captures[1]::SubString{String}), String(m.captures[2]::SubString{String})))
    end
    out = Dict{String, String}()
    for (n, (i, lev, title)) in enumerate(heads)
        stop = something(findfirst(h -> h[1] > i && h[2] <= lev, heads[(n + 1):end]), 0)
        last = stop == 0 ? length(lines) : heads[n + stop][1] - 1
        haskey(out, title) || (out[title] = join(lines[(i + 1):last], "\n"))
    end
    return out
end

# The machine values of the run blocks of the preregistration (E22 check i), by key.
function _p2_run_values(rb)
    arms = Dict(a["name"] => a for a in rb.raw["arm"])
    specs = Dict(s["id"] => s for s in get(rb.raw, "indicator", Dict{String, Any}[]))
    for (what, ok) in (("arm lyap_p", haskey(arms, "lyap_p")), ("arm win_p", haskey(arms, "win_p")), ("arm dir_tg", haskey(arms, "dir_tg")), ("indicator fq1", haskey(specs, "fq1")))
        ok || throw(ArgumentError("the run block of the preregistration has no $what"))
    end
    h4 = rb.raw["h4"]
    ij = [match(r"-i(\d+)-j(\d+)$", h["label"]) for h in h4]
    any(isnothing, ij) && throw(ArgumentError("a run block h4 label is not of the form grid-k<K index>-i<i>-j<j>"))
    return (N_1 = only(arms["lyap_p"]["horizons"]), N_c = only(arms["win_p"]["horizons"]), N_dir = only(arms["dir_tg"]["horizons"]), windows = Int.(specs["fq1"]["windows"]),
            directions = String.(arms["dir_tg"]["direction"]), drift_bound = rb.raw["health"]["drift_bound"], grid_k_values = unique(String[h["K"] for h in h4]),
            grid_i_count = length(unique(m.captures[1] for m in ij)), grid_j_count = length(unique(m.captures[2] for m in ij)), grid_orbit_count = length(h4),
            panel_build_sha256 = rb.raw["panel_prereg_sha256_accepted"])
end

function _p2_machine(path, p, rb)
    run = _p2_run_values(rb)
    for k in P2_MACHINE_KEYS
        want, v = getfield(run, k), getfield(p, k)
        ok = k === :panel_build_sha256 ? v in want : isequal(v, want)
        ok || _p2err(path, k, "is $(repr(v)) but the run block of the preregistration has $(repr(want))")
    end
    p.grid_orbit_count == length(p.grid_k_values) * p.grid_i_count * p.grid_j_count ||
        _p2err(path, :grid_orbit_count, "is not the product of the K count, grid_i_count and grid_j_count")
end

# The text that stands for a value in a quote (A2): Julia's `string` of the number, of each element of a list, or the text itself.
_p2_texts(raw, v) = v isa AbstractVector ? [string(x) for x in v] : unique([string(raw), string(v)])

function _p2_quotes(path, p, entries, text)
    secs = _p2_sections(text)
    for k in P2_PARAM_KEYS
        e = entries[k]
        e.machine && continue
        haskey(secs, e.section) || _p2err(path, k, "names the section \"$(e.section)\", which the preregistration does not have")
        occursin(replace(e.qtext, "\r\n" => "\n"), secs[e.section]) || _p2err(path, k, "has a quote that does not occur in the section \"$(e.section)\"")
        v = getfield(p, k)
        ok = v isa AbstractVector ? all(t -> occursin(t, e.qtext), _p2_texts(e.raw, v)) : any(t -> occursin(t, e.qtext), _p2_texts(e.raw, v))
        ok || _p2err(path, k, "has the value $(repr(v)), whose text does not occur in its quote")
    end
end

function _p2_links(path, p, defaults, p1_params)
    seed = defaults === nothing ? default_seed() : Int(TOML.parsefile(defaults)["seed"])
    p.bootstrap_seed == seed || _p2err(path, :bootstrap_seed, "is $(p.bootstrap_seed), but the seed of numerics-defaults.toml is $seed")
    p1_params === nothing && return
    q = read_params(p1_params)
    for (k, k1) in ((:label_lambda_min, :lambda_threshold), (:label_growth_ratio, :growth_ratio), (:label_rotation_tol, :rotation_tol), (:label_estimator_bound, :h1a_bound))
        getfield(p, k) == getfield(q, k1) || _p2err(path, k, "is $(getfield(p, k)), but the key $k1 of the P1 parameter table is $(getfield(q, k1))")
    end
end

"""
    p2_read_params(path; prereg, p1_params = nothing, defaults = nothing) -> NamedTuple

Read the parameter table of study `p2-indicators` (E22, 58 keys) and apply its checks: types and ranges, the consistency of the budgets, the machine keys against
the run blocks of the preregistration `prereg`, the quote of every other key against its section, `bootstrap_seed` against the seed of `defaults` (default: the
repository's `numerics-defaults.toml`) and, with `p1_params`, the four label keys against the P1 table. An `ArgumentError` names the key.
"""
function p2_read_params(path::AbstractString; prereg::AbstractString, p1_params = nothing, defaults = nothing)
    isfile(path) || throw(ArgumentError("parameter file $path does not exist"))
    d = _p2_parse(path)
    _p2_keys(path, d)
    entries = Dict(k => _p2_entry(path, k, d[string(k)]) for k in P2_PARAM_KEYS)
    p = NamedTuple{P2_PARAM_KEYS}(Tuple(entries[k].value for k in P2_PARAM_KEYS))
    _p2_ranges(path, p)
    _p2_budgets(path, p)
    _p2_machine(path, p, read_run_block(prereg))
    _p2_quotes(path, p, entries, read(prereg, String))
    _p2_links(path, p, defaults, p1_params)
    return p
end
