# Indicator specs of the run block (Amendment D86, R11): the repeated table `[[run.indicator]]` and the arm key `indicators`. Validation runs in two
# passes inside `run_block`, before any computation: each spec alone (`_rb_indicator`), then the specs against the arms (`_ind_cross`). The module
# holds no default for any key of a spec (I27) and no key that could hold a threshold or a class.

const RB_KINDS = ("tangent", "pair", "stretching", "frequency", "phase_advance", "birkhoff")
const IND_WINDOW_KINDS = ("frequency", "phase_advance", "birkhoff")
const IND_TANGENT_MEMBERS = ("SALI", "FLI", "MEGNO")
const IND_PAIR_MEMBERS = ("RLI", "MEGNO")
const IND_MAX_OBSERVABLES = 31

const _IND_KEYS = Dict{String, Tuple}("tangent" => ("members",), "pair" => ("members", "separation", "frame_source"),
    "stretching" => ("range", "bins", "interval", "min_per_cell", "counts"),
    "frequency" => ("observable", "windows", "order", "positions", "pad", "mean_removal", "mirror_ratio", "floor_multiplier"),
    "phase_advance" => ("observable", "windows", "birkhoff_order"), "birkhoff" => ("observables", "windows", "birkhoff_order"))

"True for the kinds read at their own evaluation windows (frequency, phase advance, Birkhoff)."
ind_window_kind(kind::AbstractString) = kind in IND_WINDOW_KINDS

# ---------------------------------------------------------------- one spec alone
function _ind_members(v, path, allowed)
    m = _rb_list(v, path, (x, p) -> _rb_str(x, p, allowed))
    allunique(m) || _rb_err(path, "a member may appear once, got $(join(m, ", "))")
    "SALI" in m && last(m) != "SALI" && _rb_err(path, "SALI replaces the columns and may only be the last member (D71), got $(join(m, ", "))")
    return m
end

function _ind_windows(v, path, minw)
    w = _rb_list(v, path, (x, p) -> _rb_int(x, p; min = minw))
    issorted(w; lt = <=) || _rb_err(path, "must be strictly increasing, got $(join(w, ", "))")
    return w
end

function _ind_order(v, path; lo, hi)
    n = _rb_int(v, path; min = lo)
    n <= hi || _rb_err(path, "must be at most $hi, got $n")
    return n
end

function _ind_range(v, path)
    v == "derived" && return "derived"
    (v isa AbstractVector && length(v) == 2) || _rb_err(path, "must be \"derived\" or a list of two numbers lo < hi, got $(repr(v))")
    lo, hi = _rb_real(v[1], "$path[1]"), _rb_real(v[2], "$path[2]")
    lo < hi || _rb_err(path, "needs lo < hi, got ($lo, $hi)")
    return [lo, hi]
end

function _ind_observable_tables(v, path)
    (v isa AbstractVector && !isempty(v)) || _rb_err(path, "must be a non-empty list of tables")
    length(v) <= IND_MAX_OBSERVABLES || _rb_err(path, "holds at most $IND_MAX_OBSERVABLES observables, got $(length(v))")
    return map(enumerate(v)) do (i, o)
        p = "$path[$i]"
        o isa AbstractDict || _rb_err(p, "must be a table")
        kind = _rb_str(get(o, "kind", nothing), "$p.kind", ("cos", "sin", "lifted"))
        _rb_keys(o, p, kind == "lifted" ? ("kind", "coordinate") : ("kind", "coordinate", "frequency"))
        out = Dict{String, Any}("kind" => kind, "coordinate" => _ind_order(o["coordinate"], "$p.coordinate"; lo = 1, hi = 2))
        kind == "lifted" || (out["frequency"] = _rb_int(o["frequency"], "$p.frequency"))
        return out
    end
end

function _ind_kind_keys(s, kind, path)
    if kind == "tangent"
        return Dict{String, Any}("members" => _ind_members(s["members"], "$path.members", IND_TANGENT_MEMBERS))
    elseif kind == "pair"
        return Dict{String, Any}("members" => _ind_members(s["members"], "$path.members", IND_PAIR_MEMBERS),
            "separation" => _rb_real(s["separation"], "$path.separation"; ok = >(0), what = "positive"),
            "frame_source" => _rb_str(s["frame_source"], "$path.frame_source", ("seeded", "direction")))
    elseif kind == "stretching"
        return Dict{String, Any}("range" => _ind_range(s["range"], "$path.range"), "bins" => _rb_int(s["bins"], "$path.bins"),
            "interval" => _rb_int(s["interval"], "$path.interval"), "min_per_cell" => _rb_int(s["min_per_cell"], "$path.min_per_cell"; min = 0),
            "counts" => _rb_str(s["counts"], "$path.counts", ("none", RB_SUBSETS...)))
    elseif kind == "frequency"
        order = _ind_order(s["order"], "$path.order"; lo = 0, hi = 5)
        s["mean_removal"] isa Bool || _rb_err("$path.mean_removal", "must be true or false, got $(repr(s["mean_removal"]))")
        return Dict{String, Any}("observable" => _rb_str(s["observable"], "$path.observable", ("exp_i_q", "q_plus_i_p")),
            "windows" => _ind_windows(s["windows"], "$path.windows", 8 * (order + 1)), "order" => order,
            "positions" => _rb_str(s["positions"], "$path.positions", ("half", "edge")), "pad" => _rb_int(s["pad"], "$path.pad"; min = 2),
            "mean_removal" => s["mean_removal"], "mirror_ratio" => _rb_real(s["mirror_ratio"], "$path.mirror_ratio"; ok = x -> 0 < x <= 1, what = "inside (0, 1]"),
            "floor_multiplier" => _rb_real(s["floor_multiplier"], "$path.floor_multiplier"; ok = >(0), what = "positive"))
    elseif kind == "phase_advance"
        return Dict{String, Any}("observable" => _rb_str(s["observable"], "$path.observable", ("planar", "lifted")),
            "windows" => _ind_windows(s["windows"], "$path.windows", 16), "birkhoff_order" => _ind_order(s["birkhoff_order"], "$path.birkhoff_order"; lo = 1, hi = 4))
    end
    return Dict{String, Any}("observables" => _ind_observable_tables(s["observables"], "$path.observables"),
        "windows" => _ind_windows(s["windows"], "$path.windows", 16), "birkhoff_order" => _ind_order(s["birkhoff_order"], "$path.birkhoff_order"; lo = 1, hi = 4))
end

# One `[[run.indicator]]` table, without its budget (checked against the precisions of the arms in the second pass).
function _rb_indicator(s, i)
    path = "indicator[$i]"
    s isa AbstractDict || _rb_err(path, "must be a table")
    haskey(s, "id") || _rb_err("$path.id", "missing key")
    id = s["id"]
    (id isa AbstractString && occursin(r"^[a-z][a-z0-9_]*$", id)) ||
        _rb_err("$path.id", "must be lowercase letters, digits and underscore, starting with a letter, got $(repr(id))")
    lowercase(id) in map(lowercase, RB_METHODS) && _rb_err("$path.id", "`$id` is the id of a catalog method")
    haskey(s, "kind") || _rb_err("$path.kind", "missing key")
    kind = _rb_str(s["kind"], "$path.kind", RB_KINDS)
    _rb_keys(s, path, ("id", "kind", "budget", _IND_KEYS[kind]...))
    s["budget"] isa AbstractDict || _rb_err("$path.budget", "must be a table")
    return merge(Dict{String, Any}("id" => String(id), "kind" => kind), _ind_kind_keys(s, kind, path))
end

# ---------------------------------------------------------------- specs against arms
function _ind_budget(b, precisions, path)
    out = Dict{String, Any}()
    for k in keys(b)
        k in RB_PRECISIONS || _rb_err("$path.$k", "unknown precision")
    end
    for p in precisions
        haskey(b, p) || _rb_err("$path.$p", "missing key for a precision of the arms that use the spec")
        t = b[p]
        _rb_keys(t, "$path.$p", ("ns_cpu_thread", "safety"), ("ns_gpu",))
        e = Dict{String, Any}("ns_cpu_thread" => _rb_real(t["ns_cpu_thread"], "$path.$p.ns_cpu_thread"; ok = >(0), what = "positive"),
                              "safety" => _rb_real(t["safety"], "$path.$p.safety"; ok = >=(1), what = "at least 1"))
        haskey(t, "ns_gpu") && (e["ns_gpu"] = _rb_real(t["ns_gpu"], "$path.$p.ns_gpu"; ok = >(0), what = "positive"))
        out[p] = e
    end
    return out
end

"Run length (handoffs) of a window spec (E22); `nothing` for the kinds read on the recording grid."
function ind_run_length(s::AbstractDict)
    s["kind"] in ("frequency", "birkhoff") && return 2 * maximum(s["windows"])
    s["kind"] == "phase_advance" && return 2 * maximum(s["windows"]) + 1
    return nothing
end

# Mode rules of E26 and the interval rule, per arm that uses the spec.
function _ind_arm_rules(s, i, arm, ai, prim)
    path, taus, modes = "indicator[$i]", _rb_taus(arm, prim), _arm_dim(arm, "coordinate_mode", prim)
    s["kind"] == "stretching" && for tau in taus
        s["interval"] % tau == 0 || _rb_err("$path.interval", "interval $(s["interval"]) is not a multiple of tau = $tau of arm[$ai]")
    end
    "lifted" in modes || return
    s["kind"] == "frequency" && s["observable"] == "exp_i_q" &&
        _rb_err("$path.observable", "exp_i_q needs the wrapped coordinate mode, arm[$ai] uses lifted")
    s["kind"] == "birkhoff" && any(o -> o["kind"] in ("cos", "sin"), s["observables"]) &&
        _rb_err("$path.observables", "cos and sin need the wrapped coordinate mode, arm[$ai] uses lifted")
end

# A window kind fixes the horizon of every arm that uses it: tau times the run length, one tau per arm (R11 item 3).
function _ind_window_horizons(s, arm, ai, prim)
    L = ind_run_length(s)
    L === nothing && return
    taus = _rb_taus(arm, prim)
    length(taus) == 1 || _rb_err("arm[$ai].tau", "an arm that uses the window spec `$(s["id"])` has a single tau, got $(join(taus, ", "))")
    want = taus[1] * L
    all(==(want), arm["horizons"]) || _rb_err("arm[$ai].horizons", "an arm that uses the window spec `$(s["id"])` has the horizon tau x run length = $want, got $(join(arm["horizons"], ", "))")
end

function _ind_cross(raw_specs, arms, prim, tables)
    ids = [s["id"] for s in raw_specs]
    for (i, id) in enumerate(ids)
        id in ids[1:(i - 1)] && _rb_err("indicator[$i].id", "duplicate id `$id`")
    end
    users = Dict(id => Int[] for id in ids)
    for (ai, a) in enumerate(arms), id in get(a, "indicators", String[])
        haskey(users, id) || _rb_err("arm[$ai].indicators", "unknown indicator id `$id`")
        push!(users[id], ai)
        "BigFloat" in a["precisions"] && _rb_err("arm[$ai].precisions", "BigFloat is not offered for indicator specs (Float64 and Double64 only)")
    end
    return map(enumerate(raw_specs)) do (i, s)
        isempty(users[s["id"]]) && _rb_err("indicator[$i]", "spec `$(s["id"])` is used by no arm")
        for ai in users[s["id"]]
            _ind_arm_rules(s, i, arms[ai], ai, prim)
            _ind_window_horizons(s, arms[ai], ai, prim)
        end
        precs = [p for p in RB_PRECISIONS if any(ai -> p in arms[ai]["precisions"], users[s["id"]])]
        return merge(s, Dict{String, Any}("budget" => _ind_budget(tables[i]["budget"], precs, "indicator[$i].budget")))
    end
end

function _rb_indicators(v, arms, prim)
    if v === nothing
        any(a -> haskey(a, "indicators"), arms) && _rb_err("indicator", "arms list indicator ids but the block has no [[run.indicator]] table")
        return Dict{String, Any}[]
    end
    (v isa AbstractVector && !isempty(v)) || _rb_err("indicator", "must be a non-empty list of tables")
    return _ind_cross([_rb_indicator(s, i) for (i, s) in enumerate(v)], arms, prim, v)
end

# ---------------------------------------------------------------- canonical text and hash of a spec (E27)
"Spec `id` of a block (`ArgumentError` if absent)."
function block_spec(rb::RunBlock, id::AbstractString)
    i = findfirst(s -> s["id"] == id, get(rb.raw, "indicator", Dict{String, Any}[]))
    i === nothing && throw(ArgumentError("run block has no indicator spec `$id`"))
    return rb.raw["indicator"][i]
end

"Canonical text of a spec without its `budget`: `indicator.<key> = <value>`, keys sorted, lists in order (E27)."
function spec_text(s::AbstractDict)
    io = IOBuffer()
    _canon!(io, Dict{String, Any}(k => v for (k, v) in s if k != "budget"), "indicator")
    return String(take!(io))
end
spec_text(rb::RunBlock, id::AbstractString) = spec_text(block_spec(rb, id))

"SHA-256 (64 hexadecimal digits) of [`spec_text`](@ref)."
spec_hash(s::AbstractDict) = bytes2hex(sha256(spec_text(s)))
spec_hash(rb::RunBlock, id::AbstractString) = spec_hash(block_spec(rb, id))
