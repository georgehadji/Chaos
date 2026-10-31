# The machine-readable run block of the preregistration (Amendment A4): validation before any computation (D53a) and its
# canonical text. Every study value is read from the block; the module holds none (D53r).

const RB_METHODS = ("L1", "L2", "L3", "L4", "L4b", "L5", "L5b", "L6", "L9")
const RB_PRECISIONS = ("Float64", "Double64", "BigFloat")
const RB_SUBSETS = ("P", "S_D", "S_L", "h4", "D", "D_B")
const RB_RELEASES = ("each_round", "closure")      # D94: when the round log releases an arm
const RB_MODES = ("lifted", "wrapped")
const RB_REQUESTS = ("auto", "gpu", "cpu")
const RB_CLASSES = ("T", "R", "C", "S", "X", "control_elliptic", "control_hyperbolic", "h4")

"A validated run block; `raw` is the normalized tree (integers as `Int`, reals as `Float64`, directions as text)."
struct RunBlock
    raw::Dict{String, Any}
end

_rb_err(key, msg) = throw(ArgumentError("run block `$key`: $msg"))

function _rb_keys(d, path, required, optional = ())
    d isa AbstractDict || _rb_err(path, "must be a table")
    for k in keys(d)
        (k in required || k in optional) || _rb_err("$path.$k", "unknown key")
    end
    for k in required
        haskey(d, k) || _rb_err("$path.$k", "missing key")
    end
end

function _rb_int(v, path; min = 1)
    (v isa Integer && !(v isa Bool)) || _rb_err(path, "must be an integer, got $(repr(v))")
    v >= min || _rb_err(path, "must be at least $min, got $v")
    return Int(v)
end

function _rb_real(v, path; ok = x -> true, what = "a finite number")
    (v isa Real && !(v isa Bool) && isfinite(v)) || _rb_err(path, "must be a finite number, got $(repr(v))")
    ok(v) || _rb_err(path, "must be $what, got $v")
    return Float64(v)
end

function _rb_str(v, path, allowed)
    (v isa AbstractString && String(v) in allowed) || _rb_err(path, "must be one of $(join(allowed, ", ")), got $(repr(v))")
    return String(v)
end

function _rb_list(v, path, f)
    (v isa AbstractVector && !isempty(v)) || _rb_err(path, "must be a non-empty list")
    return [f(x, "$path[$i]") for (i, x) in enumerate(v)]
end

# A direction is a closed-form text evaluated at 256 bits; a number is taken as its shortest decimal text.
function _rb_direction(v, path)
    s = v isa Real && !(v isa Bool) ? string(v) : v
    s isa AbstractString || _rb_err(path, "must be a closed-form text, got $(repr(v))")
    try
        setprecision(() -> _eval_exact(String(s)), BigFloat, 256)
    catch e
        e isa ArgumentError ? _rb_err(path, "cannot evaluate \"$s\"") : rethrow()
    end
    return String(s)
end

function _rb_primary(p)
    _rb_keys(p, "primary", ("tau", "eps", "direction", "coordinate_mode"))
    return Dict{String, Any}("tau" => _rb_int(p["tau"], "primary.tau"),
        "eps" => _rb_real(p["eps"], "primary.eps"; ok = >(0), what = "positive"),
        "direction" => _rb_direction(p["direction"], "primary.direction"),
        "coordinate_mode" => _rb_str(p["coordinate_mode"], "primary.coordinate_mode", RB_MODES))
end

const _ARM_REQUIRED = ("name", "orbits", "precisions", "horizons")
const _ARM_OPTIONAL = ("methods", "indicators", "tau", "eps", "direction", "coordinate_mode", "h1_check", "stop_exempt", "release")
const RB_CONTROL_SCOPES = ("all", "exact_start")

function _rb_arm(a, path)
    _rb_keys(a, path, _ARM_REQUIRED, _ARM_OPTIONAL)
    name = a["name"]
    (name isa AbstractString && !isempty(name)) || _rb_err("$path.name", "must be a non-empty text")
    out = Dict{String, Any}("name" => String(name), "orbits" => _rb_str(a["orbits"], "$path.orbits", RB_SUBSETS),
        "precisions" => _rb_list(a["precisions"], "$path.precisions", (x, p) -> _rb_str(x, p, RB_PRECISIONS)),
        "horizons" => _rb_list(a["horizons"], "$path.horizons", (x, p) -> _rb_int(x, p)))
    _rb_arm_lists!(out, a, path)
    haskey(a, "tau") && (out["tau"] = _rb_list(a["tau"], "$path.tau", (x, p) -> _rb_int(x, p)))
    haskey(a, "eps") && (out["eps"] = _rb_list(a["eps"], "$path.eps", (x, p) -> _rb_real(x, p; ok = >(0), what = "positive")))
    haskey(a, "direction") && (out["direction"] = _rb_list(a["direction"], "$path.direction", _rb_direction))
    haskey(a, "coordinate_mode") &&
        (out["coordinate_mode"] = _rb_list(a["coordinate_mode"], "$path.coordinate_mode", (x, p) -> _rb_str(x, p, RB_MODES)))
    haskey(a, "release") && (out["release"] = _rb_str(a["release"], "$path.release", RB_RELEASES))
    for k in ("h1_check", "stop_exempt")
        haskey(a, k) || continue
        a[k] isa Bool || _rb_err("$path.$k", "must be true or false, got $(repr(a[k]))")
        out[k] = a[k]
    end
    return out
end

# D86: an arm lists catalog methods, indicator spec ids or both; `methods` may be absent or empty when `indicators` is non-empty.
function _rb_arm_lists!(out, a, path)
    has_i = haskey(a, "indicators")
    (haskey(a, "methods") || has_i) || _rb_err(path, "the arm lists neither methods nor indicators")
    has_i && (out["indicators"] = _rb_list(a["indicators"], "$path.indicators", function (x, p)
        (x isa AbstractString && !isempty(x)) || _rb_err(p, "must be an indicator id, got $(repr(x))")
        return String(x)
    end))
    if haskey(a, "methods")
        out["methods"] = has_i && a["methods"] isa AbstractVector && isempty(a["methods"]) ? String[] :
                         _rb_list(a["methods"], "$path.methods", (x, p) -> _rb_str(x, p, RB_METHODS))
    end
    has_i && !allunique(out["indicators"]) && _rb_err("$path.indicators", "an id may appear once")
    return out
end

_rb_taus(arm, prim) = get(arm, "tau", [prim["tau"]])

function _rb_arms(v, prim, has_h4)
    arms = _rb_list(v, "arm", (a, p) -> _rb_arm(a, p))
    names = [a["name"] for a in arms]
    allunique(names) || _rb_err("arm.name", "duplicate arm name $(repr(first(n for n in names if count(==(n), names) > 1)))")
    for (i, a) in enumerate(arms)
        a["orbits"] == "D_B" && get(a, "release", "") == "each_round" &&
            _rb_err("arm[$i].release", "arm `$(a["name"])` has orbits D_B, which does not exist before closure: release each_round is refused")
        a["orbits"] == "h4" && !has_h4 && _rb_err("h4", "arm `$(a["name"])` uses the subset h4 but the block lists no [[run.h4]] start")
        for N in a["horizons"], tau in _rb_taus(a, prim)
            N % tau == 0 || _rb_err("arm[$i].horizons", "horizon $N is not a multiple of tau = $tau")
        end
    end
    return arms
end

function _rb_grid(g, arms, prim)
    _rb_keys(g, "grid", ("points_per_decade", "first_steps", "extra_fractions"))
    fr = _rb_list(g["extra_fractions"], "grid.extra_fractions", (x, p) -> _rb_real(x, p; ok = x -> 0 < x < 1, what = "inside (0, 1)"))
    for (i, a) in enumerate(arms), N in a["horizons"], tau in _rb_taus(a, prim), f in fr
        f == 0.1 && !isempty(get(a, "methods", String[])) && !(N % 10 == 0 && (N ÷ 10) % tau == 0) &&
            _rb_err("grid.extra_fractions", "arm[$i]: the fraction 0.1 needs N / 10 = $(N / 10) to be a multiple of tau = $tau")
    end
    return Dict{String, Any}("points_per_decade" => _rb_int(g["points_per_decade"], "grid.points_per_decade"),
                             "first_steps" => _rb_int(g["first_steps"], "grid.first_steps"), "extra_fractions" => fr)
end

function _rb_budget(b, used)
    b isa AbstractDict || _rb_err("budget", "must be a table")
    out = Dict{String, Any}()
    for k in keys(b)
        k in RB_PRECISIONS || _rb_err("budget.$k", "unknown precision")
    end
    for p in used
        haskey(b, p) || _rb_err("budget.$p", "missing key for a used precision")
        t = b[p]
        _rb_keys(t, "budget.$p", ("ns_cpu_thread", "safety"), ("ns_gpu",))
        e = Dict{String, Any}("ns_cpu_thread" => _rb_real(t["ns_cpu_thread"], "budget.$p.ns_cpu_thread"; ok = >(0), what = "positive"),
                              "safety" => _rb_real(t["safety"], "budget.$p.safety"; ok = >=(1), what = "at least 1"))
        haskey(t, "ns_gpu") && (e["ns_gpu"] = _rb_real(t["ns_gpu"], "budget.$p.ns_gpu"; ok = >(0), what = "positive"))
        out[p] = e
    end
    return out
end

function _rb_backend(b, used)
    b isa AbstractDict || _rb_err("backend", "must be a table")
    out = Dict{String, Any}()
    for k in keys(b)
        k in RB_PRECISIONS || _rb_err("backend.$k", "unknown precision")
    end
    for p in used
        haskey(b, p) || _rb_err("backend.$p", "missing key for a used precision")
        out[p] = _rb_str(b[p], "backend.$p", RB_REQUESTS)
    end
    return out
end

function _rb_health(h, arms, spec_ids = String[])
    keys_ = ("fail_fraction", "control_bound", "control_methods", "setup_check_steps", "setup_check_tol", "drift_bound",
             "rerun_attempts", "stop_failed_fraction", "first_chunk_cost_factor", "wall_cap_hours")
    _rb_keys(h, "health", keys_, ("control_scope", "control_bound_method"))
    pos(k) = _rb_real(h[k], "health.$k"; ok = >(0), what = "positive")
    frac(k) = _rb_real(h[k], "health.$k"; ok = x -> 0 <= x < 1, what = "inside [0, 1)")
    s = h["setup_check_steps"]
    _rb_keys(s, "health.setup_check_steps", RB_CLASSES)
    shortest = minimum(minimum(a["horizons"]) for a in arms)
    steps = Dict{String, Any}(c => _rb_int(s[c], "health.setup_check_steps.$c") for c in RB_CLASSES)
    for (c, n) in steps
        n <= shortest || _rb_err("health.setup_check_steps.$c", "$n steps exceed the shortest horizon $shortest")
    end
    scope = haskey(h, "control_scope") ? Dict{String, Any}("control_scope" => _rb_str(h["control_scope"], "health.control_scope", RB_CONTROL_SCOPES)) :
            Dict{String, Any}()
    cbm = Dict{String, Any}()                      # D57: optional per-method control bound; a method not listed keeps `control_bound`
    if haskey(h, "control_bound_method")
        t = h["control_bound_method"]
        t isa AbstractDict || _rb_err("health.control_bound_method", "must be a table")
        for (m, v) in t
            m in RB_METHODS || _rb_err("health.control_bound_method.$m", "unknown method")
            m in h["control_methods"] || _rb_err("health.control_bound_method.$m", "method is not in control_methods")
            cbm[m] = _rb_real(v, "health.control_bound_method.$m"; ok = >(0), what = "positive")
        end
        scope["control_bound_method"] = cbm
    end
    return merge(scope, Dict{String, Any}("fail_fraction" => frac("fail_fraction"), "control_bound" => pos("control_bound"),
        "control_methods" => _rb_list(h["control_methods"], "health.control_methods", function (x, p)
            x in spec_ids && _rb_err(p, "`$x` is an indicator spec id: the exact-exponent control rule is for the exponent methods")
            return _rb_str(x, p, RB_METHODS)
        end),
        "setup_check_steps" => steps, "setup_check_tol" => pos("setup_check_tol"), "drift_bound" => pos("drift_bound"),
        "rerun_attempts" => _rb_int(h["rerun_attempts"], "health.rerun_attempts"; min = 0),
        "stop_failed_fraction" => frac("stop_failed_fraction"), "first_chunk_cost_factor" => pos("first_chunk_cost_factor"),
        "wall_cap_hours" => pos("wall_cap_hours")))
end

# D86 (R20 item 1): the pilot runs the window kinds at their own run length, so an arm that holds window specs only is outside the pilot horizon rule.
_pilot_arm(a, specs) = !isempty(get(a, "methods", String[])) || any(id -> !ind_window_kind(only(s["kind"] for s in specs if s["id"] == id)), get(a, "indicators", String[]))

function _rb_misc(d, arms, prim, specs = Dict{String, Any}[])
    _rb_keys(d["k19"], "k19", ("threads",))
    v = d["verify"]
    _rb_keys(v, "verify", ("backend_agreement", "rerun_sample"))
    ba = v["backend_agreement"]
    _rb_keys(ba, "verify.backend_agreement", ("per_orbit_sigma", "mean_standard_errors"))
    p = d["pilot"]
    _rb_keys(p, "pilot", ("orbits_per_k", "horizon", "seed_label"))
    ph = _rb_int(p["horizon"], "pilot.horizon")
    parms = [a for a in arms if _pilot_arm(a, specs)]
    shortest = minimum((minimum(a["horizons"]) for a in parms); init = typemax(Int))
    ph <= shortest || _rb_err("pilot.horizon", "$ph exceeds the shortest arm horizon $shortest")
    for a in parms, tau in _rb_taus(a, prim)
        ph % tau == 0 || _rb_err("pilot.horizon", "$ph is not a multiple of tau = $tau")
    end
    (p["seed_label"] isa AbstractString && !isempty(p["seed_label"])) || _rb_err("pilot.seed_label", "must be a non-empty text")
    return Dict{String, Any}("k19" => Dict{String, Any}("threads" => _rb_int(d["k19"]["threads"], "k19.threads")),
        "verify" => Dict{String, Any}("backend_agreement" => Dict{String, Any}(
            "per_orbit_sigma" => _rb_real(ba["per_orbit_sigma"], "verify.backend_agreement.per_orbit_sigma"; ok = >(0), what = "positive"),
            "mean_standard_errors" => _rb_real(ba["mean_standard_errors"], "verify.backend_agreement.mean_standard_errors"; ok = >(0), what = "positive")),
            "rerun_sample" => _rb_int(v["rerun_sample"], "verify.rerun_sample")),
        "pilot" => Dict{String, Any}("orbits_per_k" => _rb_int(p["orbits_per_k"], "pilot.orbits_per_k"), "horizon" => ph,
                                     "seed_label" => String(p["seed_label"])))
end

function _rb_h4(v)
    v === nothing && return Dict{String, Any}[]
    list = _rb_list(v, "h4", function (e, p)
        _rb_keys(e, p, ("K", "x0", "p0", "label"))
        txt(k) = (x = e[k]; x isa Real && !(x isa Bool) ? string(x) : x isa AbstractString ? String(x) : _rb_err("$p.$k", "must be a text"))
        for k in ("K", "x0", "p0")
            _rb_direction(txt(k), "$p.$k")
        end
        return Dict{String, Any}("K" => txt("K"), "x0" => txt("x0"), "p0" => txt("p0"), "label" => txt("label"))
    end)
    labels = [e["label"] for e in list]
    allunique(labels) || _rb_err("h4.label", "duplicate label")
    return list
end

function _rb_accepted(v)
    list = _rb_list(v, "panel_prereg_sha256_accepted", function (x, p)
        (x isa AbstractString && occursin(r"^[0-9a-f]{64}$", x)) || _rb_err("panel_prereg_sha256_accepted", "entry $p must be 64 lowercase hexadecimal digits")
        return String(x)
    end)
    return list
end

"""
    run_block(dict) -> RunBlock

Validate the `[run]` table of Amendment A4 (string keys, as parsed from TOML). An `ArgumentError` names the key, before any
computation.
"""
function run_block(d::AbstractDict)
    top = ("panel_prereg_sha256_accepted", "primary", "arm", "grid", "backend", "health", "budget", "k19", "verify", "pilot")
    _rb_keys(d, "run", top, ("h4", "indicator"))
    h4 = _rb_h4(get(d, "h4", nothing))
    prim = _rb_primary(d["primary"])
    arms = _rb_arms(d["arm"], prim, !isempty(h4))
    used = [p for p in RB_PRECISIONS if any(a -> p in a["precisions"], arms)]
    raw = Dict{String, Any}("panel_prereg_sha256_accepted" => _rb_accepted(d["panel_prereg_sha256_accepted"]), "primary" => prim,
        "arm" => arms, "grid" => _rb_grid(d["grid"], arms, prim), "backend" => _rb_backend(d["backend"], used),
        "budget" => _rb_budget(d["budget"], used), "h4" => h4)
    specs = _rb_indicators(get(d, "indicator", nothing), arms, prim)
    raw["health"] = _rb_health(d["health"], arms, [s["id"] for s in specs])
    isempty(specs) || (raw["indicator"] = specs)
    merge!(raw, _rb_misc(d, arms, prim, specs))
    return RunBlock(raw)
end

# Canonical text: sorted keys, lists in their order, one `path = value` line each.
_cv(x::AbstractString) = repr(String(x))
_cv(x::AbstractVector) = "[" * join((_cv(e) for e in x), ", ") * "]"
_cv(x) = repr(x)

function _canon!(io, x, path)
    if x isa AbstractDict
        foreach(k -> _canon!(io, x[k], path * "." * k), sort!(collect(keys(x))))
    elseif x isa AbstractVector && !isempty(x) && first(x) isa AbstractDict
        foreach(((i, e),) -> _canon!(io, e, "$path[$i]"), enumerate(x))
    else
        println(io, path, " = ", _cv(x))
    end
end

"Canonical text of a parsed block, independent of the textual key order of the file."
function block_text(rb::RunBlock)
    io = IOBuffer()
    _canon!(io, rb.raw, "run")
    return String(take!(io))
end

"SHA-256 (64 hexadecimal digits) of [`block_text`](@ref)."
block_hash(rb::RunBlock) = bytes2hex(sha256(block_text(rb)))

# Merge of the run tables of successive amendments (A5 over A4): a later key replaces an earlier one, tables merge key by key, the
# lists of tables `arm` (by `name`) and `h4` (by `label`) are replaced entry by entry and new entries are appended in order.
const _LIST_KEYS = Dict("arm" => "name", "h4" => "label", "indicator" => "id")

function merge_run(a::AbstractDict, b::AbstractDict)
    out = Dict{String, Any}(a)
    for (k, v) in b
        old = get(out, k, nothing)
        if v isa AbstractDict && old isa AbstractDict
            out[k] = merge_run(old, v)
        elseif haskey(_LIST_KEYS, k) && v isa AbstractVector && old isa AbstractVector
            id = _LIST_KEYS[k]
            list = Any[e for e in old]
            for e in v
                i = findfirst(x -> get(x, id, nothing) == get(e, id, missing), list)
                i === nothing ? push!(list, e) : (list[i] = e)
            end
            out[k] = list
        else
            out[k] = v
        end
    end
    return out
end

"""
    read_run_block(prereg_path) -> RunBlock

Check the preregistration hash (E18), merge the `[run]` tables of all its TOML blocks in file order (Amendment A4, then A5: see
`merge_run`) and validate the result (LF and CRLF files).
"""
function read_run_block(path::AbstractString)
    check_prereg(path)
    merged = nothing
    for m in eachmatch(r"```toml\n(.*?)\n```"s, _lf(read(path, String)))
        d = try
            TOML.parse(string(m.captures[1]))
        catch
            continue
        end
        haskey(d, "run") && (merged = merged === nothing ? d["run"] : merge_run(merged, d["run"]))
    end
    merged === nothing && throw(ArgumentError("$path: no ```toml block with a [run] table (Amendment A4)"))
    return run_block(merged)
end
