# The run plan (R1): work items, chunks, hashes, step counts. A pure function of (preregistration, panel file, block, defaults,
# code version); nothing is computed here and nothing depends on threads, memory, backend, host or time (I11).

const HASH_VERSION = 1
const PRECISION_ORDER = RB_PRECISIONS

"""
The plan: `items`, `chunks`, `counts (items, chunks)`, `steps`, `hash`, plus the inputs the runner needs (`block`, `defaults`,
`panel`, `study`, `code_version`, `prereg_header`).
"""
struct RunPlan
    items::Vector{NamedTuple}
    chunks::Vector{NamedTuple}
    counts::NamedTuple{(:items, :chunks), Tuple{Int, Int}}
    steps::NamedTuple
    hash::String
    text::String
    block::RunBlock
    defaults::Dict{String, Any}
    panel::NamedTuple
    study::String
    code_version::String
    prereg_header::String
    not_applicable::Dict{String, Vector{Int}}                   # D86: orbits without an item for a spec, by spec id
    omitted_cells::Dict{String, Vector{Tuple{Int, String}}}     # the same with the panel cell of each orbit
end

const _ITEM_FIELDS = (:orbit, :system, :param, :x0, :p0, :run_seed, :orbit_seed, :method, :precision, :horizon, :tau, :eps,
                      :direction, :coordinate_mode, :vectors, :grid_sha256, :code_version, :hash_version, :kind, :spec)

# `eps` for L2 only, `direction` only for the items that have one, `kind` and `spec` only for the indicator kinds (D86, E27).
_item_prop(it, f::Symbol) = hasproperty(it, f) ? getproperty(it, f) : nothing

"LF text of the hashed fields of a work item, one `name=value` line each (`eps` for L2 only; indicator items add `kind` and `spec` last)."
function item_text(it)
    io = IOBuffer()
    for f in _ITEM_FIELDS
        (f in (:eps, :direction, :kind, :spec) && _item_prop(it, f) === nothing) && continue
        print(io, f, "=", getproperty(it, f), "\n")
    end
    return String(take!(io))
end

"SHA-256 (64 hexadecimal digits) of [`item_text`](@ref)."
item_hash(it) = bytes2hex(sha256(item_text(it)))

# ---------------------------------------------------------------- inputs
function _read_panel(path, rb)
    bytes = read(path)
    t = TOML.parse(String(copy(bytes)))
    pre = get(get(t, "manifest", Dict()), "prereg", Dict())
    built = String(get(pre, "computed_sha256", ""))
    accepted = rb.raw["panel_prereg_sha256_accepted"]
    built in accepted || throw(ArgumentError("$path: built under preregistration hash $built, which is not in panel_prereg_sha256_accepted"))
    orbits = sort!(Dict{String, Any}[o for o in t["orbit"]]; by = o -> o["id"])
    return (orbits = orbits, sha256 = bytes2hex(sha256(bytes)), build_hash = built)
end

function _read_defaults(path)
    path = path === nothing ? joinpath(_root(), "_config", "numerics-defaults.toml") : path
    return TOML.parsefile(path)
end

function _study_slug(text)
    m = match(r"(?m)^\| Study slug \| *(.*?) *\|", text)
    m === nothing && throw(ArgumentError("preregistration header has no \"Study slug\" row"))
    return string(m.captures[1])
end

# ---------------------------------------------------------------- orbit bases
# A control without a `stability` key is hyperbolic only if its analytic exponent is positive (M6); anything else is refused.
function _setup_class(o)
    o["control"] || return o["class"]
    st = get(o, "stability") do
        get(o, "exponent", 0.0) > 0 ? "hyperbolic" : throw(ArgumentError("control orbit $(o["id"]) has no stability key and no positive analytic exponent"))
    end
    return st == "hyperbolic" ? "control_hyperbolic" : "control_elliptic"
end

# D56: a control is under the control bound only if its binary64 start is an exact periodic point of the lifted Float64 map: the Y5
# map, or a point whose Float64 iterate after the known period is the start itself, bit for bit (residual exactly 0). An orbit with a
# non-zero winding never is: the lifted state moves by the binary64 value of 2 pi, which is not the period of the map. Decided here,
# from the panel, once per orbit.
function _exact_start(o, param, x0, p0)
    o["control"] || return false
    get(o, "control_kind", "") == "y5" && return true
    sys = StandardMap(exact(Float64, param))
    start = SVector(exact(Float64, x0), exact(Float64, p0))
    x = start
    for _ in 1:Int(get(o, "period", 1))
        x = forward(sys, x)
    end
    return x[1] == start[1] && x[2] == start[2]
end

function _panel_base(o)
    ctrl = o["control"]
    lin = o["system"] == "linear"
    param = lin ? string(o["mu"]) : _decimal_text(String(get(o, "K_text", string(o["K"]))))
    x0, p0 = String(get(o, "x0_dec", string(o["x0"]))), String(get(o, "p0_dec", string(o["p0"])))
    return (id = Int(o["id"]), system = String(o["system"]), param = param, x0 = x0,
            p0 = p0, run_seed = String(o["run_seed"]), orbit_seed = String(o["orbit_seed"]),
            design_class = ctrl ? "control" : String(o["class"]), setup_class = _setup_class(o),
            cell = ctrl ? "controls" : "k$(o["k_index"]):$(o["class"])", exponent = ctrl ? Float64(o["exponent"]) : nothing,
            in_SD = Bool(o["in_SD"]), in_SL = Bool(o["in_SL"]), k_index = Int(get(o, "k_index", 0)), control_exact = _exact_start(o, param, x0, p0))
end

function _h4_base(h, id, seed)
    return (id = id, system = "standard", param = _decimal_text(h["K"]), x0 = _decimal_text(h["x0"]), p0 = _decimal_text(h["p0"]),
            run_seed = string(seed), orbit_seed = string(orbit_seed(seed, id)), design_class = "h4", setup_class = "h4", cell = "h4",
            exponent = nothing, in_SD = false, in_SL = false, k_index = 0, control_exact = false)
end

# ---------------------------------------------------------------- work items
_arm_dim(arm, key, prim) = get(arm, key, [prim[key]])

function _make_item(id, b, method, prec, N, tau, eps, dir, mode, ctx; spec = nothing)
    g = get!(ctx.grids, (N, tau)) do
        gr = record_grid(N, tau; points_per_decade = ctx.grid["points_per_decade"], first_steps = ctx.grid["first_steps"],
                         extra_fractions = ctx.grid["extra_fractions"])
        (grid = gr, hash = grid_hash(gr))
    end
    label = prec == "BigFloat" ? "BigFloat($(ctx.bits))" : prec
    base = (orbit = b.id, system = b.system, param = b.param, x0 = b.x0, p0 = b.p0, run_seed = b.run_seed, orbit_seed = b.orbit_seed,
            method = method, precision = label, horizon = N, tau = tau, eps = eps === nothing ? nothing : string(eps), direction = dir,
            coordinate_mode = mode, vectors = spec === nothing ? (method in ("L1", "L2") ? 1 : 2) : (ind_window_kind(spec["kind"]) ? 1 : 2),
            grid_sha256 = g.hash, code_version = ctx.cv, hash_version = HASH_VERSION,
            kind = spec === nothing ? nothing : spec["kind"], spec = spec === nothing ? nothing : spec_hash(spec))
    return (; id = id, base..., hash = item_hash(base), arms = String[], design_class = b.design_class, setup_class = b.setup_class,
            cell = b.cell, exponent = b.exponent, control_exact = b.control_exact, in_SD = b.in_SD, in_SL = b.in_SL)
end

# Whether a spec item carries the direction of the arm (R11 item 2): the tangent and stretching kinds, the pair with the direction frame source.
ind_has_direction(spec) = spec["kind"] in ("tangent", "stretching") || (spec["kind"] == "pair" && spec["frame_source"] == "direction")

# All (orbit id, method or spec id, precision, horizon, tau, eps, direction, mode) of one arm and one precision, in the order of R1 item 4: the methods,
# then the indicator specs of the arm, over the dimensions that apply (a spec without a direction is not swept over the directions).
function _arm_products(arm, prec, prim, subset, specs = Dict{String, Any}())
    out = Tuple[]
    dirs = _arm_dim(arm, "direction", prim)
    for N in arm["horizons"], mode in _arm_dim(arm, "coordinate_mode", prim), dir in dirs, tau in _arm_dim(arm, "tau", prim)
        for m in get(arm, "methods", String[]), eps in (m == "L2" ? _arm_dim(arm, "eps", prim) : (nothing,)), id in subset
            push!(out, (id, m, prec, N, tau, eps, dir, mode))
        end
        for sid in get(arm, "indicators", String[])
            d = ind_has_direction(specs[sid]) ? dir : nothing
            (d === nothing && dir != first(dirs)) && continue
            for id in subset
                push!(out, (id, sid, prec, N, tau, nothing, d, mode))
            end
        end
    end
    return out
end

function _build_items(rb, bases, ctx)
    raw = rb.raw
    prim = raw["primary"]
    subsets = Dict("P" => [b.id for b in bases if b.design_class != "h4"], "S_D" => [b.id for b in bases if b.in_SD], "S_L" => [b.id for b in bases if b.in_SL],
                   "h4" => [b.id for b in bases if b.design_class == "h4"])
    byid = Dict(b.id => b for b in bases)
    specs = Dict{String, Any}(s["id"] => s for s in get(raw, "indicator", Dict{String, Any}[]))
    items, index = NamedTuple[], Dict{Tuple, Int}()
    before = Dict{String, Tuple{Int, Int}}(a["name"] => (0, 0) for a in raw["arm"])
    pre_steps = Dict{String, Int}(p => 0 for p in PRECISION_ORDER)
    omitted = Dict{String, Vector{Tuple{Int, String}}}()
    checked = Set{Tuple}()
    for prec in PRECISION_ORDER, arm in raw["arm"]
        prec in arm["precisions"] || continue
        for key in _arm_products(arm, prec, prim, subsets[arm["orbits"]], specs)
            spec = get(specs, key[2], nothing)
            b = byid[key[1]]
            if spec !== nothing
                if !ind_applicable(spec, b.system)
                    any(o -> o[1] == b.id, get!(omitted, key[2], Tuple{Int, String}[])) || push!(omitted[key[2]], (b.id, b.cell))
                    continue
                end
                (key[2], b.id, key[4], key[5], key[8]) in checked || (ind_check_plan(spec, b, key[5], key[4], key[8]); push!(checked, (key[2], b.id, key[4], key[5], key[8])))
            end
            n, s = before[arm["name"]]
            before[arm["name"]] = (n + 1, s + key[4])
            pre_steps[prec] += key[4]
            i = get(index, key, 0)
            if i == 0
                push!(items, _make_item(length(items) + 1, b, key[2], prec, key[4], key[5], key[6], key[7], key[8], ctx; spec = spec))
                i = index[key] = length(items)
            end
            push!(items[i].arms, arm["name"])
        end
    end
    return items, before, pre_steps, omitted
end

# ---------------------------------------------------------------- chunks
_label(p) = replace(p, r"\(\d+\)" => "")

_runtype(it) = (method = it.method, precision = it.precision, horizon = it.horizon, tau = it.tau, eps = it.eps,
                direction = it.direction, coordinate_mode = it.coordinate_mode)

_chunk_key(it) = "method=$(it.method) " * (_item_prop(it, :kind) === nothing ? "" : "kind=$(it.kind) spec=$(it.spec[1:16]) ") *
                 "precision=$(it.precision) horizon=$(it.horizon) tau=$(it.tau) eps=$(it.eps === nothing ? "-" : it.eps) " *
                 "direction=$(it.direction === nothing ? "-" : it.direction) mode=$(it.coordinate_mode) cell=$(it.cell)"

# The values that decide a chunk's status (R3 item 2): they change the hash, not the key.
function _status_text(rb, df)
    h = rb.raw["health"]
    io = IOBuffer()
    for k in ("orthogonality_tol", "pairing_residual_tol", "det_residual_tol")
        println(io, k, "=", df[k])
    end
    for k in ("fail_fraction", "control_bound", "drift_bound", "setup_check_tol")
        println(io, k, "=", repr(h[k]))
    end
    println(io, "control_methods=", join(h["control_methods"], ","))
    haskey(h, "control_bound_method") && println(io, "control_bound_method=", join(("$m:$(repr(h["control_bound_method"][m]))" for m in sort(collect(keys(h["control_bound_method"])))), ","))
    println(io, "setup_check_steps=", join(("$c:$(h["setup_check_steps"][c])" for c in sort(collect(keys(h["setup_check_steps"])))), ","))
    return String(take!(io))
end

# The values a chunk takes from the panel besides the item hashes (M6): design class, setup class and the exponent of a control.
_panel_text(it) = "item=$(it.hash) design=$(it.design_class) setup=$(it.setup_class) exponent=$(it.exponent === nothing ? "-" : repr(it.exponent)) exact=$(it.control_exact)"

function _chunk_text(key, its, status, h1)
    return "chaoslab-chunk 2\nkey=" * key * "\n" * join(sort!([_panel_text(it) for it in its]), "\n") * "\nh1_check=$(h1)\n" * status
end

function _build_chunks(items, rb, df)
    status = _status_text(rb, df)
    h1arms = Set(a["name"] for a in rb.raw["arm"] if get(a, "h1_check", false))
    armorder = Dict(a["name"] => i for (i, a) in enumerate(rb.raw["arm"]))
    index, groups = Dict{String, Int}(), Vector{Int}[]
    for it in items
        i = get!(index, _chunk_key(it)) do
            push!(groups, Int[])
            length(groups)
        end
        push!(groups[i], it.id)
    end
    return map(enumerate(groups)) do (ordinal, ids)
        its = items[ids]
        key = _chunk_key(first(its))
        arms = sort!(unique(vcat((it.arms for it in its)...)); by = a -> armorder[a])
        h1 = any(it -> any(in(h1arms), it.arms), its)
        (ordinal = ordinal, key = key, hash = bytes2hex(sha256(_chunk_text(key, its, status, h1))),
         cell = first(its).cell, run_type = _runtype(first(its)), items = ids, orbits = [it.orbit for it in its], arms = arms)
    end
end

"LF text whose SHA-256 is the hash of chunk `ordinal`: its key, the sorted item hashes, the values that decide its status."
function chunk_text(plan::RunPlan, ordinal::Integer)
    c = plan.chunks[ordinal]
    h1arms = Set(a["name"] for a in plan.block.raw["arm"] if get(a, "h1_check", false))
    h1 = any(i -> any(in(h1arms), plan.items[i].arms), c.items)
    return _chunk_text(c.key, plan.items[c.items], _status_text(plan.block, plan.defaults), h1)
end

# ---------------------------------------------------------------- step counts
function _steps(items, before, pre_steps, rb)
    first_arm = Dict(a["name"] => 0 for a in rb.raw["arm"])
    after_p = Dict(p => 0 for p in PRECISION_ORDER)
    per_spec = Dict{String, NamedTuple{(:items, :steps), Tuple{Int, Int}}}()
    for it in items
        first_arm[first(it.arms)] += it.horizon
        after_p[_label(it.precision)] += it.horizon
        it.kind === nothing || (c = get(per_spec, it.method, (items = 0, steps = 0)); per_spec[it.method] = (items = c.items + 1, steps = c.steps + it.horizon))
    end
    arms = Dict(name => (items_before = n, steps_before = s, steps_after = first_arm[name]) for (name, (n, s)) in before)
    precs = Dict(p => (steps_before = pre_steps[p], steps_after = after_p[p]) for p in PRECISION_ORDER if pre_steps[p] > 0)
    tb, ta = sum(a.steps_before for a in values(arms)), sum(a.steps_after for a in values(arms))
    return (arms = arms, precisions = precs, total_before = tb, total_after = ta, saving = tb - ta, specs = per_spec)
end

"Canonical bytes of a plan: counts, step counts, chunk table with hashes, panel and block hashes, defaults, code version."
function _plan_text(chunks, items, steps, rb, df, panel, study, cv, omitted = Dict{String, Vector{Tuple{Int, String}}}())
    io = IOBuffer()
    println(io, "chaoslab-plan 1\nstudy=", study, "\nblock_sha256=", block_hash(rb), "\npanel_sha256=", panel.sha256, "\ncode_version=", cv)
    println(io, "hash_version=", HASH_VERSION)
    foreach(k -> println(io, "default ", k, "=", repr(df[k])), sort!(collect(keys(df))))
    println(io, "items=", length(items), " chunks=", length(chunks), " steps_before=", steps.total_before, " steps_after=", steps.total_after)
    for a in rb.raw["arm"]
        s = steps.arms[a["name"]]
        println(io, "arm ", a["name"], " items_before=", s.items_before, " steps_before=", s.steps_before, " steps_after=", s.steps_after)
    end
    for id in sort!(collect(keys(steps.specs)))      # D86: per spec, only when the block has specs
        println(io, "spec ", id, " items=", steps.specs[id].items, " steps=", steps.specs[id].steps)
    end
    for id in sort!(collect(keys(omitted)))
        println(io, "not_applicable ", id, " orbits=", join(sort!(first.(omitted[id])), ","))
    end
    for c in chunks
        println(io, "chunk ", c.ordinal, " ", c.hash, " items=", length(c.items), " arms=", join(c.arms, ","), " ", c.key)
    end
    return String(take!(io))
end

"Canonical text of a plan (the bytes behind `plan.hash`)."
plan_text(plan::RunPlan) = plan.text

# ---------------------------------------------------------------- the plan
"""
    plan_run(prereg, panel; code_version = nothing, defaults = nothing, ntasks = nothing) -> RunPlan

Build the run plan. Refuses (`ArgumentError`, nothing computed) on a preregistration hash different from its header, a panel built
under a hash outside `panel_prereg_sha256_accepted`, a malformed block or a plan without items. `ntasks` is accepted for symmetry
with the panel builder; planning is sequential and never depends on it.
"""
function plan_run(prereg::AbstractString, panel::AbstractString; code_version = nothing, defaults = nothing, ntasks = nothing)
    chk = check_prereg(prereg)
    rb = read_run_block(prereg)
    pnl = _read_panel(panel, rb)
    df = _read_defaults(defaults)
    cv = code_version === nothing ? _code_version_default() : String(code_version)
    study = _study_slug(_lf(read(prereg, String)))
    seed = df["seed"]
    nP = isempty(pnl.orbits) ? 0 : maximum(Int(o["id"]) for o in pnl.orbits)
    bases = vcat([_panel_base(o) for o in pnl.orbits], [_h4_base(h, nP + i, seed) for (i, h) in enumerate(rb.raw["h4"])])
    ctx = (grids = Dict{Tuple{Int, Int}, NamedTuple}(), grid = rb.raw["grid"], bits = df["bigfloat_bits"], cv = cv)
    items, before, pre_steps, omitted = _build_items(rb, bases, ctx)
    for (id, os) in omitted        # an arm that selects no orbit for a spec is an error naming the spec (R12 item 2)
        any(it -> it.method == id, items) || throw(ArgumentError("indicator spec `$id` is undefined for every orbit its arms select (orbits $(join(sort!(first.(os)), ", ")))"))
    end
    isempty(items) && throw(ArgumentError("run plan holds zero items: the arms select no orbit of the panel"))
    for it in items
        it.method == "L6" && _label(it.precision) != "Float64" &&
            throw(ArgumentError("run plan has method L6 at $(it.precision): L6 (library QR) runs in Float64 only"))
    end
    chunks = _build_chunks(items, rb, df)
    steps = _steps(items, before, pre_steps, rb)
    text = _plan_text(chunks, items, steps, rb, df, pnl, study, cv, omitted)
    return RunPlan(items, chunks, (items = length(items), chunks = length(chunks)), steps, bytes2hex(sha256(text)), text, rb, df, pnl,
                   study, cv, chk.header_sha256, Dict{String, Vector{Int}}(id => sort!(first.(os)) for (id, os) in omitted), omitted)
end

_code_version_default() = code_version()
