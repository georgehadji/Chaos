# Pilot mode (R9) and the CPU reference rerun (O12). The pilot uses starts that are not panel orbits, writes outside the study store
# and reports costs and projections; its numbers are not study data.

const PILOT_NOTE = "pilot measurements: not study data (Amendment A4, item 2); no pilot start is a panel orbit and no file is in the study store"
const PILOT_MIN_DISTANCE = 1e-6      # torus distance below which a drawn start counts as a panel start (R9 item 2)
const PILOT_MAX_DRAWS = 1000

function _pilot_seed(label, k)
    h = sha256("chaoslab-pilot-1|$(label)|$(k)")
    return foldl((acc, b) -> (acc << 8) | UInt64(b), h[1:8]; init = UInt64(0))
end

# Panel starts as binary64 pairs, the extra starts included.
function _panel_points(plan::RunPlan)
    pts = [(Float64(o["x0"]), Float64(o["p0"])) for o in plan.panel.orbits]
    for h in plan.block.raw["h4"]
        push!(pts, setprecision(() -> (Float64(_eval_exact(h["x0"])), Float64(_eval_exact(h["p0"]))), BigFloat, 256))
    end
    return pts
end

_k_indices(plan::RunPlan) = sort!(unique(Int(o["k_index"]) for o in plan.panel.orbits if o["k_index"] > 0))
_k_texts(plan::RunPlan) = Dict(Int(o["k_index"]) => _decimal_text(String(o["K_text"])) for o in plan.panel.orbits if o["k_index"] > 0)

"Pilot starts: `orbits_per_k` per K index, drawn uniformly on the torus from a stream keyed by (seed_label, K index); none within 1e-6 of a panel start."
function pilot_starts(plan::RunPlan)
    p = plan.block.raw["pilot"]
    pts = _panel_points(plan)
    out = Dict{String, Any}[]
    for k in _k_indices(plan)
        rng = StableRNG(_pilot_seed(p["seed_label"], k))
        got, draws = 0, 0
        while got < p["orbits_per_k"]
            (draws += 1) > PILOT_MAX_DRAWS * p["orbits_per_k"] &&
                throw(ArgumentError("pilot: no start away from the panel starts for K index $k after $PILOT_MAX_DRAWS draws"))
            x, y = 2pi * rand(rng), 2pi * rand(rng)
            minimum(torus_distance((x, y), q) for q in pts) > PILOT_MIN_DISTANCE || continue
            got += 1
            push!(out, Dict{String, Any}("k_index" => k, "x0" => x, "p0" => y))
        end
    end
    return out
end

# The pilot items: every run type of the plan (horizon replaced by the pilot horizon) on every start, in groups of one run type.
function pilot_items(plan::RunPlan, starts)
    rb, df = plan.block, plan.defaults
    ktext = _k_texts(plan)
    nid = maximum(it.orbit for it in plan.items)
    ctx0 = (grids = Dict{Tuple{Int, Int}, NamedTuple}(), grid = rb.raw["grid"], bits = df["bigfloat_bits"], cv = plan.code_version)
    types = unique((it.method, _label(it.precision), it.tau, it.eps, it.direction, it.coordinate_mode) for it in plan.items)
    items, groups = NamedTuple[], Vector{Int}[]
    for (m, prec, tau, eps, dir, mode) in types
        ids = Int[]
        for (j, st) in enumerate(starts)
            id = nid + j
            base = (id = id, system = "standard", param = ktext[st["k_index"]], x0 = string(st["x0"]), p0 = string(st["p0"]), run_seed = string(df["seed"]),
                    orbit_seed = string(orbit_seed(df["seed"], id)), design_class = "pilot", setup_class = "T", cell = "pilot:k$(st["k_index"])",
                    exponent = nothing, in_SD = false, in_SL = false, k_index = st["k_index"])
            push!(items, _make_item(length(items) + 1, base, m, prec, rb.raw["pilot"]["horizon"], tau, eps === nothing ? nothing : parse(Float64, eps), dir, mode, ctx0))
            push!(ids, length(items))
        end
        push!(groups, ids)
    end
    return items, groups
end

# First step at which an extended-precision orbit leaves the Float64 orbit by more than `tol`, up to `N` steps (nothing if never).
function _departure(::Type{T}, bits, Ktext, x0, y0, tol, N) where {T}
    sys64, a = StandardMap(exact(Float64, Ktext)), SVector(x0, y0)
    return setprecision(BigFloat, bits) do
        sysT, b = StandardMap(exact(T, Ktext)), SVector(T(x0), T(y0))
        for n in 1:N
            a, b = forward(sys64, a), forward(sysT, b)
            max(abs(a[1] - Float64(b[1])), abs(a[2] - Float64(b[2]))) > tol && return n
        end
        return nothing
    end
end

function _first_departure(plan::RunPlan, starts)
    tol, N, bits = plan.block.raw["health"]["setup_check_tol"], plan.block.raw["pilot"]["horizon"], plan.defaults["bigfloat_bits"]
    ktext = _k_texts(plan)
    out = Dict{String, Any}()
    for k in _k_indices(plan)
        mine = [st for st in starts if st["k_index"] == k]
        first_of(T) = (v = [d for d in (_departure(T, bits, ktext[k], st["x0"], st["p0"], tol, N) for st in mine) if d !== nothing]; isempty(v) ? nothing : minimum(v))
        out[string(k)] = Dict{String, Any}("Double64" => first_of(Double64), "BigFloat" => first_of(BigFloat))
    end
    return out
end

# Indices of a group by cell, in order of appearance.
function _cells(its)
    d = Dict{String, Vector{Int}}()
    for (i, it) in enumerate(its)
        push!(get!(d, it.cell, Int[]), i)
    end
    return sort!(collect(values(d)); by = first)
end

# One pilot group: compute, time by the clock hook, write one chunk file per K, return what the report needs.
function _pilot_group!(plan, ctx, hooks, store, its, gi, ordinal0, o)
    steps = sum(it.horizon for it in its)
    row = o.table[(_label(first(its).precision), first(its).method)]
    info = (group = gi, chunks = Int[], arms = String[], orbits = length(its), steps = steps,
            backend = row.backend === :gpu && !o.lost[] ? :gpu : :cpu)
    ntask = tasks_for(first(its).precision, o.threads, steps)
    prime_frames!(ctx, its)
    warm_up!(ctx, its, ntask; dev = info.backend === :gpu)
    times = zeros(length(its))
    t0 = _now(hooks, :group_start, info)
    r0 = time_ns()
    recs, prov = compute_group(ctx, its, ntask, row, o.lost; times = times)
    t1 = _now(hooks, :group_end, info)
    real_wall = (time_ns() - r0) / NS_PER_S
    workers = prov.backend === :gpu ? 1 : min(ntask, length(its))
    cells = _cells(its)
    for (n, idx) in enumerate(cells)
        it, key = its[first(idx)], _chunk_key(its[first(idx)])
        hash = bytes2hex(sha256(key * join(sort(String[its[i].hash for i in idx]), ",")))
        attrs = _attrs(plan, it, o.table[(_label(it.precision), it.method)], ordinal0 + n, key, hash, String[], length(idx); git = o.git.sha,
                       threads = o.threads, capacity = o.capacity, utc = _utc(t1), status = "ok", attempt = 1, first_reasons = String[],
                       failure_reasons = String[], setup = "not applicable", pilot = true, prov = prov, dirty = o.dirty)
        write_chunk(chunk_path(store, ordinal0 + n), attrs, recs[idx])
    end
    return (its = its, recs = recs, steps = steps, wall = t1 - t0, workers = workers, cells = length(cells), backend = prov.backend,
            eff = real_wall > 0 ? clamp(sum(times) / (real_wall * workers), 0.0, 1.0) : 1.0)
end

function _pilot_costs(groups)
    acc = Dict{Tuple{String, String}, Vector{Float64}}()
    for g in groups
        it = first(g.its)
        c = get!(acc, (it.method, _label(it.precision)), zeros(3))
        c .+= (g.steps, g.wall, g.wall * g.workers)
    end
    return acc
end

function _pilot_arms(plan::RunPlan, cost, eff, threads)
    rb = plan.block
    h = rb.raw["health"]
    out = Dict{String, Any}()
    for arm in rb.raw["arm"]
        steps = plan.steps.arms[arm["name"]].steps_after
        thr = sum(it.horizon * cost[(it.method, _label(it.precision))][3] / cost[(it.method, _label(it.precision))][1]
                  for it in plan.items if first(it.arms) == arm["name"]; init = 0.0)
        line = budget_line_ns(rb, arm)
        out[arm["name"]] = Dict{String, Any}("planned_unique_steps" => steps, "projected_thread_seconds" => thr,
            "projected_wall_seconds" => thr / (threads * max(eff, 0.001)), "budget_line_ns" => line,
            "first_chunk_cap_ns" => h["first_chunk_cost_factor"] * line, "wall_cap_hours" => h["wall_cap_hours"])
    end
    return out
end

function _pilot_report(plan, starts, groups, o, dirty)
    cost = _pilot_costs(groups)
    bk = Dict((first(g.its).method, _label(first(g.its).precision)) => g.backend for g in groups)
    eff = isempty(groups) ? 1.0 : sum(g.eff for g in groups) / length(groups)
    share = Dict{String, Any}()
    for m in unique(it.method for it in plan.items)
        rs = [r for g in groups if first(g.its).method == m for r in g.recs]
        share[m] = isempty(rs) ? 0.0 : count(exceeded, rs) / length(rs)
    end
    return Dict{String, Any}("not_study_data" => PILOT_NOTE, "starts" => starts, "plan_sha256" => plan.hash,
        "cost" => [Dict{String, Any}("method" => m, "precision" => p, "backend" => String(bk[(m, p)]), "thread_seconds_per_orbit_step" => c[3] / c[1],
                                     "wall_seconds_per_orbit_step" => c[2] / c[1]) for ((m, p), c) in sort!(collect(cost); by = first)],
        "arms" => _pilot_arms(plan, cost, eff, o.threads), "first_departure" => _first_departure(plan, starts),
        "k19" => Dict{String, Any}("threads" => plan.block.raw["k19"]["threads"], "device_record" => _record_dict(o.records[1]),
                                   "prior_record" => _record_dict(o.records[2])),
        "diagnostics_share" => share, "parallel_efficiency" => eff, "threads" => o.threads,
        "dirty_tree" => Dict{String, Any}("override" => dirty !== nothing, "paths" => dirty === nothing ? String[] : dirty.paths))
end

"""
    run_pilot(prereg, panel; out, threads, thread_override, allow_dirty, repo, code_version, capacity, device_record, prior_record, hooks)

Run the pilot: starts that are not panel orbits, every run type of the plan at `[run.pilot].horizon`. Writes `<out>/pilot-report.json`
and the pilot store `<out>/store`; returns `(report, path, starts)`. Nothing is written to the study store.
"""
function run_pilot(prereg::AbstractString, panel::AbstractString; out::AbstractString, threads::Integer = Threads.nthreads(),
                   thread_override::Bool = false, allow_dirty::Bool = false, repo = nothing, code_version = nothing, defaults = nothing,
                   capacity = nothing, device_record = nothing, prior_record = nothing, hooks::RunHooks = RunHooks())
    plan = plan_run(prereg, panel; code_version = code_version, defaults = defaults)
    _check_threads(plan, threads, thread_override)
    git, dirty = _check_tree(repo === nothing ? _root() : repo, allow_dirty)
    rec = device_record === nothing ? gpu_record() : device_record
    starts = pilot_starts(plan)
    items, groups = pilot_items(plan, starts)
    ctx = make_ctx(plan, items)
    cap = capacity === nothing ? _launch_capacity(plan, Int(Sys.free_memory()), 0.5) : max(1, Int(capacity))
    o = (threads = Int(threads), capacity = cap, git = git, dirty = dirty !== nothing, records = (rec, prior_record), lost = Ref(false), table = backend_table(plan; threads = threads, record = rec))
    store = joinpath(out, "store")
    mkpath(store)
    Threads.nthreads() > 1 && LinearAlgebra.BLAS.set_num_threads(1)
    done, ordinal = [], 0
    for (gi, ids) in enumerate(groups)
        g = _pilot_group!(plan, ctx, hooks, store, NamedTuple[items[i] for i in ids], gi, ordinal, o)
        ordinal += g.cells
        push!(done, g)
    end
    report = _pilot_report(plan, starts, done, o, dirty)
    path = _write_json(joinpath(out, "pilot-report.json"), report)
    return (report = report, path = path, starts = report["starts"])
end

# ---------------------------------------------------------------- CPU reference rerun (O12)
# Items of an arm in the order of a hash keyed by (arm, item hash): the sample depends on the plan only, never on threads.
function _sample_items(plan::RunPlan, arm, n)
    ids = [it.id for it in plan.items if arm in it.arms]
    keyed = sort!(ids; by = id -> bytes2hex(sha256(arm * "|" * plan.items[id].hash)))
    return first(keyed, min(n, length(keyed)))
end

function _same_bits(a::OrbitRecord, b::OrbitRecord)
    ok = _same(a.t, b.t) && _same(a.state, b.state) && keys(a.lyap) == keys(b.lyap) && keys(a.diagnostics) == keys(b.diagnostics) &&
         (a.winding === nothing ? b.winding === nothing : b.winding !== nothing && _same(a.winding, b.winding))
    return ok && all(m -> _same(a.lyap[m], b.lyap[m]), keys(a.lyap)) && all(m -> _same(a.diagnostics[m], b.diagnostics[m]), keys(a.diagnostics))
end

_max_diff(a::OrbitRecord, b::OrbitRecord) = maximum((maximum(abs.(a.lyap[m] .- b.lyap[m])) for m in keys(a.lyap)); init = 0.0)

"""
    rerun_check(prereg, panel, store; code_version, threads, thread_override, allow_dirty, sample) -> Vector

Rerun `[run.verify].rerun_sample` items per arm (or `sample`) on the CPU and compare them with the stored chunk files: one
`(item, arm, backend, bitwise, max_abs_diff)` per item; items whose chunk file is missing are left out.
"""
function rerun_check(prereg::AbstractString, panel::AbstractString, store::AbstractString; code_version = nothing, defaults = nothing,
                     threads::Integer = Threads.nthreads(), thread_override::Bool = false, allow_dirty::Bool = false, repo = nothing, sample = nothing)
    plan = plan_run(prereg, panel; code_version = code_version, defaults = defaults)
    _check_threads(plan, threads, thread_override)
    _check_tree(repo === nothing ? _root() : repo, allow_dirty)
    ctx = make_ctx(plan)
    n = sample === nothing ? plan.block.raw["verify"]["rerun_sample"] : sample
    chunk_of = Dict(i => c.ordinal for c in plan.chunks for i in c.items)
    files = Dict{Int, Any}()
    out = NamedTuple[]
    for arm in plan.block.raw["arm"], id in _sample_items(plan, arm["name"], n)
        ord = chunk_of[id]
        path = chunk_path(store, ord)
        isfile(path) || continue
        ch = get!(() -> read_chunk(path), files, ord)
        it = plan.items[id]
        stored = only(o for o in ch.orbits if o.id == it.orbit)
        again = run_item(it, ctx)
        push!(out, (item = id, arm = arm["name"], backend = ch.attrs["backend"], bitwise = _same_bits(stored, again), max_abs_diff = _max_diff(stored, again)))
    end
    return out
end
