# Helpers for the D53 runner suite (prefix ex_). Definitions only: nothing here calls the new runner API at load time.
# Reuses ex_helpers.jl (EX, ex_sha256, ex_path, ex_block, ex_prereg_text, ex_with_device, ex_rejects, ...) without changing it.
# Every fixture is generated here from the sizes of spec D53 (<= 12 panel orbits in 3 cells, 3 controls, 2 extra starts).

ex_s53(letter::Char, sub::Integer = 0) = EX_SEED + 28000 + 50 * (letter - 'a' + 1) + sub   # testset D53<letter> uses 50 n + sub, sub < 50

const EX_RUN_SEED = 20261012
const EX_PANEL_HASH = "e"^64                    # the hash under which the fixture panel is "built" (the block accepts it)
const EX_CODE = "ex-code-1"                     # fixed code version string for plans that must not depend on the source

# ---------------------------------------------------------------- independent hashing (platform tool, one process for many texts)
function ex_sha256_many(texts::AbstractVector{<:AbstractString})
    Sys.which("sha256sum") === nothing && return [ex_sha256(t) for t in texts]
    dir = mktempdir()
    try
        names = ["t$i.txt" for i in eachindex(texts)]
        for (n, t) in zip(names, texts)
            write(joinpath(dir, n), t)
        end
        out = read(Cmd(`sha256sum $names`; dir = dir), String)
        got = Dict(String(m.captures[2]) => String(m.captures[1]) for m in eachmatch(r"(?m)^([0-9a-f]{64}) [ *](t\d+\.txt)\s*$", out))
        return [got[n] for n in names]
    finally
        rm(dir; recursive = true, force = true)
    end
end

# ---------------------------------------------------------------- fixture panel
# 12 drawn orbits in cells (1,T), (2,C), (3,S), three controls (elliptic fixed point, hyperbolic fixed point, Y5), ids in E20 order.
# S_D = the first two of each cell plus the controls; S_L = orbit 1. "big" has one cell of 40 orbits and no control.
function ex_panel_orbits(kind::Symbol = :std)
    rng = StableRNG(ex_s53('b', kind === :big ? 1 : 0))
    Ks = ["0.5", "1.5", "2.5"]
    orbs = Dict{String, Any}[]
    function drawn(ki, cls, n)
        for j in 1:n
            id = length(orbs) + 1
            x, p = 2pi * rand(rng), 2pi * rand(rng)
            push!(orbs, Dict{String, Any}("id" => id, "class" => cls, "control" => false, "control_kind" => "none", "k_index" => ki,
                "K_text" => Ks[ki], "K" => parse(Float64, Ks[ki]), "system" => "standard", "x0" => x, "p0" => p,
                "x0_dec" => string(x), "p0_dec" => string(p), "draw_index" => j, "run_seed" => string(EX_RUN_SEED),
                "orbit_seed" => string(EX.orbit_seed(EX_RUN_SEED, id)), "in_SD" => kind === :big || j <= 2, "in_SL" => id == 1))
        end
    end
    if kind === :big
        drawn(1, "T", 40)
        return orbs
    end
    drawn(1, "T", 4); drawn(2, "C", 4); drawn(3, "S", 4)
    function control(kind_, ki, x, p, stab, exponent)
        id = length(orbs) + 1
        push!(orbs, Dict{String, Any}("id" => id, "class" => "control", "control" => true, "control_kind" => kind_, "k_index" => ki,
            "K_text" => Ks[ki], "K" => parse(Float64, Ks[ki]), "system" => "standard", "x0" => x, "p0" => p, "x0_dec" => string(x),
            "p0_dec" => string(p), "run_seed" => string(EX_RUN_SEED), "orbit_seed" => string(EX.orbit_seed(EX_RUN_SEED, id)),
            "in_SD" => true, "in_SL" => false, "period" => 1, "wind_m" => 0, "wind_j" => 0, "stability" => stab, "exponent" => exponent))
    end
    control("fixed", 2, Float64(pi), 0.0, "elliptic", 0.0)
    lam(K) = setprecision(() -> Float64(ex_exponent(BigFloat(2) + K, 1)), BigFloat, 256)       # E4 at tr = 2 + K
    control("fixed", 3, 0.0, 0.0, "hyperbolic", lam(BigFloat("2.5")))
    mu = setprecision(() -> Float64((3 + sqrt(BigFloat(5))) / 2), BigFloat, 256)
    id = length(orbs) + 1
    push!(orbs, Dict{String, Any}("id" => id, "class" => "control", "control" => true, "control_kind" => "y5", "k_index" => 0,
        "system" => "linear", "mu" => mu, "x0" => 0.0, "p0" => 0.0, "x0_dec" => "0.0", "p0_dec" => "0.0",
        "run_seed" => string(EX_RUN_SEED), "orbit_seed" => string(EX.orbit_seed(EX_RUN_SEED, id)), "in_SD" => true, "in_SL" => false,
        "exponent" => setprecision(() -> Float64(log((3 + sqrt(BigFloat(5))) / 2)), BigFloat, 256)))
    return orbs
end

function ex_write_panel(path::AbstractString; kind::Symbol = :std, accepted::AbstractString = EX_PANEL_HASH, orbits = ex_panel_orbits(kind))
    man = Dict{String, Any}("format_version" => 1, "controls" => [o["id"] for o in orbits if o["control"]],
                            "prereg" => Dict{String, Any}("header_sha256" => accepted, "computed_sha256" => accepted))
    open(path, "w") do io
        TOML.print(io, Dict{String, Any}("manifest" => man, "orbit" => orbits))
    end
    return path
end

# ---------------------------------------------------------------- fixture run blocks (dictionaries of the [run] table)
ex_arm(name, orbits, methods, precisions, horizons; kw...) =
    merge(Dict{String, Any}("name" => name, "orbits" => orbits, "methods" => collect(methods), "precisions" => collect(precisions),
                            "horizons" => collect(horizons)), Dict{String, Any}(string(k) => v for (k, v) in kw))

const EX_ALL9 = ["L1", "L2", "L3", "L4", "L4b", "L5", "L5b", "L6", "L9"]
const EX_V = ["L1", "L3", "L4", "L4b", "L5", "L5b", "L6", "L9"]            # the variational methods (no L2)

function ex_rb_common()
    return Dict{String, Any}(
        "panel_prereg_sha256_accepted" => [EX_PANEL_HASH],
        "primary" => Dict{String, Any}("tau" => 1, "eps" => 1e-8, "direction" => "pi/4", "coordinate_mode" => "lifted"),
        "grid" => Dict{String, Any}("points_per_decade" => 5, "first_steps" => 10, "extra_fractions" => [0.1, 0.5]),
        "backend" => Dict{String, Any}("Float64" => "auto", "Double64" => "cpu", "BigFloat" => "cpu"),
        "health" => Dict{String, Any}("fail_fraction" => 0.05, "control_bound" => 0.05, "control_methods" => ["L1", "L4", "L5"],
            "setup_check_steps" => Dict{String, Any}("T" => 8, "R" => 8, "control_elliptic" => 8, "control_hyperbolic" => 3,
                                                     "C" => 3, "S" => 3, "X" => 3, "h4" => 3),
            "setup_check_tol" => 1e-10, "drift_bound" => 1e6, "rerun_attempts" => 1, "stop_failed_fraction" => 0.05,
            "first_chunk_cost_factor" => 2.0, "wall_cap_hours" => 0.125),
        "budget" => Dict{String, Any}(
            "Float64" => Dict{String, Any}("ns_cpu_thread" => 100.0, "ns_gpu" => 10.0, "safety" => 2),
            "Double64" => Dict{String, Any}("ns_cpu_thread" => 2000.0, "safety" => 3),
            "BigFloat" => Dict{String, Any}("ns_cpu_thread" => 20000.0, "safety" => 3)),
        "k19" => Dict{String, Any}("threads" => 4),
        "verify" => Dict{String, Any}("backend_agreement" => Dict{String, Any}("per_orbit_sigma" => 6.0, "mean_standard_errors" => 4.0),
                                      "rerun_sample" => 3),
        "pilot" => Dict{String, Any}("orbits_per_k" => 2, "horizon" => 40, "seed_label" => "ex-pilot"),
        "h4" => [Dict{String, Any}("K" => "1.5", "x0" => "0.001", "p0" => "-0.001", "label" => "ex_h4_a"),
                 Dict{String, Any}("K" => "2.5", "x0" => "0.002", "p0" => "-0.002", "label" => "ex_h4_b")])
end

function ex_rb(kind::Symbol = :small)
    d = ex_rb_common()
    if kind === :small
        d["arm"] = [ex_arm("main", "P", ["L1", "L5"], ["Float64"], [100]; h1_check = true),
                    ex_arm("tau", "S_D", ["L4"], ["Float64"], [100]; tau = [1, 2]),
                    ex_arm("h4f", "h4", ["L5"], ["Float64"], [100]),
                    ex_arm("h4d", "h4", ["L5"], ["Double64"], [40])]
    elseif kind === :noref
        d["arm"] = [ex_arm("main", "P", ["L1"], ["Float64"], [100]; h1_check = true),
                    ex_arm("h4d", "h4", ["L5"], ["Double64"], [40])]
    elseif kind === :full
        d["arm"] = [ex_arm("main", "P", EX_ALL9, ["Float64"], [200]; h1_check = true),
                    ex_arm("dup", "S_D", ["L5"], ["Float64"], [200]),
                    ex_arm("dirs", "S_D", ["L1", "L5"], ["Float64"], [200]; direction = ["pi/4", "0.3"]),
                    ex_arm("epsw", "S_D", ["L2"], ["Float64"], [200]; eps = [1e-6, 1e-8]),
                    ex_arm("tauw", "S_D", ["L4", "L5"], ["Float64"], [200]; tau = [1, 2, 4]),
                    ex_arm("wrap", "S_D", ["L1", "L5"], ["Float64"], [200]; coordinate_mode = ["wrapped"]),
                    ex_arm("h4f", "h4", ["L1", "L5"], ["Float64"], [200]),
                    ex_arm("h4d", "h4", ["L5"], ["Double64"], [100]),
                    ex_arm("bf", "S_L", ["L5"], ["BigFloat"], [200])]
    elseif kind === :alt         # differs from :small in every value of R1 to R5
        d["primary"] = Dict{String, Any}("tau" => 2, "eps" => 1e-6, "direction" => "0.3 + pi/5", "coordinate_mode" => "lifted")
        d["grid"] = Dict{String, Any}("points_per_decade" => 4, "first_steps" => 6, "extra_fractions" => [0.25, 0.75])
        d["backend"] = Dict{String, Any}("Float64" => "cpu", "Double64" => "cpu", "BigFloat" => "cpu")
        h = d["health"]
        h["fail_fraction"] = 0.1; h["control_bound"] = 0.2; h["control_methods"] = ["L3", "L6"]
        h["setup_check_steps"] = Dict{String, Any}("T" => 4, "R" => 4, "control_elliptic" => 4, "control_hyperbolic" => 2,
                                                   "C" => 2, "S" => 2, "X" => 2, "h4" => 2)
        h["setup_check_tol"] = 1e-9; h["drift_bound"] = 1e5; h["stop_failed_fraction"] = 0.1
        h["first_chunk_cost_factor"] = 3.0; h["wall_cap_hours"] = 8.0
        d["budget"]["Float64"] = Dict{String, Any}("ns_cpu_thread" => 50.0, "ns_gpu" => 5.0, "safety" => 3)
        d["k19"] = Dict{String, Any}("threads" => 2)
        d["verify"] = Dict{String, Any}("backend_agreement" => Dict{String, Any}("per_orbit_sigma" => 5.0, "mean_standard_errors" => 3.0),
                                        "rerun_sample" => 2)
        d["pilot"] = Dict{String, Any}("orbits_per_k" => 3, "horizon" => 24, "seed_label" => "ex-pilot-alt")
        d["h4"] = [Dict{String, Any}("K" => "pi/2", "x0" => "1e-3", "p0" => "-1e-3", "label" => "alt_a"),
                   Dict{String, Any}("K" => "2pi", "x0" => "1e-3", "p0" => "-1e-3", "label" => "alt_b")]
        d["arm"] = [ex_arm("alpha", "P", ["L2", "L3", "L6"], ["Float64"], [60]; h1_check = true),
                    ex_arm("beta", "S_L", ["L9"], ["Float64"], [80]; direction = ["0.3"], coordinate_mode = ["wrapped"]),
                    ex_arm("gamma", "h4", ["L5"], ["Double64"], [24]),
                    ex_arm("delta", "S_D", ["L1"], ["Float64"], [60]; tau = [2, 6])]
    elseif kind === :c40         # 40 chunks: 9 methods x 4 cells (36) + h4 x 2 methods x 2 horizons (4)
        d["arm"] = [ex_arm("main", "P", EX_ALL9, ["Float64"], [40]), ex_arm("h4f", "h4", ["L1", "L5"], ["Float64"], [40, 60])]
    elseif kind === :c20         # 20 chunks: 5 methods x 4 cells
        d["arm"] = [ex_arm("main", "P", ["L1", "L3", "L4", "L5", "L6"], ["Float64"], [40])]
    elseif kind === :big         # one cell of 40 orbits (panel :big)
        d["arm"] = [ex_arm("main", "P", ["L1"], ["Float64"], [60])]
        delete!(d, "h4")
    else
        error("ex_rb: unknown kind $kind")
    end
    return d
end

# deep copy with a path of keys changed: ex_rb_with(d, ("health", "fail_fraction") => 0.1)
function ex_rb_with(d::AbstractDict, edits::Pair...)
    d = deepcopy(d)
    for (path, v) in edits
        cur = d
        for k in path[1:(end - 1)]
            cur = cur isa AbstractVector ? cur[k] : cur[k]
        end
        v === :delete ? delete!(cur, path[end]) : (cur[path[end]] = v)
    end
    return d
end

# ---------------------------------------------------------------- TOML emitter with a chosen textual key order (TOML.print sorts)
function ex_emit_run(run::AbstractDict; rev::Bool = false)
    io = IOBuffer()
    ord(ks) = rev ? reverse(ks) : ks
    isarr(v) = v isa AbstractVector && !isempty(v) && first(v) isa AbstractDict
    function scalar(v)
        v isa AbstractString && return repr(v)
        v isa Bool && return string(v)
        v isa Integer && return string(v)
        v isa AbstractFloat && return (s = repr(v); occursin(r"[.eE]", s) ? s : s * ".0")
        error("ex_emit_run: unsupported value $(typeof(v))")
    end
    value(v) = v isa AbstractVector ? "[" * join(value.(v), ", ") * "]" : scalar(v)
    function emit(path, d, double)
        ks = sort(collect(keys(d)))
        plain = [k for k in ks if !(d[k] isa AbstractDict) && !isarr(d[k])]
        println(io, double ? "[[" : "[", path, double ? "]]" : "]")
        for k in ord(plain)
            println(io, k, " = ", value(d[k]))
        end
        println(io)
        for k in ord([k for k in ks if d[k] isa AbstractDict])
            emit(path * "." * k, d[k], false)
        end
        for k in ord([k for k in ks if isarr(d[k])]), e in d[k]
            emit(path * "." * k, e, true)
        end
    end
    emit("run", run, false)
    return String(take!(io))
end

function ex_prereg_text_run(run::AbstractDict; panel::AbstractDict = ex_block(), rev::Bool = false)
    text = ex_prereg_text(panel)
    return replace(text, "End of file." => "## Amendment A4\n\n```toml\n" * ex_emit_run(run; rev = rev) * "```\n\nEnd of file.")
end

function ex_write_prereg_run(path::AbstractString, run::AbstractDict; panel::AbstractDict = ex_block(), rev::Bool = false)
    text = ex_prereg_text_run(run; panel = panel, rev = rev)
    h = ex_sha256(ex_prepare(text))
    write(path, replace(text, "0"^64 => h; count = 1))
    return path
end

# fixtures: files in the suite's temporary directory, cached by tag
const EX_FIXES = Dict{Any, Any}()
function ex_fix(kind::Symbol; panel::Symbol = kind === :big ? :big : :std, tag = string(kind, "-", panel))
    get!(EX_FIXES, tag) do
        run = ex_rb(kind)
        pp = ex_write_panel(ex_path("fix-$tag-panel.toml"); kind = panel)
        pr = ex_write_prereg_run(ex_path("fix-$tag-prereg.md"), run)
        (; kind, run, panel = pp, prereg = pr, tag, orbits = ex_panel_orbits(panel))
    end
end
function ex_fix_with(run::AbstractDict, tag::AbstractString; panel::Symbol = :std, rev::Bool = false, accepted = EX_PANEL_HASH)
    pp = ex_write_panel(ex_path("fx-$tag-panel.toml"); kind = panel, accepted = accepted)
    pr = ex_write_prereg_run(ex_path("fx-$tag-prereg.md"), run; rev = rev)
    return (; run, panel = pp, prereg = pr, tag, orbits = ex_panel_orbits(panel))
end

ex_plan(fx; kw...) = EX.plan_run(fx.prereg, fx.panel; code_version = EX_CODE, kw...)

# campaign with the defaults every test needs: thread check overridden, dirty tree allowed (the repository under test is dirty)
function ex_campaign(fx; store, manifest = store * "-manifest.json", threads = 4, hooks = EX.RunHooks(), kw...)
    # capacity 64 orbits per launch unless a test sets it: the chunk-size attribute then does not depend on the free memory of the host
    # no clock given: a zero-cost clock, so the first-chunk cost cap (R5 item 5) cannot pause arms on a slow host; the cap tests inject their own
    hooks.clock === nothing && (hooks = EX.RunHooks(clock = ex_clock(info -> 0.0), kill = hooks.kill, mutate = hooks.mutate, observe = hooks.observe))
    return EX.run_campaign(fx.prereg, fx.panel; store = store, manifest = manifest, threads = threads, thread_override = true,
                           allow_dirty = true, code_version = EX_CODE, capacity = 64, hooks = hooks,
                           manifest_reset = true, kw...)   # D53 review M4: stores copied without their manifest still resume; ex_review.jl tests the refusal
end
ex_store(name) = (d = ex_path(name); rm(d; recursive = true, force = true); d)

# ---------------------------------------------------------------- kill and fault hooks
struct ExKill <: Exception
    point::Symbol
    ordinal::Int
end
ex_has_kill(::ExKill) = true
ex_has_kill(e::CompositeException) = any(ex_has_kill, e.exceptions)
ex_has_kill(e::TaskFailedException) = ex_has_kill(e.task.result)
ex_has_kill(e::LoadError) = ex_has_kill(e.error)
ex_has_kill(::Any) = false
# true when f was killed by an ExKill, false when it returned; anything else is rethrown
function ex_killed(f)
    try
        f()
    catch e
        ex_has_kill(e) && return true
        rethrow()
    end
    return false
end
ex_kill_at(point::Symbol, ordinal::Integer) = (p, info) -> (p === point && get(info, :ordinal, 0) == ordinal && throw(ExKill(p, ordinal)); nothing)
ex_kill_after(k::Integer) = (p, info) -> (p === :chunk_done && info.chunks_done == k && throw(ExKill(p, k)); nothing)

# stateful fake clock: the clock hook returns the time at an event; a group takes `cost(info)` seconds
function ex_clock(cost; start = 1.7e9)
    t = Ref(start)
    return (event::Symbol, info) -> (event === :group_end && (t[] += cost(info)); t[])
end

# ---------------------------------------------------------------- reading stores
ex_chunk_files(store) = sort(filter(f -> occursin(r"^chunk-\d{5}\.h5$", basename(f)), readdir(store; join = true)))
ex_tmp_files(store) = sort(filter(f -> endswith(f, ".tmp"), readdir(store; join = true)))
ex_ordinal(f) = parse(Int, match(r"chunk-(\d{5})\.h5", basename(f)).captures[1])

# canonical bytes of a read-back chunk: datasets always; attributes except host and time fields on request (I12, I14)
const EX_ROOT_SKIP = ("created_utc", "cpu", "threads", "session")      # the session index differs between a killed and an uninterrupted run
function ex_canon(ch; attrs::Bool = false)
    io = IOBuffer()
    if attrs
        for k in sort(collect(keys(ch.attrs)))
            k in EX_ROOT_SKIP || println(io, k, "=", repr(ch.attrs[k]))
        end
    end
    for o in sort(ch.orbits; by = o -> o.id)
        println(io, "orbit ", o.id)
        if attrs
            for k in sort(collect(keys(o.attrs)))
                println(io, "  ", k, "=", repr(o.attrs[k]))
            end
        end
        write(io, o.t); write(io, o.state)
        o.winding === nothing || write(io, o.winding)
        for k in sort(collect(keys(o.lyap)))
            println(io, "lyap/", k); write(io, o.lyap[k])
        end
        for k in sort(collect(keys(o.diagnostics)))
            println(io, "diag/", k); write(io, o.diagnostics[k]); println(io, o.thresholds[k], " ", o.first_exceedance[k])
        end
        for k in sort(collect(keys(o.cost)))
            c = o.cost[k]; println(io, "cost/", k); write(io, Float64[c[1], c[2], c[5]])     # wall and thread seconds (3, 4) excluded
        end
    end
    return take!(io)
end
function ex_store_canon(store; attrs::Bool = false)
    io = IOBuffer()
    for f in ex_chunk_files(store)
        println(io, basename(f)); write(io, ex_canon(EX.read_chunk(f); attrs = attrs))
    end
    return take!(io)
end

# ---------------------------------------------------------------- reference implementations written from the spec
# R2 item 1 at 1024 bits: the independent grid.
function ex_ref_grid(N::Integer, tau::Integer; ppd::Integer, first_steps::Integer, fracs)
    return setprecision(BigFloat, 1024) do
        half_up(x) = BigInt(floor(x + BigFloat(1) / 2))
        targets = Set{BigInt}(BigInt.(1:first_steps))
        j = 0
        while true
            v = BigFloat(10)^(BigFloat(j) / ppd)
            v <= N || break
            push!(targets, half_up(v))
            j += 1
        end
        for f in fracs
            push!(targets, half_up(parse(BigFloat, string(f)) * N))
        end
        push!(targets, BigInt(N))
        sort!(unique(Int(tau) * cld(Int(t), Int(tau)) for t in targets))
    end
end

ex_ref_cell(o) = o["control"] ? "controls" : "k$(o["k_index"]):$(o["class"])"

# R1 items 2 and 3: the work items of a block on a panel, merged by content (backend, arm names excluded). Key =
# (orbit id, method, precision label, horizon, tau, eps or nothing, direction text, coordinate mode); value = arms in block order.
function ex_ref_items(run::AbstractDict, orbits)
    prim = run["primary"]
    nP = isempty(orbits) ? 0 : maximum(o["id"] for o in orbits)
    h4ids = collect((nP + 1):(nP + length(get(run, "h4", []))))
    subsets = Dict("P" => [o["id"] for o in orbits], "S_D" => [o["id"] for o in orbits if o["in_SD"]],
                   "S_L" => [o["id"] for o in orbits if o["in_SL"]], "h4" => h4ids)
    items = Dict{Tuple, Vector{String}}()
    firstpos = Dict{Tuple, Int}()
    before = Dict{String, NamedTuple{(:items, :steps), Tuple{Int, Int}}}()
    for arm in run["arm"]
        n = 0; st = 0
        for prec in arm["precisions"], N in arm["horizons"], mode in get(arm, "coordinate_mode", [prim["coordinate_mode"]]),
            dir in get(arm, "direction", [prim["direction"]]), tau in get(arm, "tau", [prim["tau"]]), m in arm["methods"],
            eps in (m == "L2" ? get(arm, "eps", [prim["eps"]]) : Any[nothing]), id in subsets[arm["orbits"]]
            key = (id, m, prec, N, tau, eps, dir, mode)
            haskey(items, key) || (firstpos[key] = length(firstpos) + 1)
            push!(get!(items, key, String[]), arm["name"])
            n += 1; st += N
        end
        before[arm["name"]] = (items = n, steps = st)
    end
    return (; items, before, nunique = length(items), steps_after = sum(k[4] for k in keys(items); init = 0),
            steps_before = sum(v.steps for v in values(before); init = 0))
end

# D53h reference table (R4), written from the spec; `dev` is nothing or an ExDev, `rec` nothing or (threads, gpu_faster)
struct ExDev <: ex_F.AbstractDevice
    loadable::Bool
    functional::Bool
    compiles::Bool
end
ex_F.probe(d::ExDev, ::Type) = (loadable = d.loadable, functional = d.functional, device = "FakeGPU", compiles = d.compiles, gpu_faster = true)
ex_Tg.device_backend(::ExDev) = ex_Tg.KernelAbstractions.CPU()
function ex_expect_backend(req::Symbol, T::Type, cpu_only::Bool, dev, rec, threads::Integer)
    R = Chaoslab.FALLBACK_REASONS
    req === :cpu && return (:cpu, "")
    T === BigFloat && return (:cpu, R[3])
    cpu_only && return (:cpu, R[6])
    dev.loadable || return (:cpu, R[1])
    dev.functional || return (:cpu, R[2])
    dev.compiles || return (:cpu, R[3])
    if req === :auto
        ok = rec !== nothing && rec.threads == threads && rec.gpu_faster
        ok || return (:cpu, R[7])
    end
    return (:gpu, "")
end

# ---------------------------------------------------------------- direct estimator calls (O9), written from the lower-layer API
const EX_RED = Dict("L1" => () -> Chaoslab.Benettin(), "L3" => () -> Chaoslab.QRSpectrum(Chaoslab.CGS()),
                    "L4" => () -> Chaoslab.QRSpectrum(Chaoslab.MGS()), "L4b" => () -> Chaoslab.QRSpectrum(Chaoslab.RGS()),
                    "L5" => () -> Chaoslab.QRSpectrum(Chaoslab.HouseholderQR()), "L5b" => () -> Chaoslab.QRSpectrum(Chaoslab.AugmentedQR()),
                    "L6" => () -> Chaoslab.QRSpectrum(Chaoslab.LapackQR()), "L9" => () -> Chaoslab.SVDSpectrum())
ex_vectors(method) = method in ("L1", "L2") ? 1 : 2

# precision label of an item -> (element type, bits); runs f inside the right BigFloat precision
function ex_with_prec(f, label::AbstractString)
    T, bits = Chaoslab.parse_precision(label)
    return bits === nothing ? f(T) : setprecision(() -> f(T), BigFloat, bits)
end

# system and start of a plan item in element type T (inside the right BigFloat precision); wrapped mode gives the 4-state
function ex_sys_x(item, T)
    sys = if item.system == "linear"
        Chaoslab.y5_control(parse(Float64, item.param)).sys
    elseif item.coordinate_mode == "wrapped"
        Chaoslab.WrappedStandardMap(Chaoslab.exact(T, item.param))
    else
        Chaoslab.StandardMap(Chaoslab.exact(T, item.param))
    end
    x = SVector(Chaoslab.exact(T, item.x0), Chaoslab.exact(T, item.p0))
    return sys, (item.coordinate_mode == "wrapped" && item.system != "linear" ? Chaoslab.torus_state(x) : x)   # Y5 runs in plane coordinates
end

# states of a plan item after the steps ns (sorted), by direct iteration of the map, rounded once to binary64
function ex_iter_states(item, ns::AbstractVector{<:Integer})
    return ex_with_prec(item.precision) do T
        sys, x = ex_sys_x(item, T)
        out = Vector{Vector{Float64}}(); n = 0
        for target in ns
            while n < target
                x = Chaoslab.forward(sys, x); n += 1
            end
            push!(out, Float64.(collect(x)))
        end
        out
    end
end

# lambda_1 .. lambda_k of a plan item after n steps by the verified estimator run directly (no recording; the item's own
# deviation frame is the module input checked separately)
function ex_direct(item, n::Integer)
    return ex_with_prec(item.precision) do T
        sys, x = ex_sys_x(item, T)
        W = EX.deviation_frame(T, item.direction, ex_vectors(item.method))
        if item.method == "L2"
            u0 = SVector{2, T}(W[1, 1], W[2, 1])
            _, _, fin = Chaoslab.propagate_shadow_orbit(sys, x, u0, n, item.tau, Chaoslab.exact(T, item.eps))
            return Float64[fin.X2]
        end
        res = Chaoslab.propagate(sys, [x], n; k = ex_vectors(item.method), tau = item.tau, seed = parse(Int, item.run_seed),
                                 reducers = (EX_RED[item.method](),), W0 = [W])
        return Float64.(collect(res.finals[1][1].lambda))
    end
end

# factor c of lambda_1(n) = ln mu + ln(c) / n for a 2 x 2 symplectic hyperbolic A and a unit start v, at 256 bits (O11)
function ex_c_factor(A, v)
    return setprecision(BigFloat, 256) do
        M = BigFloat.(Matrix(A)); t = M[1, 1] + M[2, 2]
        mu = (t + sqrt(t^2 - 4)) / 2
        e = [M[1, 2], mu - M[1, 1]]; l = [M[2, 1], mu - M[1, 1]]      # right and left eigenvectors of mu
        vb = BigFloat.(collect(v))
        Float64(abs(l[1] * vb[1] + l[2] * vb[2]) * sqrt(e[1]^2 + e[2]^2) / abs(l[1] * e[1] + l[2] * e[2]))
    end
end

# 256-bit finite-time exponent of a linear map A started along unit vector v: ln ||A^n v|| / n (exact for the stored matrix)
function ex_lin_lambda(A, v, n::Integer)
    return setprecision(BigFloat, 256) do
        M = BigFloat.(Matrix(A)); w = BigFloat.(collect(v))
        for _ in 1:n
            w = M * w
        end
        Float64(log(sqrt(w[1]^2 + w[2]^2)) / n)
    end
end

# closed-form direction text (decimals exact, pi at the current BigFloat precision), as the preregistration writes them
function ex_dir_big(s::AbstractString)
    t = replace(s, r"(\d)pi" => s"\1*pi")
    t = replace(t, r"(?<![\w.])(\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)" => s"parse(BigFloat,\"\1\")")
    t = replace(t, r"\bpi\b" => "BigFloat(pi)")
    return eval(Meta.parse(t))
end

# ---------------------------------------------------------------- fault injection through the mutate hook
# f(orbits, info) runs on the chunk with ordinal `ord` at the listed attempts; every call is logged
function ex_fault(f; ord, attempts = (1,), log = nothing)
    return (info, orbits) -> begin
        log === nothing || push!(log, (ordinal = info.ordinal, attempt = info.attempt, backend = info.backend))
        info.ordinal == ord && info.attempt in attempts && f(orbits, info)
        nothing
    end
end
ex_find_chunk(plan; method = nothing, precision = nothing, cell = nothing, arm = nothing) =
    first(c for c in plan.chunks if (method === nothing || c.run_type.method == method) &&
         (precision === nothing || c.run_type.precision == precision) && (cell === nothing || c.cell == cell) &&
         (arm === nothing || arm in c.arms))
