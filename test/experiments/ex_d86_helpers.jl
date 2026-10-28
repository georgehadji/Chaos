# Helpers for the Amendment D86 suite (prefix d86_). Definitions only: nothing here calls the D86 API at load time, so a missing function
# errors the testsets that need it and nothing else. Reuses ex_helpers.jl and ex_run_helpers.jl without changing them.
# Spec: workspaces/module-build/stages/01-spec/output/experiments-spec.md, "Amendment D86". API assumed: experiments-test-plan.md, "API assumed by the D86 tests".

d86_seed(letter::Char, sub::Integer = 0) = EX_SEED + 29000 + 50 * (letter - 'a' + 1) + sub     # testset D86<letter> uses 29000 + 50 n + sub, sub < 50

# ---------------------------------------------------------------- spec tables (dictionaries of one [[run.indicator]] table)
d86_budget() = Dict{String, Any}("Float64" => Dict{String, Any}("ns_cpu_thread" => 100.0, "safety" => 2),
                                 "Double64" => Dict{String, Any}("ns_cpu_thread" => 2000.0, "safety" => 3))
d86_obs() = [Dict{String, Any}("kind" => "cos", "coordinate" => 1, "frequency" => 1), Dict{String, Any}("kind" => "sin", "coordinate" => 1, "frequency" => 1),
             Dict{String, Any}("kind" => "cos", "coordinate" => 2, "frequency" => 2), Dict{String, Any}("kind" => "lifted", "coordinate" => 2)]
const D86_KINDS = ("tangent", "pair", "stretching", "frequency", "phase_advance", "birkhoff")
const D86_IDS = Dict("tangent" => "tg", "pair" => "pr", "stretching" => "st", "frequency" => "fq", "phase_advance" => "pa", "birkhoff" => "bk")

# Fixture values are the suite's own: no study value. Every key of a kind is required (I27), so the constructor always writes all of them.
function d86_spec(kind::AbstractString; id = D86_IDS[kind], kw...)
    d = Dict{String, Any}("id" => id, "kind" => kind, "budget" => d86_budget())
    if kind == "tangent"
        d["members"] = ["FLI", "MEGNO", "SALI"]
    elseif kind == "pair"
        merge!(d, Dict{String, Any}("members" => ["RLI", "MEGNO"], "separation" => 1e-8, "frame_source" => "seeded"))
    elseif kind == "stretching"
        merge!(d, Dict{String, Any}("range" => "derived", "bins" => 20, "interval" => 4, "min_per_cell" => 2, "counts" => "S_L"))
    elseif kind == "frequency"
        merge!(d, Dict{String, Any}("observable" => "exp_i_q", "windows" => [16, 32], "order" => 1, "positions" => "half", "pad" => 4,
                                    "mean_removal" => true, "mirror_ratio" => 0.5, "floor_multiplier" => 8.0))
    elseif kind == "phase_advance"
        merge!(d, Dict{String, Any}("observable" => "planar", "windows" => [16, 32], "birkhoff_order" => 2))
    elseif kind == "birkhoff"
        merge!(d, Dict{String, Any}("observables" => d86_obs(), "windows" => [16, 32], "birkhoff_order" => 2))
    else
        error("d86_spec: unknown kind $kind")
    end
    for (k, v) in kw
        v === :delete ? delete!(d, string(k)) : (d[string(k)] = v)
    end
    return d
end

# run length (E22) and the horizon of an arm that uses the spec (R11 item 3)
function d86_runlen(s::AbstractDict)
    s["kind"] in ("frequency", "birkhoff") && return 2 * maximum(s["windows"])
    s["kind"] == "phase_advance" && return 2 * maximum(s["windows"]) + 1
    return nothing
end
d86_window_kind(s) = s["kind"] in ("frequency", "phase_advance", "birkhoff")

# an arm; `methods` absent unless given (R11 item 2: an arm may hold specs only)
function d86_arm(name, orbits, inds, precisions, horizons; methods = nothing, kw...)
    a = Dict{String, Any}("name" => name, "orbits" => orbits, "indicators" => collect(inds), "precisions" => collect(precisions),
                          "horizons" => collect(horizons))
    methods === nothing || (a["methods"] = collect(methods))
    for (k, v) in kw
        a[string(k)] = v
    end
    return a
end

# a run block: the common part of the D53 fixture; the fraction 0.1 stays in the grid: an arm with indicator specs only is exempt from the rule
# "N / 10 is a multiple of tau" (correction 1 of the D86 stage 02 checkpoint), so the run length 2 W + 1 of phase_advance is allowed
function d86_block(arms, specs)
    d = ex_rb_common()
    d["grid"] = Dict{String, Any}("points_per_decade" => 5, "first_steps" => 10, "extra_fractions" => [0.1, 0.5])
    d["arm"] = Any[arms...]
    d["indicator"] = Any[specs...]
    return d
end

# The fixture block with all six kinds (D86c): hand counts in the plan section of experiments-test-plan.md.
function d86_all_block()
    S = Dict(k => d86_spec(k) for k in D86_KINDS)
    prd = d86_spec("pair"; id = "prd", frame_source = "direction", members = ["RLI"])
    specs = [S["tangent"], S["pair"], prd, S["stretching"], S["frequency"], S["phase_advance"], S["birkhoff"]]
    arms = [d86_arm("ind1", "S_D", ["tg", "st", "pr", "prd"], ["Float64"], [100]),
            d86_arm("ind2", "S_D", ["tg", "st", "prd", "pr", "fq", "bk"], ["Float64"], [64]; direction = ["pi/4", "0.3", "0.5"], coordinate_mode = ["wrapped"]),
            d86_arm("ind1dup", "S_D", ["tg"], ["Float64"], [100]),
            d86_arm("ind3", "S_D", ["pa"], ["Float64"], [65]; coordinate_mode = ["wrapped"]),
            d86_arm("mix", "S_L", ["tg"], ["Float64"], [100]; methods = ["L1", "L5"]),
            d86_arm("dbl", "h4", ["tg", "fq", "bk"], ["Double64"], [64]; coordinate_mode = ["wrapped"])]
    return d86_block(arms, specs)
end

# a one-arm block over the spec dictionaries given; horizon of window kinds from E22
function d86_one(specs...; orbits = "S_D", precisions = ["Float64"], horizon = 100, tau = nothing, mode = nothing, name = "a", kw...)
    sp = collect(specs)
    wk = [s for s in sp if d86_window_kind(s)]
    tl = tau === nothing ? 1 : tau
    h = isempty(wk) ? horizon : tl * d86_runlen(first(wk))
    extra = Dict{Symbol, Any}(kw)
    tau === nothing || (extra[:tau] = [tau])
    mode === nothing || (extra[:coordinate_mode] = [mode])
    return d86_block([d86_arm(name, orbits, [s["id"] for s in sp], precisions, [h]; extra...)], sp)
end

d86_fx_cache = Dict{String, Any}()
function d86_fx(tag::AbstractString, run::AbstractDict; panel::Symbol = :std, rev::Bool = false)
    return get!(d86_fx_cache, tag * (rev ? "-rev" : "")) do
        ex_fix_with(run, "d86-" * tag * (rev ? "-rev" : ""); panel = panel, rev = rev)
    end
end
const D86_PLANS = Dict{Any, Any}()                                   # plans are pure and read-only: one per fixture tag
d86_plan(tag, run; kw...) = get!(() -> ex_plan(d86_fx(tag, run; kw...)), D86_PLANS, (tag, get(kw, :panel, :std), get(kw, :rev, false)))

# campaigns shared between testsets, keyed by tag (the campaign itself is the first thing that errors without the D86 code)
const D86_RUNS = Dict{String, Any}()
function d86_run(tag::AbstractString, run::AbstractDict; threads = 4, capacity = 64, order = :plan, panel::Symbol = :std, kw...)
    return get!(D86_RUNS, tag) do
        fx = d86_fx(tag, run; panel = panel)
        store = ex_store("d86-" * tag)
        res = ex_campaign(fx; store = store, threads = threads, capacity = capacity, order = order, kw...)
        (; fx, store, res, plan = ex_plan(fx))
    end
end

# ---------------------------------------------------------------- reading a chunk file without the package reader (spec fixes the layout, R15)
# path => array for datasets, "path@name" => value for attributes ("@name" for root attributes)
function _d86_raw(file::AbstractString)
    H = EX.HDF5
    out = Dict{String, Any}()
    function walk(obj, prefix)
        for a in keys(H.attributes(obj))
            out[prefix * "@" * a] = H.read_attribute(obj, a)
        end
        obj isa H.Dataset && return
        for k in keys(obj)
            x = obj[k]
            p = isempty(prefix) ? k : prefix * "/" * k
            if x isa H.Dataset
                out[p] = read(x)
                walk(x, p)
            else
                walk(x, p)
            end
        end
    end
    H.h5open(f -> walk(f, ""), String(file), "r")
    return out
end

d86_bits(x::AbstractArray{Float64}) = reinterpret(UInt64, vec(x))
d86_bits(x::AbstractArray) = vec(x)
d86_same(a, b) = a isa AbstractArray && b isa AbstractArray ? (size(a) == size(b) && d86_bits(a) == d86_bits(b)) : isequal(a, b)

# canonical bytes of a raw read-back (datasets and attributes), without wall and thread seconds, host and time attributes
# `skip` names root attributes (exact keys, "@name") that differ by design between runs that must agree on everything else
function d86_canon_raw(raw::AbstractDict; skip = ())
    io = IOBuffer()
    for k in sort(collect(keys(raw)))
        k in skip && continue
        any(s -> occursin(s, k), ("@created_utc", "@cpu", "@threads", "@session", "@chunk_size", "@git_sha", "@dirty", "@julia_version", "@packages")) && continue
        v = raw[k]
        if occursin(r"/cost/", k) && v isa AbstractVector
            println(io, k, " ", v[[1, 2, 5]])                       # entries 3, 4 (seconds) and 6 (transform seconds) are not results
        elseif v isa AbstractArray && isbitstype(eltype(v))
            println(io, k, " ", size(v)); write(io, d86_bits(v)); println(io)
        else
            println(io, k, " ", repr(v))
        end
    end
    return take!(io)
end
d86_canon_file(f; skip = ()) = d86_canon_raw(d86_raw(f); skip = skip)
d86_store_sha(store; skip = ()) = ex_sha256(reduce(vcat, [Vector{UInt8}(codeunits(basename(f))); d86_canon_file(f; skip = skip)] for f in ex_chunk_files(store); init = UInt8[]))
# the attributes that differ by design between a device run, a fallback run, a rerun and a CPU run (the backend and its reason, the device, the launch shape, the attempt):
# I24 asks every stored dataset and the other attributes to be identical
const D86_BACKEND_ATTRS = ("@backend", "@backend_fallback_reason", "@device", "@gpu_chunks", "@attempt", "@first_attempt_reasons")

# one orbit of a chunk file: the dataset of an indicator, its attribute
d86_ds(raw, orbit, spec, name) = raw["orbits/$orbit/indicators/$spec/$name"]
d86_at(raw, orbit, spec, ds, name) = raw["orbits/$orbit/indicators/$spec/$ds@$name"]

# every orbit of a campaign: f(item, orbit id, raw, chunk, file)
function d86_each(f, run; methods = nothing)
    for file in ex_chunk_files(run.store)
        c = run.plan.chunks[ex_ordinal(file)]
        methods === nothing || c.run_type.method in methods || continue
        raw = d86_raw(file)
        for i in c.items
            it = run.plan.items[i]
            f(it, it.orbit, raw, c, file)
        end
    end
end

# ---------------------------------------------------------------- independent references written from the spec
# D53 R1 item 2 with the dimensions of R11 item 2 and R12 item 2: key = (orbit, id, precision, horizon, tau, direction or nothing, mode).
# Window kinds carry no direction; stretching and tangent always do; a pair does only with frame_source "direction". A spec is undefined
# for the linear control Y5 when its observable is exp_i_q, cos or sin (E26).
function d86_applicable(spec, orbit)
    orbit["control_kind"] == "y5" || return true                                    # only Y5 is torus-free in the fixture panel
    spec["kind"] == "frequency" && return spec["observable"] != "exp_i_q"
    spec["kind"] == "birkhoff" && return !any(o -> o["kind"] in ("cos", "sin"), spec["observables"])
    return true
end
d86_has_dir(spec) = spec["kind"] in ("tangent", "stretching") || (spec["kind"] == "pair" && spec["frame_source"] == "direction")

function d86_ref_keys(run::AbstractDict, orbits)
    prim = run["primary"]
    nP = maximum(o["id"] for o in orbits)
    sub = Dict("P" => [o["id"] for o in orbits], "S_D" => [o["id"] for o in orbits if o["in_SD"]], "S_L" => [o["id"] for o in orbits if o["in_SL"]],
               "h4" => collect((nP + 1):(nP + length(get(run, "h4", [])))))
    byid = Dict(o["id"] => o for o in orbits)
    specs = Dict(s["id"] => s for s in run["indicator"])
    items = Dict{Tuple, Vector{String}}()
    before = Dict{String, Int}()
    omitted = Dict{String, Vector{Int}}()
    perspec = Dict{String, Vector{Int}}()                                           # id => [unique items, steps]
    for arm in run["arm"], prec in arm["precisions"], N in arm["horizons"], mode in get(arm, "coordinate_mode", [prim["coordinate_mode"]]),
        sid in get(arm, "indicators", String[]), id in sub[arm["orbits"]]
        s = specs[sid]
        if id <= nP && !d86_applicable(s, byid[id])
            id in get!(omitted, sid, Int[]) || push!(omitted[sid], id)
            continue
        end
        dirs = d86_has_dir(s) ? get(arm, "direction", [prim["direction"]]) : Any[nothing]
        for tau in get(arm, "tau", [prim["tau"]]), dir in dirs
            key = (id, sid, prec, N, tau, dir, mode)
            push!(get!(items, key, String[]), arm["name"])
            before[arm["name"]] = get(before, arm["name"], 0) + 1
        end
    end
    for (k, arms) in items
        v = get!(perspec, k[2], [0, 0])
        v[1] += 1; v[2] += k[4]
    end
    return (; items, before, omitted, perspec)
end

# canonical text of a spec without its budget (E27; root path `indicator`, 1-based list index: correction 3 of the D86 stage 02 checkpoint): `indicator.<key> = <value>`, keys sorted, lists in order, tables of a list as `key[i].name`
function d86_spec_text(spec::AbstractDict)
    io = IOBuffer()
    cv(x::AbstractString) = repr(String(x))
    cv(x::AbstractVector) = "[" * join((cv(e) for e in x), ", ") * "]"
    cv(x) = repr(x)
    function rec(x, path)
        if x isa AbstractDict
            foreach(k -> rec(x[k], path * "." * k), sort!(collect(keys(x))))
        elseif x isa AbstractVector && !isempty(x) && first(x) isa AbstractDict
            foreach(((i, e),) -> rec(e, "$path[$i]"), enumerate(x))
        else
            println(io, path, " = ", cv(x))
        end
    end
    rec(Dict{String, Any}(k => v for (k, v) in spec if k != "budget"), "indicator")
    return String(take!(io))
end
d86_spec_hash(spec) = (t = d86_spec_text(spec); d86_memo(() -> ex_sha256(t), (:sh, t)))

# item text of a new kind (E27): the D53 text without `eps`, with `direction` only when present, then kind and spec lines
function d86_item_text(it, spec_hash::AbstractString)
    io = IOBuffer()
    for f in (:orbit, :system, :param, :x0, :p0, :run_seed, :orbit_seed, :method, :precision, :horizon, :tau, :direction, :coordinate_mode,
              :vectors, :grid_sha256, :code_version, :hash_version)
        v = getproperty(it, f)
        f === :direction && v === nothing && continue
        print(io, f, "=", v, "\n")
    end
    print(io, "kind=", it.kind, "\n", "spec=", spec_hash, "\n")
    return String(take!(io))
end
d86_chunk_key(it, spec_hash) = "method=$(it.method) kind=$(it.kind) spec=$(spec_hash[1:16]) precision=$(it.precision) horizon=$(it.horizon) tau=$(it.tau) " *
                               "eps=- direction=$(it.direction === nothing ? "-" : it.direction) mode=$(it.coordinate_mode) cell=$(it.cell)"

# ---------------------------------------------------------------- stand-alone references (O16): the verified reducer run directly on the identical orbit
const D86_MEMBERS = Dict("SALI" => Chaoslab.SALI, "FLI" => Chaoslab.FLI, "MEGNO" => Chaoslab.MEGNO, "RLI" => Chaoslab.RLI)
d86_red(m::AbstractString; replace::Bool = true) = m == "SALI" ? Chaoslab.SALI() : D86_MEMBERS[m](; replace = replace)

# forms of a chain (E23): every member but the last non-replacing, the last replacing
d86_forms(members) = [i == length(members) for i in eachindex(members)]
d86_chain_record(members) = ["$m(replace=$(f))" for (m, f) in zip(members, d86_forms(members))]

# names and values of the columns of a member at the final of its stand-alone run (E24)
function d86_member_cols(m::AbstractString, fin)
    m == "SALI" && return ["SALI.sali" => Float64(fin.sali), "SALI.L" => Float64(fin.L)]
    m == "FLI" && return ["FLI.fli" => Float64(fin.fli), "FLI.fli1" => Float64(fin.fli1)]
    m == "MEGNO" && return reduce(vcat, [["MEGNO.Y$c" => Float64(fin.Y[c]), "MEGNO.Ybar$c" => Float64(fin.Ybar[c]), "MEGNO.Ybar20_$c" => Float64(fin.Ybar20[c])] for c in 1:2])
    m == "RLI" && return ["RLI.DL" => Float64(fin.DL), "RLI.RLI" => Float64(fin.RLI)]
    m == "stretching" && return ["stretching.D" => Float64(fin.D)]
    error("d86_member_cols: $m")
end
d86_col_names(members) = reduce(vcat, [first.(d86_member_cols(m, d86_dummy_final(m))) for m in members])
function d86_dummy_final(m)
    z = zero(Float64)
    return (sali = z, L = z, fli = z, fli1 = z, Y = (z, z), Ybar = (z, z), Ybar20 = (z, z), DL = z, RLI = z, D = z)
end

# the identical orbit of a plan item in its element type (system, start, frame by the module's own frame rule R2 item 6)
function d86_orbit_args(it, T)
    sys, x = ex_sys_x(it, T)
    W = EX.deviation_frame(T, it.direction, 2)
    return sys, x, SMatrix{2, 2, T}(W)
end

# final of a chain run directly to step n; `members` in order, only the member `which` is read; forms as the module sets them, or all replacing
# (the stand-alone replacing form of I22: the member alone in a chain)
function _d86_stand_alone(it, which::AbstractString, n::Integer)
    return ex_with_prec(it.precision) do T
        sys, x, W = d86_orbit_args(it, T)
        r = d86_red(which)
        _, _, st = ex_Tg.propagate_orbit(sys, x, W, n, it.tau, (r,))
        ex_Tg.reducer_final(r, st[1])
    end
end

# StretchingNumbers for the item's system: range derived by the verified helper, parameters from the spec dictionary
function d86_stretch_red(sys, spec; T = Float64)
    rng = spec["range"] == "derived" ? Chaoslab.stretching_range(sys) : Tuple(spec["range"])
    return Chaoslab.StretchingNumbers(; range = rng, bins = spec["bins"], interval = spec["interval"], min_per_cell = spec["min_per_cell"])
end
function _d86_stretch_final(it, spec, n::Integer)
    return ex_with_prec(it.precision) do T
        sys, x, W = d86_orbit_args(it, T)
        r = d86_stretch_red(sys, spec; T = T)
        ws = zeros(UInt32, ex_Tg.reducer_workspace(r, T, Val(2), Val(2)))
        _, _, st = ex_Tg.propagate_orbit(sys, x, W, n, it.tau, (r,); workspace = ws)
        (fin = ex_Tg.reducer_final(r, st[1], ws), range = (r.lo, r.hi), bins = r.bins)
    end
end

# series of a window kind built by plain iteration of the map (O17): `forward` on a vector of the element type and the sample convention
# of E21 (sample j is the observable at the point after (j + 1) tau steps); no code of the recorder
function _d86_complex_series(it, nsamp::Integer, observable::AbstractString)
    return ex_with_prec(it.precision) do T
        pts = d86_points_T(it, nsamp, T)
        observable == "exp_i_q" ? Complex{T}[cis(p[1]) for p in pts] : Complex{T}[Complex(p[1], p[2]) for p in pts]
    end
end
function d86_points_T(it, nsamp::Integer, ::Type{T}) where {T}
    sys, x = ex_sys_x(it, T)
    out = Vector{typeof(x)}()
    for _ in 1:nsamp
        for _ in 1:it.tau
            x = Chaoslab.forward(sys, x)
        end
        push!(out, x)
    end
    return out
end

# the orbit points of an item at the grid steps, by plain iteration (for I21 checks) are `ex_iter_states`; here the L1 twin of a plan item
function d86_l1_twin(it)
    return merge(it, (method = "L1", vectors = 1, direction = it.direction === nothing ? "pi/4" : it.direction))
end

# flag words: bit i (1-based, in the order of the `bits` attribute) is 2^(i-1)
d86_flagset(word::Integer, names) = Set(String(names[i]) for i in eachindex(names) if (word >> (i - 1)) & 1 == 1)

# StableRNGs-based seeded starts for the shear tests: x0 and p0 in [0, 2 pi)
d86_rand_starts(sub::Integer, n::Integer) = (rng = StableRNG(d86_seed('z', sub)); [(2pi * rand(rng), 2pi * rand(rng)) for _ in 1:n])

# the same block without the arm that holds methods (indicator arms only), shared by several testsets
function d86_all_ind_block()
    run = d86_all_block()
    run["arm"] = [a for a in run["arm"] if !haskey(a, "methods")]
    return run
end

# true when a store holds no chunk file and no temporary file (the directory may not exist)
d86_nofiles(store) = !isdir(store) || (isempty(ex_chunk_files(store)) && isempty(ex_tmp_files(store)))
d86_files(store) = isdir(store) ? ex_chunk_files(store) : String[]

# the panel entry of an orbit id, `nothing` for the extra starts (ids after the panel)
d86_pinfo(plan, o) = o <= length(plan.panel.orbits) ? plan.panel.orbits[o] : nothing
d86_spec_of(plan, it) = only(s for s in plan.block.raw["indicator"] if s["id"] == it.method)

# the stand-alone pair run on the identical orbit (O16): xi0 as stored, the block frame (xi_0, 0), (0, xi_0), the member in its replacing form
_d86_pair_ref(it, spec, xi0, n, members) = ex_with_prec(it.precision) do T
    sys, _ = ex_sys_x(it, T)
    x = SVector(exact(T, it.x0), exact(T, it.p0))              # pair_start takes the lifted start and forms the wrapped state itself
    u = SVector{2, T}(xi0...)
    ps = Chaoslab.pair_start(sys, x, u, spec["separation"]; orbit = it.orbit)
    W = SMatrix{4, 2, T}(u[1], u[2], zero(T), zero(T), zero(T), zero(T), u[1], u[2])
    red = d86_red(members; replace = true)
    _, _, st = ex_Tg.propagate_orbit(Chaoslab.PairMap(sys), ps.z, W, n, it.tau, (red,))
    (fin = ex_Tg.reducer_final(red, st[1]), start = ps)
end


# observables of a birkhoff spec as the verified callables (E26); lifted mode has no winding count, the lifted coordinate is the state entry itself
function d86_observables(spec, mode)
    return map(spec["observables"]) do o
        c = o["coordinate"]
        o["kind"] == "cos" ? Chaoslab.CosObservable(c, o["frequency"]) : o["kind"] == "sin" ? Chaoslab.SinObservable(c, o["frequency"]) :
        mode == "wrapped" ? Chaoslab.LiftedObservable(c, c + 2) : (x -> x[c])
    end
end

function _d86_ref_freq(it, spec)
    Ws = spec["windows"]; Wm = maximum(Ws)
    series = d86_complex_series(it, 2Wm, spec["observable"])
    return map(Ws) do W
        r = Chaoslab.frequency_map(series[1:(2W)], W; order = spec["order"], positions = Symbol(spec["positions"]), pad = spec["pad"],
                                   mean_removal = spec["mean_removal"], mirror_ratio = spec["mirror_ratio"], floor_multiplier = spec["floor_multiplier"])
        w1, w2 = r.w1[1], r.w2[1]
        vals = Float64[w1.nu, w2.nu, w1.amplitude, w2.amplitude, w1.mirror_amplitude, w2.mirror_amplitude, w1.second_ratio, w2.second_ratio, r.dist[1], r.delta, r.D]
        fl = Set{String}()
        r.below_floor && push!(fl, "below_floor")
        r.status[1] === :degenerate && push!(fl, "degenerate")
        r.status[1] === :failed && push!(fl, "failed")
        w1.tie && push!(fl, "tie1"); w2.tie && push!(fl, "tie2")
        w1.near_zero && push!(fl, "near_zero1"); w2.near_zero && push!(fl, "near_zero2")
        w1.mirror && push!(fl, "mirror1"); w2.mirror && push!(fl, "mirror2")
        (; vals, flags = fl, status = r.status[1])
    end
end
const D86_FREQ_COLS = ["nu1", "nu2", "amp1", "amp2", "mirror1", "mirror2", "second1", "second2", "dist", "delta", "D"]
const D86_FLAGS_CHECKED = ("below_floor", "degenerate", "failed", "tie1", "tie2", "near_zero1", "near_zero2", "mirror1", "mirror2")

# the E2 rows from the pure functions on the series built by plain iteration (O17)
function _d86_ref_pa(it, spec)
    Ws = spec["windows"]; Wm = maximum(Ws)
    z = spec["observable"] == "planar" ? d86_complex_series(it, 2Wm + 1, "q_plus_i_p") : nothing
    return map(Ws) do W
        if z === nothing                                      # lifted: increments of the lifted first coordinate rebuilt from the winding count / the lifted state
            pts = ex_with_prec(it.precision) do T
                d86_points_T(it, 2Wm + 1, T)
            end
            q = [it.coordinate_mode == "wrapped" ? p[1] + 2 * Chaoslab.Systems._pi(p[1]) * p[3] : p[1] for p in pts]
            inc = diff(q)
        else
            inc = Chaoslab.phase_increments(z)
        end
        nu = Chaoslab.weighted_phase_frequency(inc, W; p = spec["birkhoff_order"])
        delta = Chaoslab.circle_distance(nu[1], nu[2])
        fi = Chaoslab.frequency_index(delta)
        (; vals = Float64[nu[1], nu[2], delta, fi.D], below = fi.below_floor)
    end
end

function _d86_ref_bk(it, spec)
    Ws = spec["windows"]; Wm = maximum(Ws)
    obs = d86_observables(spec, it.coordinate_mode)
    return ex_with_prec(it.precision) do T
        pts = d86_points_T(it, 2Wm, T)
        series = Matrix{T}(undef, 2Wm, length(obs))
        for (j, p) in enumerate(pts), (k, o) in enumerate(obs)
            series[j, k] = T(o(p))
        end
        r = Chaoslab.weighted_birkhoff(series, Tuple(Ws); order = spec["birkhoff_order"])
        rows = Vector{Float64}[]; zeros_ = Bool[]; floors = Bool[]
        for (h, W) in enumerate(Ws), k in eachindex(obs)
            ix = Chaoslab.birkhoff_index(r.avg1[k, h], r.avg2[k, h], r.S[k, h])
            push!(rows, Float64[r.avg1[k, h], r.avg2[k, h], r.S[k, h], r.B[h], ix.delta, ix.absdig, ix.reldig, ix.dig])
            push!(zeros_, ix.zero); push!(floors, ix.below_floor)
        end
        (; rows, zero = zeros_, below = floors)
    end
end

d86_isequal(a::AbstractVector, b::AbstractVector) = length(a) == length(b) && all(isequal.(a, b))

# the item table of an orbit: rows x columns of `final`, the flag words
d86_final(raw, orbit, id) = (d86_ds(raw, orbit, id, "final"), d86_ds(raw, orbit, id, "flags"), d86_at(raw, orbit, id, "flags", "bits"))


# ---------------------------------------------------------------- bytes(orbit) of E28 written from the spec
const D86_E24_COLS = Dict("SALI" => 2, "FLI" => 2, "MEGNO" => 6, "RLI" => 2, "stretching" => 1)
function d86_ref_bytes(plan, it)
    spec = only(s for s in plan.block.raw["indicator"] if s["id"] == it.method)
    kind = spec["kind"]
    T = it.precision == "Double64" ? Chaoslab.Double64 : Float64
    wrapped = it.coordinate_mode == "wrapped"
    s = (wrapped && it.system != "linear") ? 4 : 2
    pair = kind == "pair"
    d, k, sK = pair ? 4 : 2, it.vectors, pair ? 2s : s
    M = ex_Tg.bytes_per_orbit(Val(d), Val(k), T, sK)
    rb(r) = ex_Tg.reducer_bytes(r, T, Val(d), Val(k))
    mon = 3 * sizeof(T)
    mem = Dict("SALI" => Chaoslab.SALI(), "FLI" => Chaoslab.FLI(replace = false), "MEGNO" => Chaoslab.MEGNO(replace = false), "RLI" => Chaoslab.RLI(replace = false))
    R, cols = if kind in ("tangent", "pair")
        sum(rb(mem[m]) for m in spec["members"]) + mon, sum(D86_E24_COLS[m] for m in spec["members"])
    elseif kind == "stretching"
        sys = it.system == "linear" ? Chaoslab.y5_control(parse(Float64, it.param)).sys : Chaoslab.StandardMap(1.0)
        rb(d86_stretch_red(sys, spec)) + mon, 1
    elseif kind == "frequency"
        mon, 0
    elseif kind == "phase_advance"
        rb(Chaoslab.WeightedPhaseAdvance(window = Tuple(spec["windows"]), observable = Chaoslab.PlanarObservable(1, 2), birkhoff_order = spec["birkhoff_order"])) + mon, 0
    else
        rb(Chaoslab.WeightedBirkhoff(window = Tuple(spec["windows"]), observables = Tuple(d86_observables(spec, it.coordinate_mode)), birkhoff_order = spec["birkhoff_order"])) + mon, 0
    end
    gr = plan.block.raw["grid"]
    ng = length(ex_ref_grid(it.horizon, it.tau; ppd = gr["points_per_decade"], first_steps = gr["first_steps"], fracs = gr["extra_fractions"]))
    c = 1 + 2 + (wrapped ? 2 : 0) + 2 + cols
    S = 8 * c * ng
    B = kind == "frequency" ? (T === Float64 ? 32 : 64) * maximum(spec["windows"]) : 0
    return (; M, R, S, B, total = M + R + S + B, ng, c)
end
# the launch capacity of the plan (R17 item 1): the most demanding chunk bounds it; f = 1/2
d86_ref_capacity(plan, free) = max(1, minimum(Int(fld(big(free), 2 * big(d86_ref_bytes(plan, plan.items[first(c.items)]).total))) for c in plan.chunks))


# ---------------------------------------------------------------- memoized references and reads (the same reference run is used by several testsets and files)
# A reference depends on the orbit (system, parameter, start, precision, mode, tau, direction), the member or the spec text and the step: nothing else, so the key is that.
const D86_MEMO = Dict{Any, Any}()
const D86_RAW = Dict{Any, Any}()
d86_okey(it) = (it.system, it.param, it.x0, it.p0, it.precision, it.coordinate_mode, it.tau, it.direction)
d86_memo(f, key) = get!(f, D86_MEMO, key)
d86_stand_alone(it, which::AbstractString, n::Integer) = d86_memo(() -> _d86_stand_alone(it, which, n), (:sa, d86_okey(it), which, Int(n)))
d86_stretch_final(it, spec, n::Integer) = d86_memo(() -> _d86_stretch_final(it, spec, n), (:st, d86_okey(it), d86_spec_text(spec), Int(n)))
d86_pair_ref(it, spec, xi0, n, members) = d86_memo(() -> _d86_pair_ref(it, spec, xi0, n, members), (:pr, d86_okey(it), d86_spec_text(spec), Tuple(xi0), Int(n), members))
d86_ref_freq(it, spec) = d86_memo(() -> _d86_ref_freq(it, spec), (:fq, d86_okey(it), d86_spec_text(spec)))
d86_ref_pa(it, spec) = d86_memo(() -> _d86_ref_pa(it, spec), (:pa, d86_okey(it), d86_spec_text(spec)))
d86_ref_bk(it, spec) = d86_memo(() -> _d86_ref_bk(it, spec), (:bk, d86_okey(it), d86_spec_text(spec)))
d86_complex_series(it, nsamp::Integer, observable::AbstractString) = d86_memo(() -> _d86_complex_series(it, nsamp, observable), (:cs, d86_okey(it), Int(nsamp), observable))
# a chunk file is read once per (path, modification time, size); the files of a store are never changed after the campaign
d86_raw(file::AbstractString) = get!(() -> _d86_raw(file), D86_RAW, (String(file), mtime(file), filesize(file)))

# the all-kinds block without the arm that holds methods, run once and shared (D86q, D86t, D86v, D86w)
d86_all_ind_run() = d86_run("all-ind", d86_all_ind_block(); threads = 4, capacity = 10^6)
