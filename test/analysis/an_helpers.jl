# Helpers for the analysis suite (prefix an_). Definitions only: nothing here calls Chaoslab.Analysis at load time.
# Spec: workspaces/module-build/stages/01-spec/output/analysis-spec.md. API under test: the test plan, section "API under test".
using Test, Chaoslab, StableRNGs, TOML

# AN.name resolves Chaoslab.Analysis.name when called, so a missing module is an UndefVarError inside the testset.
struct AnLazy end
Base.getproperty(::AnLazy, s::Symbol) = getfield(Chaoslab.Analysis, s)
const AN = AnLazy()
const AX = Chaoslab.Experiments        # store writer and reader, used by the fixture only

const AN_SEED = 20261012               # _config/numerics-defaults.toml, seed
const AN_BASE = AN_SEED + 28000        # this slug uses 20261012 + 28000 .. 28999; testset AC<n> uses 50 n + sub (sub < 50)
an_seed(ac::Integer, sub::Integer = 0) = AN_BASE + 50 * ac + sub

const AN_DIR = mktempdir()
an_path(name) = joinpath(AN_DIR, name)

# ---- independent SHA-256 (platform command-line tool); the module's value is compared against it
function an_sha256(bytes::AbstractVector{UInt8})
    mktemp() do path, io
        write(io, bytes)
        close(io)
        for cmd in (`sha256sum $path`, `certutil -hashfile $path SHA256`, `shasum -a 256 $path`)
            out = try
                read(cmd, String)
            catch
                continue
            end
            m = match(r"\b[0-9a-fA-F]{64}\b", out)
            m === nothing || return lowercase(m.match)
        end
        error("no command-line SHA-256 tool found")
    end
end
an_sha256(s::AbstractString) = an_sha256(Vector{UInt8}(codeunits(s)))

# ---- parameter table E0 (flat keys). Thresholds are the frozen study values; horizons are scaled down in the pipeline fixture.
const AN_PARAMS = (; N_1 = 10_000_000, N_long = 100_000_000, N_B = 1_000_000, n_H2 = 10_000,
    h1a_bound = 1e-8, h1a_fraction = 0.98, h1b_bound = 1e-6, h1b_fraction = 0.98, h1c_bound = 5e-3, h1c_fraction = 0.90, h1c_eps = 1e-8,
    control_bound = 1e-8, control_bound_L9 = 5e-8,
    h2_median = 0.10, h2_interval = 0.05, bootstrap_resamples = 10_000, bootstrap_level = 0.95,
    h3_band = 0.05, h3_ratio = 10.0, h3_interval = 3.0,
    h4_ref_halfpi = 0.298, h4_tol_halfpi = 0.01, h4_ref_2pi = 1.172, h4_tol_2pi = 0.02, h4_T = 1_000_000_000, h4_T_cross = 100_000_000,
    h5a_min = 0.90, h5b_max = 0.30, h5c_max = 0.30,
    lambda_threshold = 1e-4, growth_ratio = 5.0, rotation_tol = 1e-7, return_tol = 1e-10, ladder_tol = 5e-3,
    drift_bound = 1e6, min_set = 20, K_c = 0.971635406)
const AN_PARAM_KEYS = collect(keys(AN_PARAMS))
an_params(; kw...) = merge(AN_PARAMS, NamedTuple(kw))
# the pipeline fixture: horizons scaled down, regular floor of the scaled horizon below the instability threshold, small sets allowed
an_fx_params(; kw...) = an_params(; N_1 = 10_000, N_long = 100_000, N_B = 1_000, n_H2 = 100, h4_T = 5_000, h4_T_cross = 1_000,
    lambda_threshold = 1e-2, min_set = 3, bootstrap_resamples = 1000, kw...)
function an_write_params(path::AbstractString, p::NamedTuple)
    open(io -> TOML.print(io, Dict{String, Any}(string(k) => v for (k, v) in pairs(p))), path, "w")
    return path
end

# ---- closed forms
an_k0(n; a = pi / 4) = log(sqrt((cos(a) + n * sin(a))^2 + sin(a)^2)) / n      # regular orbit, preregistration A1 item 2
an_grid(N::Integer; ppd = 10, first_steps = 10, fracs = (0.1, 0.5)) = begin
    s = Set{Int}(1:min(first_steps, N))
    j = 0
    while 10.0^(j / ppd) <= N
        push!(s, Int(floor(10.0^(j / ppd) + 0.5)))
        j += 1
    end
    foreach(f -> push!(s, Int(floor(f * N + 0.5))), fracs)
    push!(s, N)
    sort!(collect(s))
end

# ---- fixture panel. P = 20 orbits (ids 1 to 20), h4 starts 21 and 22 (extra starts, ids after the panel).
const AN_K1, AN_K2 = 0.5, 2 * pi
const AN_V = ["L1", "L3", "L4", "L4b", "L5", "L5b", "L6", "L9"]
function an_orbits()
    os = Dict{String, Any}[]
    function add(id, K, kidx, class; control = false, kind = "none", exponent = nothing, stability = nothing, sd = false, sl = false,
                 system = "standard", x0 = 0.5, p0 = 0.7, period = 1)
        d = Dict{String, Any}("id" => id, "K" => K, "K_text" => string(K), "class" => class, "control" => control, "control_kind" => kind,
            "in_SD" => sd, "in_SL" => sl, "k_index" => kidx, "system" => system, "x0" => x0, "p0" => p0, "x0_dec" => string(x0),
            "p0_dec" => string(p0), "run_seed" => "20261012", "orbit_seed" => string(1000 + id), "draw_index" => id)
        if control
            d["exponent"], d["period"] = exponent, period
            stability === nothing || (d["stability"] = stability)
        end
        system == "linear" && (d["mu"] = "2.618033988749895"; delete!(d, "K"); delete!(d, "K_text"))       # the real Y5 control has no K and no K_text
        push!(os, d)
    end
    for i in 1:4;  add(i, AN_K1, 1, "C"; sd = i <= 2, sl = i == 1, p0 = 0.1 * i); end
    for i in 5:6;  add(i, AN_K1, 1, "S"; sd = i == 5, sl = i == 5, p0 = 0.1 * i); end
    for i in 7:8;  add(i, AN_K1, 1, "X"; sd = i == 7, p0 = 0.1 * i); end
    for i in 9:10; add(i, AN_K1, 1, "T"; sd = i == 9, sl = i == 9, p0 = 0.1 * i); end
    for i in 11:14; add(i, AN_K2, 2, "C"; sd = i <= 12, sl = i == 11, p0 = 0.1 * (i - 10)); end
    for i in 15:16; add(i, AN_K2, 2, "T"; sd = i == 15, p0 = 0.1 * (i - 14)); end
    add(17, AN_K1, 1, "control"; control = true, kind = "fixed", exponent = 0.69, stability = "hyperbolic", sd = true, x0 = 0.0, p0 = 0.0)
    add(18, AN_K1, 1, "control"; control = true, kind = "newton", exponent = 0.5, stability = "hyperbolic", sd = true, x0 = 1.0, p0 = 0.2, period = 2)
    add(19, AN_K1, 1, "control"; control = true, kind = "fixed", exponent = 0.0, stability = "elliptic", sd = true, x0 = Float64(pi), p0 = 0.0)
    add(20, 0.0, 0, "control"; control = true, kind = "y5", exponent = 0.96, stability = "hyperbolic", sd = true, system = "linear", x0 = 0.3, p0 = 0.1)
    return os
end
const AN_ORBITS = an_orbits()
an_orbit(id) = AN_ORBITS[id]
const AN_H4 = [(id = 21, K = pi / 2, x0 = 1e-3, p0 = -1e-3), (id = 22, K = 2 * pi, x0 = 1e-3, p0 = -1e-3)]

# profile of an orbit: the closed-form series lambda_m(n) = lref (1 + c_m / n) for the chaotic classes, the integrable closed form for T and the
# elliptic control. c_m defaults: 0 for L1 and the QR family, 0.02 for L9, 0.001 for L2.
function an_profile(id::Integer, cover::AbstractDict = Dict())
    o = an_orbit(id)
    kind = (o["class"] == "T" || id == 19) ? :regular : :const
    lref = o["control"] ? (id == 18 ? 0.8 : Float64(o["exponent"])) : (o["class"] == "C" ? 0.30 + 0.01 * id : 0.02 + 0.001 * id)
    c = Dict{String, Float64}(m => 0.0 for m in AN_V)
    c["L9"], c["L2"] = (o["control"] ? 0.0 : 0.02), (o["control"] ? 0.0 : 0.001)
    q = Dict{String, Float64}(m => 1.0 for m in vcat(AN_V, "L2"))
    for (m, v) in get(cover, id, Dict())                 # a number c (decay 1/n) or a pair (c, q) (decay 1/n^q)
        c[m], q[m] = v isa Tuple ? (Float64(v[1]), Float64(v[2])) : (Float64(v), 1.0)
    end
    return (kind = kind, lref = lref, c = c, q = q)
end
an_vals(prof, m::AbstractString, steps) = prof.kind === :regular ? an_k0.(Float64.(steps)) : prof.lref .* (1 .+ prof.c[m] ./ Float64.(steps) .^ prof.q[m])

# variant offsets of the sweep items (so that F7 numbers are known): tau, eps, direction, coordinate mode
function an_offset(spec)
    off = 1e-4 * (spec.tau - 1)
    spec.eps === nothing || (off += 1e-3 * log10(spec.eps / 1e-8))
    off += 1e-5 * (findfirst(==(spec.direction), ["pi/4", "0.3", "0.3 + pi/5", "0.3 + 2pi/5", "0.3 + 3pi/5"]) - 1)
    spec.mode == "wrapped" && spec.method != "L2" && (off += 5e-5)
    return off
end

# ---- the item list of the fixture (one chunk per method and parameter set; arms merge identical items as the runner does)
const AN_DIRS = ["0.3", "0.3 + pi/5", "0.3 + 2pi/5", "0.3 + 3pi/5"]
function an_specs(p)
    S = NamedTuple[]
    allP, sd, sl = collect(1:20), [o["id"] for o in AN_ORBITS if o["in_SD"]], [o["id"] for o in AN_ORBITS if o["in_SL"]]
    add(m, prec, hor, orbits, arms; tau = 1, eps = nothing, dir = "pi/4", mode = "lifted") =
        push!(S, (method = m, precision = prec, horizon = hor, tau = tau, eps = eps, direction = dir, mode = mode, orbits = orbits, arms = arms))
    for m in AN_V
        arms = ["main"]
        m in ("L1", "L5", "L9") && push!(arms, "directions")
        m in ("L4", "L5") && push!(arms, "tau_sweep")
        add(m, "Float64", p.N_1, allP, arms)
    end
    add("L2", "Float64", p.N_1, allP, ["main_l2", "wrapped", "eps_sweep"]; eps = 1e-8, mode = "wrapped")
    for m in AN_V; add(m, "Float64", p.N_long, sd, ["long"]); end
    add("L2", "Float64", p.N_long, sd, ["long_l2"]; eps = 1e-8, mode = "wrapped")
    for m in ("L4", "L5"), tau in (2, 4); add(m, "Float64", p.N_1, sd, ["tau_sweep"]; tau = tau); end
    for e in (1e-6, 1e-10); add("L2", "Float64", p.N_1, sd, ["eps_sweep"]; eps = e, mode = "wrapped"); end
    for m in ("L1", "L5", "L9"), d in AN_DIRS; add(m, "Float64", p.N_1, sd, ["directions"]; dir = d); end
    for m in AN_V; add(m, "Float64", p.N_1, sd, ["wrapped"]; mode = "wrapped"); end
    for m in ("L1", "L3", "L5", "L9"); add(m, "Double64", p.N_1, sl, ["ladder_double64_short"]); end
    add("L5", "Double64", p.N_long, sl, ["ladder_double64_long"])
    for m in ("L1", "L5", "L9"); add(m, "BigFloat(256)", p.N_B, sl, ["ladder_bigfloat"]); end
    for m in ("L1", "L5"); add(m, "Float64", p.h4_T, [21, 22], ["h4_float64"]); end
    add("L5", "Double64", p.h4_T_cross, [21, 22], ["h4_double64"])
    return S
end

# ---- one orbit record of one item
function an_state(id, steps, mode, hook)
    o = id <= 20 ? an_orbit(id) : nothing
    n = length(steps)
    x0, p0 = o === nothing ? (AN_H4[id - 20].x0, AN_H4[id - 20].p0) : (o["x0"], o["p0"])
    x = x0 .+ steps .* p0
    o !== nothing && id == 18 && (x = x0 .+ 1e-6 .* steps)         # the info control: its binary64 start leaves the orbit
    o !== nothing && o["control"] && o["control_kind"] == "fixed" && (x = fill(x0, n))
    st = hcat(x, fill(p0, n))
    return hook === nothing ? st : hook(id, st)
end
const AN_DIAG = Dict("L1" => ["det", "symp"], "L9" => ["det", "symp", "sum", "pair"])
an_diag(m) = m == "L2" ? String[] : get(AN_DIAG, m, ["det", "symp", "orth", "sum", "pair"])

function an_record(id, spec, steps; vals, flags, attrs, state, state3 = nothing)
    m = spec.method
    k = m in ("L1", "L2") ? 1 : 2
    lam = k == 1 ? reshape(vals, :, 1) : hcat(vals, -vals)
    st = an_state(id, steps, spec.mode, state)
    state3 === nothing || (st = state3(id, st, steps))
    n = length(steps)
    # a breakdown item (Amendment A5 item 2) stores t = 0 and NaN from the first row that was not computed
    blank = "breakdown_tau" in flags ? .!isfinite.(vals) : falses(n)
    t = Int64.(steps)
    t[blank] .= 0
    st[blank, :] .= NaN
    a = Dict{String, Any}("system" => id <= 20 && an_orbit(id)["system"] == "linear" ? "linear" : "standard", "params" => id <= 20 ? get(an_orbit(id), "K_text", "0.0") : string(AN_H4[id - 20].K),
        "ic" => "0,0", "class_label" => id <= 20 ? an_orbit(id)["class"] : "h4", "seed" => string(1000 + id), "item_hash" => "i"^64,
        "design_class" => id <= 20 ? an_orbit(id)["class"] : "h4", "flags" => String.(flags), "steps_recorded" => n)
    merge!(a, Dict{String, Any}(attrs))
    wind = spec.mode == "wrapped" ? zeros(Int64, n, 2) : nothing
    diag = Dict{String, Vector{Float64}}(nm => fill(1e-16, n) for nm in an_diag(m))
    thr = Dict{String, Float64}(nm => 1e-12 for nm in an_diag(m))
    fe = Dict{String, Int}(nm => 0 for nm in an_diag(m))
    return AX.OrbitRecord(id, a, t, spec.mode == "wrapped" ? mod2pi.(st) : st, wind, Dict(m => lam), diag, thr, fe,
                          Dict(m => Float64[spec.horizon, k, 0.01, 0.1, 1e5]), "")
end

# ---- fixture study. Hooks (all optional):
#   lam(id, spec, steps, vals) -> vals               replace the exponent series of one orbit and item
#   flags(id, spec) -> Vector{String}                stored orbit flags of that item
#   attrs(id, spec) -> Dict                          extra orbit attributes (control_scope and so on)
#   state(id, st) -> st                              replace the stored points (n x 2) of an orbit
#   cover = Dict(id => Dict(method => c))            override the profile constants
const AN_STUDIES = Dict{String, Any}()
function an_study(tag::AbstractString; p = an_fx_params(), lam = nothing, flags = nothing, attrs = nothing, state = nothing, state3 = nothing, cover = Dict(),
                  status = "complete", pilot = false)
    dir = an_path("study-$tag")
    rm(dir; recursive = true, force = true)
    store = joinpath(dir, "store")
    mkpath(store)
    # preregistration with the SELF hash rule
    prep(t) = replace(replace(t, "\r\n" => "\n"), r"(?m)^(\| SHA-256 of this file \(current version\) \|)[^\n|]*\|" => s"\1 SELF |")
    text = "# Preregistration: an-fixture\n\n## Header\n\n| Field | Value |\n|-------|-------|\n| Study slug | an-fixture |\n" *
           "| SHA-256 of this file (current version) | " * "0"^64 * " (rule: fixture) |\n| Amendment log | A1 (2026-01-01, fixture) |\n\nFixture. Not a study.\n"
    h = an_sha256(prep(text))
    text = replace(text, "0"^64 => h; count = 1)
    prereg = joinpath(dir, "an-fixture-preregistration.md")
    write(prereg, text)
    # panel file
    panel = joinpath(dir, "an-fixture-panel.toml")
    cells = [Dict{String, Any}("class" => c, "drawn" => n, "k_index" => k, "qualifying" => 100, "reason" => "", "status" => "populated")
             for (k, c, n) in ((1, "C", 4), (1, "S", 2), (1, "X", 2), (1, "T", 2), (2, "C", 4), (2, "T", 2))]
    open(panel, "w") do io
        TOML.print(io, Dict{String, Any}("manifest" => Dict{String, Any}("K_binary64" => [AN_K1, AN_K2], "controls" => [17, 18, 19, 20],
            "cells" => cells, "format_version" => 1), "orbit" => AN_ORBITS))
    end
    ph = an_sha256(read(panel))
    params_file = an_write_params(joinpath(dir, "an-fixture-parameters.toml"), p)
    specs = an_specs(p)
    chunks = Dict{String, Any}[]
    for (ord, sp) in enumerate(specs)
        steps = an_grid(sp.horizon)
        recs = AX.OrbitRecord[]
        for id in sp.orbits
            prof = id <= 20 ? an_profile(id, cover) : (kind = :const, lref = id == 21 ? 0.298 : 1.172, c = Dict{String, Float64}(m => 0.0 for m in vcat(AN_V, "L2")), q = Dict{String, Float64}(m => 1.0 for m in vcat(AN_V, "L2")))
            vals = an_vals(prof, sp.method, steps) .+ an_offset(sp)
            lam === nothing || (vals = lam(id, sp, steps, vals))
            fl = flags === nothing ? String[] : flags(id, sp)
            at = Dict{String, Any}()
            if "main" in sp.arms && sp.method in AN_V && id in (17, 18, 20)      # what the runner stores for a control of an h1_check arm
                at["control_scope"] = id == 18 ? "info" : "bound"
                at["control_deviation"] = abs(vals[end] - Float64(an_orbit(id)["exponent"]))
            end
            attrs === nothing || merge!(at, Dict{String, Any}(attrs(id, sp)))
            push!(recs, an_record(id, sp, steps; vals = vals, flags = fl, attrs = at, state = state, state3 = state3))
        end
        key = "method=$(sp.method) precision=$(sp.precision) horizon=$(sp.horizon) tau=$(sp.tau) eps=$(sp.eps) direction=$(sp.direction) mode=$(sp.mode)"
        ch = "c"^60 * lpad(string(ord), 4, '0')
        root = Dict{String, Any}("schema_version" => 2, "study" => "an-fixture", "chunk" => ord, "chunk_key" => key, "chunk_hash" => ch, "plan_sha256" => "p"^64,
            "panel_sha256" => ph, "prereg_sha256" => h, "git_sha" => "a"^40, "code_version" => "an-code-1", "backend_request" => "cpu", "backend" => "cpu",
            "device" => "cpu", "backend_fallback_reason" => "", "precision" => sp.precision, "method" => sp.method, "horizon" => sp.horizon, "tau" => sp.tau,
            "vectors" => sp.method in ("L1", "L2") ? 1 : 2, "coordinate_mode" => sp.mode, "direction" => sp.direction, "grid_sha256" => "g"^64,
            "status" => "ok", "attempt" => 1, "setup_check" => "not applicable", "first_attempt_reasons" => String[], "failure_reasons" => String[],
            "pilot" => pilot, "arms" => sp.arms, "session" => 1, "control_stop" => false, "dirty" => false, "created_utc" => "2026-10-20T00:00:00Z")
        sp.eps === nothing || (root["perturbation"] = string(sp.eps))
        AX.write_chunk(joinpath(store, "chunk-" * lpad(string(ord), 5, '0') * ".h5"), root, recs)
        push!(chunks, Dict{String, Any}("ordinal" => ord, "key" => key, "hash" => ch, "arms" => sp.arms, "status" => "ok", "attempt" => 1, "orbits" => length(sp.orbits),
            "flags" => Dict{String, Any}(f => Int[] for f in ("diagnostics exceeded", "drift", "excluded", "breakdown_tau"))))
    end
    manifest = joinpath(dir, "an-fixture-manifest.json")
    AX._write_json(manifest, Dict{String, Any}("study" => "an-fixture", "schema_version" => 2, "status" => status, "plan_sha256" => "p"^64, "chunks" => chunks,
        "prereg" => Dict("header" => h, "sha256_end" => h), "panel" => Dict("sha256_end" => ph, "build_hash" => "e"^64), "overrides" => Dict{String, Any}(),
        "reasons" => String[]))
    return (; dir, store, prereg, panel, manifest, params = params_file, p, out = joinpath(dir, "out"), specs, chunks = length(chunks), study = "an-fixture",
            prereg_hash = h, panel_hash = ph)
end
an_cached(tag; kw...) = get!(() -> an_study(tag; kw...), AN_STUDIES, tag)

# copy of a fixture study into a new directory (for malformed-input cases); returns the study with the new paths
function an_copy(st, tag::AbstractString)
    dir = an_path("copy-$tag")
    rm(dir; recursive = true, force = true)
    cp(st.dir, dir)
    rm(joinpath(dir, "out"); recursive = true, force = true)         # the copy starts without outputs
    sub(s) = replace(s, st.dir => dir)
    return merge(st, (; dir, store = sub(st.store), prereg = sub(st.prereg), panel = sub(st.panel), manifest = sub(st.manifest), params = sub(st.params),
                      out = joinpath(dir, "out")))
end

an_run(st; params = st.params, out = st.out, kw...) = AN.analyze(; prereg = st.prereg, manifest = st.manifest, panel = st.panel, params, store = st.store, out, kw...)
an_files(dir) = isdir(dir) ? sort(readdir(dir)) : String[]

# true if f() throws an ArgumentError whose message contains every word (case-insensitive): the refusal rule of E17
function an_rejects(f, words)
    try
        f()
        return false
    catch e
        e isa ArgumentError || rethrow()                    # any other error (an UndefVarError of a missing module too) is not a refusal
        msg = lowercase(sprint(showerror, e))
        return all(w -> occursin(lowercase(w), msg), words)
    end
end

# median without the Statistics library (not a test dependency)
an_median(v) = (s = sort(v); n = length(s); isodd(n) ? s[(n + 1) ÷ 2] : (s[n ÷ 2] + s[n ÷ 2 + 1]) / 2)
