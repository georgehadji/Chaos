using Random
# Helpers for the Amendment D125 suite (prefix d125_). Definitions only: nothing here calls the D125 API at load time, so a missing function
# errors the testsets that need it and nothing else. Reuses ex_helpers.jl, ex_run_helpers.jl, ex_d86_helpers.jl and ex_d94_helpers.jl unchanged.
# Spec: workspaces/module-build/stages/01-spec/output/experiments-spec.md, "Amendment D125". API assumed: the test plan, "Amendment D125 (stage 02)".
# Every reference below is written from the text of the spec (E44 to E62, R29 to R43), not from the source of the module.

# testset D125<letter> draws from StableRNG(d125_seed(letter, sub)), sub < 50; the verify-only calibration block is 65000 + ... (scratch scripts, plan)
d125_seed(letter::Char, sub::Integer = 0) = EX_SEED + 63000 + 50 * (letter - 'a' + 1) + sub

# Amendment D126 defaults, pinned from the K24 calibration (plan, "D126"): 4 x the worst realized C of the verify-only sample; not changed after the data
const D126_PAIR_CONST = 37.0   # pair_const
const D126_CTL_CONST = 56.0     # control_const

const D125_BIT = (regular = 0x01, chaotic = 0x02, weak = 0x04, strong = 0x08)

# ---------------------------------------------------------------- streams (E46 to E53, E56), hashed by the platform tool through the D94 helpers
d125_screen_text(seed, cell) = "chaoslab-screen-md-1|$(seed)|$(cell)"
d125_dev_text(seed, cell) = "chaoslab-screen-md-1|$(seed)|$(cell)|dev"
d125_panel_text(seed, cell, cls) = "chaoslab-panel-md-1|$(seed)|$(cell)|$(cls)"
d125_d_text(seed, cell, cls) = "chaoslab-d-md-1|$(seed)|$(cell)|$(cls)"
d125_pilot_text(label, cell) = "chaoslab-pilot-md-1|$(label)|$(cell)"
d125_frame_text(orbit, label) = "chaoslab-frame-1|$(orbit)|$(label)"

# E46: point j has coordinates RN(2 pi u_((j-1) d + i)), i = 1..d, RN of the 256-bit value
function d125_points(seed, cell, d, M; text = d125_screen_text(seed, cell))
    rng = d94_stream(text)
    return setprecision(BigFloat, 256) do
        tp = 2 * BigFloat(pi)
        [[Float64(tp * BigFloat(rand(rng))) for _ in 1:d] for _ in 1:M]
    end
end

# exact rational of a text such as "-0.3", "1/12", "9/10" and its BigFloat value at the current precision
function d125_rat(s::AbstractString)
    t = strip(s); neg = startswith(t, "-")
    r = d94_rat(neg ? t[2:end] : t)
    return neg ? -r : r
end
d125_big(s::AbstractString) = (r = d125_rat(s); BigFloat(numerator(r)) / BigFloat(denominator(r)))
d125_f64(s::AbstractString) = setprecision(() -> Float64(d125_big(s)), BigFloat, 256)

# E47: pairs (u1, u2) in order; q2 = q2_lo + (q2_hi - q2_lo) u1 at 256 bits, rounded once; accepted when the energy-surface start exists (status ok).
# Returns (starts, pairs consumed). p1 is evaluated in binary64 from the binary64 q2, p2 and E (reading U3 of the plan), the start is (0, q2, p1, p2).
function d125_flow_starts(seed, cell, E::AbstractString, w::NTuple{4, String}, M::Integer; maxpairs = 1000 * M, text = d125_screen_text(seed, cell))
    rng = d94_stream(text)
    Ef = d125_f64(E)
    starts = Vector{Float64}[]; pairs = 0
    while length(starts) < M && pairs < maxpairs
        u1, u2 = rand(rng), rand(rng); pairs += 1
        q2, p2 = setprecision(BigFloat, 256) do
            lo1, hi1, lo2, hi2 = (d125_big(x) for x in w)
            (Float64(lo1 + (hi1 - lo1) * BigFloat(u1)), Float64(lo2 + (hi2 - lo2) * BigFloat(u2)))
        end
        z, st = Chaoslab.energy_start(Chaoslab.HenonHeilesFlow(), Ef, q2, p2)
        st === :ok && push!(starts, Float64.(collect(z)))
    end
    return starts, pairs
end

# E49 on one record: lambda, r (nothing = no ratio), the floor phi, the rule values (a NamedTuple), Lambda_q (nothing when undefined) -> class bits
function d125_bits(lam::Float64, r, phi::Float64, rules, Lq)
    r === nothing && return 0x00
    bits = 0x00
    reg = r <= rules.ratio_regular_max && lam < phi
    cha = r >= rules.ratio_chaotic_min && lam >= phi
    reg && (bits |= D125_BIT.regular)
    if cha
        bits |= D125_BIT.chaotic
        wk = r >= rules.ratio_weak_min
        if !wk && Lq !== nothing
            thr = setprecision(() -> Float64(BigFloat(rules.weak_fraction) * BigFloat(Lq)), BigFloat, 256)
            wk = lam <= thr
        end
        bits |= wk ? D125_BIT.weak : D125_BIT.strong
    end
    return bits
end

# E49 Lambda_q: the value of rank ceil(q M_c) of the ascending list (exact integer arithmetic); nothing when M_c = 0
function d125_Lq(lams::AbstractVector, q::AbstractString)
    isempty(lams) && return nothing
    qr = d125_rat(q)
    M = length(lams)
    rank = Int(cld(numerator(qr) * M, denominator(qr)))
    return sort(Float64.(lams))[clamp(rank, 1, M)]
end

# E48 floor: max(floor_abs, floor_log ln(t) / t) at 1024 bits, rounded once
d125_phi(floor_abs::Real, floor_log::Real, t::Integer) =
    setprecision(() -> Float64(max(BigFloat(floor_abs), BigFloat(floor_log) * log(BigFloat(t)) / BigFloat(t))), BigFloat, 1024)

# E51 and E53: draw (the base partial Fisher-Yates) and subset D indices
d125_draw(seed, cell, cls, pool, size) = size <= 0 ? Int[] : ex_ref_draw(d94_stream(d125_panel_text(seed, cell, cls)), pool, size)
function d125_dindices(seed, cell, cls, n, f::AbstractString, dmin)
    n == 0 && return Int[]
    m = min(n, max(dmin, d94_dsize(d125_rat(f), n)))
    return sort(ex_ref_draw(d94_stream(d125_d_text(seed, cell, cls)), 1:n, m))
end
d125_status(size, drawn) = size == 0 ? "not drawn" : drawn == size ? "populated" : drawn == 0 ? "empty" : "short"

# E54 exact exponents at 256 bits, rounded once, descending
function d125_exact_y5(mu::AbstractVector{<:AbstractString})
    setprecision(BigFloat, 256) do
        v = BigFloat[]
        for m in mu
            x = d125_big(m)
            push!(v, log(x), -log(x))
        end
        Float64.(sort(v; rev = true))
    end
end
function d125_exact_chain_origin(N::Integer, K::AbstractString)
    setprecision(BigFloat, 256) do
        k = d125_big(K)
        x = log((2 + k + sqrt((2 + k)^2 - 4)) / 2)
        Float64.(sort(vcat(fill(x, N), fill(-x, N)); rev = true))
    end
end

# ---------------------------------------------------------------- fixtures of the panel builder (scheme "systems"); A and B differ in every value
d125_rules(; rm, cm, wm, fa, fl, wf, wq) = Dict{String, Any}("ratio_regular_max" => rm, "ratio_chaotic_min" => cm, "ratio_weak_min" => wm,
                                                              "floor_abs" => fa, "floor_log" => fl, "weak_fraction" => wf, "weak_quantile" => wq)
const D125_FA = (name = "A", seed = EX_SEED, classes = ("reg", "wk", "st"),
    rules = d125_rules(rm = 2.0, cm = 6.0, wm = 8.5, fa = 1e-4, fl = 3.0, wf = 0.25, wq = "9/10"),
    cells = Any[
        Dict{String, Any}("id" => "y3a", "kind" => "froeschle", "nu" => "0.5", "kappa" => "0.7", "mu" => "0.3", "n_start" => 48, "scr_steps" => 400,
                          "scr_tau" => 1, "scr_back" => 4, "quota" => Dict("reg" => 2, "wk" => 1, "st" => 2),
                          "d_fraction" => Dict("reg" => "1/2", "wk" => "1", "st" => "0.4"), "d_min" => 1),
        Dict{String, Any}("id" => "ch2", "kind" => "chain", "N" => 2, "K" => "1.1", "eps" => "0.2", "boundary" => "open", "n_start" => 40,
                          "scr_steps" => 400, "scr_tau" => 2, "scr_back" => 4, "quota" => Dict("reg" => 1, "wk" => 1, "st" => 2),
                          "d_fraction" => Dict("reg" => "1", "wk" => "1/2", "st" => "1/3"), "d_min" => 0),
        Dict{String, Any}("id" => "hh1", "kind" => "henon_heiles", "energy" => "1/12", "q2_lo" => "-0.3", "q2_hi" => "0.4", "p2_lo" => "-0.25",
                          "p2_hi" => "0.25", "n_start" => 24, "scr_steps" => 800, "scr_tau" => 4, "scr_back" => 4,
                          "integrator" => Dict("method" => "Leapfrog", "scheme" => "DKD", "step" => "0.0625"),
                          "quota" => Dict("reg" => 1, "wk" => 1, "st" => 2), "d_fraction" => Dict("reg" => "1", "wk" => "1", "st" => "1/2"),
                          "d_min" => 1)],
    controls = Any[Dict{String, Any}("id" => "c_y5", "kind" => "y5", "dimension" => 4, "mu" => ["3", "5/2"]),
                   Dict{String, Any}("id" => "c_ch", "kind" => "chain_origin", "N" => 2, "K" => "1.5", "boundary" => "open", "eps" => "0")])
const D125_FB = (name = "B", seed = 777, classes = ("ca", "cb", "cc"),
    rules = d125_rules(rm = 2.4, cm = 5.0, wm = 7.0, fa = 3e-4, fl = 2.0, wf = 0.5, wq = "0.8"),
    cells = Any[
        Dict{String, Any}("id" => "y3b", "kind" => "froeschle", "nu" => "0.9", "kappa" => "1.3", "mu" => "0.15", "n_start" => 36, "scr_steps" => 360,
                          "scr_tau" => 2, "scr_back" => 3, "quota" => Dict("ca" => 3, "cb" => 2, "cc" => 1),
                          "d_fraction" => Dict("ca" => "1/3", "cb" => "1", "cc" => "0.5"), "d_min" => 2),
        Dict{String, Any}("id" => "chr3", "kind" => "chain", "N" => 3, "K" => "0.8", "eps" => "0.35", "boundary" => "ring", "n_start" => 30,
                          "scr_steps" => 300, "scr_tau" => 1, "scr_back" => 5, "quota" => Dict("ca" => 2, "cb" => 2, "cc" => 3),
                          "d_fraction" => Dict("ca" => "1/2", "cb" => "1/4", "cc" => "1"), "d_min" => 0),
        Dict{String, Any}("id" => "hh2", "kind" => "henon_heiles", "energy" => "1/10", "q2_lo" => "-0.2", "q2_hi" => "0.35", "p2_lo" => "-0.3",
                          "p2_hi" => "0.3", "n_start" => 20, "scr_steps" => 640, "scr_tau" => 2, "scr_back" => 4,
                          "integrator" => Dict("method" => "KahanLi", "scheme" => "8", "step" => "0.125"),
                          "quota" => Dict("ca" => 2, "cb" => 1, "cc" => 2), "d_fraction" => Dict("ca" => "1", "cb" => "1/2", "cc" => "1/2"),
                          "d_min" => 0)],
    controls = Any[Dict{String, Any}("id" => "c_y5b", "kind" => "y5", "dimension" => 6, "mu" => ["2", "3/2", "4"])])

# the panel block (the [panel] table of the preregistration) of a fixture
function d125_block(fx)
    return Dict{String, Any}("scheme" => "systems", "seed" => fx.seed, "classes" => Dict{String, Any}("regular" => fx.classes[1], "weak" => fx.classes[2],
                             "strong" => fx.classes[3]), "rules" => fx.rules, "cell" => fx.cells, "control" => fx.controls)
end
d125_deep(x) = x isa AbstractDict ? Dict{String, Any}(k => d125_deep(v) for (k, v) in x) : x isa AbstractVector ? Any[d125_deep(v) for v in x] : x
# an edited deep copy: d125_edit(block, "rules" => "ratio_regular_max" => 7.0, ...) ; a path is a tuple of keys, an integer indexes a vector
function d125_edit(block, edits::Pair...)
    b = d125_deep(block)
    for (path, v) in edits
        node = b
        p = path isa Tuple ? path : (path,)
        for k in p[1:(end - 1)]
            node = node[k]
        end
        v === :delete ? delete!(node, p[end]) : (node[p[end]] = v)
    end
    return b
end

function d125_prereg(block::AbstractDict, tag::AbstractString)
    pre = ex_path("d125-prereg-" * tag * "-" * ex_sha256(ex_prereg_text(block))[1:12] * ".md")
    isfile(pre) || ex_write_prereg(pre, block)
    return pre
end

const D125_BUILDS = Dict{String, Any}()
function d125_build(tag::AbstractString, fx; block = d125_block(fx), kw...)
    get!(D125_BUILDS, tag) do
        pre = d125_prereg(block, tag)
        out = ex_path("d125-$tag-panel.toml")
        EX.build_panel(pre; out = out, kw...)
        (panel = TOML.parsefile(out), path = out, prereg = pre, pools = TOML.parsefile(out * ".pools.toml"), fx = fx, block = block)
    end
end

# 0 files check: a refused build writes neither the panel nor a sidecar
function d125_refused(block::AbstractDict, keys::AbstractVector; tag = "bad")
    pre = d125_prereg(block, tag)
    out = ex_path("d125-refused-$(abs(hash(block)))-panel.toml")
    ok = ex_rejects(() -> EX.build_panel(pre; out = out), keys)
    return ok && !isfile(out) && !isfile(out * ".pools.toml") && !isfile(out * ".provenance.toml")
end

# ---------------------------------------------------------------- the hand panel of the runner tests (R33): ids 1.. in (cell, class, draw index) order, controls last
# Panel orbit keys pinned by the plan (reading U1): id, cell, class, system, param, start, start_dec, control, control_kind, point_id, draw_index, in_D,
# run_seed, orbit_seed, and for controls dimension, exact_exponents.
const D125_PARAM = Dict("y3a" => "0.5;0.7;0.3", "ch2" => "2;1.1;0.2;open", "hh1" => "1/12")
const D125_KIND = Dict("y3a" => "froeschle", "ch2" => "chain", "hh1" => "henon_heiles")
const D125_CELLS = [("y3a", "reg", 2), ("y3a", "st", 1), ("ch2", "reg", 1), ("ch2", "st", 2), ("hh1", "st", 3)]
function d125_panel_orbits(cells = D125_CELLS; controls::Bool = true)
    rng = StableRNG(d125_seed('l'))
    orbs = Dict{String, Any}[]
    Ef = 1 / 12
    for (cell, cls, n) in cells, t in 1:n
        id = length(orbs) + 1
        start = if cell == "hh1"
            z = nothing
            while z === nothing
                q2 = -0.3 + 0.7 * rand(rng); p2 = -0.25 + 0.5 * rand(rng)
                z, st = Chaoslab.energy_start(Chaoslab.HenonHeilesFlow(), Ef, q2, p2)
                st === :ok || (z = nothing)
            end
            Float64.(collect(z))
        else
            [2pi * rand(rng) for _ in 1:4]
        end
        push!(orbs, Dict{String, Any}("id" => id, "cell" => cell, "class" => cls, "system" => D125_KIND[cell], "param" => D125_PARAM[cell], "start" => start,
                                      "start_dec" => string.(start), "control" => false, "control_kind" => "none", "point_id" => 10 * id, "draw_index" => t,
                                      "in_D" => t == 1, "run_seed" => string(EX_RUN_SEED), "orbit_seed" => string(EX.orbit_seed(EX_RUN_SEED, id))))
    end
    if controls
        ex = d125_exact_y5(["3", "5/2"])
        id = length(orbs) + 1
        push!(orbs, Dict{String, Any}("id" => id, "cell" => "controls", "class" => "control", "system" => "linear", "param" => "3;5/2", "start" => zeros(4),
                                      "start_dec" => fill("0", 4), "control" => true, "control_kind" => "y5", "dimension" => 4, "exact_exponents" => ex,
                                      "in_D" => false, "run_seed" => string(EX_RUN_SEED), "orbit_seed" => string(EX.orbit_seed(EX_RUN_SEED, id))))
        ex = d125_exact_chain_origin(2, "1.5")
        id += 1
        push!(orbs, Dict{String, Any}("id" => id, "cell" => "controls", "class" => "control", "system" => "chain", "param" => "2;1.5;0;open", "start" => zeros(4),
                                      "start_dec" => fill("0", 4), "control" => true, "control_kind" => "chain_origin", "dimension" => 4, "exact_exponents" => ex,
                                      "in_D" => false, "run_seed" => string(EX_RUN_SEED), "orbit_seed" => string(EX.orbit_seed(EX_RUN_SEED, id))))
    end
    return orbs
end
function d125_write_panel(path::AbstractString; orbits = d125_panel_orbits(), accepted = EX_PANEL_HASH)
    man = Dict{String, Any}("format_version" => 1, "scheme" => "systems", "controls" => [o["id"] for o in orbits if o["control"]],
                            "prereg" => Dict{String, Any}("header_sha256" => accepted, "computed_sha256" => accepted))
    open(path, "w") do io
        TOML.print(io, Dict{String, Any}("manifest" => man, "orbit" => orbits))
    end
    return path
end

# ---------------------------------------------------------------- run blocks
d125_integrators() = Any[Dict{String, Any}("id" => "lf", "method" => "Leapfrog", "scheme" => "DKD", "step" => "0.0625"),
                         Dict{String, Any}("id" => "kl", "method" => "KahanLi", "scheme" => "8", "step" => "0.125")]
d125_cqr(; kw...) = merge(Dict{String, Any}("id" => "cq", "kind" => "cqr", "budget" => d86_budget(), "order" => 2, "vectors" => 4, "form" => "auto",
                                            "penalty_factor" => "auto", "orthogonality_tol" => 1e-2), Dict{String, Any}(string(k) => v for (k, v) in kw))
d125_clv(; kw...) = merge(Dict{String, Any}("id" => "cl", "kind" => "clv", "budget" => d86_budget(), "vectors" => 4, "orthonormalizer" => "L5",
        "storage" => "store", "block" => 0, "window" => 40, "pilot_handoffs" => 60, "transient_c" => 5.0, "transient_factor" => 2, "transient_cap" => 30,
        "gap_exclude" => Int[2], "unstable_fraction" => 0.5, "store_angles" => "summary", "probabilities" => ["1/10", "1/2"], "thresholds" => [0.1, 0.01],
        "hist_lo" => -4.0, "hist_hi" => 0.0, "hist_bins" => 8), Dict{String, Any}(string(k) => v for (k, v) in kw))
d125_tangent(; kw...) = d86_spec("tangent"; kw...)

# an arm on the new kinds: methods and specs optional, the new keys through kw (cells, vectors, frames, integrators, estimators)
function d125_arm(name, orbits, methods, horizons; inds = nothing, ests = nothing, precisions = ["Float64"], kw...)
    a = Dict{String, Any}("name" => name, "orbits" => orbits, "precisions" => collect(precisions), "horizons" => collect(horizons))
    methods === nothing || (a["methods"] = collect(methods))
    inds === nothing || (a["indicators"] = collect(inds))
    ests === nothing || (a["estimators"] = collect(ests))
    for (k, v) in kw
        a[string(k)] = v
    end
    return a
end

function d125_run(arms; integrators = d125_integrators(), estimators = Any[], specs = Any[], pilot_clv = nothing)
    d = ex_rb_common()
    for (k, v) in ("reg" => 3, "wk" => 3, "st" => 3, "control" => 3)               # reading U2: setup_check_steps is keyed by the class ids of the panel
        d["health"]["setup_check_steps"][k] = v
    end
    d["health"]["control_methods"] = ["L1", "L3", "L4", "L5"]
    d["pilot"]["horizon"] = 40
    pilot_clv === nothing || (d["pilot"]["clv_window"] = pilot_clv)
    delete!(d, "h4")
    d["arm"] = arms
    isempty(integrators) || (d["integrator"] = integrators)
    isempty(estimators) || (d["estimator"] = estimators)
    isempty(specs) || (d["indicator"] = specs)
    return d
end

# the standard fixture block: spectrum methods on the map cells, on the flow cell, with controls
function d125_block_spectra()
    arms = Any[d125_arm("maps", "P", ["L1", "L3", "L4", "L4b", "L5"], [100]; cells = ["y3a", "ch2"], vectors = "d", frames = ["a"]),
               d125_arm("flow", "P", ["L1", "L3", "L5"], [160]; cells = ["hh1"], vectors = "d", frames = ["a"], integrators = ["lf", "kl"]),
               d125_arm("ctl", "controls", ["L1", "L3", "L5"], [100]; vectors = "d", frames = ["a"], h1_check = true)]
    return d125_run(arms)
end

const D125_FIXES = Dict{String, Any}()
function d125_fix(tag::AbstractString, run = d125_block_spectra(); orbits = d125_panel_orbits(), rev::Bool = false)
    return get!(D125_FIXES, tag) do
        pp = d125_write_panel(ex_path("d125-$tag-panel.toml"); orbits = orbits)
        pr = ex_write_prereg_run(ex_path("d125-$tag-prereg.md"), run; rev = rev)
        (; run, panel = pp, prereg = pr, tag, orbits)
    end
end
d125_plan(fx; kw...) = EX.plan_run(fx.prereg, fx.panel; code_version = EX_CODE, kw...)

# refusal of a run block through the reader and a plan: the error names the key path (strings in `keys`), nothing is written
function d125_block_refused(run::AbstractDict, keys::AbstractVector; tag::AbstractString)
    pr = ex_write_prereg_run(ex_path("d125-bad-$tag-prereg.md"), run)
    pp = d125_write_panel(ex_path("d125-bad-panel.toml"))
    store = ex_path("d125-bad-$tag-store")
    ok = ex_rejects(() -> EX.plan_run(pr, pp; code_version = EX_CODE), keys)
    return ok && !isdir(store)
end

# ---------------------------------------------------------------- orbits and spectra from the lower layers (references of D125c, D125k, D125n)
function d125_system(kind::AbstractString, params::AbstractVector, T::Type; wrapped::Bool = false)
    if kind == "froeschle"
        nu, ka, mu = (Chaoslab.exact(T, p) for p in params)
        return wrapped ? Chaoslab.WrappedFroeschleMap(nu, ka, mu) : Chaoslab.FroeschleMap(nu, ka, mu)
    elseif kind == "chain"
        N, K, eps, bd = params
        N isa AbstractString && (N = parse(Int, N))
        Kt, et = Chaoslab.exact(T, K), Chaoslab.exact(T, eps)
        return wrapped ? Chaoslab.WrappedStandardMapChain(N, Kt, et; boundary = Symbol(bd)) : Chaoslab.StandardMapChain(N, Kt, et; boundary = Symbol(bd))
    end
    error("d125_system: $kind")
end
d125_method(name::AbstractString, scheme::AbstractString, T::Type) =
    name == "Leapfrog" ? Chaoslab.Leapfrog(T, Symbol(scheme)) : name == "KahanLi" ? Chaoslab.KahanLi(T, parse(Int, scheme)) : error("d125_method: $name")
# the one-step map of flow cell hh1 with an integrator row of the block
d125_flow_map(integ::AbstractDict, T::Type) =
    Chaoslab.IntegratorMap(d125_method(integ["method"], integ["scheme"], T), Chaoslab.HenonHeilesFlow(), Chaoslab.exact(T, integ["step"]))
# n steps of a map or one-step map from the start (a vector of binary64 values), states rounded once to binary64 at the steps ns
function d125_states(sys, x0::AbstractVector, ns::AbstractVector{<:Integer}, T::Type)
    x = Chaoslab.SVector{length(x0), T}(T.(x0)...)
    out = Vector{Vector{Float64}}(); n = 0
    for target in ns
        while n < target
            x = Chaoslab.forward(sys, x); n += 1
        end
        push!(out, Float64.(collect(x)))
    end
    return out
end
# lambda_1..lambda_k after n steps by the verified estimator on the item's own frame (label "a": run seed and orbit id, reading U4)
function d125_lambda(sys, x0, n::Integer, method::AbstractString, k::Integer, run_seed::Integer, orbit::Integer; T::Type = Float64, tau::Integer = 1)
    d = length(x0)
    x = Chaoslab.SVector{d, T}(T.(x0)...)
    W = Chaoslab.initial_deviations(run_seed, orbit, Val(d), Val(k), T)
    res = Chaoslab.propagate(sys, [x], n; k = k, tau = tau, seed = run_seed, reducers = (EX_RED[method](),), W0 = [W])
    return Float64.(collect(res.finals[1][1].lambda))
end

# independent text of an item on a new kind (E55, reading U5): D53 R1 item 3 lines with start= in place of x0 and p0, then frame, integrator
const D125_FIELDS_OLD = (:orbit, :system, :param, :run_seed, :orbit_seed, :method, :precision, :horizon, :tau)
function d125_item_text(it; integ = nothing)
    io = IOBuffer()
    for f in D125_FIELDS_OLD
        print(io, f, "=", getproperty(it, f), "\n")
        f === :param && print(io, "start=", join(it.start, ","), "\n")
    end
    print(io, "direction=seeded\n", "frame=", it.frame, "\n")
    integ === nothing || print(io, "integrator=", integ, "\n")
    for f in (:coordinate_mode, :vectors, :grid_sha256, :code_version, :hash_version)
        print(io, f, "=", getproperty(it, f), "\n")
    end
    return String(take!(io))
end
d125_integ_text(i::AbstractDict) = join((i["method"], get(i, "scheme", ""), i["step"]), ",")

# ---------------------------------------------------------------- orbit of a plan item, from the item's own texts (references of D125c, D125k, D125m)
d125_row(it) = first(r for r in d125_integrators() if d125_integ_text(r) == it.integrator)
d125_start(it) = parse.(Float64, split(it.start, ","))
d125_item_system(it, T::Type = Float64) = it.system == "henon_heiles" ? d125_flow_map(d125_row(it), T) : d125_system(it.system, String.(split(it.param, ";")), T)
d125_item_h(it) = it.system == "henon_heiles" ? d125_f64(d125_row(it)["step"]) : 1.0
d125_each(f, run) = for file in ex_chunk_files(run.store)
    c = run.plan.chunks[ex_ordinal(file)]
    raw = d86_raw(file)
    for i in c.items
        f(run.plan.items[i], raw, c, file)
    end
end
const D125_RUNS = Dict{String, Any}()
function d125_campaign(tag::AbstractString, run = d125_block_spectra(); orbits = d125_panel_orbits(), threads = 4, kw...)
    return get!(D125_RUNS, tag) do
        fx = d125_fix(tag, run; orbits = orbits)
        store = ex_store("d125-" * tag)
        res = ex_campaign(fx; store = store, threads = threads, kw...)
        (; fx, store, res, plan = d125_plan(fx))
    end
end
d125_variant(tag, f; orbits = d125_panel_orbits()) = (r = d125_block_spectra(); f(r); d125_plan(d125_fix(tag, r; orbits = orbits)))
