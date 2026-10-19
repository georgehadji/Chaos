# Health checks (R5 items 1 and 2) of one chunk attempt. Thresholds come from the block and the defaults; nothing reads a
# hypothesis quantity. Checks flag orbits in place (`diagnostics exceeded`, `drift`) and return the failure reasons.

"Whether a running maximum of the orbit is above the threshold in force (strict)."
exceeded(o::OrbitRecord) = any(nm -> o.diagnostics[nm][end] > o.thresholds[nm], keys(o.diagnostics))

# Lifted |p| above the bound at any recorded step (wrapped mode adds the winding count).
drifted(o::OrbitRecord, bound) = _drifted(o.state, o.winding, bound)
_drifted(state, ::Nothing, bound) = any(i -> abs(state[i, 2]) > bound, axes(state, 1))
_drifted(state, winding::Matrix, bound) = any(i -> abs(state[i, 2] + 2pi * winding[i, 2]) > bound, axes(state, 1))

function _first_nonfinite(o::OrbitRecord)
    all(isfinite, o.state) || return "state"
    for (m, v) in o.lyap
        all(isfinite, v) || return "lyap/" * m
    end
    for (m, v) in o.diagnostics
        all(isfinite, v) || return "diagnostics/" * m
    end
    return ""
end

# D56: an orthonormalisation breakdown at a large reorthonormalisation interval is a limit of the method, a result of the tau sweep,
# not a defect: (kind, step) of such a failure text, else nothing.
function breakdown_of(failure::AbstractString)
    m = match(r"NumericalHealthError: (orthonormalisation breakdown) = .* at step (\d+),", failure)
    return m === nothing ? nothing : (kind = String(m.captures[1]), step = parse(Int, m.captures[2]))
end

"A chunk whose arms are all marked `stop_exempt` in the block (D56)."
chunk_exempt(plan::RunPlan, chunk) = !isempty(chunk.arms) && all(a -> any(x -> x["name"] == a && get(x, "stop_exempt", false), plan.block.raw["arm"]), chunk.arms)

# Flag the breakdown orbits of an exempt chunk (`breakdown_tau`, with step and kind) and return the other orbits.
function _flag_breakdowns!(recs)
    for o in recs
        b = breakdown_of(o.failure)
        b === nothing && continue
        push!(o.attrs["flags"], "breakdown_tau")
        o.attrs["breakdown_step"], o.attrs["breakdown_kind"] = b.step, b.kind
    end
    return OrbitRecord[o for o in recs if !("breakdown_tau" in o.attrs["flags"])]
end

function nonfinite_reasons(recs)
    out = String[]
    for o in recs
        what = _first_nonfinite(o)
        isempty(o.failure) || push!(out, "numerical health error: " * o.failure)
        isempty(what) || !isempty(o.failure) || push!(out, "non-finite value in $what of orbit $(o.id)")
    end
    return out
end

# Exact-exponent controls in an h1_check arm for a method of control_methods; exponent zero is not checked (R5 item 1c). With
# `control_scope = "exact_start"` (D56) only controls whose binary64 start is an exact periodic point are under the bound; the
# others keep their measured deviation as `control_scope = "info"` on the orbit (reported in the manifest, never a stop).
function control_reasons(plan::RunPlan, chunk, recs)
    h = plan.block.raw["health"]
    h1 = Set(a["name"] for a in plan.block.raw["arm"] if get(a, "h1_check", false))
    out = String[]
    for (i, o) in zip(chunk.items, recs)
        it = plan.items[i]
        (it.exponent === nothing || it.exponent == 0 || !(it.method in h["control_methods"]) || !any(in(h1), it.arms)) && continue
        err = abs(o.lyap[it.method][end, 1] - it.exponent)
        if get(h, "control_scope", "all") == "exact_start" && isfinite(err)
            o.attrs["control_scope"], o.attrs["control_deviation"] = it.control_exact ? "bound" : "info", err
            it.control_exact || continue
        end
        err > h["control_bound"] && push!(out, "control bound exceeded: orbit $(o.id) method $(it.method), |lambda_1 - exact| = $err above $(h["control_bound"])")
    end
    return out
end

# Float64 states after the steps `ns` (sorted) by direct iteration of the map, for the setup check.
function _f64_states(item, ns)
    sys, x = _system_start(Float64, item)
    out, n = Vector{Vector{Float64}}(), 0
    for target in ns
        while n < target
            x = forward(sys, x)
            n += 1
        end
        push!(out, Float64.(collect(x)))
    end
    return out
end

_coord_diff(wrapped, a, b) = wrapped ? max(torus_delta(a[1] - b[1]), torus_delta(a[2] - b[2])) : max(abs(a[1] - b[1]), abs(a[2] - b[2]))

# Setup check of one orbit: "ok", "failed" (with the largest difference) or "no reference".
function _setup_one(ctx, item, o::OrbitRecord)
    h = ctx.plan.block.raw["health"]
    span = h["setup_check_steps"][item.setup_class]
    ref = get(ctx.ref, (item.orbit, item.method, item.tau, item.direction, item.coordinate_mode), nothing)
    rows = [i for (i, n) in enumerate(o.t) if n <= span]
    ref === nothing && return "no reference", 0.0
    refgrid = ctx.grids[(ref.horizon, ref.tau)]
    all(i -> o.t[i] in refgrid, rows) || return "no reference", 0.0
    states = _f64_states(item, o.t[rows])
    d = maximum(_coord_diff(item.coordinate_mode == "wrapped", o.state[i, 1:2], s) for (i, s) in zip(rows, states); init = 0.0)
    return (d > h["setup_check_tol"] ? "failed" : "ok"), d
end

"Setup check of a chunk of an extended precision (R5 item 1d): `(status, reasons)`; Float64 chunks are \"not applicable\"."
function setup_check(ctx, chunk, recs)
    _label(chunk.run_type.precision) == "Float64" && return "not applicable", String[]
    res = [(_setup_one(ctx, ctx.plan.items[i], o), o.id) for (i, o) in zip(chunk.items, recs)]
    tol = ctx.plan.block.raw["health"]["setup_check_tol"]
    bad = ["setup check failed: orbit $id differs from the Float64 run by $d (tolerance $tol)" for (((st, d), id)) in res if st == "failed"]
    isempty(bad) || return "failed", bad
    return (all(r -> r[1][1] == "no reference", res) ? "no reference" : "ok"), String[]
end

"""
    check_chunk!(ctx, chunk, recs) -> (reasons, fatal, setup)

Run the checks of R5 on the records of one attempt: flag `diagnostics exceeded` and `drift` in place; `reasons` are the failure
reasons (non-finite values or a health error, the exceedance fraction, the setup check), `fatal` the control-bound reasons (a code
defect that stops the campaign), `setup` the setup-check status.
"""
function check_chunk!(ctx, chunk, recs)
    h = ctx.plan.block.raw["health"]
    reasons = nonfinite_reasons(chunk_exempt(ctx.plan, chunk) ? _flag_breakdowns!(recs) : recs)
    nex = 0
    for o in recs
        isempty(_first_nonfinite(o)) || continue
        exceeded(o) && (push!(o.attrs["flags"], "diagnostics exceeded"); nex += 1)
        drifted(o, h["drift_bound"]) && push!(o.attrs["flags"], "drift")
    end
    n = length(recs)
    nex > h["fail_fraction"] * n * (1 + STRICT_SLACK) &&
        push!(reasons, "diagnostics exceeded in $nex of $n orbits (more than fail_fraction $(h["fail_fraction"]))")
    status, sreasons = setup_check(ctx, chunk, recs)
    append!(reasons, sreasons)
    return (reasons = reasons, fatal = control_reasons(ctx.plan, chunk, recs), setup = status)
end

"Mark every orbit of a chunk that failed twice as excluded from the analysis (R5 item 3)."
function mark_excluded!(recs)
    for o in recs
        "excluded" in o.attrs["flags"] || push!(o.attrs["flags"], "excluded")
    end
    return recs
end

flagged_ids(recs, flag) = [o.id for o in recs if flag in o.attrs["flags"]]

"Measured deviations of the controls outside the control bound (D56): `(orbit, method, deviation)` per orbit marked `info`."
control_info(recs) = [Dict{String, Any}("orbit" => o.id, "method" => first(keys(o.lyap)), "deviation" => o.attrs["control_deviation"])
                      for o in recs if get(o.attrs, "control_scope", "") == "info"]
