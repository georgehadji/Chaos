# Reducer chains, observables, table layouts and memory of the indicator specs (Amendment D86: E23 to E26, E28). The block names kinds,
# members and observables; this file builds the verified reducers from them, so the recorder and the device path share one construction.

const IND_FREQ_COLUMNS = ["nu1", "nu2", "amp1", "amp2", "mirror1", "mirror2", "second1", "second2", "dist", "delta", "D"]
const IND_FREQ_BITS = ["below_floor", "degenerate", "failed", "tie1", "tie2", "near_zero1", "near_zero2", "mirror1", "mirror2", "refined1", "refined2"]
const IND_PA_COLUMNS = ["nu1", "nu2", "delta", "D"]
const IND_PA_BITS = ["below_floor", "undefined window"]
const IND_PA_SCALARS = ["handoffs", "near_pi", "complete"]
const IND_BK_COLUMNS = ["avg1", "avg2", "S", "B", "delta", "absdig", "reldig", "dig"]
const IND_BK_BITS = ["below_floor", "zero", "reldig_undefined"]
const IND_BK_SCALARS = ["handoffs", "complete"]
const IND_ST_SCALARS = ["handoffs", "sum_1", "sum_2", "comp_1", "comp_2", "sumsq_1", "sumsq_2", "compsq_1", "compsq_2", "mean_1", "mean_2", "std_1", "std_2",
                        "below_1", "below_2", "above_1", "above_2", "first_full_step"]
const IND_DIRECTION_RULE = "u = xi_0"

const IND_MEMBER_COLUMNS = Dict{String, Vector{String}}("SALI" => ["SALI.sali", "SALI.L"], "FLI" => ["FLI.fli", "FLI.fli1"],
    "MEGNO" => ["MEGNO.Y1", "MEGNO.Ybar1", "MEGNO.Ybar20_1", "MEGNO.Y2", "MEGNO.Ybar2", "MEGNO.Ybar20_2"], "RLI" => ["RLI.DL", "RLI.RLI"])

"A coordinate of the orbit point as an isbits observable (the lifted coordinate of a lifted state or of the torus-free control)."
struct IndCoordinate
    i::Int
end
(o::IndCoordinate)(x) = x[o.i]

# ---------------------------------------------------------------- chains (E23)
_ind_member(name, replace::Bool) = name == "SALI" ? SALI() : name == "FLI" ? FLI(; replace = replace) : name == "MEGNO" ? MEGNO(; replace = replace) : RLI(; replace = replace)

"Members of a tangent or pair chain in order: the module sets the forms, every member but the last non-replacing, the last replacing (E23)."
chain_members(members::AbstractVector{<:AbstractString}) = Tuple(_ind_member(m, i == length(members)) for (i, m) in enumerate(members))

"Record of a chain, member names with their forms, for the chunk file and the manifest."
function chain_record(spec::AbstractDict)
    spec["kind"] in ("tangent", "pair") && return ["$m(replace=$(i == length(spec["members"])))" for (i, m) in enumerate(spec["members"])]
    spec["kind"] == "stretching" && return ["stretching(replace=true)"]
    return [spec["kind"]]
end

"Member names of the chain record."
ind_members(spec::AbstractDict) = spec["kind"] in ("tangent", "pair") ? String.(spec["members"]) : [spec["kind"] == "stretching" ? "stretching" : spec["kind"]]

"Columns of the series of a spec (E24); empty for the window kinds."
function ind_series_columns(spec::AbstractDict)
    spec["kind"] in ("tangent", "pair") && return reduce(vcat, [IND_MEMBER_COLUMNS[m] for m in spec["members"]])
    return spec["kind"] == "stretching" ? ["stretching.D"] : String[]
end

# ---------------------------------------------------------------- applicability and observables (E26)
"Whether a spec is defined for an orbit of `system` (the torus-free linear control has no exp_i_q, cos or sin)."
function ind_applicable(spec::AbstractDict, system::AbstractString)
    system == "linear" || return true
    spec["kind"] == "frequency" && return spec["observable"] != "exp_i_q"
    spec["kind"] == "birkhoff" && return !any(o -> o["kind"] in ("cos", "sin"), spec["observables"])
    return true
end

ind_winding(item) = item.coordinate_mode == "wrapped" && item.system != "linear"

function ind_pa_observable(spec, item)
    spec["observable"] == "planar" && return PlanarObservable(1, 2)
    return ind_winding(item) ? LiftedAngleObservable(1, 3) : LiftedAngleObservable(1)
end

function _ind_bk_observable(o, wound::Bool)
    c = o["coordinate"]
    o["kind"] == "cos" && return CosObservable(c, o["frequency"])
    o["kind"] == "sin" && return SinObservable(c, o["frequency"])
    return wound ? LiftedObservable(c, c + 2) : IndCoordinate(c)
end

"The observables of a Birkhoff spec for an item, a tuple of isbits callables."
ind_bk_observables(spec, item) = Tuple(_ind_bk_observable(o, ind_winding(item)) for o in spec["observables"])

"The names of the observables of a Birkhoff spec for an item (the module's names)."
function ind_observable_names(spec, item)
    return String[o isa IndCoordinate ? "lifted(x[$(o.i)])" : observable_name(o) for o in ind_bk_observables(spec, item)]
end

@inline ind_freq_sample(::Val{:exp_i_q}, x, ::Type{T}) where {T} = Complex{T}(cis(x[1]))
@inline ind_freq_sample(::Val{:q_plus_i_p}, x, ::Type{T}) where {T} = Complex{T}(x[1], x[2])

# ---------------------------------------------------------------- reducers of an item
function _ind_stretching(::Type{T}, spec, sys) where {T}
    rng = spec["range"] == "derived" ? stretching_range(sys) : Tuple(spec["range"])
    return StretchingNumbers(; range = rng, bins = spec["bins"], interval = spec["interval"], min_per_cell = spec["min_per_cell"])
end

"Readers of an item without the monitor: the chain members, or the window reader (none for `frequency`)."
function ind_readers(::Type{T}, item, spec, sys) where {T}
    kind = spec["kind"]
    kind in ("tangent", "pair") && return chain_members(spec["members"])
    kind == "stretching" && return (_ind_stretching(T, spec, sys),)
    kind == "frequency" && return ()
    kind == "phase_advance" && return (WeightedPhaseAdvance(; window = Tuple(spec["windows"]), observable = ind_pa_observable(spec, item), birkhoff_order = spec["birkhoff_order"]),)
    return (WeightedBirkhoff(; window = Tuple(spec["windows"]), observables = ind_bk_observables(spec, item), birkhoff_order = spec["birkhoff_order"]),)
end

"The monitor of an item: the non-raising residual monitor of the orbit (of orbit A for a pair, on the pair map)."
ind_monitor(sys, kind) = RunMonitor(kind == "pair" ? PairMap(sys) : sys)

"The full reducer tuple of an item: members then the monitor (series kinds), the monitor then the reader (window kinds)."
function ind_chain(::Type{T}, item, spec, sys) where {T}
    rd = ind_readers(T, item, spec, sys)
    mon = ind_monitor(sys, spec["kind"])
    return ind_window_kind(spec["kind"]) ? (mon, rd...) : (rd..., mon)
end

# ---------------------------------------------------------------- memory per orbit (E28)
_ind_prec_type(label) = _element_type(_label(label))

"""
    orbit_bytes(plan, item) -> (M, R, S, B, total, ng, c)

Memory per orbit of an indicator item (E28): the stream `M`, the reducers and the workspace `R` (read from the reducers), the series `S = 8 c ng` and the
frequency buffer `B`. Items of the exponent methods are not covered.
"""
function orbit_bytes(plan::RunPlan, it)
    it.kind === nothing && throw(ArgumentError("orbit_bytes: item $(it.id) is not an indicator item"))
    spec = block_spec(plan.block, it.method)
    T = _ind_prec_type(it.precision)
    kind = spec["kind"]
    d, pair = kind == "pair" ? 4 : 2, kind == "pair"
    s = ind_winding(it) ? 4 : 2
    M = Tangent.bytes_per_orbit(Val(d), Val(it.vectors), T, pair ? 2s : s)
    R = sum(reducer_bytes(r, T, Val(d), Val(it.vectors)) for r in ind_chain(T, it, spec, StandardMap(one(T))))
    gr = plan.block.raw["grid"]
    ng = length(record_grid(it.horizon, it.tau; points_per_decade = gr["points_per_decade"], first_steps = gr["first_steps"], extra_fractions = gr["extra_fractions"]))
    c = 1 + 2 + (it.coordinate_mode == "wrapped" ? 2 : 0) + 2 + length(ind_series_columns(spec))
    B = kind == "frequency" ? 2 * maximum(spec["windows"]) * sizeof(Complex{T}) : 0
    return (M = M, R = R, S = 8 * c * ng, B = B, total = M + R + 8 * c * ng + B, ng = ng, c = c)
end

"""
    group_capacity(plan, item; free_bytes, f) -> Int

Orbits per launch group for an indicator item: `floor(f free_bytes / bytes(orbit))` with exact arithmetic (E28, R17).
"""
function group_capacity(plan::RunPlan, it; free_bytes::Integer, f::Real = Tangent.CHUNK_HEADROOM)
    b = orbit_bytes(plan, it)
    r = Rational{BigInt}(f)
    return Int(fld(numerator(r) * big(free_bytes), denominator(r) * b.total))
end

# ---------------------------------------------------------------- checks at plan time (R13 item 1)
# `b` is the orbit base of the panel; the run check of the reducers of the spec with the item's tau and horizon, the existence of the derived range, the
# observable at the start. A failure names the spec, the orbit and the parameter; nothing is run.
function ind_check_plan(spec::AbstractDict, b, tau::Integer, N::Integer, mode::AbstractString)
    kind, id = spec["kind"], spec["id"]
    fail(what) = throw(ArgumentError("indicator spec `$id`, orbit $(b.id): $what"))
    if kind == "stretching"
        spec["range"] == "derived" && !(b.system in ("standard", "linear")) &&
            fail("parameter `range` = \"derived\" is defined for the standard maps and 2 x 2 linear maps only, the system is $(b.system)")
        _ind_try(fail, () -> check_run(StretchingNumbers(; range = (-1.0, 1.0), bins = spec["bins"], interval = spec["interval"], min_per_cell = spec["min_per_cell"]), tau, N))
    end
    ind_window_kind(kind) || return nothing
    item = (id = b.id, system = b.system, param = b.param, x0 = b.x0, p0 = b.p0, coordinate_mode = mode)
    x = try
        last(_system_start(Float64, item))
    catch e
        fail("the start ($(b.x0), $(b.p0)) cannot be formed: " * first(sprint(showerror, e), 120))
    end
    all(isfinite, x) || fail("parameter `observable`: the start ($(b.x0), $(b.p0)) is not finite")
    kind == "frequency" && (isfinite(ind_freq_sample(Val(Symbol(spec["observable"])), x, Float64)) || fail("parameter `observable` is not finite at the start"))
    kind == "phase_advance" && _ind_try(fail, () -> check_run(first(ind_readers(Float64, item, spec, nothing)), tau, N; d = 2))
    kind == "birkhoff" && _ind_try(fail, () -> check_run(first(ind_readers(Float64, item, spec, nothing)), tau, N; x0 = x))
    return nothing
end

_ind_try(fail, f) = try
    f()
catch e
    e isa ArgumentError ? fail(e.msg) : rethrow()
end
