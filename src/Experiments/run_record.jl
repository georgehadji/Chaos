# Recording of one work item (R2): the verified estimator reducers run segment by segment between grid steps on the shared
# per-orbit body, so the value stored at a grid step is the final value of the same estimator run directly to that step. The
# residual diagnostics are kept by a non-raising monitor (P18): running maxima, never an exception.

"Series of one orbit and its record in a chunk file: attributes, grid, point, winding, exponents, running maxima, costs."
struct OrbitRecord
    id::Int
    attrs::Dict{String, Any}
    t::Vector{Int64}
    state::Matrix{Float64}
    winding::Union{Nothing, Matrix{Int64}}
    lyap::Dict{String, Matrix{Float64}}
    diagnostics::Dict{String, Vector{Float64}}
    thresholds::Dict{String, Float64}
    first_exceedance::Dict{String, Int}
    cost::Dict{String, Vector{Float64}}
    failure::String
end

# ---------------------------------------------------------------- the non-raising monitor
struct RunMonitor{S}
    sys::S
end

struct MonState{T}
    det_max::T
    symp_max::T
    M_max::T
end

Tangent.reducer_init(::RunMonitor, ::Type{T}, ::Val, ::Val) where {T} = MonState{T}(zero(T), zero(T), zero(T))

function Tangent.reducer_step(r::RunMonitor, st::MonState{T}, n::Int, x::SVector, W::SMatrix) where {T}
    M = tangent(r.sys, x)
    return MonState{T}(max(st.det_max, eps_det(M)), max(st.symp_max, eps_symp(M)), max(st.M_max, T(maximum(abs, M)))), W
end

Tangent.reducer_update(r::RunMonitor, st::MonState, n::Int, x::SVector, W::SMatrix) = Tangent.reducer_step(r, st, n, x, W)
Tangent.reducer_final(::RunMonitor, st::MonState) = (det_max = st.det_max, symp_max = st.symp_max, M_max = st.M_max)

# ---------------------------------------------------------------- methods
_alg(m::AbstractString) = m == "L3" ? CGS() : m == "L4" ? MGS() : m == "L4b" ? RGS() : m == "L5" ? HouseholderQR() :
                          m == "L5b" ? AugmentedQR() : LapackQR()

function _estimator(method::AbstractString, orth_tol)
    method == "L1" && return Benettin()
    method == "L9" && return SVDSpectrum(; orthogonality_tol = orth_tol)
    return QRSpectrum(_alg(method); orthogonality_tol = orth_tol)
end

"Diagnostics kept for a method (R2 item 2): none for L2; det and symplecticity; orthogonality for the QR family; sum and pairing for k = 2."
function diagnostic_names(method::AbstractString)
    method == "L2" && return String[]
    method == "L1" && return ["det", "symp"]
    method == "L9" && return ["det", "symp", "sum", "pair"]
    return ["det", "symp", "orth", "sum", "pair"]
end

# ---------------------------------------------------------------- buffers
mutable struct RecBuf
    grid::Vector{Int}
    t::Vector{Int64}
    state::Matrix{Float64}
    winding::Matrix{Int64}
    lyap::Matrix{Float64}
    det::Vector{Float64}
    symp::Vector{Float64}
    orth::Vector{Float64}
    sumd::Vector{Float64}
    pair::Vector{Float64}
    gmax::Float64
    mmax::Float64
    failrow::Int
    msg::String
end

function _new_buf(grid::Vector{Int}, k::Int, wrapped::Bool)
    n = length(grid)
    return RecBuf(grid, zeros(Int64, n), zeros(n, 2), zeros(Int64, wrapped ? n : 0, 2), zeros(n, k), zeros(n), zeros(n), zeros(n),
                  zeros(n), zeros(n), 1.0, 1.0, 0, "")
end

_toint(x) = Int64(x)
_wind!(b::RecBuf, i, x::SVector{4}) = (b.winding[i, 1] = _toint(x[3]); b.winding[i, 2] = _toint(x[4]); nothing)
_wind!(b::RecBuf, i, x::SVector{2}) = nothing

# Fields that only some estimators report, selected at compile time (no dead branch reads a missing field).
_orth(fin::NamedTuple{N}) where {N} = _orth(fin, Val(:orth_max in N))
_orth(fin, ::Val{true}) = Float64(fin.orth_max)
_orth(fin, ::Val{false}) = 0.0
_gmax(fin::NamedTuple{N}) where {N} = _gmax(fin, Val(:G_max in N))
_gmax(fin, ::Val{true}) = Float64(fin.G_max)
_gmax(fin, ::Val{false}) = 1.0

function _sumpair!(b::RecBuf, i, lam::SVector{2})
    s, p = Float64(eps_sum(lam)), Float64(eps_pair(lam))
    b.sumd[i] = i > 1 ? max(b.sumd[i - 1], s) : s
    b.pair[i] = i > 1 ? max(b.pair[i - 1], p) : p
    return nothing
end
_sumpair!(b::RecBuf, i, lam::SVector{1}) = nothing

# Store the record of grid step i from the final estimator and monitor summaries.
function _store!(b::RecBuf, i::Int, x, est, mon)
    b.t[i] = b.grid[i]
    b.state[i, 1] = Float64(x[1])
    b.state[i, 2] = Float64(x[2])
    _wind!(b, i, x)
    lam = est.lambda
    for j in axes(b.lyap, 2)
        b.lyap[i, j] = Float64(lam[j])
    end
    b.det[i] = Float64(mon.det_max)
    b.symp[i] = Float64(mon.symp_max)
    b.orth[i] = _orth(est)
    _sumpair!(b, i, lam)
    b.gmax = _gmax(est)
    b.mmax = Float64(mon.M_max)
    return nothing
end

function _fail!(b::RecBuf, i, code, reducers, st, x, W, step, orbit)
    err = code == Tangent.STATUS_REDUCER ? Tangent._reducer_error(reducers, st, orbit, :cpu) :
          Tangent._health_error(code, x, W, step, orbit, :cpu)
    b.failrow = i
    b.msg = sprint(showerror, err)
    return nothing
end

"""
Run one orbit from `x0` with deviation matrix `W0` along the grid in segments of the shared per-orbit body; a health failure is
recorded in `b` (`failrow`, `msg`), never thrown. Allocates nothing on the success path.
"""
function _record_loop!(b::RecBuf, sys, x0, W0, reducers::Tuple, st0::Tuple, tau::Int, N::Int, orbit::Int)
    x, W, st = x0, W0, st0
    n = 0
    for (i, target) in enumerate(b.grid)
        x, W, st, code, step = Tangent._orbit_body(sys, x, W, st, n, target - n, N, tau, reducers)
        if code != Tangent.STATUS_OK
            return _fail!(b, i, code, reducers, st, x, W, step, orbit)
        end
        n = target
        _store!(b, i, x, reducer_final(reducers[1], st[1]), reducer_final(reducers[2], st[2]))
    end
    return nothing
end

# L2: the shadow orbit carried across the segments of the grid.
function _record_shadow!(b::RecBuf, sys, x0, u0, d0, tau::Int, N::Int, orbit::Int)
    e = Lyapunov._blank_end(x0, d0)
    n = 0
    for (i, target) in enumerate(b.grid)
        e = Lyapunov._shadow_body(sys, x0, u0, e, n, target - n, N, tau, d0)
        if e.code != 0
            b.failrow = i
            b.msg = sprint(showerror, Lyapunov._shadow_error(e, d0, orbit, :cpu))
            return nothing
        end
        n = target
        b.t[i] = target
        b.state[i, 1] = Float64(e.x[1])
        b.state[i, 2] = Float64(e.x[2])
        _wind!(b, i, e.x)
        b.lyap[i, 1] = Float64(total(e.sum) / target)
    end
    return nothing
end

# ---------------------------------------------------------------- system and start of an item
# The linear control has no torus: in wrapped mode it runs in its plane coordinates and the winding counts stay zero.
function _system_start(::Type{T}, item) where {T}
    lin = item.system == "linear"
    sys = if lin
        y5_control(parse(Float64, item.param)).sys
    elseif item.coordinate_mode == "wrapped"
        WrappedStandardMap(exact(T, item.param))
    else
        StandardMap(exact(T, item.param))
    end
    x = SVector(exact(T, item.x0), exact(T, item.p0))
    return sys, (item.coordinate_mode == "wrapped" && !lin ? torus_state(x) : x)
end

# ---------------------------------------------------------------- thresholds and the record
function _thresholds(::Type{T}, item, b::RecBuf, df) where {T}
    M = SMatrix{2, 2, T}(T(b.mmax), zero(T), zero(T), zero(T))
    out = Dict{String, Float64}("det" => Float64(det_threshold(df["det_residual_tol"], T, M)),
                                "symp" => Float64(symp_threshold(df["det_residual_tol"], T, M)))
    if item.vectors == 2
        s = Float64(sum_threshold(df["pairing_residual_tol"], T, b.gmax, item.tau, b.mmax))
        out["sum"] = s
        out["pair"] = s
        item.method == "L9" || (out["orth"] = Float64(orthogonality_threshold(_alg(item.method), df["orthogonality_tol"], T, 2, b.gmax + 1 / b.gmax)))
    end
    return out
end

_series(b::RecBuf, name) = name == "det" ? b.det : name == "symp" ? b.symp : name == "orth" ? b.orth : name == "sum" ? b.sumd : b.pair

function _first_exceedance(v::AbstractVector, thr, t)
    i = findfirst(>(thr), v)
    return i === nothing ? 0 : Int(t[i])
end

_series_bytes(item, n) = 8 * (1 + 2 + (item.coordinate_mode == "wrapped" ? 2 : 0) + item.vectors + length(diagnostic_names(item.method))) * n
_state_bytes(st) = isbitstype(typeof(st)) ? sizeof(st) : Base.summarysize(st)

function _orbit_record(::Type{T}, item, b::RecBuf, st0, df) where {T}
    n = length(b.grid)
    names = diagnostic_names(item.method)
    thr = item.method == "L2" ? Dict{String, Float64}() : _thresholds(T, item, b, df)
    diag = Dict{String, Vector{Float64}}(nm => copy(_series(b, nm)) for nm in names)
    fe = Dict{String, Int}(nm => _first_exceedance(diag[nm], thr[nm], b.t) for nm in names)
    resident = Tangent.bytes_per_orbit(Val(2), Val(item.vectors), T, item.coordinate_mode == "wrapped" ? 4 : 2) + _state_bytes(st0) + _series_bytes(item, n)
    attrs = Dict{String, Any}("system" => item.system, "params" => item.param, "ic" => item.x0 * "," * item.p0,
        "class_label" => item.design_class, "seed" => item.orbit_seed, "item_hash" => item.hash, "design_class" => item.design_class,
        "flags" => String[], "steps_recorded" => n)
    wind = item.coordinate_mode == "wrapped" ? copy(b.winding) : nothing
    return OrbitRecord(item.orbit, attrs, copy(b.t), copy(b.state), wind, Dict(item.method => copy(b.lyap)), diag, thr, fe,
                       Dict(item.method => Float64[item.horizon, item.vectors, 0.0, 0.0, resident]), b.msg)
end

# NaN from the failed grid row on, so the health check names the orbit.
function _blank_failed!(b::RecBuf)
    b.failrow == 0 && return
    for v in (b.state, b.lyap), i in b.failrow:size(v, 1), j in axes(v, 2)
        v[i, j] = NaN
    end
    for v in (b.det, b.symp, b.orth, b.sumd, b.pair)
        v[b.failrow:end] .= NaN
    end
    return
end

function _run_item(::Type{T}, item, ctx) where {T}
    sys, x0 = _system_start(T, item)
    W0 = get(ctx.frames, (item.precision, item.direction, item.vectors), nothing)
    W0 === nothing && (W0 = deviation_frame(T, item.direction, item.vectors))
    b = _new_buf(ctx.grids[(item.horizon, item.tau)], item.vectors, item.coordinate_mode == "wrapped")
    if item.method == "L2"
        _record_shadow!(b, sys, x0, SVector{2, T}(W0[1, 1], W0[2, 1]), exact(T, item.eps), item.tau, item.horizon, item.orbit)
        st0 = ()
    else
        reducers = (_estimator(item.method, ctx.df["orthogonality_tol"]), RunMonitor(sys))
        st0 = Tangent._init_states(reducers, T, Val(2), Val(item.vectors))
        _record_loop!(b, sys, x0, W0, reducers, st0, item.tau, item.horizon, item.orbit)
    end
    _blank_failed!(b)
    return _orbit_record(T, item, b, st0, ctx.df)
end

# Deviation frames of the items, computed once per (precision, direction, k) before the parallel region (read-only afterwards).
function prime_frames!(ctx, items)
    for it in items
        key = (it.precision, it.direction, it.vectors)
        haskey(ctx.frames, key) && continue
        T, bits = parse_precision(it.precision)
        ctx.frames[key] = bits === nothing ? deviation_frame(T, it.direction, it.vectors) :
                          setprecision(() -> deviation_frame(T, it.direction, it.vectors), BigFloat, bits)
    end
end

const WARM_STEPS = 1000   # steps of the warm-up run of a run type (an item shorter than this is warmed by itself)

# Compile every code path of the run types once (and the parallel path when `nt > 1`), outside the timed region of the first group
# (R5 item 5 measures a warm run).
function warm_up!(ctx, items, nt::Integer)
    for it in items
        key = (typeof(it), it.precision, it.method, it.vectors, it.coordinate_mode, it.system)
        key in ctx.warmed && continue
        push!(ctx.warmed, key)
        w = merge(it, (horizon = min(it.horizon, it.tau * cld(WARM_STEPS, it.tau)),))
        get!(ctx.grids, (w.horizon, w.tau)) do
            gr = ctx.plan.block.raw["grid"]
            record_grid(w.horizon, w.tau; points_per_decade = gr["points_per_decade"], first_steps = gr["first_steps"], extra_fractions = gr["extra_fractions"])
        end
        compute_records(ctx, NamedTuple[w for _ in 1:max(1, Int(nt))], nt)
    end
end

"Record one work item of the plan in its own element type and BigFloat precision."
function run_item(item, ctx)
    T, bits = parse_precision(item.precision)
    return bits === nothing ? _run_item(T, item, ctx) : setprecision(() -> _run_item(T, item, ctx), BigFloat, bits)
end

"""
    hot_loop_bytes(method; T = Float64) -> Int

Bytes allocated by the recording loop of one orbit after warm-up (D53v); zero when the loop is allocation free.
"""
function hot_loop_bytes(method::AbstractString; T::Type = Float64)
    sys = StandardMap(exact(T, 1.5))
    x0 = SVector(exact(T, 1.0), exact(T, 2.0))
    k = method in ("L1", "L2") ? 1 : 2
    W0 = deviation_frame(T, "pi/4", k)
    grid = record_grid(64, 1; points_per_decade = 5, first_steps = 10, extra_fractions = [0.1, 0.5])
    b = _new_buf(grid, k, false)
    if method == "L2"
        u0 = SVector{2, T}(W0[1, 1], W0[2, 1])
        d0 = exact(T, 1e-8)
        _record_shadow!(b, sys, x0, u0, d0, 1, 64, 1)
        return @allocated _record_shadow!(b, sys, x0, u0, d0, 1, 64, 1)
    end
    reducers = (_estimator(method, 8.0), RunMonitor(sys))
    st0 = Tangent._init_states(reducers, T, Val(2), Val(k))
    _record_loop!(b, sys, x0, W0, reducers, st0, 1, 64, 1)
    return @allocated _record_loop!(b, sys, x0, W0, reducers, st0, 1, 64, 1)
end
