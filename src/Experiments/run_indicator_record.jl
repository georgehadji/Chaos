# Recording of the indicator items on the CPU (Amendment D86, R14): the orbit is the identical orbit of every other method, advanced by the shared
# per-orbit body of the Tangent layer; the verified reducers are called, never re-implemented (O16). The tangent, pair and stretching kinds are read at
# the recording grid and stored as series; the window kinds are stepped one handoff at a time, the deviation column reset at each handoff (it is never
# read, and an unnormalised column would overflow in a long run), and read once at the run end.

# The non-raising monitor of a pair stream reads orbit A (the first half of the pair state) on the base map.
function Tangent.reducer_step(r::RunMonitor{<:PairMap}, st::MonState{T}, n::Int, x::SVector, W::SMatrix) where {T}
    M = tangent(r.sys.sys, Systems._halves(x)[1])
    return MonState{T}(max(st.det_max, eps_det(M)), max(st.symp_max, eps_symp(M)), max(st.M_max, T(maximum(abs, M)))), W
end
Tangent.reducer_update(r::RunMonitor{<:PairMap}, st::MonState, n::Int, x::SVector, W::SMatrix) = Tangent.reducer_step(r, st, n, x, W)

# ---------------------------------------------------------------- buffers and rows
mutable struct IndBuf
    grid::Vector{Int}
    t::Vector{Int64}
    state::Matrix{Float64}
    winding::Matrix{Int64}
    det::Vector{Float64}
    symp::Vector{Float64}
    series::Matrix{Float64}
    mmax::Float64
    failrow::Int
    msg::String
end

function _new_ind_buf(grid::Vector{Int}, wrapped::Bool, ncols::Int)
    n = length(grid)
    return IndBuf(grid, zeros(Int64, n), zeros(n, 2), zeros(Int64, wrapped ? n : 0, 2), zeros(n), zeros(n), zeros(n, ncols), 1.0, 0, "")
end

struct HostRow{M}
    S::M
    i::Int
end
@inline _set!(h::HostRow, j, v) = (h.S[h.i, j] = v; nothing)

# The columns of E24 of one member at a grid step, written through `_set!` (the host series and the device record array share them).
@inline _put!(row, off, r, st, ws) = _put!(row, off, r, st)
@inline function _put!(row, off, r::SALI, st)
    f = reducer_final(r, st)
    _set!(row, off + 1, Float64(f.sali)); _set!(row, off + 2, Float64(f.L))
    return off + 2
end
@inline function _put!(row, off, r::FLI, st)
    f = reducer_final(r, st)
    _set!(row, off + 1, Float64(f.fli)); _set!(row, off + 2, Float64(f.fli1))
    return off + 2
end
@inline function _put!(row, off, r::MEGNO, st)
    f = reducer_final(r, st)
    for c in 1:2
        _set!(row, off + 3c - 2, Float64(f.Y[c])); _set!(row, off + 3c - 1, Float64(f.Ybar[c])); _set!(row, off + 3c, Float64(f.Ybar20[c]))
    end
    return off + 6
end
@inline function _put!(row, off, r::RLI, st)
    f = reducer_final(r, st)
    _set!(row, off + 1, Float64(f.DL)); _set!(row, off + 2, Float64(f.RLI))
    return off + 2
end
# D from the two count columns of the workspace, without the count matrix that `reducer_final(r, st, ws)` copies (same function, same operands, same bits).
@inline function _stretch_D(r::StretchingNumbers, st, ws)
    B, T = r.bins, typeof(st.value)
    return spectral_distance(view(ws, 1:(B + 2)), view(ws, (B + 3):(2B + 4)), (T(r.hi) - T(r.lo)) / B)
end
@inline function _put!(row, off, r::StretchingNumbers, st, ws)
    _set!(row, off + 1, Float64(_stretch_D(r, st, ws)))
    return off + 1
end
@inline _put_members!(row, off, ::Tuple{}, ::Tuple{}, ws) = off
@inline function _put_members!(row, off, rs::Tuple, ss::Tuple, ws)
    off = _put!(row, off, first(rs), first(ss), ws)
    return _put_members!(row, off, Base.tail(rs), Base.tail(ss), ws)
end

@inline function _ind_point!(b::IndBuf, i, x, wound::Bool)
    b.t[i] = b.grid[i]
    b.state[i, 1] = Float64(x[1])
    b.state[i, 2] = Float64(x[2])
    wound && (b.winding[i, 1] = _toint(x[3]); b.winding[i, 2] = _toint(x[4]))
    return nothing
end

@inline function _ind_monitor!(b::IndBuf, i, mon)
    b.det[i] = Float64(mon.det_max)
    b.symp[i] = Float64(mon.symp_max)
    b.mmax = Float64(mon.M_max)
    return nothing
end

function _ind_fail!(b::IndBuf, i, code, red, st, x, W, step, orbit)
    err = code == Tangent.STATUS_REDUCER ? Tangent._reducer_error(red, st, orbit, :cpu) : Tangent._health_error(code, x, W, step, orbit, :cpu)
    b.failrow = i
    b.msg = sprint(showerror, err)
    return nothing
end

# ---------------------------------------------------------------- the loops
"Tangent, pair and stretching kinds: advance along the grid in segments of the shared body, store the columns of E24 at each grid step."
function _ind_series_loop!(b::IndBuf, sys, x0, W0, red::Tuple, st0::Tuple, ws, tau::Int, N::Int, orbit::Int, wound::Bool)
    x, W, st = x0, W0, st0
    n = 0
    for (i, target) in enumerate(b.grid)
        x, W, st, code, step = Tangent._orbit_body(sys, x, W, st, n, target - n, N, tau, red, ws)
        code == Tangent.STATUS_OK || return _ind_fail!(b, i, code, red, st, x, W, step, orbit)
        n = target
        _ind_point!(b, i, x, wound)
        _put_members!(HostRow(b.series, i), 0, Base.front(red), Base.front(st), ws)
        _ind_monitor!(b, i, reducer_final(last(red), last(st)))
    end
    return st
end

"Window kinds: one handoff at a time with the constant deviation column, `feed(j, x)` called at handoff `j`, the orbit point and the monitor stored at the grid steps."
function _ind_window_loop!(b::IndBuf, sys, x0, W1, red::Tuple, st0::Tuple, tau::Int, N::Int, orbit::Int, wound::Bool, feed)
    x, st = x0, st0
    n, i, ng = 0, 1, length(b.grid)
    for j in 1:(N ÷ tau)
        x, W, st, code, step = Tangent._orbit_body(sys, x, W1, st, n, tau, N, tau, red)
        code == Tangent.STATUS_OK || return _ind_fail!(b, i, code, red, st, x, W, step, orbit)
        n += tau
        feed(j, x)
        if i <= ng && n == b.grid[i]
            _ind_point!(b, i, x, wound)
            _ind_monitor!(b, i, reducer_final(first(red), first(st)))
            i += 1
        end
    end
    return st
end

# NaN from the failed grid row on, so the health check names the orbit.
function _blank_ind!(b::IndBuf)
    b.failrow == 0 && return
    for v in (b.state, b.series), i in b.failrow:size(v, 1), j in axes(v, 2)
        v[i, j] = NaN
    end
    for v in (b.det, b.symp)
        v[b.failrow:end] .= NaN
    end
    return
end

# ---------------------------------------------------------------- the result tables of a kind
_ind_nt(; kw...) = (; series = nothing, columns = String[], final = nothing, rows = String[], fcolumns = String[], flags = nothing, bits = String[], scalars = nothing,
                   scolumns = String[], range = nothing, counts = nothing, pool = nothing, kw...)

function _stretch_tables(r::StretchingNumbers, fin, grid, tau, b::IndBuf, store_counts::Bool)
    spec_first = findfirst(g -> g ÷ tau >= r.min_per_cell * r.bins, grid)
    f(v) = Float64.(collect(v))
    sc = Float64[fin.handoffs; f(fin.sum); f(fin.comp); f(fin.sum2); f(fin.comp2); f(fin.mean); f(fin.std); f(fin.below); f(fin.above);
                 spec_first === nothing ? 0.0 : Float64(grid[spec_first])]
    return _ind_nt(series = b.series, columns = ["stretching.D"], scalars = sc, scolumns = IND_ST_SCALARS, range = Float64[r.lo, r.hi],
                   counts = store_counts ? fin.counts : nothing, pool = fin.counts)
end

_flag_word(bits::Vector{Bool}) = sum((Int64(1) << (i - 1) for i in eachindex(bits) if bits[i]); init = Int64(0))

function _pa_tables(r::WeightedPhaseAdvance, fin, windows)
    H = length(windows)
    final = Matrix{Float64}(undef, H, 4)
    flags = Vector{Int64}(undef, H)
    for h in 1:H
        final[h, :] .= (Float64(fin.nu1[h]), Float64(fin.nu2[h]), Float64(fin.delta[h]), Float64(fin.D[h]))
        flags[h] = _flag_word([fin.below_floor[h], isnan(fin.delta[h])])
    end
    return _ind_nt(final = final, rows = ["W=$w" for w in windows], fcolumns = IND_PA_COLUMNS, flags = flags, bits = IND_PA_BITS,
                   scalars = Float64[fin.handoffs, fin.near_pi, fin.complete ? 1.0 : 0.0], scolumns = IND_PA_SCALARS)
end

function _bk_tables(r::WeightedBirkhoff, fin, windows, names)
    H, K = length(windows), length(names)
    final = Matrix{Float64}(undef, H * K, 8)
    flags = Vector{Int64}(undef, H * K)
    rows = String[]
    for h in 1:H, k in 1:K
        i = (h - 1) * K + k
        final[i, :] .= (Float64(fin.avg1[k, h]), Float64(fin.avg2[k, h]), Float64(fin.S[k, h]), Float64(fin.B[h]), Float64(fin.delta[k, h]),
                        Float64(fin.absdig[k, h]), Float64(fin.reldig[k, h]), Float64(fin.dig[k, h]))
        flags[i] = _flag_word([fin.below_floor[k, h], fin.zero[k, h], isnan(fin.reldig[k, h]) && !fin.zero[k, h]])
        push!(rows, "W=$(windows[h]) $(names[k])")
    end
    return _ind_nt(final = final, rows = rows, fcolumns = IND_BK_COLUMNS, flags = flags, bits = IND_BK_BITS,
                   scalars = Float64[fin.handoffs, fin.complete ? 1.0 : 0.0], scolumns = IND_BK_SCALARS)
end

# E1: the pure function of naff-fma on the first 2 W samples of the series buffer, one row per window; the failure of a non-finite sample is returned for the host.
function _freq_tables(spec, buf::AbstractVector{Complex{T}}, orbit::Int) where {T}
    Ws = spec["windows"]
    final = Matrix{Float64}(undef, length(Ws), 11)
    flags = Vector{Int64}(undef, length(Ws))
    failure = ""
    for (r, W) in enumerate(Ws)
        m = frequency_map(view(buf, 1:(2W)), W; order = spec["order"], positions = Symbol(spec["positions"]), pad = spec["pad"], mean_removal = spec["mean_removal"],
                          mirror_ratio = spec["mirror_ratio"], floor_multiplier = spec["floor_multiplier"])
        a, c = m.w1[1], m.w2[1]
        final[r, :] .= Float64.((a.nu, c.nu, a.amplitude, c.amplitude, a.mirror_amplitude, c.mirror_amplitude, a.second_ratio, c.second_ratio, m.dist[1], m.delta, m.D))
        flags[r] = _flag_word([m.below_floor, m.status[1] === :degenerate, m.status[1] === :failed, a.tie, c.tie, a.near_zero, c.near_zero, a.mirror, c.mirror, false, false])
        if m.failure !== nothing && isempty(failure)
            failure = "NumericalHealthError: series = NaN, non-finite sample $(m.failure.sample) of window $(m.failure.window) (plane $(m.failure.plane)) of spec " *
                      "`$(spec["id"])` at orbit $orbit, backend cpu"
        end
    end
    return _ind_nt(final = final, rows = ["W=$w" for w in Ws], fcolumns = IND_FREQ_COLUMNS, flags = flags, bits = IND_FREQ_BITS), failure
end

# ---------------------------------------------------------------- pair
function _ind_xi0(::Type{T}, item, spec, ctx) where {T}
    spec["frame_source"] == "direction" && return SVector{2, T}(deviation_frame(T, item.direction, 2)[1, 1], deviation_frame(T, item.direction, 2)[2, 1])
    W = Tangent.initial_deviations(parse(Int, item.run_seed), item.orbit, Val(2), Val(1), T)
    return SVector{2, T}(W[1, 1], W[2, 1])
end

function _pair_failed(b, msg, spec, item, xi0)
    b.failrow, b.msg = 1, msg
    return nothing, Dict{String, Any}("Dx" => spec["separation"], "Dx_r" => NaN, "mode" => item.coordinate_mode, "xi0" => Float64.(collect(xi0)))
end

# Pair start and block frame (E30); a refusal is a health error of the orbit, recorded with the pair id (R13 item 2, R18).
function _pair_setup(::Type{T}, item, spec, ctx, b, sys) where {T}
    xi0 = _ind_xi0(T, item, spec, ctx)
    x = SVector(exact(T, item.x0), exact(T, item.p0))
    ps = try
        pair_start(sys, x, xi0, spec["separation"]; orbit = item.orbit)
    catch e
        e isa NumericalHealthError || rethrow()
        return _pair_failed(b, sprint(showerror, e), spec, item, xi0)..., xi0
    end
    z = zero(T)
    W = SMatrix{4, 2, T}(xi0[1], xi0[2], z, z, z, z, xi0[1], xi0[2])
    ctx.hooks.frame === nothing || (W = ctx.hooks.frame(:pair_frame, (precision = item.precision, orbit = item.orbit), W))
    try
        check_frame(RLI(), W; orbit = item.orbit)
    catch e
        e isa NumericalHealthError || rethrow()
        return _pair_failed(b, sprint(showerror, e), spec, item, xi0)..., xi0
    end
    attrs = Dict{String, Any}("Dx" => Float64(ps.Dx), "Dx_r" => Float64(ps.Dx_r), "mode" => item.coordinate_mode, "xi0" => Float64.(collect(xi0)))
    return (z = ps.z, W = W), attrs, xi0
end

# ---------------------------------------------------------------- one item
function _ind_frame(ctx, item, ::Type{T}) where {T}
    W = get(ctx.frames, (item.precision, item.direction, 2), nothing)
    return W === nothing ? deviation_frame(T, item.direction, 2) : W
end

function _ind_tables(::Type{T}, item, spec, ctx, b, sys, x0) where {T}
    kind, N, tau, wound = spec["kind"], item.horizon, item.tau, ind_winding(item)
    red = ind_chain(T, item, spec, sys)
    attrs, secs = Dict{String, Any}(), 0.0
    if kind in ("tangent", "stretching")
        st0 = Tangent._init_states(red, T, Val(2), Val(2))
        nws = Tangent._chain_workspace(red, T, Val(2), Val(2))
        ws = nws == 0 ? Tangent._NO_WS : zeros(UInt32, nws)
        st = _ind_series_loop!(b, sys, x0, _ind_frame(ctx, item, T), red, st0, ws, tau, N, item.orbit, wound)
        kind == "tangent" && return _ind_nt(series = b.series, columns = ind_series_columns(spec)), attrs, secs
        b.failrow == 0 || return _ind_nt(series = b.series, columns = ["stretching.D"], scalars = fill(NaN, 18), scolumns = IND_ST_SCALARS, range = Float64[first(red).lo, first(red).hi]), attrs, secs
        sub = spec["counts"]
        store = sub == "P" ? item.design_class != "h4" : sub == "S_D" ? item.in_SD : sub == "S_L" ? item.in_SL : sub == "h4" ? item.design_class == "h4" : false
        return _stretch_tables(first(red), reducer_final(first(red), first(st), ws), b.grid, tau, b, store), attrs, secs
    end
    kind == "pair" && return _ind_pair(T, item, spec, ctx, b, sys)
    return _ind_window(T, item, spec, ctx, b, sys, x0, red)
end

function _ind_pair(::Type{T}, item, spec, ctx, b, sys) where {T}
    start, attrs, _ = _pair_setup(T, item, spec, ctx, b, sys)
    start === nothing && return _ind_nt(series = b.series, columns = ind_series_columns(spec)), attrs, 0.0
    red = ind_chain(T, item, spec, sys)
    st0 = Tangent._init_states(red, T, Val(4), Val(2))
    _ind_series_loop!(b, PairMap(sys), start.z, start.W, red, st0, Tangent._NO_WS, item.tau, item.horizon, item.orbit, ind_winding(item))
    return _ind_nt(series = b.series, columns = ind_series_columns(spec)), attrs, 0.0
end

function _ind_window(::Type{T}, item, spec, ctx, b, sys, x0, red) where {T}
    kind, tau, N, wound = spec["kind"], item.tau, item.horizon, ind_winding(item)
    W1 = SMatrix{2, 1, T}(one(T), zero(T))
    st0 = Tangent._init_states(red, T, Val(2), Val(1))
    attrs = Dict{String, Any}()
    if kind == "frequency"
        buf = Vector{Complex{T}}(undef, 2 * maximum(spec["windows"]))
        obs = Val(Symbol(spec["observable"]))
        _ind_window_loop!(b, sys, x0, W1, red, st0, tau, N, item.orbit, wound, (j, x) -> (buf[j] = ind_freq_sample(obs, x, T); nothing))
        b.failrow == 0 || return _failed_freq_tables(spec), attrs, 0.0
        ctx.hooks.series === nothing || ctx.hooks.series((precision = item.precision, spec = item.method), item.orbit, buf)
        t0 = time_ns()
        tab, failure = _freq_tables(spec, buf, item.orbit)
        isempty(failure) || (b.failrow = max(b.failrow, 1); b.msg = failure)
        return tab, attrs, (time_ns() - t0) / NS_PER_S
    end
    st = _ind_window_loop!(b, sys, x0, W1, red, st0, tau, N, item.orbit, wound, (j, x) -> nothing)
    b.failrow == 0 || return _failed_window_tables(spec, item), attrs, 0.0
    fin = reducer_final(last(red), last(st))
    kind == "phase_advance" && return _pa_tables(last(red), fin, spec["windows"]), attrs, 0.0
    return _bk_tables(last(red), fin, spec["windows"], ind_observable_names(spec, item)), attrs, 0.0
end

_failed_freq_tables(spec) = _ind_nt(final = fill(NaN, length(spec["windows"]), 11), rows = ["W=$w" for w in spec["windows"]], fcolumns = IND_FREQ_COLUMNS,
                                    flags = zeros(Int64, length(spec["windows"])), bits = IND_FREQ_BITS)

# Tables of a window item that failed before its run end: NaN values, no flag, `complete` false.
function _failed_window_tables(spec, item)
    H = length(spec["windows"])
    if spec["kind"] == "phase_advance"
        return _ind_nt(final = fill(NaN, H, 4), rows = ["W=$w" for w in spec["windows"]], fcolumns = IND_PA_COLUMNS, flags = zeros(Int64, H), bits = IND_PA_BITS,
                       scalars = Float64[0, 0, 0], scolumns = IND_PA_SCALARS)
    end
    names = ind_observable_names(spec, item)
    rows = ["W=$w $nm" for w in spec["windows"] for nm in names]
    return _ind_nt(final = fill(NaN, length(rows), 8), rows = rows, fcolumns = IND_BK_COLUMNS, flags = zeros(Int64, length(rows)), bits = IND_BK_BITS,
                   scalars = Float64[0, 0], scolumns = IND_BK_SCALARS)
end

# The record of one indicator item (E25, R15): grid rows of the orbit, the monitor maxima, the tables, the cost entry.
function _ind_record(::Type{T}, item, spec, b::IndBuf, tab, attrs0, secs, ctx) where {T}
    n = length(b.grid)
    M = SMatrix{2, 2, T}(T(b.mmax), zero(T), zero(T), zero(T))
    df = ctx.df
    thr = Dict{String, Float64}("det" => Float64(det_threshold(df["det_residual_tol"], T, M)), "symp" => Float64(symp_threshold(df["det_residual_tol"], T, M)))
    diag = Dict{String, Vector{Float64}}("det" => copy(b.det), "symp" => copy(b.symp))
    fe = Dict{String, Int}(nm => _first_exceedance(diag[nm], thr[nm], b.t) for nm in ("det", "symp"))
    attrs = Dict{String, Any}("system" => item.system, "params" => item.param, "ic" => item.x0 * "," * item.p0, "class_label" => item.design_class,
        "seed" => item.orbit_seed, "item_hash" => item.hash, "design_class" => item.design_class, "flags" => String[], "steps_recorded" => n)
    merge!(attrs, attrs0)
    cost = Float64[item.horizon, item.vectors, 0.0, 0.0, ind_item_bytes(ctx, item)]
    spec["kind"] == "frequency" && push!(cost, secs)
    wind = item.coordinate_mode == "wrapped" ? copy(b.winding) : nothing
    return OrbitRecord(item.orbit, attrs, copy(b.t), copy(b.state), wind, Dict{String, Matrix{Float64}}(), diag, thr, fe, Dict(item.method => cost), b.msg,
                       Dict{String, NamedTuple}(item.method => tab))
end

"Record one indicator item of the plan in its element type (called by `run_item`)."
function run_indicator_item(::Type{T}, item, ctx) where {T}
    spec = ctx.specs[item.method]
    sys, x0 = _system_start(T, item)
    b = _new_ind_buf(ctx.grids[(item.horizon, item.tau)], item.coordinate_mode == "wrapped", length(ind_series_columns(spec)))
    tab, attrs, secs = _ind_tables(T, item, spec, ctx, b, sys, x0)
    _blank_ind!(b)
    return _ind_record(T, item, spec, b, tab, attrs, secs, ctx)
end

# bytes(orbit) of an item (E28), computed once per run type in the sequential part (`prime_frames!`); a miss in a parallel region only recomputes.
_ind_bytes_key(it) = (it.method, it.precision, it.horizon, it.tau, it.coordinate_mode, it.system == "linear")
ind_item_bytes(ctx, it) = get(() -> orbit_bytes(ctx.plan, it).total, ctx.bytes, _ind_bytes_key(it))

# The frame checks of the members that have one (D67, R13 item 2), once per (frame, spec): a failure refuses the start naming the spec.
function check_item_frame!(ctx, it, W)
    spec = ctx.specs[it.method]
    tag = (:framecheck, it.precision, it.direction, it.vectors, it.method)
    tag in ctx.warmed && return nothing
    push!(ctx.warmed, tag)
    kind = spec["kind"]
    names = kind == "tangent" ? spec["members"] : kind == "stretching" ? ["stretching"] : kind == "pair" ? [m for m in spec["members"] if m == "MEGNO"] : String[]
    for m in names
        r = m == "stretching" ? StretchingNumbers(; range = (-1.0, 1.0), bins = 1, interval = 1, min_per_cell = 0) : _ind_member(m, true)
        try
            check_frame(r, W; orbit = 0)
        catch e
            e isa NumericalHealthError || rethrow()
            throw(ArgumentError("indicator spec `$(spec["id"])`, member $m, direction $(it.direction): " * sprint(showerror, e)))
        end
    end
    return nothing
end

"""
    hot_loop_bytes(spec; T = Float64) -> Int

Bytes allocated by the recording loop of one orbit of an indicator spec after warm-up (R14): zero for the tangent, pair, phase-advance, Birkhoff and
frequency kinds (the frequency transform is outside the loop); the stretching final builds its count matrices at every grid step and allocates.
"""
function hot_loop_bytes(spec::AbstractDict; T::Type = Float64)
    kind = spec["kind"]
    item = (system = "standard", param = "1.5", x0 = "1.0", p0 = "2.0", coordinate_mode = "lifted")
    sys, x0 = StandardMap(exact(T, 1.5)), SVector(exact(T, 1.0), exact(T, 2.0))
    red = ind_chain(T, item, spec, sys)
    L = ind_run_length(spec)
    tau, N = 1, L === nothing ? 64 : L
    grid = record_grid(N, tau; points_per_decade = 5, first_steps = 10, extra_fractions = [0.1, 0.5])
    b = _new_ind_buf(grid, false, length(ind_series_columns(spec)))
    if kind == "pair"
        ps = pair_start(sys, x0, SVector{2, T}(exact(T, 0.6), exact(T, 0.8)), 1e-8)
        z = zero(T)
        u1, u2 = exact(T, 0.6), exact(T, 0.8)
        run = () -> _ind_series_loop!(b, PairMap(sys), ps.z, SMatrix{4, 2, T}(u1, u2, z, z, z, z, u1, u2), red, Tangent._init_states(red, T, Val(4), Val(2)), Tangent._NO_WS, tau, N, 1, false)
    elseif ind_window_kind(kind)
        W1, buf = SMatrix{2, 1, T}(one(T), zero(T)), Vector{Complex{T}}(undef, 2 * (kind == "frequency" ? maximum(spec["windows"]) : 1))
        feed = kind == "frequency" ? (j, x) -> (buf[j] = ind_freq_sample(Val(Symbol(spec["observable"])), x, T); nothing) : (j, x) -> nothing
        run = () -> _ind_window_loop!(b, sys, x0, W1, red, Tangent._init_states(red, T, Val(2), Val(1)), tau, N, 1, false, feed)
    else
        nws = Tangent._chain_workspace(red, T, Val(2), Val(2))
        ws = nws == 0 ? Tangent._NO_WS : zeros(UInt32, nws)
        W0 = deviation_frame(T, "pi/4", 2)
        run = () -> _ind_series_loop!(b, sys, x0, W0, red, Tangent._init_states(red, T, Val(2), Val(2)), ws, tau, N, 1, false)
    end
    run()
    return @allocated run()
end
