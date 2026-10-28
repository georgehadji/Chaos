# Device path of the indicator recorder (Amendment D86, R16): the tangent, pair, phase-advance and Birkhoff kinds run the shared per-orbit body in
# step-chunked launches of the one KernelAbstractions kernel per family and write the grid rows into the record array R[col, grid step, orbit]; the window
# tables are read on the host from the carried reducer state. Float64 per record; stretching and frequency are CPU only. Kernel-path rule: nothing
# reachable from the kernels throws with a built string, allocates or has unstable inference; every argument check is in the block, the plan and the host.

struct DevRow{A}
    R::A
    g::Int
    i::Int
end
@inline _set!(d::DevRow, j, v) = (d.R[j, d.g, d.i] = v; nothing)

# Row of a series kind: the columns of E24, then det, symp, M_max of the monitor, then the first `ns` components of the orbit point.
@inline function _dev_ind_row!(R, i, g, x, red::Tuple, st::Tuple, ns)
    off = _put_members!(DevRow(R, g, i), 0, Base.front(red), Base.front(st), nothing)
    _dev_mon_state!(R, i, g, off, x, last(red), last(st), ns)
    return nothing
end

@inline function _dev_mon_state!(R, i, g, off, x, mon, mst, ns)
    f = reducer_final(mon, mst)
    R[off + 1, g, i] = Float64(f.det_max)
    R[off + 2, g, i] = Float64(f.symp_max)
    R[off + 3, g, i] = Float64(f.M_max)
    for j in 1:ns
        R[off + 3 + j, g, i] = Float64(x[j])
    end
    return nothing
end

@kernel function _ind_series_kernel!(R, fail, x, W, st, code, hstep, sys, mem, grid, glo, ghi, n0, steps, N, tau, ns)
    i = @index(Global)
    if n0 == 0 || code[i] == Tangent.STATUS_OK
        s = sys[i]
        red = (mem..., RunMonitor(s))
        xi, Wi = x[i], W[i]
        st0 = n0 == 0 ? Tangent._init_states(red, eltype(xi), Val(size(Wi, 1)), Val(size(Wi, 2))) : st[i]
        n, c, h = n0, Tangent.STATUS_OK, n0
        for g in glo:ghi
            target = grid[g]
            xi, Wi, st0, c, h = Tangent._orbit_body(s, xi, Wi, st0, n, target - n, N, tau, red)
            if c != Tangent.STATUS_OK
                fail[i] = g
                break
            end
            n = target
            _dev_ind_row!(R, i, g, xi, red, st0, ns)
        end
        if c == Tangent.STATUS_OK && n < n0 + steps
            xi, Wi, st0, c, h = Tangent._orbit_body(s, xi, Wi, st0, n, n0 + steps - n, N, tau, red)
            c == Tangent.STATUS_OK || (fail[i] = ghi + 1)
        end
        x[i] = xi
        W[i] = Wi
        st[i] = st0
        code[i] = c
        hstep[i] = h
    end
end

# The window kinds: one handoff at a time with the constant deviation column (as the host loop), the monitor and the orbit point at the grid steps.
@kernel function _ind_window_kernel!(R, fail, x, st, code, hstep, sys, rd, grid, glo, ghi, n0, steps, N, tau, ns)
    i = @index(Global)
    if n0 == 0 || code[i] == Tangent.STATUS_OK
        s = sys[i]
        red = (RunMonitor(s), rd...)
        xi = x[i]
        T = eltype(xi)
        W1 = SMatrix{2, 1, T}(one(T), zero(T))
        st0 = n0 == 0 ? Tangent._init_states(red, T, Val(2), Val(1)) : st[i]
        n, g, c, h = n0, glo, Tangent.STATUS_OK, n0
        for _ in 1:(steps ÷ tau)
            xi, Wd, st0, c, h = Tangent._orbit_body(s, xi, W1, st0, n, tau, N, tau, red)
            if c != Tangent.STATUS_OK
                fail[i] = g <= ghi ? g : ghi + 1
                break
            end
            n += tau
            if g <= ghi && n == grid[g]
                _dev_mon_state!(R, i, g, 0, xi, first(red), first(st0), ns)
                g += 1
            end
        end
        x[i] = xi
        st[i] = st0
        code[i] = c
        hstep[i] = h
    end
end

# ---------------------------------------------------------------- host side
# Start data of an item for the device: the stream map, the start state and the deviation frame; `failed` is the record of an item refused at its start (pair).
function _ind_setup(::Type{T}, item, spec, ctx) where {T}
    sys, x0 = _system_start(T, item)
    kind = spec["kind"]
    if kind == "pair"
        b = _new_ind_buf(ctx.grids[(item.horizon, item.tau)], item.coordinate_mode == "wrapped", length(ind_series_columns(spec)))
        start, attrs, _ = _pair_setup(T, item, spec, ctx, b, sys)
        start === nothing && return (failed = run_indicator_item(T, item, ctx),)
        return (sys = PairMap(sys), base = sys, x = start.z, W = start.W, attrs = attrs, failed = nothing)
    end
    W = ind_window_kind(kind) ? SMatrix{2, 1, T}(one(T), zero(T)) : _ind_frame(ctx, item, T)
    return (sys = sys, base = sys, x = x0, W = W, attrs = Dict{String, Any}(), failed = nothing)
end

# Orbits per device batch (R3 item 3, R17 item 1): the real reducer bytes of the chain, not an allowance.
function _ind_device_capacity(ctx::RunCtx, it, ::Type{T}, free::Integer) where {T}
    b = orbit_bytes(ctx.plan, it)
    return launch_capacity(free; f = Tangent.CHUNK_HEADROOM, M = b.M, reducer = b.R, columns = b.c, grid_steps = b.ng)
end

# Host buffer of orbit j from the record array: the rows before the failing row.
function _ind_unpack!(b::IndBuf, Rh, j, nser, ns, wound, fail)
    b.failrow = fail
    last = fail == 0 ? length(b.grid) : fail - 1
    for g in 1:last
        b.t[g] = b.grid[g]
        for c in 1:nser
            b.series[g, c] = Rh[c, g, j]
        end
        b.det[g], b.symp[g] = Rh[nser + 1, g, j], Rh[nser + 2, g, j]
        b.state[g, 1], b.state[g, 2] = Rh[nser + 4, g, j], Rh[nser + 5, g, j]
        wound && (b.winding[g, 1] = _toint(Rh[nser + 6, g, j]); b.winding[g, 2] = _toint(Rh[nser + 7, g, j]))
    end
    last >= 1 && (b.mmax = Rh[nser + 3, last, j])
    return b
end

function _ind_tables_device(::Type{T}, item, spec, b, rd, sth_j) where {T}
    kind = spec["kind"]
    kind in ("tangent", "pair") && return _ind_nt(series = b.series, columns = ind_series_columns(spec))
    b.failrow == 0 || return _failed_window_tables(spec, item)
    fin = reducer_final(last(rd), last(sth_j))
    kind == "phase_advance" && return _pa_tables(last(rd), fin, spec["windows"])
    return _bk_tables(last(rd), fin, spec["windows"], ind_observable_names(spec, item))
end

# The launch closures of the two families, each with its own kernel (a kernel variable shared by two closures is flagged by JET).
function _ind_window_call(backend, R, fail, x, st, code, hstep, dsys, rds, dgrid, N, tau, ns, m)
    kern = _ind_window_kernel!(backend, Tangent.workgroup(backend))
    return (glo, ghi, n0, steps) -> kern(R, fail, x, st, code, hstep, dsys, rds, dgrid, glo, ghi, n0, steps, N, tau, ns; ndrange = m), (x = x, states = st)
end

function _ind_series_call(backend, R, fail, x, W0, st, code, hstep, dsys, rds, dgrid, N, tau, ns, m)
    W = _todev(backend, W0)
    kern = _ind_series_kernel!(backend, Tangent.workgroup(backend))
    return (glo, ghi, n0, steps) -> kern(R, fail, x, W, st, code, hstep, dsys, rds, dgrid, glo, ghi, n0, steps, N, tau, ns; ndrange = m), (x = x, W = W, states = st)
end

function _ind_batch!(recs, ctx, ::Type{T}, items, spec, setups, r, backend, L, ci, hook) where {T}
    it = items[first(r)]
    kind, m, N, tau = spec["kind"], length(r), it.horizon, it.tau
    grid, wound, wrapped = ctx.grids[(N, tau)], ind_winding(it), it.coordinate_mode == "wrapped"
    sysv, x0, W0 = [setups[j].sys for j in r], [setups[j].x for j in r], [setups[j].W for j in r]
    base = setups[first(r)].base
    rds = ind_readers(T, it, spec, base)                       # members of the chain, or the one reader of a window kind
    red0 = ind_window_kind(kind) ? (ind_monitor(base, kind), rds...) : (rds..., ind_monitor(base, kind))
    nser, ns = length(ind_series_columns(spec)), wound ? 4 : 2
    s0 = Tangent._init_states(red0, T, Val(size(W0[1], 1)), Val(size(W0[1], 2)))
    x, st = _todev(backend, x0), KernelAbstractions.allocate(backend, typeof(s0), m)
    code, hstep, fail = KernelAbstractions.zeros(backend, Int8, m), KernelAbstractions.zeros(backend, Int, m), KernelAbstractions.zeros(backend, Int, m)
    dsys, dgrid = _todev(backend, sysv), _todev(backend, grid)
    R = KernelAbstractions.zeros(backend, Float64, nser + 3 + ns, length(grid), m)
    call, carried = ind_window_kind(kind) ? _ind_window_call(backend, R, fail, x, st, code, hstep, dsys, rds, dgrid, N, tau, ns, m) :
                    _ind_series_call(backend, R, fail, x, W0, st, code, hstep, dsys, rds, dgrid, N, tau, ns, m)
    _launches!(call, carried, backend, grid, N, L, ci, hook)
    Rh, fh, ch, hh, xh, sth = Tangent._host(R), Tangent._host(fail), Tangent._host(code), Tangent._host(hstep), Tangent._host(x), Tangent._host(st)
    for (jj, j) in enumerate(r)
        b = _ind_unpack!(_new_ind_buf(grid, wrapped, nser), Rh, jj, nser, ns, wound, fh[jj])
        red = ind_window_kind(kind) ? (ind_monitor(setups[j].base, kind), rds...) : (rds..., ind_monitor(setups[j].base, kind))
        fh[jj] == 0 || (b.msg = _failure_text(red, ch[jj], sth[jj], xh[jj], W0[jj], hh[jj], items[j].orbit))
        _blank_ind!(b)
        recs[j] = _ind_record(T, items[j], spec, b, _ind_tables_device(T, items[j], spec, b, rds, sth[jj]), setups[j].attrs, 0.0, ctx)
    end
    return nothing
end

"Device records of the items of one indicator run type (called by `device_records`); same contract: `recs` is `nothing` when the group falls back."
function ind_device_records(ctx::RunCtx, items; hook::Bool = true)
    prime_frames!(ctx, items)
    it = first(items)
    spec = ctx.specs[it.method]
    cpu_only_method(it.method, spec["kind"]) && return (recs = nothing, prov = cpu_prov(FALLBACK_REASONS[6]))
    T, bits = parse_precision(it.precision)
    bits === nothing || return (recs = nothing, prov = cpu_prov(FALLBACK_REASONS[3]))
    L, cap = Tangent._device_plan(Tangent._forced_launch(it.tau), T, it.horizon, it.tau)
    L == 0 && return (recs = nothing, prov = cpu_prov(FALLBACK_REASONS[8]))
    recs, nbatch, csize = Vector{OrbitRecord}(undef, length(items)), 0, 0
    try
        dev = Foundation.DEVICE[]
        backend = Tangent.device_backend(dev)
        csize = min(_ind_device_capacity(ctx, it, T, Foundation.free_memory(dev)), cap)
        csize >= 1 || return (recs = nothing, prov = cpu_prov(FALLBACK_REASONS[4]))
        setups = [_ind_setup(T, i, spec, ctx) for i in items]
        for (j, st) in enumerate(setups)
            st.failed === nothing || (recs[j] = st.failed)
        end
        ok = [j for j in eachindex(items) if setups[j].failed === nothing]
        starts = [(setups[j].base, nothing) for j in ok]
        for idx in _by_system(starts), r in Iterators.partition(idx, csize)
            nbatch += 1
            _ind_batch!(recs, ctx, T, items, spec, setups, [ok[k] for k in r], backend, L, nbatch, hook)
        end
        return (recs = recs, prov = (backend = :gpu, device = Foundation.probe(dev, T).device, reason = "", gpu_chunks = nbatch, chunk_size = csize,
                                     launch_steps = L, launches = cld(it.horizon, L)))
    catch e
        why = Tangent._fallback_reason(e)
        why == 0 && rethrow()
        why == 9 && (Foundation.DEVICE[] = Foundation.NoDevice())
        return (recs = nothing, prov = cpu_prov(FALLBACK_REASONS[why]))
    end
end
