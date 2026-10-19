# Device path of the recorder (return of D53 to stage 03): the kernel-path methods run the shared per-orbit body on the device in
# step-chunked launches (D33, D44) and write the grid records into preallocated per-orbit arrays, so a launch boundary may fall
# between grid steps. One KernelAbstractions kernel serves the CPU backend (the package tests) and CUDA. L6 (library QR, not
# isbits) stays on the CPU path of run_record.jl, which is also the reference of every device record (backend_agreement).

# Columns of the record array R[col, grid step, orbit]: lambda (k), det, symp, orth, sum, pair, G_max, M_max, then the state.
_ncol(k::Integer, s::Integer) = k + 7 + s

# One grid row of an orbit: the same quantities as `_store!`, written into R. The running maxima of sum and pairing read the row before.
@inline function _dev_row!(R, i, g, k, x, est, mon)
    lam = est.lambda
    for j in 1:k
        R[j, g, i] = Float64(lam[j])
    end
    R[k + 1, g, i] = Float64(mon.det_max)
    R[k + 2, g, i] = Float64(mon.symp_max)
    R[k + 3, g, i] = _orth(est)
    _dev_sumpair!(R, k, i, g, lam)
    R[k + 6, g, i] = _gmax(est)
    R[k + 7, g, i] = Float64(mon.M_max)
    _dev_state!(R, k, i, g, x)
    return nothing
end

@inline function _dev_state!(R, k, i, g, x)
    for j in 1:length(x)
        R[k + 7 + j, g, i] = Float64(x[j])
    end
    return nothing
end

@inline function _dev_sumpair!(R, k, i, g, lam::SVector{2})
    s, p = Float64(eps_sum(lam)), Float64(eps_pair(lam))
    R[k + 4, g, i] = g > 1 ? max(R[k + 4, g - 1, i], s) : s
    R[k + 5, g, i] = g > 1 ? max(R[k + 5, g - 1, i], p) : p
    return nothing
end
@inline _dev_sumpair!(R, k, i, g, lam::SVector{1}) = nothing

# One launch of the variational methods: for each orbit the grid steps in (n0, n0 + steps] in segments of the shared body, then the
# tail up to the end of the launch. A failure is recorded as the grid row whose segment holds it (`fail`) and the orbit is not
# advanced in later launches; the state, matrix, reducer states, status and failing step are carried in place.
@kernel function _record_kernel!(R, fail, x, W, st, code, hstep, sys, est, grid, glo, ghi, n0, steps, N, tau, k)
    i = @index(Global)
    if n0 == 0 || code[i] == Tangent.STATUS_OK
        s = sys[i]
        red = (est, RunMonitor(s))
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
            _dev_row!(R, i, g, k, xi, reducer_final(est, st0[1]), reducer_final(red[2], st0[2]))
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

# The same for L2 (two-trajectory): the carried end record replaces the state, the matrix and the reducer states.
@kernel function _shadow_record_kernel!(R, fail, ends, x0, u0, sys, d0, grid, glo, ghi, n0, steps, N, tau, k)
    i = @index(Global)
    if n0 == 0 || ends[i].code == 0
        s = sys[i]
        e = n0 == 0 ? Lyapunov._blank_end(x0[i], d0) : ends[i]
        n = n0
        for g in glo:ghi
            target = grid[g]
            e = Lyapunov._shadow_body(s, x0[i], u0[i], e, n, target - n, N, tau, d0)
            if e.code != 0
                fail[i] = g
                break
            end
            n = target
            R[1, g, i] = Float64(total(e.sum) / target)
            _dev_state!(R, k, i, g, e.x)
        end
        if e.code == 0 && n < n0 + steps
            e = Lyapunov._shadow_body(s, x0[i], u0[i], e, n, n0 + steps - n, N, tau, d0)
            e.code == 0 || (fail[i] = ghi + 1)
        end
        ends[i] = e
    end
end

# ---------------------------------------------------------------- host side
"Provenance of a group that ran on the CPU; `reason` is why it did not run on the device (empty when the CPU was asked for)."
cpu_prov(reason::AbstractString = "") = (backend = :cpu, device = "cpu", reason = String(reason), gpu_chunks = 0, chunk_size = 0, launch_steps = 0, launches = 0)

_todev(backend, v::AbstractVector) = (a = KernelAbstractions.allocate(backend, eltype(v), length(v)); copyto!(a, v); a)

# Orbits per device batch from the free device memory (R3 item 3 with the device figures): resident bytes, reducer allowance, record array.
function _device_capacity(ctx::RunCtx, it, ::Type{T}, free::Integer) where {T}
    k, s = it.vectors, it.coordinate_mode == "wrapped" ? 4 : 2
    M = Tangent.bytes_per_orbit(Val(2), Val(k), T, s)
    return launch_capacity(free; f = Tangent.CHUNK_HEADROOM, M = M, reducer = 256, columns = _ncol(k, s), grid_steps = length(ctx.grids[(it.horizon, it.tau)]))
end

# Launches of one batch: `call(glo, ghi, n0, steps)` starts the kernel on the grid steps inside the launch; the hook of tangent-core runs before each.
function _launches!(call, carried, backend, grid, N::Int, L::Int, ci::Int, hook::Bool)
    n0, launch = 0, 0
    while n0 < N
        steps = min(L, N - n0)
        launch += 1
        hook && Tangent._before_launch(Foundation.DEVICE[], ci, launch, carried)
        call(searchsortedfirst(grid, n0 + 1), searchsortedlast(grid, n0 + steps), n0, steps)
        KernelAbstractions.synchronize(backend)
        n0 += steps
    end
    return launch
end

# Host record buffer of orbit j from the record array: rows before the failing row, the running G_max and M_max of the last good row.
function _unpack!(b::RecBuf, Rh, j, k, s, fail)
    b.failrow = fail
    last = fail == 0 ? length(b.grid) : fail - 1
    for g in 1:last
        b.t[g] = b.grid[g]
        for c in 1:k
            b.lyap[g, c] = Rh[c, g, j]
        end
        b.det[g], b.symp[g], b.orth[g], b.sumd[g], b.pair[g] = Rh[k + 1, g, j], Rh[k + 2, g, j], Rh[k + 3, g, j], Rh[k + 4, g, j], Rh[k + 5, g, j]
        b.state[g, 1], b.state[g, 2] = Rh[k + 8, g, j], Rh[k + 9, g, j]
        s == 4 && size(b.winding, 1) > 0 && (b.winding[g, 1] = _toint(Rh[k + 10, g, j]); b.winding[g, 2] = _toint(Rh[k + 11, g, j]))
    end
    last >= 1 && (b.gmax = Rh[k + 6, last, j]; b.mmax = Rh[k + 7, last, j])
    return b
end

function _failure_text(red, code, st, x, W, step, orbit)
    err = code == Tangent.STATUS_REDUCER ? Tangent._reducer_error(red, st, orbit, :gpu) : Tangent._health_error(code, x, W, step, orbit, :gpu)
    return sprint(showerror, err)
end

# State arrays of a batch on the device, the host start data and the launch loop of the variational kernel.
function _variational_batch!(recs, ctx, ::Type{T}, items, starts, r, backend, L, ci, hook) where {T}
    it = items[first(r)]
    m, k, N, tau = length(r), it.vectors, it.horizon, it.tau
    grid, wrapped = ctx.grids[(N, tau)], it.coordinate_mode == "wrapped"
    sysv, x0 = [starts[j][1] for j in r], [starts[j][2] for j in r]
    W0 = ctx.frames[(it.precision, it.direction, k)]
    est = _estimator(it.method, ctx.df["orthogonality_tol"])
    s0 = Tangent._init_states((est, RunMonitor(sysv[1])), T, Val(2), Val(k))
    x, W = _todev(backend, x0), _todev(backend, fill(W0, m))
    st, code, hstep = KernelAbstractions.allocate(backend, typeof(s0), m), KernelAbstractions.zeros(backend, Int8, m), KernelAbstractions.zeros(backend, Int, m)
    fail, dsys, dgrid = KernelAbstractions.zeros(backend, Int, m), _todev(backend, sysv), _todev(backend, grid)
    R = KernelAbstractions.zeros(backend, Float64, _ncol(k, length(x0[1])), length(grid), m)
    kern = _record_kernel!(backend, Tangent.workgroup(backend))
    _launches!((glo, ghi, n0, steps) -> kern(R, fail, x, W, st, code, hstep, dsys, est, dgrid, glo, ghi, n0, steps, N, tau, k; ndrange = m),
               (x = x, W = W, states = st), backend, grid, N, L, ci, hook)
    Rh, fh, ch, hh, xh, Wh, sth = Tangent._host(R), Tangent._host(fail), Tangent._host(code), Tangent._host(hstep), Tangent._host(x), Tangent._host(W), Tangent._host(st)
    for (jj, j) in enumerate(r)
        red = (est, RunMonitor(sysv[jj]))
        b = _unpack!(_new_buf(grid, k, wrapped), Rh, jj, k, length(x0[1]), fh[jj])
        fh[jj] == 0 || (b.msg = _failure_text(red, ch[jj], sth[jj], xh[jj], Wh[jj], hh[jj], items[j].orbit))
        _blank_failed!(b)
        recs[j] = _orbit_record(T, items[j], b, Tangent._init_states(red, T, Val(2), Val(k)), ctx.df)
    end
    return nothing
end

function _shadow_batch!(recs, ctx, ::Type{T}, items, starts, r, backend, L, ci, hook) where {T}
    it = items[first(r)]
    m, N, tau = length(r), it.horizon, it.tau
    grid, wrapped = ctx.grids[(N, tau)], it.coordinate_mode == "wrapped"
    sysv, x0 = [starts[j][1] for j in r], [starts[j][2] for j in r]
    W0 = ctx.frames[(it.precision, it.direction, 1)]
    u0, d0 = SVector{2, T}(W0[1, 1], W0[2, 1]), exact(T, it.eps)
    ends = KernelAbstractions.allocate(backend, Lyapunov.ShadowEnd{length(x0[1]), T}, m)
    fail, dsys, dgrid = KernelAbstractions.zeros(backend, Int, m), _todev(backend, sysv), _todev(backend, grid)
    dx0, du0 = _todev(backend, x0), _todev(backend, fill(u0, m))
    R = KernelAbstractions.zeros(backend, Float64, _ncol(1, length(x0[1])), length(grid), m)
    kern = _shadow_record_kernel!(backend, Tangent.workgroup(backend))
    _launches!((glo, ghi, n0, steps) -> kern(R, fail, ends, dx0, du0, dsys, d0, dgrid, glo, ghi, n0, steps, N, tau, 1; ndrange = m),
               (ends = ends,), backend, grid, N, L, ci, hook)
    Rh, fh, eh = Tangent._host(R), Tangent._host(fail), Tangent._host(ends)
    for (jj, j) in enumerate(r)
        b = _unpack!(_new_buf(grid, 1, wrapped), Rh, jj, 1, length(x0[1]), fh[jj])
        fh[jj] == 0 || (b.msg = sprint(showerror, Lyapunov._shadow_error(eh[jj], d0, items[j].orbit, :gpu)))
        _blank_failed!(b)
        recs[j] = _orbit_record(T, items[j], b, (), ctx.df)
    end
    return nothing
end

# Start data of the items, grouped by system type (the Y5 control is a linear map next to the standard maps), in first-seen order.
function _by_system(starts)
    types = unique(typeof(s[1]) for s in starts)
    return [findall(s -> typeof(s[1]) === S, starts) for S in types]
end

"""
    device_records(ctx, items; hook = true) -> (recs, prov)

Record the work items of one run type on the device in force (`Foundation.DEVICE[]`), in batches of the capacity that the free device memory allows
and in launches of the step-chunked budget (D44). `recs` is `nothing` when the group falls back to the CPU; `prov` then holds the reason (one of the
nine) and `backend = :cpu`. Otherwise `prov` holds the tangent-core provenance fields: `device`, `gpu_chunks`, `chunk_size`, `launch_steps`, `launches`.
`hook = false` skips the test hook of tangent-core (the warm-up launch).
"""
function device_records(ctx::RunCtx, items; hook::Bool = true)
    prime_frames!(ctx, items)
    it = first(items)
    T, bits = parse_precision(it.precision)
    bits === nothing || return (recs = nothing, prov = cpu_prov(FALLBACK_REASONS[3]))
    L, cap = Tangent._device_plan(Tangent._forced_launch(it.tau), T, it.horizon, it.tau)
    L == 0 && return (recs = nothing, prov = cpu_prov(FALLBACK_REASONS[8]))
    recs, nbatch, csize = Vector{OrbitRecord}(undef, length(items)), 0, 0
    try
        dev = Foundation.DEVICE[]
        backend = Tangent.device_backend(dev)
        csize = min(_device_capacity(ctx, it, T, Foundation.free_memory(dev)), cap)
        csize >= 1 || return (recs = nothing, prov = cpu_prov(FALLBACK_REASONS[4]))   # not one orbit fits: out of memory, not an error
        starts = [_system_start(T, i) for i in items]
        for idx in _by_system(starts), r in Iterators.partition(idx, csize)
            nbatch += 1
            (it.method == "L2" ? _shadow_batch! : _variational_batch!)(recs, ctx, T, items, starts, collect(r), backend, L, nbatch, hook)
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

"""
    compute_group(ctx, items, ntask, row, lost; times = nothing) -> (recs, prov)

Records of one launch group on the backend its (precision, method) row resolved to: the device when the row says so and the device is still
usable, else the CPU threads. A device failure (out of memory, lost, compile error, launch over the budget) sends this group to the CPU with
the reason recorded; a lost device (`lost[]`) keeps every later group there. `times` receives the per-item seconds as in `compute_records`.
"""
function compute_group(ctx::RunCtx, items, ntask::Integer, row, lost::Base.RefValue{Bool}; times = nothing)
    row.backend === :gpu || return compute_records(ctx, items, ntask; times = times), cpu_prov(row.reason)
    if lost[] || Foundation.DEVICE[] isa Foundation.NoDevice
        lost[] = true
        return compute_records(ctx, items, ntask; times = times), cpu_prov(FALLBACK_REASONS[9])
    end
    t0 = time_ns()
    g = device_records(ctx, items)
    if g.recs === nothing
        g.prov.reason == FALLBACK_REASONS[9] && (lost[] = true)
        return compute_records(ctx, items, ntask; times = times), g.prov
    end
    times === nothing || (times .= (time_ns() - t0) / NS_PER_S / length(items))
    return g.recs, g.prov
end
