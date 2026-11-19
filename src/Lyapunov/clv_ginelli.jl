# Covariant Lyapunov vectors by the method of Ginelli (catalog L12): a forward pass of the delivered discrete QR that records the R factors,
# then a backward pass of triangular solves with column normalisation (V6, V8). Spec: clv-ginelli V1 to V23, Algorithm 1 to 7.
# The module owns a per-orbit loop (not a reducer of the tangent stream): the orbit is advanced by `forward(map, z)` and the tangent of each
# step is `tangent(map, z)`, one call of each per step, tau steps per handoff. CPU only (P9). The backward half is in clv_backward.jl.

const CLV_STREAM_LABEL = 0x434c56434f454646          # distinguishes the stream of the starting coefficients from the frame stream (P13)

"""
    ClvGinelli(map; tau = 1, h = 1.0, orthonormalizer = HouseholderQR(), storage = :store, block = 0)

Covariant-Lyapunov-vector estimator (Ginelli's method). `tau >= 1` map steps per handoff; `h` the time per map step (1 for a map), used only to give
the exponents per unit time (`Dt = tau h`); `orthonormalizer` any delivered [`Orthonormalizer`](@ref) of the forward phase; `storage` is `:store`
(all R factors in memory) or `:checkpoint` (state and frame at the start of each block of `block >= 1` handoffs, the blocks recomputed backwards).
`tau < 1`, an `h` that is not finite and positive, an unknown `storage` or `:checkpoint` with `block < 1` throw `ArgumentError`. `h` is kept
as given and converted to the element type of the run, so a `BigFloat` or `Rational` step keeps its precision (Amendment D124).
"""
struct ClvGinelli{M, A <: Orthonormalizer, H <: Real}
    map::M
    tau::Int
    h::H                                          # kept as given, converted to the element type at use (Amendment D124, M2)
    alg::A
    storage::Symbol
    block::Int
end

function ClvGinelli(map; tau::Integer = 1, h::Real = 1.0, orthonormalizer::Orthonormalizer = HouseholderQR(), storage::Symbol = :store, block::Integer = 0)
    tau >= 1 || throw(ArgumentError("tau must be at least 1, got $tau"))
    (isfinite(h) && h > 0) || throw(ArgumentError("h must be finite and positive, got $h"))   # Amendment D124 (M1)
    storage in (:store, :checkpoint) || throw(ArgumentError("storage must be :store or :checkpoint, got :$storage"))
    storage === :checkpoint && block < 1 && throw(ArgumentError("the checkpoint block must be at least 1, got $block"))
    return ClvGinelli(map, Int(tau), h, orthonormalizer, storage, Int(block))
end

"""
    ClvResult

Result of [`propagate_clv_orbit`](@ref); every array holds the samples n = 0 .. n_s at index n + 1. `C` (k, k, n_s + 1) the upper-triangular coefficient
matrices with unit columns; `lambda` (k, n_s + 1) the local exponents ln Lambda / (tau h), the last column NaN when `n_b = 0`; `cosines` (k, k, n_s + 1)
the matrices `|C' C|`; `V`, `G` (d, k, n_s + 1) the phase-space CLVs and the Gram-Schmidt frames, `R` (k, k, n_s + 1) the factors (last slice NaN when
`n_b = 0`), each `nothing` unless requested in `outputs`; `angles` (min(m1, m2), n_s + 1) the principal angles of the two `angle_sets`, ascending, else
`nothing`. With a `sink` all of these are `nothing`. `lambda_fwd` (k) the forward exponents of the window per unit time, bitwise those of `QRSpectrum`
divided by `h`; `x0` the state at the window start, `x` the state after all handoffs; `bytes` the memory formula [`clv_bytes`](@ref) of the run;
`record` the run parameters.
"""
struct ClvResult{T, X <: SVector}
    C::Union{Nothing, Array{T, 3}}
    lambda::Union{Nothing, Matrix{T}}
    cosines::Union{Nothing, Array{T, 3}}
    V::Union{Nothing, Array{T, 3}}
    G::Union{Nothing, Array{T, 3}}
    R::Union{Nothing, Array{T, 3}}
    angles::Union{Nothing, Matrix{T}}
    lambda_fwd::Vector{T}
    x0::X
    x::X
    bytes::Int
    record::NamedTuple
end

"""
    clv_initial_coefficients(seed, orbit, Val(k), T) -> SMatrix{k,k,T}

Starting backward matrix C: upper triangular, positive diagonal, unit columns (to a few eps of `T`), column 1 exactly `e_1`, from a stream keyed by
`(seed, orbit)` and labelled differently from the frame stream; column j does not depend on k (spec P13).
"""
function clv_initial_coefficients(seed::Integer, orbit::Integer, ::Val{k}, ::Type{T}) where {k, T}
    rng = StableRNG(_mix(_mix(seed % UInt64) + orbit % UInt64) ⊻ CLV_STREAM_LABEL)
    M = zeros(T, k, k)
    for j in 1:k
        v = randn(rng, j)
        v[j] = abs(v[j]) + 0.1
        col = T.(v)
        M[1:j, j] = col / norm(col)
    end
    M[1, 1] = one(T)
    return SMatrix{k, k, T}(M)
end

"""
    clv_bytes(Val(d), Val(k), T; n_s, n_b, storage = :store, block = 0, frames = true) -> Int

Memory of one orbit (spec V17): store-all `sizeof(T) (N k(k+1)/2 + (n_s + 1) d k) + bytes_work` with `N = n_s + n_b` (`frames = false` drops the frame
term); checkpointed `sizeof(T) (n_c (d + d k) + h_b (k(k+1)/2 + d k) + k(k+1)/2) + bytes_work`, `n_c = ceil(N / h_b)`;
`bytes_work = sizeof(T) (d + 2 d k + k(k+1)/2 + 2 k) + 16`. The result arrays handed to the caller are not counted.
"""
function clv_bytes(::Val{d}, ::Val{k}, ::Type{T}; n_s::Integer, n_b::Integer, storage::Symbol = :store, block::Integer = 0, frames::Bool = true) where {d, k, T}
    sT = sizeof(T); N = n_s + n_b; kk = k * (k + 1) ÷ 2; work = sT * (d + 2d * k + kk + 2k) + 16
    storage === :store && return sT * (N * kk + (frames ? (n_s + 1) * d * k : 0)) + work
    return sT * (cld(N, block) * (d + d * k) + block * (kk + d * k) + kk) + work
end

# ---------------------------------------------------------------------------------------------------------------- forward half
function _factor_error(R, status, column, step, orbit)
    status == QR_PASSES && return NumericalHealthError("orthonormalisation breakdown", Float64(MAX_PASSES + 1), Float64(MAX_PASSES), step, orbit, :cpu, column)
    pivot = Float64(R[column, column])
    return NumericalHealthError("orthonormalisation breakdown", pivot, isfinite(pivot) ? 0.0 : Inf, step, orbit, :cpu, column)
end

# One handoff: tau steps of the tangent and the orbit (tangent at the old point, as the delivered stream), the health checks, the forward QR.
@inline function _clv_handoff(clv::ClvGinelli, z, G::SMatrix{d, k, T}, step::Int, orbit::Int) where {d, k, T}
    W = G
    for _ in 1:clv.tau
        W = fused_product(tangent(clv.map, z), W)
        z = forward(clv.map, z)
    end
    all(isfinite, z) || throw(NumericalHealthError("state", NaN, 0.0, step, orbit, :cpu))
    all(isfinite, W) || throw(NumericalHealthError("deviation matrix", NaN, 0.0, step, orbit, :cpu))
    Q, R, _, status, column = orthonormalize(clv.alg, W)
    status == QR_OK || throw(_factor_error(R, status, column, step, orbit))
    if checks_conditioning(clv.alg)
        gw = Float64(window_Gw(R) * effective_eps(T))
        gw < 0.5 || throw(NumericalHealthError("orthonormalisation breakdown", gw, 0.5, step, orbit, :cpu))
    end
    return z, Q, R
end

@inline _addlog(acc::NTuple{k}, R::SMatrix{k, k}) where {k} = ntuple(i -> add(acc[i], log(R[i, i])), Val(k))

# The forward transient: n_f handoffs, nothing stored.
function _transient(clv::ClvGinelli, z, G, nf::Int, orbit::Int)
    for i in 1:nf
        z, G, _ = _clv_handoff(clv, z, G, i, orbit)
    end
    return z, G
end

# Store-all (V17): R of all N intervals packed, the frames n = 0 .. n_s when wanted; then one backward sweep. Returns (x_window, x_end, acc).
function _run_store!(o::ClvOut{T}, clv::ClvGinelli, z0, G0::SMatrix{d, k, T}, C0, nf::Int, ns::Int, nb::Int, orbit::Int) where {T, d, k}
    N = ns + nb; frames = o.wantV || o.wantG
    P = Matrix{T}(undef, k * (k + 1) ÷ 2, N)
    Gf = o.wantG ? o.G : Array{T, 3}(undef, d, k, frames ? ns + 1 : 0)
    z, G = _transient(clv, z0, G0, nf, orbit)
    xw = z; acc = ntuple(_ -> CompensatedSum{T}(), Val(k))
    for n in 0:N
        frames && n <= ns && _putframe!(Gf, n + 1, G)
        n == N && break
        z, G, R = _clv_handoff(clv, z, G, nf + n + 1, orbit)
        _pack!(P, n + 1, R)
        n < ns && (acc = _addlog(acc, R))
    end
    _start_backward!(o, C0, nb, ns, G, Val(k))
    _bwd_range!(o, P, Gf, 0, 0, N, ns, nf, Val(k), Val(d), frames, orbit)
    return xw, z, acc
end

# The backward matrix starts as C0; for n_b = 0 the sample n_s is C0 itself with NaN exponents (test plan U2).
function _start_backward!(o::ClvOut{T}, C0, nb::Int, ns::Int, Gfinal::SMatrix{d, k, T}, ::Val{k}) where {T, d, k}
    o.C .= C0
    if nb == 0
        fill!(o.lam, T(NaN))
        _deliver!(o, ns, Gfinal, nothing, 0)
    end
    return nothing
end

# Checkpointed (V17): a first pass through blocks 0 .. n_c - 2 stores (z, G) at each block start; then, from the last block, each block is rerun
# forward storing R and the frames and swept backward, the carried C passing from block to block. Last block first, so the window sum is in forward order.
function _run_ckpt!(o::ClvOut{T}, clv::ClvGinelli, z0, G0::SMatrix{d, k, T}, C0, nf::Int, ns::Int, nb::Int, orbit::Int) where {T, d, k}
    N = ns + nb; hb = clv.block; nc = cld(N, hb); frames = o.wantV || o.wantG; blk = min(hb, N)
    Zs = Vector{typeof(z0)}(undef, nc); Gs = Vector{typeof(G0)}(undef, nc)
    P = Matrix{T}(undef, k * (k + 1) ÷ 2, blk); Gf = Array{T, 3}(undef, d, k, frames ? blk : 0)
    z, G = _transient(clv, z0, G0, nf, orbit)
    xw = z; acc = ntuple(_ -> CompensatedSum{T}(), Val(k))
    for b in 0:(nc - 2), n in (b * hb):((b + 1) * hb - 1)
        n == b * hb && (Zs[b + 1] = z; Gs[b + 1] = G)
        z, G, R = _clv_handoff(clv, z, G, nf + n + 1, orbit)
        n < ns && (acc = _addlog(acc, R))
    end
    xend = z
    for b in (nc - 1):-1:0
        b < nc - 1 && ((z, G) = (Zs[b + 1], Gs[b + 1]))
        a = b * hb; e = min(N, a + hb)
        for n in a:(e - 1)
            frames && _putframe!(Gf, n - a + 1, G)
            z, G, R = _clv_handoff(clv, z, G, nf + n + 1, orbit)
            _pack!(P, n - a + 1, R)
            b == nc - 1 && n < ns && (acc = _addlog(acc, R))
        end
        b == nc - 1 && (xend = z; _start_backward!(o, C0, nb, ns, G, Val(k)))
        _bwd_range!(o, P, Gf, a, a, e, ns, nf, Val(k), Val(d), frames, orbit)
    end
    return xw, xend, acc
end

# ---------------------------------------------------------------------------------------------------------------- API
function _check_clv(clv::ClvGinelli, d::Int, k::Int, nf::Int, ns::Int, nb::Int, outputs, angle_sets, W0, C0)
    1 <= k <= d || throw(ArgumentError("k must be in 1:$d, got $k"))
    ns >= 1 || throw(ArgumentError("n_s must be at least 1, got $ns"))
    nf >= 0 || throw(ArgumentError("n_f must be nonnegative, got $nf"))
    nb >= 0 || throw(ArgumentError("n_b must be nonnegative, got $nb"))
    all(in((:V, :G, :R)), outputs) || throw(ArgumentError("outputs must be a tuple of :V, :G, :R, got $outputs"))
    angle_sets === nothing || all(r -> !isempty(r) && all(in(1:k), r), angle_sets) || throw(ArgumentError("angle_sets must be two nonempty column ranges within 1:$k"))
    W0 === nothing || size(W0) == (d, k) || throw(ArgumentError("W0 must be d x k = $d x $k, got $(size(W0))"))
    C0 === nothing || size(C0) == (k, k) || throw(ArgumentError("C0 must be k x k = $k x $k, got $(size(C0))"))
    return nothing
end

_clv_integrator(map) = hasproperty(map, :method) ? (method_name(map.method), nominal_order(map.method)) : ("", 0)

function _clv_record(clv::ClvGinelli, ::Type{T}, k, nf, ns, nb, seed, cseed, outputs, choice) where {T}
    integ, order = _clv_integrator(clv.map)
    return (system = sprint(show, clv.map), h = clv.h, tau = clv.tau, k = k, n_f = nf, n_s = ns, n_b = nb, orthonormalizer = orthonormalizer_name(clv.alg),
            storage = clv.storage, block = clv.block, precision = precision_label(T), seed = seed, cseed = cseed, outputs = outputs,
            backend = choice.backend, backend_fallback_reason = choice.reason, integrator = integ, integrator_order = order)
end

function _clv_out(clv::ClvGinelli, ::Type{T}, d::Int, k::Int, ns::Int, outputs, angle_sets, sink) where {T}
    nan(dims...) = fill(T(NaN), dims...)
    sets = angle_sets === nothing ? nothing : (collect(Int, angle_sets[1]), collect(Int, angle_sets[2]))
    m = sets === nothing ? 0 : min(length(sets[1]), length(sets[2]))
    big = sink === nothing
    return ClvOut(zeros(T, k, k), zeros(T, k, k), zeros(T, k), T(clv.tau) * T(clv.h), sink, sets,
                  big ? zeros(T, k, k, ns + 1) : zeros(T, 0, 0, 0), big ? nan(k, ns + 1) : zeros(T, 0, 0), big ? zeros(T, k, k, ns + 1) : zeros(T, 0, 0, 0),
                  big && :V in outputs ? zeros(T, d, k, ns + 1) : zeros(T, 0, 0, 0), big && :G in outputs ? zeros(T, d, k, ns + 1) : zeros(T, 0, 0, 0),
                  big && :R in outputs ? nan(k, k, ns + 1) : zeros(T, 0, 0, 0), sets === nothing || !big ? zeros(T, 0, 0) : nan(m, ns + 1),
                  big && :V in outputs, big && :G in outputs, big && :R in outputs)
end

"""
    propagate_clv_orbit(clv, x0, n_f, n_s, n_b; k, seed, orbit = 1, W0 = nothing, C0 = nothing, cseed = seed, outputs = (), angle_sets = nothing,
                        sink = nothing, backend = :cpu) -> ClvResult

One orbit: `n_f` forward handoffs (transient), the recording of `n_s + n_b` handoffs, the backward sweep from `C0` delivering the samples n = n_s .. 0
(`n_b` backward handoffs are the backward transient). `k` in `1:d`; `W0` the starting d x k frame (default `initial_deviations(seed, orbit, ...)`); `C0` the
starting k x k upper-triangular matrix (default [`clv_initial_coefficients`](@ref)`(cseed, orbit, ...)`). `outputs` is a tuple of `:V`, `:G`, `:R`;
`angle_sets = (r1, r2)` two column ranges whose principal angles are returned. `sink(n, C_n, lambda_n)` is called once per sample for n = n_s down to 0
with the matrix and vector of the run (copy to keep). `ArgumentError` before any work for `k` outside `1:d`, `n_s < 1`, `n_f < 0`, `n_b < 0`;
`NumericalHealthError` kinds "state", "deviation matrix", "orthonormalisation breakdown" (forward) and "backward breakdown"; `step` counts handoffs from 1
over the transient and the window in forward order. CPU only: `:gpu` and `:auto` run on the CPU and record the reason of `FALLBACK_REASONS[6]`.
"""
function propagate_clv_orbit(clv::ClvGinelli, x0::SVector{s, T}, nf::Integer, ns::Integer, nb::Integer; k::Integer, seed::Integer, orbit::Integer = 1,
                             W0 = nothing, C0 = nothing, cseed::Integer = seed, outputs::Tuple = (), angle_sets = nothing, sink = nothing,
                             backend::Symbol = :cpu) where {s, T}
    d = dimension(clv.map)
    _check_clv(clv, d, Int(k), Int(nf), Int(ns), Int(nb), outputs, angle_sets, W0, C0)
    choice = resolve_backend(backend, T; cpu_only_reducer = true)
    G0 = W0 === nothing ? initial_deviations(seed, orbit, Val(d), Val(Int(k)), T) : SMatrix{d, Int(k), T}(W0)
    Cm = C0 === nothing ? clv_initial_coefficients(cseed, orbit, Val(Int(k)), T) : SMatrix{Int(k), Int(k), T}(C0)
    all(isfinite, x0) || throw(NumericalHealthError("state", NaN, 0.0, 0, orbit, :cpu))
    all(isfinite, G0) || throw(NumericalHealthError("deviation matrix", NaN, 0.0, 0, orbit, :cpu))
    return _clv_orbit(clv, x0, G0, Cm, Int(nf), Int(ns), Int(nb), Int(orbit), outputs, angle_sets, sink,
                      _clv_record(clv, T, Int(k), Int(nf), Int(ns), Int(nb), seed, cseed, outputs, choice))
end

function _clv_orbit(clv::ClvGinelli, x0::SVector, G0::SMatrix{d, k, T}, Cm, nf, ns, nb, orbit, outputs, angle_sets, sink, record) where {d, k, T}
    o = _clv_out(clv, T, d, k, ns, outputs, angle_sets, sink)
    xw, xend, acc = clv.storage === :store ? _run_store!(o, clv, x0, G0, Cm, nf, ns, nb, orbit) : _run_ckpt!(o, clv, x0, G0, Cm, nf, ns, nb, orbit)
    lf = T[total(acc[i]) / (ns * clv.tau) / T(clv.h) for i in 1:k]
    frames = sink === nothing && (:V in outputs || :G in outputs)
    bytes = clv_bytes(Val(d), Val(k), T; n_s = ns, n_b = nb, storage = clv.storage, block = clv.block, frames = frames)
    sink === nothing && return ClvResult(o.Cs, o.lams, o.cosines, o.wantV ? o.V : nothing, o.wantG ? o.G : nothing, o.wantR ? o.R : nothing,
                                         o.sets === nothing ? nothing : o.angles, lf, xw, xend, bytes, record)
    return ClvResult{T, typeof(xw)}(nothing, nothing, nothing, nothing, nothing, nothing, nothing, lf, xw, xend, bytes, record)
end

"""
    ClvEnsemble

Result of [`propagate_clv`](@ref): `results` (a [`ClvResult`](@ref) per start, member j the run `propagate_clv_orbit(...; orbit = j)`, bitwise) and
`provenance` (`backend`, `device`, `backend_fallback_reason`, `precision`, `seed`, `k`, `chunk`, `ntasks`, `bytes_per_orbit`).
"""
struct ClvEnsemble{R <: ClvResult}
    results::Vector{R}
    provenance::NamedTuple
end

# A task failure reaches the caller as the original exception (the first orbit that fails, as the tangent propagate).
_unwrap(e) = e
_unwrap(e::CompositeException) = _unwrap(first(e.exceptions))
_unwrap(e::TaskFailedException) = _unwrap(e.task.result)

"""
    propagate_clv(clv, xs, n_f, n_s, n_b; k, seed, chunk = nothing, ntasks = nothing, backend = :cpu, outputs = (), angle_sets = nothing) -> ClvEnsemble

Many orbits, orbit id j the index in `xs`; the result of orbit j does not depend on `chunk`, `ntasks` or the position (spec I6). `chunk` defaults to the
orbits that fit in a quarter of the free host memory by [`clv_bytes`](@ref), `ntasks` to `Threads.nthreads()`. CPU only: a `:gpu` or `:auto` request runs on
the CPU and records "reducer is CPU-only" (`FALLBACK_REASONS[6]`). The first orbit that fails throws its `NumericalHealthError`.
"""
function propagate_clv(clv::ClvGinelli, xs::AbstractVector{<:SVector{s, T}}, nf::Integer, ns::Integer, nb::Integer; k::Integer, seed::Integer, chunk = nothing,
                       ntasks = nothing, backend::Symbol = :cpu, outputs::Tuple = (), angle_sets = nothing) where {s, T}
    isempty(xs) && throw(ArgumentError("xs must hold at least one start"))
    d = dimension(clv.map)
    choice = resolve_backend(backend, T; cpu_only_reducer = true)
    per = clv_bytes(Val(d), Val(Int(k)), T; n_s = ns, n_b = nb, storage = clv.storage, block = clv.block)
    ch = chunk === nothing ? clamp(Int(Sys.free_memory() ÷ (4 * per)), 1, length(xs)) : Int(chunk)
    nt = ntasks === nothing ? Threads.nthreads() : Int(ntasks)
    (ch >= 1 && nt >= 1) || throw(ArgumentError("chunk and ntasks must be at least 1"))
    one_orbit(j) = propagate_clv_orbit(clv, xs[j], nf, ns, nb; k = k, seed = seed, orbit = j, outputs = outputs, angle_sets = angle_sets, backend = :cpu)
    first_res = one_orbit(1)
    results = Vector{typeof(first_res)}(undef, length(xs)); results[1] = first_res
    try
        for part in Iterators.partition(2:length(xs), ch)
            @sync for sub in Iterators.partition(part, cld(length(part), nt))
                Threads.@spawn for j in sub
                    results[j] = one_orbit(j)
                end
            end
        end
    catch e
        throw(_unwrap(e))
    end
    prov = (backend = choice.backend, device = choice.device, backend_fallback_reason = choice.reason, precision = precision_label(T), seed = seed, k = Int(k),
            chunk = ch, ntasks = nt, bytes_per_orbit = per)
    return ClvEnsemble(results, prov)
end
