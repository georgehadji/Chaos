# The weighted Birkhoff reader reducer (spec weighted-birkhoff: P1, P6 to P8, WB14, Amendment D84): K observables of the orbit point over H horizons in
# one pass, O(1) isbits state. The weight, the sums and the index are those of `wbirkhoff.jl`; nothing here may throw a message built at run time
# (the device kernel would not compile): arguments are checked in the constructor and in `check_run`.

"""
    WeightedBirkhoff(; window, observables, birkhoff_order = 1)

Weighted Birkhoff average reducer (I13, Das2017, DasYorke2018) for `propagate` with any `k >= 1`: a reader of the orbit point that returns the deviation
matrix unchanged and replaces no column (any position in a chain; the stream still carries one deviation column). `window` is `N`, or a tuple of `N`
for several horizons of one pass; each `N` is at least 16 (`ArgumentError` naming `window`). `observables` is a nonempty tuple of `K` isbits callables
of the orbit point (the wrapped state of `WrappedStandardMap`, a lifted point, a plane of a block map), for example [`CosObservable`](@ref),
[`SinObservable`](@ref), [`LiftedObservable`](@ref) or any function; the value is converted to the element type of the orbit (`ArgumentError` naming
`observables`). The caller guarantees smoothness. The first handoff is the sample 0 of the series; a run has `2 max(window)` handoffs
([`check_run`](@ref)). Per observable and horizon: `avg = A / B` of the two windows with the weight `birkhoff_weight` of exponent `birkhoff_order`
(`>= 1`, `ArgumentError` naming `birkhoff_order`), Neumaier sums, and the largest `abs(f)` of the `2N` samples. The state is isbits with no workspace:
`16 + sizeof(T) H (2 + 5 K)` bytes. A non-finite observable is recorded in the state as its index (kind "observable"); the host raises with the orbit id.

`reducer_final` gives `(n, handoffs, complete, avg1, avg2, S, B, delta, absdig, reldig, dig, below_floor, zero)`: `n` the last handoff (map step),
`complete` true once the largest horizon has its `2N` handoffs; `avg1`, `avg2`, `S`, `delta`, `absdig`, `reldig`, `dig` are `SMatrix{K, H}` (row the
observable, column the horizon in the order given), `below_floor` and `zero` the `Bool` matrices of [`birkhoff_index`](@ref), `B` the `SVector{H}` of
weight totals. `avg1` is defined from the handoff `N_h`, `avg2`, `delta` and the digits from `2 N_h` (`NaN` before); `S` is the running maximum; a
completed horizon never changes.
"""
struct WeightedBirkhoff{H, K, O <: Tuple}
    window::NTuple{H, Int}
    observables::O
    birkhoff_order::Int
end

function WeightedBirkhoff(; window, observables, birkhoff_order::Integer = 1)
    ws = _wb_windows(window)
    p = _wb_check_order(birkhoff_order)
    observables isa Tuple && 1 <= length(observables) <= typemax(Int16) && all(o -> isbitstype(typeof(o)), observables) ||
        throw(ArgumentError("observables must be a nonempty tuple of at most $(typemax(Int16)) isbits callables"))
    return WeightedBirkhoff{length(ws), length(observables), typeof(observables)}(ws, observables, p)
end

# The sums of one horizon: the weight total B (both windows share it) and a cell per observable.
struct WBHorizon{T, K}
    B::CompensatedSum{T}
    cells::NTuple{K, WBCell{T}}
end

# `n`: the last handoff (the failing one if `code` is not 0), a map step as `Int`; `handoffs`: successful handoffs, saturating; `code`: the index of
# the first non-finite observable of the failing handoff. The three fields take 16 bytes, so the size is 16 + sizeof(T) H (2 + 5 K).
struct WBState{T, H, K}
    hor::NTuple{H, WBHorizon{T, K}}
    n::Int
    handoffs::Int32
    code::Int16
end

function reducer_init(::WeightedBirkhoff{H, K}, ::Type{T}, ::Val, ::Val) where {H, K, T}
    hor = WBHorizon{T, K}(CompensatedSum{T}(), ntuple(_ -> WBCell{T}(), Val(K)))
    return WBState{T, H, K}(ntuple(_ -> hor, Val(H)), 0, Int32(0), Int16(0))
end

# Sample j (0-based) into the window of the horizon `N` that holds it; a horizon is frozen after its 2N handoffs. B takes the weights of window 1 only.
@inline function _wb_advance(h::WBHorizon{T, K}, N::Int, p::Int, j::Int, vals::NTuple{K, T}) where {T, K}
    j < 2N || return h
    i = j >= N ? 1 : 0
    w = _wb(j - i * N, N, p, T)
    B = (i == 0 && !iszero(w)) ? add(h.B, w) : h.B
    return WBHorizon{T, K}(B, ntuple(k -> _wb_sample(h.cells[k], i, w, vals[k]), Val(K)))
end

# Index of the first non-finite value, 0 if none.
@inline function _first_nonfinite(vals::Tuple)
    for k in 1:length(vals)
        isfinite(vals[k]) || return k
    end
    return 0
end

function reducer_step(r::WeightedBirkhoff{H, K}, st::WBState{T, H, K}, n::Int, x::SVector, W::SMatrix) where {H, K, T}
    st.code == 0 || return st, W
    vals = map(o -> T(o(x)), r.observables)::NTuple{K, T}
    bad = _first_nonfinite(vals)
    bad == 0 || return WBState{T, H, K}(st.hor, n, st.handoffs, Int16(bad)), W
    j = Int(st.handoffs)
    hor = ntuple(h -> _wb_advance(st.hor[h], r.window[h], r.birkhoff_order, j, vals), Val(H))
    return WBState{T, H, K}(hor, n, j < typemax(Int32) ? Int32(j + 1) : st.handoffs, Int16(0)), W
end

# Direct use throws at once (orbit unknown, 0); inside `propagate` the body uses `reducer_step`.
function reducer_update(r::WeightedBirkhoff, st::WBState, n::Int, x::SVector, W::SMatrix)
    st1, W1 = reducer_step(r, st, n, x, W)
    reducer_failed(r, st1) && throw(reducer_error(r, st1, 0, :cpu))
    return st1, W1
end

reducer_failed(::WeightedBirkhoff, st::WBState) = st.code != 0
reducer_error(::WeightedBirkhoff, st::WBState, orbit, backend) = NumericalHealthError("observable", NaN, Inf, st.n, orbit, backend, Int(st.code))

# The index fields of one (observable, horizon): NaN digits and clear flags until the horizon has its 2N handoffs.
function _wb_cell_final(c::WBCell{T}, B::T, s::Int, N::Int) where {T}
    a1 = s >= N ? total(c.A1) / B : T(NaN)
    s >= 2N || return (avg1 = a1, avg2 = T(NaN), S = c.S, ix = (delta = T(NaN), absdig = T(NaN), reldig = T(NaN), dig = T(NaN), below_floor = false, zero = false))
    a2 = total(c.A2) / B
    return (avg1 = a1, avg2 = a2, S = c.S, ix = birkhoff_index(a1, a2, c.S))
end

function reducer_final(r::WeightedBirkhoff{H, K}, st::WBState{T, H, K}) where {H, K, T}
    s = Int(st.handoffs)
    cells = ntuple(Val(K * H)) do i
        h, k = (i - 1) ÷ K + 1, (i - 1) % K + 1
        _wb_cell_final(st.hor[h].cells[k], total(st.hor[h].B), s, r.window[h])
    end
    M(f) = SMatrix{K, H}(ntuple(i -> f(cells[i]), Val(K * H)))
    return (n = st.n, handoffs = s, complete = s >= 2 * maximum(r.window), avg1 = M(c -> c.avg1), avg2 = M(c -> c.avg2), S = M(c -> c.S),
            B = SVector{H, T}(ntuple(h -> total(st.hor[h].B), Val(H))), delta = M(c -> c.ix.delta), absdig = M(c -> c.ix.absdig), reldig = M(c -> c.ix.reldig),
            dig = M(c -> c.ix.dig), below_floor = M(c -> c.ix.below_floor), zero = M(c -> c.ix.zero))
end

"""
    check_run(r::WeightedBirkhoff, tau, N; x0 = nothing) -> nothing

The caller's check before the stream (the protocol gives `reducer_init` neither `tau` nor `N`, and the kernel cannot build error messages):
`ArgumentError` naming `window` unless `tau` divides `N` and the run has `N / tau = 2 max(window)` handoffs. With `x0`, the initial point, every observable
is evaluated once: `ArgumentError` naming `observable` for one that throws, is not a real number or is not finite there. A chain that holds
[`WeightedPhaseAdvance`](@ref) takes the longer run its members need; the frozen horizons make the extra handoff harmless, so run the check on the
stand-alone reducer only (spec D84 item 1).
"""
function check_run(r::WeightedBirkhoff, tau::Integer, N::Integer; x0 = nothing)
    need = 2 * maximum(r.window)
    (N % tau == 0 && N ÷ tau == need) || throw(ArgumentError("the run of $N steps at tau = $tau does not have the $need handoffs (2 window) that window = $(r.window) needs"))
    x0 === nothing && return nothing
    for (k, o) in enumerate(r.observables)
        v = try
            o(x0)
        catch e
            throw(ArgumentError("observable $k throws at the initial point: $(sprint(showerror, e))"))
        end
        (v isa Real && isfinite(v)) || throw(ArgumentError("observable $k is not a finite real number at the initial point: $v"))
    end
    return nothing
end

function reducer_provenance(r::WeightedBirkhoff)
    return (indicator = "weighted-birkhoff", estimator = "weighted-birkhoff-average", weight = "wb_exp", birkhoff_order = r.birkhoff_order, window = r.window,
            observables = map(observable_name, r.observables), handoffs = 2 * maximum(r.window), floor_multiplier = WB_FLOOR_MULTIPLIER,
            min_window = WB_MIN_WINDOW)
end
