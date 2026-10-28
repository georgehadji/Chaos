# Weighted mean phase advance (spec naff-fma: E2, NF4 to NF6, Amendment D81 items 1, 10, 12 to 14): the pure functions and the O(1) reader
# reducer. The pure function and the reducer share one update (`_accumulate`), one weight (`_wb`) and one increment (`_phase_inc`), so the
# stream result is bitwise the function on the series sampled at the handoffs.

const NF_E2_MIN_WINDOW = 16

"""
    birkhoff_weight(t, p = 1)

Weighted Birkhoff weight `exp(-1 / (t^p (1 - t)^p))` for `0 < t < 1`, `zero(t)` elsewhere (NF4, Das2017 Eqs. (3), (16)). `p` is Das's exponent,
not the window order of [`window_weights`](@ref).
"""
birkhoff_weight(t::T, p::Integer = 1) where {T <: AbstractFloat} = (t <= 0 || t >= 1) ? zero(T) : exp(-1 / (t^p * (1 - t)^p))

# The weight of sample j of a window of N: evaluated at the mirrored argument min(j, N - j) / N, so the weights are bitwise symmetric
# under j <-> N - j and the weighted mean of j is N/2 exactly (D81 item 10).
_wb(j::Integer, N::Integer, p::Integer, ::Type{T}) where {T} = birkhoff_weight(T(min(j, N - j)) / T(N), p)

"""
    birkhoff_weights(T, N; p = 1) -> Vector{T}

Normalised weights `w_hat_j = w(j/N) / sum`, `j = 0..N-1` (`w_hat_0 = 0`), the sum Neumaier-compensated (NF4, Das2017 Eq. (4)).
"""
function birkhoff_weights(::Type{T}, N::Integer; p::Integer = 1) where {T <: AbstractFloat}
    w = T[_wb(j, N, p, T) for j in 0:(N - 1)]
    s = CompensatedSum{T}()
    for x in w
        s = add(s, x)
    end
    return w ./ total(s)
end

# `z` scaled by an exact power of two so that its larger component lies in [1, 2): the direction is unchanged, the products cannot overflow.
@inline function _unit_scale(z::Complex{T}) where {T}
    m = Float64(max(abs(real(z)), abs(imag(z))))
    e = (isfinite(m) && m > 0) ? clamp(Int((reinterpret(UInt64, m) >> 52) & 0x7ff) - 1023, -1000, 1000) : 0   # the binary exponent, no throwing `exponent`
    return Complex(ldexp(real(z), -e), ldexp(imag(z), -e))
end

# Principal argument of z' / z in (-pi, pi]; NaN (undefined increment) if either sample is exactly 0. The product is formed directly while it is
# well inside the floating-point range, else from the power-of-two normalised samples (samples near 1e150 or 1e-160 and beyond).
function _phase_inc(zprev::Complex{T}, z::Complex{T}) where {T}
    (iszero(zprev) || iszero(z)) && return T(NaN)
    r = z * conj(zprev)
    s = max(abs(real(r)), abs(imag(r)))
    if !(s > T(1e-280) && s < T(1e280))
        r = _unit_scale(z) * conj(_unit_scale(zprev))
    end
    return atan(imag(r), real(r))
end

"""
    phase_increments(z) -> Vector

`Delta_j`, the principal argument of `z_{j+1} / z_j` in `(-pi, pi]` (NF6), as `atan(imag(r), real(r))` with `r = z_{j+1} conj(z_j)` (formed from
power-of-two normalised samples when the product leaves the safe range, so any finite magnitude works); `NaN` where one of the two samples is exactly 0. The reducer's planar observable uses the same function.
"""
phase_increments(z::AbstractVector{Complex{T}}) where {T <: AbstractFloat} = T[_phase_inc(z[j], z[j + 1]) for j in 1:(length(z) - 1)]

# One accumulation of a window: A += w d, B += w, both Neumaier. A weight that is exactly 0 adds nothing, so an undefined (NaN) increment with
# no weight does not spoil the window (the first sample of every window has the weight 0).
@inline _accumulate(A::CompensatedSum{T}, B::CompensatedSum{T}, w::T, d::T) where {T} = iszero(w) ? (A, B) : (add(A, w * d), add(B, w))
@inline _frequency(A::CompensatedSum{T}, B::CompensatedSum{T}, tp::T) where {T} = total(A) / (tp * total(B))

"""
    weighted_phase_frequency(increments, N; weights = nothing, p = 1) -> (nu1, nu2)

NF6 for the two windows of `N` increments (radians) of any float type, BigFloat included: `nu = A / (2 pi B)`, `A = sum w_j Delta_j`, `B = sum w_j`
(Neumaier sums over the unnormalised weights). The default weights are `birkhoff_weight` (exponent `p`) at the mirrored index, window 2 restarting
at `t = 0`; `weights` is a vector of `N` weights used as given (equal, `sin^2`, ... for the recorded line R2). `pi` is taken at the precision of the
data. `ArgumentError` naming `window` (`N` below 16), `increments` (fewer than `2N`), `weights` (length not `N`) or `birkhoff_order` (`p` not an
integer from 1 to 4 with the default weights: the centre weight underflows from `p = 5`).
"""
function weighted_phase_frequency(increments::AbstractVector{T}, N::Integer; weights = nothing, p::Integer = 1) where {T <: AbstractFloat}
    N >= NF_E2_MIN_WINDOW || throw(ArgumentError("window length $N is below the minimum window $NF_E2_MIN_WINDOW"))
    length(increments) >= 2N || throw(ArgumentError("increments: $(length(increments)) given, 2 window = $(2N) needed"))
    weights === nothing || length(weights) == N || throw(ArgumentError("weights: $(length(weights)) given, window = $N needed"))
    weights === nothing && _wb_check_order(p)
    tp = 2 * Systems._pi(first(increments))
    nu = ntuple(Val(2)) do h
        A = B = CompensatedSum{T}()
        for j in 0:(N - 1)
            w = weights === nothing ? _wb(j, N, p, T) : T(weights[j + 1])
            A, B = _accumulate(A, B, w, increments[(h - 1) * N + j + 1])
        end
        _frequency(A, B, tp)
    end
    return nu
end

abstract type NaffObservable end

"""
    PlanarObservable(iq, ip)

Observable of [`WeightedPhaseAdvance`](@ref): the sample is `z = x[iq] + i x[ip]` and the increment its [`phase_increments`](@ref) value
(planar elliptic maps, Y5 rotation blocks, the standard map with the wrapped `(q, p)`). A sample exactly at the origin has no phase: the window
that holds it has an undefined frequency (`NaN`), no failure is recorded.
"""
struct PlanarObservable <: NaffObservable
    iq::Int
    ip::Int
end

"""
    LiftedAngleObservable(iq, im = 0)

Observable of [`WeightedPhaseAdvance`](@ref): the sample is the pair `(x[iq], x[im])` of a wrapped coordinate and its winding count (`im = 0`: a
lifted coordinate with no count) and the increment is `(q' - q) + 2 pi (m' - m)` with `2 pi` in the element type as in `lift` (D28); for the
standard map this is the lifted post-kick momentum, so the weighted rotation number of NF5 in radians.
"""
struct LiftedAngleObservable <: NaffObservable
    iq::Int
    im::Int
end
LiftedAngleObservable(iq::Integer) = LiftedAngleObservable(iq, 0)

@inline _sample(o::PlanarObservable, x, ::Type{T}) where {T} = (T(x[o.iq]), T(x[o.ip]))
@inline _sample(o::LiftedAngleObservable, x, ::Type{T}) where {T} = (T(x[o.iq]), o.im == 0 ? zero(T) : T(x[o.im]))
@inline _increment(::PlanarObservable, a, b) = _phase_inc(Complex(a[1], a[2]), Complex(b[1], b[2]))
@inline _increment(::LiftedAngleObservable, a::Tuple{T, T}, b::Tuple{T, T}) where {T} = (b[1] - a[1]) + 2 * Systems._pi(one(T)) * (b[2] - a[2])
_observable_name(::PlanarObservable) = "planar"
_observable_name(::LiftedAngleObservable) = "lifted-angle"

"""
    WeightedPhaseAdvance(; window, observable, birkhoff_order = 1)

Weighted mean phase advance reducer (E2, Bazzani2023 pp.3, 19; spec P1 option (c)) for `propagate` with any `k >= 1`: a reader of the orbit point
that ignores the deviation matrix and returns it unchanged (the stream still carries one column, spec P1 fact (iii)). `window` is `N`, or a tuple of
`N` for several horizons of one run (one set of sums each); each `N` is at least 16 (`ArgumentError` naming `window`), `birkhoff_order` an integer
from 1 to 4 (`ArgumentError` naming `birkhoff_order`: the centre weight `exp(-4^p)` underflows from `p = 5`, so `B = 0` and every frequency is NaN). The first handoff is the sample 0 of the series, the increment `j` runs from sample `j` to `j + 1`, and a
run has `2 max(window) + 1` handoffs ([`check_run`](@ref)). Per horizon and window `nu = A / (2 pi B)` with the Neumaier sums `A = sum w_j Delta_j`,
`B = sum w_j` and the weight `birkhoff_weight` (NF6, D81 item 10); the state is isbits with no workspace: 96 bytes in Float64 and 176 in
Double64 for one horizon, plus 64 and 128 per extra horizon. A non-finite sample is recorded in the state (kind "series", the host raises).

`reducer_final` gives `(n, handoffs, complete, near_pi, nu1, nu2, delta, D, below_floor)`: `n` the last handoff (map step), `complete` true once the
largest horizon has its `2N + 1` handoffs, `near_pi` the number of increments that entered a window with `|Delta| > 0.9 pi` (planar observable only, 0 for
the lifted one; saturating at 32767): the principal argument is then near its branch cut and the estimate is biased, a property of the data to be read
before `nu` is used (spec I9); the rest `SVector{H}`: `nu1[h]` from the handoff `N_h + 1`, `nu2[h]`, `delta[h]` (circle distance), `D[h]`,
`below_floor[h]` from `2 N_h + 1`, `NaN` before and for an undefined window; after a horizon is complete further handoffs change nothing.
"""
struct WeightedPhaseAdvance{H, O <: NaffObservable}
    window::NTuple{H, Int}
    observable::O
    birkhoff_order::Int
end

function WeightedPhaseAdvance(; window, observable::NaffObservable, birkhoff_order::Integer = 1)
    ws = window isa Integer ? (Int(window),) : Tuple(Int.(window))
    (!isempty(ws) && all(>=(NF_E2_MIN_WINDOW), ws)) || throw(ArgumentError("window must be at least $NF_E2_MIN_WINDOW, got $window"))
    return WeightedPhaseAdvance{length(ws), typeof(observable)}(ws, observable, _wb_check_order(birkhoff_order))
end

struct WPAHorizon{T}
    A1::CompensatedSum{T}
    B1::CompensatedSum{T}
    A2::CompensatedSum{T}
    B2::CompensatedSum{T}
end
WPAHorizon{T}() where {T} = WPAHorizon{T}(CompensatedSum{T}(), CompensatedSum{T}(), CompensatedSum{T}(), CompensatedSum{T}())

# `prev1`, `prev2`: the previous sample; `n`: the last handoff (the failing one if `code` is 1), a map step as `Int`; `handoffs`: successful
# handoffs, saturating at the `Int32` maximum; `nwrap`: increments that entered a window with |Delta| > 0.9 pi (planar observable only, saturating).
# The field order keeps the size at 96 bytes (Float64) and 176 (Double64) for one horizon.
struct WPAState{T, H}
    hor::NTuple{H, WPAHorizon{T}}
    prev1::T
    prev2::T
    n::Int
    handoffs::Int32
    nwrap::Int16
    code::Int8
end

_observable_indices(o::PlanarObservable) = (o.iq, o.ip)
_observable_indices(o::LiftedAngleObservable) = o.im == 0 ? (o.iq,) : (o.iq, o.im)

# No argument check here: this runs in the device kernel, where a message built at run time cannot compile ([`check_run`](@ref) validates the observable).
function reducer_init(::WeightedPhaseAdvance{H}, ::Type{T}, ::Val, ::Val) where {H, T}
    return WPAState{T, H}(ntuple(_ -> WPAHorizon{T}(), Val(H)), zero(T), zero(T), 0, Int32(0), Int16(0), Int8(0))
end

# Increment j (from sample j to j + 1) into the window of the horizon `N` that holds it.
@inline function _advance(h::WPAHorizon{T}, N::Int, p::Int, j::Int, d::T) where {T}
    if j < N
        A, B = _accumulate(h.A1, h.B1, _wb(j, N, p, T), d)
        return WPAHorizon{T}(A, B, h.A2, h.B2)
    elseif j < 2N
        A, B = _accumulate(h.A2, h.B2, _wb(j - N, N, p, T), d)
        return WPAHorizon{T}(h.A1, h.B1, A, B)
    end
    return h
end

# Increments with |Delta| near pi are near the branch cut of the principal argument: the planar estimator is then biased (spec I9, a recorded flag).
@inline _near_pi(::PlanarObservable, d::T) where {T} = abs(d) > T(9) * Systems._pi(one(T)) / 10
@inline _near_pi(::LiftedAngleObservable, d) = false

function reducer_step(r::WeightedPhaseAdvance{H}, st::WPAState{T, H}, n::Int, x::SVector, W::SMatrix) where {H, T}
    st.code == 0 || return st, W
    a = _sample(r.observable, x, T)
    if !(isfinite(a[1]) && isfinite(a[2]))
        return WPAState{T, H}(st.hor, st.prev1, st.prev2, n, st.handoffs, st.nwrap, Int8(1)), W
    end
    s = Int(st.handoffs)
    h1 = s < typemax(Int32) ? Int32(s + 1) : st.handoffs
    s == 0 && return WPAState{T, H}(st.hor, a[1], a[2], n, h1, st.nwrap, Int8(0)), W
    d = _increment(r.observable, (st.prev1, st.prev2), a)
    hor = ntuple(h -> _advance(st.hor[h], r.window[h], r.birkhoff_order, s - 1, d), Val(H))
    flagged = _near_pi(r.observable, d) && s <= 2 * maximum(r.window) && st.nwrap < typemax(Int16)
    return WPAState{T, H}(hor, a[1], a[2], n, h1, flagged ? st.nwrap + Int16(1) : st.nwrap, Int8(0)), W
end

# Direct use throws at once (orbit unknown, 0); inside `propagate` the body uses `reducer_step`.
function reducer_update(r::WeightedPhaseAdvance, st::WPAState, n::Int, x::SVector, W::SMatrix)
    st1, W1 = reducer_step(r, st, n, x, W)
    reducer_failed(r, st1) && throw(reducer_error(r, st1, 0, :cpu))
    return st1, W1
end

reducer_failed(::WeightedPhaseAdvance, st::WPAState) = st.code != 0
reducer_error(::WeightedPhaseAdvance, st::WPAState, orbit, backend) = NumericalHealthError("series", NaN, Inf, st.n, orbit, backend, 1)

function reducer_final(r::WeightedPhaseAdvance{H}, st::WPAState{T, H}) where {H, T}
    tp = 2 * Systems._pi(one(T))
    s = Int(st.handoffs)
    nu1 = SVector{H, T}(ntuple(h -> s >= r.window[h] + 1 ? _frequency(st.hor[h].A1, st.hor[h].B1, tp) : T(NaN), Val(H)))
    nu2 = SVector{H, T}(ntuple(h -> s >= 2 * r.window[h] + 1 ? _frequency(st.hor[h].A2, st.hor[h].B2, tp) : T(NaN), Val(H)))
    delta = circle_distance.(nu1, nu2)
    idx = ntuple(h -> isnan(delta[h]) ? (D = T(NaN), below_floor = false) : frequency_index(delta[h]), Val(H))
    return (n = st.n, handoffs = s, complete = s >= 2 * maximum(r.window) + 1, near_pi = Int(st.nwrap), nu1 = nu1, nu2 = nu2, delta = delta,
            D = SVector{H, T}(ntuple(h -> idx[h].D, Val(H))), below_floor = SVector{H, Bool}(ntuple(h -> idx[h].below_floor, Val(H))))
end

"""
    check_run(r::WeightedPhaseAdvance, tau, N; d = nothing) -> nothing

The caller's check before the stream (the protocol gives `reducer_init` neither `tau` nor `N`, and the kernel cannot build error messages): `ArgumentError`
naming `window` unless the run of `N` map steps has `N / tau = 2 max(window) + 1` handoffs and `tau` divides `N` (D81 item 1). With `d` (the dimension
of the map) it also checks the observable indices against the orbit point, which has `d` entries (`2d` for the wrapped maps): `ArgumentError` naming `observable`.
"""
function check_run(r::WeightedPhaseAdvance, tau::Integer, N::Integer; d::Union{Nothing, Integer} = nothing)
    need = 2 * maximum(r.window) + 1
    (N % tau == 0 && N ÷ tau == need) || throw(ArgumentError("the run of $N steps at tau = $tau does not have the $need handoffs (2 window + 1) that window = $(r.window) needs"))
    idx = _observable_indices(r.observable)
    d === nothing || all(i -> 1 <= i <= 2d, idx) || throw(ArgumentError("observable indices $idx must lie in 1:$(2d) for a map of dimension $d"))
    return nothing
end

function reducer_provenance(r::WeightedPhaseAdvance)
    return (indicator = "naff-fma", estimator = "weighted-phase-advance", weight = "wb_exp", birkhoff_order = r.birkhoff_order, window = r.window,
            observable = _observable_name(r.observable), sample_position = "j/N", handoffs = 2 * maximum(r.window) + 1, floor_multiplier = 32,
            min_window = NF_E2_MIN_WINDOW)
end
