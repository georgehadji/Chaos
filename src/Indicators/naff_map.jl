# Frequency map of a pair of windows and its index (spec naff-fma: NF1, NF9, NF10, P3, P4, P8, P11): the series path, pure functions.

"""
    circle_distance(a, b)

`abs(d - round(d))` with `d = a - b`, in the type of the arguments (NF9): the distance of two frequencies on the circle of circumference 1.
"""
circle_distance(a, b) = (d = a - b; abs(d - round(d)))

"""
    frequency_distance(nu1, nu2)

Euclidean norm of the plane circle distances (NF9) of two frequency vectors of equal length; `DimensionMismatch` otherwise.
"""
function frequency_distance(nu1::AbstractVector, nu2::AbstractVector)
    length(nu1) == length(nu2) || throw(DimensionMismatch("frequency vectors of lengths $(length(nu1)) and $(length(nu2))"))
    return sqrt(sum(abs2, circle_distance.(nu1, nu2)))
end

"""
    frequency_index(delta; floor_multiplier = 32) -> (D, below_floor, delta_floor)

NF10: `delta_floor = floor_multiplier * effective_eps(T)`, `D = log10(max(delta, delta_floor))`, `below_floor = delta <= delta_floor`.
`delta` is never clipped by the function (it only reports); the type of `delta` is kept.
"""
function frequency_index(delta::T; floor_multiplier::Real = 32) where {T <: AbstractFloat}
    fl = T(floor_multiplier) * effective_eps(T)
    return (D = log10(max(delta, fl)), below_floor = delta <= fl, delta_floor = fl)
end

"""
    FrequencyMap{T}

Result of [`frequency_map`](@ref) for one orbit: `w1`, `w2` (one [`WindowFrequency`](@ref) per plane), `dist` (per plane, the circle
distance after the mirror rule, `NaN` where undefined), `delta` (NF9, never clipped, `NaN` unless every plane is `:ok`), `D` (NF10),
`below_floor`, `status` (per plane: `:ok`, `:degenerate`, `:failed`), `failure` (`nothing`, or `(plane, window, sample)` of the first
non-finite sample, `sample` 1-based in the series) and `settings` (the identifiers of P11).
"""
struct FrequencyMap{T}
    w1::Vector{WindowFrequency{T}}
    w2::Vector{WindowFrequency{T}}
    dist::Vector{T}
    delta::T
    D::T
    below_floor::Bool
    status::Vector{Symbol}
    failure::Union{Nothing, NamedTuple{(:plane, :window, :sample), Tuple{Int, Int, Int}}}
    settings::NamedTuple
end

# Status of a plane from its two windows, and its distance with the mirror rule of P4 (the smaller of the direct and the sign-flipped
# distance when either window carries the flag).
function _plane_result(a::WindowFrequency{T}, b::WindowFrequency{T}) where {T}
    (a.status === :failed || b.status === :failed) && return :failed, T(NaN)
    (a.status === :degenerate || b.status === :degenerate) && return :degenerate, T(NaN)
    d = circle_distance(a.nu, b.nu)
    return :ok, (a.mirror || b.mirror) ? min(d, circle_distance(a.nu, -b.nu)) : d
end

function _check_map_args(N::Integer, p::Integer, pad::Integer, mirror_ratio::Real, floor_multiplier::Real)
    N >= 8 * (p + 1) || throw(ArgumentError("window length $N is below the minimum window $(8 * (p + 1)) for order $p"))
    pad >= 1 || throw(ArgumentError("pad must be positive, got $pad"))
    0 < mirror_ratio <= 1 || throw(ArgumentError("mirror_ratio must be in (0, 1], got $mirror_ratio"))
    floor_multiplier > 0 || throw(ArgumentError("floor_multiplier must be positive, got $floor_multiplier"))
    return nothing
end

"""
    frequency_map(series, N; order = 1, positions = :half, pad = 8, mean_removal = false, mirror_ratio = 0.9, floor_multiplier = 32) -> FrequencyMap

Frequency map distance of the pair of windows `[0, N)`, `[N, 2N)` of one plane (`series::AbstractVector{Complex{T}}`) or several (a vector
of such series, one per plane): the [`window_frequency`](@ref) of each window, the plane distances (NF9 with the mirror rule of P4),
`delta = sqrt(sum(dist.^2))` (NF1, unclipped) and `D = log10 max(delta, floor_multiplier eps)` with the `below_floor` flag (NF10). A series
longer than `2N` uses its first `2N` samples. A plane with a failed or degenerate window gives no `delta` for the pair and the plane status
names which; a non-finite sample is recorded in `failure` (the host raises, [`check_frequency_map`](@ref)). `ArgumentError` naming `window`,
`series` (fewer than `2N` samples), `pad`, `order`, `positions`, `mirror_ratio` or `floor_multiplier`, before any work.
"""
function frequency_map(series::AbstractVector{<:AbstractVector{Complex{T}}}, N::Integer; order::Integer = 1, positions::Symbol = :half,
                       pad::Integer = 8, mean_removal::Bool = false, mirror_ratio::Real = 0.9, floor_multiplier::Real = 32) where {T <: AbstractFloat}
    p = _check_order_positions(order, positions)
    _check_map_args(N, p, pad, mirror_ratio, floor_multiplier)
    all(z -> length(z) >= 2N, series) || throw(ArgumentError("series: every plane needs at least 2 window = $(2N) samples"))
    kw = (order = p, positions = positions, pad = pad, mean_removal = mean_removal, mirror_ratio = mirror_ratio)
    w1 = [window_frequency(Vector(view(z, 1:N)); kw...) for z in series]
    w2 = [window_frequency(Vector(view(z, (N + 1):(2N))); kw...) for z in series]
    res = [_plane_result(a, b) for (a, b) in zip(w1, w2)]
    status, dist = first.(res), T[last(r) for r in res]
    ok = all(==(:ok), status)
    delta = ok ? sqrt(sum(abs2, dist)) : T(NaN)
    fi = ok ? frequency_index(delta; floor_multiplier = floor_multiplier) : (D = T(NaN), below_floor = false)
    settings = (estimator = "windowed-fourier", weight = "cosine", window_order = p, positions = positions, pad = Int(pad),
                mean_removal = mean_removal, mirror_ratio = Float64(mirror_ratio), floor_multiplier = floor_multiplier, window = Int(N),
                eltype = precision_label(T), planes = length(series), min_window = 8 * (p + 1))
    return FrequencyMap{T}(w1, w2, dist, delta, fi.D, fi.below_floor, status, _first_failure(series, N), settings)
end
frequency_map(z::AbstractVector{Complex{T}}, N::Integer; kw...) where {T <: AbstractFloat} = frequency_map([z], N; kw...)

# `(plane, window, sample)` of the first non-finite sample of the first plane that has one.
function _first_failure(series, N::Integer)
    for (k, z) in enumerate(series)
        i = findfirst(!isfinite, view(z, 1:(2N)))
        i === nothing || return (plane = k, window = i <= N ? 1 : 2, sample = i)
    end
    return nothing
end

"""
    check_frequency_map(r; orbit = 0) -> r

The host check of a [`FrequencyMap`](@ref): returns `r`, or throws `NumericalHealthError("series", NaN, Inf, step = sample, orbit, :cpu, column = plane)`
for a recorded non-finite sample.
"""
function check_frequency_map(r::FrequencyMap; orbit::Integer = 0)
    f = r.failure
    f === nothing || throw(NumericalHealthError("series", NaN, Inf, f.sample, orbit, :cpu, f.plane))
    return r
end
