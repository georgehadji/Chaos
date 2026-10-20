# Bootstrap (E12): percentile interval of a median or of a ratio of medians, streams keyed by seed and label.

"UInt64 seed of a stream: first 8 bytes, big-endian, of SHA-256 over `chaoslab-analysis-1|<seed>|<label>` (the rule of the panel streams)."
stream_seed(label::AbstractString, seed::Integer) = reinterpret(UInt64, reverse(sha256(codeunits("chaoslab-analysis-1|$(seed)|$(label)"))[1:8]))[1]

bootstrap_rng(label::AbstractString, seed::Integer) = StableRNG(stream_seed(label, seed))

# Median of a sorted vector.
function _median_sorted(s::AbstractVector)
    n = length(s)
    return isodd(n) ? s[(n + 1) ÷ 2] : (s[n ÷ 2] + s[n ÷ 2 + 1]) / 2
end

median_of(v::AbstractVector{<:Real}) = _median_sorted(sort(Float64.(v)))

# Percentile of a sorted vector, linear interpolation between order statistics (type 7).
function _quantile_sorted(s::AbstractVector, p::Real)
    n = length(s)
    h = (n - 1) * p + 1
    lo = floor(Int, h)
    hi = min(lo + 1, n)
    return s[lo] + (h - lo) * (s[hi] - s[lo])
end

# `resamples` medians of resamples of `values` (m draws with replacement), in the stream of `rng`.
function _resample_medians(values::Vector{Float64}, rng, resamples::Integer)
    m = length(values)
    buf = Vector{Float64}(undef, m)
    out = Vector{Float64}(undef, resamples)
    for b in 1:resamples
        for i in 1:m
            buf[i] = values[rand(rng, 1:m)]
        end
        sort!(buf)
        out[b] = _median_sorted(buf)
    end
    return out
end

function _interval(stats::Vector{Float64}, level::Real)
    s = sort(stats)
    a = (1 - level) / 2
    return (lo = _quantile_sorted(s, a), hi = _quantile_sorted(s, 1 - a))
end

"""
    bootstrap_interval(values, label; seed, resamples, level) -> (point, lo, hi)

Median of `values` and the percentile interval of its bootstrap distribution (`resamples` resamples of the orbits, with replacement).
"""
function bootstrap_interval(values::AbstractVector{<:Real}, label::AbstractString; seed::Integer, resamples::Integer, level::Real)
    v = Float64.(values)
    iv = _interval(_resample_medians(v, bootstrap_rng(label, seed), resamples), level)
    return (point = median_of(v), lo = iv.lo, hi = iv.hi)
end

"""
    bootstrap_ratio(w, s, label_w, label_s; seed, resamples, level) -> (point, lo, hi)

Median of `w` over median of `s`; W and S are resampled independently, each in the stream of its own label, and the ratio is formed pairwise.
"""
function bootstrap_ratio(w::AbstractVector{<:Real}, s::AbstractVector{<:Real}, label_w::AbstractString, label_s::AbstractString;
                         seed::Integer, resamples::Integer, level::Real)
    wv, sv = Float64.(w), Float64.(s)
    q = _resample_medians(wv, bootstrap_rng(label_w, seed), resamples) ./ _resample_medians(sv, bootstrap_rng(label_s, seed), resamples)
    iv = _interval(q, level)
    return (point = median_of(wv) / median_of(sv), lo = iv.lo, hi = iv.hi)
end

"""
    bootstrap_band(M, label; seed, resamples, level) -> (median, lo, hi)

Pointwise band of the median of the rows of `M` (orbits x steps): one set of resampled orbit indices (stream of `label`) serves every column.
"""
function bootstrap_band(M::AbstractMatrix{<:Real}, label::AbstractString; seed::Integer, resamples::Integer, level::Real)
    m, ns = size(M)
    rng = bootstrap_rng(label, seed)
    ix = [rand(rng, 1:m) for _ in 1:m, _ in 1:resamples]
    med, lo, hi = zeros(ns), zeros(ns), zeros(ns)
    buf = Vector{Float64}(undef, m)
    stats = Vector{Float64}(undef, resamples)
    for j in 1:ns
        med[j] = median_of(M[:, j])
        for b in 1:resamples
            for i in 1:m
                buf[i] = M[ix[i, b], j]
            end
            sort!(buf)
            stats[b] = _median_sorted(buf)
        end
        iv = _interval(stats, level)
        lo[j], hi[j] = iv.lo, iv.hi
    end
    return (median = med, lo = lo, hi = hi)
end
