# Study p2-indicators (Amendment D89): agreement of the indicators (E30), AUC, regular range and leaving budget (E31).

_decided(v) = v === :chaotic || v === :regular

# Kahan sum of a vector of reals (compensated summation for accumulated sums).
function _ksum(v)
    s = 0.0
    c = 0.0
    for x in v
        y = x - c
        t = s + y
        c = (t - s) - y
        s = t
    end
    return s
end

"""
    p2_contingency(vx, vy) -> (cc, cr, rc, rr, n, A, kappa)

2 x 2 counts of the orbits decided for both indicators (first index the verdict of X), agreement `A = (cc + rr) / n` and Cohen's kappa (E30);
`A` and `kappa` are `nothing` if `n = 0`, and kappa is 0 if the expected agreement is 1.
"""
function p2_contingency(vx::AbstractVector, vy::AbstractVector)
    cc = cr = rc = rr = 0
    for (a, b) in zip(vx, vy)
        (_decided(a) && _decided(b)) || continue
        cc += (a === :chaotic && b === :chaotic)
        cr += (a === :chaotic && b === :regular)
        rc += (a === :regular && b === :chaotic)
        rr += (a === :regular && b === :regular)
    end
    n = cc + cr + rc + rr
    n == 0 && return (cc = cc, cr = cr, rc = rc, rr = rr, n = n, A = nothing, kappa = nothing)
    A = (cc + rr) / n
    pe = ((cc + cr) * (cc + rc) + (rc + rr) * (cr + rr)) / n^2
    return (cc = cc, cr = cr, rc = rc, rr = rr, n = n, A = A, kappa = pe == 1 ? 0.0 : (A - pe) / (1 - pe))
end

"""
    p2_unanimity(V, I) -> (U, used, dropped, A_min, A, kappa, kappa_mean)

`V` maps an indicator to its verdicts over the orbits of a set (`:na` where the orbit has no item), `I` is the indicator list. `U` is the fraction of the
orbits decided for every indicator of `I` whose verdicts are not all equal; `A` and `kappa` map each pair of `I` (in its order) to the agreement and kappa.
"""
function p2_unanimity(V::AbstractDict, I)
    n = isempty(I) ? 0 : length(V[first(I)])
    used = split = 0
    for o in 1:n
        vs = [V[X][o] for X in I]
        all(_decided, vs) || continue
        used += 1
        split += !all(==(vs[1]), vs)
    end
    pairs = [(I[a], I[b]) for a in eachindex(I) for b in (a + 1):length(I)]
    tabs = Dict(pr => p2_contingency(V[pr[1]], V[pr[2]]) for pr in pairs)
    As = [tabs[pr].A for pr in pairs]
    ks = [tabs[pr].kappa for pr in pairs]
    ok = all(!isnothing, As)
    return (U = used == 0 ? nothing : split / used, used = used, dropped = n - used, A_min = ok && !isempty(As) ? minimum(As) : nothing,
            A = Dict(pr => tabs[pr].A for pr in pairs), kappa = Dict(pr => tabs[pr].kappa for pr in pairs),
            kappa_mean = all(!isnothing, ks) && !isempty(ks) ? _ksum(ks) / length(ks) : nothing)
end

# Midranks (ties get the mean of the ranks they occupy).
function _midranks(x::AbstractVector)
    ix = sortperm(x)
    r = zeros(length(x))
    i = 1
    while i <= length(x)
        j = i
        while j < length(x) && x[ix[j + 1]] == x[ix[i]]
            j += 1
        end
        r[ix[i:j]] .= (i + j) / 2
        i = j + 1
    end
    return r
end

"Pearson correlation of the midranks (E30); `nothing` for fewer than 3 values or a constant rank vector."
function p2_spearman(x::AbstractVector{<:Real}, y::AbstractVector{<:Real})
    n = length(x)
    n < 3 && return nothing
    rx, ry = _midranks(Float64.(x)), _midranks(Float64.(y))
    (all(==(rx[1]), rx) || all(==(ry[1]), ry)) && return nothing
    mx, my = _ksum(rx) / n, _ksum(ry) / n
    return _ksum((rx .- mx) .* (ry .- my)) / sqrt(_ksum((rx .- mx) .^ 2) * _ksum((ry .- my) .^ 2))
end

"""
    p2_auc(c, r; min_set = 0)

Probability that a chaotic score exceeds a regular one, ties one half (rank-sum form, E31); `nothing` if either class has fewer than `min_set` scores.
"""
function p2_auc(c::AbstractVector{<:Real}, r::AbstractVector{<:Real}; min_set::Integer = 0)
    (length(c) < min_set || length(r) < min_set || isempty(c) || isempty(r)) && return nothing
    rk = _midranks(vcat(Float64.(c), Float64.(r)))
    nc, nr = length(c), length(r)
    return (_ksum(rk[1:nc]) - nc * (nc + 1) / 2) / (nc * nr)
end

# Percentile of a sorted vector, type 7; an exact order statistic is returned as it is (an infinite neighbour is not interpolated).
function _p2_percentile(s::AbstractVector, q::Real)
    n = length(s)
    h = (n - 1) * q + 1
    lo = floor(Int, h)
    hi = min(lo + 1, n)
    return h == lo ? s[lo] : s[lo] + (h - lo) * (s[hi] - s[lo])
end

"Regular range of E31: `(lo, hi, median, n)` from the lower and upper percentile (in percent, type 7) of the scores; `nothing` below `min_set` scores."
function p2_regular_range(scores::AbstractVector{<:Real}, p)
    length(scores) < p.min_set && return nothing
    s = sort(Float64.(scores))
    return (lo = _p2_percentile(s, p.range_lower_percentile / 100), hi = _p2_percentile(s, p.range_upper_percentile / 100), median = _median_sorted(s), n = length(s))
end

"The budget b_k of the smallest k from which the score is strictly above the upper percentile of the same budget up to k = 9; `nothing` if there is none."
function p2_leaving_budget(scores::AbstractVector, upper::AbstractVector, p)
    above(j) = scores[j] !== nothing && upper[j] !== nothing && scores[j] > upper[j]
    for k in eachindex(scores)
        all(above, k:length(scores)) && return p.budgets[k]
    end
    return nothing
end

"The share of the scores below the median `m` of the regular class (M10)."
p2_share_below(scores::AbstractVector{<:Real}, m::Real) = count(<(m), scores) / length(scores)
