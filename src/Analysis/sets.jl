# Orbit sets (E7), finite-time spread and H2 (E10), convergence time and H3 (E11).

# S of the preregistration (H3): chaotic-sea orbits at K = 2 pi and K = 10. The parameter table has no key for these two values (spec gap, see
# the implementation notes), so they are named here as part of the definition of the set.
const STRONG_K = (2 * pi, 10.0)
_strong(K) = any(k -> abs(K - k) <= 1e-9, STRONG_K)

const POSITIVE = "positive growth supported"

"""
    orbit_sets(orbits, p) -> (E_ch, W, S, S_nonpositive, dropped)

`orbits` carry `id`, `K`, `class`, `control`, `exponent`, `flags` and the final `instability` label. Orbits flagged "excluded" or "drift" are
removed from all sets and listed with the reason.
"""
function orbit_sets(orbits, p)
    E, W, S, dropped = Int[], Int[], Int[], NamedTuple[]
    npos = 0
    for o in orbits
        reason = "excluded" in o.flags ? "excluded" : ("drift" in o.flags ? "drift" : "")
        if !isempty(reason)
            push!(dropped, (id = o.id, reason = reason))
            continue
        end
        o.control && o.exponent > 0 && push!(E, o.id)
        o.control && continue
        o.class == "C" && o.K > 0 && push!(E, o.id)
        pos = o.instability == POSITIVE
        (o.class in ("S", "X") && o.K > 0 && pos || o.class == "C" && 0 < o.K <= p.K_c && pos) && push!(W, o.id)
        if o.class == "C" && _strong(o.K)
            push!(S, o.id)
            npos += !pos
        end
    end
    return (E_ch = E, W = W, S = S, S_nonpositive = npos, dropped = dropped)
end

"""
    relative_spread(lam, lam_ref) -> Vector

R(n) = (max over V of lambda_m(n) - min over V) / lambda_ref at every recorded step; L2 is not in V.
"""
function relative_spread(lam::AbstractDict, lam_ref::Real)
    cols = [lam[m] for m in V_METHODS if haskey(lam, m)]
    return [(maximum(c[i] for c in cols) - minimum(c[i] for c in cols)) / lam_ref for i in eachindex(first(cols))]
end

"""
    convergence_time(steps, lam, lam_ref, band) -> (n, defined)

Smallest recorded step n with |lambda(k) - lambda_ref| <= band lambda_ref for every recorded k >= n. Outside the band at the last step, or
lambda_ref <= 0: undefined, n is the last step.
"""
function convergence_time(steps::AbstractVector{<:Integer}, lam::AbstractVector{<:Real}, lam_ref::Real, band::Real)
    N = Int(steps[end])
    lam_ref > 0 || return (n = N, defined = false)
    tol = band * lam_ref
    inside(i) = abs(lam[i] - lam_ref) <= tol
    inside(length(steps)) || return (n = N, defined = false)
    j = findlast(i -> !inside(i), eachindex(steps))
    return (n = Int(steps[j === nothing ? 1 : j + 1]), defined = true)
end

"Convergence time of an orbit: the maximum over V; undefined if any method is."
function orbit_convergence(steps, lam::AbstractDict, lam_ref::Real, band::Real)
    r = [convergence_time(steps, lam[m], lam_ref, band) for m in V_METHODS if haskey(lam, m)]
    return (n = maximum(x -> x.n, r), defined = all(x -> x.defined, r))
end

_undecided(n, p) = n < p.min_set

"""
    decide_h2(R, p; seed) -> (verdict, median, lo, hi, n)

R is R(n_H2) over W. Supported if the median reaches `h2_median` and the lower end of its bootstrap interval reaches `h2_interval`.
"""
function decide_h2(R::AbstractVector{<:Real}, p; seed::Integer)
    n = length(R)
    n == 0 && return (verdict = :not_decidable, median = NaN, lo = NaN, hi = NaN, n = 0)
    b = bootstrap_interval(R, "H2/W"; seed = seed, resamples = p.bootstrap_resamples, level = p.bootstrap_level)
    verdict = _undecided(n, p) ? :not_decidable : (b.point >= p.h2_median && b.lo >= p.h2_interval ? :supported : :not_supported)
    return (verdict = verdict, median = b.point, lo = b.lo, hi = b.hi, n = n)
end

"""
    decide_h3(nconv_W, nconv_S, p; seed) -> (verdict, ratio, lo, hi, n_w, n_s)

Supported if median(W) / median(S) reaches `h3_ratio` and the lower end of the bootstrap interval of the ratio is strictly above `h3_interval`.
"""
function decide_h3(w::AbstractVector{<:Real}, s::AbstractVector{<:Real}, p; seed::Integer)
    nw, ns = length(w), length(s)
    (nw == 0 || ns == 0) && return (verdict = :not_decidable, ratio = NaN, lo = NaN, hi = NaN, n_w = nw, n_s = ns)
    q = bootstrap_ratio(w, s, "H3/W", "H3/S"; seed = seed, resamples = p.bootstrap_resamples, level = p.bootstrap_level)
    verdict = (_undecided(nw, p) || _undecided(ns, p)) ? :not_decidable : (q.point >= p.h3_ratio && q.lo > p.h3_interval ? :supported : :not_supported)
    return (verdict = verdict, ratio = q.point, lo = q.lo, hi = q.hi, n_w = nw, n_s = ns)
end
