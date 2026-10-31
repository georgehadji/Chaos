# Study p3-abstention (Amendment D93): the standard results of E56 and E57, defined in full in the amendment (P61): the exact one-sided tail of the McNemar test and its critical
# values, rounding half away from zero, the multiple of the step, the Hanley and McNeil variance of an AUC, the Mann-Whitney AUC, the binormal shift, the votes of HC.

"The exact tail `P(Bin(m, 1/2) >= n01)` with `m = n10 + n01` as a rational (1 for `m = 0`)."
function p3_mcnemar_p(n10::Integer, n01::Integer)
    m = n10 + n01
    m == 0 && return Rational{BigInt}(1)
    return sum(binomial(big(m), j) for j in n01:m) // big(2)^m
end

"""
    p3_mcnemar_critical(m, alpha) -> Union{Nothing, Int}

The smallest `n01` whose tail `P(Bin(m, 1/2) >= n01)` is at most the rational `alpha` (the critical-value table of PP, exact integer arithmetic), or `nothing`.
"""
function p3_mcnemar_critical(m::Integer, alpha::Rational)
    lim = alpha * big(2)^m
    tail, c, best = big(0), big(1), nothing
    for j in m:-1:0
        tail += c
        tail <= lim || break
        best = j
        c = c * j ÷ (m - j + 1)
    end
    return best
end

"Round half away from zero of a rational (an exact product of keys, so the tie cases are exact)."
p3_round_half_away(x::Rational) = x >= 0 ? Int(floor(x + 1 // 2)) : -Int(floor(-x + 1 // 2))

"The rational `x` rounded up to a multiple of the integer `step`."
p3_ceil_step(x::Rational, step::Integer) = step * Int(ceil(x / step))

"The Hanley and McNeil variance of an AUC `A` with `na` failures and `nn` non-failures (E57)."
function p3_hanley_mcneil(A::Real, na::Integer, nn::Integer)
    Q1 = A / (2 - A)
    Q2 = 2 * A^2 / (1 + A)
    return (A * (1 - A) + (na - 1) * (Q1 - A^2) + (nn - 1) * (Q2 - A^2)) / (na * nn)
end

# The Mann-Whitney count of E57: twice the number of pairs (failure above non-failure) plus the tied pairs, by ranks (sort of the non-failures).
function _p3_auc_count(f, n)
    ns = sort(collect(n))
    s = 0
    for x in f
        lo = searchsortedfirst(ns, x) - 1
        s += 2 * lo + (searchsortedlast(ns, x) - lo)
    end
    return s
end

"The AUC of the scores of the failures `f` against the non-failures `n` (Mann-Whitney, ties one half) as an exact rational, or `nothing` if a class is empty."
function p3_auc(f, n)
    (isempty(f) || isempty(n)) && return nothing
    return _p3_auc_count(f, n) // (2 * length(f) * length(n))
end

"The binormal shift `d(A) = sqrt(2) Phi^-1(A)` (E57); an `ArgumentError` naming the AUC outside the open unit interval."
function p3_binormal_shift(A::Real)
    (0 < A < 1) || throw(ArgumentError("the AUC $A must lie in the open unit interval for the binormal model"))
    return sqrt(2.0) * p3_normal_quantile(A)
end

"""
    p3_vote_summary(votes, rho) -> (nC, nR, label, f, s)

The votes of one orbit in HC (E57): the plain label is chaotic if `n_C > n_R`, regular if `n_R > n_C`, else a tie (a failure); `f` iff the label differs from `rho`
(`:C` or `:R`); the score `s = 1 - |n_C - n_R| / (n_C + n_R)`, and 0 if there is no decided vote.
"""
function p3_vote_summary(votes, rho::Symbol)
    nC, nR = count(==(:chaotic), votes), count(==(:regular), votes)
    label = nC > nR ? :chaotic : nR > nC ? :regular : :tie
    f = label === :tie || (label === :chaotic) != (rho === :C)
    return (nC = nC, nR = nR, label = label, f = f, s = nC + nR == 0 ? 0.0 : 1 - abs(nC - nR) / (nC + nR))
end
