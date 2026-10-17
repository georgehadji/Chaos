# Class predicates on screening records (E11). Comparison directions as preregistered: "<=" for T, ">=" for C and S and the
# distances, strict "<" for the block test and the growth ratio.

"Median of the block exponents; the mean of the two middle values when their count is even."
function block_median(blocks)
    s = sort(collect(blocks))
    n = length(s)
    return isodd(n) ? s[(n + 1) ÷ 2] : (s[n ÷ 2] + s[n ÷ 2 + 1]) / 2
end

"Length in steps of the longest run of consecutive blocks with `e_i < f_S e_med` (number of blocks times `w`)."
function stuck_length(blocks, f_S, w)
    thr = f_S * block_median(blocks)
    best = run = 0
    for e in blocks
        run = e < thr ? run + 1 : 0
        best = max(best, run)
    end
    return best * w
end

"Predicate T: `lam <= theta_T`, `S(N/10) > 0`, growth ratio `S(N) / S(N/10) < g_T` and distance `d >= r_T` to the periodic set."
pred_T(lam, S10, SN, d; theta_T, g_T, r_T) = lam <= theta_T && S10 > 0 && SN / S10 < g_T && d >= r_T

"Predicate C: `lam >= theta_C` and distance `d >= r_C` to the elliptic subset."
pred_C(lam, d; theta_C, r_C) = lam >= theta_C && d >= r_C

"Predicate S: `lam >= theta_S`, a positive median block exponent and a stuck stretch of at least `L_S` steps."
function pred_S(lam, blocks; theta_S, w, f_S, L_S)
    lam >= theta_S || return false
    block_median(blocks) > 0 || return false
    return stuck_length(blocks, f_S, w) >= L_S
end
