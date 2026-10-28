# Included by verify/experiments-indicators.jl (section FILES): the orbit of an item by plain iteration of the wrapped Chirikov-Taylor map in Float64, written here from
# the equations of the map (decision D28, equations W-E2 and W-E5: p~ = p + K sin x and x~ = x + p', each reduced into [0, 2 pi) with its integer count, m' = m + n' + b),
# with no code of the module in the iteration (O17: the independent series that feeds the pure functions). `exact` only parses the decimal text of K, x0 and p0.
# Returns the states after j * tau steps, j = 1 .. nsamp (sample j - 1 of E21) as (q, p, m, n).
function plain_series(it, nsamp)
    K = Chaoslab.exact(Float64, it.param)
    tp = 2 * (3.141592653589793 + 1.2246467991473532e-16)      # 2 pi as the system forms it: hi + lo of pi, in the element type
    function red(u)
        a = floor(u / tp)
        r = u - tp * a
        if r < 0
            r += tp; a -= 1
        end
        if r >= tp
            r -= tp; a += 1
        end
        return r, a
    end
    q, m = red(Chaoslab.exact(Float64, it.x0))
    p, n = red(Chaoslab.exact(Float64, it.p0))
    out = NamedTuple{(:q, :p, :m, :n), NTuple{4, Float64}}[]
    for _ in 1:nsamp
        for _ in 1:it.tau
            pt, a = red(p + K * sin(q))
            n2 = n + a
            qt, b = red(q + pt)
            q, p, m, n = qt, pt, m + n2 + b, n2
        end
        push!(out, (q = q, p = p, m = m, n = n))
    end
    return out
end
