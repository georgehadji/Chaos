# Study p2-indicators (Amendment D89), pure rules: budget steps, growth ratio, verdicts R1 and R2, chaos scores, detection time (E26 to E28, E31).
# Nothing here reads a file; every threshold comes from the parameter table `p` (invariant I14).

const P2_I7 = [:LE, :SALI, :FLI, :MEGNO, :E1, :E2, :WBA]
const P2_I8 = vcat(P2_I7, :RLI)

_fin(x) = x isa Real && isfinite(x)
_bit2(w) = (Int(w) & 2) != 0                # the bit of the second flag name (degenerate, undefined window, zero) has the value 2^1 (P34)

"""
    p2_budget_step(t, b) -> Int

The first recorded step of the sorted vector `t` that is not before the budget `b` (E26). An `ArgumentError` naming the budget if there is none.
"""
function p2_budget_step(t::AbstractVector{<:Integer}, b::Real)
    i = searchsortedfirst(t, b)
    i <= length(t) || throw(ArgumentError("budget $b is above the last recorded step $(isempty(t) ? "(none)" : t[end])"))
    return Int(t[i])
end

"The recorded step nearest to `x`; a tie goes to the larger step (E26)."
function p2_nearest_step(t::AbstractVector{<:Integer}, x::Real)
    i = searchsortedfirst(t, x)
    i > length(t) && return Int(t[end])
    i == 1 && return Int(t[1])
    return Int(x - t[i - 1] < t[i] - x ? t[i - 1] : t[i])
end

"Position of the window W = b / 2 in `windows`; an `ArgumentError` naming the window if b is odd or W is not an element (E26)."
function p2_window_row(windows::AbstractVector{<:Integer}, b::Integer)
    iseven(b) || throw(ArgumentError("budget $b has no window: b / 2 is not an integer"))
    i = findfirst(==(b ÷ 2), windows)
    i === nothing && throw(ArgumentError("window $(b ÷ 2) (budget $b / 2) is not in the windows of the parameter table"))
    return i
end

"Growth function g_X(n) of E26: g_LE = n lambda_1(n), g_SALI = -L (L = ln GALI_2), g_FLI = fli, g_MEGNO = Ybar20_1 (P16)."
p2_growth(X::Symbol, n::Real, value::Real) = X === :LE ? n * value : X === :SALI ? -value : value

"r = g(n) / g(n'); `nothing` unless both are finite and g(n') is positive."
p2_ratio(g::Real, gp::Real) = (isfinite(g) && isfinite(gp) && gp > 0) ? g / gp : nothing

# Table R2 (sensitivity): the four rules that change; every other indicator keeps its R1 rule.
function _verdict_r2(X, v, p)
    X === :LE && return (_fin(v.lambda) && _fin(v.n)) ? (v.n * v.lambda >= p.r2_le_factor * log(v.n) ? :chaotic : :regular) : :undecided
    X === :SALI && return _fin(v.sali) ? (v.sali <= p.r2_sali_threshold ? :chaotic : :regular) : :undecided
    X === :FLI && return (_fin(v.fli) && _fin(v.n)) ? (v.fli >= p.r2_fli_factor * log(v.n) ? :chaotic : :regular) : :undecided
    X === :MEGNO && return _fin(v.ybar) ? (v.ybar >= p.r2_megno_threshold ? :chaotic : :regular) : :undecided
    return nothing
end

"""
    p2_verdict(X, rule, v, p) -> :chaotic | :regular | :undecided

Verdict of indicator `X` by table R1 or R2 (E27) from the value record `v`: LE `(n, lambda, r)`, SALI `(r, sali)`, FLI `(n, r, fli)`, MEGNO `(r, ybar)`,
RLI `(rli,)`, E1 and E2 `(D, flags)`, WBA `(dig, flags)` (the four rows of the window). `:STR` has no verdict.
"""
function p2_verdict(X::Symbol, rule::Symbol, v, p)
    rule in (:R1, :R2) || throw(ArgumentError("unknown verdict rule $rule"))
    X in P2_I8 || throw(ArgumentError("indicator $X has no verdict"))
    if rule === :R2
        r2 = _verdict_r2(X, v, p)
        r2 === nothing || return r2
    end
    if X === :LE
        (v.r === nothing || !_fin(v.lambda)) && return :undecided
        return (v.lambda >= p.le_floor && v.r >= p.ratio_threshold) ? :chaotic : :regular
    elseif X in (:SALI, :FLI, :MEGNO)
        v.r === nothing && return :undecided
        return v.r >= p.ratio_threshold ? :chaotic : :regular
    elseif X === :RLI
        return !_fin(v.rli) ? :undecided : (v.rli >= p.rli_threshold ? :chaotic : :regular)
    elseif X === :WBA
        (any(_bit2, v.flags) || !all(_fin, v.dig)) && return :undecided
        return minimum(v.dig) < p.wba_digits_threshold ? :chaotic : :regular
    end
    (_bit2(v.flags) || !_fin(v.D)) && return :undecided
    return v.D >= (X === :E1 ? p.e1_d_threshold : p.e2_d_threshold) ? :chaotic : :regular
end

"""
    p2_score(X, v) -> Union{Nothing, Float64}

Chaos score of E31, oriented so that a larger score means more chaotic (`nothing` if unavailable; minus infinity is kept for RLI = 0). `:STR` takes `(D,)`.
"""
function p2_score(X::Symbol, v)
    if X in (:LE, :SALI, :FLI, :MEGNO)
        return (v.r === nothing || !isfinite(v.r)) ? nothing : Float64(v.r)       # a ratio that overflowed to Inf is chaotic but has no score
    elseif X === :RLI
        (v.rli isa Real && !isnan(v.rli) && v.rli >= 0 && v.rli < Inf) || return nothing
        return v.rli == 0 ? -Inf : log10(Float64(v.rli))
    elseif X in (:E1, :E2)
        return (_bit2(v.flags) || !_fin(v.D)) ? nothing : Float64(v.D)
    elseif X === :WBA
        return (any(_bit2, v.flags) || !all(_fin, v.dig)) ? nothing : -Float64(minimum(v.dig))
    elseif X === :STR
        return (_fin(v.D) && v.D > 0) ? -log10(Float64(v.D)) : nothing
    end
    throw(ArgumentError("unknown indicator $X"))
end

"""
    p2_detection(V, rho, p) -> (t_det, outcome, censored)

Stable verdict, detection time, outcome and censored detection time (E28) from the ten verdicts `V` and the reference class `rho` (`:C` or `:R`).
"""
function p2_detection(V::AbstractVector{Symbol}, rho::Symbol, p)
    nb = length(V)
    last = V[nb]
    tdet = nothing
    if last !== :undecided
        for k in 1:nb
            p.budgets[k] * p.persist_divisor <= p.N_c || break
            if all(j -> V[j] === last, k:nb)
                tdet = p.budgets[k]
                break
            end
        end
    end
    outcome = tdet === nothing ? :undetected : ((last === :chaotic && rho === :C) || (last === :regular && rho === :R)) ? :correct : :wrong
    return (t_det = tdet, outcome = outcome, censored = outcome === :correct ? tdet : p.censor_factor * p.N_c)
end

"Median of m values: the middle value, or the mean of the two middle values for even m (Notation of the amendment)."
p2_median(v::AbstractVector{<:Real}) = _median_sorted(sort(Float64.(v)))
