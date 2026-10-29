# Amendment D89, stage 02. Independent reference implementation of the P2 decision layer (spec E26 to E34), written from the amendment text
# and used as the oracle of the tests (spec O6, in Julia instead of Python). It never calls Chaoslab.Analysis, so it works before the module exists.
# Names: p2ref_*. Verdict symbols :chaotic, :regular, :undecided; :na marks an orbit without an item of the indicator.

const c_, r_, u_ = :chaotic, :regular, :undecided          # verdict symbols of the tests
const P2_X7 = [:LE, :SALI, :FLI, :MEGNO, :E1, :E2, :WBA]
const P2_X8 = vcat(P2_X7, :RLI)

p2ref_median(v) = (s = sort(Float64.(v)); n = length(s); isodd(n) ? s[(n + 1) ÷ 2] : (s[n ÷ 2] + s[n ÷ 2 + 1]) / 2)

# type 7 percentile (q in [0, 1]); infinities are kept, so an equal pair of neighbours is not interpolated
function p2ref_percentile(v, q)
    s = sort(Float64.(v)); n = length(s)
    h = (n - 1) * q + 1
    lo = floor(Int, h); hi = min(lo + 1, n)
    return h == lo ? s[lo] : s[lo] + (h - lo) * (s[hi] - s[lo])
end

# budget rule (E26): first recorded step not before b; the recorded step nearest to x, a tie going to the larger step
p2ref_n(t, b) = (i = findfirst(>=(b), t); i === nothing ? nothing : t[i])
function p2ref_nearest(t, x)
    best = t[1]
    for s in t
        d, db = abs(s - x), abs(best - x)
        (d < db || (d == db && s > best)) && (best = s)
    end
    return best
end

# ---------------------------------------------------------------- verdicts, E27. v: see the test plan, "API proposed by the tests".
function p2ref_verdict(X::Symbol, rule::Symbol, v, p)
    fin(x) = x isa Real && isfinite(x)
    bit2(w) = (Int(w) & 2) != 0                      # the bit values are 2^j: degenerate, undefined window and zero are the second name (value 2)
    if rule === :R2
        X === :LE && return (fin(v.lambda) && fin(v.n)) ? (v.n * v.lambda >= p.r2_le_factor * log(v.n) ? :chaotic : :regular) : :undecided
        X === :SALI && return fin(v.sali) ? (v.sali <= p.r2_sali_threshold ? :chaotic : :regular) : :undecided
        X === :FLI && return (fin(v.fli) && fin(v.n)) ? (v.fli >= p.r2_fli_factor * log(v.n) ? :chaotic : :regular) : :undecided
        X === :MEGNO && return fin(v.ybar) ? (v.ybar >= p.r2_megno_threshold ? :chaotic : :regular) : :undecided
    end
    if X === :LE
        (v.r === nothing || !fin(v.lambda)) && return :undecided
        return (v.lambda >= p.le_floor && v.r >= p.ratio_threshold) ? :chaotic : :regular
    elseif X in (:SALI, :FLI, :MEGNO)
        v.r === nothing && return :undecided
        return v.r >= p.ratio_threshold ? :chaotic : :regular
    elseif X === :RLI
        fin(v.rli) || return :undecided
        return v.rli >= p.rli_threshold ? :chaotic : :regular
    elseif X === :E1
        (bit2(v.flags) || !fin(v.D)) && return :undecided
        return v.D >= p.e1_d_threshold ? :chaotic : :regular
    elseif X === :E2
        (bit2(v.flags) || !fin(v.D)) && return :undecided
        return v.D >= p.e2_d_threshold ? :chaotic : :regular
    elseif X === :WBA
        (any(bit2, v.flags) || !all(fin, v.dig)) && return :undecided
        return minimum(v.dig) < p.wba_digits_threshold ? :chaotic : :regular
    end
    error("p2ref_verdict: unknown indicator $X")
end

# growth ratio (E26): g_LE = n lambda, g_SALI = -L, g_FLI = fli, g_MEGNO = Ybar20_1; undefined unless both are finite and g(n') > 0
p2ref_g(X::Symbol, n, value) = X === :LE ? n * value : X === :SALI ? -value : value
p2ref_ratio(g, gp) = (isfinite(g) && isfinite(gp) && gp > 0) ? g / gp : nothing

function p2ref_score(X::Symbol, v)
    fin(x) = x isa Real && isfinite(x)
    if X in (:LE, :SALI, :FLI, :MEGNO)
        return v.r === nothing ? nothing : Float64(v.r)
    elseif X === :RLI
        (v.rli isa Real && !isnan(v.rli) && v.rli >= 0 && v.rli < Inf) || return nothing
        return v.rli == 0 ? -Inf : log10(v.rli)
    elseif X in (:E1, :E2)
        return ((Int(v.flags) & 2) != 0 || !fin(v.D)) ? nothing : Float64(v.D)
    elseif X === :WBA
        return (any(w -> (Int(w) & 2) != 0, v.flags) || !all(fin, v.dig)) ? nothing : -minimum(v.dig)
    elseif X === :STR
        return (fin(v.D) && v.D > 0) ? -log10(v.D) : nothing
    end
    error("p2ref_score: unknown indicator $X")
end

# ---------------------------------------------------------------- detection, E28. V: ten verdicts, index j = budget index j - 1
function p2ref_detection(V::AbstractVector{Symbol}, rho::Symbol, p)
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
    cens = outcome === :correct ? tdet : p.censor_factor * p.N_c
    return (t_det = tdet, outcome = outcome, censored = cens)
end

# ---------------------------------------------------------------- agreement, E30
function p2ref_table(vx, vy)
    cc = cr = rc = rr = 0
    for (a, b) in zip(vx, vy)
        (a in (:chaotic, :regular) && b in (:chaotic, :regular)) || continue
        cc += (a === :chaotic && b === :chaotic); cr += (a === :chaotic && b === :regular)
        rc += (a === :regular && b === :chaotic); rr += (a === :regular && b === :regular)
    end
    n = cc + cr + rc + rr
    A = n == 0 ? nothing : (cc + rr) / n
    kappa = nothing
    if n > 0
        pe = ((cc + cr) * (cc + rc) + (rc + rr) * (cr + rr)) / n^2
        kappa = pe == 1 ? 0.0 : (A - pe) / (1 - pe)
    end
    return (cc = cc, cr = cr, rc = rc, rr = rr, n = n, A = A, kappa = kappa)
end

# V: Dict indicator => vector of verdicts over the orbits of the set (:na where the orbit has no item), I: list of indicators
function p2ref_unanimity(V::AbstractDict, I)
    n = length(first(values(V)))
    used = 0; split = 0
    for o in 1:n
        vs = [V[X][o] for X in I]
        all(v -> v in (:chaotic, :regular), vs) || continue
        used += 1
        split += !all(==(vs[1]), vs)
    end
    pairs = [(I[a], I[b]) for a in 1:length(I) for b in (a + 1):length(I)]
    tabs = Dict(pr => p2ref_table(V[pr[1]], V[pr[2]]) for pr in pairs)
    As = [tabs[pr].A for pr in pairs]
    ks = [tabs[pr].kappa for pr in pairs]
    ok = all(!isnothing, As)
    return (U = used == 0 ? nothing : split / used, used = used, dropped = n - used, A_min = ok ? minimum(As) : nothing,
            A = Dict(pr => tabs[pr].A for pr in pairs), kappa = Dict(pr => tabs[pr].kappa for pr in pairs),
            kappa_mean = all(!isnothing, ks) ? sum(ks) / length(ks) : nothing)
end

function p2ref_midranks(x)
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

function p2ref_spearman(x, y)
    n = length(x)
    n < 3 && return nothing
    rx, ry = p2ref_midranks(x), p2ref_midranks(y)
    (all(==(rx[1]), rx) || all(==(ry[1]), ry)) && return nothing
    mx, my = sum(rx) / n, sum(ry) / n
    return sum((rx .- mx) .* (ry .- my)) / sqrt(sum((rx .- mx) .^ 2) * sum((ry .- my) .^ 2))
end

# AUC (E31): probability that a chaotic score exceeds a regular one, ties one half; nothing below min_set in either class
function p2ref_auc(c, r; min_set = 0)
    (length(c) < min_set || length(r) < min_set || isempty(c) || isempty(r)) && return nothing
    s = 0.0
    for a in c, b in r
        s += a > b ? 1.0 : a == b ? 0.5 : 0.0
    end
    return s / (length(c) * length(r))
end

# ---------------------------------------------------------------- hypothesis statistics, E32. rec = (id, V::Dict{Symbol,Vector{Symbol}}, T::Dict{Symbol,Int})
function p2ref_stat(h::Symbol, sets::Tuple, p; I = P2_X7, kp = findfirst(==(p.b_primary), p.budgets))
    nb = length(p.budgets)
    if h === :H1
        recs = sets[1]
        V = Dict(X => [haskey(r.V, X) ? r.V[X][kp] : :na for r in recs] for X in I)
        u = p2ref_unanimity(V, I)
        return (n = length(recs), used = u.used, dropped = u.dropped, U = u.U, A_min = u.A_min, A = u.A, kappa_mean = u.kappa_mean)
    elseif h === :H2
        recs = sets[1]
        M = Dict(X => p2ref_median([r.T[X] for r in recs if haskey(r.T, X)]) for X in I)
        return (n = length(recs), M = M, Q = maximum(values(M)) / minimum(values(M)))
    elseif h === :H3
        w, s = sets
        Mw = Dict(X => p2ref_median([r.T[X] for r in w if haskey(r.T, X)]) for X in I)
        Ms = Dict(X => p2ref_median([r.T[X] for r in s if haskey(r.T, X)]) for X in I)
        return (n_w = length(w), n_s = length(s), M_w = Mw, M_s = Ms, R = Dict(X => Mw[X] / Ms[X] for X in I))
    elseif h === :H4
        recs = sets[1]
        J = [X for X in I if X !== :LE]
        dec = Dict(X => count(r -> haskey(r.V, X) && r.V[X][nb] in (:chaotic, :regular), recs) for X in J)
        phi = Dict(X => dec[X] == 0 ? nothing : count(r -> haskey(r.V, X) && r.V[X][nb] === :regular, recs) / dec[X] for X in J)
        v = collect(values(phi))
        return (n = length(recs), phi = phi, decided = dec, dphi = any(isnothing, v) ? nothing : maximum(v) - minimum(v))
    elseif h === :H5
        recs = sets[1]
        dec = Dict{Symbol, Int}(); F = Dict{Symbol, Union{Nothing, Float64}}()
        for X in I
            sel = [r for r in recs if haskey(r.V, X) && r.V[X][kp] in (:chaotic, :regular) && r.V[X][nb] in (:chaotic, :regular)]
            dec[X] = length(sel)
            F[X] = isempty(sel) ? nothing : count(r -> r.V[X][kp] !== r.V[X][nb], sel) / length(sel)
        end
        v = [F[X] for X in I]
        return (n = length(recs), F = F, decided = dec, F_med = any(isnothing, v) ? nothing : p2ref_median(Float64.(v)))
    end
    error("p2ref_stat: unknown hypothesis $h")
end

# the scalar statistics that carry an interval, by hypothesis
function p2ref_scalars(h::Symbol, s)
    h === :H1 && return Dict{Symbol, Any}(:U => s.U)
    h === :H2 && return Dict{Symbol, Any}(:Q => s.Q)
    h === :H3 && return Dict{Symbol, Any}(Symbol(:R_, X) => s.R[X] for X in keys(s.R))
    h === :H4 && return Dict{Symbol, Any}(:dphi => s.dphi)
    h === :H5 && return Dict{Symbol, Any}(:F_med => s.F_med)
end

# ---------------------------------------------------------------- bootstrap stream (E12 rule, written again from the amendment text)
function p2_stream_seed(label::AbstractString, seed::Integer)
    bytes = hex2bytes(an_sha256("chaoslab-analysis-1|$(seed)|$(label)"))
    return reinterpret(UInt64, reverse(bytes[1:8]))[1]
end
p2_rng(label, seed) = StableRNG(p2_stream_seed(label, seed))

# replicate r draws m indices uniformly with replacement (E34); the draws of a label depend on the label and the replicate only
p2_draws(label, m, B, seed) = (rng = p2_rng(label, seed); [rand(rng, 1:m) for _ in 1:m, _ in 1:B])

# Interval of every scalar statistic of H1..H5. labels: one string, or two (W then S) for H3.
function p2ref_bootstrap(h::Symbol, sets::Tuple, labels, p; I = P2_X7, seed = p.bootstrap_seed)
    B = p.bootstrap_resamples
    draws = h === :H3 ? (p2_draws(labels[1], length(sets[1]), B, seed), p2_draws(labels[2], length(sets[2]), B, seed)) :
                        (p2_draws(labels, length(sets[1]), B, seed),)
    reps = Dict{Symbol, Vector{Float64}}()
    dropped = 0
    for r in 1:B
        rs = Tuple(sets[j][draws[j][:, r]] for j in eachindex(sets))
        sc = p2ref_scalars(h, p2ref_stat(h, rs, p; I = I))
        if any(isnothing, values(sc))
            dropped += 1
            continue
        end
        for (k, v) in sc
            push!(get!(reps, k, Float64[]), Float64(v))
        end
    end
    a = (1 - p.bootstrap_level) / 2
    ok = dropped <= 0.05 * B
    iv = Dict(k => (lo = p2ref_percentile(v, a), hi = p2ref_percentile(v, 1 - a)) for (k, v) in reps)
    return (interval = ok ? iv : nothing, dropped = dropped, defined = ok)
end

# ---------------------------------------------------------------- decisions, E32
function p2ref_decide(h::Symbol, s, p)
    sup(c) = c ? :supported : :not_supported
    if h === :H1
        (s.n < p.min_set || s.U === nothing || s.U_lo === nothing || s.A_min === nothing) && return :not_decidable
        return sup(s.U >= p.h1_u_min && s.U_lo >= p.h1_u_lower_min && s.A_min <= p.h1_amin_max)
    elseif h === :H2
        (s.n < p.min_set || s.Q_lo === nothing) && return :not_decidable
        return sup(s.Q >= p.h2_q_min && s.Q_lo > p.h2_q_lower_gt)
    elseif h === :H3
        (s.n_w < p.min_set || s.n_s < p.min_set) && return :not_decidable
        any(isnothing, s.R_lo) && return :not_decidable
        return sup(count(i -> s.R[i] >= p.h3_r_min && s.R_lo[i] > p.h3_r_lower_gt, eachindex(s.R)) >= p.h3_count_min)
    elseif h === :H4
        (s.n < p.min_set || s.dphi === nothing || s.dphi_lo === nothing) && return :not_decidable
        return sup(s.dphi >= p.h4_dphi_min && s.dphi_lo >= p.h4_dphi_lower_min)
    elseif h === :H5
        (s.n < p.min_set || s.F_med === nothing || s.F_med_lo === nothing) && return :not_decidable
        return sup(s.F_med >= p.h5_fmed_min && s.F_med_lo >= p.h5_fmed_lower_min)
    elseif h === :H6
        da = s.dir.n >= p.min_set; db = s.prec.n >= p.min_set
        a = da && !isempty(s.dir.delta) && maximum(values(s.dir.delta)) >= p.h6_flip_min
        b = db && !isempty(s.prec.delta) && maximum(values(s.prec.delta)) >= p.h6_flip_min
        (a || b) && return :supported
        return (da && db) ? :not_supported : :not_decidable
    end
end

# flip count Phi_X (M5, E32 H5): changes between consecutive indices k = 1..9 where both verdicts are decided
p2ref_flips(V) = count(k -> V[k] in (:chaotic, :regular) && V[k + 1] in (:chaotic, :regular) && V[k] !== V[k + 1], 1:(length(V) - 1))

# ---------------------------------------------------------------- the digest rule (E19, E37): 12 significant digits, LF line ends
p2_fmt(x::AbstractFloat) = isfinite(x) ? string(round(x; sigdigits = 12)) : string(x)
p2_fmt(x) = string(x)
p2_table_lines(cols, rows) = vcat("| " * join(cols, " | ") * " |", "|" * repeat("---|", length(cols)), ["| " * join(p2_fmt.(r), " | ") * " |" for r in rows])
p2_digest(cols, rows) = an_sha256(join(p2_table_lines(cols, rows), "\n"))

# deep comparison with a relative tolerance for reals (dictionaries by key, collections by position, nothing equal to nothing)
function p2_same(u, v; rtol = 1e-12)
    if u isa AbstractDict
        return v isa AbstractDict && keys(u) == keys(v) && all(k -> p2_same(u[k], v[k]; rtol), keys(u))
    elseif u isa NamedTuple
        return v isa NamedTuple && propertynames(u) == propertynames(v) && all(k -> p2_same(getproperty(u, k), getproperty(v, k); rtol), propertynames(u))
    elseif u isa Union{Tuple, AbstractVector}
        return (v isa Union{Tuple, AbstractVector}) && length(u) == length(v) && all(i -> p2_same(u[i], v[i]; rtol), eachindex(u))
    elseif u isa Real && v isa Real
        return isequal(u, v) || isapprox(u, v; rtol = rtol, atol = 1e-15)
    end
    return isequal(u, v)
end
# compare the fields of the reference tuple `b` with the same fields of the module's tuple `a`; returns the list of differing fields
p2_diff(a, b; skip = ()) = [k for k in propertynames(b) if k ∉ skip && !(hasproperty(a, k) && p2_same(getproperty(a, k), getproperty(b, k)))]
