# Amendment D93, stage 02. Independent reference implementation of the p3-abstention calibration layer (spec E50 to E57), written from the amendment text and used as
# the oracle of the tests (spec O15 in Julia instead of Python). It never calls Chaoslab.Analysis, so it works before the module exists.
# Names: p3ref_*. Exact rationals (BigInt) for every count-based quantity, as I32 requires. Verdict symbols :chaotic :regular :undecided (an_p2_ref.jl).
# Output symbols of the rule: :chaotic_supported :regular_supported :unresolved.

# ---------------------------------------------------------------- the rational of a decimal text (I32): 207.6 is 2076 // 10, never the nearest binary fraction
p3ref_rat(x::Integer) = Rational{BigInt}(x)
function p3ref_rat(x::Real)
    s = string(Float64(x))
    m = match(r"^(-?)(\d+)(?:\.(\d+))?(?:e([+-]?\d+))?$", s)
    m === nothing && error("p3ref_rat: cannot read $s")
    frac = m.captures[3] === nothing ? "" : m.captures[3]
    num = parse(BigInt, m.captures[2] * frac) // big(10)^length(frac)
    m.captures[1] == "-" && (num = -num)
    e = m.captures[4] === nothing ? 0 : parse(Int, m.captures[4])
    return Rational{BigInt}(e >= 0 ? num * big(10)^e : num // big(10)^(-e))
end

# ---------------------------------------------------------------- the search space (E52)
p3ref_names(p) = Symbol.(p.voter_bits)
p3ref_count(p) = sum(binomial(length(p.voter_bits), k) * k * (k - k ÷ 2) * length(p.w_range) * length(p.v_range) for k in p.min_voters:length(p.voter_bits))
function p3ref_theta(p)
    names = p3ref_names(p); n = length(names)
    out = NamedTuple[]
    for sz in p.min_voters:n
        for m in sort!([m for m in 0:(2^n - 1) if count_ones(m) == sz])
            V = Symbol[names[i] for i in 1:n if (m >> (i - 1)) & 1 == 1]
            for qC in 1:sz, qR in (sz ÷ 2 + 1):sz, w in p.w_range, v in p.v_range
                push!(out, (index = length(out) + 1, V = V, mask = m, qC = qC, qR = qR, w = w, v = v))
            end
        end
    end
    return out
end

# ---------------------------------------------------------------- the rule (E51, E52). V: Dict indicator => verdicts at the ten budgets; k = budget index 0 to 9
function p3ref_stable(vs, k, w)
    k - w >= 0 || error("budget b_k / 2^$w is not an element of B")
    c = vs[k + 1]
    return c !== :undecided && all(j -> vs[k + 1 - j] == c, 1:w)
end
function p3ref_evidence(V, k, th)
    nC = nR = nA = 0
    for X in th.V
        c = V[X][k + 1]
        s = p3ref_stable(V[X], k, th.w)
        c === :chaotic && (nA += 1; s && (nC += 1))
        c === :regular && s && (nR += 1)
    end
    return (nC = nC, nR = nR, nCany = nA)
end
function p3ref_out(V, k, th)
    e = p3ref_evidence(V, k, th)
    e.nC >= th.qC && return :chaotic_supported
    (e.nCany <= th.v && e.nR >= th.qR) && return :regular_supported
    return :unresolved
end
# R_ab: unresolved if some X in V has decided, different Float64 and Double64 verdicts at b_k; Vd = nothing is not defined
function p3ref_out_ab(V, Vd, k, th)
    Vd === nothing && error("R_ab is not defined without Double64 items")
    out = p3ref_out(V, k, th)
    any(X -> (a = V[X][k + 1]; b = Vd[X][k + 1]; a !== :undecided && b !== :undecided && a != b), th.V) && return :unresolved
    return out
end

# ---------------------------------------------------------------- cost model (E53)
p3ref_spec(p, X) = p.spec_of[string(X)]
p3ref_c(p, spec) = p3ref_rat(getfield(p, Symbol("cost_", spec)))
p3ref_CR(p, V, b) = b * sum(p3ref_c(p, s) for s in unique(p3ref_spec(p, X) for X in V))
p3ref_cap(p) = minimum(p3ref_c(p, p3ref_spec(p, X)) for X in p3ref_names(p)) * p.N_c
function p3ref_matched(p, V)
    CR = p3ref_CR(p, V, p.b_primary)
    out = Dict{Symbol, NamedTuple}()
    for X in vcat(p3ref_names(p), :LE)
        c = p3ref_c(p, p3ref_spec(p, X))
        i = findfirst(b -> b * c >= CR, p.budgets)
        j = findlast(b -> b * c <= CR, p.budgets)
        out[X] = (b = i === nothing ? nothing : p.budgets[i], capped = i === nothing, round_down = j === nothing ? nothing : p.budgets[j], round_down_none = j === nothing, cost = c)
    end
    le = out[:LE]
    b_read = le.capped ? p.N_c : min(le.b, p.N_c)
    return (CR = CR, X = out, le_read = b_read, le_unmatched = b_read * le.cost < CR)
end

# ---------------------------------------------------------------- C3 to C6 (E54). orbits: rows (id, source, kindex, in_cal, rho, S2, Ts, W2, gate, V)
p3ref_kP(p) = findfirst(==(p.b_primary), p.budgets) - 1
function p3ref_stats(rows, p, th)
    k = p3ref_kP(p)
    cal = [o for o in rows if o.in_cal]
    gate = [o for o in rows if o.gate]
    out(o) = p3ref_out(o.V, k, th)
    dec = 0; contra = 0; unres = 0
    for o in cal
        r = out(o)
        r === :unresolved && (unres += 1; continue)
        dec += 1
        (r === :chaotic_supported && o.rho === :R || r === :regular_supported && o.rho === :C) && (contra += 1)
    end
    function share(set, f)
        s = [o for o in cal if set(o)]
        isempty(s) && return nothing
        return count(f, s) // length(s)
    end
    s2 = [o for o in cal if o.S2]
    return (e = dec == 0 ? nothing : contra // dec, decided = dec, contradicting = contra, u_cal = isempty(cal) ? nothing : unres // length(cal),
            u_s2 = isempty(s2) ? nothing : count(o -> out(o) === :unresolved, s2) // length(s2),
            fr_ts = share(o -> o.Ts, o -> out(o) === :regular_supported), fr_w2 = share(o -> o.W2, o -> out(o) === :regular_supported),
            dec_ts = share(o -> o.Ts, o -> out(o) !== :unresolved), dec_w2 = share(o -> o.W2, o -> out(o) !== :unresolved),
            g = count(o -> out(o) === :chaotic_supported, gate), cost = p3ref_CR(p, th.V, p.b_primary))
end

function p3ref_feasible(st, p)
    cap = p3ref_cap(p)
    le(a, b) = a !== nothing && a <= p3ref_rat(b)
    ge(a, b) = a !== nothing && a >= p3ref_rat(b)
    return (st.cost <= cap, le(st.e, p.e_cal), le(st.u_s2, p.s_cal), st.g == 0, le(st.u_cal, p.a_cal), ge(st.dec_w2, p.c_min) && ge(st.dec_ts, p.c_min))
end

# the order of C6: lowest C_R, smallest |V|, w, v, bitmask, q_C, q_R
p3ref_key(th, cost) = (cost, length(th.V), th.w, th.v, th.mask, th.qC, th.qR)

function p3ref_select(rows, p)
    th = p3ref_theta(p)
    cal = [o for o in rows if o.in_cal]
    sizes = (Ts = count(o -> o.Ts, cal), W2 = count(o -> o.W2, cal), S2 = count(o -> o.S2, cal))
    stats = [p3ref_stats(rows, p, t) for t in th]
    flags = [p3ref_feasible(s, p) for s in stats]
    first_fail = [(i = findfirst(!, f); i === nothing ? 0 : i) for f in flags]
    alone = ntuple(i -> count(f -> f[i], flags), 6)
    F = [i for i in eachindex(th) if first_fail[i] == 0]
    base = (thetas = th, stats = stats, flags = flags, first_fail = first_fail, alone = alone, feasible = F, sizes = sizes)
    if any(<(p.min_set), values(sizes))
        return merge(base, (outcome = :sets_too_small, J = Dict{Int, Rational{BigInt}}(), J_min = nothing, tol = Int[], star = nothing))
    end
    if isempty(F)
        return merge(base, (outcome = :no_feasible_rule, J = Dict{Int, Rational{BigInt}}(), J_min = nothing, tol = Int[], star = nothing))
    end
    J = Dict{Int, Rational{BigInt}}(i => (stats[i].fr_ts + stats[i].fr_w2) // 2 for i in F)
    Jmin = minimum(values(J))
    tol = [i for i in F if J[i] <= Jmin + p3ref_rat(p.delta_tol)]
    star = tol[argmin([p3ref_key(th[i], stats[i].cost) for i in tol])]
    return merge(base, (outcome = :selected, J = J, J_min = Jmin, tol = tol, star = star))
end
# ties of the key tuple are impossible (the tuple identifies the element)
p3ref_order(th, stats, ids) = sort(ids; by = i -> p3ref_key(th[i], stats[i].cost))

# ---------------------------------------------------------------- relabelling (E49 c). d: (lam5, lam5_tenth, lam1, rho1, rho2, diag5, diag1) -> rho
function p3ref_label(d, q)
    (d.diag5 || d.diag1) && return (instability = :unresolved, structure = :unresolved, rho = :U)
    r = d.lam5_tenth > 0 ? 10 * d.lam5 / d.lam5_tenth : nothing
    inst = r === nothing ? :unresolved : (d.lam5 >= q.label_lambda_min && r >= q.label_growth_ratio) ? :positive : (d.lam5 < q.label_lambda_min && r < q.label_growth_ratio) ? :none : :unresolved
    str = (inst === :none && abs(d.rho1 - d.rho2) <= q.label_rotation_tol) ? :quasiperiodic : inst === :positive ? :irregular : :unresolved
    if abs(d.lam1 - d.lam5) > q.label_estimator_bound
        inst = str = :unresolved
    end
    return (instability = inst, structure = str, rho = (inst === :positive && str === :irregular) ? :C : (inst === :none && str === :quasiperiodic) ? :R : :U)
end

# ---------------------------------------------------------------- exact tail of the McNemar test, rounding, Hanley and McNeil (E56, E57)
p3ref_mcnemar(n10, n01) = (m = n10 + n01; m == 0 ? Rational{BigInt}(1) : sum(binomial(big(m), j) for j in n01:m; init = big(0)) // big(2)^m)
p3ref_round_away(x::Rational) = x >= 0 ? floor(Int, x + 1 // 2) : -floor(Int, -x + 1 // 2)
p3ref_ceil_step(x::Rational, s::Integer) = s * Int(ceil(x / s))
function p3ref_hm(A, na, nn)
    A = Rational{BigInt}(A); Q1 = A / (2 - A); Q2 = 2 * A^2 / (1 + A)
    return (A * (1 - A) + (na - 1) * (Q1 - A^2) + (nn - 1) * (Q2 - A^2)) / (na * nn)
end
# Mann-Whitney: failures against non-failures, a tie counts one half
function p3ref_auc(f, n)
    (isempty(f) || isempty(n)) && return nothing
    s = sum(x -> sum(y -> x > y ? 2 : x == y ? 1 : 0, n), f)
    return Rational{Int}(s, 2 * length(f) * length(n))
end

# ---------------------------------------------------------------- the draw path (E56 d)
function p3ref_path(p)
    step = p3ref_rat(p.rho_path_second) - p3ref_rat(p.rho_path_start)
    out = NamedTuple[]
    j = 0
    while true
        rho = p3ref_rat(p.rho_path_start) + j * step
        rho > p3ref_rat(p.r_max) && break
        nT = p3ref_round_away(rho * p.y_t); nRW = p3ref_round_away(rho * (p.y_w - p.y_t))
        push!(out, (j = j, rho = rho, nT = nT, nRW = nRW, nW = nT + nRW, allowed = nT >= p.n_t_min && nT + nRW >= p.n_w_min))
        j += 1
    end
    return out
end

# ---------------------------------------------------------------- generators (E55). SplitMix64 fills the state of xoshiro256++; the seed is the first 8 bytes of SHA-256, big-endian
mutable struct P3RefGen
    s::NTuple{4, UInt64}
end
function p3ref_splitmix(seed::UInt64)
    s = seed; out = UInt64[]
    for _ in 1:4
        s += 0x9e3779b97f4a7c15
        z = s
        z = (z ⊻ (z >> 30)) * 0xbf58476d1ce4e5b9
        z = (z ⊻ (z >> 27)) * 0x94d049bb133111eb
        push!(out, z ⊻ (z >> 31))
    end
    return Tuple(out)
end
function p3ref_next!(g::P3RefGen)
    s0, s1, s2, s3 = g.s
    r = bitrotate(s0 + s3, 23) + s0
    t = s1 << 17
    s2 ⊻= s0; s3 ⊻= s1; s1 ⊻= s2; s0 ⊻= s3; s2 ⊻= t; s3 = bitrotate(s3, 45)
    g.s = (s0, s1, s2, s3)
    return r
end
# SHA-256 in plain Julia (FIPS 180-4), so that a stream seed costs microseconds instead of a process; the fixture self-test compares it with the platform tool
const P3_K = UInt32[0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74,
    0x80deb1fe, 0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da, 0x983e5152, 0xa831c66d, 0xb00327c8,
    0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1,
    0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070, 0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f,
    0x682e6ff3, 0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2]
function p3ref_sha256(s::AbstractString)
    msg = Vector{UInt8}(codeunits(s))
    L = length(msg)
    push!(msg, 0x80)
    while length(msg) % 64 != 56; push!(msg, 0x00); end
    append!(msg, reverse(reinterpret(UInt8, [UInt64(8 * L)])))
    h = UInt32[0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19]
    rotr(x, n) = (x >> n) | (x << (32 - n))
    for blk in 1:64:length(msg)
        w = zeros(UInt32, 64)
        for i in 1:16
            w[i] = reduce((a, b) -> (a << 8) | b, UInt32.(msg[blk + 4 * (i - 1):blk + 4 * (i - 1) + 3]))
        end
        for i in 17:64
            s0 = rotr(w[i - 15], 7) ⊻ rotr(w[i - 15], 18) ⊻ (w[i - 15] >> 3)
            s1 = rotr(w[i - 2], 17) ⊻ rotr(w[i - 2], 19) ⊻ (w[i - 2] >> 10)
            w[i] = w[i - 16] + s0 + w[i - 7] + s1
        end
        a, b, c, d, e, f, g, hh = h
        for i in 1:64
            S1 = rotr(e, 6) ⊻ rotr(e, 11) ⊻ rotr(e, 25)
            ch = (e & f) ⊻ (~e & g)
            t1 = hh + S1 + ch + P3_K[i] + w[i]
            S0 = rotr(a, 2) ⊻ rotr(a, 13) ⊻ rotr(a, 22)
            maj = (a & b) ⊻ (a & c) ⊻ (b & c)
            t2 = S0 + maj
            hh, g, f, e, d, c, b, a = g, f, e, d + t1, c, b, a, t1 + t2
        end
        h .+= UInt32[a, b, c, d, e, f, g, hh]
    end
    return join(lpad.(string.(h; base = 16), 8, '0'))
end

p3ref_unit(u) = u isa Tuple ? join(u, ",") : string(u)
p3ref_seed(test_seed, purpose, unit) = parse(UInt64, p3ref_sha256("chaoslab-analysis-1|$(test_seed)|p3|$(purpose)|$(p3ref_unit(unit))")[1:16]; base = 16)
p3ref_stream(test_seed, purpose, unit) = P3RefGen(p3ref_splitmix(p3ref_seed(test_seed, purpose, unit)))
p3ref_uniform!(g) = Float64(p3ref_next!(g) >> 11 + 0.5) / 2.0^53
# an integer in 1..m: modulo with rejection (a reading pinned by the tests, plan "Readings pinned"): draw until x < m floor(2^64 / m), return x mod m + 1
function p3ref_randint!(g, m::Integer)
    lim = typemax(UInt64) - (typemax(UInt64) % UInt64(m) + 1) % UInt64(m)
    while true
        x = p3ref_next!(g)
        x <= lim && return Int(x % UInt64(m)) + 1
    end
end

# ---------------------------------------------------------------- PP for HB: one simulated run (E56 c, e). model: see an_p3_helpers.jl (p3_hb_model)
# the event "false-regular" is u < probability (a reading pinned by the tests)
function p3ref_hb_orbit!(g, st, model)
    R = model.rule[st]
    u = [p3ref_uniform!(g) for _ in 1:9]
    rfr = u[1] < R.p
    comp = Dict{Symbol, Bool}()
    for (i, X) in enumerate(P2_X7)
        c = model.comp[X][st]
        q = rfr ? c.s : (R.p >= 1 ? 0.0 : max(c.p - c.s * R.p, 0.0) / (1 - R.p))
        comp[X] = u[1 + i] < q
    end
    dprob = R.p >= 1 ? 1.0 : clamp((R.d - R.p) / (1 - R.p), 0.0, 1.0)
    return (rfr = rfr, comp = comp, decided = rfr || u[9] < dprob)
end
function p3ref_hb_run(test_seed, point, i, model, nT, nRW, p; alpha = p3ref_rat(p.alpha_hb), rho_hb = p3ref_rat(p.rho_hb), c_min = p3ref_rat(p.c_min))
    g = p3ref_stream(test_seed, "pp-hb", (point, i))
    sets = Dict(X => Dict(:T => zeros(Int, 4), :RW => zeros(Int, 4)) for X in P2_X7)   # n11, n10, n01, n00
    dec = Dict(:T => 0, :RW => 0)
    for (st, n) in ((:T, nT), (:RW, nRW)), _ in 1:n
        o = p3ref_hb_orbit!(g, st, model)
        o.decided && (dec[st] += 1)
        for X in P2_X7
            idx = o.rfr ? (o.comp[X] ? 1 : 2) : (o.comp[X] ? 3 : 4)
            sets[X][st][idx] += 1
        end
    end
    function ok(cnt)
        n11, n10, n01, _ = cnt
        m = n10 + n01; a = n10 + n11; c = n01 + n11
        return m >= 1 && p3ref_mcnemar(n10, n01) <= alpha && c >= 1 && a <= rho_hb * c
    end
    good = all(X -> ok(sets[X][:T]), model.ts) && all(X -> ok(sets[X][:T] .+ sets[X][:RW]), model.w2)
    cov = dec[:T] >= c_min * nT && dec[:T] + dec[:RW] >= c_min * (nT + nRW)
    return good && cov
end
p3ref_hb_power(test_seed, point, model, nT, nRW, p, runs) = count(i -> p3ref_hb_run(test_seed, point, i, model, nT, nRW, p), 1:runs)

# exact power of a single comparison on a single stratum with every orbit decided (O16 a): the sum over (n11, n10, n01, n00) of the multinomial probability
function p3ref_exact_power(n::Integer, p_rule, pX, s, alpha, rho_hb)
    p, pX, s = Rational{BigInt}(p_rule), Rational{BigInt}(pX), Rational{BigInt}(s)
    q = p >= 1 ? 0 // 1 : max(pX - s * p, 0 // 1) / (1 - p)
    P = (p * s, p * (1 - s), (1 - p) * q, (1 - p) * (1 - q))
    total = Rational{BigInt}(0)
    for n11 in 0:n, n10 in 0:(n - n11), n01 in 0:(n - n11 - n10)
        n00 = n - n11 - n10 - n01
        m = n10 + n01; a = n10 + n11; c = n01 + n11
        (m >= 1 && c >= 1 && p3ref_mcnemar(n10, n01) <= alpha && a <= rho_hb * c) || continue
        coef = factorial(big(n)) ÷ (factorial(big(n11)) * factorial(big(n10)) * factorial(big(n01)) * factorial(big(n00)))
        total += coef * P[1]^n11 * P[2]^n10 * P[3]^n01 * P[4]^n00
    end
    return total
end


# the optimism of E56 a by the reference: resamples of Cal (listed by increasing id) from the streams (optimism, b); C3 to C6 on each with the floor min_set not applied (P66);
# an empty F skips the resample; the mean of FR(theta_b; Cal set) - FR(theta_b; resample set) summed in increasing b
function p3ref_optimism(rows, q)
    cal = sort([o for o in rows if o.in_cal]; by = o -> o.id)
    n = length(cal)
    th = p3ref_theta(q)
    thetas = Union{Int, Nothing}[]; sizes = NamedTuple[]
    sT = 0.0; sW = 0.0; used = 0
    for b in 1:q.b_opt
        g = p3ref_stream(q.test_seed, "optimism", b)
        res = [cal[p3ref_randint!(g, n)] for _ in 1:n]
        push!(sizes, (Ts = count(o -> o.Ts, res), W2 = count(o -> o.W2, res)))
        sel = p3ref_select(res, merge(q, (; min_set = 0)))
        if sel.star === nothing
            push!(thetas, nothing); continue
        end
        push!(thetas, sel.star)
        full = p3ref_stats(cal, q, th[sel.star])
        sT += Float64(full.fr_ts) - Float64(sel.stats[sel.star].fr_ts); sW += Float64(full.fr_w2) - Float64(sel.stats[sel.star].fr_w2)
        used += 1
    end
    return (used = used, skipped = q.b_opt - used, thetas = thetas, sizes = sizes, optimism = used == 0 ? nothing : (T_s = sT / used, W2 = sW / used))
end

