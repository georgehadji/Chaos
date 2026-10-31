# Stage 04 verification of the analysis module, Amendment D93 (study p3-abstention), oracle O16 (criteria D93k, D93l). Run from the repo root:
#   julia --project=verify -t 23 verify/analysis-p3-o16.jl
# Spec: workspaces/module-build/stages/01-spec/output/analysis-spec.md, "Amendment D93", E56 c to f, E57, O16 a to d.
# The exact powers are summed in this script (BigInt integers over a common denominator, i.e. exact rationals), from the model of E56 c and the decision of E56 e; the module is
# called only for the simulations (p3_simulate_hb, p3_hb_type1, p3_simulate_hc). The HC reference simulation of O16 d uses Random.Xoshiro and `randn`, not the module's streams.
# One printed line per check:  <id> | quantity | value | bound | PASS/FAIL
# Nothing in data/, src/, the studies or the specs is written.
using Printf, Random, Chaoslab

const AN = Chaoslab.Analysis
const COMPS = AN.P2_I7                       # LE, SALI, FLI, MEGNO, E1, E2, WBA
const I6 = COMPS[2:end]
const ALPHA, RHO, CMIN = 5 // 100, 1 // 2, 1 // 2
const RUNS = parse(Int, get(ENV, "O16_RUNS", "400000"))   # E56 e / O16 a (the env override is for smoke runs only)
const HCRUNS = parse(Int, get(ENV, "O16_HCRUNS", "200000"))
const TH = Threads.nthreads()
# the thresholds of the table (alpha_hb, rho_hb, c_min, alpha_ci, delta_c, z_one_sided, nonfail_per_fail, test_seed): the values of the register
const Q = (; test_seed = 20261030, alpha_hb = 0.05, rho_hb = 0.5, c_min = 0.5, alpha_ci = 0.05, delta_c = 0.1, z_one_sided = 1.645, nonfail_per_fail = 20)
const DELTA_C = 1 // 10

nfail = Ref(0)
function emit(id, what, value, bound, ok)
    ok || (nfail[] += 1)
    println(id, " | ", what, " | ", value, " | ", bound, " | ", ok ? "PASS" : "FAIL"); flush(stdout)
end
f4(x) = @sprintf("%.5f", x)

# ------------------------------------------------------------------------------------------------- exact power (O16 a)
const FACT = [factorial(big(i)) for i in 0:130]
const TAIL = Dict{Tuple{Int, Int}, Bool}()
# exact one-sided tail P(Bin(m, 1/2) >= n01) <= alpha, in integers
okmc(m, n01) = get!(TAIL, (m, n01)) do
    t = sum(binomial(big(m), j) for j in n01:m; init = big(0))
    t * denominator(ALPHA) <= numerator(ALPHA) * big(2)^m
end
# n orbits of one stratum, every orbit decided: the rule is false-regular with probability p; X with probability s given that, and so with probability max(pX - s p, 0) jointly
# with "rule not"; cell probabilities (n11, n10, n01, n00) = (p s, p (1 - s), max(pX - s p, 0), the rest). Success iff m = n10 + n01 >= 1, tail <= alpha, c = n01 + n11 >= 1, a = n10 + n11 <= rho c.
function exactpow(n, p, pX, s)
    p, pX, s = Rational{BigInt}(p), Rational{BigInt}(pX), Rational{BigInt}(s)
    P3 = max(pX - s * p, 0 // 1)
    P = [p * s, p * (1 - s), P3, 1 - p - P3]
    D = lcm(denominator.(P)...)
    N = [BigInt(numerator(x * D)) for x in P]
    pw = [[N[k]^e for e in 0:n] for k in 1:4]
    tot = big(0)
    for n11 in 0:n, n10 in 0:(n - n11), n01 in 0:(n - n11 - n10)
        n00 = n - n11 - n10 - n01
        m = n10 + n01; a = n10 + n11; c = n01 + n11
        (m >= 1 && c >= 1 && okmc(m, n01) && a * denominator(RHO) <= numerator(RHO) * c) || continue
        tot += FACT[n + 1] ÷ (FACT[n11 + 1] * FACT[n10 + 1] * FACT[n01 + 1] * FACT[n00 + 1]) * pw[1][n11 + 1] * pw[2][n10 + 1] * pw[3][n01 + 1] * pw[4][n00 + 1]
    end
    return tot // D^n
end

model(; rule_T, rule_RW = rule_T, comp, ts, w2) = (rule = Dict(:T => rule_T, :RW => rule_RW), comp = comp, ts = ts, w2 = w2)
zero_c = (T = (p = 0.0, s = 0.0), RW = (p = 0.0, s = 0.0))

println("# O16 a: exact power of one comparison (one stratum, every orbit decided), simulated power at $RUNS runs, $TH threads")
sets = [(12, 1 // 20, 3 // 5, 1 // 2), (16, 1 // 10, 1 // 2, 1 // 2), (20, 0 // 1, 1 // 4, 1 // 2), (25, 1 // 20, 2 // 5, 1 // 2), (30, 1 // 10, 1 // 2, 1 // 4), (40, 1 // 10, 1 // 4, 1 // 2),
        (40, 1 // 10, 3 // 10, 1 // 2), (50, 1 // 5, 1 // 2, 1 // 2), (70, 1 // 5, 2 // 5, 1 // 2), (80, 1 // 10, 3 // 10, 0 // 1), (100, 1 // 5, 3 // 10, 1 // 2), (110, 1 // 5, 2 // 5, 1 // 4),
        (120, 3 // 10, 1 // 2, 1 // 2), (120, 1 // 10, 1 // 5, 1 // 2), (120, 1 // 10, 1 // 4, 1 // 2), (120, 1 // 20, 1 // 5, 1 // 2)]
inrange = 0
for (i, (n, p, pX, s)) in enumerate(sets)
    ex = Float64(exactpow(n, p, pX, s))
    global inrange += 0.1 <= ex <= 0.99
    m = model(rule_T = (p = Float64(p), d = 1.0), comp = Dict(X => (X == :SALI ? (T = (p = Float64(pX), s = Float64(s)), RW = (p = 0.0, s = 0.0)) : zero_c) for X in COMPS), ts = [:SALI], w2 = Symbol[])
    r = AN.p3_simulate_hb(m, n, 0, Q; point = i, runs = RUNS, threads = TH)
    se = sqrt(ex * (1 - ex) / RUNS)
    z = (r.power - ex) / se
    emit("O16a-$i", "n=$n p=$p pX=$pX s=$s: exact $(f4(ex)), simulated $(f4(r.power)), z = (sim - exact) / SE", @sprintf("%.2f", z), "|z| <= 3", abs(z) <= 3 && 0.1 <= ex <= 0.99)
end
emit("O16a", "parameter sets with exact power in [0.1, 0.99], n up to 120", "$inrange of $(length(sets))", ">= 12", inrange >= 12)

# ------------------------------------------------------------------------------------------------- joint power of the 13 comparisons (O16 b)
println("# O16 b: joint power of the 13 comparisons (six of I_6 on T_s, seven of I_7 on W2) against min single and Bonferroni, same models")
# strata T_s (nT orbits) and R_W (nRW orbits) with the same rates, so W2 is nT + nRW orbits with those rates; every orbit decided
joint_models = [
    ("A", 1 // 20, 30, 40, Dict(:LE => (3 // 5, 1 // 2), :SALI => (1 // 2, 1 // 2), :FLI => (3 // 5, 1 // 4), :MEGNO => (1 // 2, 1 // 2), :E1 => (2 // 5, 1 // 2), :E2 => (3 // 5, 1 // 2), :WBA => (1 // 2, 1 // 4))),
    ("B", 1 // 10, 40, 40, Dict(:LE => (1 // 2, 1 // 2), :SALI => (2 // 5, 1 // 2), :FLI => (1 // 2, 1 // 4), :MEGNO => (2 // 5, 1 // 2), :E1 => (3 // 10, 1 // 2), :E2 => (1 // 2, 1 // 2), :WBA => (2 // 5, 1 // 4))),
    ("C", 1 // 20, 24, 24, Dict(:LE => (7 // 10, 1 // 2), :SALI => (3 // 5, 1 // 2), :FLI => (1 // 2, 1 // 2), :MEGNO => (3 // 5, 1 // 4), :E1 => (7 // 10, 1 // 2), :E2 => (3 // 5, 1 // 2), :WBA => (1 // 2, 1 // 2))),
]
for (k, (name, p, nT, nRW, par)) in enumerate(joint_models)
    comp = Dict(X => (T = (p = Float64(par[X][1]), s = Float64(par[X][2])), RW = (p = Float64(par[X][1]), s = Float64(par[X][2]))) for X in COMPS)
    m = model(rule_T = (p = Float64(p), d = 1.0), comp = comp, ts = I6, w2 = COMPS)
    singles = Float64[]
    for X in I6; push!(singles, Float64(exactpow(nT, p, par[X]...))); end              # the six on T_s
    for X in COMPS; push!(singles, Float64(exactpow(nT + nRW, p, par[X]...))); end    # the seven on W2
    r = AN.p3_simulate_hb(m, nT, nRW, Q; point = 100 + k, runs = RUNS, threads = TH)
    upper = minimum(singles) + 3 * r.se
    lower = 1 - sum(1 .- singles) - 3 * r.se
    emit("O16b-$name", "13 comparisons, nT=$nT nRW=$nRW: joint simulated $(f4(r.power)) (SE $(f4(r.se))), min single $(f4(minimum(singles)))", f4(r.power), "<= $(f4(upper))", r.power <= upper)
    emit("O16b-$name", "Bonferroni 1 - sum(1 - single) = $(f4(1 - sum(1 .- singles))) over $(length(singles)) comparisons", f4(r.power), ">= $(f4(lower))", r.power >= lower)
end

# ------------------------------------------------------------------------------------------------- type I (O16 c)
println("# O16 c: type I, rule rate equal to the comparator rate")
for (i, (n, p, s)) in enumerate([(20, 1 // 5, 1 // 2), (40, 1 // 10, 1 // 2), (60, 3 // 10, 1 // 4), (120, 1 // 20, 1 // 2), (120, 1 // 5, 0 // 1), (12, 1 // 2, 1 // 2)])
    pf, sf = Float64(p), Float64(s)
    m = model(rule_T = (p = pf, d = 1.0), comp = Dict(X => (T = (p = pf, s = sf), RW = (p = pf, s = sf)) for X in COMPS), ts = [:SALI], w2 = Symbol[])
    t = AN.p3_hb_type1(m, n, 0, Q; point = 200 + i, runs = RUNS, threads = TH)
    ex = Float64(exactpow(n, p, p, s))
    se = sqrt(ex * (1 - ex) / RUNS)
    r = t.rates[:SALI]
    emit("O16c-$i", "n=$n p=pX=$p s=$s: single-comparison rejection rate $(f4(r.single)) (exact $(f4(ex)))", f4(r.single), "< alpha_hb = 0.05", r.single < Q.alpha_hb)
    emit("O16c-$i", "same rate against the exact rejection probability, in SE", @sprintf("%.2f", (r.single - ex) / se), "|z| <= 3", abs(r.single - ex) <= 3 * se)
end
# the joint type-I check of E56 f on a 13-comparison model with two strata: the largest joint rate over X is at most alpha_hb plus 3 SE
let (name, p, nT, nRW, par) = joint_models[2]
    comp = Dict(X => (T = (p = Float64(par[X][1]), s = Float64(par[X][2])), RW = (p = Float64(par[X][1]), s = Float64(par[X][2]))) for X in COMPS)
    m = model(rule_T = (p = Float64(p), d = 1.0), comp = comp, ts = I6, w2 = COMPS)
    t = AN.p3_hb_type1(m, nT, nRW, Q; point = 300, runs = RUNS, threads = TH)
    se = sqrt(Q.alpha_hb * (1 - Q.alpha_hb) / RUNS)
    emit("O16c-J", "13-comparison model B, rule rate set to each X in turn: max joint rejection rate $(f4(t.max_joint))", f4(t.max_joint), "<= $(f4(Q.alpha_hb + 3 * se))", t.max_joint <= Q.alpha_hb + 3 * se)
end

# ------------------------------------------------------------------------------------------------- PP for HC (O16 d)
println("# O16 d: PP for HC, $HCRUNS runs per point")
# own normal quantile (Newton on the Taylor series of erf; |x| < 2 only): the shift d(A) = sqrt(2) Phi^{-1}(A)
function Phi(x)

    erf = 0.0; z = x / sqrt(2)
    for k in 0:80
        erf += (-1)^k * z^(2k + 1) / (factorial(big(k)) * (2k + 1))
    end
    return 0.5 * (1 + 2 / sqrt(pi) * Float64(erf))
end
function qnorm(A)
    x = 0.0
    for _ in 1:60
        x -= (Phi(x) - A) / (exp(-x^2 / 2) / sqrt(2pi))
    end
    return x
end
shift(A) = sqrt(2) * qnorm(A)

function hm(A, na, nn)
    Q1 = A / (2 - A); Q2 = 2A^2 / (1 + A)
    return (A * (1 - A) + (na - 1) * (Q1 - A^2) + (nn - 1) * (Q2 - A^2)) / (na * nn)
end
# Mann-Whitney U (failures above non-failures; continuous scores, no ties)
function ucount(f, nonf)
    nf = length(f)
    sp = sortperm(vcat(f, nonf))
    R = 0
    for (pos, idx) in enumerate(sp); idx <= nf && (R += pos); end
    return R - nf * (nf + 1) ÷ 2
end
# independent simulation: run r has its own Xoshiro; the decision is Delta >= delta_c (exact, integers) and Delta - z SE > 0
function indep_hc(Ai, Ap, nf, runs; seed = 20261031)
    nn = Q.nonfail_per_fail * nf
    di, dp = shift(Ai), shift(Ap)
    hits = zeros(Int, TH)
    Threads.@threads for t in 1:TH
        for r in t:TH:runs
            rng = Random.Xoshiro(seed + 1_000_003 * nf + r)
            fI = di .+ randn(rng, nf); nI = randn(rng, nn); fP = dp .+ randn(rng, nf); nP = randn(rng, nn)
            UI, UP = ucount(fI, nI), ucount(fP, nP)
            AI, AP = UI / (nf * nn), UP / (nf * nn)
            se = sqrt(hm(AI, nf, nn) + hm(AP, nf, nn))
            ((UI - UP) * denominator(DELTA_C) >= numerator(DELTA_C) * nf * nn && (AI - AP) - Q.z_one_sided * se > 0) && (hits[t] += 1)
        end
    end
    return sum(hits)
end

for (i, (nf, a)) in enumerate([(20, 0.55), (49, 0.55), (65, 0.55), (30, 0.5)])
    r = AN.p3_simulate_hc(Q, a, a; nf = nf, runs = HCRUNS, threads = TH)
    bound = Q.alpha_ci + 3 * sqrt(Q.alpha_ci * (1 - Q.alpha_ci) / HCRUNS)
    emit("O16d-0-$i", "true AUC difference 0 (AUC $a both), n_f=$nf: rate of the decision", f4(r.power), "<= alpha_ci + 3 SE = $(f4(bound))", r.power <= bound)
    pj = indep_hc(a, a, nf, HCRUNS) / HCRUNS
    se = sqrt(r.power * (1 - r.power) / HCRUNS + pj * (1 - pj) / HCRUNS)
    emit("O16d-0x-$i", "same point, module rate $(f4(r.power)) against the independent simulation $(f4(pj)) (Xoshiro, randn)", @sprintf("%.5f", abs(r.power - pj)), "<= 3 SE = $(f4(3 * se))", abs(r.power - pj) <= 3 * se)
end
for (i, (Ai, Ap, nf)) in enumerate([(0.8, 0.5, 6), (0.8, 0.5, 10), (0.8, 0.5, 15), (0.7, 0.55, 49), (0.7, 0.55, 65), (0.75, 0.55, 30)])
    r = AN.p3_simulate_hc(Q, Ai, Ap; nf = nf, runs = HCRUNS, threads = TH)
    h = indep_hc(Ai, Ap, nf, HCRUNS)
    pj = h / HCRUNS
    se = sqrt(r.power * (1 - r.power) / HCRUNS + pj * (1 - pj) / HCRUNS)
    emit("O16d-p-$i", "AUC $Ai against $Ap, n_f=$nf: module power $(f4(r.power)), independent simulation $(f4(pj)) (Xoshiro, randn)", @sprintf("%.5f", abs(r.power - pj)), "<= 3 SE = $(f4(3 * se))", abs(r.power - pj) <= 3 * se && 0 < r.power < 1)
end

println("# O16 done: ", nfail[] == 0 ? "all checks PASS" : "$(nfail[]) FAIL")
exit(nfail[] == 0 ? 0 : 1)
