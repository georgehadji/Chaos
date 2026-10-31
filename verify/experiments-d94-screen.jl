# Included by verify/experiments-d94.jl: O24 (screening classes by the L1 estimator on the reduced real panel) and O25 (overlap of the starts).

const T_BIT, C_BIT, S_BIT = 0x01, 0x02, 0x04
function side_bits(side, i)          # the class bits of K index i from the run-length lists of the pool sidecar
    k = side["k"][string(i)]
    bits = zeros(UInt8, Int(k["points"]))
    for (c, bit) in (("T", T_BIT), ("C", C_BIT), ("S", S_BIT)), r in k[c], j in Int(r[1]):Int(r[2])
        bits[j] |= bit
    end
    return bits
end
const EPS53 = eps(1.0)
Gfirst(K, x) = (c = K * cos(x); (1 + c)^2 + 1 + c^2 + 1)
BS(n, K, G) = 8 * EPS53 * (n + G) + 16 * abs(K) * n^2 * EPS53            # B_S of lyap-benettin S3 as amended by D35 (the O5 bound, as test AC25)

# L1 (Benettin reducer of propagate, the verified estimator) at the horizons w, 2w, ..., N: S(j w) = sum + comp of the final state of each run
function l1_series(K, pts, N, w, alpha, seed)
    sys = StandardMap(K)
    x0 = [SVector{2, Float64}(p[1], p[2]) for p in pts]
    dir = SMatrix{2, 1, Float64}(cos(alpha), sin(alpha))
    nb = N ÷ w
    S = zeros(length(pts), nb)
    lam = zeros(length(pts))
    for j in 1:nb
        r = propagate(sys, x0, j * w; k = 1, tau = 1, seed = seed, reducers = (Benettin(),), W0 = fill(dir, length(pts)))
        for m in eachindex(pts)
            f = r.finals[m][1]
            S[m, j] = f.sum[1] + f.comp[1]
            j == nb && (lam[m] = f.lambda[1])
        end
    end
    return S, lam
end

# class bits from (lambda, S(N/10), S(N), block exponents) by the predicates of E11; dT, dE the distances to the periodic set and its elliptic subset
function bits_from(lam, s10, sN, blocks, dT, dE, inp)
    f = 0x00
    XP.pred_T(lam, s10, sN, dT; theta_T = inp.theta_T, g_T = inp.g_T, r_T = inp.r_T) && (f |= T_BIT)
    XP.pred_C(lam, dE; theta_C = inp.theta_C, r_C = inp.r_C) && (f |= C_BIT)
    XP.pred_S(lam, blocks; theta_S = inp.theta_S, w = inp.w, f_S = inp.f_S, L_S = inp.L_S) && (f |= S_BIT)
    return f
end

# a mismatch of the bits is "near a threshold" when some quantity of the predicates lies within the O5 error of its threshold
function near_threshold(lam, s10, sN, blocks, dT, dE, inp, tolS)
    tl = tolS / inp.N_scr
    near = abs(lam - inp.theta_T) <= tl || abs(lam - inp.theta_C) <= tl || abs(lam - inp.theta_S) <= tl || abs(s10) <= tolS
    near |= abs(sN / s10 - inp.g_T) <= inp.g_T * (tolS / abs(s10) + tolS / abs(sN))
    tb = 2 * tolS / inp.w
    med = XP.block_median(blocks)
    near |= abs(med) <= tb || any(e -> abs(e - inp.f_S * med) <= tb * (1 + inp.f_S), blocks)
    near |= abs(dT - inp.r_T) <= 1e-12 || abs(dE - inp.r_C) <= 1e-12
    return near
end

# one point: L5 against L1 under the O5 bound, then the class bits from L5 and from L1 against the sidecar bit; returns (lambda_ok, L5 bits ok, L1 verdict)
function o24_point(i, K, id, x, p, Sm, l1, bit, ps, inp)
    N = inp.N_scr
    r = XP.screen_point(K, x, p; N_scr = N, w = inp.w, alpha = inp.alpha, tau = 1)
    tolS = 2 * BS(N, K, Gfirst(K, x))
    dT, dE = XP.periodic_distance((x, p), ps), XP.periodic_distance((x, p), ps; elliptic_only = true)
    s10 = Sm[N ÷ 10 ÷ inp.w]
    blocks = diff(vcat(0.0, Sm)) ./ inp.w
    b1 = bits_from(l1, s10, Sm[end], blocks, dT, dE, inp)
    b5 = bits_from(r.lambda, r.S_tenth, r.S_N, r.blocks, dT, dE, inp)
    verdict = b1 == bit ? :equal : near_threshold(l1, s10, Sm[end], blocks, dT, dE, inp, tolS) ? :near : :failure
    verdict == :equal || say("O24 |   $(verdict == :near ? "near-threshold case" : "FAILURE") | K index $i point $id lambda_L1 $l1 lambda_L5 $(r.lambda) bits L1 $b1 L5 $b5 sidecar $bit | $(verdict == :near ? "listed, not counted" : "counted")")
    return abs(r.lambda - l1) / (tolS / N), b5 == bit, verdict
end

section("O24") do
    dir = joinpath(V, "reduced")
    inp = XP.read_prereg(joinpath(dir, "p3-reduced-prereg.md")).inputs
    side = TOML.parsefile(joinpath(dir, "p3-reduced-panel.toml.pools.toml"))
    rng = StableRNG(20261031)
    nK = length(inp.K)
    ratio, nl5, nnear, nfail, ntot = 0.0, 0, 0, 0, 0
    for i in 1:nK
        K = inp.K[i]
        pts = XP.random_screen_points(inp.seed, i, Int(side["k"][string(i)]["points"]))
        ids = sort!(randperm(rng, length(pts))[1:200])                           # 200 random point ids per K
        bits = side_bits(side, i)
        ps = setprecision(BigFloat, inp.bigfloat_bits) do
            XP._search_k(inp, i, Threads.nthreads())
        end
        t = @elapsed S, l1 = l1_series(K, pts[ids], inp.N_scr, inp.w, inp.alpha, 1)
        res = [o24_point(i, K, id, pts[id][1], pts[id][2], S[m, :], l1[m], bits[id], ps, inp) for (m, id) in enumerate(ids)]
        ratio = max(ratio, maximum(first, res)); nl5 += count(r -> !r[2], res); ntot += length(res)
        nnear += count(r -> r[3] == :near, res); nfail += count(r -> r[3] == :failure, res)
        cnt(b) = count(id -> bits[id] & b != 0, ids)
        say("O24 |   K index $i (K = $K): $(length(ids)) points, L1 series in $(round(t; digits = 1)) s, class counts in the sample T/C/S = $(cnt(T_BIT))/$(cnt(C_BIT))/$(cnt(S_BIT)) | -")
    end
    emit("O24", "points checked (200 random point ids per K index, $nK K)", ntot, "$(200 * nK)", pf(ntot == 200 * nK))
    emit("O24", "largest |lambda_L5 - lambda_L1| as a multiple of the O5 bound 2 B_S / N_scr", ratio, "<= 1", pf(ratio <= 1))
    emit("O24", "class bits recomputed from L5 (screening values, predicates of E11) against the sidecar, differing points", nl5, "0", pf(nl5 == 0))
    emit("O24", "class bits recomputed from the L1 values against the sidecar: differing points near a threshold (listed, not counted)", nnear, "reported", "INFO")
    emit("O24", "class bits recomputed from the L1 values against the sidecar: differing points not near a threshold (failures)", nfail, "0", pf(nfail == 0))
end

# ---------------------------------------------------------------------------------------------------
# O25: no start of the reduced panel equals a start of the P1 panel (= the p2 panel, reused) or a p2 grid start (h4 tables)
pair(x0, p0) = (Float64(x0) + 0.0, Float64(p0) + 0.0)         # + 0.0 maps -0.0 to 0.0: equality as binary64 values
section("O25") do
    red = TOML.parsefile(joinpath(V, "reduced", "p3-reduced-panel.toml"))["orbit"]
    p1 = TOML.parsefile(P1_PANEL)["orbit"]
    txt = replace(read(P2_PREREG, String), "\r\n" => "\n")
    blk = [m.captures[1] for m in eachmatch(r"```toml\n(.*?)\n```"s, txt) if occursin("[[run.h4]]", m.captures[1])]
    h4 = TOML.parse(only(blk))["run"]["h4"]
    emit("O25", "p2 grid starts read from the h4 tables of the p2 preregistration", length(h4), "2304", pf(length(h4) == 2304))
    S1 = [pair(o["x0"], o["p0"]) for o in p1]
    Sg = [pair(parse(Float64, h["x0"]), parse(Float64, h["p0"])) for h in h4]
    Sr = [pair(o["x0"], o["p0"]) for o in red]
    emit("O25", "starts of the reduced panel (orbits) / distinct pairs among them (I29)", "$(length(Sr)) / $(length(unique(Sr)))", "equal", pf(length(Sr) == length(unique(Sr))))
    e1, eg = count(in(Set(S1)), Sr), count(in(Set(Sg)), Sr)
    emit("O25", "reduced-panel starts equal (binary64 pair) to a start of the P1 panel ($(length(S1)) orbits, the p2 panel is this file)", e1, "0", pf(e1 == 0))
    emit("O25", "reduced-panel starts equal (binary64 pair) to a p2 grid start ($(length(Sg)))", eg, "0", pf(eg == 0))
    d1 = [minimum(tdist(s, q) for q in S1) for s in Sr]
    dg = [minimum(tdist(s, q) for q in Sg) for s in Sr]
    emit("O25", "smallest torus distance of the reduced-panel starts to the P1 panel starts (any K), at orbit id $(red[argmin(d1)]["id"])", minimum(d1), "reported", "INFO")
    emit("O25", "same, restricted to the same K (K index of both)", minimum(minimum(tdist(pair(o["x0"], o["p0"]), q) for (q, oq) in zip(S1, p1) if oq["k_index"] == o["k_index"]; init = Inf) for o in red if !o["control"]), "reported", "INFO")
    emit("O25", "smallest torus distance of the reduced-panel starts to the p2 grid starts", minimum(dg), "reported", "INFO")
    info("O25", "p2 panel starts: the p2 panel is the P1 panel file (p2 preregistration, Panel: 'the frozen P1 panel P (reused)'); no other p2 panel file exists under workspaces/, nothing under data/ was read")
end
