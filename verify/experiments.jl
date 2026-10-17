# Stage 04 verification of module experiments (study panel builder, decision D49, amendment D51).
# Run from the repo root:  julia --project=verify -t 23 verify/experiments.jl
# Spec: workspaces/module-build/stages/01-spec/output/experiments-spec.md ("Oracle" O6, O7; criteria 2, 12, 19, 21, 23, 27; D51c).
# One printed line per check:  <id> | quantity | value | bound | PASS/FAIL/INFO
# Environment variables (debugging aids): EXV_ONLY="AC27,O7" runs only the named sections (default all),
# EXV_REUSE=1 reuses the full panel already in the scratch folder instead of rebuilding it (about 7 min on 23 threads).
using Printf, Random, LinearAlgebra, Statistics, SHA
using Chaoslab, StaticArrays, ForwardDiff, StableRNGs
import TOML
import ChaosTools as CT
import CairoMakie as CM

const EX = Chaoslab.Experiments
const ROOT = dirname(@__DIR__)
const OUT = joinpath(ROOT, "workspaces", "module-build", "stages", "04-verify", "output")
const STUDY_OUT = joinpath(ROOT, "workspaces", "study", "stages", "02-panel", "output")
const PREREG = joinpath(ROOT, "workspaces", "study", "stages", "01-design", "output", "p1-le-benchmark-preregistration.md")
const SCR = get(ENV, "EXV_SCRATCH", raw"C:\Users\tesse\AppData\Local\Temp\claude\C--Users-tesse-Documents-Vibe-Coding-Chaos\5df15ee1-47df-4875-85bd-1bd4d7c665a6\scratchpad\verify-panel")
const PANEL = joinpath(SCR, "panel.toml")
mkpath(OUT); mkpath(SCR)
const ONLY = filter(!isempty, split(get(ENV, "EXV_ONLY", ""), ","))
const REUSE = get(ENV, "EXV_REUSE", "") == "1"
const LOG = open(get(ENV, "EXV_LOG", joinpath(OUT, "experiments-verify-log.txt")), "w")
const TWO_PI = 2 * pi

fmt(x::AbstractFloat) = @sprintf("%.4g", x)
fmt(x) = string(x)
const COUNTS = Dict("PASS" => 0, "FAIL" => 0, "INFO" => 0)
function say(s)
    println(s); println(LOG, s); flush(stdout); flush(LOG)
end
function emit(ac, quantity, value, bound, verdict)
    v = startswith(verdict, "FAIL") ? "FAIL" : startswith(verdict, "PASS") ? "PASS" : "INFO"
    COUNTS[v] += 1
    say(string(ac, " | ", quantity, " | ", fmt(value), " | ", bound, " | ", verdict))
end
pf(ok::Bool) = ok ? "PASS" : "FAIL"
info(ac, msg) = (COUNTS["INFO"] += 1; say("$ac | INFO | $msg"))

function section(f, name)
    isempty(ONLY) || name in ONLY || return
    t0 = time()
    try
        f()
    catch e
        emit(name, "ERROR", sprint(showerror, e, catch_backtrace())[1:min(end, 600)], "-", "FAIL")
    end
    say(@sprintf("# section %s took %.1f s", name, time() - t0))
end

tdelta(s) = abs(s - TWO_PI * round(s / TWO_PI))
tdist(u, v) = hypot(tdelta(u[1] - v[1]), tdelta(u[2] - v[2]))
sha_file(path) = bytes2hex(open(sha256, path))

say("# experiments verification, julia $(VERSION), threads $(Threads.nthreads()), ChaosTools 3.5.4 / DynamicalSystemsBase 3.19.3 (verify/Manifest.toml)")
say("# prereg $(relpath(PREREG, ROOT)), scratch $SCR")
const T0 = time()

# ---------------------------------------------------------------------------------------------------
# Panel state shared by the sections
const PR = EX.read_prereg(PREREG)           # throws on a hash mismatch: the whole script depends on the frozen file
const INP = PR.inputs
const CLASS_ORDER = ["T", "R", "C", "S", "X"]
const KIDX_2PI = findfirst(==(2pi), INP.K)
PANEL_DATA = Ref{Any}(nothing)
panel() = (PANEL_DATA[] === nothing && (PANEL_DATA[] = TOML.parsefile(PANEL)); PANEL_DATA[])
BUILD_SECONDS = Ref(NaN)
FULL_OK = Ref(false)

# ---------------------------------------------------------------------------------------------------
section("BUILD") do
    if REUSE && isfile(PANEL)
        info("BUILD", "reusing the panel in the scratch folder (EXV_REUSE=1), build time not measured")
    else
        isfile(PANEL) && rm(PANEL)
        t = @elapsed EX.build_panel(PREREG; out = PANEL)
        BUILD_SECONDS[] = t
        emit("D51c", "build_panel on the real preregistration completes, threads $(Threads.nthreads()), wall seconds (compile included)", t, "completes", "PASS")
    end
    FULL_OK[] = isfile(PANEL)
    P = panel()
    man = P["manifest"]
    emit("BUILD", "orbits (drawn, controls)", "$(man["totals"]["orbits"]) ($(man["totals"]["drawn"]), $(man["totals"]["controls"]))", "-", "INFO")
    emit("BUILD", "panel file bytes and SHA-256", "$(filesize(PANEL)) $(sha_file(PANEL))", "-", "INFO")
    emit("BUILD", "manifest git sha / dirty flag / package version", "$(man["git_sha"][1:10]) / $(man["git_dirty"]) / $(man["package_version"])", "-", "INFO")
    emit("BUILD", "Newton statistics (tried, accepted, singular, diverged, cap, ill_conditioned)",
         join((man["newton"][k] for k in ("tried", "accepted", "singular", "diverged", "cap", "ill_conditioned")), ", "), "-", "INFO")
end

# ---------------------------------------------------------------------------------------------------
# AC2 [full]: hash of the real frozen file; O3 independent implementation (platform tool, own SELF rule)
function selfrule_bytes(text)
    lines = split(replace(text, "\r\n" => "\n"), "\n")
    idx = findall(l -> startswith(l, "| SHA-256 of this file"), lines)
    length(idx) == 1 || error("expected one hash row, found $(length(idx))")
    cells = split(lines[idx[1]], "|")
    cells[3] = " SELF "                        # the whole value cell
    lines[idx[1]] = join(cells, "|")
    return join(lines, "\n")
end
function platform_sha(text)
    tmp = replace(joinpath(SCR, "selfrule.bin"), "\\" => "/")
    write(tmp, selfrule_bytes(text))
    try
        return lowercase(lstrip(split(read(`sha256sum $tmp`, String))[1], '\\')), "sha256sum"
    catch
        out = read(`certutil -hashfile $tmp SHA256`, String)
        return lowercase(replace(split(out, "\n")[2], " " => "", "\r" => "")), "certutil"
    end
end
section("AC2") do
    text = read(PREREG, String)
    chk = EX.check_prereg(PREREG)
    plat, tool = platform_sha(text)
    emit("AC2a", "frozen file: header SHA-256 equals module value", chk.header_sha256 == chk.computed_sha256, "true", pf(chk.header_sha256 == chk.computed_sha256))
    emit("AC2a/e", "module value equals the independent value ($tool, own SELF rule)", chk.computed_sha256 == plat, "true", pf(chk.computed_sha256 == plat))
    info("AC2", "header $(chk.header_sha256), amendments logged $(chk.amendments), panel block parsed: K = $(INP.K), grid_side $(INP.grid_side), N_scr $(INP.N_scr)")
    P = panel()
    man = P["manifest"]
    pm = man["prereg"]
    emit("AC2", "panel manifest records the same header and computed hash", pm["header_sha256"] == chk.header_sha256 && pm["computed_sha256"] == chk.computed_sha256, "true",
         pf(pm["header_sha256"] == chk.header_sha256 && pm["computed_sha256"] == chk.computed_sha256))
    # CRLF copy passes (AC2d) on the real file
    crlf = joinpath(SCR, "prereg-crlf.md")
    write(crlf, replace(replace(text, "\r\n" => "\n"), "\n" => "\r\n"))
    c2 = EX.check_prereg(crlf)
    emit("AC2d", "CRLF copy of the real file passes with the same hash", c2.computed_sha256 == chk.computed_sha256, "true", pf(c2.computed_sha256 == chk.computed_sha256))
end

# ---------------------------------------------------------------------------------------------------
# D51c [full]
section("D51c") do
    P = panel(); man = P["manifest"]
    k5 = KIDX_2PI
    Kb = man["K_binary64"][k5]
    emit("D51c", "K index $k5 is K = 2pi in binary64", Kb, "== 2pi (Float64)", pf(Kb == 2pi))
    ill = [r for r in man["ill_conditioned"] if r["k_index"] == k5]
    total = sum(r["count"] for r in ill; init = 0)
    emit("D51c", "ill-conditioned count at K = 2pi (all periods)", total, "> 0", pf(total > 0))
    for r in sort(ill; by = r -> r["q"])
        info("D51c", "  period $(r["q"]): $(r["count"]) rejected, smallest |2 - t| $(fmt(r["min_abs_2_minus_t"]))")
    end
    others = [r for r in man["ill_conditioned"] if r["k_index"] != k5]
    emit("D51c", "ill-conditioned rows at the other K values (reported)", length(others), "reported", "INFO")
    rows1 = [r for r in man["periodic"] if r["k_index"] == k5 && r["q"] == 1]
    best = Inf
    bad = 0
    setprecision(BigFloat, 256) do
        for r in rows1, pt in r["points"]
            z = (Float64(parse(BigFloat, pt[1])), Float64(parse(BigFloat, pt[2])))
            d = tdist(z, (pi / 2, 0.0))
            best = min(best, d)
            if d < 1e-5          # the module's own bound would have to pass: recompute e = ||F|| ||J - I|| / |2 - t| in 256 bits
                x, p = parse(BigFloat, pt[1]), parse(BigFloat, pt[2])
                Kbig = BigFloat(Kb)
                p2 = p + Kbig * sin(x); x2 = x + p2
                F = (x2 - x - TWO_PI * round(Float64(x2 - x) / TWO_PI), p2 - p - TWO_PI * round(Float64(p2 - p) / TWO_PI))
                c = Kbig * cos(x); t = 2 + c
                e = hypot(Float64(F[1]), Float64(F[2])) * hypot(Float64(c), 1.0, Float64(c), 0.0) / abs(Float64(2 - t))
                e <= INP.delta_dup / 10 || (bad += 1)
            end
        end
    end
    emit("D51c", "accepted period-1 rows at K = 2pi / nearest accepted period-1 point to (pi/2, 0) (torus distance)", "$(length(rows1)) / $(fmt(best))", ">= 1e-5 unless bound passes", pf(bad == 0))
    emit("D51c", "smallest separation of distinct orbit points (all K) against 10 delta_dup", man["min_separation"], ">= $(10 * INP.delta_dup)", pf(man["min_separation"] >= 10 * INP.delta_dup))
    ok = !haskey(man, "merged") && !haskey(man["inputs"], "sn_trace") && all(r["class"] in ("elliptic", "parabolic", "hyperbolic") for r in man["periodic"])
    emit("D51d", "no D50 output remains (merged table, sn_trace input, near-saddle-node class)", ok, "true", pf(ok))
end

# ---------------------------------------------------------------------------------------------------
# AC12 [full]: K = 0 screening exponent over the full grid, closed form (O2), T cell at K = 0
section("AC12") do
    N, w, tau, alpha = INP.N_scr, INP.w, INP.tau, INP.alpha
    cen = EX.grid_centres(INP.grid_side)
    n = length(cen)
    lam = zeros(n); s10 = zeros(n); sN = zeros(n)
    t = @elapsed Threads.@threads :greedy for tid in 1:Threads.nthreads()
        buf = EX.screen_buffer(N, w)
        for i in tid:Threads.nthreads():n
            EX.screen_point!(buf, 0.0, cen[i][1], cen[i][2], N, w, alpha, tau)
            lam[i] = buf.lambda; s10[i] = buf.S_tenth; sN[i] = buf.S_N
        end
    end
    ref = setprecision(BigFloat, 256) do
        a = BigFloat(alpha); c, s = cos(a), sin(a)
        log(sqrt((c + N * s)^2 + s^2)) / N
    end
    reff = Float64(ref)
    err = maximum(abs, lam .- reff)
    emit("AC12", "K = 0: $(n) grid cells screened in $(round(t; digits = 1)) s; closed-form lambda_scr (256 bit) $(fmt(reff)); max |lambda_scr - closed form|", err, "<= $(fmt(8 * eps(1.0)))", pf(err <= 8 * eps(1.0)))
    emit("AC12", "same, relative to the closed form", err / reff, "reported", "INFO")
    emit("AC12", "spread of lambda_scr over the grid (max - min)", maximum(lam) - minimum(lam), "reported", "INFO")
    # T cell at K = 0, recomputed from the closed-form screening values and the lines p = 2 pi a / b
    lines = [(a, b) for b in 1:INP.q_max for a in 0:(b - 1) if gcd(a, b) == 1]
    dline(p) = minimum(tdelta(p - TWO_PI * a / b) for (a, b) in lines)
    q = count(i -> lam[i] <= INP.theta_T && s10[i] > 0 && sN[i] / s10[i] < INP.g_T && dline(cen[i][2]) >= INP.r_T, 1:n)
    cell = only(c for c in panel()["manifest"]["cells"] if c["k_index"] == 1 && c["class"] == "T")
    emit("AC12", "K = 0 T cell: qualifying count recomputed here / in the panel manifest", "$q / $(cell["qualifying"])", "equal", pf(q == cell["qualifying"]))
    g = sN[1] / s10[1]
    info("AC12", "K = 0 T cell in the panel: status $(cell["status"]), drawn $(cell["drawn"]), reason \"$(cell["reason"])\"; lambda_scr $(fmt(lam[1])) <= theta_T $(INP.theta_T), growth ratio $(fmt(g)) < g_T $(INP.g_T)")
    emit("AC12", "K = 0 other classes recorded undefined (C, S, X)", all(c["status"] == "undefined" for c in panel()["manifest"]["cells"] if c["k_index"] == 1 && c["class"] in ("C", "S", "X")), "true",
         pf(all(c["status"] == "undefined" for c in panel()["manifest"]["cells"] if c["k_index"] == 1 && c["class"] in ("C", "S", "X"))))
end

# ---------------------------------------------------------------------------------------------------
# AC21 [full], cell table, AC22 / AC23 file checks
section("AC21") do
    P = panel(); man = P["manifest"]; orbs = P["orbit"]
    inp = INP
    cellno(o) = (o["k_index"] - 1) * 5 + findfirst(==(o["class"]), CLASS_ORDER)
    members = [(id = o["id"], cell = o["control"] ? 0 : cellno(o), rank = o["control"] ? 0 : o["draw_index"], control = o["control"]) for o in orbs]
    ctrl = count(o -> o["control"], orbs)
    drawn_by_cell = Dict{Int, Int}()
    for o in orbs
        o["control"] || (drawn_by_cell[cellno(o)] = get(drawn_by_cell, cellno(o), 0) + 1)
    end
    expD = sum(min(n, inp.s_D) for n in values(drawn_by_cell)) + ctrl
    expL = sum(min(n, inp.s_L) for n in values(drawn_by_cell))
    SD = sort([o["id"] for o in orbs if o["in_SD"]]); SL = sort([o["id"] for o in orbs if o["in_SL"]])
    ids = Set(o["id"] for o in orbs)
    emit("AC21", "S_L within S_D within P", issubset(SL, SD) && issubset(SD, ids), "true", pf(issubset(SL, SD) && issubset(SD, ids)))
    emit("AC21", "|S_D| / cap_D, |S_L| / cap_L", "$(length(SD)) / $(inp.cap_D), $(length(SL)) / $(inp.cap_L)", "<= caps", pf(length(SD) <= inp.cap_D && length(SL) <= inp.cap_L))
    emit("AC21", "|S_D| = sum over cells of min(drawn, s_D) plus controls, before caps (expected $expD)", length(SD), "== $expD (or caps)", pf(length(SD) == expD || length(SD) == inp.cap_D))
    emit("AC21", "|S_L| = sum over cells of min(drawn, s_L) (expected $expL)", length(SL), "== $expL (or cap)", pf(length(SL) == expL || length(SL) == inp.cap_L))
    percell_ok = all(filter(o -> !o["control"], orbs)) do o
        o["in_SD"] == (o["draw_index"] <= inp.s_D) && o["in_SL"] == (o["draw_index"] <= inp.s_L)
    end
    emit("AC21", "no cap applied: every drawn orbit of rank <= s_D is in S_D, of rank <= s_L in S_L", percell_ok, "true", pf(percell_ok || man["subsets"]["capped_D"] || man["subsets"]["capped_L"]))
    emit("AC21", "controls all in S_D, none in S_L, count", "$(all(o -> o["in_SD"] && !o["in_SL"], filter(o -> o["control"], orbs))) / $ctrl", "true", pf(all(o -> o["in_SD"] && !o["in_SL"], filter(o -> o["control"], orbs))))
    r = EX.select_subsets(members; s_D = inp.s_D, cap_D = inp.cap_D, s_L = inp.s_L, cap_L = inp.cap_L)
    same = r.SD == SD && r.SL == SL && r.removed_D == man["subsets"]["removed_D"] && r.removed_L == man["subsets"]["removed_L"]
    emit("AC21", "recomputed subsets and removed ids equal the file (removed D / L: $(length(r.removed_D)) / $(length(r.removed_L)))", same, "true", pf(same))
end

section("CELLS") do
    P = panel(); man = P["manifest"]; orbs = P["orbit"]
    say("# cell table (drawn / qualifying, status):")
    say("#   K index, K | " * join(CLASS_ORDER, " | "))
    ok_all = true
    for i in eachindex(INP.K)
        row = String[]
        for cls in CLASS_ORDER
            c = only(c for c in man["cells"] if c["k_index"] == i && c["class"] == cls)
            push!(row, "$(c["drawn"]) / $(c["qualifying"]) $(c["status"])" * (isempty(c["reason"]) ? "" : " [$(c["reason"])]"))
            n = count(o -> !o["control"] && o["k_index"] == i && o["class"] == cls, orbs)
            ok = n == c["drawn"] && (c["status"] in ("populated",) ? n == INP.n_c : c["status"] == "short" ? 0 < n < INP.n_c && !isempty(c["reason"]) :
                                     c["status"] == "empty" ? n == 0 && !isempty(c["reason"]) : c["status"] == "undefined" && n == 0 && INP.K[i] == 0 && cls in ("C", "S", "X"))
            ok_all &= ok
        end
        say("#   $i, $(INP.K[i]) | " * join(row, " | "))
    end
    emit("CELLS", "every cell populated (n_c), short or empty with a reason, or undefined at K = 0 (C, S, X); counts in the file equal the orbits present", ok_all, "true", pf(ok_all))
    for cls in CLASS_ORDER
        n = count(o -> !o["control"] && o["class"] == cls, orbs)
        emit("AC23", "class $cls: orbits in the panel (over all K)", n, ">= 1", pf(n >= 1))
    end
    r2 = only(c for c in man["cells"] if c["k_index"] == KIDX_2PI && c["class"] == "R")
    emit("CELLS", "K = 2pi, class R: status $(r2["status"]), drawn $(r2["drawn"]), reason \"$(r2["reason"])\" (ACCEPTED by the owner: the four elliptic points of the D50 build were ghosts)", r2["status"], "empty accepted", "INFO")
    # AC22-type file checks
    ids_ok = [o["id"] for o in orbs] == collect(1:length(orbs))
    emit("AC22", "orbit ids run 1..N in file order", ids_ok, "true", pf(ids_ok))
    tot = man["totals"]
    emit("AC22", "totals: orbits = drawn + controls = sum of cell counts + controls", tot["orbits"] == tot["drawn"] + tot["controls"] == (sum(c["drawn"] for c in man["cells"]) + tot["controls"]) && tot["orbits"] == length(orbs), "true",
         pf(tot["orbits"] == tot["drawn"] + tot["controls"] && tot["drawn"] == sum(c["drawn"] for c in man["cells"]) && tot["orbits"] == length(orbs)))
    req = ("id", "class", "control", "control_kind", "system", "x0", "p0", "x0_dec", "p0_dec", "run_seed", "orbit_seed", "in_SD", "in_SL")
    missing_f = count(o -> !all(haskey(o, k) for k in req), orbs)
    drawn_missing = count(o -> !o["control"] && !(haskey(o, "cell_seed") && haskey(o, "draw_index") && haskey(o, "k_index")), orbs)
    emit("AC23", "orbits missing a required field / drawn orbits missing cell seed or draw index", "$missing_f / $drawn_missing", "0 / 0", pf(missing_f == 0 && drawn_missing == 0))
    # cell stream seeds re-derived with the platform-independent rule of the manifest (E15)
    seed = parse(Int, string(man["inputs"]["seed"]["value"]))
    badseed = count(orbs) do o
        o["control"] && return false
        h = sha256("chaoslab-panel-1|$(seed)|$(o["k_index"])|$(o["class"])")
        s = foldl((a, b) -> (a << 8) | UInt64(b), h[1:8]; init = UInt64(0))
        string(s) != o["cell_seed"]
    end
    emit("AC17", "cell stream seed of every drawn orbit re-derived from E15 (SHA-256 of the text)", badseed, "0 mismatches", pf(badseed == 0))
    emit("AC17", "orbit seed of every orbit equals tangent-core rule (EX.orbit_seed)", count(o -> string(EX.orbit_seed(seed, o["id"])) != o["orbit_seed"], orbs), "0 mismatches", pf(all(o -> string(EX.orbit_seed(seed, o["id"])) == o["orbit_seed"], orbs)))
end

# ---------------------------------------------------------------------------------------------------
# O7: separate process, independent re-evaluation of the panel file (no Chaoslab loaded in the child)
const O7_CODE = raw"""
using TOML, Printf
setprecision(BigFloat, 1024)
const TWO_PI = 2 * BigFloat(pi)
panel = TOML.parsefile(ARGS[1])
man, orbs = panel["manifest"], panel["orbit"]
inp = man["inputs"]
val(k) = inp[k]["value"]
rho, ddup, dR = Float64(val("rho")), Float64(val("delta_dup")), Float64(val("delta_R"))
lo, hi = Float64(val("log10_r")[1]), Float64(val("log10_r")[2])
Ng, qmax = Int(val("grid_side")), Int(val("q_max"))
Ks = [Float64(k) for k in man["K_binary64"]]
fails = Ref(0)
function line(id, what, value, bound, ok)
    ok || (fails[] += 1)
    @printf("%s | %s | %s | %s | %s\n", id, what, value, bound, ok ? "PASS" : "FAIL"); flush(stdout)
end
d256(s) = setprecision(() -> parse(BigFloat, s), BigFloat, 256)
tdel(s) = abs(s - TWO_PI * round(s / TWO_PI))
smap(K, x, p) = (p2 = p + K * sin(x); (x + p2, p2))
# 1. starts rebuilt from the decimal text at 256 bits
bad1 = 0; bad_exact = 0
for o in orbs
    x, p = d256(o["x0_dec"]), d256(o["p0_dec"])
    (Float64(x) == o["x0"] && Float64(p) == o["p0"]) || (global bad1 += 1)
    if o["class"] in ("T", "C", "S", "X") || (o["class"] == "R" && o["k_index"] >= 1 && haskey(o, "a"))
        (x == BigFloat(o["x0"]) && p == BigFloat(o["p0"])) || (global bad_exact += 1)
    end
end
line("O7", "decimal text rebuilt at 256 bits rounds to the stored binary64 start, orbits $(length(orbs))", bad1, "0 mismatches", bad1 == 0)
line("O7", "starts formed in binary64 (T, C, S, X, R at K = 0): decimal text equals the binary64 value exactly", bad_exact, "0 mismatches", bad_exact == 0)
# 2. periodic orbits at 1024 bits
prows = man["periodic"]
maxres = 0.0; nbad = 0; maxtr = 0.0; maxex = 0.0; bad_prim = 0; bad_pts = 0; bad_class = 0; bad_wind = 0; maxdet = 0.0
for r in prows
    K = BigFloat(Ks[r["k_index"]]); q = r["q"]
    pts = [(parse(BigFloat, s[1]), parse(BigFloat, s[2])) for s in r["points"]]
    z0 = pts[1]
    x, p = z0; a, b, c, d = BigFloat(1), BigFloat(0), BigFloat(0), BigFloat(1)
    iter = Tuple{BigFloat, BigFloat}[]
    zs = [(x, p)]
    for t in 1:q
        cc = K * cos(x)
        a, b, c, d = (1 + cc) * a + c, (1 + cc) * b + d, cc * a + c, cc * b + d
        x, p = smap(K, x, p)
        push!(zs, (x, p))
    end
    m, j = round(Int, Float64((x - z0[1]) / TWO_PI)), round(Int, Float64((p - z0[2]) / TWO_PI))
    res = Float64(max(abs(x - z0[1] - TWO_PI * m), abs(p - z0[2] - TWO_PI * j)))
    global maxres = max(maxres, res)
    res <= rho || (global nbad += 1)
    (m == r["m"] && j == r["j"]) || (global bad_wind += 1)
    t = a + d
    global maxdet = max(maxdet, Float64(abs(a * d - b * c - 1)))
    global maxtr = max(maxtr, Float64(abs(t - r["trace"]) / max(1, abs(t))))
    cls = abs(abs(Float64(t)) - 2) <= rho ? "parabolic" : (abs(Float64(t)) < 2 ? "elliptic" : "hyperbolic")
    cls == r["class"] || (global bad_class += 1)
    ex = abs(t) > 2 ? Float64(log((abs(t) + sqrt(t * t - 4)) / 2) / q) : 0.0
    global maxex = max(maxex, abs(ex - r["exponent"]) / max(1e-300, abs(ex) + (ex == 0 ? 1.0 : 0.0)))
    # every orbit point listed (wrapped) is an iterate of the representative
    for t2 in 0:(q - 1)
        zt = zs[t2 + 1]
        zw = (zt[1] - TWO_PI * floor(zt[1] / TWO_PI), zt[2] - TWO_PI * floor(zt[2] / TWO_PI))
        any(pt -> Float64(hypot(tdel(pt[1] - zw[1]), tdel(pt[2] - zw[2]))) < 1e-8, pts) || (global bad_pts += 1)
    end
    # primitive period: no divisor d < q returns within rho
    for dv in 1:(q - 1)
        q % dv == 0 || continue
        zd = zs[dv + 1]
        rd = max(abs(tdel(zd[1] - z0[1])), abs(tdel(zd[2] - z0[2])))
        Float64(rd) <= rho && (global bad_prim += 1)
    end
end
line("O7", "periodic rows $(length(prows)): max return residual at 1024 bits (windings of the stored representative)", maxres, "<= rho = $rho", nbad == 0)
line("O7", "windings (m, j) recomputed equal the stored ones, mismatches", bad_wind, "0", bad_wind == 0)
line("O7", "trace recomputed at 1024 bits against the stored trace, max relative", maxtr, "<= 1e-12", maxtr <= 1e-12)
line("O7", "max |det J_q - 1| at 1024 bits", maxdet, "<= 1e-200", maxdet <= 1e-200)
line("O7", "class re-derived from the trace (parabolic band rho), mismatches", bad_class, "0", bad_class == 0)
line("O7", "stored exponent against ln|mu_+| / q at 1024 bits, max error", maxex, "<= 1e-12", maxex <= 1e-12)
line("O7", "every listed orbit point is an iterate of the representative (tolerance 1e-8: points are polished to residual ~rho^2, not 256 bits), mismatches", bad_pts, "0", bad_pts == 0)
line("O7", "primitive period: a divisor of q passes the return test, rows", bad_prim, "0", bad_prim == 0)
# 3. separation of distinct orbit points per K (Float64 wrapped), sweep over x
allsep = Inf
for k in eachindex(Ks)
    pts = Tuple{Float64, Float64, Int}[]
    for r in prows
        r["k_index"] == k || continue
        for s in r["points"]
            push!(pts, (Float64(parse(BigFloat, s[1])), Float64(parse(BigFloat, s[2])), r["id"]))
        end
    end
    ext = copy(pts)
    for (x, p, id) in pts
        x < 1e-2 && push!(ext, (x + Float64(TWO_PI), p, id))
        x > Float64(TWO_PI) - 1e-2 && push!(ext, (x - Float64(TWO_PI), p, id))
    end
    sort!(ext; by = t -> t[1])
    best = Inf
    for i in eachindex(ext)
        for j in (i + 1):length(ext)
            ext[j][1] - ext[i][1] > 1e-2 && break
            ext[i][3] == ext[j][3] && continue
            best = min(best, hypot(tdel(ext[i][1] - ext[j][1]), tdel(ext[i][2] - ext[j][2])))
        end
    end
    global allsep = min(allsep, Float64(best))
end
line("O7", "smallest distance between points of distinct periodic orbits (all K, recomputed); stored value $(man["min_separation"])", allsep, ">= 10 delta_dup = $(10 * ddup)", allsep >= 10 * ddup)
line("O7", "recomputed separation equals the stored one (relative)", abs(allsep - man["min_separation"]) / man["min_separation"], "<= 1e-6", abs(allsep - man["min_separation"]) / man["min_separation"] <= 1e-6)
# 4. torus distances of the starts
col = Set(Float64(TWO_PI * (2i + 1) / (2Ng)) for i in 0:(Ng - 1))
prow = Dict(r["id"] => r for r in prows)
nX = 0; nXr = 0; nXl = 0; nR = 0; nRd = 0; nG = 0; nR0 = 0; nCtl = 0; maxXr = 0.0; maxRd = 0.0
for o in orbs
    o["control"] && continue
    x, p = d256(o["x0_dec"]), d256(o["p0_dec"])
    if o["class"] in ("T", "C", "S")
        (o["x0"] in col && o["p0"] in col) || (global nG += 1)
    elseif o["class"] == "X"
        global nX += 1
        r = Float64(hypot(tdel(x), tdel(p)))
        e = abs(r - o["r"]) / o["r"]
        global maxXr = max(maxXr, e); e <= 1e-12 || (global nXr += 1)
        (lo - 1e-12 <= log10(o["r"]) <= hi + 1e-12) || (global nXl += 1)
    elseif o["class"] == "R" && haskey(o, "a")
        global nR0 += 1
        a, b = o["a"], o["b"]
        (gcd(a, b) == 1 && 1 <= b <= qmax && o["p0"] == Float64(TWO_PI * a / b)) || (global nR += 1)
    elseif o["class"] == "R"
        global nR += 0
        src = prow[o["source_orbit"]]["points"][o["source_point"]]
        sx, sp = parse(BigFloat, src[1]), parse(BigFloat, src[2])
        d = Float64(hypot(tdel(x - sx), tdel(p - sp)))
        e = abs(d - dR) / dR
        global maxRd = max(maxRd, e); e <= 1e-12 || (global nRd += 1)
    end
end
line("O7", "T, C, S starts are screening grid centres of the N_g grid (1024-bit nearest binary64), mismatches", nG, "0", nG == 0)
line("O7", "X starts ($nX): torus distance from (0, 0) equals stored r, max relative error", maxXr, "<= 1e-12", nXr == 0)
line("O7", "X starts: log10 r in [$lo, $hi], violations", nXl, "0", nXl == 0)
line("O7", "R starts at K > 0: torus distance from the recorded periodic point equals delta_R, max relative error", maxRd, "<= 1e-12", nRd == 0)
line("O7", "R starts at K = 0 ($nR0): gcd(a, b) = 1, b <= q_max, p0 = nearest binary64 of 2 pi a / b, violations", nR, "0", nR == 0)
# 5. controls
for o in orbs
    o["control"] || continue
    o["control_kind"] == "fixed" || continue
    K = BigFloat(Ks[o["k_index"]])
    t = o["x0"] == 0 ? 2 + K : 2 - K
    cls = abs(abs(Float64(t)) - 2) <= rho ? "parabolic" : (abs(Float64(t)) < 2 ? "elliptic" : "hyperbolic")
    ex = abs(t) > 2 ? Float64(log((abs(t) + sqrt(t * t - 4)) / 2)) : 0.0
    ok = cls == o["stability"] && abs(ex - o["exponent"]) <= 1e-12 * max(1, ex)
    line("O7", "fixed-point control id $(o["id"]) (K index $(o["k_index"]), x0 = $(o["x0"] == 0 ? "0" : "pi")): class and exponent against the closed form", "$(o["stability"]) $(o["exponent"])", "$cls $(Float64(ex))", ok)
end
for o in orbs
    o["control_kind"] == "y5" || continue
    mu = (3 + sqrt(BigFloat(5))) / 2
    ok = abs(o["exponent"] - Float64(log(mu))) <= 1e-14 * Float64(log(mu))
    line("O7", "Y5 control id $(o["id"]): exponent against ln mu at 1024 bits", o["exponent"], Float64(log(mu)), ok)
end
println("O7-SUMMARY | failures | $(fails[])")
"""
section("O7") do
    f = joinpath(SCR, "o7-child.jl")
    write(f, O7_CODE)
    t = @elapsed out = read(`$(Base.julia_cmd()) --startup-file=no -t 1 $f $PANEL`, String)
    for l in split(strip(out), "\n")
        startswith(l, "O7-SUMMARY") && continue
        say(l)
        endswith(l, "| FAIL") ? (COUNTS["FAIL"] += 1) : endswith(l, "| PASS") ? (COUNTS["PASS"] += 1) : nothing
    end
    info("O7", "separate process (julia --startup-file=no, TOML and Base only, no Chaoslab), $(round(t; digits = 1)) s")
end

# ---------------------------------------------------------------------------------------------------
# AC27 / O6: ChaosTools periodic orbits against the module's periodic set
ct_rule(x, par, n) = (p = x[2] + par[1] * sin(x[1]); th = x[1] + p; mod2pi.(SVector(th, p)))
function lifted_o(z, K, o)
    x, p = z[1], z[2]
    for _ in 1:o
        p = p + K * sin(x); x = x + p
    end
    return SVector(x, p)
end
ad_trace(z, K, o) = (J = ForwardDiff.jacobian(w -> lifted_o(w, K, o), SVector{2, Float64}(z)); J[1, 1] + J[2, 2])
function trace_power(t, k)
    k == 1 && return t
    a, b = 2.0, t
    for _ in 2:k
        a, b = b, t * b - a
    end
    return b
end
cls_of(t, band) = abs(abs(t) - 2) <= band ? "parabolic" : (abs(t) < 2 ? "elliptic" : "hyperbolic")
function dedup(pts, tol)
    s = sort(pts; by = p -> p[1])
    keep = SVector{2, Float64}[]
    for p in s
        j = length(keep)
        dup = false
        while j >= 1 && p[1] - keep[j][1] <= tol
            tdist(p, keep[j]) <= tol && (dup = true; break)
            j -= 1
        end
        dup || push!(keep, p)
    end
    return keep
end
bdist(z, K, o) = (w = SVector{2, Float64}(z); m = Inf; for _ in 1:o; m = min(m, w[1], TWO_PI - w[1], w[2], TWO_PI - w[2]); w = ct_rule(w, (K,), 1); end; m)

section("AC27") do
    P = panel(); man = P["manifest"]
    ics = EX.grid_centres(64)
    inds, signs = CT.lambdaperms(2)
    lams = [0.001, 0.01, 0.1]
    Kidx = [i for i in eachindex(INP.K) if INP.K[i] > 0]
    say("# ChaosTools settings: periodicorbits, orders 1..3, ICs = centres of the 64 x 64 sub-grid ($(length(ics))), lambda in $lams, all $(length(inds)) index permutations x $(length(signs)) sign patterns, maxiters 1e5, disttol 1e-10, inftol 10.0, abstol 1e-8; davidchacklai n = 3, m = 2, package defaults")
    # task list: (K index, order, lambda, perm, sign, chunk)
    chunks = collect(Iterators.partition(1:length(ics), 256))
    tasks = [(ki, o, li, a, b, c) for ki in Kidx for o in 1:3 for li in eachindex(lams) for a in eachindex(inds) for b in eachindex(signs) for c in eachindex(chunks)]
    res = Vector{Vector{SVector{2, Float64}}}(undef, length(tasks))
    t_pc = @elapsed Threads.@threads :greedy for n in eachindex(tasks)
        ki, o, li, a, b, c = tasks[n]
        ds = CT.DeterministicIteratedMap(ct_rule, SVector(1.0, 1.0), [INP.K[ki]])
        res[n] = try        # ChaosTools throws a BoundsError when it finds no point (StateSpaceSet of an empty list)
            SVector{2, Float64}[SVector{2, Float64}(p) for p in CT.periodicorbits(ds, o, ics[chunks[c]], [lams[li]], [inds[a]], [signs[b]])]
        catch e
            e isa BoundsError ? SVector{2, Float64}[] : rethrow()
        end
    end
    info("AC27", "periodicorbits: $(length(tasks)) tasks in $(round(t_pc; digits = 1)) s")
    dl = Dict{Int, Vector{Vector{SVector{2, Float64}}}}()
    t_dl = @elapsed Threads.@threads :greedy for ki in Kidx
        ds = CT.DeterministicIteratedMap(ct_rule, SVector(1.0, 1.0), [INP.K[ki]])
        dl[ki] = [SVector{2, Float64}[SVector{2, Float64}(p) for p in v] for v in CT.davidchacklai(ds, 3, ics, 2)]
    end
    info("AC27", "davidchacklai: $(length(Kidx)) K values in $(round(t_dl; digits = 1)) s")
    total_unmatched = 0; total_mismatch = 0; total_trace_bad = 0
    for ki in Kidx
        K = INP.K[ki]
        mod_all = [(q = r["q"], id = r["id"], cls = r["class"], tr = r["trace"], pt = SVector(Float64(parse(BigFloat, s[1])), Float64(parse(BigFloat, s[2]))))
                   for r in man["periodic"] if r["k_index"] == ki for s in r["points"]]
        ill_q = Set(r["q"] for r in man["ill_conditioned"] if r["k_index"] == ki && r["count"] > 0)
        for o in 1:3
            hits_by = Dict{String, Vector{SVector{2, Float64}}}()
            for li in eachindex(lams)
                hits_by["pc_$(lams[li])"] = dedup(reduce(vcat, [res[n] for n in eachindex(tasks) if tasks[n][1] == ki && tasks[n][2] == o && tasks[n][3] == li]; init = SVector{2, Float64}[]), 1e-7)
            end
            hits_by["dl"] = dedup(dl[ki][o], 1e-7)
            hits = dedup(reduce(vcat, values(hits_by)), 1e-7)
            cand = [m for m in mod_all if o % m.q == 0]
            matched = 0; boundary = 0; ghost = 0; other = 0; mism = 0; trbad = 0; maxtr = 0.0
            unmatched_list = String[]
            mism_list = String[]
            for h in hits
                best = Inf; bm = nothing
                for m in cand
                    d = tdist(h, m.pt)
                    d < best && (best = d; bm = m)
                end
                tad = ad_trace(h, K, o)
                if best <= 1e-6
                    matched += 1
                    texp = trace_power(bm.tr, o ÷ bm.q)
                    rel = abs(tad - texp) / max(1.0, abs(texp))
                    maxtr = max(maxtr, rel)
                    rel <= 1e-5 || (trbad += 1)
                    cad = cls_of(tad, 1e-5)
                    cok = bm.cls == "parabolic" ? abs(abs(tad) - 2) <= 1e-5 : cad == bm.cls
                    cok || (mism += 1; push!(mism_list, "($(fmt(h[1])), $(fmt(h[2]))) module $(bm.cls) t_expected $(fmt(texp)), AD t $(fmt(tad))"))
                else
                    near_b = bdist(h, K, o) < 1e-3
                    near_t = abs(tad - 2) < 1e-2
                    if ki == KIDX_2PI && near_t && any(q -> o % q == 0, ill_q)
                        ghost += 1
                        length(unmatched_list) < 6 && push!(unmatched_list, "ghost ($(fmt(h[1])), $(fmt(h[2]))) nearest module point $(fmt(best)), AD t - 2 = $(fmt(tad - 2))")
                    elseif near_b
                        boundary += 1
                        push!(unmatched_list, "boundary ($(fmt(h[1])), $(fmt(h[2]))) nearest module point $(fmt(best)), AD t $(fmt(tad))")
                    else
                        other += 1
                        push!(unmatched_list, "OTHER ($(fmt(h[1])), $(fmt(h[2]))) nearest module point $(fmt(best)), AD t $(fmt(tad))")
                    end
                end
            end
            total_unmatched += other; total_mismatch += mism; total_trace_bad += trbad
            # reverse coverage and analytic share (reported, not scored)
            ownp = [m.pt for m in cand]
            found_mod = count(p -> any(h -> tdist(h, p) <= 1e-6, hits), ownp)
            ana = SVector{2, Float64}[]
            if o == 2
                append!(ana, [SVector(0.0, Float64(pi)), SVector(Float64(pi), Float64(pi))])          # family A
                if setprecision(() -> BigFloat(K) > 2 * BigFloat(pi), BigFloat, 256)
                    xa = asin(TWO_PI / K)
                    append!(ana, [SVector(xa, Float64(pi)), SVector(mod2pi(xa + pi), Float64(pi)), SVector(mod2pi(pi - xa), Float64(pi)), SVector(mod2pi(2pi - xa), Float64(pi))])   # family B
                end
            end
            anafound = isempty(ana) ? "-" : "$(count(p -> any(h -> tdist(h, p) <= 1e-6, hits), ana)) / $(length(ana))"
            counts = join(("$(k)=$(length(v))" for (k, v) in sort(collect(hits_by))), " ")
            emit("AC27", "K index $ki (K = $(fmt(K))), o = $o: ChaosTools points (union) $(length(hits)) [$counts]; matched $matched, unmatched near wrap boundary $boundary, ghost (D51) $ghost, unmatched other $other, class mismatches $mism",
                 "$matched/$(length(hits))", "100% matched, 0 mismatches", pf(other == 0 && mism == 0))
            emit("AC27", "  trace agreement of the matched points (AD of the wrapped rule vs module trace power), max relative; points above 1e-5", "$(fmt(maxtr)); $trbad", "<= 1e-5", pf(trbad == 0))
            info("AC27", "  reverse coverage (module points with period | $o found by ChaosTools) $found_mod / $(length(ownp)); analytic E6 families A, B found $anafound")
            for s in unmatched_list
                info("AC27", "    unmatched: " * s)
            end
            for s in mism_list
                info("AC27", "    class mismatch: " * s)
            end
        end
    end
    emit("AC27", "totals over K and o: unmatched other / class mismatches / trace outliers", "$total_unmatched / $total_mismatch / $total_trace_bad", "0 / 0 / 0", pf(total_unmatched == 0 && total_mismatch == 0 && total_trace_bad == 0))
end

# ---------------------------------------------------------------------------------------------------
# AC19: real thread runs on a reduced non-trivial input (copy of the real block with grid_side = 64)
section("AC19") do
    txt = replace(read(PREREG, String), "\r\n" => "\n")
    occursin("grid_side = 256", txt) || error("grid_side = 256 not found in the panel block")
    txt2 = replace(txt, "grid_side = 256" => "grid_side = 64"; count = 1)
    h = EX.prereg_hash(txt2)
    txt3 = replace(txt2, r"(\| SHA-256 of this file[^|\n]*\|[ \t]*)[0-9a-f]{64}" => SubstitutionString("\\g<1>" * h))
    red = joinpath(SCR, "prereg-reduced.md")
    write(red, txt3)
    chk = EX.check_prereg(red)
    plat, tool = platform_sha(txt3)
    emit("AC19", "reduced preregistration (grid_side 64) passes the hash rule; module hash equals independent ($tool)", chk.header_sha256 == chk.computed_sha256 == plat, "true", pf(chk.header_sha256 == chk.computed_sha256 == plat))
    code = "using Chaoslab; t = @elapsed Chaoslab.build_panel(ARGS[1]; out = ARGS[2]); println(\"BUILD_SECONDS \", t, \" threads \", Threads.nthreads())"
    shas = Dict{Int, String}()
    for nt in (1, 4, 23)
        out = joinpath(SCR, "reduced-t$nt.toml")
        isfile(out) && rm(out)
        o = read(`$(Base.julia_cmd()) --project=$(joinpath(ROOT, "verify")) --startup-file=no -t $nt -e $code $red $out`, String)
        shas[nt] = sha_file(out)
        info("AC19", "julia -t $nt: $(strip(o)), SHA-256 $(shas[nt]), $(filesize(out)) bytes")
    end
    same = length(unique(values(shas))) == 1
    emit("AC19", "panel files from julia -t 1, 4, 23 byte-identical (SHA-256)", same, "true", pf(same))
    for nt in (1, 4, 23)
        say("#   AC19 SHA-256 -t $nt: $(shas[nt])")
    end
    forbidden = ("lambda", "block", "screen", "wall", "host", "thread", "device", "backend", "elapsed")
    txtp = read(joinpath(SCR, "reduced-t23.toml"), String)
    hits = [k for k in forbidden if occursin(Regex("^\\s*\\[?\\[?[\\w.]*$k", "mi"), txtp)]
    emit("AC19", "forbidden key names in the reduced panel file (lambda, block, screen, wall, host, thread, device, backend)", join(hits, ","), "none", pf(isempty(hits)))
    # the full-input panel, built in this process at $(Threads.nthreads()) threads, against the reduced one is not comparable; report its SHA-256
    info("AC19", "full-input panel (this process, $(Threads.nthreads()) threads) SHA-256 $(sha_file(PANEL))")
end

# ---------------------------------------------------------------------------------------------------
# Portraits from portrait_content of the full-input build
section("PORTRAITS") do
    pc = EX.portrait_content(PANEL; n_plot = 1000)
    colors = Dict("T" => "#0072B2", "R" => "#E69F00", "C" => "#D55E00", "S" => "#009E73", "X" => "#CC79A7", "control" => "#000000")
    fig = CM.Figure(size = (2100, 1500), fontsize = 15)
    for (n, c) in enumerate(pc)
        r, col = fldmod1(n, 3)
        ax = CM.Axis(fig[r, col]; title = "k index $(c.k_index), K = $(fmt(c.K))", subtitle = (isempty(c.notes) ? "all cells populated" : join(c.notes, "; ")) * "; " * string(count(q -> q.symbol == :cross, c.periodic)) * " hyperbolic periodic points not drawn",
                     subtitlesize = 11, aspect = CM.DataAspect(), limits = (0, TWO_PI, 0, TWO_PI), xlabel = "x", ylabel = "p",
                     xticks = ([0, pi, TWO_PI], ["0", "π", "2π"]), yticks = ([0, pi, TWO_PI], ["0", "π", "2π"]))
        for cls in vcat(CLASS_ORDER, ["control"])
            its = [it for it in c.iterates if it.class == cls]
            isempty(its) && continue
            xs = Float64[x[1] for it in its for x in it.xy]; ps = Float64[x[2] for it in its for x in it.xy]
            CM.scatter!(ax, xs, ps; markersize = 1.6, color = (colors[cls], 0.55))
        end
        for cls in CLASS_ORDER
            st = [s for s in c.starts if s.class == cls]
            isempty(st) || CM.scatter!(ax, [s.x for s in st], [s.p for s in st]; marker = :star5, markersize = 9, color = colors[cls], strokecolor = :black, strokewidth = 0.6)
        end
        ct = [s for s in c.starts if s.class == "control"]
        isempty(ct) || CM.scatter!(ax, [s.x for s in ct], [s.p for s in ct]; marker = :rect, markersize = 11, color = :transparent, strokecolor = :black, strokewidth = 1.6)
        for (sym, mk) in ((:circle, :circle), (:diamond, :diamond))      # hyperbolic periodic points are counted in the subtitle, not drawn (controls are marked by their start)
            pp = [q for q in c.periodic if q.symbol == sym]
            isempty(pp) && continue
            CM.scatter!(ax, [q.x for q in pp], [q.p for q in pp]; marker = mk, markersize = sym == :circle ? 6 : 7, color = :transparent, strokecolor = (:black, 0.8), strokewidth = 0.8)
        end
    end
    els = [CM.MarkerElement(color = colors[k], marker = :circle, markersize = 10) for k in CLASS_ORDER]
    labs = ["T torus", "R resonant", "C chaotic sea", "S sticky", "X near separatrix"]
    push!(els, CM.MarkerElement(color = :transparent, strokecolor = :black, strokewidth = 1.6, marker = :rect, markersize = 11)); push!(labs, "control start")
    push!(els, CM.MarkerElement(color = :black, marker = :star5, markersize = 11)); push!(labs, "orbit start (star)")
    push!(els, CM.MarkerElement(color = :transparent, strokecolor = :black, strokewidth = 1, marker = :circle, markersize = 10)); push!(labs, "periodic point, elliptic")
    push!(els, CM.MarkerElement(color = :transparent, strokecolor = :black, strokewidth = 1, marker = :diamond, markersize = 10)); push!(labs, "parabolic")
    CM.Legend(fig[3, 1:3], els, labs; orientation = :horizontal, nbanks = 2, tellheight = true, framevisible = false)
    CM.Label(fig[0, 1:3], "p1-le-benchmark panel: first 1000 iterates (wrapped) of every panel orbit, 870 drawn + 13 controls", fontsize = 18)
    png = joinpath(OUT, "experiments-portraits.png")
    CM.save(png, fig; px_per_unit = 1)
    emit("PORTRAITS", "figure written", relpath(png, ROOT), "-", "INFO")
    mkpath(STUDY_OUT)
    cp(png, joinpath(STUDY_OUT, "p1-le-benchmark-portraits.png"); force = true)
    info("PORTRAITS", "copied to $(relpath(joinpath(STUDY_OUT, "p1-le-benchmark-portraits.png"), ROOT)); portrait content: $(length(pc)) K values, starts per K $(join((length(c.starts) for c in pc), ", ")), periodic markers per K $(join((length(c.periodic) for c in pc), ", "))")
end

# ---------------------------------------------------------------------------------------------------
section("FINAL") do
    say(@sprintf("# summary: PASS %d, FAIL %d, INFO %d, total %.1f s", COUNTS["PASS"], COUNTS["FAIL"], COUNTS["INFO"], time() - T0))
    if COUNTS["FAIL"] == 0 && isempty(ONLY)
        mkpath(STUDY_OUT)
        dst = joinpath(STUDY_OUT, "p1-le-benchmark-panel.toml")
        cp(PANEL, dst; force = true)
        emit("FINAL", "all checks passed: study panel written", relpath(dst, ROOT), "-", "INFO")
    else
        emit("FINAL", "study panel NOT written (a check failed or a partial run)", COUNTS["FAIL"], "0 failures", "INFO")
    end
end
close(LOG)
