# Stage 04 verification of module experiments, Amendment D96 (indicator specs held only by closure arms), check O27 (library flow and one count step).
# Run from the repo root (default project, the environment of `bin/chaoslab.jl run`):  julia --project -t 4 verify/experiments-d96.jl <mode> <fixture dir> <out dir> [tag]
#   mode lib:   the library flow: run_campaign by rounds with count_round between the runs, the closure run, and a second closure run (resume, nothing to do)
#   mode count: one count_round on the store and round log of <tag> (used by verify/experiments-d96-cli.sh between the CLI runs)
#   mode analyze: Analysis.analyze (study p3-abstention test side, Amendment D97) on the closed store, round log and manifest of <tag>; output folder <out dir>/analysis-<tag> (D97 stage 04, O28 and O29)
# Environment FAULTS=1 (lib mode): chunks failed twice by a mutate hook in the arms ref, ind, ind_pa1, d64, d64_pa1, chk (O29 of D97); default none (the D96 flow).
# The runner is not patched. At D96 count_round (Analysis of D95) did not read the runner's cell form (M1); the scratch patches p1_cells, p3a_missing, p3b_find of the discovery were then
# loaded for count_round, only when <fixture dir>/patches exists. Since D97 the fixture has no patches/ folder and no Analysis function is redefined (D97 stage 04, O28).
# Fixture: <fixture dir> holds the inputs of the discovery d96 (inputs/, p1/, p2/, p3-abstention-test-panel.toml, patches/); read only. Nothing under data/ is touched.
using Chaoslab, SHA
const AN = Chaoslab.Analysis
const XP = Chaoslab.Experiments
const MODE, FIX, OUTD = ARGS[1], ARGS[2], ARGS[3]
const TAG = length(ARGS) >= 4 ? ARGS[4] : "lib"
const IN = joinpath(FIX, "inputs")
const PRE = joinpath(IN, "p3-abstention-preregistration.sealed.md")
const PANEL = joinpath(FIX, "p3-abstention-test-panel.toml")
const STORE, MAN, LOG = joinpath(OUTD, "store-$TAG"), joinpath(OUTD, "manifest-$TAG.json"), joinpath(OUTD, "round-log-$TAG.txt")
for p in (OUTD, STORE, MAN, LOG); occursin("/data/", replace(abspath(p), "\\" => "/")) && error("refused: $p is inside data/"); end
for pf in ("p1_cells", "p3a_missing", "p3b_find"); isfile(joinpath(FIX, "patches", pf * ".jl")) && include(joinpath(FIX, "patches", pf * ".jl")); end
cr_args() = (; prereg = PRE, params = joinpath(IN, "p3-abstention-calibration-parameters.toml"), test_params = joinpath(IN, "p3-abstention-test-parameters.toml"),
             store = STORE, panel = PANEL, round_log = LOG, p2_manifest = joinpath(FIX, "p2", "p2-indicators-manifest.json"),
             p2_prereg = joinpath(FIX, "p2", "p2-indicators-preregistration.md"), p2_params = joinpath(FIX, "p2", "p2-indicators-analysis-parameters.toml"),
             calibration_record = joinpath(IN, "p3-abstention-calibration-record.json"))
ind_keys() = sort!(collect(keys(get(XP.read_json(MAN), "indicators", Dict{String, Any}()))))
function count_step()
    r = AN.count_round(; cr_args()...)
    println("O27 | count_round | action=$(r.action) round=$(r.round) decision=$(r.decision) orbits=$(r.counts.orbits) | - | INFO")
    return r
end
# FAULTS=1: chunks failed twice (attempts 1 and 2), one per arm: ref L1 k3:T:r2, ind tg k2:C:r2, ind_pa1 k3:C:r2, d64 tg k3:C:r1, d64_pa1 k2:C:r1 (the two cells hold orbits of E_test and D, so that orbits leave D, E77), chk L5 Double64 k2:C:r2, chk L5 BigFloat k3:C:r2
const FAULT_CELLS = (("ref", "method=L1 ", "cell=k3:T:r2"), ("ind", "method=tg ", "cell=k2:C:r2"), ("ind_pa1", "method=pa1 ", "cell=k3:C:r2"), ("d64", "method=tg ", "cell=k3:C:r1"),
                     ("d64_pa1", "method=pa1 ", "cell=k2:C:r1"), ("chk", "method=L5 precision=Double64 ", "cell=k2:C:r2"), ("chk", "method=L5 precision=BigFloat", "cell=k3:C:r2"))
function fault_hook(info, orbits)
    key = info.chunk.key
    any(f -> info.chunk.arms == [f[1]] && occursin(f[2], key) && endswith(key, f[3]), FAULT_CELLS) && info.attempt in (1, 2) && (for o in orbits; o.diagnostics["det"][end] = 2 * o.thresholds["det"]; end)
    nothing
end
const HOOKS = get(ENV, "FAULTS", "0") == "1" ? XP.RunHooks(mutate = fault_hook) : XP.RunHooks()
function run_step(label)
    res = XP.run_campaign(PRE, PANEL; store = STORE, manifest = MAN, threads = 4, allow_dirty = true, round_log = LOG, hooks = HOOKS)
    want = any(c["status"] == "failed" for c in XP.read_json(MAN)["chunks"]) ? "complete with failed chunks" : "complete"      # E76 of D97: the status follows the failed entries of the manifest
    println("O27 | run $label | status=$(res.status) ran=$(length(res.ran)) skipped=$(length(res.skipped)) failed=$(length(res.failed)) indicators=$(ind_keys()) | status ", want, " | ", res.status == want ? "PASS" : "FAIL")
    return res
end
if MODE == "count"
    count_step()
elseif MODE == "analyze"
    out = joinpath(OUTD, "analysis-$TAG"); occursin("/data/", replace(abspath(out), "\\" => "/")) && error("refused: $out is inside data/")
    res = AN.count_round(; cr_args()...); println("D97 | count_round on the closed log | action=$(res.action) round=$(res.round) decision=$(res.decision) | action closed | ", string(res.action) == "closed" ? "PASS" : "FAIL")
    res = AN.analyze(; cr_args()..., manifest = MAN, p1_prereg = joinpath(FIX, "p1", "p1-le-benchmark-preregistration.md"), p1_params = joinpath(FIX, "p1", "p1-le-benchmark-parameters.toml"), out, threads = 4)
    println("D97 | analyze $TAG | status=$(res.status) decisions=$(res.decisions) | no refusal | PASS")
elseif MODE == "lib"
    for p in (STORE, MAN, LOG); ispath(p) && error("refused: $p exists (use a fresh out dir)"); end
    mkpath(STORE); write(LOG, "")
    closed = false
    for r in 1:4
        run_step("round $r")
        r == 1 && println("O27 | indicator table before closure | ", ind_keys(), " | empty (R28) | ", isempty(ind_keys()) ? "PASS" : "FAIL")
        count_step().decision == "close" && (global closed = true; break)
    end
    println("O27 | log closed | $closed | true | ", closed ? "PASS" : "FAIL")
    run_step("closure")
    println("O27 | indicator table at closure | ", ind_keys(), " | the specs of the closure arms | INFO")
    run_step("closure again (resume)")
    println("O27 | round log sha256 | ", bytes2hex(sha256(read(LOG))), " | - | INFO")
else
    error("mode must be lib, count or analyze")
end
