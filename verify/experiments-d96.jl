# Stage 04 verification of module experiments, Amendment D96 (indicator specs held only by closure arms), check O27 (library flow and one count step).
# Run from the repo root (default project, the environment of `bin/chaoslab.jl run`):  julia --project -t 4 verify/experiments-d96.jl <mode> <fixture dir> <out dir> [tag]
#   mode lib:   the library flow: run_campaign by rounds with count_round between the runs, the closure run, and a second closure run (resume, nothing to do)
#   mode count: one count_round on the store and round log of <tag> (used by verify/experiments-d96-cli.sh between the CLI runs)
# The runner is not patched. count_round is the Analysis of D95, which does not yet read the runner's cell form (mismatch M1, Amendment D97);
# the scratch patches p1_cells, p3a_missing, p3b_find of the discovery are loaded for count_round only (Analysis code, process-local; no Experiments function is redefined).
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
for pf in ("p1_cells", "p3a_missing", "p3b_find"); include(joinpath(FIX, "patches", pf * ".jl")); end
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
function run_step(label)
    res = XP.run_campaign(PRE, PANEL; store = STORE, manifest = MAN, threads = 4, allow_dirty = true, round_log = LOG)
    println("O27 | run $label | status=$(res.status) ran=$(length(res.ran)) skipped=$(length(res.skipped)) failed=$(length(res.failed)) indicators=$(ind_keys()) | status complete | ",
            res.status == "complete" ? "PASS" : "FAIL")
    return res
end
if MODE == "count"
    count_step()
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
    error("mode must be lib or count")
end
