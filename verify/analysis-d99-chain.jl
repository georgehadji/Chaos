# Stage 04 verification of module analysis, Amendment D99, criterion D99g (the chain calibrate, count-round, analyze with the default paths of a study layout).
# Run from the repo root (default project):  julia --project -t 4 verify/analysis-d99-chain.jl <stages dir> <stores dir> <out dir> [fixture record]
#   <stages dir>: a study layout (01-design, 02-panel, 03-run, 04-analyze) holding the p3, p2 and P1 inputs at their default names; the 04-analyze output folder is empty.
#   <stores dir>: holds p2/ and p1/, the stores of the p2 and P1 manifests. <out dir>: must not exist; it receives the p3 run store.
#   [fixture record]: optional; the record that `calibrate` writes is compared with it byte for byte.
# Options are given only where the default is outside the layout: the stores (default under data/) and the manifest and round log of `run` (default the study folders of the
# repository). calibrate, count-round and analyze take no output, record, panel, manifest or parameter option. The commands are called in process through the CLI parsers.
# Nothing under data/ is read or written.
using Chaoslab, SHA
const AN, XP = Chaoslab.Analysis, Chaoslab.Experiments
const L, ST, OUTD = abspath(ARGS[1]), abspath(ARGS[2]), abspath(ARGS[3])
const FIXREC = length(ARGS) >= 4 ? ARGS[4] : nothing
for p in (L, ST, OUTD); occursin("/data/", replace(p, "\\" => "/") * "/") && error("refused: $p is inside data/"); end
ispath(OUTD) && error("refused: $OUTD exists (use a fresh out dir)")
lay(sub, name) = normpath(joinpath(L, sub, name))
const PRE, PAR, TP = lay("01-design/output", "p3-abstention-preregistration.md"), lay("01-design/references", "p3-abstention-calibration-parameters.toml"), lay("01-design/references", "p3-abstention-test-parameters.toml")
const OUT4, REC = lay("04-analyze", "output"), lay("04-analyze/output", "p3-abstention-calibration-record.json")
const LOG, MAN, PANEL, S3 = lay("03-run/output", "p3-abstention-round-log.txt"), lay("03-run/output", "p3-abstention-manifest.json"), lay("02-panel/output", "p3-abstention-panel.toml"), joinpath(OUTD, "store")
sh(f) = bytes2hex(sha256(read(f)))
say(what, got, want, ok) = println("D99g | ", what, " | ", got, " | ", want, " | ", ok ? "PASS" : "FAIL")
function cmd(mod, args)
    io = IOBuffer()
    code = mod.cli(args; io)
    txt = String(take!(io))
    print(txt)
    return code, txt
end

isempty(readdir(OUT4)) || error("refused: $OUT4 is not empty")
cal_out = AN._p3_inputs(PRE, PAR, nothing, Dict{Symbol, Any}()).out    # the folder calibrate writes to when no output folder is given
say("default folder of calibrate equals the record folder of count-round and analyze", cal_out, OUT4, normpath(cal_out) == normpath(OUT4))
if FIXREC === nothing
    before01 = sort(readdir(lay("01-design", "output")))
    code, _ = cmd(AN, ["calibrate", PRE, PAR, "--p2-store", joinpath(ST, "p2"), "--p1-store", joinpath(ST, "p1"), "--threads", "4"])
    say("calibrate exit", code, 0, code == 0)
    say("stage 01 output folder unchanged", sort(readdir(lay("01-design", "output"))) == before01, true, sort(readdir(lay("01-design", "output"))) == before01)
else
    # shortcut: the record of the D96 fixture is built by the test helpers and its hash is pinned in the test-side table, so calibrate cannot produce it on these inputs;
    # it is placed in the default folder of calibrate (resolved above), and the chain from there on uses the defaults only
    cp(FIXREC, joinpath(cal_out, basename(REC)))
    println("D99g | record source | fixture record placed in the default folder of calibrate | - | INFO")
end
say("record at the default folder of count-round and analyze", isfile(REC), true, isfile(REC))
println("D99g | files in the stage 04 output folder before the rounds | ", join(sort(readdir(OUT4)), ", "), " | - | INFO")
rec_sha = sh(REC)

mkpath(S3); write(LOG, "")
runargs = ["run", PRE, PANEL, "--store", S3, "--manifest", MAN, "--round-log", LOG, "--threads", "4", "--allow-dirty"]
closed = false
for r in 1:10
    c, _ = cmd(XP, runargs)
    say("run round $r exit", c, 0, c == 0)
    c, txt = cmd(AN, ["count-round", PRE, PAR, TP, "--store", S3])
    say("count-round after round $r exit", c, 0, c == 0)
    occursin(r"\bclose\b", txt) && (global closed = true; break)
end
say("round log closed", closed, true, closed)
c, _ = cmd(XP, runargs)
say("closure run exit", c, 0, c == 0)
c, _ = cmd(AN, ["count-round", PRE, PAR, TP, "--store", S3])
say("count-round on the closed log exit", c, 0, c == 0)

c, _ = cmd(AN, ["analyze", PRE, PAR, TP, "--store", S3])
say("analyze exit", c, 0, c == 0)
res = lay("04-analyze/output", "p3-abstention-test-results.json")
say("analyze outputs at the default folder", isfile(res) && isfile(lay("04-analyze/output", "p3-abstention-tables.md")), true, isfile(res))
fp = isfile(res) ? XP.read_json(res)["fingerprints"]["calibration_record"] : nothing
say("record read by analyze (fingerprint) is the record in the default folder of calibrate", fp === nothing ? "-" : fp[1:12], rec_sha[1:12], fp == rec_sha)
say("record unchanged after the chain", sh(REC)[1:12], rec_sha[1:12], sh(REC) == rec_sha)
println("D99g | round log sha256 | ", sh(LOG), " | - | INFO")
