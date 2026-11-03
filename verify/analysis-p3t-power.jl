# D95j (V): the power at the achieved sizes (E69 c, E70) gives identical counts at 1, 4 and 23 threads, over a grid of sizes and run counts with non-trivial success counts.
# Uses the module's own readers on the toy campaign inputs (<root>/campaign) and the internal functions p3t_power_hb and p3t_power_hc. Run: julia --project=. -t 23 verify/analysis-p3t-power.jl [<root>]
using Chaoslab, JSON
const AN = Chaoslab.Analysis
const ROOT = length(ARGS) >= 1 ? ARGS[1] : "C:/Users/tesse/AppData/Local/Temp/claude/C--Users-tesse-Documents-Vibe-Coding-Chaos/5df15ee1-47df-4875-85bd-1bd4d7c665a6/scratchpad/d95v"
C(x) = joinpath(ROOT, "campaign", x)
qp = AN._p3_read(C("inputs/p3-abstention-calibration-parameters.toml"), C("inputs/p3-abstention-preregistration.md"), C("p2/p2-parameters.toml"), C("p2/p2-preregistration.md"), C("p2/p2-manifest.json"), nothing, nothing)
t = AN.p3t_read_params(C("inputs/p3-abstention-test-parameters.toml"); prereg = C("inputs/p3-abstention-preregistration.md"), record = C("inputs/p3-abstention-calibration-record.json"), params = C("inputs/p3-abstention-calibration-parameters.toml"))
rec = JSON.parsefile(C("inputs/p3-abstention-calibration-record.json"))
qt0 = merge(qp.q, t)
model = AN.p3t_record_model(rec, qt0)
nbad = 0; ncase = 0; maxsucc = 0; nonzero = 0
for runs in (40, 400, 2000), (nT, nRW) in ((1, 1), (3, 5), (5, 6), (10, 12), (20, 25), (40, 50), (80, 100))
    qt = merge(qt0, (; pp_runs = runs))
    r = [AN.p3t_power_hb(model, nT, nRW, qt, th) for th in (1, 4, 23)]
    global ncase += 1; same = r[1] == r[2] == r[3]; global nbad += !same; global maxsucc = max(maxsucc, r[1].successes); global nonzero += r[1].successes > 0
    println("D95j | HB power runs=$runs nT=$nT nRW=$nRW | successes at 1/4/23 threads: ", join(getproperty.(r, :successes), "/"), " of ", runs, " | identical incl. power and se: ", same, " | ", same ? "PASS" : "FAIL")
end
for runs in (40, 400, 2000), (nf, nn) in ((2, 4), (6, 12), (12, 40), (30, 100), (60, 300))
    qt = merge(qt0, (; pp_runs = runs))
    r = [AN.p3t_power_hc(rec, nf, nn, qt, th) for th in (1, 4, 23)]
    global ncase += 1; same = r[1] == r[2] == r[3]; global nbad += !same; global maxsucc = max(maxsucc, r[1].successes); global nonzero += r[1].successes > 0
    println("D95j | HC power runs=$runs nf=$nf nn=$nn | successes at 1/4/23 threads: ", join(getproperty.(r, :successes), "/"), " of ", runs, " | identical incl. power and se: ", same, " | ", same ? "PASS" : "FAIL")
end
println("D95j | summary: $ncase cases, $nbad differ between thread counts, $nonzero cases with successes > 0, largest success count $maxsucc | ", nbad == 0 && nonzero > 0 ? "PASS" : "FAIL")
