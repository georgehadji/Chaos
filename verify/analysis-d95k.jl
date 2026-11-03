# D95k (V): the HC power at achieved failures (E70) against an independent simulation. Writes the module's success counts at a grid of
# (n_f, n_n) to <out>/module.jsonl; verify/analysis-d95k.py then compares them with sim_hc of verify/analysis-p3-oracle.py within 3 SE.
# Run: julia --project=. -t 4 verify/analysis-d95k.jl <campaign folder> <out folder>
using Chaoslab, JSON
const AN = Chaoslab.Analysis
const CAMP, OUT = ARGS[1], ARGS[2]
C(x) = joinpath(CAMP, x)
qp = AN._p3_read(C("inputs/p3-abstention-calibration-parameters.toml"), C("inputs/p3-abstention-preregistration.md"), C("p2/p2-parameters.toml"), C("p2/p2-preregistration.md"), C("p2/p2-manifest.json"), nothing, nothing)
t = AN.p3t_read_params(C("inputs/p3-abstention-test-parameters.toml"); prereg = C("inputs/p3-abstention-preregistration.md"), record = C("inputs/p3-abstention-calibration-record.json"), params = C("inputs/p3-abstention-calibration-parameters.toml"))
rec = JSON.parsefile(C("inputs/p3-abstention-calibration-record.json"))
qt = merge(qp.q, t, (; pp_runs = 4000))
hc = rec["pp_hc"]
open(joinpath(OUT, "module.jsonl"), "w") do io
    for (nf, nn) in ((2, 4), (3, 6), (4, 8), (5, 10), (6, 12), (8, 16), (10, 30), (20, 60))
        r = AN.p3t_power_hc(rec, nf, nn, qt, Threads.nthreads())
        JSON.print(io, Dict("nf" => nf, "nn" => nn, "successes" => r.successes, "runs" => r.runs, "auc_i_used" => hc["auc_i_used"], "auc_p" => hc["auc_p"],
                            "delta_c" => string(qt.delta_c), "z_one_sided" => string(qt.z_one_sided)))
        println(io)
        println("module nf=$nf nn=$nn successes $(r.successes) of $(r.runs)")
    end
end
