# Stage 04 verification of the analysis module, Amendment D95 (test side of study p3-abstention), part 1: materialize the toy campaign of the test helpers
# (31 orbits, 3 rounds, 50 chunk files) as plain files, write the round log round by round with count_round, run analyze. Run from the repo root:
#   julia --project=verify -t 4 verify/analysis-p3t.jl [<scratch root, default the d95v folder of the session>]
# Output folder: <root>/campaign/, or <root>/$D95_CAMPAIGN/ when that variable is set (self-describing, see INPUTS.md there). Nothing in data/, src/, test/ is written.
using Printf, Test, SHA
using Chaoslab

const REPO = dirname(@__DIR__)
const TESTDIR = joinpath(REPO, "test", "analysis")
const ROOT = length(ARGS) >= 1 ? ARGS[1] : "C:/Users/tesse/AppData/Local/Temp/claude/C--Users-tesse-Documents-Vibe-Coding-Chaos/5df15ee1-47df-4875-85bd-1bd4d7c665a6/scratchpad/d95v"

function include_defs(f)
    ex = Meta.parseall(read(f, String); filename = f)
    for e in ex.args
        e isa LineNumberNode && continue
        (e isa Expr && e.head == :macrocall && e.args[1] == Symbol("@testset")) && continue
        Base.eval(Main, e)
    end
end
include(joinpath(TESTDIR, "an_helpers.jl"))
for f in ("an_fixture", "an_p2_ref", "an_p2_helpers", "an_p2_study", "an_p2_build", "an_p2_refstore", "an_p2_fixture", "an_p2_dispatch", "an_p2_refusals")
    include_defs(joinpath(TESTDIR, f * ".jl"))
end
include(joinpath(TESTDIR, "an_p3_ref.jl")); include(joinpath(TESTDIR, "an_p3_helpers.jl"))
include(joinpath(TESTDIR, "an_p3t_ref.jl")); include(joinpath(TESTDIR, "an_p3t_helpers.jl"))

const INPUTS_MD = raw"""# Toy campaign of study p3-abstention (Amendment D95 verification)

A synthetic campaign (not the study): 31 orbits, 3 rounds, 50 chunk files written by the real chunk writer. Every path below is relative to this folder.
Commands run from the repository root with the Julia project `.` (CLI) or `verify` (library). Nothing here is a module output except `round-log.txt` and `module-out/`.

## Files

| File | Role |
|------|------|
| inputs/p3-abstention-preregistration.md | preregistration of the toy campaign with the section "Amendment A1 (freeze of the calibration, 2026-11-01)" (A1 values, quotes of the test-side table) |
| inputs/p3-abstention-calibration-parameters.toml | the 71-key parameter table (D93 E47), r_max 3, scaled floors and yields, 100 bootstrap resamples |
| inputs/p3-abstention-test-parameters.toml | the test-side table (E60): pilot costs, targets, theta*, planned powers, each with section, quote, text |
| inputs/p3-abstention-calibration-record.json | the calibration record (D93 output form) the test side reads |
| inputs/p3-abstention-test-panel.toml | the panel of the test run (carries the rounds: `round`, `k_index`, `in_D`, `r_max`, `gate_rounds`) |
| inputs/p3-abstention-test-manifest.json | the manifest of the completed test run (D94 form) |
| store/chunk-NNNNN.h5 | the store, 50 chunk files: arm ref (L1 and L5, Float64, one pair per cell (round, K index), 24 files), arm ind (tg, fq1, bk1, pa1, 20 files), arm d64 (Double64, 4 files), arm chk (L5 Double64 and BigFloat, 2 files). The cell is in the root attribute chunk_key (`cell=k<K index>:r<round>` for ref, `cell=k<K index>` for ind, `cell=D`, `cell=DB`). |
| p2/p2-manifest.json, p2/p2-preregistration.md, p2/p2-parameters.toml | the P2 study inputs (D89) |
| p1/p1-preregistration.md, p1/p1-parameters.toml | the P1 study inputs |
| round-log.txt | the round log written by `count-round`, round by round (3 entries) |
| module-out/ | the two files written by `analyze` (E73) |

Rounds were counted with a staged store: for round r the store held exactly the ref chunk files of the rounds 1 to r (a subset of `store/`, same bytes), because `count-round` refuses a store with closure-arm chunks before closure. The staged stores are not included; the oracle reads `store/` and selects the ref chunks by chunk_key.

## Commands (library form, keywords)

    count_round(; prereg, params, test_params, store, panel, round_log, p2_manifest, p2_prereg, p2_params, calibration_record, threads)
    analyze(; prereg, params, test_params, store, panel, round_log, manifest, p2_manifest, p2_prereg, p2_params, calibration_record, p1_prereg, p1_params, out, threads)

## Commands (CLI form)

    julia --project bin/chaoslab.jl count-round inputs/p3-abstention-preregistration.md inputs/p3-abstention-calibration-parameters.toml inputs/p3-abstention-test-parameters.toml \
        --store <store> --panel inputs/p3-abstention-test-panel.toml --round-log round-log.txt --p2-manifest p2/p2-manifest.json \
        --p2-preregistration p2/p2-preregistration.md --p2-parameters p2/p2-parameters.toml --calibration-record inputs/p3-abstention-calibration-record.json [--threads N]
    julia --project bin/chaoslab.jl analyze inputs/p3-abstention-preregistration.md inputs/p3-abstention-calibration-parameters.toml inputs/p3-abstention-test-parameters.toml \
        --store store --panel inputs/p3-abstention-test-panel.toml --round-log round-log.txt --manifest inputs/p3-abstention-test-manifest.json \
        --p2-manifest p2/p2-manifest.json --p2-preregistration p2/p2-preregistration.md --p2-parameters p2/p2-parameters.toml \
        --calibration-record inputs/p3-abstention-calibration-record.json --p1-preregistration p1/p1-preregistration.md --p1-parameters p1/p1-parameters.toml --out module-out [--threads N]

Run for this folder: `count-round` three times with `--store` the staged store of round 1, 2, 3 (the log is created by the first call), then `analyze` once with `--store store`.
"""
S, I, st = p3t_cached(), p3t_cached_inputs(), p3_cached("d95")
const CAMPNAME = get(ENV, "D95_CAMPAIGN", "campaign")      # a folder other than "campaign" keeps the earlier campaign; its staged stores go to stage-<name>
const STAGE = CAMPNAME == "campaign" ? "stage" : "stage-" * CAMPNAME
camp = joinpath(ROOT, CAMPNAME); rm(camp; recursive = true, force = true)
for d in ("inputs", "p2", "p1", "store"); mkpath(joinpath(camp, d)); end
full = p3t_instance(S, I, "campaign"; rounds = 3, closure = true)
cpx(src, dst) = cp(src, joinpath(camp, dst))
cpx(I.prereg, "inputs/p3-abstention-preregistration.md"); cpx(I.params, "inputs/p3-abstention-calibration-parameters.toml")
cpx(I.tparams, "inputs/p3-abstention-test-parameters.toml"); cpx(I.record, "inputs/p3-abstention-calibration-record.json")
cpx(full.panel, "inputs/p3-abstention-test-panel.toml"); cpx(full.manifest, "inputs/p3-abstention-test-manifest.json")
cpx(st.s2.manifest, "p2/p2-manifest.json"); cpx(st.s2.prereg, "p2/p2-preregistration.md"); cpx(st.s2.params, "p2/p2-parameters.toml")
cpx(st.p1.prereg, "p1/p1-preregistration.md"); cpx(st.p1.params, "p1/p1-parameters.toml")
for f in readdir(full.store; join = true); cp(f, joinpath(camp, "store", basename(f))); end
P(x) = joinpath(camp, x)
inp = (prereg = P("inputs/p3-abstention-preregistration.md"), params = P("inputs/p3-abstention-calibration-parameters.toml"), tparams = P("inputs/p3-abstention-test-parameters.toml"),
       record = P("inputs/p3-abstention-calibration-record.json"), panel = P("inputs/p3-abstention-test-panel.toml"), manifest = P("inputs/p3-abstention-test-manifest.json"),
       p2m = P("p2/p2-manifest.json"), p2p = P("p2/p2-preregistration.md"), p2q = P("p2/p2-parameters.toml"), p1p = P("p1/p1-preregistration.md"), p1q = P("p1/p1-parameters.toml"),
       store = P("store"), log = P("round-log.txt"), out = P("module-out"))
cr(store; kw...) = Chaoslab.Analysis.count_round(; prereg = inp.prereg, params = inp.params, test_params = inp.tparams, store, panel = inp.panel, round_log = inp.log,
    p2_manifest = inp.p2m, p2_prereg = inp.p2p, p2_params = inp.p2q, calibration_record = inp.record, kw...)
for r in 1:3
    # the staged store of round r: the ref chunks of the rounds up to r only (count-round refuses closure chunks before closure)
    sdir = joinpath(ROOT, STAGE, "r$r"); rm(sdir; recursive = true, force = true); mkpath(sdir)
    keep = Set(c.ordinal for c in S.chunks if c.arm == "ref" && c.round <= r)
    for f in readdir(full.store; join = true)
        m = match(r"chunk-(\d{5})\.h5$", f); m !== nothing && parse(Int, m.captures[1]) in keep && cp(f, joinpath(sdir, basename(f)))
    end
    x = cr(sdir)
    println("count_round round $r: action=$(x.action) round=$(x.round) decision=$(x.decision) counts=$(x.counts)")
    x.decision == "close" && break
end
a = Chaoslab.Analysis.analyze(; prereg = inp.prereg, params = inp.params, test_params = inp.tparams, store = inp.store, panel = inp.panel, round_log = inp.log, manifest = inp.manifest,
    p2_manifest = inp.p2m, p2_prereg = inp.p2p, p2_params = inp.p2q, calibration_record = inp.record, p1_prereg = inp.p1p, p1_params = inp.p1q, out = inp.out)
println("analyze returned: ", typeof(a))
println("module-out: ", sort(readdir(inp.out)))
open(P("INPUTS.md"), "w") do io
    print(io, INPUTS_MD)
end
write(P("READY"), "")
println("campaign written: ", camp, "; files in store: ", length(readdir(inp.store)))
