# Stage 04 verification of the analysis module, Amendment D93 (study p3-abstention). Run from the repo root:
#   julia --project=verify -t 23 verify/analysis-p3.jl
# Spec: workspaces/module-build/stages/01-spec/output/analysis-spec.md, "Amendment D93" (criteria D93a to D93r; oracles O15 to O19).
# Section SUITE (stage 03, decision D38): the D93 files of the package suite at full size. The package suite (`Pkg.test()`) carries a representative subset of the 71-key malformed parameter
# files, of the 52 refusal cases of D93c and of the five missing-data fixtures of C0; this script sets P3_FULL_RUN = 1 and runs the same test files with every case and the same
# assertions and tolerances. The (V) parts of the criteria (the real stores, the independent layer O15, O19, the wall time and the memory at 23 threads) are added at stage 04.
# One printed line per check:  <id> | quantity | value | bound | PASS/FAIL
# Nothing in data/, src/, the studies or the specs is written; the test fixtures are written under a temporary folder.
ENV["P3_FULL_RUN"] = "1"
using Printf, Test
using Chaoslab

const TESTDIR = joinpath(dirname(@__DIR__), "test", "analysis")

# The definitions of a test file (its helpers), without its testsets: the P2 helper files mix both.
function include_defs(f)
    ex = Meta.parseall(read(f, String); filename = f)
    for e in ex.args
        e isa LineNumberNode && continue
        (e isa Expr && e.head == :macrocall && e.args[1] == Symbol("@testset")) && continue
        Base.eval(Main, e)
    end
end

emit(id, what, value, bound, status) = println(id, " | ", what, " | ", value, " | ", bound, " | ", status)

include(joinpath(TESTDIR, "an_helpers.jl"))
for f in ("an_fixture", "an_p2_ref", "an_p2_helpers", "an_p2_study", "an_p2_build", "an_p2_refstore", "an_p2_fixture", "an_p2_dispatch", "an_p2_refusals")
    include_defs(joinpath(TESTDIR, f * ".jl"))
end
include(joinpath(TESTDIR, "an_p3_ref.jl"))
include(joinpath(TESTDIR, "an_p3_helpers.jl"))

FILES = ["an_p3_fixture", "an_p3_rule", "an_p3_sets", "an_p3_select", "an_p3_pp", "an_p3_params", "an_p3_dispatch", "an_p3_refusals", "an_p3_pipeline", "an_p3_gate"]
for f in FILES
    local ts
    t = @elapsed try
        ts = @testset verbose = true "$f (full size)" begin
            include(joinpath(TESTDIR, f * ".jl"))
        end
    catch e
        e isa Test.TestSetException || rethrow()
        ts = nothing
        println("# $f: TestSetException: $(e.pass) passed, $(e.fail) failed, $(e.error) errors")
    end
    c = ts === nothing ? nothing : Test.get_test_counts(ts)
    npass = c === nothing ? -1 : c.passes + c.cumulative_passes
    nbad = c === nothing ? -1 : c.fails + c.errors + c.cumulative_fails + c.cumulative_errors
    emit("SUITE", "$f at full size (P3_FULL_RUN), assertions passed (fail + error = $nbad); wall time $(round(t; digits = 1)) s", npass, "0 fail, 0 error", c !== nothing && nbad == 0 && npass > 0 ? "PASS" : "FAIL")
end
