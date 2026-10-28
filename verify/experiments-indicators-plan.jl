# Child script of verify/experiments-indicators.jl, criterion D86b [full]: the plan of a preregistration and a panel file, written as plain text so that
# the plan of the code before Amendment D86 (a git worktree of commit ba21bcc, whose src/ is the pre-D86 source) and the plan of the code now can be compared.
# Usage: julia --project=<env whose Chaoslab is the code to test> experiments-indicators-plan.jl <prereg> <panel> <code_version> <defaults file> <out file>
# Output lines: plan_sha256, plan_text_sha256, block_sha256, counts, then one line per work item (id, hash) and per chunk (ordinal, key, hash). Nothing else is read or written.
using Chaoslab, SHA
const EX = Chaoslab.Experiments
prereg, panel, cv, defaults, out = ARGS
plan = EX.plan_run(prereg, panel; code_version = cv, defaults = defaults)
open(out, "w") do io
    println(io, "plan_sha256 ", plan.hash)
    println(io, "plan_text_sha256 ", bytes2hex(sha256(plan.text)))
    println(io, "block_sha256 ", EX.block_hash(plan.block))
    println(io, "counts ", plan.counts.items, " ", plan.counts.chunks)
    for it in plan.items
        println(io, "item\t", it.id, "\t", it.hash)
    end
    for c in plan.chunks
        println(io, "chunk\t", c.ordinal, "\t", c.key, "\t", c.hash)
    end
end
println("PLAN WRITTEN ", plan.hash)
