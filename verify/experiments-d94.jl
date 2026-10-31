# Stage 04 verification of module experiments, Amendment D94 (panel builder in rounds, round log).
# Run from the repo root, one section per process or all:  EXV_ONLY="P1BUILD" julia --project=verify -t 23 verify/experiments-d94.jl
# Spec: workspaces/module-build/stages/01-spec/output/experiments-spec.md, "Amendment D94" (D94b [full], D94q, D94r, D94s, D94t V, O23 to O26).
# One printed line per check:  <id> | quantity | value | bound | PASS/FAIL/INFO
# Sections (EXV_ONLY, comma list; default = all but the two helper sections):
#   P1BUILD P1PLAN  (D94b full)   REDUCED FIXA O26  (inputs of the independent oracle, README.txt)   GENERATOR  O24 O25  SRC  THREADS  COST
#   P1PLAN_OLD (plan of the commit before D94, run under a worktree project)   THREADCHILD (child process of THREADS)
# Environment: EXV_V = output folder (default: the scratchpad folder d94v), EXV_LOG = log file (appended), EXV_HDIR = folder of the test helpers' files.
using Printf, Random, SHA, LinearAlgebra, Statistics
using Chaoslab, StaticArrays, StableRNGs
import TOML

# `EX` is the name the test helpers define (they are included below for FIXA, O26, THREADCHILD); this script says XP.
const XP = Chaoslab.Experiments
const ROOT = dirname(@__DIR__)
const OUT = joinpath(ROOT, "workspaces", "module-build", "stages", "04-verify", "output")
const V = get(ENV, "EXV_V", raw"C:\Users\tesse\AppData\Local\Temp\claude\C--Users-tesse-Documents-Vibe-Coding-Chaos\5df15ee1-47df-4875-85bd-1bd4d7c665a6\scratchpad\d94v")
const STUDY = joinpath(ROOT, "workspaces", "study", "stages")
const P1_PREREG = joinpath(STUDY, "01-design", "output", "p1-le-benchmark-preregistration.md")
const P2_PREREG = joinpath(STUDY, "01-design", "output", "p2-indicators-preregistration.md")
const P1_PANEL = joinpath(STUDY, "02-panel", "output", "p1-le-benchmark-panel.toml")
const P1_MANIFEST = joinpath(STUDY, "03-run", "output", "p1-le-benchmark-manifest.json")
const FIXED_CODE = "d94-verify-fixed-code-version"          # the plan hash depends on the code version: both commits plan with this one
const DEFAULT_SECTIONS = ["P1BUILD", "P1BUILD_AT", "P1PLAN", "REDUCED", "FIXA", "GENERATOR", "O24", "O25", "SRC", "THREADS", "O26", "REV", "O26R", "COST", "README"]
const SECS = filter(!isempty, split(get(ENV, "EXV_ONLY", ""), ","))
wants(s) = isempty(SECS) ? s in DEFAULT_SECTIONS : s in SECS
const HDIR = Ref(get(ENV, "EXV_HDIR", joinpath(V, "helper")))
mkpath(V); mkpath(OUT)
const LOG = open(get(ENV, "EXV_LOG", joinpath(OUT, "experiments-d94-verify-log.txt")), "a")

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
    wants(name) || return
    t0 = time()
    try
        f()
    catch e
        emit(name, "ERROR", sprint(showerror, e, catch_backtrace())[1:min(end, 900)], "-", "FAIL")
    end
    say(@sprintf("# section %s took %.1f s", name, time() - t0))
end
sha_file(path) = bytes2hex(open(sha256, path))
tdelta(s) = abs(s - 2pi * round(s / 2pi))
tdist(u, v) = hypot(tdelta(u[1] - v[1]), tdelta(u[2] - v[2]))

say("# experiments D94 verification, julia $(VERSION), threads $(Threads.nthreads()), sections $(isempty(SECS) ? "default" : join(SECS, ","))")

# The test helpers (fixtures A and B, the hand panel and block, the round log writer) serve FIXA, O26 and THREADCHILD; their files go to HDIR[].
if wants("FIXA") || wants("O26") || wants("THREADCHILD") || wants("REV") || wants("O26R")
    for f in ("ex_helpers.jl", "ex_run_helpers.jl", "ex_d86_helpers.jl", "ex_d94_helpers.jl")
        include(joinpath(ROOT, "test", "experiments", f))
    end
    ex_path(name) = joinpath(HDIR[], name)
end

# ---------------------------------------------------------------------------------------------------
# D94b [full]: the real P1 panel rebuilt, and the real P1 plan
const HOST_KEYS = ("git_sha", "git_dirty", "package_version", "julia_version", "stablerngs_version")
strip_host(txt) = [l for l in split(txt, '\n') if !any(k -> startswith(l, k * " = "), HOST_KEYS)]

# lines of a panel text outside the ignored fields: host and time fields of [manifest] (HOST_KEYS) and the given keys of [manifest.prereg]
function strip_fields(txt, prereg_keys)
    sec, keep = "", String[]
    for l in split(txt, '\n')
        startswith(lstrip(l), "[") && (sec = strip(strip(l), ['[', ']']))
        sec == "manifest" && any(k -> startswith(l, k * " = "), HOST_KEYS) && continue
        sec == "manifest.prereg" && any(k -> startswith(lstrip(l), k * " = "), prereg_keys) && continue
        push!(keep, l)
    end
    return keep
end

function p1_compare(out, tag, prereg_keys)
    a, b = read(P1_PANEL), read(out)
    ncr(x) = count(==(UInt8(13)), x)             # carriage returns
    info("D94b", "[$tag] committed panel (working tree) $(length(a)) bytes SHA-256 $(bytes2hex(sha256(a))); rebuilt $(length(b)) bytes SHA-256 $(bytes2hex(sha256(b))); carriage returns committed/rebuilt $(ncr(a)) / $(ncr(b)) ($(strip(read(`git -C $ROOT ls-files --eol $(relpath(P1_PANEL, ROOT))`, String))): the working tree has the CRLF form of the checkout, the blob has LF)")
    lf(x) = replace(String(x), "\r\n" => "\n")
    la, lb = strip_fields(lf(a), prereg_keys), strip_fields(lf(b), prereg_keys)
    nd = length(la) == length(lb) ? count(i -> la[i] != lb[i], eachindex(la)) : max(length(la), length(lb))
    ign = "git_sha, git_dirty, package_version, julia_version, stablerngs_version" * (isempty(prereg_keys) ? "" : "; [manifest.prereg] " * join(prereg_keys, ", "))
    emit("D94b", "[$tag] rebuilt panel against the committed file (CRLF read as LF, as the module reads files), lines differing outside the ignored fields ($ign); lines $(length(la)) / $(length(lb))", nd, "0", pf(nd == 0 && length(la) == length(lb)))
    same = join(la, "\n") == join(lb, "\n")
    emit("D94b", "[$tag] same, the two texts without those lines are equal (byte for byte)", same, "true", pf(same))
    pa, pb = TOML.parsefile(P1_PANEL), TOML.parsefile(out)
    info("D94b", "[$tag] ignored values: " * join(["$k committed $(pa["manifest"][k]) / rebuilt $(pb["manifest"][k])" for k in HOST_KEYS], "; ") *
         "; prereg " * join(["$k $(pa["manifest"]["prereg"][k]) / $(pb["manifest"]["prereg"][k])" for k in ("path", "header_sha256", "amendments")], "; "))
    emit("D94b", "[$tag] manifest key lists and orbit counts equal", sort(collect(keys(pa["manifest"]))) == sort(collect(keys(pb["manifest"]))) && length(pa["orbit"]) == length(pb["orbit"]), "true",
         pf(sort(collect(keys(pa["manifest"]))) == sort(collect(keys(pb["manifest"]))) && length(pa["orbit"]) == length(pb["orbit"])))
end

function p1_build(prereg, out, tag)
    mkpath(dirname(out))
    if get(ENV, "EXV_REUSE", "") == "1" && isfile(out)
        info("D94b", "[$tag] EXV_REUSE=1: the rebuilt panel already in the scratch folder is compared (its build time is in the first lines of this log)")
    else
        foreach(f -> rm(f; force = true), (out, out * ".provenance.toml"))
        t = @elapsed XP.build_panel(prereg; out = out)
        emit("D94b", "[$tag] real P1 panel rebuilt with the code of HEAD, threads $(Threads.nthreads()), wall seconds", t, "completes", "PASS")
    end
end

# the P1 preregistration as it is now (amendments A4 to A6 changed its hash after the panel was built; the identity fields of the manifest differ and are ignored here)
section("P1BUILD") do
    out = joinpath(V, "p1", "p1-rebuilt-panel.toml")
    p1_build(P1_PREREG, out, "current preregistration")
    p1_compare(out, "current preregistration", ("header_sha256", "computed_sha256", "amendments"))
end
# the preregistration of the commit under which the panel was built (header hash 650b3558..., commit 3991155): nothing may differ but the host fields and the path
section("P1BUILD_AT") do
    pre = joinpath(V, "p1", "p1-preregistration-at-3991155.md")
    write(pre, read(`git -C $ROOT show 3991155:workspaces/study/stages/01-design/output/p1-le-benchmark-preregistration.md`))
    out = joinpath(V, "p1", "p1-rebuilt-panel-at-3991155.toml")
    p1_build(pre, out, "preregistration of commit 3991155")
    p1_compare(out, "preregistration of commit 3991155", ("path",))
end

function plan_dump(path; panel = P1_PANEL)
    plan = XP.plan_run(P1_PREREG, panel; code_version = FIXED_CODE)
    open(path, "w") do io
        println(io, "plan_hash ", plan.hash)
        println(io, "plan_text_sha256 ", bytes2hex(sha256(XP.plan_text(plan))))
        println(io, "counts chunks=", length(plan.chunks), " items=", length(plan.items))
        for c in plan.chunks
            println(io, c.ordinal, " ", c.hash, " ", c.key)
        end
    end
    return plan
end
const REC_PY = raw"""
import json, sys
m = json.load(open(sys.argv[1]))
with open(sys.argv[2], "w") as f:
    for c in m["chunks"]:
        f.write("%d %s\n" % (c["ordinal"], c["key"]))
    f.write("counts %d %d\n" % (m["plan_counts"]["chunks"], m["plan_counts"]["items"]))
"""
section("P1PLAN_OLD") do
    mkpath(joinpath(V, "p1"))
    plan = plan_dump(joinpath(V, "p1", "plan-7f70aa7.txt"))
    info("D94b", "plan of the commit before the amendment (Chaoslab loaded from $(dirname(dirname(pathof(Chaoslab)))), HEAD there $(strip(read(`git -C $(dirname(dirname(pathof(Chaoslab)))) rev-parse HEAD`, String)))), code_version fixed to \"$FIXED_CODE\": hash $(plan.hash), $(length(plan.chunks)) chunks")
end
section("P1PLAN") do
    mkpath(joinpath(V, "p1"))
    newp = joinpath(V, "p1", "plan-new.txt")
    plan = plan_dump(newp)
    oldp = joinpath(V, "p1", "plan-7f70aa7.txt")
    emit("D94b", "P1 plan with the code of HEAD (code_version fixed), chunks / items", "$(length(plan.chunks)) / $(length(plan.items))", "1363 / 15569", pf(length(plan.chunks) == 1363 && length(plan.items) == 15569))
    if isfile(oldp)
        lo, ln = readlines(oldp), readlines(newp)
        emit("D94b", "P1 plan hash, code before D94 (commit 7f70aa7, git worktree under the scratch folder, same fixed code_version) / HEAD", "$(split(lo[1])[2][1:16]) / $(split(ln[1])[2][1:16])", "equal", pf(lo[1] == ln[1]))
        emit("D94b", "P1 plan text SHA-256 equal", lo[2] == ln[2], "true", pf(lo[2] == ln[2]))
        nd = length(lo) == length(ln) ? count(i -> lo[i] != ln[i], 4:length(lo)) : -1
        emit("D94b", "ordinal, key and hash of all chunks, lines differing between 7f70aa7 and HEAD", nd, "0 (of $(length(ln) - 3))", pf(nd == 0))
    else
        emit("D94b", "plan of 7f70aa7 not found (run P1PLAN_OLD under the worktree project first)", oldp, "present", "FAIL")
    end
    # the values recorded by the accepted P1 run: ordinals and keys are independent of code_version; plan and chunk hashes are not
    recp = joinpath(V, "p1", "recorded-keys.txt")
    py = joinpath(V, "p1", "rec.py"); write(py, REC_PY)
    run(`python $py $P1_MANIFEST $recp`)
    rec = readlines(recp)
    mine = ["$(c.ordinal) $(c.key)" for c in plan.chunks]
    nk = count(i -> rec[i] != mine[i], 1:min(length(mine), length(rec) - 1))
    emit("D94b", "ordinal and key of the 1363 chunks against the manifest of the accepted P1 run (recorded plan_sha256 81b04cdf..., made with another code_version, so only keys compare), differing", nk, "0", pf(nk == 0 && rec[end] == "counts 1363 15569"))
    info("D94b", "recorded plan and chunk hashes in the repo (p1 manifest) used code_version of the run commit; the pre-D94 values above are computed here with 7f70aa7")
    rb = joinpath(V, "p1", "p1-rebuilt-panel-at-3991155.toml")
    if isfile(rb)
        p2 = plan_dump(joinpath(V, "p1", "plan-rebuilt-panel.txt"); panel = rb)
        lo, ln = readlines(joinpath(V, "p1", "plan-new.txt")), readlines(joinpath(V, "p1", "plan-rebuilt-panel.txt"))
        info("D94b", "plan of the rebuilt panel (prereg of commit 3991155; host fields differ, so the file SHA-256 differs): plan hash $(p2.hash[1:16]) against $(plan.hash[1:16]) of the committed panel; chunk lines (ordinal, hash, key) differing: $(count(i -> lo[i] != ln[i], 4:length(lo)))")
    end
end

# ---------------------------------------------------------------------------------------------------
# Generator check data for O23: the installed StableRNGs
section("GENERATOR") do
    path = joinpath(V, "generator.txt")
    seeds = (0, 1, 20261012, UInt64(2)^63 + 5)
    mask = (UInt64(1) << 52) - 1
    nbad = 0
    open(path, "w") do io
        println(io, "# StableRNGs ", pkgversion(StableRNGs), ", julia ", VERSION, ". Per seed: 250 words rand(StableRNG(seed), UInt64) (decimal) then 250 values rand(StableRNG(seed)) (Float64, 17 significant digits, %.17g).")
        println(io, "# line: <seed> <kind> <index> <value>; kind u64 or f64; each kind starts from a fresh StableRNG(seed)")
        for s in seeds
            r = StableRNG(s); words = [rand(r, UInt64) for _ in 1:250]
            r = StableRNG(s); fl = [rand(r) for _ in 1:250]
            for (i, w) in enumerate(words)
                println(io, s, " u64 ", i, " ", w)
            end
            for (i, f) in enumerate(fl)
                println(io, s, " f64 ", i, " ", @sprintf("%.17g", f))
            end
            nbad += count(i -> fl[i] != Float64(words[i] & mask) * 2.0^-52, 1:250)
        end
    end
    emit("O23", "generator check data written: 4 seeds x (250 words + 250 floats) = 2000 numbers in $(basename(path)); StableRNGs $(pkgversion(StableRNGs))", filesize(path), "written", "PASS")
    emit("O23", "reading of the spec 'uniform from the low 52 bits of a word': float i equals (word i & (2^52 - 1)) * 2^-52 for the same stream, all 1000 pairs, differing", nbad, "0", pf(nbad == 0))
    info("O23", "seed 2^63 + 5 = $(UInt64(2)^63 + 5) is accepted by StableRNG; the module itself refuses seeds above typemax(Int) (impl note, reading 4)")
end

# ---------------------------------------------------------------------------------------------------
# D94t V: the numeric values of the p3 preregistration in the module source
const P3_VALUES = [("seed 20261030", "20261030", true), ("screen points 262144", "262144", true), ("screen side 512", "512", true), ("orbits per round 1470", "1470", true),
                   ("orbit total 7430", "7430", true), ("D total 1309", "1309", true), ("D per round 261", "261", true), ("gate orbits 80", "80", false),
                   ("size 600/320/550", "(600|320|550)", true), ("D parts 165/36", "(165|36)", false), ("D fraction 0.30", "0\\.30?", true), ("D fraction 0.05", "0\\.05", true),
                   ("quota 250", "250", false), ("quota 100", "100", false), ("quota 40", "40", false), ("quota 25", "25", false), ("round count 5", "5", false)]
# hits read by hand and found to be coincidences (file, line) => reason; they are listed and not scored
const REVIEWED = Dict(("run_grid.jl", 59) => "docstring example of a direction text (0.3 + pi/5), not a D fraction")
section("SRC") do
    files = sort([joinpath(d, f) for (d, _, fs) in walkdir(joinpath(ROOT, "src", "Experiments")) for f in fs if endswith(f, ".jl")])
    emit("D94t", "source files searched (src/Experiments, *.jl, comments included)", length(files), ">= 40", pf(length(files) >= 40))
    nstrict = 0
    for (name, pat, strict) in P3_VALUES
        re = Regex("(?<![\\w.])" * pat * "(?![\\w]|\\.\\d)")
        hits, nrev = String[], 0
        for f in files, (ln, l) in enumerate(eachline(f))
            match(re, l) === nothing && continue
            code = first(split(l, '#'; limit = 2))
            why = get(REVIEWED, (basename(f), ln), "")
            nrev += isempty(why) ? 0 : 1
            push!(hits, "$(basename(f)):$ln [$(occursin(re, code) ? "code" : "comment")$(isempty(why) ? "" : ", reviewed: " * why)] $(strip(l)[1:min(end, 110)])")
        end
        nopen = length(hits) - nrev
        strict && (nstrict += nopen)
        emit("D94t", "hits of the p3 value ($name) in src/Experiments" * (strict ? "" : " (a small or common number: hits listed for reading, not scored)") * (nrev > 0 ? "; $nrev reviewed coincidence" : ""), nopen, strict ? "0" : "reported", strict ? pf(nopen == 0) : "INFO")
        foreach(h -> say("D94t |   hit | " * h), hits[1:min(end, 12)])
    end
    emit("D94t", "hits of the distinctive p3 values (seed, sizes, totals, fractions) in the module source", nstrict, "0", pf(nstrict == 0))
end

include("experiments-d94-panels.jl")
include("experiments-d94-screen.jl")
include("experiments-d94-review.jl")
say("# done: $(COUNTS["PASS"]) PASS, $(COUNTS["FAIL"]) FAIL, $(COUNTS["INFO"]) INFO")
