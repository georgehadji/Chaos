# Stage 04 review fixes of Amendment D94 (included from experiments-d94.jl): section REV, one probe per fixed finding (C1, H1, H2, M1, M2, L4 to L7, L9),
# and section O26R, the runner side of O26 (the independent oracle's log and its 24 altered copies against the oracle's store).
# Same output form as the other sections: emit(id, quantity, value, bound, verdict).

# the message of the ArgumentError that f throws, "" when it returns
refusal(f) = try; f(); "" catch e; e isa ArgumentError ? e.msg : rethrow(); end
has(msg, parts...) = !isempty(msg) && all(p -> occursin(p, msg), parts)
snapshot(dirs...) = sort!([(relpath(joinpath(d, f), d), filesize(joinpath(d, f)), mtime(joinpath(d, f))) for d in dirs if isdir(d) for f in readdir(d)])

# H2: a campaign through round 1, the log with entry 1, then the ref chunk of round 1 missing, incomplete or with another chunk hash, in copies of the store
function h2_store()
    fx = d94_std(); store = ex_store("h2-camp"); logp = ex_path("h2.log"); write(logp, "")
    ex_campaign(fx; store = store, round_log = logp)
    digest = r -> d94_digest(d94_ref_hashes(store, r))
    d94_write_log(logp, fx, [(round = 1, decision = "continue")], digest)
    plan = d94_plan(fx; log = logp)
    ord = first(sort([c.ordinal for c in plan.chunks if "ref" in c.arms && XP._cell_round(c.cell) == 1]))
    H = XP.HDF5
    set(f, name, v) = H.h5open(h -> (H.delete_attribute(h, name); H.write_attribute(h, name, v)), f, "r+")
    cases = (("missing", f -> rm(f)), ("incomplete", f -> set(f, "complete", false)), ("chunk hash changed", f -> set(f, "chunk_hash", "0"^64)))
    for (name, damage) in cases
        st = ex_store("h2-$(replace(name, ' ' => '-'))"); cp(store, st)
        damage(XP.chunk_path(st, ord))
        before = snapshot(st, HDIR[])
        m = refusal(() -> ex_campaign(fx; store = st, round_log = logp))
        same = before == snapshot(st, HDIR[])
        emit("REV-H2", "run_campaign, ref chunk $ord of round 1 $name: refused naming entry 1 and the file; no file written: $(same)", first(m, 140), "entry 1, chunk-$(lpad(ord, 5, '0')).h5", pf(has(m, "entry 1", "chunk-$(lpad(ord, 5, '0')).h5") && same))
    end
    st = ex_store("h2-intact"); cp(store, st)
    r = refusal(() -> ex_campaign(fx; store = st, round_log = logp))
    emit("REV-H2", "same store intact: the run is accepted (no refusal)", r, "no refusal", pf(isempty(r)))
end

section("REV") do
    HDIR[] = joinpath(V, "fix", "rev"); mkpath(HDIR[])
    ref(; kw...) = ex_arm("ref", "P", ["L1"], ["Float64"], [20]; kw...)
    emptylog() = (p = ex_path("empty.log"); write(p, ""); p)
    # C1: an indicator spec held only by a closure arm: the plan before closure and the dry plan build (before: "selected by no orbit")
    arms = [ref(; release = "each_round"), d86_arm("ind", "D", ["tg"], ["Float64"], [20]; release = "closure")]
    blk = d86_block(arms, [d86_spec("tangent")]); blk["pilot"]["horizon"] = 20
    fx = d94_fix("rev-c1", blk)
    p1 = try d94_plan(fx; log = emptylog()) catch e; sprint(showerror, e) end
    p2 = try d94_plan(fx; check = true) catch e; sprint(showerror, e) end
    emit("REV-C1", "plan before closure with an indicator spec held only by a closure arm builds; dry check plan builds", "$(p1 isa XP.RunPlan ? p1.counts.chunks : p1) / $(p2 isa XP.RunPlan ? p2.rounds.R_plan : p2)", "builds", pf(p1 isa XP.RunPlan && p2 isa XP.RunPlan))
    # H1: a D_B arm without a release key in a block with release keys
    arms = [ref(; release = "each_round"), ex_arm("dd", "D", ["L5"], ["Float64"], [20]; release = "closure"), ex_arm("db", "D_B", ["L4"], ["Float64"], [20])]
    m = refusal(() -> d94_plan(d94_fix("rev-h1", d94_block(arms)); log = emptylog()))
    emit("REV-H1", "D_B arm without release key refused naming the key path", first(m, 110), "arm[3].release", pf(has(m, "arm[3].release", "closure")))
    # M1: a panel with r_max > 1 and a block without release keys; a log given to a block without release keys
    m = refusal(() -> d94_plan(d94_fix("rev-m1a", d94_block([ref()])); log = emptylog()))
    emit("REV-M1", "panel with 3 rounds, block without release keys refused", first(m, 110), "release keys", pf(has(m, "release keys")))
    plain = [Dict{String, Any}(k => v for (k, v) in o if !(k in ("round", "in_D"))) for o in d94_panel_orbits()]
    f1 = d94_fix("rev-m1b", d94_block([ref()]); orbits = plain, r_max = nothing)
    m = refusal(() -> d94_plan(f1; log = emptylog()))
    emit("REV-M1", "one-round panel, block without release keys, round log given: refused", first(m, 110), "no release key", pf(has(m, "no release key")))
    emit("REV-M1", "the same without the log still plans", d94_plan(f1) isa XP.RunPlan, "true", pf(d94_plan(f1) isa XP.RunPlan))
    # M2: an arm without release in a block with release keys: planning, plan text and digest agree (each_round)
    arms = [ref(; release = "each_round"), ex_arm("fr", "P", ["L1"], ["Float64"], [40]),
            ex_arm("dd", "D", ["L5"], ["Float64"], [20]; release = "closure"), ex_arm("db", "D_B", ["L4"], ["Float64"], [20]; release = "closure")]
    pl = d94_plan(d94_fix("rev-m2", d94_block(arms)); log = emptylog())
    ours = [c.hash for c in pl.chunks if any(a -> a in ("ref", "fr"), c.arms) && occursin(r":r1$", c.cell)]
    emit("REV-M2", "ref digest of round 1 equals SHA-256 of the sorted round 1 chunk hashes of the arms ref and fr (free arm counted)", XP._ref_digest(pl.chunks, pl.block, 1)[1:16],
         bytes2hex(sha256(join(sort(ours) .* "\n")))[1:16], pf(XP._ref_digest(pl.chunks, pl.block, 1) == bytes2hex(sha256(join(sort(ours) .* "\n"))) && any(c -> "fr" in c.arms, pl.chunks) && pl.rounds.release["fr"] == "each_round"))
    # H2: the store check of run_campaign
    h2_store()
    # L4 to L7
    K = [0.0, 1.5]
    full = Dict("1" => Dict("T" => "0.5", "R" => "0.5"), "2" => Dict("T" => "0.5", "R" => "0.5", "C" => "0.5", "S" => "0.5", "X" => "0.5"))
    for (name, key) in (("01", "01"), ("space", " 1"))
        t = Dict(key => Dict("T" => "0.2"), "1" => full["1"], "2" => full["2"])
        m = refusal(() -> XP._check_fractions(Dict("d_fraction" => t), K)); m2 = refusal(() -> XP._check_quota(Dict("quota" => Dict{String, Any}(k => Dict{String, Any}(c => 1 for c in keys(r)) for (k, r) in t)), K, 1))
        emit("REV-L4", "duplicate K index key $(repr(key)) refused in d_fraction and in quota", "$(has(m, "twice")) / $(has(m2, "twice"))", "true / true", pf(has(m, "twice") && has(m2, "twice")))
    end
    ok = XP._check_fractions(Dict("d_fraction" => full), K)
    m1 = refusal(() -> XP._check_fractions(Dict("d_fraction" => Dict("1" => full["1"])), K))
    m2 = refusal(() -> XP._check_fractions(Dict("d_fraction" => Dict("1" => full["1"], "2" => Dict("T" => "0.5"))), K))
    emit("REV-L5", "d_fraction: a complete table (K = 0 classes C, S, X absent) read; a missing K index and a missing class refused naming them", "$(length(ok)) / $(first(m1, 50)) / $(first(m2, 50))", "2 / misses K index 2 / misses class R",
         pf(length(ok) == 2 && has(m1, "misses K index 2") && has(m2, "misses class R at K index 2")))
    t = @elapsed (m1 = refusal(() -> XP._ratio("x", "1e-5000000")); m2 = refusal(() -> XP._ratio("x", "1e99999999999999999999")))
    m3 = refusal(() -> XP._check_rounds_keys(Dict("r_max" => 10^9))); m4 = refusal(() -> XP._check_source(Dict("screen_source" => "random", "screen_points" => 2^40)))
    emit("REV-L6", "huge exponent (1e-5000000, 1e99999999999999999999), r_max 10^9 and screen_points 2^40 refused; seconds for the two fractions", "$(has(m1, "exponent")) $(has(m2, "exponent")) $(has(m3, "at most")) $(has(m4, "at most")); $(round(t; digits = 3))",
         "true true true true; < 1", pf(has(m1, "exponent") && has(m2, "exponent") && has(m3, "at most") && has(m4, "at most") && t < 1))
    f7 = d94_fix("rev-l7", d94_block([ref(; release = "each_round")]); orbits = plain, r_max = 3)
    m = refusal(() -> d94_plan(f7; check = true))
    emit("REV-L7", "release keys on a panel with no orbit of a round: named refusal (before: BoundsError or empty collection)", first(m, 110), "no orbit with a round", pf(has(m, "no orbit with a round")))
    lp = ex_path("utf8.log"); write(lp, UInt8[0x72, 0x6f, 0x75, 0x6e, 0x64, 0x3d, 0x31, 0x0a, 0xff, 0xfe, 0x0a])
    m = refusal(() -> XP.read_round_log(lp))
    emit("REV-L9", "round log with invalid UTF-8 refused naming the file", first(m, 90), "UTF-8", pf(has(m, "UTF-8", "utf8.log")))
end

# ---------------------------------------------------------------------------------------------------
# O26R: the oracle's log (entries 1 to 3) and its 24 altered copies, against the oracle's store (a copy of it per run), through the library guard and the entry point
section("O26R") do
    OR = joinpath(V, "oracle", "o26"); SRC = joinpath(V, "o26")
    HDIR[] = joinpath(V, "fix", "o26r"); mkpath(HDIR[])
    fx = (; prereg = joinpath(SRC, "d94-std-prereg.md"), panel = joinpath(SRC, "d94-std-panel.toml"))
    store = ex_store("store"); cp(joinpath(SRC, "camp"), store)
    nfiles = length(readdir(store))
    mf(tag) = ex_path("manifest-$tag.json")
    cli(log, tag) = (io = IOBuffer(); code = XP.cli(["run", fx.prereg, fx.panel, "--store", store, "--manifest", mf(tag), "--round-log", log, "--thread-override", "--allow-dirty",
                                                   "--code-version", EX_CODE, "--capacity", "64", "--threads", "4", "--manifest-reset"]; io = io); (code, String(take!(io))))
    good = joinpath(OR, "round-log.txt")
    res = ex_campaign(fx; store = store, manifest = mf("good"), round_log = good)
    (code, out) = cli(good, "goodcli")
    emit("O26", "oracle's valid log through run_campaign: status / ran / skipped, store files", "$(res.status) / $(length(res.ran)) / $(length(res.skipped)); $(length(readdir(store)))", "complete / 0 / $nfiles; $nfiles",
         pf(res.status == "complete" && isempty(res.ran) && length(res.skipped) == nfiles && length(readdir(store)) == nfiles))
    emit("O26", "oracle's valid log through the entry point (--round-log): exit code", code, "0", pf(code == 0))
    bad = sort(filter(f -> startswith(f, "round-log-bad-"), readdir(OR)))
    nlib = ncli = nname = nclean = 0; undetected = String[]
    for f in bad
        field, n = match(r"round-log-bad-(.*)-(\d)\.txt", f).captures; n = parse(Int, n)
        before = snapshot(store, HDIR[])
        m = refusal(() -> ex_campaign(fx; store = store, manifest = mf("bad"), round_log = joinpath(OR, f)))
        (code, out) = cli(joinpath(OR, f), "badcli")
        after = snapshot(store, HDIR[])
        # a changed count is in no check of the entry itself: it is caught at the next entry through the chain hash; in the last entry it is not caught at all (counts are Analysis's, P81)
        want = startswith(field, "count_") ? n + 1 : n
        nlib += !isempty(m); ncli += code == 1; nname += has(m, "entry $want"); nclean += !isempty(m) && before == after
        isempty(m) && push!(undetected, f)
        say("O26 |   $f | library: $(first(m, 100)) | exit $code | files unchanged $(before == after)")
    end
    emit("O26", "altered copies refused by the library guard / by the entry point (exit 1) / naming the expected entry / store and manifest folder unchanged, of $(length(bad))", "$nlib / $ncli / $nname / $nclean",
         "$(length(bad) - length(undetected)) / $(length(bad) - length(undetected)) / $(length(bad) - length(undetected)) / $(length(bad) - length(undetected))", pf(nlib == length(bad) - length(undetected) && ncli == nlib && nname == nlib && nclean == nlib))
    isempty(undetected) || info("O26", "accepted, not covered by any check of the log (a count of the last entry, which no later entry chains): " * join(undetected, ", "))
end
