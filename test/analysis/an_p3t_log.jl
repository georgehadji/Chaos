# D95d the round log writer and D95e the verification of the log (V1 of D94). Spec: analysis-spec.md, Amendment D95, E62, E63, I38 to I40, O22 (e) (f), O23, criteria D95d, D95e.
# The expected entries come from the independent reference (an_p3t_ref.jl, an_p3t_helpers.jl: p3tref_entry, p3t_log_entries), written from E62 a, b and the design of the campaign.
# API proposed by the tests: count_round(; prereg, params, test_params, store, panel, round_log, p2_manifest, p2_prereg, p2_params, calibration_record, threads, on_read) returning
# (status = :ok, action in :appended :nothing :closed, round, entry, counts, decision, targets); ArgumentError for every refusal. See the test plan.

# one run of the campaign round by round: after round r the instance holds the ref chunks of the rounds up to r; returns the results and the final log text
function p3t_drive(S, I, st, tag; threads = 1, rounds = P3T_ROUNDS)
    text = ""; res = NamedTuple[]
    for r in 1:rounds
        inst = p3t_instance(S, I, "$tag-r$r"; rounds = r, closure = false)
        isempty(text) || write(inst.log, text)             # the first call creates the file
        x = AN.count_round(; p3t_cr_args(I, st, inst)..., threads)
        push!(res, x)
        text = String(read(inst.log))
        x.decision == "close" && break
    end
    return res, text
end

@testset "D95d the log is written round by round: entry text, counts, decisions, chain" begin
    S, I, st = p3t_cached(), p3t_cached_inputs(), p3_cached("d95")
    res, text = p3t_drive(S, I, st, "main")
    ents, decs = p3t_log_entries(S, I)
    @test length(res) == 3 && text == join(ents) && [x.decision for x in res] == decs == ["continue", "continue", "close"]
    @test [x.round for x in res] == [1, 2, 3] && all(x -> x.status === :ok && x.action === :appended, res)
    for r in 1:3
        c = p3t_counts(S, I.qt, r)
        @test Dict(string(k) => v for (k, v) in pairs(res[r].counts)) == c && res[r].entry == ents[r]
    end
    @test res[1].targets == (t_s = 5, w2 = 11) && occursin("prev_sha256=e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855\n", ents[1])
    # 13 lines, canonical integers, the chain of the previous entries
    @test all(e -> count(==('\n'), e) == 13 && !occursin(r"=0\d", e) && !occursin("\n\n", e), ents) && occursin("prev_sha256=$(p3ref_sha256(ents[1]))\n", ents[2]) && occursin("prev_sha256=$(p3ref_sha256(ents[2]))\n", ents[3])
    # the count lines are in ascending byte order of their names (E62 a)
    names = [match(r"^count\.([a-z_0-9]+)=", l).captures[1] for l in split(ents[3], '\n') if startswith(l, "count.")]
    @test names == sort(names) && names == P3T_COUNT_KEYS
    # drop reasons by hand (I46): round 2 has the drift orbit 14 and the excluded orbit 24 inside count.orbits; the U orbit 23 is counted among the orbits only
    @test res[2].counts.orbits == 22 && res[2].counts.drift == 1 && res[2].counts.excluded == 1 && res[2].counts.t_s == 4 && res[2].counts.w2 == 7
    # the digest of round r covers the ref chunks of the rounds up to r only (a failed chunk included)
    @test occursin("ref_digest=$(p3t_digest(S, 2))\n", ents[2]) && p3t_digest(S, 1) != p3t_digest(S, 2) && length(p3t_ref_hashes(S, 3)) == 24
end

@testset "D95d the decision at the boundaries of the targets and at r_max" begin
    S, st = p3t_cached(), p3_cached("d95")
    # (targets, decisions of the rounds): both reached at round 2; one short in w2; one short in t_s; below the targets at r_max; the floors of the table reached at round 2
    for (tt, tw, want) in ((4, 7, ["continue", "close"]), (4, 8, ["continue", "continue", "close"]), (5, 7, ["continue", "continue", "close"]), (6, 12, ["continue", "continue", "close"]), (3, 5, ["continue", "close"]))
        I = p3t_inputs(S, st, "b$tt-$tw"; target_t_s = tt, target_w2 = tw)
        res, text = p3t_drive(S, I, st, "b$tt-$tw")
        ents, decs = p3t_log_entries(S, I)
        @test [x.decision for x in res] == want == decs && text == join(ents)
        @test all(x -> x.targets == (t_s = tt, w2 = tw), res) && occursin("count.target_t_s=$tt\n", text) && occursin("count.target_w2=$tw\n", text)
    end
end

@testset "D95d a repeated call and a closed log append nothing; a closed log tolerates later chunks" begin
    S, I, st = p3t_cached(), p3t_cached_inputs(), p3_cached("d95")
    ents, _ = p3t_log_entries(S, I)
    # repeated call: the log has the entries of the rounds 1 and 2, the store has the ref chunks of the rounds 1 and 2
    inst = p3t_instance(S, I, "rep"; rounds = 2, closure = false)
    write(inst.log, join(ents[1:2]))
    b0 = read(inst.log); fp = p3t_fp_all(I, st, inst)
    x = AN.count_round(; p3t_cr_args(I, st, inst)...)
    @test x.status === :ok && x.action === :nothing && read(inst.log) == b0 && p3t_fp_all(I, st, inst) == fp
    io = IOBuffer(); @test AN.cli(p3t_cli_cr(I, st, inst); io = io) == 0 && occursin("nothing to append", lowercase(String(take!(io)))) && read(inst.log) == b0
    # a closed log: nothing is appended, the log is verified, exit code 0
    full = p3t_instance(S, I, "closed"; closure = false)
    write(full.log, join(ents))
    b1 = read(full.log)
    x = AN.count_round(; p3t_cr_args(I, st, full)...)
    @test x.action === :closed && read(full.log) == b1
    io = IOBuffer(); @test AN.cli(p3t_cli_cr(I, st, full); io = io) == 0 && occursin("closed", lowercase(String(take!(io)))) && read(full.log) == b1
    # a log closed at round 2 (targets 4 and 7) while the store holds the ref chunks of round 3: still closed, nothing appended
    I2 = p3t_inputs(S, st, "closed2"; target_t_s = 4, target_w2 = 7)
    e2, d2 = p3t_log_entries(S, I2)
    inst2 = p3t_instance(S, I2, "closed2"; closure = false)
    write(inst2.log, join(e2)); b2 = read(inst2.log)
    @test d2 == ["continue", "close"] && AN.count_round(; p3t_cr_args(I2, st, inst2)...).action === :closed && read(inst2.log) == b2
end

@testset "D95d refusals leave the log with its old bytes (E62 c, I39)" begin
    S, I, st = p3t_cached(), p3t_cached_inputs(), p3_cached("d95")
    ents, _ = p3t_log_entries(S, I)
    # a ref chunk of round 1 missing: the round is not finished, the missing cell is named, nothing is appended, the file is not created
    inst = p3t_instance(S, I, "miss"; rounds = 1, closure = false)
    rm(p3t_chunk(S, inst; arm = "ref", role = :L5, cell = "k3:r1"))
    @test an_rejects(() -> AN.count_round(; p3t_cr_args(I, st, inst)...), ["k3:r1"]) && !isfile(inst.log)
    # a ref chunk of round 1 incomplete (attribute complete false)
    inst = p3t_instance(S, I, "incomplete"; rounds = 1, closure = false)
    p3t_set_attr(p3t_chunk(S, inst; arm = "ref", role = :L1, cell = "k2:r1"), "complete", false)
    @test an_rejects(() -> AN.count_round(; p3t_cr_args(I, st, inst)...), ["k2:r1", "complete"]) && !isfile(inst.log)
    # ref chunks of a round above r (the guard of the runner was bypassed): refused, naming the round
    inst = p3t_instance(S, I, "above"; rounds = 2, closure = false)
    @test an_rejects(() -> AN.count_round(; p3t_cr_args(I, st, inst)...), ["round 2"]) && !isfile(inst.log)
    inst3 = p3t_instance(S, I, "above3"; rounds = 3, closure = false); write(inst3.log, ents[1]); b = read(inst3.log)
    @test an_rejects(() -> AN.count_round(; p3t_cr_args(I, st, inst3)...), ["round 3"]) && read(inst3.log) == b
    # a closure chunk before closure: refused (the existence is itself a refusal), the log unchanged, no dataset of a closure arm opened
    inst = p3t_instance(S, I, "closure"; rounds = 1, closure = true)
    reads = NamedTuple[]
    @test an_rejects(() -> AN.count_round(; p3t_cr_args(I, st, inst)..., on_read = (kind, path) -> push!(reads, (kind = kind, path = String(path)))), ["closure"]) && !isfile(inst.log)
    clo = Set(basename(f) for f in p3t_files(inst.store) if !("ref" in AX.read_chunk_attrs(f)["arms"]))
    @test length(clo) == 26 && all(r -> !(basename(r.path) in clo && r.kind === :datasets), reads)
    # a failed write leaves the old bytes: a read-only log file (a defect that is not a refusal of the inputs)
    inst = p3t_instance(S, I, "readonly"; rounds = 2, closure = false)
    write(inst.log, ents[1]); chmod(inst.log, 0o444)
    b = read(inst.log)
    failed = try
        AN.count_round(; p3t_cr_args(I, st, inst)...); false
    catch e
        !(e isa UndefVarError || e isa MethodError)
    end
    chmod(inst.log, 0o666)
    @test failed && read(inst.log) == b
    # the store changed between two calls: entry 1 says t_s = 2; the L1 item of orbit 3 now gives lambda(n_s) above 0.2 lambda_ref, so T_s has one orbit less and the verification refuses entry 1
    inst = p3t_instance(S, I, "changed"; rounds = 1, closure = false)
    write(inst.log, ents[1])
    @test AN.count_round(; p3t_cr_args(I, st, inst)...).action === :nothing
    lns2 = copy(S.design.f.lns); lns2[3] = 0.5
    p3t_alter_chunk(S, inst, I, first(c.ordinal for c in S.chunks if c.arm == "ref" && c.role === :L1 && c.round == 1 && 3 in c.ids); design = merge(S.design.f, (; lns = lns2)))
    @test an_rejects(() -> AN.count_round(; p3t_cr_args(I, st, inst)...), ["entry 1", "count.t_s"]) && read(inst.log, String) == ents[1]
end

@testset "D95d the instrumented reads: root attributes and ref datasets only (I39); the bytes do not depend on the thread count or the run" begin
    S, I, st = p3t_cached(), p3t_cached_inputs(), p3_cached("d95")
    ents, _ = p3t_log_entries(S, I)
    inst = p3t_instance(S, I, "reads"; rounds = 2, closure = false)
    write(inst.log, ents[1])
    reads = NamedTuple[]
    x = AN.count_round(; p3t_cr_args(I, st, inst)..., on_read = (kind, path) -> push!(reads, (kind = kind, path = String(path))))
    ref = Set(basename(f) for f in p3t_files(inst.store) if "ref" in AX.read_chunk_attrs(f)["arms"])
    @test x.action === :appended && !isempty(reads) && all(r -> r.kind in (:attrs, :datasets), reads) && all(r -> r.kind === :attrs || basename(r.path) in ref, reads)
    @test any(r -> r.kind === :datasets, reads) && read(inst.log, String) == join(ents[1:2])
    # threads 1 and 4, and a second run: identical bytes
    texts = String[]
    for th in (1, 4, 1)
        i2 = p3t_instance(S, I, "thr$th-$(length(texts))"; rounds = 2, closure = false)
        write(i2.log, ents[1])
        AN.count_round(; p3t_cr_args(I, st, i2)..., threads = th)
        push!(texts, read(i2.log, String))
    end
    @test allequal(texts) && texts[1] == join(ents[1:2])
end

# ---------------------------------------------------------------- D95e: the alterations of a valid log
function p3t_alterations(ents)
    full = join(ents); out = Pair{String, Tuple{String, Int}}[]               # name => (altered text, entry that must be named)
    lines(e) = split(e, '\n')[1:13]
    rebuild(i, ls) = join(vcat(ents[1:(i - 1)], [join(ls, '\n') * "\n"], ents[(i + 1):end]))
    for i in (1, 3)
        L = lines(ents[i])
        for (j, l) in enumerate(L)
            k, v = split(l, '='; limit = 2)
            nv = (startswith(k, "count.") || k == "round") ? string(parse(Int, v) + 1) : k == "decision" ? (v == "close" ? "continue" : "close") : string(v[1] == '0' ? '1' : '0') * v[2:end]
            M = copy(L); M[j] = k * "=" * nv
            push!(out, "entry $i line $k" => (rebuild(i, M), i))
        end
    end
    L = lines(ents[2])
    push!(out, "reorder two lines" => (rebuild(2, vcat(L[1:2], L[4], L[3], L[5:13])), 2))
    push!(out, "reorder count and panel" => (rebuild(2, vcat(L[1:8], L[10], L[9], L[11:13])), 2))
    push!(out, "count name added" => (rebuild(2, vcat(L[1:8], ["count.extra=1"], L[9:13])), 2))
    push!(out, "count name removed" => (rebuild(2, vcat(L[1:4], L[6:13])), 2))
    push!(out, "leading zero" => (rebuild(2, replace.(L, "count.t_s=4" => "count.t_s=04")), 2))
    push!(out, "plus sign" => (rebuild(2, replace.(L, "count.t_s=4" => "count.t_s=+4")), 2))
    push!(out, "blank line" => (rebuild(2, vcat(L[1:6], [""], L[7:13])), 2))
    push!(out, "comment line" => (rebuild(2, vcat(L[1:6], ["# note"], L[7:13])), 2))
    push!(out, "trailing space" => (rebuild(2, vcat(L[1:12], [L[13] * " "])), 2))
    push!(out, "missing final LF" => (full[1:(end - 1)], 3))
    push!(out, "CRLF" => (replace(full, "\n" => "\r\n"), 1))
    push!(out, "entry after close" => (full * ents[3], 4))
    push!(out, "decision continue at r_max" => (rebuild(3, replace.(lines(ents[3]), "decision=close" => "decision=continue")), 3))
    push!(out, "entry 2 repeated as 3" => (ents[1] * ents[2] * ents[2], 3))
    push!(out, "entry 1 missing" => (ents[2] * ents[3], 1))
    return out
end

@testset "D95e every alteration of a valid log is refused by count_round, naming the entry; the valid logs of 1, 2 and 3 entries pass" begin
    S, I, st = p3t_cached(), p3t_cached_inputs(), p3_cached("d95")
    ents, _ = p3t_log_entries(S, I)
    inst = p3t_instance(S, I, "alter"; closure = false)
    # valid logs of 1, 2 and 3 entries pass the verification (the call does nothing or reports the closed log)
    for n in 1:3
        i = p3t_instance(S, I, "valid$n"; rounds = n, closure = false); write(i.log, join(ents[1:n]))
        @test AN.count_round(; p3t_cr_args(I, st, i)...).action === (n == 3 ? :closed : :nothing)
    end
    alts = p3t_alterations(ents)
    @test length(alts) >= 40 && allunique(first.(alts))
    bad = String[]
    for (name, (text, entry)) in alts
        write(inst.log, text)
        b = read(inst.log)
        ok = an_rejects(() -> AN.count_round(; p3t_cr_args(I, st, inst)..., threads = 4), ["entry $entry"]) && read(inst.log) == b
        ok || push!(bad, name)
    end
    @test isempty(bad)
    # a log whose entry 1 carries the hash of another panel, another preregistration or another digest
    for f in (t -> replace(t, "panel_sha256=" * S.panel_hash => "panel_sha256=" * "f"^64; count = 1), t -> replace(t, "prereg_sha256=" * I.prereg_hash => "prereg_sha256=" * "f"^64; count = 1),
              t -> replace(t, "ref_digest=" * p3t_digest(S, 1) => "ref_digest=" * "0"^64; count = 1))
        write(inst.log, f(join(ents))); @test an_rejects(() -> AN.count_round(; p3t_cr_args(I, st, inst)...), ["entry 1"])
    end
end
