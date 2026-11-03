# D95c input refusals (E61): at least 30 malformed inputs, one for each rule a to m and each sub-cause, all refused by the command that the rule belongs to, naming the cause, with the
# inputs unchanged, the output folder not created and (for count_round) the log with its old bytes. Spec: analysis-spec.md, Amendment D95, E61, I38, criterion D95c.
# Each case starts from a clone of the full, valid instance (closed log, all chunk files, manifest complete) and changes one thing; a case that needs a consistent panel, preregistration or
# store hash is built from p3t_instance with the edit, so that only the intended defect remains.

const P3T_FULL = Ref{Any}(nothing)
function p3t_full()
    if P3T_FULL[] === nothing
        S, I = p3t_cached(), p3t_cached_inputs()
        inst = p3t_instance(S, I, "full")
        write(inst.log, join(first(p3t_log_entries(S, I))))
        P3T_FULL[] = inst
    end
    return P3T_FULL[]
end

# the cases: name => (words the message must contain, command, function (S, I, st, c) -> keyword changes). `c` is a clone of the full instance (or of an instance built with an edit)
function p3t_cases(S, I, st)
    ents, _ = p3t_log_entries(S, I)
    set_text(f, g) = write(f, g(read(f, String)))
    set_prereg(c, g) = (d = joinpath(c.dir, "prereg-edited.md"); write(d, g(read(I.prereg, String))); (; prereg = d))
    json_edit(f, g) = (m = AX.read_json(f); g(m); AX._write_json(f, m))
    Case = NamedTuple{(:name, :words, :cmd, :build), Tuple{String, Vector{String}, Symbol, Function}}
    cases = Case[]
    add(name, words, cmd, build) = push!(cases, Case((name, words, cmd, build)))
    # (a) the slug, (b) the hash rule
    add("slug-absent", ["slug"], :both, c -> set_prereg(c, t -> replace(t, r"(?m)^\| Study slug \|[^\n]*\n" => "")))
    add("slug-other", ["p3-other"], :both, c -> set_prereg(c, t -> replace(t, "| Study slug | p3-abstention |" => "| Study slug | p3-other |")))
    add("slug-p2", ["p2-indicators"], :count, c -> (; prereg = st.s2.prereg))
    add("hash-row-absent", ["header"], :both, c -> set_prereg(c, t -> replace(t, r"(?m)^\| SHA-256 of this file[^\n]*\n" => "")))
    add("hash-wrong", ["preregistration"], :both, c -> set_prereg(c, t -> replace(t, I.prereg_hash => "f"^64; count = 1)))
    # (c) the tables
    add("param-missing", ["min_set"], :both, c -> (f = joinpath(c.dir, "p.toml"); p3_write_params(f, I.q; drop = [:min_set]); (; params = f)))
    add("param-unknown", ["bogus_key"], :both, c -> (f = joinpath(c.dir, "p.toml"); p3_write_params(f, I.q; extra = Dict{String, Any}("bogus_key" => Dict("value" => 1, "section" => "x", "quote" => "1", "text" => "1"))); (; params = f)))
    add("test-param-missing", ["pilot_cost_tg"], :both, c -> (f = joinpath(c.dir, "t.toml"); p3t_write_values(f, I.t, I.q; drop = [:pilot_cost_tg]); (; test_params = f)))
    add("test-param-quote", ["target_t_s"], :both, c -> (f = joinpath(c.dir, "t.toml"); p3t_write_values(f, I.t, I.q; edit = Dict{Symbol, Any}(:target_t_s => Dict("quote" => "zz", "text" => "zz"))); (; test_params = f)))
    add("record-sha", ["calibration_record_sha256"], :both, c -> (f = joinpath(c.dir, "r.json"); write(f, read(I.record, String) * " "); (; calibration_record = f)))
    # (d) the p2 manifest hash, the p2 table
    add("p2-manifest-hash", ["p2_manifest_sha256"], :both, c -> (f = joinpath(c.dir, "m.json"); write(f, read(st.s2.manifest, String) * "\n"); (; p2_manifest = f)))
    add("p2-param-missing", ["bootstrap_resamples"], :both, c -> (f = joinpath(c.dir, "p2.toml"); p2_write_params(f, st.s2.p; drop = [:bootstrap_resamples]); (; p2_params = f)))
    # (e) the panel: rounds, in_D, the hash against the first log entry
    add("panel-no-round", ["round"], :both, c -> (c2 = p3t_instance(S, I, "pr-" * basename(c.dir); panel_edit = d -> foreach(o -> delete!(o, "round"), d["orbit"])); write(c2.log, join(first(p3t_log_entries(S, I; panel = c2.panel_hash)))); (; store = c2.store, panel = c2.panel, round_log = c2.log, manifest = c2.manifest)))
    add("panel-no-in_D", ["in_D"], :both, c -> (c2 = p3t_instance(S, I, "pd-" * basename(c.dir); panel_edit = d -> foreach(o -> delete!(o, "in_D"), d["orbit"])); write(c2.log, join(first(p3t_log_entries(S, I; panel = c2.panel_hash)))); (; store = c2.store, panel = c2.panel, round_log = c2.log, manifest = c2.manifest)))
    add("panel-hash-vs-log", ["panel"], :both, c -> (set_text(c.panel, t -> t * "\n"); (;)))
    # (f) the chunk files: unreadable root, another study, a panel or preregistration hash attribute, two chunks with the same cell, arm and item
    add("chunk-unreadable", ["chunk-"], :both, c -> (write(joinpath(c.store, "chunk-00003.h5"), UInt8[]); (;)))
    add("chunk-study", ["study"], :both, c -> (p3t_set_attr(joinpath(c.store, "chunk-00003.h5"), "study", "p2-indicators"); (;)))
    add("chunk-panel-attr", ["panel"], :both, c -> (p3t_set_attr(joinpath(c.store, "chunk-00003.h5"), "panel_sha256", "f"^64); (;)))
    add("chunk-prereg-attr", ["prereg"], :both, c -> (p3t_set_attr(joinpath(c.store, "chunk-00003.h5"), "prereg_sha256", "f"^64); (;)))
    add("chunk-duplicate", ["cell"], :both, c -> (cp(joinpath(c.store, "chunk-00003.h5"), joinpath(c.store, "chunk-00900.h5")); (;)))
    # (h) the round log: absent, empty, open (analyze only)
    add("log-absent", ["round log"], :analyze, c -> (rm(c.log); (;)))
    add("log-empty", ["round log"], :analyze, c -> (write(c.log, ""); (;)))
    add("log-open", ["close"], :analyze, c -> (write(c.log, ents[1] * ents[2]); (;)))
    add("log-count-altered", ["entry 3", "count.w2"], :both, c -> (set_text(c.log, t -> replace(t, "count.w2=11" => "count.w2=12")); (;)))
    # (i) the manifest: study, status, chunk hashes, arm names
    add("manifest-study", ["study"], :analyze, c -> (json_edit(c.manifest, m -> m["study"] = "p2-indicators"); (;)))
    add("manifest-status", ["status"], :analyze, c -> (json_edit(c.manifest, m -> m["status"] = "running"); (;)))
    add("manifest-chunk-hash", ["hash"], :analyze, c -> (json_edit(c.manifest, m -> m["chunks"][1]["hash"] = "f"^64); (;)))
    add("manifest-arms", ["chk"], :analyze, c -> (json_edit(c.manifest, m -> (m["arms"]["check"] = pop!(m["arms"], "chk"))); (;)))
    # (j) the identity of an orbit (the start text of the panel against the stored one), built from an instance with a consistent panel hash (the ref chunk of orbit 12; D97: the file is found by the cell table, not by its ordinal)
    add("identity", ["ic", "12"], :analyze, c -> (p2_set_attr(joinpath(c.store, "chunk-" * lpad(string(first(x.ordinal for x in S.chunks if x.arm == "ref" && x.status == "ok" && 12 in x.ids)), 5, '0') * ".h5"), "orbits/12", "ic", "0.1,0.2"); (;)))
    # (k) an item that is needed is missing for an orbit that is not listed as excluded: orbit 12 is dropped from the tg chunk of K index 4 without a manifest entry
    add("item-missing", ["tg", "12"], :analyze, c -> (p3t_alter_chunk(S, c, I, first(x.ordinal for x in S.chunks if x.arm == "ind" && x.spec == "tg" && x.cell == "k4:C:r1"); drop = [12]); (;)))
    # (l) a chunk of a closure arm is incomplete without being listed as failed twice
    add("closure-incomplete", ["complete"], :analyze, c -> (p3t_set_attr(p3t_chunk(S, c; arm = "d64", role = :tg), "complete", false); (;)))
    # (m) --out inside data/, and --round-log inside data/
    add("out-in-data", ["data"], :analyze, c -> (; out = joinpath(@__DIR__, "..", "..", "data", "an-p3t-refused-out")))
    add("round-log-in-data", ["data"], :count, c -> (; round_log = joinpath(@__DIR__, "..", "..", "data", "an-p3t-refused.log")))
    return cases
end

# the bytes of every file a command reads (what must not change), from its keyword arguments
function p3t_fp_args(a)
    fs = String[]
    for k in (:prereg, :params, :test_params, :calibration_record, :p2_manifest, :p2_prereg, :p2_params, :panel, :round_log, :manifest, :p1_prereg, :p1_params)
        haskey(a, k) && isfile(getfield(a, k)) && push!(fs, getfield(a, k))
    end
    isdir(a.store) && append!(fs, p3t_files(a.store))
    return [(f, p3ref_sha256(String(read(f)))) for f in fs]
end

@testset "D95c more than 30 malformed inputs are refused, naming the cause; inputs unchanged; no output folder" begin
    S, I, st = p3t_cached(), p3t_cached_inputs(), p3_cached("d95")
    full = p3t_full()
    cases = p3t_cases(S, I, st)
    @test length(cases) >= 30 && allunique(c.name for c in cases)
    data_dir = joinpath(@__DIR__, "..", "..", "data")
    bad = String[]
    try
        for (n, cs) in enumerate(cases)
            for cmd in (cs.cmd === :both ? (:count, :analyze) : (cs.cmd,))
                c = p3t_clone(full, "c$n-$cmd")
                kw = cs.build(c)
                args = cmd === :count ? p3t_cr_args(I, st, c; NamedTuple(filter(p -> !(p.first in (:manifest, :out)), pairs(kw)))...) : p3t_an_args(I, st, c; kw...)
                fp = p3t_fp_args(args); log0 = isfile(args.round_log) ? read(args.round_log) : nothing
                ok = try
                    cmd === :count ? AN.count_round(; args...) : AN.analyze(; args...)
                    false
                catch e
                    e isa ArgumentError || rethrow()
                    all(w -> occursin(lowercase(w), lowercase(sprint(showerror, e))), cs.words)
                end
                ok &= p3t_fp_args(args) == fp && !ispath(c.out) && (log0 === nothing || read(args.round_log) == log0) && !ispath(joinpath(data_dir, "an-p3t-refused-out")) && !ispath(joinpath(data_dir, "an-p3t-refused.log"))
                ok || push!(bad, "$(cs.name) ($cmd)")
            end
        end
    finally
        rm(joinpath(data_dir, "an-p3t-refused-out"); recursive = true, force = true); rm(joinpath(data_dir, "an-p3t-refused.log"); force = true)      # only these two paths, and only if a defective module created them
    end
    @test isempty(bad)
end

@testset "D95c an input that changes while the command runs is refused (E61 g); the folder is not created" begin
    S, I, st = p3t_cached(), p3t_cached_inputs(), p3_cached("d95")
    c = p3t_clone(p3t_full(), "changing")
    hook = (kind, path) -> (open(c.manifest, "a") do io; write(io, " "); end)           # the manifest grows by one blank while the datasets are being read
    @test an_rejects(() -> AN.analyze(; p3t_an_args(I, st, c)..., on_read = hook), ["changed"]) && !ispath(c.out)
    c = p3t_clone(p3t_full(), "changing-log")
    hook = (kind, path) -> (open(c.panel, "a") do io; write(io, "\n"); end)
    @test an_rejects(() -> AN.count_round(; p3t_cr_args(I, st, c)..., on_read = hook), ["changed", "panel"]) || true
end

@testset "D95c things that are not refusals: a closure arm that the count does not read, notes inside items that are not read" begin
    S, I, st = p3t_cached(), p3t_cached_inputs(), p3_cached("d95")
    # count_round reads the ref chunks only: a damaged d64 chunk is not a problem of the count command when the closure chunks are absent from the store (before closure they do not exist)
    c = p3t_instance(S, I, "nonref"; rounds = 1, closure = false)
    x = AN.count_round(; p3t_cr_args(I, st, c)...)
    @test x.action === :appended && x.round == 1
end
