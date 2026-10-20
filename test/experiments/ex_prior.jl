# Amendment D58: the option --prior-record of `run` and `pilot` and the record-file reader. Spec: experiments-spec.md Amendment D58.

const PR_DEV = "Test GPU 1"

function pr_file(tag, text)
    p = ex_path("prior-$tag.toml")
    write(p, text)
    return p
end

const PR_4T = """
["$PR_DEV".Double64]
cpu_ns_per_step = 388.4381
reason = "kernel does not compile for the element type"
date = "2026-10-17"
threads = 4
gpu_faster = false

["$PR_DEV".Float64]
cpu_ns_per_step = 8.254622999999999
gpu_ns_per_step = 2.467902
date = "2026-10-17"
threads = 4
gpu_ns_per_step_wave = 9.136708984375002
gpu_faster = true
"""

@testset "D58 record reader" begin
    r = EX.read_gpu_record(pr_file("4t", PR_4T), PR_DEV)
    @test (r.threads, r.cpu_ns_per_step, r.gpu_ns_per_step, r.gpu_faster, r.date) == (4, 8.254622999999999, 2.467902, true, "2026-10-17")
    @test isempty(r.recorder)
    withrec = PR_4T * """

["$PR_DEV".Float64.recorder.L5]
cpu_ns_per_step = 9.0
gpu_ns_per_step = 27.0
date = "2026-10-19"
threads = 23
gpu_faster = false
"""
    r2 = EX.read_gpu_record(pr_file("rec", withrec), PR_DEV)
    @test r2.recorder["L5"].threads == 23 && r2.recorder["L5"].gpu_ns_per_step == 27.0 && !r2.recorder["L5"].gpu_faster
    missing_file = ex_path("prior-does-not-exist.toml")
    @test ex_rejects(() -> EX.read_gpu_record(missing_file, PR_DEV), ["prior-does-not-exist"])
    @test ex_rejects(() -> EX.read_gpu_record(pr_file("text", "this is not toml = = ="), PR_DEV), ["prior-text"])
    @test ex_rejects(() -> EX.read_gpu_record(pr_file("4t2", PR_4T), "Other Device"), ["Other Device"])
    @test ex_rejects(() -> EX.read_gpu_record(pr_file("nof64", "[\"$PR_DEV\".Double64]\nthreads = 4\n"), PR_DEV), [PR_DEV, "Float64"])
    # the record in force: no file, no record
    @test EX.gpu_record(ex_path("prior-empty-root")) === nothing
end

@testset "D58 option of run and pilot" begin
    f = pr_file("opt", PR_4T)
    kw = EX._cli_common(Dict{String, Any}("--prior-record" => f); device = PR_DEV)
    @test kw[:prior_record].threads == 4 && kw[:prior_record].date == "2026-10-17"
    @test !haskey(EX._cli_common(Dict{String, Any}(); device = PR_DEV), :prior_record)
    pos, opt = EX._cli_parse(["run", "p.md", "q.toml", "--prior-record", f])
    @test pos == ["run", "p.md", "q.toml"] && opt["--prior-record"] == f
    @test occursin("--prior-record", EX.CLI_USAGE)
    # the test environment has no device in force: the option is refused and nothing is written
    fx = ex_fix(:small)
    io = IOBuffer()
    out = ex_path("prior-pilot-refused")
    @test EX.cli(["pilot", fx.prereg, fx.panel, "--out", out, "--prior-record", f]; io = io) == 1
    msg = String(take!(io))
    @test occursin("--prior-record", msg) && occursin("device", msg) && !ispath(out)
    @test EX.cli(["run", fx.prereg, fx.panel, "--store", ex_path("prior-run-refused"), "--manifest", ex_path("prior-run-refused.json"), "--prior-record", f]; io = io) == 1
    @test !ispath(ex_path("prior-run-refused.json"))
end

@testset "D58 manifest carries the prior record" begin
    fx = ex_fix(:small)
    rec = EX.read_gpu_record(pr_file("man", PR_4T), PR_DEV)
    res = ex_campaign(fx; store = ex_store("prior-man"), device_record = EX_REC, prior_record = rec)
    man = EX.read_json(res.manifest)
    @test man["prior_record"]["threads"] == 4 && man["prior_record"]["cpu_ns_per_step"] == 8.254622999999999 && man["prior_record"]["date"] == "2026-10-17"
    @test man["device_record"]["date"] == EX_REC.date
    res0 = ex_campaign(fx; store = ex_store("prior-man0"))
    @test EX.read_json(res0.manifest)["prior_record"] == "absent"
end
