# Amendment D54: `auto` resolves per method from the record taken on the recording chain. Spec: experiments-spec.md "Amendment D54".
# Layout in data/gpu-benchmark.toml: [<device>."Float64".recorder.<method>] with threads, cpu_ns_per_step, gpu_ns_per_step, gpu_faster, date,
# launch_steps, orbits, steps. Interface: a record NamedTuple with a `recorder` field (method => (threads, gpu_faster, ...)) is read per method;
# a record without the field keeps the D53h form (one record for all methods). `Experiments.gpu_record(root)` returns the field from the file.

ex_mrec(threads, faster) = (threads = threads, cpu_ns_per_step = 5.0, gpu_ns_per_step = 9.0, gpu_faster = faster, date = "2026-10-19")

@testset "D54 resolution reads the per-method record" begin
    R = Chaoslab.FALLBACK_REASONS
    rec = merge(EX_DREC, (recorder = Dict("L5" => ex_mrec(4, true), "L1" => ex_mrec(4, false), "L3" => ex_mrec(8, true)),))   # EX_DREC says faster: ignored
    dev = ExDev(true, true, true)
    res(m; req = :auto, T = Float64, r = rec) = ex_with_device(() -> EX.resolve_run_backend(req, T, m; threads = 4, record = r), dev)
    @test (res("L5").backend, res("L5").reason) == (:gpu, "")
    @test (res("L1").backend, res("L1").reason) == (:cpu, R[7])                       # measured not faster on the recording chain
    x = res("L3")
    @test (x.backend, x.reason) == (:cpu, R[7]) && occursin("thread", x.note)           # record at another thread count
    @test (res("L9").backend, res("L9").reason) == (:cpu, R[7])                       # no record for the method
    @test (res("L6").backend, res("L6").reason) == (:cpu, R[6])
    @test res("L9"; req = :gpu).backend === :gpu && res("L1"; req = :cpu).reason == ""   # explicit requests need no record
    @test res("L5"; T = BigFloat).reason == R[3] && res("L5"; T = Double64, req = :cpu).reason == ""
    # an empty recorder table: every method absent
    @test all(m -> res(m; r = merge(EX_DREC, (recorder = Dict{String, Any}(),))).reason == R[7], ("L1", "L2", "L5", "L9"))
end

@testset "D54 matrix of the resolution with per-method records" begin
    devs = Dict("not loadable" => ExDev(false, false, false), "not functional" => ExDev(true, false, true),
                "not compiling" => ExDev(true, true, false), "device" => ExDev(true, true, true))
    recs = Dict("absent" => nothing, "equal, faster" => (threads = 4, gpu_faster = true), "equal, not faster" => (threads = 4, gpu_faster = false),
                "other threads" => (threads = 8, gpu_faster = true))
    bad = String[]
    for T in (Float64, Double64, BigFloat), req in (:cpu, :gpu, :auto), (dn, dev) in devs, (rn, r) in recs, method in ("L1", "L5", "L6")
        table = r === nothing ? Dict{String, Any}() : Dict{String, Any}(method => r)
        full = (threads = 4, gpu_faster = true, recorder = table)
        got = ex_with_device(() -> EX.resolve_run_backend(req, T, method; threads = 4, record = full), dev)
        want = ex_expect_backend(req, T, method == "L6", dev, r, 4)
        (got.backend, got.reason) == want || push!(bad, "$T $req $dn $rn $method: got $((got.backend, got.reason)), want $want")
    end
    @test isempty(bad)
end

@testset "D54 gpu_record reads the layout of the file" begin
    root = mktempdir()
    mkpath(joinpath(root, "data"))
    open(joinpath(root, "data", "gpu-benchmark.toml"), "w") do io
        TOML.print(io, Dict{String, Any}("FakeGPU" => Dict{String, Any}("Float64" => Dict{String, Any}("threads" => 23, "cpu_ns_per_step" => 3.4,
            "gpu_ns_per_step" => 2.5, "gpu_faster" => true, "date" => "2026-10-19", "recorder" => Dict{String, Any}(
                "L5" => Dict{String, Any}("threads" => 23, "cpu_ns_per_step" => 9.9, "gpu_ns_per_step" => 27.2, "gpu_faster" => false, "date" => "2026-10-19",
                                          "launch_steps" => 106883, "orbits" => 2048, "steps" => 200000))))))
    end
    r = ex_with_device(() -> EX.gpu_record(root), ex_hdev())
    @test r.threads == 23 && r.gpu_faster === true && r.recorder["L5"].gpu_faster === false && r.recorder["L5"].threads == 23
    @test r.recorder["L5"].cpu_ns_per_step == 9.9 && r.recorder["L5"].gpu_ns_per_step == 27.2 && r.recorder["L5"].launch_steps == 106883
    @test !haskey(r.recorder, "L1")
    # a file without the sub-table: an empty recorder table (every method absent)
    open(joinpath(root, "data", "gpu-benchmark.toml"), "w") do io
        TOML.print(io, Dict{String, Any}("FakeGPU" => Dict{String, Any}("Float64" => Dict{String, Any}("threads" => 23, "gpu_faster" => true))))
    end
    @test isempty(ex_with_device(() -> EX.gpu_record(root), ex_hdev()).recorder)
end

@testset "D54 campaign: backend table and chunks follow the per-method record" begin
    ref = ex_small_run("ref")
    R = Chaoslab.FALLBACK_REASONS
    rec = merge(EX_DREC, (recorder = Dict("L5" => ex_mrec(4, true), "L4" => ex_mrec(4, true), "L1" => ex_mrec(4, false)),))
    store = ex_store("d54-campaign")
    res = ex_with_device(() -> ex_campaign(ref.fx; store = store, device_record = rec, hooks = ex_zero_clock()), ex_hdev())
    @test res.status == "complete" && ex_store_canon(store) == ex_store_canon(ref.store)
    man = EX.read_json(res.manifest)
    rows = Dict((r["precision"], r["method"]) => r for r in man["backend_table"])
    @test rows[("Float64", "L5")]["backend"] == "gpu" && rows[("Float64", "L4")]["backend"] == "gpu"
    @test rows[("Float64", "L1")]["backend"] == "cpu" && rows[("Float64", "L1")]["reason"] == R[7]
    @test haskey(man["device_record"], "recorder") && man["device_record"]["recorder"]["L5"]["gpu_faster"] === true
    for c in ref.plan.chunks
        it = ref.plan.items[first(c.items)]
        a = ex_dev_attrs(store, c.ordinal)
        ex_label(it.precision) == "Float64" || continue
        @test a["backend"] == (it.method == "L1" ? "cpu" : "gpu") && a["backend_fallback_reason"] == (it.method == "L1" ? R[7] : "")
    end
end
