# AC19 determinism (I8), AC20 backend. Thread-count invariance at the process level (julia -t 1/4/N) is stage 04; here the
# concurrency cap `ntasks` of build_panel stands in for it (plan, ambiguity A6).

ex_bytes(b) = read(b.path)

@testset "AC19 determinism" begin
    base = ex_build("base")
    h0 = ex_sha256(ex_bytes(base))
    rerun = ex_build("rerun")                                         # same inputs, same source, second run
    @test ex_sha256(ex_bytes(rerun)) == h0
    @test ex_bytes(rerun) == ex_bytes(base)
    for nt in (1, 2, 4)
        b = ex_build("ntasks$nt", Dict{String, Any}(); ntasks = nt)
        @test ex_sha256(ex_bytes(b)) == h0
    end
    for ch in (16, 4096)                                              # grid chunk sizes (64 cells in the fixture)
        b = ex_build("chunk$ch", Dict{String, Any}(); chunk = ch)
        @test ex_sha256(ex_bytes(b)) == h0
    end
    # no wall-clock time, host name, thread count, backend or device, and no screening value in the file
    forbidden = Set(["wall_time", "timestamp", "time", "date", "hostname", "host", "threads", "nthreads", "backend", "device",
                     "lambda", "lambda_scr", "block_exponents", "blocks", "e_med", "screening"])
    @test isdisjoint(Set(lowercase.(ex_all_keys(base.panel))), forbidden)
    @test !occursin(gethostname(), String(ex_bytes(base)))
    @test !occursin(r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}", String(ex_bytes(base)))        # no timestamp
end

@testset "AC20 backend" begin
    h0 = ex_sha256(ex_bytes(ex_build("base")))
    none, present = ex_F.NoDevice(), ExFakeDevice(true, true)
    for (n, (req, dev)) in enumerate(((:cpu, none), (:gpu, present), (:gpu, none), (:auto, present), (:auto, none)))
        out = ex_path("backend$n-panel.toml")
        pre = ex_prereg_for()
        ex_with_device(() -> EX.build_panel(pre; out = out, backend = req), dev)
        @test ex_sha256(read(out)) == h0
        side = TOML.parsefile(out * ".provenance.toml")
        @test side["backend"] == "cpu" && side["backend_fallback_reason"] == ""
    end
end
