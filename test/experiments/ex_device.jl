# D53 device path (stage 03 return). Spec: Amendment D53 R2 to R5; tangent-core Compute, D29, D33, D44. The package suite has no GPU:
# the device is simulated by a type whose backend is the KernelAbstractions CPU backend, so the very kernel that CUDA would run
# executes here and is checked against the CPU reference loop. Interface under test (new, not present before the return):
#   Chaoslab.Experiments.device_records(ctx, items) -> (recs, prov), recs === nothing when the group falls back; prov holds
#       backend, device, reason, gpu_chunks, chunk_size, launch_steps, launches (the tangent-core provenance fields)
#   campaigns and pilots on a present device run the resolved-GPU run types on it and record the backend per chunk.

# A device with a hook before every launch (Tangent._before_launch), a one-wave record value, a free-memory figure and a compile flag.
struct ExHookDev <: ex_F.AbstractDevice
    hook::Function
    ns::Union{Nothing, Float64}
    free::Int
    compiles::Bool
end
ex_F.probe(d::ExHookDev, ::Type) = (loadable = true, functional = true, device = "FakeGPU", compiles = d.compiles, gpu_faster = true)
ex_F.free_memory(d::ExHookDev) = d.free
ex_F.gpu_step_ns(d::ExHookDev, ::Type) = d.ns
ex_Tg.device_backend(::ExHookDev) = ex_Tg.KernelAbstractions.CPU()
ex_Tg._before_launch(d::ExHookDev, ci::Int, li::Int, carried) = d.hook(ci, li, carried)
struct ExOOM <: Exception end
struct ExLost <: Exception end
struct ExCompile <: Exception end
ex_F.is_oom(::ExOOM) = true
ex_F.is_device_lost(::ExLost) = true
ex_F.is_compile_error(::ExCompile) = true

const EX_DREC = (threads = 4, cpu_ns_per_step = 8.25, gpu_ns_per_step = 2.47, gpu_faster = true, date = "2026-10-19")
ex_hdev(hook = (ci, li, c) -> nothing; ns = 2.47, free = 2^30, compiles = true) = ExHookDev(hook, ns, free, compiles)
ex_zero_clock() = EX.RunHooks(clock = ex_clock(info -> 0.0; start = 0.0))

function ex_with_launch(f, L)
    old = ex_Tg._LAUNCH_STEPS[]
    ex_Tg._LAUNCH_STEPS[] = L
    try
        return f()
    finally
        ex_Tg._LAUNCH_STEPS[] = old
    end
end

# everything a chunk file stores for an orbit except wall and thread seconds, compared bit for bit (isequal: NaN equals NaN)
function ex_rec_same(a, b)
    return a.id == b.id && isequal(a.t, b.t) && isequal(a.state, b.state) && isequal(a.winding, b.winding) && isequal(a.lyap, b.lyap) &&
           isequal(a.diagnostics, b.diagnostics) && isequal(a.thresholds, b.thresholds) && isequal(a.first_exceedance, b.first_exceedance) &&
           a.failure == b.failure && isequal(a.attrs, b.attrs) && all(m -> a.cost[m][[1, 2, 5]] == b.cost[m][[1, 2, 5]], keys(a.cost))
end

ex_dev_items(plan, rt) = NamedTuple[plan.items[i] for c in plan.chunks if c.run_type == rt for i in c.items]
ex_dev_types(plan) = [rt for rt in unique(c.run_type for c in plan.chunks) if rt.precision == "Float64" && !EX.cpu_only_method(rt.method)]
ex_dev_attrs(store, ord) = EX.read_chunk(ex_chunk_files(store)[ord]).attrs

@testset "D53 device path: methods with a kernel path" begin
    # the recorder state is isbits exactly for the methods that resolve to the device; L6 (library QR) is the only CPU-only method
    for m in EX_ALL9
        m == "L2" && continue
        k = m == "L1" ? 1 : 2
        red = (EX._estimator(m, 8.0), EX.RunMonitor(Chaoslab.StandardMap(1.0)))
        st = ex_Tg._init_states(red, Float64, Val(2), Val(k))
        @test isbitstype(typeof(st)) == !EX.cpu_only_method(m)
    end
    @test isbitstype(Chaoslab.Lyapunov.ShadowEnd{2, Float64}) && !EX.cpu_only_method("L2")
    @test [m for m in EX_ALL9 if EX.cpu_only_method(m)] == ["L6"]
end

@testset "D53 device path: records equal the CPU loop on the KA CPU backend" begin
    fx = ex_fix(:full)
    plan = ex_plan(fx)
    ctx = EX.make_ctx(plan)
    bad, seen = String[], Set{String}()
    for rt in ex_dev_types(plan)
        items = ex_dev_items(plan, rt)
        cpu = EX.compute_records(ctx, items, 1)
        got = ex_with_device(() -> EX.device_records(ctx, items), ex_hdev())
        push!(seen, rt.method)
        ok = got.recs !== nothing && got.prov.backend === :gpu && got.prov.device == "FakeGPU" && got.prov.reason == "" &&
             got.prov.launches == 1 && got.prov.launch_steps == items[1].horizon && length(got.recs) == length(cpu) && all(ex_rec_same.(got.recs, cpu))
        ok || push!(bad, string(rt))
    end
    @test isempty(bad)
    @test seen == Set(setdiff(EX_ALL9, ["L6"]))         # eight methods have a kernel path, in lifted and wrapped mode, with the Y5 control
    @test any(it -> it.system == "linear", ex_dev_items(plan, first(rt for rt in ex_dev_types(plan) if rt.method == "L5")))
    @test any(rt -> rt.coordinate_mode == "wrapped", ex_dev_types(plan))
end

@testset "D53 device path: relaunch boundaries not aligned with grid steps" begin
    fx = ex_fix(:full)
    plan = ex_plan(fx)
    ctx = EX.make_ctx(plan)
    firsts = Dict{Tuple, Any}()
    for rt in ex_dev_types(plan)
        rt.tau == 1 && get!(firsts, (rt.method, rt.coordinate_mode), rt)
    end
    @test length(firsts) >= 9
    for (L, tau_ok) in ((1, 1), (7, 1), (13, 1))
        bad = String[]
        for rt in values(firsts)
            items = ex_dev_items(plan, rt)
            N = items[1].horizon
            cpu = EX.compute_records(ctx, items, 1)
            hooks = Ref(0)
            got = ex_with_launch(L) do
                ex_with_device(() -> EX.device_records(ctx, items), ex_hdev((ci, li, c) -> (hooks[] += 1; nothing)))
            end
            ok = got.recs !== nothing && all(ex_rec_same.(got.recs, cpu)) && got.prov.launch_steps == L && got.prov.launches == cld(N, L) &&
                 hooks[] == cld(N, L) * got.prov.gpu_chunks
            ok || push!(bad, "L=$L $(rt.method) $(rt.coordinate_mode)")
        end
        @test isempty(bad)
    end
    # a launch length of 28 serves tau 1, 2 and 4 together: some boundaries fall on grid steps of tau 4, most do not
    bad = String[]
    for rt in ex_dev_types(plan)
        rt.method in ("L1", "L5", "L2") || continue
        items = ex_dev_items(plan, rt)
        got = ex_with_launch(28) do
            ex_with_device(() -> EX.device_records(ctx, items), ex_hdev())
        end
        all(ex_rec_same.(got.recs, EX.compute_records(ctx, items, 1))) || push!(bad, string(rt))
    end
    @test isempty(bad)
    # a forced length that is not a multiple of tau is refused as in tangent-core
    rt = first(rt for rt in ex_dev_types(plan) if rt.tau == 2)
    @test_throws ArgumentError ex_with_launch(7) do
        ex_with_device(() -> EX.device_records(ctx, ex_dev_items(plan, rt)), ex_hdev())
    end
end

@testset "D53 device path: chunk size from free device memory" begin
    fx = ex_fix(:full)
    plan = ex_plan(fx)
    ctx = EX.make_ctx(plan)
    rt = first(rt for rt in ex_dev_types(plan) if rt.method == "L5" && rt.tau == 1 && rt.coordinate_mode == "lifted")
    items = ex_dev_items(plan, rt)
    cpu = EX.compute_records(ctx, items, 1)
    big = ex_with_device(() -> EX.device_records(ctx, items), ex_hdev(; free = 2^30))
    nsys = length(unique(it.system for it in items))          # the Y5 control is a linear map: a batch of its own
    @test big.prov.gpu_chunks == nsys && big.prov.chunk_size >= length(items)
    small = ex_with_device(() -> EX.device_records(ctx, items), ex_hdev(; free = 8000))
    cs = small.prov.chunk_size
    @test 1 <= cs < length(items) && small.prov.gpu_chunks == sum(cld(count(it -> it.system == sy, items), cs) for sy in unique(it.system for it in items))
    @test all(ex_rec_same.(small.recs, cpu))
end

@testset "D53 device path: a failing orbit is recorded, not thrown" begin
    fx = ex_fix(:full)
    plan = ex_plan(fx)
    ctx = EX.make_ctx(plan)
    rt = first(rt for rt in ex_dev_types(plan) if rt.method == "L5" && rt.tau == 1 && rt.coordinate_mode == "lifted")
    items = ex_dev_items(plan, rt)
    N = items[1].horizon
    grid = ctx.grids[(N, 1)]
    cpu = EX.compute_records(ctx, items, 1)
    corrupt(at) = (ci, li, c) -> (li == at && ci == 1 && (c.x[1] = typeof(c.x[1])(NaN, NaN)); nothing)
    # (launch 2 starts at step 7: the failure is in a segment that ends on a grid step; launch 16 starts at step 105 and holds no grid
    # step: the failure is in the tail of the launch and belongs to the next grid step, 158)
    for (at, first_bad, step) in ((2, 8, 8), (16, 158, 106))
        got = ex_with_launch(7) do
            ex_with_device(() -> EX.device_records(ctx, items), ex_hdev(corrupt(at)))
        end
        r = got.recs[1]
        i0 = findfirst(==(first_bad), grid)
        lam = r.lyap["L5"]
        @test !isempty(r.failure) && occursin("state", r.failure) && occursin(string(step), r.failure)
        @test all(isnan, lam[i0:end, :]) && all(isfinite, lam[1:(i0 - 1), :]) && isequal(lam[1:(i0 - 1), :], cpu[1].lyap["L5"][1:(i0 - 1), :])
        @test all(isnan, r.diagnostics["det"][i0:end]) && isequal(r.state[1:(i0 - 1), :], cpu[1].state[1:(i0 - 1), :])
        @test all(ex_rec_same.(got.recs[2:end], cpu[2:end]))      # the other orbits of the launch are untouched
    end
end

@testset "D53 device path: fallback reasons of a group" begin
    fx = ex_fix(:full)
    plan = ex_plan(fx)
    ctx = EX.make_ctx(plan)
    rt = first(rt for rt in ex_dev_types(plan) if rt.method == "L5" && rt.tau == 1 && rt.coordinate_mode == "lifted")
    items = ex_dev_items(plan, rt)
    R = Chaoslab.FALLBACK_REASONS
    thrower(e) = (ci, li, c) -> throw(e)
    got(dev) = ex_with_device(() -> (g = EX.device_records(ctx, items); (g, ex_F.DEVICE[])), dev)
    (g, after) = got(ex_hdev(thrower(ExOOM())))
    @test g.recs === nothing && g.prov.backend === :cpu && g.prov.reason == R[4] && g.prov.gpu_chunks == 0
    (g, after) = got(ex_hdev(thrower(ExCompile())))
    @test g.recs === nothing && g.prov.backend === :cpu && g.prov.reason == R[3]
    (g, after) = got(ex_hdev(thrower(ExLost())))
    @test g.recs === nothing && g.prov.backend === :cpu && g.prov.reason == R[9] && after isa ex_F.NoDevice   # the device is dropped
    (g, after) = got(ex_hdev(; ns = 1e9))                                  # even one handoff interval of one wave exceeds the 1 s budget
    @test g.recs === nothing && g.prov.backend === :cpu && g.prov.reason == R[8] && g.prov.gpu_chunks == 0
    # an exception that is not ours is not swallowed
    @test_throws ErrorException got(ex_hdev(thrower(ErrorException("boom"))))
end

@testset "D53 device path: mixed backend table, records and provenance of a campaign" begin
    ref = ex_full_run("a"; threads = 4)
    fx, plan = ref.fx, ref.plan
    store = ex_store("dev-full")
    res = ex_with_device(() -> ex_campaign(fx; store = store, device_record = EX_DREC, hooks = ex_zero_clock()), ex_hdev())
    @test res.status == "complete"
    man = EX.read_json(res.manifest)
    rows = Dict((r["precision"], r["method"]) => r for r in man["backend_table"])
    for m in EX_ALL9
        r = rows[("Float64", m)]
        @test (r["backend"], r["reason"]) == (m == "L6" ? ("cpu", "reducer is CPU-only") : ("gpu", ""))
    end
    @test all(kr -> kr[2]["backend"] == "cpu" && kr[2]["reason"] == "", (kr for kr in rows if kr[1][1] != "Float64"))
    # the datasets of every chunk equal the CPU campaign bit for bit (KA CPU backend, the same body)
    @test ex_store_canon(store) == ex_store_canon(ref.store)
    ngpu, ncpu = 0, 0
    for c in plan.chunks
        a = ex_dev_attrs(store, c.ordinal)
        it = plan.items[first(c.items)]
        p = ex_label(it.precision)
        if p == "Float64" && it.method != "L6"
            ngpu += 1
            @test a["backend"] == "gpu" && a["device"] == "FakeGPU" && a["backend_fallback_reason"] == "" && a["backend_request"] == "auto"
            @test a["gpu_chunks"] >= 1 && a["steps_per_launch"] == it.horizon && a["launches_per_chunk"] == 1 && a["chunk_size"] >= length(c.items)
            @test man["chunks"][c.ordinal]["backend"] == "gpu"
        else
            ncpu += 1
            @test a["backend"] == "cpu" && a["device"] == "cpu" && a["gpu_chunks"] == 0
            @test a["backend_fallback_reason"] == (it.method == "L6" && p == "Float64" ? "reducer is CPU-only" : "")
        end
    end
    @test ngpu > 0 && ncpu > 0
    # O12 on the device: items that ran on the device reproduce within backend_agreement; on the KA CPU backend bit for bit
    chk = EX.rerun_check(fx.prereg, fx.panel, store; code_version = EX_CODE, threads = 4, thread_override = true, allow_dirty = true, sample = 4)
    @test !isempty(chk) && any(r -> r.backend == "gpu", chk) && all(r -> r.bitwise && r.max_abs_diff == 0.0, chk)
end

@testset "D53 device path: relaunch inside a campaign" begin
    ref = ex_small_run("ref")
    store = ex_store("dev-relaunch")
    res = ex_with_launch(14) do
        ex_with_device(() -> ex_campaign(ref.fx; store = store, device_record = EX_DREC, hooks = ex_zero_clock()), ex_hdev())
    end
    @test res.status == "complete" && ex_store_canon(store) == ex_store_canon(ref.store)
    for c in ref.plan.chunks
        it = ref.plan.items[first(c.items)]
        a = ex_dev_attrs(store, c.ordinal)
        if ex_label(it.precision) == "Float64"
            @test a["backend"] == "gpu" && a["steps_per_launch"] == 14 && a["launches_per_chunk"] == cld(it.horizon, 14)
        else
            @test a["backend"] == "cpu"
        end
    end
end

@testset "D53 device path: fallback reasons per pair in a campaign" begin
    ref = ex_small_run("ref")
    R = Chaoslab.FALLBACK_REASONS
    function run_with(dev, tag)
        store = ex_store("dev-fb-$tag")
        res = ex_with_device(() -> ex_campaign(ref.fx; store = store, device_record = EX_DREC, hooks = ex_zero_clock()), dev)
        @test res.status == "complete" && ex_store_canon(store) == ex_store_canon(ref.store)         # the CPU rerun gives the same records
        return store, EX.read_json(res.manifest)
    end
    f64(store) = [ex_dev_attrs(store, c.ordinal) for c in ref.plan.chunks if ex_label(ref.plan.items[first(c.items)].precision) == "Float64"]
    # out of memory in the first launch only: that group runs on the CPU with the fourth reason, later groups return to the device
    fired = Ref(false)
    store, man = run_with(ex_hdev((ci, li, c) -> (fired[] || (fired[] = true; throw(ExOOM())); nothing)), "oom")
    a = f64(store)
    @test any(x -> x["backend"] == "cpu" && x["backend_fallback_reason"] == R[4], a) && any(x -> x["backend"] == "gpu", a)
    @test all(x -> x["backend"] == "gpu" || x["backend_fallback_reason"] == R[4], a)
    # a lost device: the first group falls back with the ninth reason and every later group stays on the CPU
    store, man = run_with(ex_hdev((ci, li, c) -> throw(ExLost())), "lost")
    @test all(x -> x["backend"] == "cpu" && x["backend_fallback_reason"] == R[9], f64(store))
    # a kernel that does not compile, and a launch over the display budget
    store, man = run_with(ex_hdev((ci, li, c) -> throw(ExCompile())), "compile")
    @test all(x -> x["backend"] == "cpu" && x["backend_fallback_reason"] == R[3], f64(store))
    store, man = run_with(ex_hdev(; ns = 1e9), "watchdog")
    @test all(x -> x["backend"] == "cpu" && x["backend_fallback_reason"] == R[8], f64(store))
end

@testset "D53 device path: K19 record thread count against the run thread count" begin
    ref = ex_small_run("ref")
    R = Chaoslab.FALLBACK_REASONS
    function table(rec)
        store = ex_store("dev-k19")
        res = ex_with_device(() -> ex_campaign(ref.fx; store = store, device_record = rec, hooks = ex_zero_clock()), ex_hdev())
        man = EX.read_json(res.manifest)
        return store, man, Dict((r["precision"], r["method"]) => r for r in man["backend_table"])
    end
    # taken at 8 threads, run at 4: the device is not used, the existing seventh reason, the mismatch noted
    store, man, rows = table(merge(EX_DREC, (threads = 8,)))
    @test rows[("Float64", "L5")]["backend"] == "cpu" && rows[("Float64", "L5")]["reason"] == R[7] && occursin("thread", rows[("Float64", "L5")]["note"])
    @test all(c -> c["backend"] == "cpu", man["chunks"]) && ex_store_canon(store) == ex_store_canon(ref.store)
    # not faster at the run's thread count: the same reason
    store, man, rows = table(merge(EX_DREC, (gpu_faster = false,)))
    @test rows[("Float64", "L5")]["backend"] == "cpu" && rows[("Float64", "L5")]["reason"] == R[7]
    # record at the thread count of the run: the device
    store, man, rows = table(EX_DREC)
    @test rows[("Float64", "L5")]["backend"] == "gpu" && rows[("Float64", "L5")]["reason"] == ""
    # an explicit gpu request needs no record (the K19 gate is for auto)
    gp = ex_fix_with(ex_rb_with(ex_rb(:small), ("backend", "Float64") => "gpu"), "dev-gpu-request")
    store = ex_store("dev-gpu-request")
    res = ex_with_device(() -> ex_campaign(gp; store = store, device_record = merge(EX_DREC, (threads = 8, gpu_faster = false)), hooks = ex_zero_clock()), ex_hdev())
    @test res.status == "complete" && ex_dev_attrs(store, 1)["backend"] == "gpu" && ex_dev_attrs(store, 1)["backend_request"] == "gpu"
end

@testset "D53 device path: failed chunk reruns on the CPU, cost cap uses the device line" begin
    fx = ex_fix(:small)
    ref = ex_small_run("ref")
    R = Chaoslab.FALLBACK_REASONS
    nan!(orbits, info) = (first(values(orbits[1].lyap))[2, 1] = NaN)
    log = NamedTuple[]
    store = ex_store("dev-rerun")
    hooks = EX.RunHooks(clock = ex_clock(info -> 0.0; start = 0.0), mutate = ex_fault(nan!; ord = 2, log = log))
    res = ex_with_device(() -> ex_campaign(fx; store = store, device_record = EX_DREC, hooks = hooks), ex_hdev())
    @test res.status == "complete"
    @test [(l.attempt, l.backend) for l in log if l.ordinal == 2] == [(1, :gpu), (2, :cpu)]
    a = ex_dev_attrs(store, 2)
    @test a["status"] == "ok" && a["attempt"] == 2 && a["backend"] == "cpu" && a["backend_fallback_reason"] == R[5] && !isempty(a["first_attempt_reasons"])
    # first-chunk cost cap: 100 ns per orbit-step is under the CPU line (100 x 2 x factor 2) and over the device line (10 x 2 x 2)
    slow = EX.RunHooks(clock = ex_clock(info -> info.steps * 100e-9; start = 0.0))
    cpu_run = ex_campaign(fx; store = ex_store("dev-cap-cpu"), hooks = slow)
    @test cpu_run.status == "complete" && isempty(cpu_run.paused)
    gpu_run = ex_with_device(() -> ex_campaign(fx; store = ex_store("dev-cap-gpu"), device_record = EX_DREC, hooks = slow), ex_hdev())
    @test gpu_run.status == "paused arms" && "main" in gpu_run.paused && "h4d" ∉ gpu_run.paused
    @test any(r -> occursin("budget line 20.0 ns", r), gpu_run.reasons)           # ns_gpu 10 x safety 2 (the CPU line is 200)
end

@testset "D53 device path: pilot on the device" begin
    fx = ex_fix(:small)
    out = ex_path("dev-pilot"); rm(out; recursive = true, force = true)
    res = ex_with_device(() -> EX.run_pilot(fx.prereg, fx.panel; out = out, threads = 4, thread_override = true, allow_dirty = true, code_version = EX_CODE,
        capacity = 64, device_record = EX_DREC, hooks = EX.RunHooks(clock = ex_clock(info -> info.steps * 150e-9; start = 0.0))), ex_hdev())
    cost = res.report["cost"]
    @test !isempty(cost) && all(c -> c["backend"] == (c["precision"] == "Float64" ? "gpu" : "cpu"), cost)
    @test all(c -> isapprox(c["thread_seconds_per_orbit_step"], 150e-9; rtol = 1e-9), cost)      # one worker: the device wall time
    @test all(f -> (a = EX.read_chunk(f).attrs; a["pilot"] === true && a["backend"] in ("gpu", "cpu")), ex_chunk_files(joinpath(out, "store")))
    @test any(f -> EX.read_chunk(f).attrs["backend"] == "gpu", ex_chunk_files(joinpath(out, "store")))
end
