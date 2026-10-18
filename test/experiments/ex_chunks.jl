# D53e series content, D53f direct call (O9), D53g result independence, D53i thread independence, D53s Y5 and exact exponents (O11),
# D53t CPU reference rerun (O12), D53v time and memory, the deviation frame (R2 item 6). Full fixture: all nine methods at Float64,
# Double64 on 2 orbits, BigFloat on 1; horizons 100 and 200. Three campaigns (threads, capacity and order differ) are shared.

const EX_RUNS = Dict{Any, Any}()
function ex_full_run(tag::AbstractString; threads::Integer, capacity = 10^6, order::Symbol = :plan)
    return get!(EX_RUNS, tag) do
        fx = ex_fix(:full)
        store = ex_store("full-$tag")
        res = ex_campaign(fx; store = store, threads = threads, capacity = capacity, order = order)
        (; fx, store, res, plan = ex_plan(fx))
    end
end
ex_full_runs() = (ex_full_run("a"; threads = 4), ex_full_run("b"; threads = 1, capacity = 1, order = :reverse),
                  ex_full_run("c"; threads = 2, capacity = 7))

# the plan item of an orbit record in a chunk
ex_item_of(plan, c, o) = plan.items[c.items[findfirst(i -> plan.items[i].orbit == o.id, c.items)]]
function ex_each_orbit(f, run)
    for file in ex_chunk_files(run.store)
        c = run.plan.chunks[ex_ordinal(file)]
        ch = EX.read_chunk(file)
        for o in ch.orbits
            f(ex_item_of(run.plan, c, o), o, ch, c)
        end
    end
end

# diagnostics present per method (R2 item 2): det and symplecticity for L1, L3 to L6, L9; orthogonality for the QR family;
# sum and pairing for the full-spectrum methods (k = d = 2); none for L2
function ex_diag_names(m)
    m == "L2" && return String[]
    m == "L1" && return ["det", "symp"]
    m == "L9" && return ["det", "symp", "sum", "pair"]
    return ["det", "symp", "orth", "sum", "pair"]
end

@testset "D53 deviation frame (R2 item 6)" begin
    for T in (Float64, Double64, BigFloat), text in ("pi/4", "0.3", "0.3 + pi/5")
        setprecision(BigFloat, 256) do
            W = EX.deviation_frame(T, text, 2)
            @test size(W) == (2, 2) && eltype(W) === T
            a = setprecision(() -> ex_dir_big(text), BigFloat, 1024)
            c, s = setprecision(() -> (cos(a), sin(a)), BigFloat, 1024)
            tol = 4 * Chaoslab.effective_eps(T)
            @test abs(BigFloat(W[1, 1]) - c) <= tol && abs(BigFloat(W[2, 1]) - s) <= tol      # first vector: the direction
            @test W[1, 2] == -W[2, 1] && W[2, 2] == W[1, 1]                                   # second: counterclockwise quarter turn, exact
            W1 = EX.deviation_frame(T, text, 1)
            @test size(W1) == (2, 1) && W1[:, 1] == W[:, 1]
        end
    end
end

@testset "D53e series content" begin
    run = ex_full_run("a"; threads = 4)
    @test run.res.status == "complete"
    files = ex_chunk_files(run.store)
    @test length(files) == run.plan.counts.chunks
    grids = Dict{Tuple{Int, Int}, Vector{Int}}()
    bad = Dict{String, Int}()
    miss!(tag, ok) = ok || (bad[tag] = get(bad, tag, 0) + 1)
    n_items = Ref(0)
    ex_each_orbit(run) do it, o, ch, c
        n_items[] += 1
        n = length(o.t)
        g = get!(grids, (it.horizon, it.tau)) do
            ex_ref_grid(it.horizon, it.tau; ppd = 5, first_steps = 10, fracs = [0.1, 0.5])
        end
        miss!("grid", o.t == g)
        miss!("dtype", o.t isa Vector{Int64} && o.state isa Matrix{Float64} && all(v -> v isa Matrix{Float64}, values(o.lyap)) &&
                       all(v -> v isa Vector{Float64}, values(o.diagnostics)))
        miss!("lyap", Set(keys(o.lyap)) == Set([it.method]) && size(o.lyap[it.method]) == (n, ex_vectors(it.method)))
        miss!("diag names", Set(keys(o.diagnostics)) == Set(ex_diag_names(it.method)))
        miss!("diag monotone", all(v -> issorted(v), values(o.diagnostics)) && all(v -> length(v) == n, values(o.diagnostics)))
        miss!("threshold", all(k -> isfinite(o.thresholds[k]) && o.thresholds[k] > 0 && o.first_exceedance[k] >= 0, keys(o.diagnostics)))
        miss!("finite", all(isfinite, o.state) && all(v -> all(isfinite, v), values(o.lyap)) && all(v -> all(isfinite, v), values(o.diagnostics)))
        miss!("cost", haskey(o.cost, it.method) && length(o.cost[it.method]) == 5 && o.cost[it.method][1] == it.horizon &&
                      o.cost[it.method][2] == ex_vectors(it.method))
        miss!("attrs", o.attrs["item_hash"] == it.hash && o.attrs["steps_recorded"] == n)
        ref = ex_iter_states(it, o.t)
        if it.coordinate_mode == "wrapped"
            miss!("wrapped state", o.winding isa Matrix{Int64} && size(o.winding) == (n, 2) &&
                                   all(i -> o.state[i, :] == ref[i][1:2] && o.winding[i, :] == Int.(ref[i][3:4]), 1:n))
            # wrapped point plus winding counts reconstructs the lifted point to 1 ulp of the lifted value
            ok = all(1:n) do i
                lifted = Chaoslab.lift(SVector(ref[i]...))
                rec = o.state[i, :] .+ (2 * pi) .* Float64.(o.winding[i, :])
                all(j -> abs(rec[j] - lifted[j]) <= eps(abs(lifted[j])), 1:2)
            end
            miss!("wrapped lift", ok)
        else
            miss!("lifted state", o.winding === nothing && all(i -> o.state[i, :] == ref[i][1:2], 1:n))     # 0 ulp
        end
    end
    @test n_items[] == run.plan.counts.items
    @test isempty(bad)
end

@testset "D53f direct call (O9)" begin
    runs = ex_full_runs()
    for (ri, run) in enumerate(runs)
        rng = StableRNG(ex_s53('f', ri))
        diffs = String[]
        count = Ref(0)
        ex_each_orbit(run) do it, o, ch, c
            # run a: every item; runs b, c: a seeded third of the items (their datasets equal run a bit for bit, D53g)
            ri == 1 || rand(rng) < 1 / 3 || return
            for idx in sort(unique(rand(rng, 1:length(o.t), 5)))
                n = o.t[idx]
                direct = ex_direct(it, n)
                stored = o.lyap[it.method][idx, :]
                count[] += 1
                direct == stored || push!(diffs, "item $(it.id) $(it.method) $(it.precision) n=$n: $direct vs $stored")
            end
        end
        @test count[] >= (ri == 1 ? 5 * run.plan.counts.items ÷ 2 : 20)
        @test isempty(diffs)
    end
end

@testset "D53g stored datasets bit-identical for any capacity, order and thread count (I12); D53i threads" begin
    a, b, c = ex_full_runs()
    ref = ex_store_canon(a.store)
    @test !isempty(ref)
    @test [basename(f) for f in ex_chunk_files(b.store)] == [basename(f) for f in ex_chunk_files(a.store)]
    @test ex_store_canon(b.store) == ref                                 # threads 1, capacity 1, chunks in reverse order
    @test ex_store_canon(c.store) == ref                                 # threads 2, capacity 7
    used(run) = EX.read_json(run.res.manifest)["threads"]["used"]
    @test (used(a), used(b), used(c)) == (4, 1, 2)                       # the manifest records the thread count used
    @test all(run -> EX.read_json(run.res.manifest)["threads"]["available"] == Threads.nthreads(), (a, b, c))
    # launch capacity as recorded (inputs of R3 item 3) and chunk files carry the grouping attributes
    capof(run) = EX.read_json(run.res.manifest)["launch"]["capacity"]
    @test capof(b) == 1 && capof(c) == 7 && capof(a) >= 1
    # reverse order: the run processed the last ordinal first
    @test b.res.ran == reverse(sort(a.res.ran))
    # thread count of the CPU backend: all available threads as well
    all_threads = ex_campaign(ex_fix(:small); store = ex_store("g-allthreads"), threads = Threads.nthreads())
    one_thread = ex_campaign(ex_fix(:small); store = ex_store("g-onethread"), threads = 1)
    @test ex_store_canon(ex_path("g-allthreads")) == ex_store_canon(ex_path("g-onethread"))
    @test all_threads.status == one_thread.status == "complete"
end

@testset "D53r the runner with the second fixture block" begin
    fx = ex_fix(:alt)
    plan = ex_plan(fx)
    res = ex_campaign(fx; store = ex_store("r-alt"), threads = 2)
    @test res.status == "complete"
    files = ex_chunk_files(ex_path("r-alt"))
    @test length(files) == plan.counts.chunks == 22
    bad = String[]
    for f in files
        c = plan.chunks[ex_ordinal(f)]
        ch = EX.read_chunk(f)
        for o in ch.orbits
            it = ex_item_of(plan, c, o)
            g = ex_ref_grid(it.horizon, it.tau; ppd = 4, first_steps = 6, fracs = [0.25, 0.75])
            o.t == g || push!(bad, "grid $(it.id)")
            n = o.t[end ÷ 2]
            ex_direct(it, n) == o.lyap[it.method][end ÷ 2, :] || push!(bad, "direct $(it.id)")
            Set(keys(o.diagnostics)) == Set(ex_diag_names(it.method)) || push!(bad, "diag $(it.id)")
        end
    end
    @test isempty(bad)
    @test ex_store_canon(ex_path("r-alt")) == ex_store_canon(let s = ex_store("r-alt-1"); ex_campaign(fx; store = s, threads = 1, capacity = 3); s end)
end

@testset "D53s Y5 and exact exponents (O11)" begin
    run = ex_full_run("a"; threads = 4)
    mu = setprecision(() -> Float64((3 + sqrt(BigFloat(5))) / 2), BigFloat, 256)
    A_y5 = Chaoslab.y5_control(mu).sys.A
    A_fx = [3.5 1.0; 2.5 1.0]                                             # tangent of the standard map at (0, 0), K = 2.5
    v = Float64.(EX.deviation_frame(Float64, "pi/4", 2)[:, 1])
    lnmu(A) = setprecision(() -> Float64(log((BigFloat(A[1, 1]) + BigFloat(A[2, 2]) + sqrt((BigFloat(A[1, 1]) + BigFloat(A[2, 2]))^2 - 4)) / 2)), BigFloat, 256)
    cases = ((15, A_y5, 0.0, (mu^2 + mu^-2 + 2), EX_V), (14, A_fx, 2.5, ex_Gfirst(2.5, 0.0), filter(!=("L9"), EX_V)))
    seen = Set{Tuple{Int, String}}()
    ex_each_orbit(run) do it, o, ch, c
        for (id, A, K, G, methods) in cases
            o.id == id && it.method in methods && it.tau == 1 && it.direction == "pi/4" && it.coordinate_mode == "lifted" &&
                ex_label(it.precision) == "Float64" && it.horizon == 200 || continue
            push!(seen, (id, it.method))
            lam = o.lyap[it.method][:, 1]
            bound(n) = ex_BS(n, K, G) / n
            c_fac = ex_c_factor(A, v)
            ln_mu = lnmu(A)
            for (i, n) in enumerate(o.t)
                exact = id == 15 && it.method == "L9" ? ln_mu : ex_lin_lambda(A, v, n)         # L9 is the singular exponent: ln mu for the normal Y5
                @test abs(lam[i] - exact) <= bound(n)
                n >= 20 && @test abs(lam[i] - (ln_mu + log(c_fac) / n)) <= bound(n) + 1e-12   # ln mu + ln(c) / n
            end
            id == 15 && @test abs(c_fac - 1) < 1e-15                                  # unstable direction pi/4 of the symmetric block: c = 1
        end
    end
    @test length(seen) == length(EX_V) + length(EX_V) - 1
end

@testset "D53t CPU reference rerun (O12)" begin
    run = ex_full_run("a"; threads = 4)
    r1 = EX.rerun_check(run.fx.prereg, run.fx.panel, run.store; code_version = EX_CODE, threads = 2, thread_override = true, allow_dirty = true)
    r2 = EX.rerun_check(run.fx.prereg, run.fx.panel, run.store; code_version = EX_CODE, threads = 1, thread_override = true, allow_dirty = true)
    @test [r.item for r in r1] == [r.item for r in r2]                  # the sample is keyed by the item hash, not by the thread count
    @test !isempty(r1) && all(r -> r.bitwise && r.backend == "cpu", r1)
    for arm in run.fx.run["arm"]
        name = arm["name"]
        n_arm = count(it -> name in it.arms, run.plan.items)
        got = [r for r in r1 if r.arm == name]
        @test length(got) == min(3, n_arm)                              # rerun_sample items per arm (3 in the fixture block)
        @test allunique(r.item for r in got) && all(r -> name in run.plan.items[r.item].arms, got)
    end
end

@testset "D53v time and memory of the run itself" begin
    for m in EX_ALL9
        @test EX.hot_loop_bytes(m) == 0                                 # the recording loop allocates nothing after warm-up
    end
    run = ex_full_run("a"; threads = 4)
    file = argmax(f -> filesize(f), ex_chunk_files(run.store))
    ch = EX.read_chunk(file)
    series = sum(sizeof(o.t) + sizeof(o.state) + sum(sizeof, values(o.lyap)) + sum(sizeof, values(o.diagnostics)) for o in ch.orbits)
    EX.write_chunk(ex_path("v-warm.h5"), ch.attrs, ch.orbits)
    bytes = @allocated EX.write_chunk(ex_path("v-write.h5"), ch.attrs, ch.orbits)
    @test bytes < series + 64 * 2^20                                    # writer memory below the series buffers plus 64 MB
end
