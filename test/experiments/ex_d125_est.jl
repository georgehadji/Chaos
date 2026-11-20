# Amendment D125, stage 02: D125n (continuous QR items, R38, E59) and D125o (covariant-vector items, R39, E60, E61).
# Seeds: d125_seed('n' | 'o', sub). References are the library calls themselves (propagate_cqr_series, propagate_clv_orbit) on the item's own texts and the
# independent recomputation of the angle summary of O34 from the series; the module is the only missing part.

const D125_HH_RUN = (energy = "1/12", w = ("-0.3", "0.4", "-0.25", "0.25"))

# E59 targets: the grid of D53 (the fixture grid: points per decade 5, first 10 steps, extra fractions 0.1 and 0.5) mapped to the first multiple of m not before
d125_cqr_targets(N, m) = sort(unique(m * cld(t, m) for t in ex_ref_grid(N, 1; ppd = 5, first_steps = 10, fracs = (0.1, 0.5))))

function d125_cqr_ref(it, spec)
    fl = Chaoslab.HenonHeilesFlow(); sys = d125_item_system(it)
    cq = Chaoslab.ContinuousQR(sys, fl; order = spec["order"], orthogonality_tol = spec["orthogonality_tol"])
    k = spec["vectors"]
    Q0 = Chaoslab.initial_deviations(parse(Int, it.run_seed), it.orbit, Val(4), Val(k), Float64)
    _, _, fin, ser = Chaoslab.propagate_cqr_series(cq, Chaoslab.SVector{4, Float64}(d125_start(it)...), Q0, it.horizon; every = 1, orbit = it.orbit)
    return Dict(r.n => r for r in ser)
end

@testset "D125n continuous QR items: prefix against the library loop, targets, refusals, orbit, guard (R38, E59)" begin
    specs = Any[d125_cqr(id = "cq2", order = 2, vectors = 4), d125_cqr(id = "cq4", order = 4, vectors = 2), d125_cqr(id = "cq1", order = 4, vectors = 1)]
    arms = Any[d125_arm("cq", "P", ["L5"], [160]; ests = ["cq2", "cq4", "cq1"], cells = ["hh1"], vectors = "d", frames = ["a"], integrators = ["lf"])]
    run = d125_campaign("n-cqr", d125_run(arms; estimators = specs))
    sp = Dict(s["id"] => s for s in specs)
    got = Dict{Tuple{Int, String}, Any}(); nrows = 0
    d125_each(run) do it, raw, c, file
        it.system == "henon_heiles" || return
        o = it.orbit; pre = "orbits/$o/"
        if it.method == "L5"
            got[(o, "L5")] = raw[pre * "state"]; return
        end
        spec = sp[it.method]; k = spec["vectors"]; m = spec["order"] ÷ 2
        ser = raw[pre * "estimators/$(it.method)/series"]; t = raw[pre * "t"]
        @test size(ser) == (length(t), k + 1) && raw[pre * "estimators/$(it.method)/series@columns"] == [["lambda$i" for i in 1:k]..., "orth"]
        @test Int.(t) == d125_cqr_targets(it.horizon, m)                                                            # first multiples of m not before the D53 targets
        ref = d125_cqr_ref(it, spec); h = d125_item_h(it)
        @test all(r -> ser[r, 1:k] == ref[Int(t[r])].l && ser[r, k + 1] == ref[Int(t[r])].E, eachindex(t))            # bitwise at every target
        nrows += length(t)
        got[(o, it.method)] = raw[pre * "state"]
        @test raw["@exponent_unit"] == "time" && !haskey(raw, pre * "lyap/" * it.method)
        @test raw["@kind"] == "cqr" && raw["@order"] == spec["order"] && raw["@orthogonality_tol"] == spec["orthogonality_tol"]
    end
    @test nrows >= 3 * 3 * 10
    # the orbit state of a cqr item equals that of the L5 item of the same orbit and integrator at every common step (grids are subsets: compare the last row)
    @test all(o -> got[(o, "cq4")][end, :] == got[(o, "L5")][end, :], [o for o in 1:3 if haskey(got, (o, "cq4")) && haskey(got, (o, "L5"))]) && length(got) >= 12
    # refusals: tau other than 1, a horizon that is not a multiple of m (order 4: m = 2), a map cell
    bad(f, keys, tag) = (r = d125_run(Any[f(d125_arm("cq", "P", nothing, [160]; ests = ["cq4"], cells = ["hh1"], vectors = "d", frames = ["a"], integrators = ["lf"]))];
                                      estimators = Any[d125_cqr(id = "cq4", order = 4, vectors = 2)]); d125_block_refused(r, keys; tag = tag))
    @test bad(a -> (a["tau"] = [2]; a), ["tau", "cq4", "cqr"], "cqr-tau")
    @test bad(a -> (a["horizons"] = [161]; a), ["horizon", "multiple", "cq4"], "cqr-horizon")
    @test bad(a -> (a["cells"] = ["y3a"]; a), ["cqr", "flow", "y3a", "henon_heiles"], "cqr-map")
    # a guard below the realized defect fails the chunk with the library's health error, deterministically on the rerun; CPU only with the sixth reason
    tight = d125_run(Any[d125_arm("cq", "P", nothing, [160]; ests = ["cq2"], cells = ["hh1"], vectors = "d", frames = ["a"], integrators = ["lf"])];
                     estimators = Any[d125_cqr(id = "cq2", order = 2, vectors = 4, orthogonality_tol = 1e-9)])
    rt = d125_campaign("n-tight", tight)
    @test rt.res.status == "complete with failed chunks" && !isempty(rt.res.failed)
    ch = EX.read_chunk(first(ex_chunk_files(rt.store)))
    @test ch.attrs["status"] == "failed" && ch.attrs["attempt"] == 2 && any(r -> occursin("orthogonality", r) && occursin("orbit", r), ch.attrs["failure_reasons"])
    for req in ("auto", "gpu")
        rr = ex_with_device(ExFakeDevice(true, true)) do
            EX.resolve_run_backend(req, "Float64", d125_plan(d125_fix("std")).items[end - 8]; kind = "cqr")                 # assumed: kind keyword, plan
        end
        @test rr.backend == :cpu && rr.reason == Chaoslab.Foundation.FALLBACK_REASONS[6]
    end
end

# E60 from the definition: n*, n_f = n_b and U of a pilot spectrum, or the error kind (exact rational arithmetic where the spec gives integers)
function d125_clv_ref(lam::AbstractVector, spec::AbstractDict, tau::Integer, h::Real)
    d = length(lam)
    gaps = [lam[i] - lam[i + 1] for i in 1:(d - 1) if !(i in spec["gap_exclude"])]
    (isempty(gaps) || any(g -> !isfinite(g) || g <= 0, gaps)) && return :pilot
    gmin = minimum(gaps)
    nstar = Int(setprecision(() -> ceil(BigFloat(spec["transient_c"]) / (BigFloat(gmin) * tau * BigFloat(h))), BigFloat, 256))
    nf = min(spec["transient_cap"], spec["transient_factor"] * nstar)
    U = [i for i in 1:d if lam[i] >= spec["unstable_fraction"] * lam[1]]
    (isempty(U) || length(U) > d ÷ 2) && return :unstable_set
    return (nstar = nstar, nf = nf, nb = nf, U = U, gmin = gmin)
end

@testset "D125o covariant vectors: transients of synthetic pilot spectra (E60)" begin
    rng = StableRNG(d125_seed('o'))
    spec = d125_clv(gap_exclude = [2], transient_c = 5.0, transient_factor = 2, transient_cap = 30, unstable_fraction = 0.5)
    spec_all = d125_clv(gap_exclude = Int[])
    cases = Any[]
    for t in 1:8                                                                                    # regular spectra: Hamiltonian pairs with a zero block in the middle
        a = 0.05 + 0.2 * rand(rng); b = 0.01 + 0.5 * a * rand(rng)
        push!(cases, (spec, [a, b, -b, -a], 1, 0.0625))
    end
    for t in 1:4; push!(cases, (spec_all, [0.4, 0.1 + 0.05t, -0.1 - 0.05t, -0.4], 1, 0.125)); end       # no exclusion: gap 2 enters
    push!(cases, (spec, [0.3, 0.0, 0.0, -0.3], 1, 0.0625))                                          # the gap across the zero block is excluded: gaps 1 and 3 remain
    push!(cases, (spec_all, [0.3, 0.0, 0.0, -0.3], 1, 0.0625))                                      # not excluded: the gap 0 is non-positive: pilot error
    push!(cases, (spec, [0.3, 0.3, -0.3, -0.3], 1, 0.0625))                                         # a non-positive remaining gap: pilot error
    push!(cases, (spec, [0.3, 0.1, 0.1, -0.5], 1, 0.0625))                                          # (this spectrum has remaining gaps 0.2 and 0.6)
    push!(cases, (spec, [1e-9, 5e-10, -5e-10, -1e-9], 1, 0.0625))                                   # a cap that binds: n* huge
    push!(cases, (spec, [0.5, 0.2, -0.2, -0.5], 2, 0.125))                                          # tau = 2
    push!(cases, (spec, [-0.1, -0.2, -0.3, -0.4], 1, 0.0625))                                       # lam_1 < 0: U = {i : lam_i >= f lam_1}
    push!(cases, (spec, [0.0, 0.0, 0.0, 0.0], 1, 0.0625))                                           # all zero: gaps zero, pilot error
    push!(cases, (spec, [0.4, 0.3, 0.1, -0.2, -0.5, -0.6], 1, 0.0625))                              # d = 6, U has 2 members of 3 allowed
    push!(cases, (d125_clv(gap_exclude = Int[], unstable_fraction = 0.05), [0.4, 0.3, 0.1, -0.2, -0.5, -0.6], 1, 0.0625))   # U has more than d/2 members
    push!(cases, (d125_clv(gap_exclude = Int[], unstable_fraction = 3.0), [0.4, 0.3, 0.1, -0.2, -0.5, -0.6], 1, 0.0625))   # U empty
    push!(cases, (spec, [NaN, 0.1, -0.1, -0.2], 1, 0.0625))                                         # not finite: pilot error
    push!(cases, (spec_all, [0.6, 0.45, 0.3, 0.1, -0.2, -0.5, -0.6, -0.7], 1, 0.03125))
    @test length(cases) >= 20
    for (sp, lam, tau, h) in cases
        want = d125_clv_ref(lam, sp, tau, h)
        if want isa Symbol
            @test ex_rejects(() -> EX.clv_transients(lam, sp, tau, h), [string(want) == "pilot" ? "pilot" : "unstable set"])
        else
            got = EX.clv_transients(lam, sp, tau, h)
            @test (got.nstar, got.nf, got.nb, collect(got.U)) == (want.nstar, want.nf, want.nb, want.U)
        end
    end
    @test count(c -> d125_clv_ref(c...) isa Symbol, cases) >= 5 && any(c -> (r = d125_clv_ref(c...); !(r isa Symbol) && r.nf == c[1]["transient_cap"]), cases)
end

@testset "D125o covariant vectors: items equal the library, summary from the series, checkpoint, threads, memory, failure (R39, E61, E62)" begin
    store_spec = d125_clv(id = "cl", storage = "store", store_angles = "summary")
    ckpt_spec = d125_clv(id = "cl", storage = "checkpoint", block = 8, store_angles = "summary")
    ser_spec = d125_clv(id = "cl", storage = "store", store_angles = "series")
    arm = Any[d125_arm("clv", "P", nothing, [100]; ests = ["cl"], cells = ["hh1"], vectors = "d", frames = ["a"], integrators = ["lf"])]
    runs = Dict(nm => d125_campaign("o-" * nm, d125_run(arm; estimators = Any[sp]); threads = th)
                for (nm, sp, th) in (("store", store_spec, 4), ("ckpt", ckpt_spec, 4), ("series", ser_spec, 4), ("t1", store_spec, 1)))
    # the item's results equal the library run with the same transients, bitwise
    d125_each(runs["series"]) do it, raw, c, file
        it.system == "henon_heiles" || return
        o = it.orbit; pre = "orbits/$o/"; w = raw[pre * "estimators/cl/theta_min"]
        nf = Int(raw[pre * "@n_f"]); nb = Int(raw[pre * "@n_b"]); u = Int(raw[pre * "@unstable_count"])
        clv = Chaoslab.ClvGinelli(d125_item_system(it); tau = it.tau, h = d125_item_h(it), orthonormalizer = Chaoslab.HouseholderQR(), storage = :store, block = 0)
        r = Chaoslab.propagate_clv_orbit(clv, Chaoslab.SVector{4, Float64}(d125_start(it)...), nf, ser_spec["window"], nb; k = 4, seed = parse(Int, it.run_seed),
                                         orbit = o, outputs = (), angle_sets = (collect(1:u), collect((5 - u):4)))
        @test raw[pre * "estimators/cl/lambda_fwd"] == r.lambda_fwd && vec(w) == r.angles[1, :] && length(w) == ser_spec["window"] + 1
        @test raw[pre * "@window_start_state"] == Float64.(collect(r.x0)) && raw[pre * "@window_end_state"] == Float64.(collect(r.x))
        @test u == 1 && length(raw[pre * "@pilot_lambda"]) == 4                                                       # F3 chaotic start: U = {1}, the stable set {4}
        # the summary equals the recomputation from the series of the same orbit, 0 difference (O34): count, min, max, nearest-rank quantiles, below counts, histogram
        sm = runs["store"]
        rs = d86_raw(ex_chunk_files(sm.store)[ex_ordinal(file)]); s = rs[pre * "estimators/cl/theta_summary"]
        th = sort(vec(w)); n = length(th); probs = ser_spec["probabilities"]
        want = Float64[n, th[1], th[end]]
        for p in probs
            push!(want, th[clamp(Int(cld(numerator(d125_rat(p)) * n, denominator(d125_rat(p)))), 1, n)])
        end
        for b in ser_spec["thresholds"]
            push!(want, count(<(b), th))
        end
        @test vec(s) == want
        lg = log10.(th); nbin = ser_spec["hist_bins"]; lo = ser_spec["hist_lo"]; hi = ser_spec["hist_hi"]
        hist = zeros(Int, nbin + 2)
        for v in lg
            hist[v < lo ? 1 : v >= hi ? nbin + 2 : 2 + floor(Int, (v - lo) / (hi - lo) * nbin)] += 1
        end
        @test vec(rs[pre * "estimators/cl/theta_hist"]) == hist && sum(hist) == n
    end
    # storage checkpoint equals store bitwise; 1, 4 and all threads give bitwise equal files
    @test d86_store_sha(runs["ckpt"].store; skip = ("@block", "@storage")) == d86_store_sha(runs["store"].store; skip = ("@block", "@storage"))
    @test d86_store_sha(runs["t1"].store) == d86_store_sha(runs["store"].store)
    # E62: the memory formula of the module equals clv_bytes of the library; a start whose one orbit does not fit is refused naming kind, bytes and free memory
    plan = runs["store"].plan; it = first(i for i in plan.items if i.kind == "clv")
    nb = EX.orbit_bytes(plan, it)
    @test nb.B == Chaoslab.clv_bytes(Val(4), Val(4), Float64; n_s = 40, n_b = 30, storage = :store, block = 0)
    free = nb.total - 1
    @test ex_rejects(() -> EX.group_capacity(plan, it; free_bytes = free, f = 1.0), ["clv", string(nb.total), string(free)])
    # a failing orbit fails its chunk naming the orbit; deterministic on the rerun (the pilot gap guard: a spectrum with a zero gap, gap_exclude empty)
    badspec = d125_clv(id = "cl", gap_exclude = Int[], pilot_handoffs = 60)
    rb = d125_campaign("o-bad", d125_run(arm; estimators = Any[badspec]))
    @test rb.res.status in ("complete with failed chunks", "stopped: failed chunks") && !isempty(rb.res.failed)
    cb = EX.read_chunk(first(ex_chunk_files(rb.store)))
    @test cb.attrs["status"] == "failed" && any(r -> occursin("orbit", r) && occursin("pilot", r), cb.attrs["failure_reasons"])
end
