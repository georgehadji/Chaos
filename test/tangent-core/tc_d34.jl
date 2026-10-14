# Amendment D34 (M2 Batch B API contract, A1a and A2), included last from runtests.jl. Existing testsets are untouched.
#   A1a  NumericalHealthError gets a seventh field `column::Int` (default 0); the six-argument constructor keeps working.
#   A2   Tangent.reducer_provenance(r) -> NamedTuple (default NamedTuple()), exported, merged into the provenance of propagate
#        over the chain in order (later reducers win on a shared key).
# RED for this file is mostly Fail or MethodError on existing names (no `column` field yet) plus UndefVarError `reducer_provenance`.

@testset "D34a NumericalHealthError column" begin
    # six-argument constructor: column 0; seven-argument constructor: column set; field order and types; showerror unchanged in
    # content; every error Batch A raises (propagate state / deviation matrix, check_residual, monitor) has column 0.
    e6 = NumericalHealthError("det", 1.5, 0.5, 3, 2, :cpu)
    @test e6 isa Exception
    @test e6.what == "det" && e6.value == 1.5 && e6.threshold == 0.5 && e6.step == 3 && e6.orbit == 2 && e6.backend === :cpu
    @test e6.column == 0
    e7 = NumericalHealthError("orthonormalisation breakdown", 0.0, 0.0, 4, 5, :cpu, 3)
    @test e7.what == "orthonormalisation breakdown" && e7.step == 4 && e7.orbit == 5 && e7.backend === :cpu
    @test e7.column == 3
    @test fieldnames(NumericalHealthError) == (:what, :value, :threshold, :step, :orbit, :backend, :column)
    @test fieldtype(NumericalHealthError, :column) === Int
    msg = sprint(showerror, e7)
    @test occursin("NumericalHealthError", msg) && occursin("orthonormalisation breakdown", msg)
    @test occursin("step 4", msg) && occursin("orbit 5", msg)
    @test occursin("step 3", sprint(showerror, e6)) && occursin("orbit 2", sprint(showerror, e6))
    # errors raised by the Batch A paths carry column 0
    err = tc_health_error(() -> propagate(tc_Fault(:inf_state, 7), [SVector(0.0, 0.0, 1.0)], 20; k = 1, tau = 5, seed = TC_SEED + 3400))
    @test err isa NumericalHealthError && err.what == "state" && err.column == 0
    err = tc_health_error(() -> propagate(tc_Fault(:nan_W, 3), [SVector(0.0, 0.0, 1.0)], 20; k = 1, tau = 5, seed = TC_SEED + 3401))
    @test err isa NumericalHealthError && err.what == "deviation matrix" && err.column == 0
    err = tc_health_error(() -> check_residual("det", 1.0, 0.5; step = 2, orbit = 3))
    @test err isa NumericalHealthError && err.column == 0 && err.what == "det" && err.step == 2 && err.orbit == 3
    sys = LinearSymplecticMap(SMatrix{2, 2, Float64}(2, 1, 1, 1 + 2.0^-20))
    mon = ResidualMonitor(sys, 4.0, Inf)
    err = tc_health_error(() -> propagate_orbit(sys, SVector(0.0, 0.0), SMatrix{2, 1, Float64}(1, 0), 4, 1, (mon,); orbit = 6))
    @test err isa NumericalHealthError && err.what == "det" && err.orbit == 6 && err.column == 0
    # the error is a value type usable in the usual ways (equality of fields, rethrow)
    @test_throws NumericalHealthError throw(e7)
end

@testset "D34b reducer_provenance" begin
    # default NamedTuple() for reducers without a method; exported; merged into propagate's provenance in chain order (later wins);
    # the Batch A keys are unchanged by the merge; reducers without provenance add nothing.
    @test Base.isexported(Chaoslab, :reducer_provenance) && Base.isexported(Chaoslab.Tangent, :reducer_provenance)
    @test reducer_provenance(tc_LogNorm()) === NamedTuple()
    @test reducer_provenance(tc_Reader(1)) === NamedTuple()
    @test reducer_provenance(tc_ProvA()) == (alpha = 1, shared = "a")
    sys = tc_Fault(:inf_state, 0)   # D38: no orbit is flagged (z[3] = 0), so this is a plain run; the map only shares D34a's kernel specialization
    x0 = [SVector(0.0, 0.0, 0.0) for j in 1:4]
    run(reducers) = propagate(sys, x0, 10; k = 1, tau = 2, seed = TC_SEED + 3410, reducers = reducers, chunk = 16).provenance
    base = run(())
    @test haskey(base, :backend) && base.seed == TC_SEED + 3410 && base.tau == 2 && base.k == 1 && base.chunk_size == 16
    ab = run((tc_ProvA(), tc_Reader(2), tc_ProvB()))   # a reducer without a method between two with one adds nothing
    @test ab.alpha == 1 && ab.beta == 2.5
    @test ab.shared == "b"   # D38: the reversed chain (later reducer wins, ba.shared == "a") moved to verify/tangent-core.jl
    @test Set(keys(ab)) == union(Set(keys(base)), Set((:alpha, :beta, :shared)))
    @test all(k -> isequal(ab[k], base[k]), keys(base))
    # D38: the chains without provenance and the mixed chain (reducers with and without a method) moved to verify/tangent-core.jl.
end
