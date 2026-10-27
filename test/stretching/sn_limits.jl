# Stage 04 review finding R1 (stretching-numbers): the reducer refuses the handoff that would wrap a 4-byte count (more than 10^9 handoffs), as a recorded failure.
@testset "AC7e handoff limit inside the reducer" begin
    S = Chaoslab.Indicators
    r = StretchingNumbers(; range = (-1.0, 1.0), bins = 10, interval = 100)
    st0 = reducer_init(r, Float64, Val(2), Val(2))
    ws = zeros(UInt32, reducer_workspace(r, Float64, Val(2), Val(2)))
    W = SMatrix{2, 2}(1.0, 0.0, 0.0, 1.0)
    x = SVector(0.1, 0.2)
    st = S.StretchState{Float64}(st0.s1, st0.s2, st0.q1, st0.q2, 0.0, 10^9, 10^9, 0, Int8(0), Int8(0))
    st1, _ = reducer_step(r, st, 10^9 + 1, x, W, ws)                # handoff number 10^9 + 1 (tau = 1)
    @test reducer_failed(r, st1) && all(iszero, ws)                  # nothing counted
    e = reducer_error(r, st1, 7, :cpu)
    @test e isa NumericalHealthError && e.what == "handoffs" && e.orbit == 7 && e.step == 10^9 + 1
    st2 = S.StretchState{Float64}(st0.s1, st0.s2, st0.q1, st0.q2, 0.0, 10^9 - 1, 10^9 - 1, 0, Int8(0), Int8(0))
    st3, _ = reducer_step(r, st2, 10^9, x, W, ws)                    # handoff number 10^9 is the last allowed
    @test !reducer_failed(r, st3) && sum(ws) == 2
end
