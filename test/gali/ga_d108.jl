# Rework D108j (2026-11-05): a column with a raw norm below floatmin(T) (subnormal) is a health failure like a zero column (review finding 1).
# Tests first (RED against the stage 03 code, where `gali_values([1,2,3], [5e-324, 5e-324, 0])` gave 1.1649, true 0.8238).
# BigFloat has no reachable subnormal range (floatmin 2e-323228497), so the subnormal cases cover Float64 and Double64.
@testset "D108j subnormal column norm is a health failure" begin
    _ = (GALI, gali_values, check_frame)
    for T in GA_TYPES
        x = SVector{3, T}(0, 0, 0)
        fm = floatmin(T)
        Wsub = SMatrix{3, 2, T}(T(1), T(2), T(3), T(5e-324), T(5e-324), T(0))
        Wsub3 = SMatrix{3, 3, T}(T(1), T(2), T(3), T(5e-324), T(5e-324), T(0), T(0), T(1), T(0))
        st = reducer_init(GALI(), T, Val(3), Val(3))
        st1, W1 = reducer_step(GALI(), st, 4, x, Wsub3)
        @test reducer_failed(GALI(), st1) && W1 === Wsub3
        e = reducer_error(GALI(), st1, 6, :cpu)
        @test e isa NumericalHealthError && e.what == "deviation matrix" && e.column == 2 && e.step == 4 && e.orbit == 6 && e.backend === :cpu
        @test e.threshold == Float64(fm) && isfinite(e.value) && 0 < e.value < Float64(fm)
        e = ga_health_error(() -> reducer_update(GALI(), st, 4, x, Wsub3))
        @test e isa NumericalHealthError && e.column == 2 && e.threshold == Float64(fm)
        g = gali_values(Wsub)
        @test !any(isfinite, g) && !any(v -> isfinite(v) && v > 1, g)
        @test !any(v -> isfinite(v) && v > 1, gali_values(Wsub3))
        e = ga_health_error(() -> check_frame(GALI(), Matrix(Wsub); orbit = 3))
        @test e isa NumericalHealthError && e.what == "initial frame" && e.orbit == 3 && e.step == 0
        # just above floatmin(T): accepted, and the value is that of the scaled frame
        Wok = SMatrix{3, 3, T}(T(1), T(2), T(3), 2 * fm, 4 * fm, T(0), T(0), T(1), T(0))
        stok, _ = reducer_step(GALI(), st, 1, x, Wok)
        @test !reducer_failed(GALI(), stok)
        Wbig = SMatrix{3, 3, T}(T(1), T(2), T(3), T(2), T(4), T(0), T(0), T(1), T(0))
        @test maximum(abs, Float64.(gali_values(Wok)) .- Float64.(gali_values(Wbig))) <= 1e-12
        @test check_frame(GALI(), Matrix(Wok)) > 0
    end
end
