# Rework D109 b to d (2026-11-05), review findings 2 to 4: the d >= K check, type-stable heap path bitwise equal to the SMatrix path,
# integer and other non-float element types. Tests first (RED against the D108 code). Helpers: ga_helpers.jl.
using Test: @inferred
ga_same(a, b) = length(a) == length(b) && all(i -> (a[i] isa BigFloat ? a[i] == b[i] : a[i] === b[i]), eachindex(a))

@testset "D109b gali_values and gali_value need d >= K" begin
    _ = (GALI, gali_values, gali_value, check_frame)
    for T in (Float64, Double64, BigFloat)
        M = T[1 2 3; 0 1 1]                   # 2 x 3
        @test_throws ArgumentError gali_values(M)
        @test_throws ArgumentError gali_values(SMatrix{2, 3, T}(M))
        @test_throws ArgumentError gali_value(M, 2)
        @test_throws ArgumentError gali_value(M, 3)
        @test_throws ArgumentError check_frame(GALI(), M)
        @test_throws ArgumentError gali_values(T[1 0; 0 1; 0 0][:, 1:1])   # one column stays rejected
        @test length(gali_values(T[1 0 0; 0 1 0; 0 0 1])) == 2             # d = K stays accepted
    end
    # constant message on the SMatrix path (it may be reached from device code), D109b
    e = try gali_values(SMatrix{2, 3, Float64}(1.0, 0, 2, 1, 3, 1)); nothing catch err; err end
    @test e isa ArgumentError && !occursin("2", e.msg) && !occursin("3", e.msg)
end

@testset "D109c AbstractMatrix path: inferred, bitwise equal to the SMatrix path" begin
    _ = (gali_values, check_frame)
    rng = StableRNG(GA_BASE + 109)
    for T in (Float64, Double64, BigFloat), d in (4, 8, 200), K in (2, 4)
        (T === BigFloat && d == 200) && (d = 30)    # BigFloat at d = 200 is slow; the path is the same loops
        S = ga_randset(rng, T === BigFloat ? Float64 : T, d, K)
        S = T === BigFloat ? SMatrix{d, K, BigFloat}(BigFloat.(S)) : S
        M = Matrix(S)
        gm = @inferred gali_values(M)
        gs = gali_values(S)
        @test ga_same(gm, gs)
        @test eltype(gm) === T
        @test ga_same([gali_value(M, k) for k in 2:K], gs)
        c = @inferred check_frame(GALI(), M)
        @test c isa T && c == gs[end]
    end
end

@testset "D109d non-float element types are converted with float" begin
    _ = (gali_values, gali_value, check_frame)
    Mi = [1 0 0; 0 1 0; 0 0 1; 0 0 0]       # exact 0/1 entries
    g = gali_values(Mi)
    @test eltype(g) === Float64 && all(==(1.0), g)
    @test gali_values(SMatrix{4, 3}(Mi)) == g && eltype(gali_values(SMatrix{4, 3}(Mi))) === Float64
    Mg = [3 1 4; 1 5 9; 2 6 5; 3 5 8]       # general integers: no InexactError, equals the Float64 evaluation
    @test ga_same(gali_values(Mg), gali_values(Float64.(Mg)))
    @test ga_same(gali_values(Int32.(Mg)), gali_values(Float64.(Mg)))
    @test gali_value(Mg, 3) === gali_values(Float64.(Mg))[2]
    c = check_frame(GALI(), Mi)
    @test c isa Float64 && c == 1.0
    @test check_frame(GALI(), Mg) == gali_values(Float64.(Mg))[end]
    @test_throws NumericalHealthError check_frame(GALI(), [1 2; 1 2; 1 2])    # dependent integer columns
    Mr = Rational{Int}[1 0; 0 1; 1 1]
    @test eltype(gali_values(Mr)) === Float64
end
