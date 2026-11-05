# Stage 04 review findings of the coupled-maps module (owner decision 2026-11-05: fix all, tests first). Finding 2: the vector forms of the wrapped chain
# return the element type promoted against K and eps, as the lifted chain and the fixed-size forms do. Finding 3: the positional constructor of the chain
# validates N and the ring as the keyword constructor does. Finding D105d (the GPU kernel) is checked at stage 04 (verify/gpu.jl), not here.

@testset "review 2: wrapped chain vector forms promote against K and eps" begin
    rng = StableRNG(CM_SEED + 2101)
    for ring in (false, true)
        w = WrappedStandardMapChain(3, Double64(3), Double64(0.2); boundary = ring ? :ring : :open)
        s = torus_state(w, 4pi .* (2 .* rand(rng, 6) .- 1))
        @test eltype(s) == Float64                                                     # guard: the state is binary64
        @test eltype(forward(w, s)) == Double64
        @test eltype(inverse(w, s)) == Double64
        @test eltype(forward(w, cm105_sv(s))) == Double64                              # guard: the fixed-size form promotes (passes now)
    end
end

@testset "review 3: positional chain constructor validates N and the ring" begin
    @test_throws ArgumentError StandardMapChain(1, 3.0, 0.2, false)
    @test_throws ArgumentError StandardMapChain(101, 3.0, 0.2, false)
    @test_throws ArgumentError StandardMapChain(2, 3.0, 0.2, true)
    @test StandardMapChain(3, 3.0, 0.2, true).ring                                    # guard: a valid ring is accepted
    @test StandardMapChain(2, 3.0, 0.2, false).N == 2
end
