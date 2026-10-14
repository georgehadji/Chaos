using Test
using Chaoslab

@testset "Chaoslab" verbose = true begin   # verbose: per-suite times for the D35 ceilings
    @testset "quality" begin
        include("quality.jl")
    end
    @testset "layering" begin
        include("layering.jl")
    end
    # Module suites register themselves by existing: every test/<module>/runtests.jl runs here.
    for dir in sort(readdir(@__DIR__; join = true))
        suite = joinpath(dir, "runtests.jl")
        if isdir(dir) && isfile(suite)
            @testset "$(basename(dir))" begin
                include(suite)
            end
        end
    end
end
