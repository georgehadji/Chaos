# AC1 (parameter part), E0, E17 last rule, invariant I1. Spec: analysis-spec.md, Equations E0 and E17.

@testset "AC1 parameter table E0" begin
    f = an_write_params(an_path("params-ok.toml"), AN_PARAMS)
    p = AN.read_params(f)
    @test p isa NamedTuple && Set(keys(p)) == Set(keys(AN_PARAMS))
    @test all(getfield(p, k) == getfield(AN_PARAMS, k) for k in keys(AN_PARAMS))
    # each key missing: refused, the key named
    base = Dict{String, Any}(string(k) => v for (k, v) in pairs(AN_PARAMS))
    for k in AN_PARAM_KEYS
        d = delete!(copy(base), string(k))
        g = joinpath(an_path("params-missing-$k.toml"))
        open(io -> TOML.print(io, d), g, "w")
        @test an_rejects(() -> AN.read_params(g), [string(k)])
    end
    # an unknown key is refused, named
    g = an_path("params-unknown.toml")
    open(io -> TOML.print(io, merge(base, Dict{String, Any}("bogus_key" => 1))), g, "w")
    @test an_rejects(() -> AN.read_params(g), ["bogus_key"])
    # a value of the wrong type is refused, the key named
    g = an_path("params-type.toml")
    open(io -> TOML.print(io, merge(base, Dict{String, Any}("n_H2" => "ten thousand"))), g, "w")
    @test an_rejects(() -> AN.read_params(g), ["n_H2"])
    # a missing file
    @test an_rejects(() -> AN.read_params(an_path("params-does-not-exist.toml")), ["params-does-not-exist"])
end
