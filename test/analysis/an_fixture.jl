# Fixture self-test. Not a test of the module: it does not call Chaoslab.Analysis and passes before the module exists.
# It shows that a RED result elsewhere in this suite is not caused by a broken fixture (stage 02 audit, RED).

@testset "fixture: the synthetic study is what its closed forms say" begin
    st = an_cached("default")
    @test st.chunks == length(st.specs) && st.chunks > 40
    man = AX.read_json(st.manifest)
    @test man["status"] == "complete" && length(man["chunks"]) == st.chunks && man["prereg"]["header"] == st.prereg_hash
    @test AX.check_prereg(st.prereg).header_sha256 == st.prereg_hash
    @test an_sha256(read(st.panel)) == st.panel_hash == man["panel"]["sha256_end"]
    files = filter(f -> endswith(f, ".h5"), readdir(st.store))
    @test length(files) == st.chunks
    # one chunk read back: shapes, grid, closed-form values
    i = findfirst(s -> s.method == "L5" && s.horizon == st.p.N_1 && s.mode == "lifted" && s.tau == 1 && s.precision == "Float64", st.specs)
    ch = AX.read_chunk(joinpath(st.store, "chunk-" * lpad(string(i), 5, '0') * ".h5"))
    @test ch.attrs["method"] == "L5" && ch.attrs["horizon"] == st.p.N_1 && ch.attrs["complete"] === true && ch.attrs["status"] == "ok"
    @test length(ch.orbits) == 20 && ch.orbits[1].id == 1
    o = ch.orbits[1]
    @test o.t == an_grid(st.p.N_1) && size(o.lyap["L5"]) == (length(o.t), 2) && size(o.state) == (length(o.t), 2)
    @test o.t[end] == st.p.N_1 && (st.p.N_1 ÷ 10) in o.t && (st.p.N_1 ÷ 2) in o.t && all(k in o.t for k in 1:6)
    @test o.lyap["L5"][:, 1] ≈ an_profile(1).lref .* ones(length(o.t))
    t9 = ch.orbits[9]                                      # a T orbit: the integrable closed form, preregistration A1 item 2
    @test t9.lyap["L5"][:, 1] ≈ an_k0.(Float64.(t9.t))
    @test ch.orbits[17].state[:, 1] == fill(0.0, length(ch.orbits[17].t))      # fixed point control stays
    # the grid helper
    g = an_grid(10_000)
    @test g[1:10] == collect(1:10) && 100 in g && 1000 in g && g[end] == 10_000 && issorted(g) && allunique(g)
    # parameters round trip
    p = an_fx_params()
    f = an_write_params(an_path("fx-params.toml"), p)
    @test Symbol.(keys(TOML.parsefile(f))) |> Set == Set(AN_PARAM_KEYS) && length(AN_PARAM_KEYS) == 37
    # the L2 chunk carries the perturbation text, the V chunks do not
    j = findfirst(s -> s.method == "L2" && s.horizon == st.p.N_1 && s.eps == 1e-8, st.specs)
    @test haskey(AX.read_chunk(joinpath(st.store, "chunk-" * lpad(string(j), 5, '0') * ".h5")).attrs, "perturbation")
    @test !haskey(ch.attrs, "perturbation")
end
