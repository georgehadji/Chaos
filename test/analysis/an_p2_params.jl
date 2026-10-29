# D89b parameter table of study p2-indicators (E22, P9). Spec: analysis-spec.md, Amendment D89, E22, criterion D89b.
# Rule A2 (owner decision of the D89 stage 03 checkpoint): every key other than the 11 machine keys carries `quote` and `text`; `text` is the value as the section writes it
# (1e-4, 0.20, 10^3, 95 percent), must occur in `quote`, and `quote` in the named section; the string form of `value` is not required to occur anywhere.

p2_readp(st; params = st.params, kw...) = AN.p2_read_params(params; prereg = st.prereg, p1_params = st.p1.params, kw...)
# a consistent prereg and table pair for the parameter set p (the quotes follow the values, the run blocks follow the windows)
function p2_pair(tag, p; style = :section, texts = Dict{Symbol, String}(), kw...)
    dir = an_path("p2par-$tag"); mkpath(dir)
    text, _ = p2_prereg_text(p; style, texts)
    pre = joinpath(dir, "prereg.md"); write(pre, text)
    return pre, p2_write_params(joinpath(dir, "params.toml"), p; style, texts, kw...)
end
wrongval(v) = v isa AbstractString ? 7 : v isa AbstractVector ? "scalar" : "text"

@testset "D89b a valid table is read" begin
    st = p2_cached("default")
    got = p2_readp(st)
    @test Set(keys(got)) == Set(P2_PARAM_KEYS) && Set(AN.P2_PARAM_KEYS) == Set(P2_PARAM_KEYS) && length(AN.P2_PARAM_KEYS) == 58
    @test all(k -> isequal(getfield(got, k), getfield(st.p, k)), P2_PARAM_KEYS)
    @test all(k -> getfield(got, k) isa Integer, P2_INT_KEYS) && got.budgets isa Vector{Int} && got.directions isa Vector{String}
    # without the P1 table the label keys are not compared; the default numerics-defaults.toml of the repository gives the seed 20261012
    @test AN.p2_read_params(st.params; prereg = st.prereg).bootstrap_seed == AN_SEED
end

@testset "D89b each key missing and each key mistyped is refused naming the key" begin
    st = p2_cached("default")
    dir = an_path("p2par-each"); mkpath(dir)
    f = joinpath(dir, "t.toml")
    miss = String[]; typ = String[]; notext = String[]; badtext = String[]; wrongtext = String[]
    for k in P2_PARAM_KEYS
        p2_write_params(f, st.p; drop = [k])
        an_rejects(() -> p2_readp(st; params = f), [string(k)]) || push!(miss, string(k))
        p2_write_params(f, st.p; edit = Dict(k => Dict("value" => wrongval(getfield(st.p, k)))))
        an_rejects(() -> p2_readp(st; params = f), [string(k)]) || push!(typ, string(k))
        k in P2_MACHINE && continue
        # A2: a missing `text`, a `text` of the wrong type and a `text` that is not in the quote are refused naming the key (47 keys each)
        p2_write_params(f, st.p; edit = Dict(k => Dict("text" => nothing)))
        an_rejects(() -> p2_readp(st; params = f), [string(k)]) || push!(notext, string(k))
        p2_write_params(f, st.p; edit = Dict(k => Dict("text" => 7)))
        an_rejects(() -> p2_readp(st; params = f), [string(k)]) || push!(wrongtext, string(k))
        p2_write_params(f, st.p; edit = Dict(k => Dict("text" => "a text that is in no quote")))
        an_rejects(() -> p2_readp(st; params = f), [string(k)]) || push!(badtext, string(k))
    end
    @test isempty(miss) && isempty(typ) && isempty(notext) && isempty(wrongtext) && isempty(badtext) && length(P2_PARAM_KEYS) - length(P2_MACHINE) == 47
    # an unknown key, and a repeated key (a table defined twice), are refused naming the key
    p2_write_params(f, st.p; extra = Dict{String, Any}("bogus_key" => Dict("value" => 1, "section" => "x", "quote" => "1")))
    @test an_rejects(() -> p2_readp(st; params = f), ["bogus_key"])
    p2_write_params(f, st.p)
    write(f, read(f, String) * "\n[budgets]\nvalue = [1]\nsection = \"x\"\nquote = \"1\"\n")
    @test an_rejects(() -> p2_readp(st; params = f), ["budgets"])
    @test an_rejects(() -> p2_readp(st; params = f * ".none"), ["does not exist"])
end

@testset "D89b values out of range, an empty list and the consistency of budgets" begin
    st = p2_cached("default")
    f = an_path("p2par-range.toml")
    for (k, v) in ((:bootstrap_level, 1.5), (:bootstrap_level, 0.0), (:bootstrap_level, 1.0), (:min_set, 0), (:bootstrap_resamples, 0), (:persist_divisor, 0), (:censor_factor, 0),
                   (:h3_count_min, 0), (:ratio_threshold, NaN), (:h1_u_min, Inf), (:gate_max_share, NaN), (:f_s, Inf))
        p2_write_params(f, st.p; edit = Dict(k => Dict("value" => v)))
        @test an_rejects(() -> p2_readp(st; params = f), [string(k)])
    end
    p2_write_params(f, st.p; edit = Dict(:budgets => Dict("value" => Int[])))
    @test an_rejects(() -> p2_readp(st; params = f), ["budgets"])
    # budgets[k] must be twice windows[k], b_primary an element of budgets, N_c the last budget: a consistent preregistration and table, one inconsistency
    for (tag, bad, key) in (("b-last", p2_params(budgets = vcat(P2_BUDGETS[1:9], 20_482)), "budgets"),
                            ("b-mid", p2_params(budgets = vcat(P2_BUDGETS[1:7], 5_124, P2_BUDGETS[9:10])), "budgets"),
                            ("b-primary", p2_params(b_primary = 2561), "b_primary"),
                            ("b-primary-out", p2_params(b_primary = 40960), "b_primary"))
        pre, par = p2_pair(tag, merge(bad, (; panel_sha256 = st.p.panel_sha256)))
        @test an_rejects(() -> AN.p2_read_params(par; prereg = pre), [key])
    end
end

@testset "D89b the 11 machine keys are compared with the run blocks" begin
    st = p2_cached("default")
    f = an_path("p2par-machine.toml")
    alt = Dict(:N_1 => st.p.N_1 + 1, :N_c => st.p.N_c + 2, :N_dir => st.p.N_dir + 1, :windows => vcat(21, st.p.windows[2:end]),
               :directions => ["pi/4", "0.3", "0.3 + pi/5", "0.3 + 2pi/5", "0.3 + 4pi/5"], :drift_bound => 2e6, :grid_k_values => ["0.5", "6.2831853"],
               :grid_i_count => 4, :grid_j_count => 5, :grid_orbit_count => 25, :panel_build_sha256 => "f"^64)
    @test Set(keys(alt)) == Set(P2_MACHINE)
    for k in P2_MACHINE
        p2_write_params(f, st.p; edit = Dict(k => Dict("value" => alt[k])))
        # windows and N_c are tied to the budgets by the consistency checks, which may report first
        words = k in (:windows, :N_c) ? [[string(k)], ["budgets"]] : [[string(k)]]
        @test any(w -> an_rejects(() -> p2_readp(st; params = f), w), words)
    end
end

@testset "D89b quotes, text, sections, the seed and the label keys" begin
    st = p2_cached("default")
    f = an_path("p2par-quote.toml")
    sample = (:ratio_threshold, :min_set, :budgets, :h2_q_min, :censor_factor, :panel_sha256, :label_lambda_min, :label_rotation_tol)
    for k in sample
        # quote still checked: a quote (with a text that occurs in it) that is not in the section, and a section that does not exist
        p2_write_params(f, st.p; edit = Dict(k => Dict("quote" => "text that is not in the section", "text" => "not in the")))
        @test an_rejects(() -> p2_readp(st; params = f), [string(k)])
        p2_write_params(f, st.p; edit = Dict(k => Dict("section" => "No such section")))
        @test an_rejects(() -> p2_readp(st; params = f), [string(k)])
        # the quote alone changed: the text is no longer in it
        p2_write_params(f, st.p; edit = Dict(k => Dict("quote" => "text that is not in the section")))
        @test an_rejects(() -> p2_readp(st; params = f), [string(k)])
    end
    # the text is the value as the section writes it: the fixture writes 1e-4, 0.20, -4, 10^3, 95 percent, 5th; Julia's string(value) is accepted as well when it is the text
    got = p2_readp(st)
    @test got.le_floor == 1e-4 && got.h1_u_min == 0.2 && got.e1_d_threshold == -4.0 && got.bootstrap_resamples == 1000 && got.bootstrap_level == 0.95 && got.range_lower_percentile == 5.0
    tab = TOML.parsefile(st.params)
    @test tab["le_floor"]["text"] == "1e-4" && occursin("1e-4", tab["le_floor"]["quote"]) && !occursin("0.0001", tab["le_floor"]["quote"])
    pre, par = p2_pair("julia-text", merge(p2_params(), (; panel_sha256 = st.p.panel_sha256)); style = :julia)
    tj = TOML.parsefile(par)
    @test tj["le_floor"]["text"] == "0.0001" && tj["bootstrap_level"]["text"] == "0.95" && tj["bootstrap_resamples"]["text"] == "1000"
    gj = AN.p2_read_params(par; prereg = pre, p1_params = st.p1.params)
    @test all(k -> isequal(getfield(gj, k), getfield(got, k)), P2_PARAM_KEYS)
    # a numeric text must equal its value (A2 with review finding 4): another value with the same text and quote is refused
    p2_write_params(f, st.p; edit = Dict(:ratio_threshold => Dict("value" => 7.5)))
    @test an_rejects(() -> p2_readp(st; params = f), ["ratio_threshold"])
    # bootstrap_seed must equal the seed of numerics-defaults.toml
    pre, par = p2_pair("seed", p2_params(bootstrap_seed = 1, panel_sha256 = st.p.panel_sha256))
    @test an_rejects(() -> AN.p2_read_params(par; prereg = pre), ["bootstrap_seed"])
    dflt = an_path("p2par-defaults.toml"); write(dflt, "seed = 1
")
    @test an_rejects(() -> p2_readp(st; defaults = dflt), ["bootstrap_seed"])
    write(dflt, "seed = $AN_SEED
")
    @test p2_readp(st; defaults = dflt).bootstrap_seed == AN_SEED
    # a label key different from the P1 table key (the four label keys carry a quote and a text of the P2 section that cites them)
    pre, par = p2_pair("label", p2_params(label_rotation_tol = 2e-7, panel_sha256 = st.p.panel_sha256))
    @test an_rejects(() -> AN.p2_read_params(par; prereg = pre, p1_params = st.p1.params), ["label_rotation_tol"])
    @test AN.p2_read_params(par; prereg = pre).label_rotation_tol == 2e-7
    @test all(k -> haskey(tab[string(k)], "quote") && haskey(tab[string(k)], "text"), (:label_lambda_min, :label_growth_ratio, :label_rotation_tol, :label_estimator_bound))
end
