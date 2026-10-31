# D93b parameter table of study p3-abstention (E47, P44, P47). Spec: analysis-spec.md, Amendment D93, E47, criterion D93b.
# The table follows the rule of D89 A2 for every key: `value`, `section`, `quote` (a fragment of the section that contains the text) and `text` (the value as the section writes it, which must
# occur in the quote); there is no machine key. 71 keys; checks: quote and text, the mirrors against the p2 table, the SHA-256 of the p2 manifest, the p2 table against the p2 preregistration,
# the value ranges, voter_bits against indicators_i7.

p3_readp(st; params = st.params, prereg = st.prereg, kw...) = AN.p3_read_params(params; prereg, p2_params = st.s2.params, p2_prereg = st.s2.prereg, p2_manifest = st.s2.manifest, kw...)
# a consistent preregistration and table for the parameter set q (the quotes follow the values)
function p3_pair(tag, q; texts = Dict{Symbol, String}(), kw...)
    dir = an_path("p3par-$tag"); mkpath(dir)
    text, _ = p3_prereg_text(q; texts)
    pre = joinpath(dir, "prereg.md"); write(pre, text)
    return pre, p3_write_params(joinpath(dir, "params.toml"), q; texts, kw...)
end
p3_wrongval(v) = v isa AbstractString ? 7 : (v isa AbstractVector || v isa AbstractDict) ? "scalar" : "text"

@testset "D93b a valid table is read: 71 keys, types, values" begin
    st = p3_cached("default")
    got = p3_readp(st)
    @test Set(keys(got)) == Set(P3_PARAM_KEYS) && Set(AN.P3_PARAM_KEYS) == Set(P3_PARAM_KEYS) && length(AN.P3_PARAM_KEYS) == 71
    @test all(k -> isequal(getfield(got, k), getfield(st.q, k)), P3_PARAM_KEYS)
    @test all(k -> getfield(got, k) isa Integer, P3_INT_KEYS) && got.budgets isa Vector{Int} && got.w_range isa Vector{Int} && got.voter_bits isa Vector{String} && got.indicators_i7 isa Vector{String} &&
          got.spec_of isa AbstractDict && got.spec_of["E1"] == "fq1" && got.cost_tg isa Float64 && got.p2_manifest_sha256 isa String && got.z_one_sided == 1.645
    # the exact decimals are available from the text of the key (the amendment: a threshold is the rational of its decimal text)
    @test AN.p3_rat(got.cost_fq1) == 168504 // 10 && AN.p3_rat(got.delta_tol) == 5 // 100 && AN.p3_rat(got.e_cal) == 2 // 10
    # a table written with Julia's string(value) as the text is read as well
    pre, par = p3_pair("julia-text", st.q; texts = Dict(:a_cal => "0.3", :p_hb => "0.8", :delta_tol => "0.05"))
    @test AN.p3_read_params(par; prereg = pre, p2_params = st.s2.params, p2_prereg = st.s2.prereg, p2_manifest = st.s2.manifest).a_cal == 0.3
end

@testset "D93b each key missing and each key mistyped is refused naming the key; `text` missing, mistyped and not in its quote (71 keys each)" begin
    st = p3_cached("default")
    dir = an_path("p3par-each"); mkpath(dir)
    f = joinpath(dir, "t.toml")
    miss = String[]; typ = String[]; notext = String[]; wrongtext = String[]; badtext = String[]
    for k in P3_PARAM_KEYS
        p3_write_params(f, st.q; drop = [k])
        an_rejects(() -> p3_readp(st; params = f), [string(k)]) || push!(miss, string(k))
        p3_write_params(f, st.q; edit = Dict(k => Dict("value" => p3_wrongval(getfield(st.q, k)))))
        an_rejects(() -> p3_readp(st; params = f), [string(k)]) || push!(typ, string(k))
        p3_write_params(f, st.q; edit = Dict(k => Dict("text" => nothing)))
        an_rejects(() -> p3_readp(st; params = f), [string(k)]) || push!(notext, string(k))
        p3_write_params(f, st.q; edit = Dict(k => Dict("text" => 7)))
        an_rejects(() -> p3_readp(st; params = f), [string(k)]) || push!(wrongtext, string(k))
        p3_write_params(f, st.q; edit = Dict(k => Dict("text" => "a text that is in no quote")))
        an_rejects(() -> p3_readp(st; params = f), [string(k)]) || push!(badtext, string(k))
    end
    @test isempty(miss) && isempty(typ) && isempty(notext) && isempty(wrongtext) && isempty(badtext) && length(P3_PARAM_KEYS) == 71
    # an unknown key and a repeated key (a table defined twice) are refused naming the key; a missing file is refused
    p3_write_params(f, st.q; extra = Dict{String, Any}("bogus_key" => Dict("value" => 1, "section" => "x", "quote" => "1", "text" => "1")))
    @test an_rejects(() -> p3_readp(st; params = f), ["bogus_key"])
    p3_write_params(f, st.q)
    write(f, read(f, String) * "\n[budgets]\nvalue = [1]\nsection = \"x\"\nquote = \"1\"\ntext = \"1\"\n")
    @test an_rejects(() -> p3_readp(st; params = f), ["budgets"])
    @test an_rejects(() -> p3_readp(st; params = f * ".none"), ["does not exist"])
end

@testset "D93b a section not found, a quote not in its section, a quote changed alone; numeric text against value" begin
    st = p3_cached("default")
    f = an_path("p3par-quote.toml")
    for k in (:e_cal, :min_set, :budgets, :cost_tg, :spec_of, :z_one_sided, :p2_manifest_sha256, :test_seed, :y_t, :d_expected_total, :voter_bits, :bootstrap_resamples)
        p3_write_params(f, st.q; edit = Dict(k => Dict("quote" => "text that is not in the section", "text" => "not in the")))
        @test an_rejects(() -> p3_readp(st; params = f), [string(k)])
        p3_write_params(f, st.q; edit = Dict(k => Dict("section" => "No such section")))
        @test an_rejects(() -> p3_readp(st; params = f), [string(k)])
        p3_write_params(f, st.q; edit = Dict(k => Dict("quote" => "text that is not in the section")))
        @test an_rejects(() -> p3_readp(st; params = f), [string(k)])
    end
    # a numeric text must equal its value: another value with the same text and quote is refused (D89 A2, review finding 4)
    for (k, v) in ((:a_cal, 0.4), (:cost_tg, 2658.5), (:rho_hb, 0.6), (:z_one_sided, 1.64), (:pp_runs, 61))
        p3_write_params(f, st.q; edit = Dict(k => Dict("value" => v isa Int ? v : v)))
        @test an_rejects(() -> p3_readp(st; params = f), [string(k)])
    end
    # the section of a key must exist in the preregistration given: a preregistration without the cost section refuses the cost keys
    text, _ = p3_prereg_text(st.q; drop_sections = [:COST])
    pre = an_path("p3par-nocost.md"); write(pre, text)
    @test an_rejects(() -> p3_readp(st; prereg = pre), ["cost"])
end

@testset "D93b values out of range, lists and the consistency of budgets, b_primary, N_c, voter_bits" begin
    st = p3_cached("default")
    q = st.q
    for (tag, bad, key) in (("budgets-empty", (; budgets = Int[]), "budgets"), ("budgets-nine", (; budgets = P2_BUDGETS[1:9]), "budgets"), ("b-primary", (; b_primary = 2561), "b_primary"),
                            ("b-primary-out", (; b_primary = 40960), "b_primary"), ("n-c-last", (; budgets = vcat(P2_BUDGETS[1:9], 20_482), N_c = 20_482), "budgets"),
                            ("e-cal-high", (; e_cal = 1.5), "e_cal"), ("e-cal-neg", (; e_cal = -0.1), "e_cal"), ("p-hb-high", (; p_hb = 1.2), "p_hb"), ("shared-high", (; shared_default = 2.0), "shared_default"),
                            ("alpha-nan", (; alpha_hb = NaN), "alpha_hb"), ("c-min-inf", (; c_min = Inf), "c_min"), ("min-voters-0", (; min_voters = 0), "min_voters"),
                            ("pp-runs-0", (; pp_runs = 0), "pp_runs"), ("b-opt-0", (; b_opt = 0), "b_opt"),
                            ("voter-bits-short", (; voter_bits = q.voter_bits[1:4]), "voter_bits"), ("voter-bits-order", (; voter_bits = reverse(q.voter_bits)), "voter_bits"),
                            ("voter-bits-le", (; voter_bits = q.indicators_i7), "voter_bits"), ("indicators-order", (; indicators_i7 = vcat(q.indicators_i7[2:end], "LE")), "voter_bits"),
                            ("cost-neg", (; cost_tg = -1.0), "cost_tg"), ("ratio-neg", (; ratio_d64_tg = -1.0), "ratio_d64_tg"), ("auc-p-high", (; auc_p = 1.5), "auc_p"),
                            ("auc-shrink-high", (; auc_shrink = 1.5), "auc_shrink"), ("step-zero", (; rho_path_second = 1.0), "rho_path"), ("ceil-0", (; ceil_step = 0), "ceil_step"))
        qq = merge(q, bad)
        pre, par = p3_pair(tag, qq)
        @test an_rejects(() -> AN.p3_read_params(par; prereg = pre, p2_params = st.s2.params, p2_prereg = st.s2.prereg, p2_manifest = st.s2.manifest), [key])
    end
    # the lists of integers: a list of text in w_range and v_range, a non-integer budget
    f = an_path("p3par-lists.toml")
    for (k, v) in ((:w_range, ["a"]), (:w_range, [0.5]), (:v_range, ["b"]), (:budgets, [1.5, 3.0]), (:voter_bits, [1, 2]), (:indicators_i7, "LE"))
        p3_write_params(f, q; edit = Dict(k => Dict("value" => v)))
        @test an_rejects(() -> p3_readp(st; params = f), [string(k)])
    end
end

@testset "D93b the mirrors against the p2 table, the SHA-256 of the p2 manifest, the p2 table against the p2 preregistration" begin
    st = p3_cached("default")
    q = st.q
    # each mirror key altered (the p3 preregistration and table consistent with the altered value) is refused naming the key
    alt = Dict(:budgets => vcat(P2_BUDGETS[1:2], 121, P2_BUDGETS[4:end]), :b_primary => 1280, :N_c => 20_480 * 2, :N_1 => 80_001, :n_s => 801, :f_s => 0.25, :lam_w2 => 0.2, :lam_s2 => 1.5, :min_set => 6, :gate_size => 3)
    for (k, v) in alt
        qq = merge(q, NamedTuple{(k,)}((v,)))
        pre, par = p3_pair("mirror-$k", qq)
        words = k in (:N_c,) ? [[string(k)], ["budgets"]] : [[string(k)]]
        @test any(w -> an_rejects(() -> AN.p3_read_params(par; prereg = pre, p2_params = st.s2.params, p2_prereg = st.s2.prereg, p2_manifest = st.s2.manifest), w), words)
    end
    # the hash of the p2 manifest: altered in the table (and its quote) or the manifest edited after the table was written
    pre, par = p3_pair("manifest-sha", merge(q, (; p2_manifest_sha256 = "f"^64)))
    @test an_rejects(() -> AN.p3_read_params(par; prereg = pre, p2_params = st.s2.params, p2_prereg = st.s2.prereg, p2_manifest = st.s2.manifest), ["p2_manifest_sha256"])
    c = p3_copy(st, "params-manifest")
    p3_edit_manifest(c, m -> m["reasons"] = String["edited after the table"])
    @test an_rejects(() -> p3_readp(c), ["p2_manifest_sha256"]) && !ispath(c.out)
    # a p2 table that fails D89 E22 against the p2 preregistration (a key missing, a quote not in its section) is refused naming the p2 key
    d = p3_copy(st, "params-p2table")
    t = TOML.parsefile(d.s2.params); delete!(t, "ratio_threshold"); open(io -> TOML.print(io, t), d.s2.params, "w")
    @test an_rejects(() -> p3_readp(d), ["ratio_threshold"])
    d2 = p3_copy(st, "params-p2quote")
    t = TOML.parsefile(d2.s2.params); t["min_set"]["quote"] = "text that is not in the section"; open(io -> TOML.print(io, t), d2.s2.params, "w")
    @test an_rejects(() -> p3_readp(d2), ["min_set"])
    # a p2 table of another study (the p3 table given as the p2 table) is refused
    @test an_rejects(() -> AN.p3_read_params(st.params; prereg = st.prereg, p2_params = st.params, p2_prereg = st.s2.prereg, p2_manifest = st.s2.manifest), String[])
end
