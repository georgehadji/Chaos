# D95b the test-side table (E60): a valid file is read; malformed files and each check (a) to (f) are refused naming the key. Spec: analysis-spec.md, Amendment D95, E60, I41, criterion D95b.
# API proposed by the tests: p3t_read_params(path; prereg, record, params) -> NamedTuple of the 18 keys (P3T_PARAM_KEYS); ArgumentError naming the key (or the check) otherwise.

p3t_rd(I, path = I.tparams; prereg = I.prereg, record = I.record, params = I.params) = AN.p3t_read_params(path; prereg, record, params)

@testset "D95b a valid test-side table is read: the 18 keys, their types and values" begin
    S, I = p3t_cached(), p3t_cached_inputs()
    @test AN.P3T_PARAM_KEYS isa Union{Tuple, Vector} && Set(Symbol.(AN.P3T_PARAM_KEYS)) == Set(P3T_KEYS) && length(AN.P3T_PARAM_KEYS) == 18
    t = p3t_rd(I)
    @test Set(keys(t)) == Set(P3T_KEYS)
    @test t.pilot_cost_l1 == 210.3 && t.pilot_cost_tg == 2701.0 && t.pilot_cost_pa1 == 1790.2 && t.pilot_cost_bk1 == 3050.8 && t.pilot_cost_fq1 == 16990.4
    @test t.target_t_s === 5 && t.target_w2 === 11 && t.target_f === 6 && t.theta_q_c === 2 && t.theta_q_r === 3 && t.theta_w === 1 && t.theta_v === 0 && t.cap_flag === false
    @test t.theta_voters == ["SALI", "FLI", "MEGNO", "E2"] && t.power_hb_planned == 0.31 && t.power_hc_planned == 0.82
    @test t.calibration_record_sha256 == p3ref_sha256_file(I.record) && t.p3_table_sha256 == p3ref_sha256_file(I.params)
    # the module holds none of the values (I41): another pair of targets is read as written
    I2 = p3t_inputs(S, p3_cached("d95"), "params-b"; target_t_s = 6, target_w2 = 12, target_f = 7, values = (; pilot_cost_l1 = 211.5))
    t2 = p3t_rd(I2)
    @test t2.target_t_s == 6 && t2.target_w2 == 12 && t2.target_f == 7 && t2.pilot_cost_l1 == 211.5
end

# one fresh copy of the table per case, with the field edit of the case
function p3t_param_case(I, tag; edit = Dict{Symbol, Any}(), drop = Symbol[], extra = Dict{String, Any}(), prereg = I.prereg)
    f = joinpath(dirname(I.tparams), "case-$tag.toml")
    p3t_write_values(f, I.t, I.q; drop, extra, edit)
    return f, prereg
end

@testset "D95b more than 2 x 18 + 10 malformed files are refused naming the key" begin
    S, I = p3t_cached(), p3t_cached_inputs()
    bad_value = Dict{Symbol, Any}(:pilot_cost_l1 => "x", :pilot_cost_tg => "x", :pilot_cost_pa1 => "x", :pilot_cost_bk1 => "x", :pilot_cost_fq1 => "x", :calibration_record_sha256 => 5,
        :p3_table_sha256 => 5, :target_t_s => "x", :target_w2 => "x", :target_f => 6.5, :theta_voters => 5, :theta_q_c => "x", :theta_q_r => "x", :theta_w => "x", :theta_v => "x",
        :power_hb_planned => "x", :power_hc_planned => "x", :cap_flag => "false")
    bad = String[]; n = 0
    for k in P3T_KEYS
        f, pre = p3t_param_case(I, "miss-$k"; drop = [k]); n += 1
        an_rejects(() -> p3t_rd(I, f), [string(k)]) || push!(bad, "missing $k")
        f, pre = p3t_param_case(I, "type-$k"; edit = Dict{Symbol, Any}(k => Dict("value" => bad_value[k]))); n += 1
        an_rejects(() -> p3t_rd(I, f), [string(k)]) || push!(bad, "mistyped $k")
        f, pre = p3t_param_case(I, "text-$k"; edit = Dict{Symbol, Any}(k => Dict("text" => "zz-not-in-the-quote"))); n += 1
        an_rejects(() -> p3t_rd(I, f), [string(k)]) || push!(bad, "text $k")
        f, pre = p3t_param_case(I, "quote-$k"; edit = Dict{Symbol, Any}(k => Dict("quote" => "$(k) = zz-not-in-the-section", "text" => "zz-not-in-the-section"))); n += 1
        an_rejects(() -> p3t_rd(I, f), [string(k)]) || push!(bad, "quote $k")
    end
    @test n == 72 && isempty(bad)
    # ten more: an unknown key, a repeated key, a value out of range (five), the section of A1 not found (the preregistration without Amendment A1, and a changed section name), a text field missing
    f, _ = p3t_param_case(I, "unknown"; extra = Dict{String, Any}("bogus_key" => Dict("value" => 1, "section" => "x", "quote" => "x", "text" => "x")))
    @test an_rejects(() -> p3t_rd(I, f), ["bogus_key"])
    f, _ = p3t_param_case(I, "repeat"); write(f, read(f, String) * "
[target_f]
value = 6
section = \"x\"
quote = \"1\"
text = \"1\"
")
    @test an_rejects(() -> p3t_rd(I, f), ["target_f"])
    for (k, v) in ((:target_t_s, -1), (:target_f, 0), (:theta_q_c, 0), (:power_hb_planned, 1.5), (:power_hc_planned, -0.1))
        f, _ = p3t_param_case(I, "range-$k"; edit = Dict{Symbol, Any}(k => Dict("value" => v)))
        @test an_rejects(() -> p3t_rd(I, f), [string(k)])
    end
    text, _ = p3t_prereg_text(I.q, I.t; a1 = false)
    pre = joinpath(dirname(I.tparams), "no-a1.md"); write(pre, text)
    @test an_rejects(() -> p3t_rd(I; prereg = pre), ["Amendment A1"])
    f, _ = p3t_param_case(I, "section"; edit = Dict{Symbol, Any}(:target_f => Dict("section" => "Amendment A2")))
    @test an_rejects(() -> p3t_rd(I, f), ["target_f"])
    f, _ = p3t_param_case(I, "pilot-section"; edit = Dict{Symbol, Any}(:pilot_cost_l1 => Dict("section" => P3T_A1)))
    @test an_rejects(() -> p3t_rd(I, f), ["pilot_cost_l1"])
    # the preregistration without the pilot sentence (the cost section has no such line)
    text, _ = p3t_prereg_text(I.q, I.t); text2 = replace(text, "pilot cost of the p2 pilot, pilot_cost_tg = 2701.0 ns per orbit step" => "no pilot sentence")
    pre = joinpath(dirname(I.tparams), "no-pilot.md"); write(pre, text2)
    @test an_rejects(() -> p3t_rd(I; prereg = pre), ["pilot_cost_tg"])
end

@testset "D95b each check (a) to (f) of E60 fails alone, on a prepared record, and is refused" begin
    S, st = p3t_cached(), p3_cached("d95")
    mk(tag; kw...) = p3t_inputs(S, st, "e60-$tag"; kw...)
    ed(f) = rec -> f(rec)
    # (a) the SHA-256 of the record file differs from calibration_record_sha256; the SHA-256 of the 71-key table differs from p3_table_sha256
    I = mk("a1"); write(I.record, read(I.record, String) * " ")
    @test an_rejects(() -> p3t_rd(I), ["calibration_record_sha256"])
    I = mk("a2"); write(I.params, read(I.params, String) * "\n")
    @test an_rejects(() -> p3t_rd(I), ["p3_table_sha256"])
    # (a) the value written in the record differs from the table's hash
    I = mk("a3"; record_edit = r -> merge(r, (; fingerprints = merge(r.fingerprints, (; p3_parameters = "f"^64)))))
    @test an_rejects(() -> p3t_rd(I), ["p3_table_sha256"])
    # (b) the outcome of the record is not "selected": the study stopped at calibration
    for out in ("no feasible rule", "sets too small")
        I = mk("b-$(replace(out, ' ' => '-'))"; record_edit = r -> merge(r, (; outcome = merge(r.outcome, (; selection = out)))))
        @test an_rejects(() -> p3t_rd(I), [out])
    end
    # (c) theta*, the targets, the cap flag and the planned powers of the record equal the keys of the table
    th(f) = r -> merge(r, (; selection = merge(r.selection, (; theta = merge(r.selection.theta, f)))))
    for (tag, edit, key) in (("qc", th((; qC = 3)), "theta_q_c"), ("qr", th((; qR = 4)), "theta_q_r"), ("w", th((; w = 0)), "theta_w"), ("v", th((; v = 1)), "theta_v"), ("voters", th((; V = ["SALI", "FLI", "MEGNO", "WBA"])), "theta_voters"),
                             ("tt", r -> merge(r, (; pp_hb = merge(r.pp_hb, (; targets = merge(r.pp_hb.targets, (; n_T = 6)))))), "target_t_s"),
                             ("tw", r -> merge(r, (; pp_hb = merge(r.pp_hb, (; targets = merge(r.pp_hb.targets, (; n_W = 12)))))), "target_w2"),
                             ("tf", r -> merge(r, (; pp_hc = merge(r.pp_hc, (; target = 7)))), "target_f"),
                             ("cap", r -> merge(r, (; pp_hb = merge(r.pp_hb, (; targets = merge(r.pp_hb.targets, (; cap = true)))))), "cap_flag"),
                             ("phb", r -> merge(r, (; pp_hb = merge(r.pp_hb, (; path = [merge(r.pp_hb.path[1], (; power = 0.5))])))), "power_hb_planned"),
                             ("phc", r -> merge(r, (; pp_hc = merge(r.pp_hc, (; power_expected = merge(r.pp_hc.power_expected, (; power = 0.5)))))), "power_hc_planned"))
        I = mk("c-$tag"; record_edit = edit)
        @test an_rejects(() -> p3t_rd(I), [key])
    end
    # (d) theta_voters not a subset of I_6 or fewer than min_voters; q_C, q_R, w, v out of their ranges; b_primary / 2^w not in B (the record and the table agree, so only (d) fails)
    for (tag, vals, key) in (("vsub", (; theta_voters = ["SALI", "FLI", "MEGNO", "LE"]), "theta_voters"), ("vmin", (; theta_voters = ["SALI", "FLI"]), "theta_voters"),
                             ("qcrange", (; theta_q_c = 5), "theta_q_c"), ("qrrange", (; theta_q_r = 2), "theta_q_r"), ("wrange", (; theta_w = 7), "theta_w"), ("vrange", (; theta_v = 5), "theta_v"))
        I = mk("d-$tag"; values = vals, record_edit = r -> merge(r, (; selection = merge(r.selection, (; theta = (V = String.(collect(get(vals, :theta_voters, ["SALI", "FLI", "MEGNO", "E2"]))),
              qC = get(vals, :theta_q_c, 2), qR = get(vals, :theta_q_r, 3), w = get(vals, :theta_w, 1), v = get(vals, :theta_v, 0)))))))
        @test an_rejects(() -> p3t_rd(I), [key])
    end
    # (e) the targets against the floors and the caps: below n_t_min and n_w_min, above y_t x r_max = 9 and y_w x r_max = 15 (exact), and at the limits (accepted)
    for (tag, vals, key) in (("tlow", (; target_t_s = 2), "target_t_s"), ("wlow", (; target_w2 = 4), "target_w2"), ("thigh", (; target_t_s = 10), "target_t_s"), ("whigh", (; target_w2 = 16), "target_w2"))
        I = mk("e-$tag"; target_t_s = get(vals, :target_t_s, 5), target_w2 = get(vals, :target_w2, 11))
        @test an_rejects(() -> p3t_rd(I), [key])
    end
    for (tt, tw) in ((3, 5), (9, 15))
        I = mk("e-ok-$tt"; target_t_s = tt, target_w2 = tw)
        @test p3t_rd(I).target_t_s == tt
    end
    # (f) the record's test_seed and r_max differ from those of the 71-key table
    for (tag, k, v) in (("seed", :test_seed, 1), ("rmax", :r_max, 4.0))
        I = mk("f-$tag"; record_edit = r -> merge(r, (; parameters = merge(r.parameters, NamedTuple{(k,)}((v,))))))
        @test an_rejects(() -> p3t_rd(I), [string(k)])
    end
end
