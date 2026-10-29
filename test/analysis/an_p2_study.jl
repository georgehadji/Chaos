# Amendment D89, stage 02: the designed orbits and the writer of the synthetic study "p2-indicators" (and its P1 reference study). Definitions only.
# Every number written into the store follows from a designed verdict sequence by the rules of E27 (the self-test an_p2_fixture.jl checks this with the
# reference rule of an_p2_ref.jl). Expected results of the pipeline tests come from the design, never from the module.

# ---------------------------------------------------------------- designed verdict sequences: 'c' chaotic, 'r' regular, 'u' undecided, index j = budget index j - 1
const P2_ONSET = Dict(
    :S => (LE = 0, SALI = 0, FLI = 0, MEGNO = 1, E1 = 4, E2 = 2, WBA = 3, RLI = 0),
    :S2 => (LE = 0, SALI = 1, FLI = 0, MEGNO = 1, E1 = 6, E2 = 4, WBA = 5, RLI = 1),
    :W => (LE = 5, SALI = 4, FLI = 3, MEGNO = 6, E1 = 7, E2 = 6, WBA = 8, RLI = 5),
    :W2 => (LE = 6, SALI = 4, FLI = 4, MEGNO = 5, E1 = 7, E2 = 5, WBA = 8, RLI = 6),
    :N => (LE = 2, SALI = 2, FLI = 2, MEGNO = 2, E1 = 3, E2 = 3, WBA = 3, RLI = 2))
const P2_TSTR = Dict(:LE => "rrrrrrrccc", :SALI => "rrrrrrrrrr", :FLI => "rrrrrccccc", :MEGNO => "rrrrrrrrrr", :E1 => "rrrrrrrrrr", :E2 => "rrrrrrcccc",
                     :WBA => "rrrrrrrrrr", :RLI => "rrrrcccccc")
p2_sym(str::AbstractString) = [c == 'c' ? :chaotic : c == 'r' ? :regular : :undecided for c in str]
p2_str(v) = String([x === :chaotic ? 'c' : x === :regular ? 'r' : 'u' for x in v])

function p2_design(variant::Symbol = :default)
    prof = Dict{Int, Symbol}()
    for i in (1, 3, 7); prof[i] = :W; end
    for i in (8, 21, 22, 26); prof[i] = :W2; end
    for i in (5, 23, 24, 39); prof[i] = :T; end
    for i in (11, 12, 33, 34, 35, 17, 18, 20); prof[i] = :S; end
    for i in (36, 37, 38); prof[i] = :S2; end
    for i in (4, 6, 14, 25, 40, 2, 13, 32, 44); prof[i] = :N; end
    for i in (9, 10, 15, 16, 19, 27:31..., 41, 42, 43); prof[i] = :R; end
    variant === :var && (prof[25] = :T)
    ov = Dict{Tuple{Int, Symbol}, String}((9, :E1) => "cccrrrrrrr", (10, :SALI) => "uurrrrrrrr", (15, :E1) => "uuuuuuuuuu", (15, :WBA) => "uuuuuuuuuu",
        (16, :LE) => "rrrrrrrrrc", (28, :MEGNO) => "rrrrrrrrcc", (29, :FLI) => "rrrrrrrurr", (41, :SALI) => "ccrrrrrrrr")
    variant === :var && (ov[(10, :E1)] = "rrrrrrrrrc")
    variant === :var && foreach(i -> ov[(i, :RLI)] = "rrrrrrrrcc", (1, 21, 22, 26))                  # late RLI onset on W2 orbits of E8: the extension differs
    rho = Dict{Int, Symbol}()
    for i in (1, 3, 4, 5, 6, 7, 8, 11, 12, 14, 18, 20, 21:26..., 33:40...); rho[i] = :C; end
    for i in (9, 10, 15, 16, 27:31..., 41, 42, 43); rho[i] = :R; end
    for i in (2, 13, 17, 19, 32, 44); rho[i] = :U; end
    lref = Dict{Int, Float64}(i => (rho[i] === :R ? 1.4e-4 : 0.3) for i in 1:44)
    lns = copy(lref)
    for (i, a, b) in ((1, 0.1, 0.1), (3, 0.08, 0.08), (4, nextfloat(0.1), nextfloat(0.1)), (5, 0.5, 0.1), (6, 0.5, nextfloat(0.1)), (7, 0.05, 0.05), (8, 0.06, 0.06),
                      (11, 1.0, 1.0), (12, 1.2, 1.2), (14, prevfloat(1.0), prevfloat(1.0)), (21, 0.09, 0.09), (22, 0.02, 0.02), (23, 0.5, 0.05), (24, 0.3, 0.03),
                      (25, 0.3, variant === :var ? 0.03 : 0.3), (26, 0.04, 0.04), (39, 0.5, 0.05), (40, 0.7, 0.7), (17, 0.69, 0.69), (18, 0.8, 0.8), (20, 0.96, 0.96))
        lref[i], lns[i] = a, b
    end
    for (j, i) in enumerate(33:38); lref[i] = lns[i] = 1.1 + 0.02 * (j - 1); end
    lref[32] = lns[32] = 1.4e-4                                    # orbit 32: regular numbers, the structure axis fails (rotation numbers differ)
    seqs = Dict{Int, Dict{Symbol, Vector{Symbol}}}()
    for id in 1:44
        pr = prof[id]
        d = Dict{Symbol, Vector{Symbol}}()
        for X in P2_X8
            s = pr === :R ? "rrrrrrrrrr" : pr === :T ? P2_TSTR[X] : (o = getfield(P2_ONSET[pr], X); repeat("r", o) * repeat("c", 10 - o))
            d[X] = p2_sym(get(ov, (id, X), s))
        end
        seqs[id] = d
    end
    haveRLI(id) = id in P2_SL || id >= 21
    for id in 1:44
        haveRLI(id) || delete!(seqs[id], :RLI)
        id == 20 && (delete!(seqs[id], :E1); delete!(seqs[id], :WBA))     # the Y5 control: no torus
    end
    d64 = Dict{Tuple{Int, Symbol}, Symbol}((5, :E2) => :regular, (11, :LE) => :regular)               # Double64 verdict at k = 9 that differs from Float64
    dir = Dict{Tuple{Int, Symbol, Int}, Symbol}((7, :SALI, 3) => :regular, (12, :LE, 5) => :regular, (15, :FLI, 2) => :undecided)
    drift = variant === :var ? [8] : Int[]
    return (variant = variant, prof = prof, rho = rho, lref = lref, lns = lns, V = seqs, d64 = d64, dir = dir, drift = drift)
end

# ---------------------------------------------------------------- stored numbers that realise a verdict (E26, E27)
p2_rv(v::Symbol, id, k) = v === :chaotic ? 8.0 + 0.1 * (id % 7) + 0.05 * k : 1.2 + 0.05 * (id % 5) + 0.01 * k       # the growth ratio r_X(n) of the item
const P2_TGCOLS = ["FLI.fli", "FLI.fli1", "MEGNO.Y1", "MEGNO.Ybar1", "MEGNO.Ybar20_1", "MEGNO.Y2", "MEGNO.Ybar2", "MEGNO.Ybar20_2", "SALI.sali", "SALI.L"]
const P2_PRCOLS = ["RLI.DL", "RLI.RLI", "MEGNO.Y1", "MEGNO.Ybar1", "MEGNO.Ybar20_1", "MEGNO.Y2", "MEGNO.Ybar2", "MEGNO.Ybar20_2"]
const P2_BKROWS = ["cos(1 x[1])", "sin(1 x[1])", "cos(1 x[2])", "sin(1 x[2])"]
p2_ix(t, n) = (i = searchsortedfirst(t, n); (i <= length(t) && t[i] == n) ? i : error("step $n is not on the grid"))

function p2_le(id, t, V, lref, lns; lam5 = nothing)
    lam = fill(lref, length(t))
    for k in 0:9
        v = V[k + 1]; n, np = P2_NK[k + 1], P2_NP[k + 1]
        lp = v === :undecided ? 0.0 : 0.2
        lam[p2_ix(t, np)] = lp
        lam[p2_ix(t, n)] = v === :undecided ? 0.5 : p2_rv(v, id, k) * lp * np / n       # g(n) / g(n') = r
    end
    lam[p2_ix(t, P2_NS)] = lns
    lam[end] = lref
    return lam
end

function p2_tg(id, t, V)
    M = repeat([5.0 5.0 2.0 2.0 2.0 2.0 2.0 2.0 0.1 -1.0], length(t))
    for k in 0:9
        n, np = P2_NK[k + 1], P2_NP[k + 1]
        i, ip = p2_ix(t, n), p2_ix(t, np)
        vs, vf, vm = V[:SALI][k + 1], V[:FLI][k + 1], V[:MEGNO][k + 1]
        M[ip, 10], M[i, 10] = vs === :undecided ? (0.5, -3.0) : (-1.0, -p2_rv(vs, id, k))        # SALI.L, g = -L
        M[i, 9] = vs === :chaotic ? 1e-10 : (vs === :regular && id % 3 == 0) ? 1e-9 : 0.1          # R2 (SALI(n) <= 1e-8) differs from R1 on these regular rows
        M[ip, 1], M[i, 1] = vf === :undecided ? (0.0, 3.0) : vf === :chaotic ? (5.0, 5.0 * p2_rv(vf, id, k)) : (2.0, 2.0 * p2_rv(vf, id, k))
        M[ip, 5], M[i, 5] = vm === :undecided ? (0.0, 2.0) : vm === :chaotic ? (5.0, 5.0 * p2_rv(vm, id, k)) : (2.0 / p2_rv(vm, id, k), 2.0)
    end
    return M
end

function p2_pair(id, t, V)
    M = repeat([1.0 1e-12 2.0 2.0 2.0 2.0 2.0 2.0], length(t))
    for k in 0:9
        v = V[:RLI][k + 1]
        M[p2_ix(t, P2_NKB[k + 1]), 2] = v === :chaotic ? 1e-3 * (1 + 0.1 * (id % 5)) : 1e-12 * (1 + 0.1 * (id % 3))
    end
    return M
end

function p2_strD(id, t, rho)
    D = fill(0.3, length(t))
    for k in 0:9
        D[p2_ix(t, P2_NKB[k + 1])] = rho === :C ? 0.02 * (1 + 0.1 * (id % 4)) * (1 + 0.01 * k) : 0.5 * (1 + 0.1 * (id % 3)) * (1 + 0.01 * k)
    end
    return D
end

p2_nu0(p0) = mod(p0 / (2 * pi), 1.0)
# frequency table (E1): 10 rows x 11 columns, flags with every kind of bit
function p2_fq(id, V, p0; dshift = 0.0, odd = 0.0)
    F = zeros(10, 11); fl = zeros(Int64, 10)
    for k in 0:9
        v = V[:E1][k + 1]
        F[k + 1, 1] = F[k + 1, 2] = p2_nu0(p0) + (id in (9, 10) ? (id - 8) * k * 1e-10 : 0.0)
        F[k + 1, 3:4] .= 0.5; F[k + 1, 5:6] .= 0.01; F[k + 1, 7:8] .= 0.1; F[k + 1, 9] = 1e-3
        D = v === :chaotic ? -1.0 + 0.05 * (id % 4) : -9.0 - 0.1 * (id % 3)
        F[k + 1, 10] = D - 1.0
        F[k + 1, 11] = D + dshift + (isodd(k) ? odd : 0.0)
        fl[k + 1] = v === :undecided ? 2 + 1024 * (id % 2) : v === :regular ? (id % 3 == 0 ? 1 + 512 : 0) : (id % 2 == 0 ? 8 + 128 : 0)
        v === :undecided && (F[k + 1, 3:11] .= NaN)
    end
    return F, fl
end

function p2_pa(id, V, p0; dshift = 0.0, odd = 0.0)
    F = zeros(10, 4); fl = zeros(Int64, 10)
    for k in 0:9
        v = V[:E2][k + 1]
        F[k + 1, 1] = F[k + 1, 2] = p2_nu0(p0) + (id in (9, 10) ? (id - 8) * k * 3e-10 : 0.0)
        D = v === :chaotic ? -1.0 + 0.03 * (id % 5) : -9.0 - 0.2 * (id % 3)
        F[k + 1, 3] = D - 0.5
        F[k + 1, 4] = D + dshift + (isodd(k) ? odd : 0.0)
        fl[k + 1] = v === :undecided ? 2 : v === :regular && id % 2 == 1 ? 1 : 0
        v === :undecided && (F[k + 1, 3:4] .= NaN)
    end
    return F, fl
end

function p2_bk(id, V)
    F = zeros(40, 8); fl = zeros(Int64, 40)
    for k in 0:9
        v = V[:WBA][k + 1]
        dig = v === :chaotic ? (id % 3 == 2 ? [9.0, 3.0 + id % 3, 8.0, 7.0] : [3.0 + id % 3, 9.0, 8.0, 7.0]) : [12.0 + id % 2, 11.0, 13.0, 14.0]
        for r in 1:4
            i = 4k + r
            F[i, 1:7] .= 0.1 * r
            F[i, 8] = dig[r]
            r == 3 && v === :regular && id % 4 == 1 && (fl[i] = 4; F[i, 7] = NaN)      # reldig_undefined: reldig NaN, dig finite
            r == 2 && v === :regular && id % 3 == 0 && (fl[i] = 1)
        end
        if v === :undecided
            F[4k + 2, 3:8] .= NaN; fl[4k + 2] = 2                                        # zero: the row has no digits
        end
    end
    return F, fl
end
p2_bk_rows() = ["W=$(w) $(nm)" for w in P2_WINDOWS for nm in P2_BKROWS]

# ---------------------------------------------------------------- orbit records
function p2_info(id, h4)
    if id <= 20
        o = AN_ORBITS[id]
        return (x0 = Float64(o["x0"]), p0 = Float64(o["p0"]), ic = "$(o["x0_dec"]),$(o["p0_dec"])", params = get(o, "K_text", "2.618033988749895"),
                cls = o["class"], system = o["system"], kind = o["control_kind"], K = get(o, "K", nothing), id = id)
    end
    g = h4[id - 20]
    return (x0 = parse(Float64, g.x0), p0 = parse(Float64, g.p0), ic = "$(g.x0),$(g.p0)", params = g.K, cls = "h4", system = "standard", kind = "none", K = parse(Float64, g.K), id = id)
end

function p2_points(info, t, design)
    n = length(t)
    x = info.id <= 20 && info.kind == "fixed" ? fill(info.x0, n) : info.x0 .+ Float64.(t) .* info.p0
    p = fill(info.p0, n)
    info.id in design.drift && (p .= 2.0e6)
    info.id == 32 && (P2_N1 ÷ 2) in t && (x[searchsortedfirst(t, P2_N1 ÷ 2)] += 0.1)       # grid orbit 32: the second half rotation differs (structure axis unresolved)
    raw = hcat(x, p)
    wind = Int64.(floor.(raw ./ (2 * pi)))
    return raw .- (2 * pi) .* wind, wind
end

function p2_orbit_record(sp, info, design, p; flags = String[], extra = Dict{String, Any}())
    id, t = info.id, collect(sp.grid)
    V = design.V[id]
    lref, lns = design.lref[id], design.lns[id]
    state, wind = p2_points(info, t, design)
    at = Dict{String, Any}("system" => info.system, "params" => info.params, "ic" => info.ic, "class_label" => info.cls, "seed" => string(1000 + id), "item_hash" => "i"^64,
        "design_class" => info.cls, "flags" => String.(flags), "steps_recorded" => length(t))
    merge!(at, extra)
    lyap = Dict{String, Matrix{Float64}}()
    ind = Dict{String, NamedTuple}()
    role = sp.role
    rho = design.rho[id]
    if role === :L1
        Vd = Dict{Symbol, Vector{Symbol}}(X => copy(V[X]) for X in keys(V))
        if sp.horizon == P2_NDIR && sp.dindex > 1 && haskey(design.dir, (id, :LE, sp.dindex))
            Vd[:LE][10] = design.dir[(id, :LE, sp.dindex)]
        end
        if sp.precision == "Double64" && haskey(design.d64, (id, :LE))
            Vd[:LE][10] = design.d64[(id, :LE)]
        end
        lam = p2_le(id, t, Vd[:LE], lref, lns)
        lyap["L1"] = reshape(lam, :, 1)
    elseif role === :L5
        lam = fill(lref, length(t))
        if rho === :R || id == 32
            lam[end] = 1.4e-4; lam[p2_ix(t, P2_N1 ÷ 10)] = 1.1e-3
        end
        id == 44 && (lam[end] += 1e-6)
        lyap["L5"] = hcat(lam, -lam)
    elseif role === :tg
        Vd = Dict{Symbol, Vector{Symbol}}(X => copy(V[X]) for X in keys(V))
        for X in (:SALI, :FLI, :MEGNO)
            sp.horizon == P2_NDIR && sp.dindex > 1 && haskey(design.dir, (id, X, sp.dindex)) && (Vd[X][10] = design.dir[(id, X, sp.dindex)])
            sp.precision == "Double64" && haskey(design.d64, (id, X)) && (Vd[X][10] = design.d64[(id, X)])
        end
        ind["tg"] = AX._ind_nt(series = p2_tg(id, t, Vd), columns = P2_TGCOLS)
    elseif role === :pr
        ind["pr"] = AX._ind_nt(series = p2_pair(id, t, V), columns = P2_PRCOLS)
    elseif role === :st
        sc = zeros(18); sc[18] = 1.0
        ind["st"] = AX._ind_nt(series = reshape(p2_strD(id, t, rho), :, 1), columns = ["stretching.D"], scalars = sc, scolumns = AX.IND_ST_SCALARS, range = [-0.48, 0.48])
    elseif role in (:fq1, :fq1_raw, :fq1_qp)
        F, fl = p2_fq(id, sp.precision == "Double64" ? merge(V, Dict(:E1 => (haskey(design.d64, (id, :E1)) ? vcat(V[:E1][1:9], design.d64[(id, :E1)]) : V[:E1]))) : V, info.p0;
                      dshift = role === :fq1_raw ? 0.5 : role === :fq1_qp ? 0.0 : 0.0, odd = role === :fq1_qp ? -0.25 : 0.0)
        ind[sp.id] = AX._ind_nt(final = F, rows = ["W=$w" for w in P2_WINDOWS], fcolumns = AX.IND_FREQ_COLUMNS, flags = fl, bits = AX.IND_FREQ_BITS)
    elseif role in (:pa1, :pa_planar)
        Vp = sp.precision == "Double64" && haskey(design.d64, (id, :E2)) ? merge(V, Dict(:E2 => vcat(V[:E2][1:9], design.d64[(id, :E2)]))) : V
        F, fl = p2_pa(id, Vp, info.p0; odd = role === :pa_planar ? 0.1 : 0.0)
        ind[sp.id] = AX._ind_nt(final = F, rows = ["W=$w" for w in P2_WINDOWS], fcolumns = AX.IND_PA_COLUMNS, flags = fl, bits = AX.IND_PA_BITS,
                                scalars = [1.0, 0.0, 1.0], scolumns = AX.IND_PA_SCALARS)
    elseif role === :bk1
        F, fl = p2_bk(id, V)
        ind["bk1"] = AX._ind_nt(final = F, rows = p2_bk_rows(), fcolumns = AX.IND_BK_COLUMNS, flags = fl, bits = AX.IND_BK_BITS, scalars = [1.0, 1.0], scolumns = AX.IND_BK_SCALARS)
    end
    n = length(t)
    dn = role in (:L1, :L5) ? (role === :L1 ? ["det", "symp"] : ["det", "symp", "orth", "sum", "pair"]) : ["det", "symp"]
    diag = Dict{String, Vector{Float64}}(nm => fill(1e-16, n) for nm in dn)
    thr = Dict{String, Float64}(nm => 1e-12 for nm in dn)
    fe = Dict{String, Int}(nm => 0 for nm in dn)
    cost = Dict(sp.id => Float64[sp.horizon, 1, 0.01, 0.1, 1e5])
    return AX.OrbitRecord(id, at, Int64.(t), state, wind, lyap, diag, thr, fe, cost, "", ind)
end

# ---------------------------------------------------------------- chunk list (33 chunks)
function p2_specs(p)
    S = NamedTuple[]
    P, G = collect(1:20), collect(21:44)
    add(role, id, prec, hor, grid, orbits, arms; dindex = 1, kind = "") =
        push!(S, (; role, id, precision = prec, horizon = hor, grid, orbits, arms, dindex, direction = P2_DIRS[dindex], kind))
    add(:L1, "L1", "Float64", P2_N1, P2_GA, P, ["lyap_p"])
    add(:tg, "tg", "Float64", P2_N1, P2_GA, P, ["lyap_p"]; kind = "tangent")
    add(:L1, "L1", "Float64", P2_N1, P2_GA, G, ["lyap_g"]); add(:L5, "L5", "Float64", P2_N1, P2_GA, G, ["lyap_g"])
    add(:tg, "tg", "Float64", P2_N1, P2_GA, G, ["lyap_g"]; kind = "tangent")
    add(:st, "st", "Float64", P2_NC, P2_GB, P, ["stretch_p"]; kind = "stretching"); add(:st, "st", "Float64", P2_NC, P2_GB, G, ["stretch_g"]; kind = "stretching")
    add(:pr, "pr", "Float64", P2_NC, P2_GB, P2_SL, ["pair_l"]; kind = "pair"); add(:pr, "pr", "Float64", P2_NC, P2_GB, G, ["pair_g"]; kind = "pair")
    for (orbs, arm) in ((filter(!=(20), P), "win_p"), (G, "win_g"))
        add(:fq1, "fq1", "Float64", P2_NC, P2_GB, orbs, [arm]; kind = "frequency"); add(:bk1, "bk1", "Float64", P2_NC, P2_GB, orbs, [arm]; kind = "birkhoff")
    end
    add(:pa1, "pa1", "Float64", P2_NC + 1, P2_GB1, P, ["phase_p"]; kind = "phase_advance"); add(:pa1, "pa1", "Float64", P2_NC + 1, P2_GB1, G, ["phase_g"]; kind = "phase_advance")
    for d in 1:5; add(:L1, "L1", "Float64", P2_NDIR, P2_GAD, P2_SD, ["dir_tg"]; dindex = d); end
    for d in 1:5; add(:tg, "tg", "Float64", P2_NDIR, P2_GAD, P2_SD, ["dir_tg"]; dindex = d, kind = "tangent"); end
    add(:fq1_raw, "fq1_raw", "Float64", P2_NC, P2_GB, P2_SL, ["sens_win"]; kind = "frequency"); add(:fq1_qp, "fq1_qp", "Float64", P2_NC, P2_GB, P2_SL, ["sens_win"]; kind = "frequency")
    add(:pa_planar, "pa_planar", "Float64", P2_NC + 1, P2_GB1, P2_SL, ["sens_phase"]; kind = "phase_advance")
    add(:L1, "L1", "Double64", P2_NDIR, P2_GAD, P2_SL, ["d64_lyap"]); add(:tg, "tg", "Double64", P2_NDIR, P2_GAD, P2_SL, ["d64_lyap"]; kind = "tangent")
    add(:fq1, "fq1", "Double64", P2_NC, P2_GB, P2_SL, ["d64_win"]; kind = "frequency"); add(:bk1, "bk1", "Double64", P2_NC, P2_GB, P2_SL, ["d64_win"]; kind = "birkhoff")
    add(:pa1, "pa1", "Double64", P2_NC + 1, P2_GB1, P2_SL, ["d64_phase"]; kind = "phase_advance")
    return S
end

# bit counts per flag name: the number of orbits that carry the flag in at least one row (what the manifest's `indicator_flags` lists)
function p2_flag_counts(recs)
    out = Dict{String, Int}()
    for r in recs, (id, t) in r.indicators
        t.flags === nothing && continue
        for (j, nm) in enumerate(t.bits)
            any(w -> (w & (Int64(1) << (j - 1))) != 0, t.flags) && (out[nm] = get(out, nm, 0) + 1)
        end
    end
    return out
end

const P2_COST = Dict("tg" => (2000.0, 237.5), "pr" => (26000.0, 4100.0), "st" => (2200.0, 1702.25), "fq1" => (26800.0, 16850.5), "fq1_raw" => (26800.0, 9000.0),
                     "fq1_qp" => (26800.0, 9100.0), "pa1" => (2000.0, 1810.5), "pa_planar" => (2000.0, 1500.0), "bk1" => (2000.0, 3000.25))

# ---------------------------------------------------------------- the P1 reference study (slug p1-le-benchmark): the P1 fixture with a slug parameter
function p2_p1_study(dir::AbstractString; slug = P1_SLUG, p = an_fx_params(), lam = nothing, flags = nothing, attrs = nothing, state = nothing, cover = Dict())
    rm(dir; recursive = true, force = true)
    store = joinpath(dir, "store"); mkpath(store)
    prep(t) = replace(replace(t, "\r\n" => "\n"), r"(?m)^(\| SHA-256 of this file \(current version\) \|)[^\n|]*\|" => s"\1 SELF |")
    text = "# Preregistration: $slug\n\n## Header\n\n| Field | Value |\n|-------|-------|\n| Study slug | $slug |\n" *
           "| SHA-256 of this file (current version) | " * "0"^64 * " (rule: fixture) |\n| Amendment log | A1 (2026-01-01, fixture) |\n\nFixture. Not a study.\n"
    h = an_sha256(prep(text))
    text = replace(text, "0"^64 => h; count = 1)
    prereg = joinpath(dir, "$slug-preregistration.md"); write(prereg, text)
    panel = joinpath(dir, "$slug-panel.toml")
    cells = [Dict{String, Any}("class" => c, "drawn" => n, "k_index" => k, "qualifying" => 100, "reason" => "", "status" => "populated")
             for (k, c, n) in ((1, "C", 4), (1, "S", 2), (1, "X", 2), (1, "T", 2), (2, "C", 4), (2, "T", 2))]
    open(panel, "w") do io
        TOML.print(io, Dict{String, Any}("manifest" => Dict{String, Any}("K_binary64" => [AN_K1, AN_K2], "controls" => [17, 18, 19, 20], "cells" => cells, "format_version" => 1),
                                         "orbit" => AN_ORBITS))
    end
    ph = an_sha256(read(panel))
    params_file = an_write_params(joinpath(dir, "$slug-parameters.toml"), p)
    specs = an_specs(p)
    chunks = Dict{String, Any}[]
    for (ord, sp) in enumerate(specs)
        steps = an_grid(sp.horizon)
        recs = AX.OrbitRecord[]
        for id in sp.orbits
            prof = id <= 20 ? an_profile(id, cover) : (kind = :const, lref = id == 21 ? 0.298 : 1.172, c = Dict{String, Float64}(m => 0.0 for m in vcat(AN_V, "L2")),
                                                       q = Dict{String, Float64}(m => 1.0 for m in vcat(AN_V, "L2")))
            vals = an_vals(prof, sp.method, steps) .+ an_offset(sp)
            lam === nothing || (vals = lam(id, sp, steps, vals))
            fl = flags === nothing ? String[] : flags(id, sp)
            at = Dict{String, Any}()
            if "main" in sp.arms && sp.method in AN_V && id in (17, 18, 20)
                at["control_scope"] = id == 18 ? "info" : "bound"
                at["control_deviation"] = abs(vals[end] - Float64(an_orbit(id)["exponent"]))
            end
            attrs === nothing || merge!(at, Dict{String, Any}(attrs(id, sp)))
            push!(recs, an_record(id, sp, steps; vals = vals, flags = fl, attrs = at, state = state))
        end
        key = "method=$(sp.method) precision=$(sp.precision) horizon=$(sp.horizon) tau=$(sp.tau) eps=$(sp.eps) direction=$(sp.direction) mode=$(sp.mode)"
        ch = "c"^60 * lpad(string(ord), 4, '0')
        root = Dict{String, Any}("schema_version" => 2, "study" => slug, "chunk" => ord, "chunk_key" => key, "chunk_hash" => ch, "plan_sha256" => "p"^64,
            "panel_sha256" => ph, "prereg_sha256" => h, "git_sha" => "a"^40, "code_version" => "an-code-1", "backend_request" => "cpu", "backend" => "cpu",
            "device" => "cpu", "backend_fallback_reason" => "", "precision" => sp.precision, "method" => sp.method, "horizon" => sp.horizon, "tau" => sp.tau,
            "vectors" => sp.method in ("L1", "L2") ? 1 : 2, "coordinate_mode" => sp.mode, "direction" => sp.direction, "grid_sha256" => "g"^64,
            "status" => "ok", "attempt" => 1, "setup_check" => "not applicable", "first_attempt_reasons" => String[], "failure_reasons" => String[],
            "pilot" => false, "arms" => sp.arms, "session" => 1, "control_stop" => false, "dirty" => false, "created_utc" => "2026-10-20T00:00:00Z")
        sp.eps === nothing || (root["perturbation"] = string(sp.eps))
        AX.write_chunk(joinpath(store, "chunk-" * lpad(string(ord), 5, '0') * ".h5"), root, recs)
        push!(chunks, Dict{String, Any}("ordinal" => ord, "key" => key, "hash" => ch, "arms" => sp.arms, "status" => "ok", "attempt" => 1, "orbits" => length(sp.orbits),
            "flags" => Dict{String, Any}(f => Int[] for f in ("diagnostics exceeded", "drift", "excluded", "breakdown_tau"))))
    end
    manifest = joinpath(dir, "$slug-manifest.json")
    AX._write_json(manifest, Dict{String, Any}("study" => slug, "schema_version" => 2, "status" => "complete", "plan_sha256" => "p"^64, "chunks" => chunks,
        "prereg" => Dict("header" => h, "sha256_end" => h), "panel" => Dict("sha256_end" => ph, "build_hash" => "e"^64), "overrides" => Dict{String, Any}(),
        "reasons" => String[]))
    return (; dir, store, prereg, panel, manifest, params = params_file, p, study = slug, prereg_hash = h, panel_hash = ph, chunks = length(chunks))
end
