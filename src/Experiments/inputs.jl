# Panel inputs (E0): parse closed-form strings at BigFloat precision, validate before any computation (acceptance criterion 1).

const _REQUIRED = ("K", "n_c", "grid_side", "N_scr", "alpha", "theta_T", "g_T", "r_T", "theta_C", "r_C", "theta_S", "w", "f_S",
                   "L_S", "log10_r", "q_max", "rho", "delta_R", "newton_max", "delta_dup", "n_ctrl_max", "mu", "s_D", "cap_D",
                   "s_L", "cap_L")
const _OPTIONAL = ("tau", "seed", "bigfloat_bits")
const Q_MAX_SUPPORTED = 12          # Newton search range: the conditioning of a period-q solve grows like exp(q lambda)

_bad(keys, msg) = throw(ArgumentError("panel input " * join(("`$k`" for k in keys), " / ") * ": " * msg))

const _FUNCS = Dict{Symbol, Any}(:+ => +, :- => -, :* => *, :/ => /, :^ => ^, :sqrt => sqrt, :exp => exp, :log => log,
                                 :sin => sin, :cos => cos)

# Evaluate a parsed expression at the ambient BigFloat precision.
function _ev(x)
    x isa Number && return BigFloat(x)
    x === :pi && return BigFloat(pi)
    x === :NaN && return BigFloat(NaN)
    x === :Inf && return BigFloat(Inf)
    if x isa Expr && x.head === :call && haskey(_FUNCS, x.args[1])
        return _FUNCS[x.args[1]](map(_ev, x.args[2:end])...)
    elseif x isa Expr && x.head === :block
        return _ev(last(filter(a -> !(a isa LineNumberNode), x.args)))
    end
    throw(ArgumentError("unsupported expression $x"))
end

"Value of a decimal or closed-form entry (`\"pi/2\"`, `\"(3 + sqrt(5))/2\"`, `0.5`) at the ambient BigFloat precision."
function _real(key, v)
    v isa Bool && _bad((key,), "must be a number or a closed-form string")
    v isa Real && return BigFloat(v)
    v isa AbstractString || _bad((key,), "must be a number or a closed-form string, got $(typeof(v))")
    plain = tryparse(BigFloat, String(v))
    plain === nothing || return plain                       # a decimal: the exact decimal, never a Float64 first
    try
        return _ev(Meta.parse(replace(String(v), r"(\d)pi" => s"\1*pi")))
    catch
        _bad((key,), "cannot evaluate \"$v\"")
    end
end

_f64(key, v, bits) = setprecision(() -> Float64(_real(key, v)), BigFloat, bits)

function _finite(key, v, bits)
    x = _f64(key, v, bits)
    isfinite(x) || _bad((key,), "must be finite, got $x")
    return x
end

function _int(key, v; min = 1)
    (v isa Integer && !(v isa Bool)) || _bad((key,), "must be an integer, got $(repr(v))")
    v >= min || _bad((key,), "must be at least $min, got $v")
    return Int(v)
end

function _positive(key, v, bits; upper = Inf)
    x = _finite(key, v, bits)
    (x > 0 && x < upper) || _bad((key,), "must lie in (0, $upper), got $x")
    return x
end

function _defaults()
    path = joinpath(_root(), "_config", "numerics-defaults.toml")
    return TOML.parsefile(path)
end

function _check_keys(d)
    for k in keys(d)
        k in _REQUIRED || k in _OPTIONAL || _bad((k,), "unknown key")
    end
    for k in _REQUIRED
        haskey(d, k) || _bad((k,), "missing key")
    end
end

function _check_integers(d, tau)
    n_c, ng, N, w = _int("n_c", d["n_c"]), _int("grid_side", d["grid_side"]), _int("N_scr", d["N_scr"]), _int("w", d["w"])
    qm = _int("q_max", d["q_max"])
    qm <= Q_MAX_SUPPORTED || _bad(("q_max",), "beyond the supported range 1 to $Q_MAX_SUPPORTED, got $qm")
    nm = _int("newton_max", d["newton_max"])
    w % tau == 0 || _bad(("w", "tau"), "w = $w is not a multiple of tau = $tau")
    N % w == 0 || _bad(("N_scr", "w"), "N_scr = $N is not a multiple of w = $w")
    LS = _int("L_S", d["L_S"])
    (LS >= w && LS <= N) || _bad(("L_S",), "must lie in [w, N_scr] = [$w, $N], got $LS")
    return (n_c = n_c, grid_side = ng, N_scr = N, w = w, q_max = qm, newton_max = nm, L_S = LS)
end

function _check_subsets(d)
    sD, sL = _int("s_D", d["s_D"]), _int("s_L", d["s_L"])
    cD, cL = _int("cap_D", d["cap_D"]), _int("cap_L", d["cap_L"])
    nc = _int("n_ctrl_max", d["n_ctrl_max"]; min = 0)
    sL <= sD || _bad(("s_L", "s_D"), "s_L = $sL exceeds s_D = $sD")
    cD >= nc + 1 || _bad(("cap_D", "n_ctrl_max"), "cap_D = $cD is below the number of controls $(nc + 1)")
    return (s_D = sD, s_L = sL, cap_D = cD, cap_L = cL, n_ctrl_max = nc)
end

function _check_thresholds(d, bits)
    thT, thS, thC = (_finite(k, d[k], bits) for k in ("theta_T", "theta_S", "theta_C"))
    thT < thS || _bad(("theta_T", "theta_S"), "theta_T = $thT must be below theta_S = $thS")
    thT < thC || _bad(("theta_T", "theta_C"), "theta_T = $thT must be below theta_C = $thC")
    fS = _finite("f_S", d["f_S"], bits)
    (0 < fS < 1) || _bad(("f_S",), "must lie in (0, 1), got $fS")
    rho = _positive("rho", d["rho"], bits; upper = 1.0)
    pos = (g_T = _positive("g_T", d["g_T"], bits), r_T = _positive("r_T", d["r_T"], bits), r_C = _positive("r_C", d["r_C"], bits),
           delta_R = _positive("delta_R", d["delta_R"], bits), delta_dup = _positive("delta_dup", d["delta_dup"], bits))
    return (; theta_T = thT, theta_S = thS, theta_C = thC, f_S = fS, rho, pos...)
end

function _check_special(d, bits)
    lr = d["log10_r"]
    (lr isa AbstractVector && length(lr) == 2) || _bad(("log10_r",), "must be a pair [l_lo, l_hi]")
    lo, hi = _finite("log10_r", lr[1], bits), _finite("log10_r", lr[2], bits)
    lo < hi || _bad(("log10_r",), "needs l_lo < l_hi, got [$lo, $hi]")
    10.0^hi < pi || _bad(("log10_r",), "10^l_hi = $(10.0^hi) must stay below pi")
    mu = _finite("mu", d["mu"], bits)
    mu > 1 || _bad(("mu",), "must exceed 1, got $mu")
    alpha = _finite("alpha", d["alpha"], bits)
    Ks = d["K"]
    (Ks isa AbstractVector && !isempty(Ks)) || _bad(("K",), "must be a non-empty list")
    K = [_finite("K", k, bits) for k in Ks]
    allunique(K) || _bad(("K",), "has a repeated value (compared as binary64)")
    return (K = K, log10_r = [lo, hi], mu = mu, alpha = alpha)
end

"""
    panel_inputs(dict)

Validate the E0 inputs (string keys, as in the TOML panel block) and return them as a NamedTuple with the key names as
fields, closed-form entries evaluated in BigFloat and rounded once to binary64. Optional `tau`, `seed`, `bigfloat_bits`
default to `_config/numerics-defaults.toml`. The raw values are kept in `raw`, the origin of `tau`, `seed` and
`bigfloat_bits` in `sources`. A bad input is an `ArgumentError` naming the key.
"""
function panel_inputs(d::AbstractDict)
    _check_keys(d)
    df = _defaults()
    pick(k, dk) = haskey(d, k) ? d[k] : df[dk]
    bits = _int("bigfloat_bits", pick("bigfloat_bits", "bigfloat_bits"); min = 64)
    tau = _int("tau", pick("tau", "renorm_interval_tau"))
    seed = _int("seed", pick("seed", "seed"); min = 0)
    ints = _check_integers(d, tau)
    sub = _check_subsets(d)
    th = _check_thresholds(d, bits)
    sp = _check_special(d, bits)
    src(k) = haskey(d, k) ? "preregistration" : "default"
    sources = Dict("tau" => src("tau"), "seed" => src("seed"), "bigfloat_bits" => src("bigfloat_bits"))
    return (; sp.K, ints.n_c, ints.grid_side, ints.N_scr, sp.alpha, th.theta_T, th.g_T, th.r_T, th.theta_C, th.r_C, th.theta_S, ints.w,
            th.f_S, ints.L_S, sp.log10_r, ints.q_max, th.rho, th.delta_R, ints.newton_max, th.delta_dup, sub.n_ctrl_max, sp.mu,
            sub.s_D, sub.cap_D, sub.s_L, sub.cap_L, tau, seed, bigfloat_bits = bits, raw = Dict{String, Any}(d), sources)
end
