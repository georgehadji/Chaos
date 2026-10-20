# The parameter table E0: a flat TOML file that transcribes the frozen values of the preregistration (spec open point P1).

const PARAM_INT = (:N_1, :N_long, :N_B, :n_H2, :bootstrap_resamples, :h4_T, :h4_T_cross, :min_set)
const PARAM_REAL = (:h1a_bound, :h1a_fraction, :h1b_bound, :h1b_fraction, :h1c_bound, :h1c_fraction, :h1c_eps, :control_bound, :control_bound_L9,
                    :h2_median, :h2_interval, :bootstrap_level, :h3_band, :h3_ratio, :h3_interval, :h4_ref_halfpi, :h4_tol_halfpi, :h4_ref_2pi,
                    :h4_tol_2pi, :h5a_min, :h5b_max, :h5c_max, :lambda_threshold, :growth_ratio, :rotation_tol, :return_tol, :ladder_tol,
                    :drift_bound, :K_c)
const PARAM_KEYS = (:N_1, :N_long, :N_B, :n_H2, :h1a_bound, :h1a_fraction, :h1b_bound, :h1b_fraction, :h1c_bound, :h1c_fraction, :h1c_eps,
                    :control_bound, :control_bound_L9, :h2_median, :h2_interval, :bootstrap_resamples, :bootstrap_level, :h3_band, :h3_ratio,
                    :h3_interval, :h4_ref_halfpi, :h4_tol_halfpi, :h4_ref_2pi, :h4_tol_2pi, :h4_T, :h4_T_cross, :h5a_min, :h5b_max, :h5c_max,
                    :lambda_threshold, :growth_ratio, :rotation_tol, :return_tol, :ladder_tol, :drift_bound, :min_set, :K_c)

function _param_value(path, k::Symbol, v)
    if k in PARAM_INT
        v isa Integer && !(v isa Bool) || throw(ArgumentError("parameter file $path: key $k needs an integer, got $(repr(v))"))
        return Int(v)
    end
    v isa Real && !(v isa Bool) || throw(ArgumentError("parameter file $path: key $k needs a number, got $(repr(v))"))
    return Float64(v)
end

"""
    read_params(path) -> NamedTuple

Read the parameter table. An `ArgumentError` names the path for a missing or invalid file and the key for a missing, unknown or mistyped one.
"""
function read_params(path::AbstractString)
    isfile(path) || throw(ArgumentError("parameter file $path does not exist"))
    d = try
        TOML.parsefile(path)
    catch e
        throw(ArgumentError("parameter file $path is not valid TOML: $(first(sprint(showerror, e), 200))"))
    end
    for k in keys(d)
        Symbol(k) in PARAM_KEYS || throw(ArgumentError("parameter file $path: unknown key $k"))
    end
    for k in PARAM_KEYS
        haskey(d, string(k)) || throw(ArgumentError("parameter file $path: missing key $k"))
    end
    p = NamedTuple{PARAM_KEYS}(Tuple(_param_value(path, k, d[string(k)]) for k in PARAM_KEYS))
    _check_ranges(path, p)
    return p
end

# Values must be finite, counts and horizons at least 1, the level strictly between 0 and 1.
function _check_ranges(path, p)
    for k in PARAM_KEYS
        v = getfield(p, k)
        isfinite(v) || throw(ArgumentError("parameter file $path: key $k is not finite"))
        k in PARAM_INT && v < 1 && throw(ArgumentError("parameter file $path: key $k must be at least 1, got $v"))
    end
    0 < p.bootstrap_level < 1 || throw(ArgumentError("parameter file $path: key bootstrap_level must lie strictly between 0 and 1, got $(p.bootstrap_level)"))
end
