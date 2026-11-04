# Dispatch on the study slug (Amendment D89, E21): `analyze` reads the single header row "Study slug" of the preregistration and runs the pipeline of that study. An unknown slug,
# a manifest of another study, the parameter table of the other study and the P1 reference options with a P1 preregistration are refused before anything is written.

"""
    study_slug(prereg) -> String

The value of the single header row "Study slug" of the preregistration (E21): lowercase letters, digits and hyphens. An `ArgumentError` naming the slug otherwise.
"""
function study_slug(prereg::AbstractString)
    isfile(prereg) || _fail("preregistration $prereg does not exist")
    rows = [l for l in split(replace(read(prereg, String), "\r\n" => "\n"), '\n') if startswith(l, "| Study slug |")]
    length(rows) == 1 || _fail("preregistration $prereg has $(length(rows)) header rows \"Study slug\", expected exactly one")
    slug = String(strip(split(rows[1], '|')[3]))
    occursin(r"^[a-z0-9-]+$", slug) || _fail("the study slug \"$slug\" of $prereg is not made of lowercase letters, digits and hyphens")
    return slug
end

# The `study` value of the manifest must equal the slug (E21).
function _check_study(manifest, slug)
    isfile(manifest) || return
    man = try
        read_json(manifest)
    catch
        return
    end
    st = get(man, "study", nothing)
    (st === nothing || st == slug) || _fail("manifest study \"$st\" differs from the preregistration slug \"$slug\"")
end

# The parameter table of study p2-indicators given with a P1 preregistration: its missing and unknown keys are named (E21).
function _check_other_table(params)
    isfile(params) || return
    d = try
        TOML.parsefile(params)
    catch
        return
    end
    unknown = sort!([k for k in keys(d) if !(Symbol(k) in PARAM_KEYS)])
    any(k -> Symbol(k) in P2_PARAM_KEYS, unknown) || return
    missing_ = [string(k) for k in PARAM_KEYS if !haskey(d, string(k))]
    _fail("parameter file $params is not a table of study $P1_SLUG: unknown key $(join(unknown, ", ")); missing key $(join(missing_, ", "))")
end

"""
    analyze(; prereg, manifest, panel, params, store, out, figures = false, threads, seed, p1_prereg, p1_manifest, p1_params, p1_store) -> NamedTuple

Run the analysis of the study named by the preregistration (E21): `p1-le-benchmark` runs `analyze_p1` unchanged, `p2-indicators` runs `analyze_p2` with the P1
reference inputs `p1_*` (default: the files of study `p1-le-benchmark` by the study layout). The `p1_*` keywords are an error for `p1-le-benchmark`.
"""
function analyze(; prereg, manifest, panel, params, store, out, figures::Bool = false, threads::Integer = Threads.nthreads(), seed::Union{Nothing, Integer} = nothing,
                 p1_prereg = nothing, p1_manifest = nothing, p1_params = nothing, p1_store = nothing, kw...)
    slug = study_slug(prereg)
    slug == P3_SLUG && (figures || seed !== nothing) && _fail("the keywords figures and seed are not options of the test side of study $P3_SLUG")
    slug == P3_SLUG && return analyze_p3t(; prereg, manifest, panel, params, store, out, threads, p1_prereg, p1_params, kw...)
    kwerr() = _fail("unknown keyword $(first(keys(kw))): the keywords of the test side are for study $P3_SLUG only")
    if slug == P2_SLUG
        isempty(kw) || kwerr()
        return analyze_p2(; prereg, manifest, panel, params, store, out, figures, threads, p1_prereg, p1_manifest, p1_params, p1_store)
    elseif slug == P1_SLUG
        isempty(kw) || kwerr()
        any(!isnothing, (p1_prereg, p1_manifest, p1_params, p1_store)) &&
            _fail("the keywords p1_prereg, p1_manifest, p1_params and p1_store are for study $P2_SLUG; this preregistration is of study $slug")
        _check_study(manifest, slug)
        _check_other_table(params)
        return analyze_p1(; prereg, manifest, panel, params, store, out, figures, threads, seed)
    end
    _fail("unknown study slug \"$slug\": the known studies are $P1_SLUG and $P2_SLUG")
end
