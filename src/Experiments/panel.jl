# build_panel: preregistration in, panel file and provenance sidecar out (spec Algorithm 0 to 11).

"""
Path of the preregistration as stored in the manifest: relative to the repository root (with `/`), or, when the file lies outside it,
the marker `<outside repository>/` and the file name only, so no drive letter, host or user directory is recorded (D52).
"""
function _repo_path(path, root)
    rel = relpath(abspath(path), root)
    inside = !isabspath(rel) && !startswith(rel, "..")
    return inside ? replace(rel, '\\' => '/') : "<outside repository>/" * basename(path)
end

function _manifest(path, pr, psets, rows, cell_rows, orbs, sub, controls, root)
    inp = pr.inputs
    g = git_state(root)
    prereg_path = _repo_path(path, root)
    drawn = length(orbs) - length(controls)
    return Dict{String, Any}(
        "format_version" => FORMAT_VERSION, "package_version" => string(package_version()), "git_sha" => g.sha, "git_dirty" => g.dirty,
        "stream_rule" => STREAM_RULE, "stablerngs_version" => string(pkgversion(StableRNGs)),
        "julia_version" => string(VERSION), "draw_order" => collect(DRAW_ORDER), "K_binary64" => inp.K,
        "min_separation" => minimum(ps.min_separation for ps in psets),
        "prereg" => Dict{String, Any}("path" => prereg_path, "header_sha256" => pr.header_sha256,
                                       "computed_sha256" => pr.computed_sha256, "amendments" => pr.amendments),
        "inputs" => _inputs_table(inp), "newton" => _sum_stats(psets), "cells" => cell_rows, "periodic" => rows,
        "counts" => _counts(psets), "ill_conditioned" => _ill_rows(psets), "controls" => [d["id"] for d in controls],
        "subsets" => Dict{String, Any}("capped_D" => !isempty(sub.removed_D), "capped_L" => !isempty(sub.removed_L),
                                       "removed_D" => sub.removed_D, "removed_L" => sub.removed_L),
        "totals" => Dict{String, Any}("orbits" => length(orbs), "drawn" => drawn, "controls" => length(controls)))
end

# Periodic set of K index i; a refusal (E9) names the K index.
function _search_k(inp, i, nt)
    try
        ps = periodic_set(inp.K[i]; q_max = inp.q_max, N_g = inp.newton_grid_side, rho = inp.rho, newton_max = inp.newton_max,
                          delta_dup = inp.delta_dup, bits = inp.bigfloat_bits, ntasks = nt)
        _check_separation(ps.min_separation, inp.delta_dup)
        return ps
    catch e
        e isa ArgumentError && throw(ArgumentError("K index $i (K = $(inp.K[i])): " * e.msg))
        rethrow()
    end
end

# Everything computed at the working precision; nothing is written here.
function _assemble(path, pr, nt, chunk)
    inp = pr.inputs
    inp.screen_source == "random" && return _assemble_rounds(path, pr, nt, chunk)
    centres = grid_centres(inp.grid_side)
    psets = [_search_k(inp, i, nt) for i in eachindex(inp.K)]
    rows, pids = _periodic_rows(psets)
    cells = map(eachindex(inp.K)) do i
        flags = _screen_flags(inp.K[i], inp, psets[i], centres, chunk, nt)
        _draw_k(inp, i, psets[i], pids[i], flags, centres)
    end
    controls = _controls(inp, psets)
    orbs, sub = _number_orbits!(inp, vcat(cells...), controls)
    root = _root()
    man = _manifest(path, pr, psets, rows, _cell_rows(inp, cells), orbs, sub, filter(d -> d["control"], orbs), root)
    return Dict{String, Any}("manifest" => man, "orbit" => orbs)
end

function _write_toml(path, d)
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        TOML.print(io, d; sorted = true)
    end
end

"""
    build_panel(prereg_path; out, chunk = nothing, ntasks = nothing, backend = :auto) -> out

Build the study panel of the preregistration at `prereg_path` and write it to `out` (TOML: a `[manifest]` and one `[[orbit]]` per
orbit) and a provenance sidecar `out * ".provenance.toml"` with the backend and its fallback reason. The input is validated and the
header hash checked before any computation; nothing is written on error. With `screen_source = "random"` (Amendment D94) the panel is
drawn in rounds and the class bits of the screening are also written to the pool sidecar `out * ".pools.toml"`. Screening and the Newton search run on the CPU whatever
`backend` says (spec Compute, P6), so the panel never depends on the backend, the thread count `ntasks` or the grid `chunk`.
"""
function build_panel(prereg_path::AbstractString; out::AbstractString, chunk = nothing, ntasks = nothing, backend::Symbol = :auto)
    backend in (:auto, :gpu, :cpu) || throw(ArgumentError("backend must be :auto, :gpu or :cpu, got :$backend"))
    pr = read_prereg(prereg_path)
    nt = _ntasks(ntasks)
    panel = setprecision(BigFloat, pr.inputs.bigfloat_bits) do
        _assemble(prereg_path, pr, nt, chunk)
    end
    pools = pop!(panel, "pools_sidecar", nothing)           # D94: the class bits of the random screening, a verification artifact
    _write_toml(out, panel)
    pools === nothing || _write_toml(out * ".pools.toml", pools)
    _write_toml(out * ".provenance.toml", Dict{String, Any}("backend" => "cpu", "backend_fallback_reason" => "",
                                                            "backend_requested" => String(backend), "threads" => Threads.nthreads(),
                                                            "ntasks" => nt, "julia_version" => string(VERSION)))
    return out
end
