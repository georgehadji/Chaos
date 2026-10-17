# build_panel: preregistration in, panel file and provenance sidecar out (spec Algorithm 0 to 11).

function _manifest(path, pr, psets, rows, cells, orbs, sub, controls, root)
    inp = pr.inputs
    g = git_state(root)
    prereg_path = replace(relpath(abspath(path), root), '\\' => '/')
    drawn = sum(length(c.orbits) for cs in cells for c in cs)
    return Dict{String, Any}(
        "format_version" => FORMAT_VERSION, "package_version" => string(package_version()), "git_sha" => g.sha, "git_dirty" => g.dirty,
        "stream_rule" => STREAM_RULE, "draw_order" => collect(DRAW_ORDER), "K_binary64" => inp.K,
        "min_separation" => minimum(ps.min_separation for ps in psets),
        "prereg" => Dict{String, Any}("path" => prereg_path, "header_sha256" => pr.header_sha256,
                                       "computed_sha256" => pr.computed_sha256, "amendments" => pr.amendments),
        "inputs" => _inputs_table(inp), "newton" => _sum_stats(psets), "cells" => _cell_rows(inp, cells), "periodic" => rows,
        "counts" => _counts(psets), "ill_conditioned" => _ill_rows(psets), "controls" => [d["id"] for d in controls],
        "subsets" => Dict{String, Any}("capped_D" => !isempty(sub.removed_D), "capped_L" => !isempty(sub.removed_L),
                                       "removed_D" => sub.removed_D, "removed_L" => sub.removed_L),
        "totals" => Dict{String, Any}("orbits" => length(orbs), "drawn" => drawn, "controls" => length(controls)))
end

# Everything computed at the working precision; nothing is written here.
function _assemble(path, pr, nt, chunk)
    inp = pr.inputs
    centres = grid_centres(inp.grid_side)
    psets = [periodic_set(K; q_max = inp.q_max, N_g = inp.grid_side, rho = inp.rho, newton_max = inp.newton_max,
                          delta_dup = inp.delta_dup, bits = inp.bigfloat_bits, ntasks = nt) for K in inp.K]
    foreach(ps -> _check_separation(ps.min_separation, inp.delta_dup), psets)
    rows, pids = _periodic_rows(psets)
    cells = map(eachindex(inp.K)) do i
        flags = _screen_flags(inp.K[i], inp, psets[i], centres, chunk, nt)
        _draw_k(inp, i, psets[i], pids[i], flags, centres)
    end
    controls = _controls(inp, psets)
    orbs, sub = _number_orbits!(inp, vcat(cells...), controls)
    root = _root()
    man = _manifest(path, pr, psets, rows, cells, orbs, sub, filter(d -> d["control"], orbs), root)
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
header hash checked before any computation; nothing is written on error. Screening and the Newton search run on the CPU whatever
`backend` says (spec Compute, P6), so the panel never depends on the backend, the thread count `ntasks` or the grid `chunk`.
"""
function build_panel(prereg_path::AbstractString; out::AbstractString, chunk = nothing, ntasks = nothing, backend::Symbol = :auto)
    backend in (:auto, :gpu, :cpu) || throw(ArgumentError("backend must be :auto, :gpu or :cpu, got :$backend"))
    pr = read_prereg(prereg_path)
    nt = _ntasks(ntasks)
    panel = setprecision(BigFloat, pr.inputs.bigfloat_bits) do
        _assemble(prereg_path, pr, nt, chunk)
    end
    _write_toml(out, panel)
    _write_toml(out * ".provenance.toml", Dict{String, Any}("backend" => "cpu", "backend_fallback_reason" => "",
                                                            "backend_requested" => String(backend), "threads" => Threads.nthreads(),
                                                            "ntasks" => nt, "julia_version" => string(VERSION)))
    return out
end
