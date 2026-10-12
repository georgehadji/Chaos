# One-way layering of the submodules. The order is canonical in _config/architecture.md,
# section "Layering"; this list implements it. Peers in the same layer may not import each other.
const LAYERS = [["Foundation"], ["Systems"], ["Integrators"], ["Tangent"],
                ["Lyapunov", "Indicators", "Diagnostics"], ["Shadow", "TimeSeries"],
                ["Experiments"], ["Figures"]]

layer(name) = findfirst(l -> name in l, LAYERS)

"""Return "file: Owner imports Dep" for every import that is not from a strictly lower layer."""
function layering_violations(srcdir)
    bad = String[]
    for (dir, _, files) in walkdir(srcdir), f in files
        endswith(f, ".jl") || continue
        rel = relpath(dir, srcdir)
        rel == "." && continue
        owner = first(splitpath(rel))
        rank = layer(owner)
        rank === nothing && continue
        for line in eachline(joinpath(dir, f))
            for m in eachmatch(r"^\s*(?:using|import)\s+.*?(?:\.\.|Chaoslab\.)(\w+)", line)
                dep = m.captures[1]
                dr = layer(dep)
                if dr !== nothing && dep != owner && dr >= rank
                    push!(bad, "$(replace(rel, '\\' => '/'))/$f: $owner imports $dep")
                end
            end
        end
    end
    return bad
end

@test isempty(layering_violations(joinpath(@__DIR__, "..", "src")))

# The checker itself must catch an upward import.
mktempdir() do tmp
    mkpath(joinpath(tmp, "Systems"))
    write(joinpath(tmp, "Systems", "bad.jl"), "using ..Lyapunov\n")
    @test layering_violations(tmp) == ["Systems/bad.jl: Systems imports Lyapunov"]
end
