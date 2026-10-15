"""
    Chaoslab

Lyapunov exponents, chaos indicators and symplectic integration for Hamiltonian maps and flows.
Architecture, module map and layering: `_config/architecture.md`. Each submodule lives in
`src/<Module>/` and is added by a module-build run.
"""
module Chaoslab

include("Foundation/Foundation.jl")
include("Systems/Systems.jl")
include("Tangent/Tangent.jl")
include("Diagnostics/Diagnostics.jl")
include("Lyapunov/Lyapunov.jl")

using .Foundation, .Systems, .Tangent, .Diagnostics, .Lyapunov

# Re-export every name a submodule exports, so users and tests need only `using Chaoslab`.
for m in (Foundation, Systems, Tangent, Diagnostics, Lyapunov), n in names(m)
    n === nameof(m) || @eval export $n
end
include("precompile.jl")
end
