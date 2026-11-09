# Ensembles of independent orbits: one KernelAbstractions kernel over orbits (CPU backend in the package tests, the installed device otherwise); no method throws
# inside the kernel, a failed implicit step is recorded in the per-orbit status. The result of orbit j does not depend on the chunking.

_run_orbit(m::AbstractIntegrator, flow, z, tau, n) = (advance_n(m, flow, z, tau, n), true)
function _run_orbit(m::Gauss, flow, z, tau, n)
    for _ in 1:n
        z, st = advance_checked(m, flow, z, tau)
        st.ok || return z, false
    end
    return z, true
end

@kernel function _ensemble_kernel!(out, ok, @Const(x0), m, flow, tau, n)
    i = @index(Global)
    z, good = _run_orbit(m, flow, x0[i], tau, n)
    out[i] = z
    ok[i] = good
end

# Device hook, as `Tangent.device_backend`: the package knows the CPU only; a device extension adds a method for its device type (the GPU runs are stage 04).
ensemble_backend(::NoDevice) = CPU()
_ka_backend(b::Symbol) = b === :cpu ? CPU() : b === :gpu ? ensemble_backend(DEVICE[]) : throw(ArgumentError("backend must be :cpu or :gpu, got :$b"))

"""
    advance_ensemble(m, flow, x0, tau, n; backend = :cpu, chunk = nothing) -> (states, ok)

`n` steps of `tau` from every state of `x0` (a vector of `SVector`s), one orbit per work item; `ok[j]` is false when an implicit step of orbit `j` failed
(its state is then all NaN). `chunk` bounds the orbits per launch (default all); `backend` is `:cpu` or `:gpu` (the device installed in `Foundation.DEVICE`, by `ensemble_backend`). The result of an orbit
does not depend on `chunk`.
"""
function advance_ensemble(m::AbstractIntegrator, flow, x0::AbstractVector{S}, tau, n::Integer; backend = :cpu, chunk::Union{Nothing, Integer} = nothing) where {S <: SVector}
    _check_types(m, eltype(S), tau)
    be = _ka_backend(backend); N = length(x0)
    csz = chunk === nothing ? max(N, 1) : max(Int(chunk), 1)
    states = Vector{S}(undef, N); ok = Vector{Bool}(undef, N)
    for lo in 1:csz:N
        hi = min(lo + csz - 1, N); len = hi - lo + 1
        dx = KernelAbstractions.allocate(be, S, len); copyto!(dx, view(x0, lo:hi))
        dout = KernelAbstractions.allocate(be, S, len); dok = KernelAbstractions.allocate(be, Bool, len)
        _ensemble_kernel!(be, 16)(dout, dok, dx, m, flow, tau, n; ndrange = len)
        KernelAbstractions.synchronize(be)
        copyto!(view(states, lo:hi), dout); copyto!(view(ok, lo:hi), dok)
    end
    return states, ok
end
