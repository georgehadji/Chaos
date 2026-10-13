# Test-only maps and reducer for the diagnostics suite. Loaded lazily by dg_load_faults()
# so that a missing Chaoslab.AbstractMap errors the testsets that need it, not the whole file.

# Constant-tangent map with a fixed state: tangent is `A`, or `bad` at states whose first
# coordinate equals t1 or t2 (NaN never matches). Used to inject an exact residual.
struct dg_FaultMap{T,D,L} <: Chaoslab.AbstractMap
    A::SMatrix{D,D,T,L}
    bad::SMatrix{D,D,T,L}
    t1::T
    t2::T
end
dg_FaultMap(A::SMatrix{D,D,T,L}) where {D,T,L} = dg_FaultMap{T,D,L}(A, A, T(NaN), T(NaN))
dg_FaultMap(A::SMatrix{D,D,T,L}, bad::SMatrix{D,D,T,L}, t1, t2) where {D,T,L} =
    dg_FaultMap{T,D,L}(A, bad, T(t1), T(t2))

Chaoslab.forward(::dg_FaultMap, z) = z
Chaoslab.tangent(m::dg_FaultMap, z) = (z[1] == m.t1 || z[1] == m.t2) ? m.bad : m.A
Chaoslab.dimension(::dg_FaultMap{T,D}) where {T,D} = D

# Reducer that rescales W to unit Frobenius norm (k = 1 growth control, E7 replacement).
struct dg_Renorm end
Chaoslab.reducer_init(::dg_Renorm, ::Type, ::Val, ::Val) = nothing
Chaoslab.reducer_update(::dg_Renorm, st, n::Int, x::SVector, W::SMatrix) = (st, W / sqrt(sum(abs2, W)))
Chaoslab.reducer_final(::dg_Renorm, st) = (;)
