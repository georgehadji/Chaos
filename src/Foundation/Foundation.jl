"""
    Foundation

Value objects, errors, compensated sums, precision and backend selection (tangent-core, precision).
"""
module Foundation

using DoubleFloats: Double64
using StaticArrays
using LinearAlgebra

export NumericalHealthError, CompensatedSum, add, total, FALLBACK_REASONS, BackendChoice,
       resolve_backend, precision_label, parse_precision, exact, Double64

"Thrown when a numerical health check fails (non-finite value or residual above its threshold)."
struct NumericalHealthError <: Exception
    what::String
    value::Float64
    threshold::Float64
    step::Int
    orbit::Int
    backend::Symbol
end

function Base.showerror(io::IO, e::NumericalHealthError)
    print(io, "NumericalHealthError: ", e.what, " = ", e.value, " (threshold ", e.threshold,
          ") at step ", e.step, ", orbit ", e.orbit, ", backend ", e.backend)
end

"Neumaier compensated running sum: `s` is the sum, `c` the accumulated rounding correction."
struct CompensatedSum{T}
    s::T
    c::T
end
CompensatedSum{T}() where {T} = CompensatedSum{T}(zero(T), zero(T))

"""
Add `x` to a [`CompensatedSum`](@ref) and return the new accumulator (Neumaier update); `x` must be an
`Integer` or a float that widens exactly into `T`, anything else throws `ArgumentError`.
"""
function add(acc::CompensatedSum{T}, x) where {T}
    x isa T || x isa Integer || (x isa AbstractFloat && promote_type(T, typeof(x)) === T) ||
        throw(ArgumentError("add: x of type $(typeof(x)) would be narrowed to $T"))
    y = T(x)
    t = acc.s + y
    c = abs(acc.s) >= abs(y) ? (acc.s - t) + y : (y - t) + acc.s
    return CompensatedSum{T}(t, acc.c + c)
end

"Value of a [`CompensatedSum`](@ref): `s + c`."
total(acc::CompensatedSum) = acc.s + acc.c

"The seven recorded reasons for a CPU fallback, in `_config/architecture.md` \"Compute\" order."
const FALLBACK_REASONS = [
    "CUDA.jl not loadable", "no functional device", "kernel does not compile for the element type",
    "device out of memory", "NumericalHealthError raised only on the GPU", "reducer is CPU-only",
    "GPU not faster for the element type",
]

"Resolved compute backend: `backend` is `:cpu` or `:gpu`, `reason` is empty unless a GPU was passed over."
struct BackendChoice
    backend::Symbol
    device::String
    reason::String
end

# No CUDA dependency in the library: only report whether the package could be found, never load it.
default_probe() = (loadable = Base.find_package("CUDA") !== nothing, functional = false,
                   device = "", compiles = false, gpu_faster = false)

cpu_choice(reason) = BackendChoice(:cpu, "cpu", reason)

"""
    resolve_backend(request, T; cpu_only_reducer=false, probe=default_probe) -> BackendChoice

Resolve `:auto`, `:gpu` or `:cpu` to a concrete backend; order of checks: element type, CPU-only
reducer, loadable, functional, compiles, GPU faster (spec tangent-core E10).
"""
function resolve_backend(request::Symbol, ::Type{T}; cpu_only_reducer::Bool = false,
                         probe = default_probe) where {T}
    request in (:auto, :gpu, :cpu) || throw(ArgumentError("backend must be :auto, :gpu or :cpu, got :$request"))
    request === :cpu && return cpu_choice("")
    T === BigFloat && return cpu_choice(FALLBACK_REASONS[3])
    cpu_only_reducer && return cpu_choice(FALLBACK_REASONS[6])
    p = probe()
    p.loadable || return cpu_choice(FALLBACK_REASONS[1])
    p.functional || return cpu_choice(FALLBACK_REASONS[2])
    p.compiles || return cpu_choice(FALLBACK_REASONS[3])
    request === :auto && !p.gpu_faster && return cpu_choice(FALLBACK_REASONS[7])
    return BackendChoice(:gpu, p.device, "")
end

"Precision string recorded with a run: `Float64`, `Double64` or `BigFloat(bits)` (bits in force at the call)."
precision_label(::Type{Float64}) = "Float64"
precision_label(::Type{Double64}) = "Double64"
precision_label(::Type{BigFloat}) = "BigFloat($(precision(BigFloat)))"

"Inverse of [`precision_label`](@ref): returns `(T, bits)` with `bits === nothing` unless `T` is `BigFloat`."
function parse_precision(s::AbstractString)
    s == "Float64" && return (Float64, nothing)
    s == "Double64" && return (Double64, nothing)
    m = match(r"^BigFloat\((\d+)\)\z", s)
    m === nothing && throw(ArgumentError("unknown precision string \"$s\""))
    bits = tryparse(Int, m[1])
    (bits === nothing || bits < 64) && throw(ArgumentError("BigFloat bits must be an integer >= 64 in \"$s\""))
    return (BigFloat, bits)
end

# Build a value in T from its definition, never through a Float64 detour (spec precision E2).
"""
    exact(T, x) -> T

Convert `x` into `T` without a Float64 detour: integers, floats (exact widening), rationals, irrationals
and decimal strings (parsed in `T`).
"""
exact(::Type{T}, x::Union{Integer,AbstractFloat,Irrational}) where {T} = T(x)
exact(::Type{T}, x::Rational) where {T} = T(numerator(x)) / T(denominator(x))
exact(::Type{T}, s::String) where {T} = parse(T, s)

end
