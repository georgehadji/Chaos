"""
    Foundation

Value objects, errors, compensated sums, precision and backend selection (tangent-core, precision).
"""
module Foundation

using DoubleFloats: Double64
using StaticArrays
using LinearAlgebra

export NumericalHealthError, CompensatedSum, add, total, FALLBACK_REASONS, BackendChoice,
       resolve_backend, precision_label, parse_precision, exact, Double64, effective_eps

"""
Thrown when a numerical health check fails (non-finite value or residual above its threshold). `column` is the
deviation-vector column concerned, 0 when none (D34 A1a); the six-argument constructor sets it to 0. Kinds in use (`what`):
"state", "deviation matrix", "det", "symplecticity", "orthogonality", "sum", "pairing",
"orthonormalisation breakdown", "separation".
"""
struct NumericalHealthError <: Exception
    what::String
    value::Float64
    threshold::Float64
    step::Int
    orbit::Int
    backend::Symbol
    column::Int
end
NumericalHealthError(what, value, threshold, step, orbit, backend) =
    NumericalHealthError(what, value, threshold, step, orbit, backend, 0)

function Base.showerror(io::IO, e::NumericalHealthError)
    print(io, "NumericalHealthError: ", e.what, " = ", e.value, " (threshold ", e.threshold,
          ") at step ", e.step, ", orbit ", e.orbit, ", backend ", e.backend)
end

"""
    effective_eps(T) -> T

The unit roundoff the tolerances are written in: `eps(T)`, except Double64 where the rounding model of its elementary
functions gives `2^-100` (precision spec, Amendments; D34 A1b). Single home for every threshold that scales with eps.
"""
effective_eps(::Type{T}) where {T} = eps(T)
effective_eps(::Type{Double64}) = Double64(ldexp(1.0, -100))

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

"The eight recorded reasons for a CPU fallback, in `_config/architecture.md` \"Compute\" order."
const FALLBACK_REASONS = [
    "CUDA.jl not loadable", "no functional device", "kernel does not compile for the element type",
    "device out of memory", "NumericalHealthError raised only on the GPU", "reducer is CPU-only",
    "GPU not faster for the element type", "launch exceeds the display watchdog",
]

"Resolved compute backend: `backend` is `:cpu` or `:gpu`, `reason` is empty unless a GPU was passed over."
struct BackendChoice
    backend::Symbol
    device::String
    reason::String
end

"Compute device a package extension can install; `NoDevice` (no extension, or no functional GPU) is the CPU."
abstract type AbstractDevice end
struct NoDevice <: AbstractDevice end

"The device used when a run resolves to the GPU; the CUDA extension sets it on load (decision D29)."
const DEVICE = Ref{AbstractDevice}(NoDevice())

# No CUDA dependency in the library: only report whether the package could be found, never load it.
probe(::NoDevice, ::Type) = (loadable = Base.find_package("CUDA") !== nothing, functional = false,
                             device = "", compiles = false, gpu_faster = false)

"Free memory in bytes of the device a chunk runs on; the CUDA extension adds a method."
free_memory(::NoDevice) = Sys.free_memory()
"Whether an exception thrown by a device chunk is the device running out of memory."
is_oom(::Exception) = false
"Whether an exception thrown by a device chunk is the kernel failing to compile for the element type."
is_compile_error(::Exception) = false
"Whether an exception thrown by a device call means the device or its context is lost (e.g. a launch timeout of the display watchdog)."
is_device_lost(::Exception) = false
"Measured GPU time per orbit-step in nanoseconds for `T` on `device` (the per-machine record), or `nothing` when unknown."
gpu_step_ns(::AbstractDevice, ::Type) = nothing

cpu_choice(reason) = BackendChoice(:cpu, "cpu", reason)

"""
    resolve_backend(request, T; cpu_only_reducer=false, probe=nothing) -> BackendChoice

Resolve `:auto`, `:gpu` or `:cpu` to a concrete backend; order of checks: element type, CPU-only
reducer, loadable, functional, compiles, GPU faster (spec tangent-core E10). `probe` is a zero-argument function
returning the probe record; by default the record of `DEVICE[]` for `T`.
"""
function resolve_backend(request::Symbol, ::Type{T}; cpu_only_reducer::Bool = false,
                         probe = nothing) where {T}
    request in (:auto, :gpu, :cpu) || throw(ArgumentError("backend must be :auto, :gpu or :cpu, got :$request"))
    request === :cpu && return cpu_choice("")
    T === BigFloat && return cpu_choice(FALLBACK_REASONS[3])
    cpu_only_reducer && return cpu_choice(FALLBACK_REASONS[6])
    p = probe === nothing ? Foundation.probe(DEVICE[], T) : probe()
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
exact(::Type{T}, x::Union{AbstractFloat,Irrational}) where {T} = T(x)
# Integers pass through an exact BigFloat, so the only rounding is the final one into T (Double64(::Int64) rounds via Float64).
exact(::Type{T}, x::Integer) where {T} = T(BigFloat(x; precision = max(64, ndigits(x; base = 2) + 1)))
exact(::Type{T}, x::Rational) where {T} = exact(T, numerator(x)) / exact(T, denominator(x))
exact(::Type{T}, s::String) where {T} = parse(T, s)

end
