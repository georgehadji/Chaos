"""
    ChaoslabCUDAExt

CUDA device for `Chaoslab.propagate` (decision D29): installs `CUDADevice` when a functional GPU is found, answers the
backend probe per element type, and reads the per-machine throughput record `data/gpu-benchmark.toml`.
"""
module ChaoslabCUDAExt

using Chaoslab, CUDA, TOML, StaticArrays
using Chaoslab.Foundation: Foundation, precision_label
using Chaoslab.Tangent: Tangent, _orbit_kernel!
using Chaoslab.Systems: StandardMap
using KernelAbstractions: KernelAbstractions

struct CUDADevice <: Foundation.AbstractDevice end

function __init__()
    CUDA.functional() && (Foundation.DEVICE[] = CUDADevice())
end

Tangent.device_backend(::CUDADevice) = CUDABackend()
Tangent.workgroup(::CUDABackend) = 256
Foundation.free_memory(::CUDADevice) = CUDA.free_memory()
Foundation.is_oom(::CUDA.OutOfGPUMemoryError) = true
Foundation.is_compile_error(::Union{CUDA.KernelError, CUDA.InvalidIRError}) = true

# Launch timeout of the display watchdog (702) and the other errors after which the context cannot be used again.
const LOST_CODES = (CUDA.ERROR_LAUNCH_TIMEOUT, CUDA.ERROR_LAUNCH_FAILED, CUDA.ERROR_ILLEGAL_ADDRESS, CUDA.ERROR_CONTEXT_IS_DESTROYED,
                    CUDA.ERROR_DEINITIALIZED, CUDA.ERROR_NO_DEVICE, CUDA.ERROR_ECC_UNCORRECTABLE, CUDA.ERROR_HARDWARE_STACK_ERROR,
                    CUDA.ERROR_ILLEGAL_INSTRUCTION, CUDA.ERROR_MISALIGNED_ADDRESS, CUDA.ERROR_INVALID_CONTEXT, CUDA.ERROR_ASSERT,
                    CUDA.ERROR_UNKNOWN)
Foundation.is_device_lost(e::CUDA.CuError) = e.code in LOST_CODES

"Path of the per-machine record (git-ignored): `[<device name>.\"<precision label>\"]` with `gpu_faster`, `cpu_ns_per_step`, ..."
record_path() = joinpath(pkgdir(Chaoslab), "data", "gpu-benchmark.toml")

# The record entry of (device, T); any read, parse or shape problem gives an empty entry (gpu_faster false, no ns value).
function record_entry(path::AbstractString, device::AbstractString, ::Type{T}) where {T}
    try
        isfile(path) || return Dict{String, Any}()
        rec = TOML.parsefile(path)[device][precision_label(T)]
        return rec isa AbstractDict ? Dict{String, Any}(rec) : Dict{String, Any}()
    catch
        return Dict{String, Any}()
    end
end

gpu_faster(device, T) = get(record_entry(record_path(), device, T), "gpu_faster", false) === true

function Foundation.gpu_step_ns(::CUDADevice, ::Type{T}) where {T}
    v = get(record_entry(record_path(), CUDA.name(CUDA.device()), T), "gpu_ns_per_step", nothing)
    return v isa Real && isfinite(v) && v > 0 ? Float64(v) : nothing
end

# One-orbit, one-step launch of the real kernel on the standard map; cached per element type.
const COMPILES = Dict{DataType, Bool}()
const COMPILES_LOCK = ReentrantLock()
function compiles(::Type{T}) where {T}
    lock(COMPILES_LOCK) do
        get!(COMPILES, T) do
            try
                b = CUDABackend()
                A(x) = KernelAbstractions.allocate(b, typeof(x), 1)
                x0, W0 = [SVector(T(1) / 3, T(0))], [SMatrix{2, 1, T}(1, 0)]
                xo, Wo, so, code, h = A(x0[1]), A(W0[1]), A(()), A(Int8(0)), A(0)
                _orbit_kernel!(b, 1)(xo, Wo, so, code, h, StandardMap(T(1)), CuArray(x0), CuArray(W0), 1, 1, (); ndrange = 1)
                KernelAbstractions.synchronize(b)
                true
            catch e
                Foundation.is_compile_error(e) || rethrow()
                false
            end
        end
    end
end

function Foundation.probe(::CUDADevice, ::Type{T}) where {T}
    name = CUDA.name(CUDA.device())
    return (loadable = true, functional = CUDA.functional(), device = name, compiles = compiles(T),
            gpu_faster = gpu_faster(name, T))
end

end
