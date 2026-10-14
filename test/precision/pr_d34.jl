# Amendment D34 (M2 Batch B API contract, A1b), included last from runtests.jl. Existing testsets are untouched.
#   A1b  Foundation.effective_eps(::Type{T}) = eps(T), except effective_eps(Double64) = 2^-100 (precision spec, Amendments); exported.
# The thresholds that use it are tested in test/diagnostics/dg_d34.jl.

@testset "D34 effective_eps" begin
    # Float64 and BigFloat keep eps(T); Double64 uses 2^-100, sixteen times the nominal eps(Double64) = 2^-104.
    # Only values are asserted, not the return type (a Float64 2^-100 would be a legal choice for Double64).
    @test Base.isexported(Chaoslab.Foundation, :effective_eps)
    @test Base.isexported(Chaoslab, :effective_eps)
    @test effective_eps(Float64) == eps(Float64)
    @test effective_eps(Float64) == ldexp(1.0, -52)
    @test eps(Double64) == ldexp(Double64(1), -104)               # nominal spacing the spec starts from
    @test Float64(effective_eps(Double64)) == ldexp(1.0, -100)
    @test effective_eps(Double64) == ldexp(Double64(1), -100)
    @test effective_eps(Double64) == 16 * eps(Double64)
    for bits in (128, 256)
        setprecision(BigFloat, bits) do
            @test effective_eps(BigFloat) == eps(BigFloat)
            @test effective_eps(BigFloat) == ldexp(BigFloat(1), 1 - bits)
        end
    end
    # a threshold built with it is strictly looser than the nominal one for Double64 and equal for Float64
    @test effective_eps(Double64) > eps(Double64)
    @test !(effective_eps(Float64) > eps(Float64))
end
