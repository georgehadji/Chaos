@testset "AC2 exact integers above 2^53" begin
    # exact(T, n) must not round an integer through Float64 (spec E2): every integer of up to 106 bits is exact in Double64
    for n in (2^53 + 1, 2^60 + 1, typemax(Int64), -(2^62 + 3), Int128(2)^100 + 1, big(2)^105 - 1)
        @test BigInt(BigFloat(exact(Double64, n))) == n
    end
    @test exact(Float64, 2^60 + 1) === Float64(2^60 + 1)   # Float64 still rounds once, to nearest
    @test exact(BigFloat, big(2)^200 + 1) == big(2)^200 + 1
    @test BigInt(BigFloat(exact(Double64, (2^60 + 1) // 1))) == 2^60 + 1   # rationals take the same exact path
end
