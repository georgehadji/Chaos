# Unit tests of the new Systems construct (spec P1(c), P11, RL5): the product of a verified map with itself. State layout fixed by the API proposal: the pair state is the
# concatenation (state A; state B) of two states of the base map, the dimension doubles, the one-step tangent is block diagonal, nothing couples the blocks. These tests are
# the "small spec and tests" of P11; criteria 8 and 12 rely on them end to end. Every testset starts with a line naming the new symbols it uses.

rl_randstate(::StandardMap, rng, ::Type{T}) where {T} = SVector{2, T}(T(20 * rand(rng) - 10), T(20 * rand(rng) - 10))
rl_randstate(::WrappedStandardMap, rng, ::Type{T}) where {T} = torus_state(SVector{2, T}(T(20 * rand(rng) - 10), T(20 * rand(rng) - 10)))
rl_randstate(::LinearSymplecticMap{T, D}, rng, ::Type{T}) where {T, D} = SVector{D, T}(ntuple(_ -> T(4 * rand(rng) - 2), Val(D)))

rl_alloc_forward(P, z) = (forward(P, z); @allocated forward(P, z))
rl_alloc_tangent(P, z) = (tangent(P, z); @allocated tangent(P, z))

@testset "PAIR product map protocol" begin
    _ = (PairMap,)
    # PairMap(sys) is an AbstractMap of dimension 2 d with state (z_A; z_B): forward, inverse, tangent (block diagonal [M(z_A) 0; 0 M(z_B)], exact zeros) and inverse_tangent are
    # bitwise the blocks of the base map; nothing couples the blocks (changing z_B leaves the A half of forward unchanged); the orbit of a pair is bitwise the two single orbits
    # (1000 steps); no allocation. StandardMap, WrappedStandardMap (state length 4), LinearSymplecticMap (d = 2 and d = 4), both types, 20 random pairs of states each.
    rng = StableRNG(rl_seed(14, 1))
    for T in RL_TYPES
        systems = (StandardMap(exact(T, 2.5)), WrappedStandardMap(exact(T, 2.5)), LinearSymplecticMap(SMatrix{2, 2, T}(2, 1, 1, 1)),
                   LinearSymplecticMap(block_diagonal(rotation_block(T, 0.3), hyperbolic_block(T, 0.5))))
        for sys in systems
            P = PairMap(sys)
            d = dimension(sys)
            @test P isa AbstractMap && dimension(P) == 2d
            nbad = 0
            for _ in 1:20
                za, zb, zc = rl_randstate(sys, rng, T), rl_randstate(sys, rng, T), rl_randstate(sys, rng, T)
                z = vcat(za, zb)
                tA, tB = tangent(sys, za), tangent(sys, zb)
                bd(a, b) = [Matrix(a) zeros(T, d, d); zeros(T, d, d) Matrix(b)]
                ok = isequal(forward(P, z), vcat(forward(sys, za), forward(sys, zb)))
                ok &= isequal(inverse(P, z), vcat(inverse(sys, za), inverse(sys, zb)))
                ok &= size(tangent(P, z)) == (2d, 2d) && isequal(Matrix(tangent(P, z)), bd(tA, tB))
                ok &= isequal(Matrix(inverse_tangent(P, z)), bd(inverse_tangent(sys, za), inverse_tangent(sys, zb)))
                ok &= isequal(forward(P, vcat(za, zc))[1:length(za)], forward(sys, za)) && isequal(Matrix(tangent(P, vcat(za, zc)))[1:d, :], bd(tA, tangent(sys, zc))[1:d, :])
                nbad += !ok
            end
            @test nbad == 0
            za, zb = rl_randstate(sys, rng, T), rl_randstate(sys, rng, T)
            z, a, b = vcat(za, zb), za, zb
            ok = true
            for _ in 1:1000
                z, a, b = forward(P, z), forward(sys, a), forward(sys, b)
                ok &= z == vcat(a, b)
            end
            @test ok
            @test rl_alloc_forward(P, z) == 0 && rl_alloc_tangent(P, z) == 0
        end
    end
end

@testset "PAIR block action of the tangent" begin
    _ = (PairMap,)
    # The stream forms tangent(P, z) * W. That product must be, bitwise, the product of each single tangent with its own block rows of W: for the block frame, rows 1..d of column 1 are
    # exactly tangent(sys, z_A) * xi_0 (the product the single-orbit stream forms: a d x d matrix times a d x 1 matrix), rows d+1..2d of column 2 exactly tangent(sys, z_B) * xi_0, the
    # other two blocks exactly zero; for a general W each block row of the result depends only on its own block rows of W. A dense 2d x 2d times 2d x 2 product is not enough in Double64
    # (its accumulation, fused multiply-add chains over the zero blocks, rounds differently from the d x d product and breaks the bitwise symmetry of RL5(d), I3, criterion 9(a)):
    # the tangent action of the pair must be block aware (Deviations 5). StandardMap and WrappedStandardMap, both types, 200 random pair states.
    rng = StableRNG(rl_seed(14, 3))
    for T in RL_TYPES, sys in (StandardMap(exact(T, 2.5)), WrappedStandardMap(exact(T, 2.5)))
        P = PairMap(sys)
        nbad = 0
        for j in 1:200
            za, zb = rl_randstate(sys, rng, T), rl_randstate(sys, rng, T)
            z = vcat(za, zb)
            xi = rl_xi0(rl_seed(14, 4), j, T)
            R = tangent(P, z) * rl_frame(xi)
            RA, RB = tangent(sys, za) * rl_col(xi), tangent(sys, zb) * rl_col(xi)
            ok = R[1, 1] == RA[1, 1] && R[2, 1] == RA[2, 1] && R[3, 2] == RB[1, 1] && R[4, 2] == RB[2, 1]
            ok &= R[3, 1] == 0 && R[4, 1] == 0 && R[1, 2] == 0 && R[2, 2] == 0
            Wg = SMatrix{4, 2, T}(ntuple(_ -> T(2 * rand(rng) - 1), Val(8)))
            Rg = tangent(P, z) * Wg
            for c in 1:2
                ua = tangent(sys, za) * SMatrix{2, 1, T}(Wg[1, c], Wg[2, c])
                ub = tangent(sys, zb) * SMatrix{2, 1, T}(Wg[3, c], Wg[4, c])
                ok &= Rg[1, c] == ua[1, 1] && Rg[2, c] == ua[2, 1] && Rg[3, c] == ub[1, 1] && Rg[4, c] == ub[2, 1]
            end
            nbad += !ok
        end
        @test nbad == 0
    end
end

@testset "PAIR start and distance for other maps" begin
    _ = (PairMap, pair_start, pair_distance, RLI)
    # The generic start (any map that is not a standard map: plain difference, mode "lifted"), the distance of a Y5 pair (plain Euclidean norm of the difference, no periodicity), and
    # RLI on a pair of a d = 4 map (RL1 to RL4 hold for any dimension): constant tangent gives DL = RLI = 0 exactly and s_A = s_B bitwise, with the block frame of a 4-vector.
    for T in RL_TYPES
        sys = LinearSymplecticMap(SMatrix{2, 2, T}(2, 1, 1, 1))
        u = rl_xi0(rl_seed(14, 2), 1, T)
        rec = pair_start(sys, SVector{2, T}(0.25, -0.5), u, 1e-6; orbit = 2)
        @test rec.z == vcat(SVector{2, T}(0.25, -0.5), SVector{2, T}(0.25, -0.5) + T(1e-6) * u) && rec.mode == "lifted" && rec.Dx == T(1e-6)
        @test abs(rec.Dx_r - T(1e-6)) <= 4 * rl_eps(T) && pair_distance(sys, rec.z) == rec.Dx_r
        far = vcat(SVector{2, T}(0, 0), SVector{2, T}(100, 7))
        @test abs(pair_distance(sys, far) - sqrt(T(100)^2 + T(7)^2)) <= 4 * rl_eps(T) * 100          # no minimal image on a linear map
        A4 = block_diagonal(rotation_block(T, 0.3), hyperbolic_block(T, 0.5))
        sys4 = LinearSymplecticMap(A4)
        xi4 = SVector{4, T}(T(0.5), T(0.5), T(-0.5), T(0.5))
        W4 = rl_frame(xi4)
        @test size(W4) == (8, 2)
        z4 = vcat(SVector{4, T}(0, 0, 0, 0), SVector{4, T}(T(1e-3), 0, T(2e-3), 0))
        for tau in (1, 5)
            _, _, st = propagate_orbit(PairMap(sys4), z4, W4, 40, tau, (RLI(),))
            f = reducer_final(RLI(), st[1])
            @test f.n == 40 && f.handoffs == 40 ÷ tau && f.DL == 0 && f.RLI == 0 && f.sum[1] == f.sum[2] && f.comp[1] == f.comp[2] && f.sum[1] > 0
        end
        @test rl_health_error(() -> check_frame(RLI(), W4)) === nothing
    end
end
