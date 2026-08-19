using LinearAlgebra
using Random

@testset "Couplings" begin
    @testset "Long-range clean" begin
        α = 2.0
        N = 4

        coupling = LongRangeCouplingClean(α, N)
        J = get_matrix(coupling)

        @test size(J) == (N, N)
        @test issymmetric(J)
        @test all(diag(J) .== 0.0)

        @test J[1, 2] ≈ 1.0
        @test J[1, 3] ≈ 1 / 2^α
        @test J[1, 4] ≈ 1 / 3^α

        @test get_N(coupling) == N
    end

    @testset "Nearest-neighbor" begin
        δ = 0.2
        N = 4

        coupling = NearestNeighborCoupling(δ, N)
        J = get_matrix(coupling)

        @test size(J) == (N, N)
        @test issymmetric(J)
        @test all(diag(J) .== 0.0)

        @test J[1, 2] ≈ 1.0 - δ
        @test J[2, 3] ≈ 1.0 + δ
        @test J[3, 4] ≈ 1.0 - δ

        @test J[1, 3] == 0.0
        @test J[1, 4] == 0.0

        @test get_N(coupling) == N
    end

    @testset "Long-range disorder" begin
        α = 2.0
        L = 100
        N = 10

        coupling = LongRangeCouplingDisorder(
            α,
            L,
            N;
            reordered=false,
        )

        J = get_matrix(coupling)

        @test size(J) == (N, N)
        @test issymmetric(J)
        @test all(diag(J) .== 0.0)
        @test all(J .>= 0.0)

        @test get_N(coupling) == N



        rng1 = MersenneTwister(1234)
        rng2 = MersenneTwister(1234)

        coupling1 = LongRangeCouplingDisorder(
            rng1,
            α,
            L,
            N;
            reordered=false,
        )

        coupling2 = LongRangeCouplingDisorder(
            rng2,
            α,
            L,
            N;
            reordered=false,
        )

        @test get_matrix(coupling1) == get_matrix(coupling2)
    end
end