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





    @testset "Odd number of spins" begin
        α = 2.0
        L = 50
        N = 5

        rng1 = MersenneTwister(42)
        rng2 = MersenneTwister(42)

        coupling1 = LongRangeCouplingDisorder(
            rng1,
            α,
            L,
            N;
            reordered=true,
        )

        coupling2 = LongRangeCouplingDisorder(
            rng2,
            α,
            L,
            N;
            reordered=true,
        )

        J1 = get_matrix(coupling1)
        J2 = get_matrix(coupling2)

        @test size(J1) == (N, N)
        @test issymmetric(J1)
        @test all(diag(J1) .== 0.0)
        @test J1 == J2
    end



    @testset "Coupling input validation" begin
        @test_throws ArgumentError LongRangeCouplingDisorder(
            MersenneTwister(1),
            2.0,
            4,
            5,
        )

        @test_throws ArgumentError LongRangeCouplingDisorder(
            MersenneTwister(1),
            2.0,
            0,
            0,
        )

        @test_throws ArgumentError LongRangeCouplingDisorder(
            MersenneTwister(1),
            2.0,
            10,
            -1,
        )
    end

