using LinearAlgebra
using Random

@testset "Couplings" begin
    @testset "Long-range clean" begin
        α = 2.0
        N = 4

        coupling = LongRangeCouplingClean(α, N)
        J = coupling.matrix

        @test size(J) == (N, N)
        @test issymmetric(J)
        @test all(diag(J) .== 0.0)

        @test J[1, 2] ≈ 1.0
        @test J[1, 3] ≈ 1 / 2^α
        @test J[1, 4] ≈ 1 / 3^α

        @test coupling.N == N
    end

    @testset "Nearest-neighbor" begin
        δ = 0.2
        N = 4

        coupling = NearestNeighborCoupling(δ, N)
        J = coupling.matrix

        @test size(J) == (N, N)
        @test issymmetric(J)
        @test all(diag(J) .== 0.0)

        @test J[1, 2] ≈ 1.0 - δ
        @test J[2, 3] ≈ 1.0 + δ
        @test J[3, 4] ≈ 1.0 - δ

        @test J[1, 3] == 0.0
        @test J[1, 4] == 0.0

        @test coupling.N == N
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

        J = coupling.matrix

        @test size(J) == (N, N)
        @test issymmetric(J)
        @test all(diag(J) .== 0.0)
        @test all(J .>= 0.0)

        @test coupling.N == N

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

        @test coupling1.matrix == coupling2.matrix
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

        J1 = coupling1.matrix
        J2 = coupling2.matrix

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


    @testset "Coupling object constructor" begin
        N = 4
        coupling = NearestNeighborCoupling(0.0, N)

        m = SpinModel(
            N;
            Jxy=1.0,
            Jz=1.0,
            coupling=coupling,
        )

        @test m.Jmn == coupling.matrix
    end


    @testset "Coupling and matrix agree" begin
        N = 4
        coupling = LongRangeCouplingClean(2.0, N)

        m_coupling = SpinModel(
            N;
            Jxy=1.0,
            Jz=0.5,
            coupling=coupling,
        )

        m_matrix = SpinModel(
            N;
            Jxy=1.0,
            Jz=0.5,
            Jmn=coupling.matrix,
        )

        @test m_coupling.hamiltonian ≈ m_matrix.hamiltonian
    end


    @testset "Coupling validation" begin
        @test_throws DimensionMismatch SpinModel(
            2;
            coupling=NearestNeighborCoupling(0.0, 3),
        )

        @test_throws ArgumentError SpinModel(
            2;
            Jmn=zeros(2, 2),
            coupling=NearestNeighborCoupling(0.0, 2),
        )
    end


    @testset "Update model coupling validation" begin
        m = SpinModel(
            4;
            Jxy=1.0,
            Jz=1.0,
            coupling=NearestNeighborCoupling(0.0, 4),
        )

        @test_throws DimensionMismatch update_model!(
            m;
            coupling=NearestNeighborCoupling(0.0, 3),
        )

        @test_throws ArgumentError update_model!(
            m;
            Jmn=zeros(4, 4),
            coupling=NearestNeighborCoupling(0.0, 4),
        )
    end

end