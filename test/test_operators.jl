using LinearAlgebra
using SparseArrays



@testset "Pauli operators" begin
    
    @testset "Input validation" begin
        @test_throws ArgumentError spin_operators(0)
        @test_throws ArgumentError spin_operators(-1)
    end

    @testset "Dimensions" begin
        N = 3
        ops = spin_operators(N)

        @test size(ops.x[1]) == (2^N, 2^N)
        @test size(ops.z[2]) == (2^N, 2^N)
        @test size(ops.plus[3]) == (2^N, 2^N)
        @test size(ops.minus[1]) == (2^N, 2^N)
    end

    @testset "Single-spin algebra" begin
        ops = spin_operators(1)

        X = Matrix(ops.x[1])
        Z = Matrix(ops.z[1])
        P = Matrix(ops.plus[1])
        N = Matrix(ops.minus[1])
        I2 = Matrix{Float64}(I, 2, 2)

        @test X * X ≈ I2
        @test Z * Z ≈ I2
        @test P' ≈ N
        @test P * N + N * P ≈ I2
        @test Z * P - P * Z ≈ 2P
        @test Z * N - N * Z ≈ -2N
    end

    @testset "Different sites commute" begin
        ops = spin_operators(2)

        X1 = ops.x[1]
        Z2 = ops.z[2]

        @test X1 * Z2 ≈ Z2 * X1
    end
end