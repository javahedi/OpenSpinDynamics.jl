using LinearAlgebra
using SparseArrays



@testset "Pauli operators" begin
    
    @testset "Input validation" begin
        @test_throws ArgumentError generate_operators(0)
        @test_throws ArgumentError generate_operators(-1)
    end

    @testset "Dimensions" begin
        N = 3
        paulis = generate_operators(N)

        @test size(paulis.Xop[1]) == (2^N, 2^N)
        @test size(paulis.Zop[2]) == (2^N, 2^N)
        @test size(paulis.Pop[3]) == (2^N, 2^N)
        @test size(paulis.Nop[1]) == (2^N, 2^N)
    end

    @testset "Single-spin algebra" begin
        paulis = generate_operators(1)

        X = Matrix(paulis.Xop[1])
        Z = Matrix(paulis.Zop[1])
        P = Matrix(paulis.Pop[1])
        N = Matrix(paulis.Nop[1])
        I2 = Matrix{Float64}(I, 2, 2)

        @test X * X ≈ I2
        @test Z * Z ≈ I2
        @test P' ≈ N
        @test P * N + N * P ≈ I2
        @test Z * P - P * Z ≈ 2P
        @test Z * N - N * Z ≈ -2N
    end

    @testset "Different sites commute" begin
        paulis = generate_operators(2)

        X1 = paulis.Xop[1]
        Z2 = paulis.Zop[2]

        @test X1 * Z2 ≈ Z2 * X1
    end
end