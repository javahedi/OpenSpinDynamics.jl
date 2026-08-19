using LinearAlgebra
using SparseArrays

@testset "Spin models" begin
    @testset "Single-spin fields" begin
        N = 1
        hx = [2.0]
        hz = [3.0]

        m = model(
            N,
            0.0,
            0.0;
            hx=hx,
            hz=hz,
        )

        expected = sparse([
             3.0  2.0
             2.0 -3.0
        ])

        @test size(m.hamiltonian) == (2, 2)
        @test m.hamiltonian ≈ expected
        @test ishermitian(Matrix(m.hamiltonian))
    end

    @testset "Two-spin interaction" begin
        N = 2
        Jmn = [
            0.0 1.0
            1.0 0.0
        ]

        m = model(
            N,
            1.0,
            1.0;
            Jmn=Jmn,
        )

        @test size(m.hamiltonian) == (4, 4)
        @test ishermitian(Matrix(m.hamiltonian))
    end

    @testset "Input validation" begin
        @test_throws ErrorException model(
            2,
            1.0,
            1.0;
            hx=[1.0],
        )

        @test_throws ErrorException model(
            2,
            1.0,
            1.0;
            hz=[1.0],
        )

        @test_throws ErrorException model(
            2,
            1.0,
            1.0;
            Jmn=zeros(3, 3),
        )
    end

    @testset "Stored parameters" begin
        N = 2
        hx = [0.1, 0.2]
        hz = [0.3, 0.4]
        Jmn = [
            0.0 1.0
            1.0 0.0
        ]

        m = model(
            N,
            1.0,
            2.0;
            hx=hx,
            hz=hz,
            Jmn=Jmn,
        )

        @test m.N == N
        @test m.Jxy == 1.0
        @test m.Jz == 2.0
        @test m.hx == hx
        @test m.hz == hz
        @test m.Jmn == Jmn
    end

    @testset "Update model parameters" begin
        N = 2
        m = model(N, 1.0, 1.0)

        new_hx = [0.5, 0.25]
        new_hz = [0.1, 0.2]
        new_Jmn = [
            0.0 0.7
            0.7 0.0
        ]

        update_model!(
            m,
            2.0,
            3.0;
            hx=new_hx,
            hz=new_hz,
            Jmn=new_Jmn,
        )

        @test m.Jxy == 2.0
        @test m.Jz == 3.0
        @test m.hx == new_hx
        @test m.hz == new_hz
        @test m.Jmn == new_Jmn
    end
end