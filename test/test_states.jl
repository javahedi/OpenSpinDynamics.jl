using LinearAlgebra

@testset "Quantum states" begin
    @testset "Polarized z state" begin
        ψ = construct_state(PolarizedState(3, "z"))

        expected = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

        @test length(ψ) == 2^3
        @test norm(ψ) ≈ 1.0
        @test Vector(ψ) ≈ expected
    end

    @testset "Neel z state" begin
        ψ = construct_state(NeelState(2, "z"))

        # |↑↓⟩
        expected = [0.0, 1.0, 0.0, 0.0]

        @test length(ψ) == 2^2
        @test norm(ψ) ≈ 1.0
        @test Vector(ψ) ≈ expected
    end

    @testset "Polarized x state" begin
        ψ = construct_state(PolarizedState(2, "x"))

        # |+x,+x⟩
        expected = fill(0.5, 4)

        @test norm(ψ) ≈ 1.0
        @test Vector(ψ) ≈ expected
    end

    @testset "Neel x state" begin
        ψ = construct_state(NeelState(2, "x"))

        # |+x,-x⟩
        expected = 0.5 .* [1.0, -1.0, 1.0, -1.0]

        @test norm(ψ) ≈ 1.0
        @test Vector(ψ) ≈ expected
    end

    @testset "Invalid direction" begin
        @test_throws ErrorException construct_state(NeelState(2, "y"))
        @test_throws ErrorException construct_state(PolarizedState(2, "y"))
    end
end