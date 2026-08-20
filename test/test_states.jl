using LinearAlgebra

@testset "Quantum states" begin
    @testset "Polarized z state" begin
        ψ = polarized_state(3; direction=:z)

        expected = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

        @test length(ψ) == 2^3
        @test norm(ψ) ≈ 1.0
        @test Vector(ψ) ≈ expected
    end

    @testset "Neel z state" begin
        ψ = neel_state(2; direction=:z)

        # |↑↓⟩
        expected = [0.0, 1.0, 0.0, 0.0]

        @test length(ψ) == 2^2
        @test norm(ψ) ≈ 1.0
        @test Vector(ψ) ≈ expected
    end

    @testset "Polarized x state" begin
        ψ = polarized_state(2; direction=:x)


        # |+x,+x⟩
        expected = fill(0.5, 4)

        @test norm(ψ) ≈ 1.0
        @test Vector(ψ) ≈ expected
    end

    @testset "Neel x state" begin
        ψ = neel_state(2; direction=:x)

        # |+x,-x⟩
        expected = 0.5 .* [1.0, -1.0, 1.0, -1.0]

        @test norm(ψ) ≈ 1.0
        @test Vector(ψ) ≈ expected
    end

    @testset "Invalid direction" begin
        @test_throws ArgumentError neel_state(2; direction=:y)
        @test_throws ArgumentError polarized_state(2; direction=:y)
    end
end

