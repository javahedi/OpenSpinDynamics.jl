using SparseArrays

@testset "Lindblad evolution" begin

    @testset "Static system" begin
        H = spzeros(Float64, 2, 2)
        L = spzeros(Float64, 2, 2)

        ρ0 = sparse([
            1.0 0.0
            0.0 0.0
        ])

        Z = sparse([
             1.0  0.0
             0.0 -1.0
        ])

        times = [0.0, 0.1, 0.5, 1.0]
        system = LindbladSystem(H, [L])

        result = evolve(
            system,
            ρ0,
            times,
            [Z];
            method=:expm,
        )

        @test size(result) == (length(times), 1)
        @test result[:, 1] ≈ ones(length(times))
    end

    @testset "ODE and exponential methods agree" begin
        H = spzeros(Float64, 2, 2)
        L = spzeros(Float64, 2, 2)

        ρ0 = sparse([
            1.0 0.0
            0.0 0.0
        ])

        Z = sparse([
             1.0  0.0
             0.0 -1.0
        ])

        times = [0.0, 0.1, 0.2]
        system = LindbladSystem(H, [L])

        result_expm = evolve(
            system,
            ρ0,
            times,
            [Z];
            method=:expm,
        )

        result_ode = evolve(
            system,
            ρ0,
            times,
            [Z];
            method=:ode,
        )

        @test result_expm ≈ result_ode atol=1e-8
    end

    @testset "Invalid method" begin
        H = spzeros(Float64, 2, 2)
        L = spzeros(Float64, 2, 2)

        ρ0 = sparse([
            1.0 0.0
            0.0 0.0
        ])

        Z = sparse([
             1.0  0.0
             0.0 -1.0
        ])

        system = LindbladSystem(H, [L])

        @test_throws ArgumentError evolve(
            system,
            ρ0,
            [0.0, 0.1],
            [Z];
            method=:invalid,
        )
    end

    @testset "Time points are not mutated" begin
        H = spzeros(Float64, 2, 2)
        L = spzeros(Float64, 2, 2)

        ρ0 = sparse([
            1.0 0.0
            0.0 0.0
        ])

        Z = sparse([
             1.0  0.0
             0.0 -1.0
        ])

        system = LindbladSystem(H, [L])

        times = [1.0, 1.5, 2.0]
        original_times = copy(times)

        evolve(
            system,
            ρ0,
            times,
            [Z];
            method=:expm,
        )

        @test times == original_times
    end

    @testset "Time point validation" begin
        H = spzeros(Float64, 2, 2)
        L = spzeros(Float64, 2, 2)

        ρ0 = sparse([
            1.0 0.0
            0.0 0.0
        ])

        Z = sparse([
             1.0  0.0
             0.0 -1.0
        ])

        system = LindbladSystem(H, [L])

        @test_throws ArgumentError evolve(
            system,
            ρ0,
            Float64[],
            [Z],
        )

        @test_throws ArgumentError evolve(
            system,
            ρ0,
            [0.0, 1.0, 0.5],
            [Z],
        )
    end

end