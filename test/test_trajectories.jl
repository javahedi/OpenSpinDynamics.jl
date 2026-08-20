using SparseArrays
using Random

@testset "Stochastic trajectories" begin

    @testset "Static system" begin
        H = spzeros(Float64, 2, 2)
        L = spzeros(Float64, 2, 2)

        ψ0 = sparsevec([1], [1.0], 2)

        Z = sparse([
             1.0  0.0
             0.0 -1.0
        ])

        times = [0.0, 0.1, 0.5, 1.0]
        system = TrajectorySystem(H, [L])

        mean_values, std_values = evolve(
            system,
            ψ0,
            times,
            [Z];
            num_samples=4,
        )

        @test size(mean_values) == (length(times), 1)
        @test size(std_values) == (length(times), 1)

        @test mean_values[:, 1] ≈ ones(length(times))
        @test std_values[:, 1] ≈ zeros(length(times))
    end

    @testset "Single sample standard deviation" begin
        H = spzeros(Float64, 2, 2)
        L = spzeros(Float64, 2, 2)

        ψ0 = sparsevec([1], [1.0], 2)
        Z = sparse([1.0 0.0; 0.0 -1.0])

        system = TrajectorySystem(H, [L])

        mean_values, std_values = evolve(
            system,
            ψ0,
            [0.0, 0.1],
            [Z];
            num_samples=1,
        )

        @test mean_values[:, 1] ≈ [1.0, 1.0]
        @test std_values[:, 1] == [0.0, 0.0]
    end

    @testset "Input validation" begin
        H = spzeros(Float64, 2, 2)
        L = spzeros(Float64, 2, 2)

        ψ0 = sparsevec([1], [1.0], 2)
        Z = sparse([1.0 0.0; 0.0 -1.0])

        system = TrajectorySystem(H, [L])

        @test_throws ArgumentError evolve(
            system,
            ψ0,
            Float64[],
            [Z],
        )

        @test_throws ArgumentError evolve(
            system,
            ψ0,
            [0.0, 1.0, 0.5],
            [Z],
        )

        @test_throws ArgumentError evolve(
            system,
            ψ0,
            [0.0, 0.1],
            [Z];
            num_samples=0,
        )
    end

    @testset "Constructor validation" begin
        H = spzeros(Float64, 2, 2)
        bad_L = spzeros(Float64, 3, 3)

        @test_throws ArgumentError TrajectorySystem(
            H,
            [bad_L],
        )
    end



    @testset "Construct from SpinModel" begin
        model = SpinModel(
            1;
            Jxy=0.0,
            Jz=0.0,
        )

        L = spzeros(Float64, 2, 2)

        system = TrajectorySystem(
            model,
            [L],
        )

        @test system.hamiltonian == complex.(model.hamiltonian)
        @test length(system.lindblad_ops) == 1
    end
end