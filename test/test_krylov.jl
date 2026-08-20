using SparseArrays

@testset "Krylov evolution" begin
    @testset "Zero Hamiltonian" begin
        ψ0 = sparsevec([1], [1.0], 2)

        Z = sparse([
             1.0  0.0
             0.0 -1.0
        ])

        times = [0.0, 0.1, 0.5, 1.0]

        model = SpinModel(
            1;
            Jxy=0.0,
            Jz=0.0,
            hx=[0.0],
            hz=[0.0],
        )

        result = evolve(
            model,
            ψ0,
            times,
            [Z];
            method=:krylov,
        )

        @test size(result) == (length(times), 1)
        @test result[:, 1] ≈ ones(length(times))
    end

    @testset "Invalid method" begin
        ψ0 = sparsevec([1], [1.0], 2)
        Z = sparse([1.0 0.0; 0.0 -1.0])

        model = SpinModel(
            1;
            Jxy=0.0,
            Jz=0.0,
        )

        @test_throws ArgumentError evolve(
            model,
            ψ0,
            [0.0, 0.1],
            [Z];
            method=:invalid,
        )
    end
end


@testset "Arnoldi evolution" begin
    ψ0 = sparsevec([1], [1.0], 2)

    Z = sparse([
         1.0  0.0
         0.0 -1.0
    ])

    times = [0.0, 0.1, 0.5, 1.0]

    model = SpinModel(
        1;
        Jxy=0.0,
        Jz=0.0,
        hx=[0.0],
        hz=[0.0],
    )

    result = evolve(
        model,
        ψ0,
        times,
        [Z];
        method=:arnoldi,
    )

    @test size(result) == (length(times), 1)
    @test result[:, 1] ≈ ones(length(times))
end


@testset "Arnoldi agrees with Krylov" begin
    ψ0 = sparsevec([1], [1.0], 2)

    Z = sparse([
         1.0  0.0
         0.0 -1.0
    ])

    times = [0.0, 0.1, 0.25, 0.5, 1.0]

    # For N = 1 and hx = 1:
    # H = σx
    model = SpinModel(
        1;
        Jxy=0.0,
        Jz=0.0,
        hx=[1.0],
        hz=[0.0],
    )

    result_krylov = evolve(
        model,
        ψ0,
        times,
        [Z];
        method=:krylov,
    )

    result_arnoldi = evolve(
        model,
        ψ0,
        times,
        [Z];
        method=:arnoldi,
    )

    @test size(result_krylov) == size(result_arnoldi)
    @test result_arnoldi ≈ result_krylov atol=1e-10

    expected = cos.(2 .* times)

    @test result_krylov[:, 1] ≈ expected atol=1e-10
    @test result_arnoldi[:, 1] ≈ expected atol=1e-10
end