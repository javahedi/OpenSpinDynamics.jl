using SparseArrays

@testset "Krylov evolution" begin
    @testset "Zero Hamiltonian" begin
        H = spzeros(Float64, 2, 2)
        ψ0 = sparsevec([1], [1.0], 2)

        Z = sparse([
             1.0  0.0
             0.0 -1.0
        ])

        times = [0.0, 0.1, 0.5, 1.0]
        system = KrylovArnoldiSystem(H)

        result = evolve(
            system,
            ψ0,
            times,
            [Z];
            method=:krylov,
        )

        @test size(result) == (length(times), 1)
        @test result[:, 1] ≈ ones(length(times))


       
    end

    @testset "Invalid method" begin
        H = spzeros(Float64, 2, 2)
        ψ0 = sparsevec([1], [1.0], 2)
        Z = sparse([1.0 0.0; 0.0 -1.0])

        system = KrylovArnoldiSystem(H)

        @test_throws ArgumentError evolve(
            system,
            ψ0,
            [0.0, 0.1],
            [Z];
            method=:invalid,
        )




    end
end


@testset "Arnoldi evolution" begin
    H = spzeros(Float64, 2, 2)
    ψ0 = sparsevec([1], [1.0], 2)

    Z = sparse([
         1.0  0.0
         0.0 -1.0
    ])

    times = [0.0, 0.1, 0.5, 1.0]
    system = KrylovArnoldiSystem(H)

    result = evolve(
        system,
        ψ0,
        times,
        [Z];
        method=:arnoldi,
    )

    @test size(result) == (length(times), 1)
    @test result[:, 1] ≈ ones(length(times))
end


@testset "Arnoldi agrees with Krylov" begin
    H = sparse([
        0.0 1.0
        1.0 0.0
    ])

    ψ0 = sparsevec([1], [1.0], 2)

    Z = sparse([
         1.0  0.0
         0.0 -1.0
    ])

    times = [0.0, 0.1, 0.25, 0.5, 1.0]
    system = KrylovArnoldiSystem(H)

    result_krylov = evolve(
        system,
        ψ0,
        times,
        [Z];
        method=:krylov,
    )

    result_arnoldi = evolve(
        system,
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