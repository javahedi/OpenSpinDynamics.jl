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