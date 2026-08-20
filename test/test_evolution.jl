@testset "Evolution validation" begin
    model = SpinModel(
        1;
        Jxy=0.0,
        Jz=0.0,
    )

    ψ0 = sparsevec([1], [1.0], 2)

    ρ0 = sparse([
        1.0 0.0
        0.0 0.0
    ])

    Z = sparse([
         1.0  0.0
         0.0 -1.0
    ])

    times = [0.0, 0.1, 0.2]

    @test_throws ArgumentError evolve(
        model,
        ψ0,
        times,
        [Z];
        method=:trajectories,
    )

    @test_throws ArgumentError evolve(
        model,
        ρ0,
        times,
        [Z];
        method=:expm,
    )
end