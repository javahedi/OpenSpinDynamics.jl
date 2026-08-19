@testset "Pauli operators" begin
    N = 3
    paulis = generate_operators(N)

    @test size(paulis.Xop[1]) == (2^N, 2^N)
    @test size(paulis.Zop[2]) == (2^N, 2^N)
    @test size(paulis.Pop[3]) == (2^N, 2^N)
    @test size(paulis.Nop[1]) == (2^N, 2^N)
end
