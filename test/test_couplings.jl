@testset "Couplings" begin
    α = 2.0
    L = 100
    N = 10

    coupling_disorder = LongRangeCouplingDisorder(α, L, N)
    @test size(get_matrix(coupling_disorder)) == (N, N)

    coupling_clean = LongRangeCouplingClean(α, N)
    J = get_matrix(coupling_clean)

    @test size(J) == (N, N)
    @test isapprox(J, transpose(J); atol=1e-8)
end