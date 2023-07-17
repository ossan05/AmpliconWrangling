using Test
using .GeneralAlignment

@testset "Basic tests" begin

    @test general_pairwise_aligner((A, B, 0, 0.5, ((1, 1), 0, (0, 1), 2, (0, 3), 2, (1, 0), 2, (3, 0), 2))) 
end