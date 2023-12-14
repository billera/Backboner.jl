@testset "distance.jl" begin

    protein = read_pdb("data/1ASS.pdb")
    backbone = protein["A"].backbone
    cn_distances = carbonyl_nitrogen_distances(backbone)
    @test length(cn_distances) == length(backbone) - 1
    @test 1.32 <= sum(cn_distances) / length(cn_distances) <= 1.34

end