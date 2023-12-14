@testset "assign.jl" begin
    
    @testset "assign_secondary_structure!" begin
        protein = read_pdb("data/1ASS.pdb")
        @test !has_complete_ss(protein)
        assign_secondary_structure!(protein)
        @test has_complete_ss(protein)
    end

    @testset "assign_secondary_structure" begin
        protein = read_pdb("data/1ASS.pdb")
        @test !has_complete_ss(protein)
        new_protein = assign_secondary_structure(protein)
        @test !has_complete_ss(protein)
        @test has_complete_ss(new_protein)
    end

end