@testset "io.jl" begin

    @testset "PDB" begin

        @testset "read" begin
            protein = read_pdb("data/1ASS.pdb")
            @test length.(protein) == [152]
        end

        @testset "write" begin
            out = "temp.pdb"
            try
                protein = read_pdb("data/1ASS.pdb")
                write_pdb(protein, out)
                protein2 = read_pdb(out)
                @test protein == protein2
            finally
                rm(out, force=true)
            end
        end

    end

end