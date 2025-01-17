using Backboner
using Test

using LinearAlgebra

@testset "Backboner.jl" begin

    include("backbone/backbone.jl")
    include("residue.jl")
    include("chain.jl")
    include("protein.jl")
    include("assign.jl")
    include("io.jl")

end
