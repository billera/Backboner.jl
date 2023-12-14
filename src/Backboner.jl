module Backboner

using LinearAlgebra

import Rotations
import PDBTools

export has_complete_ss

has_complete_ss(ssvector::Vector{Char}) = all(!=(' '), ssvector)

include("backbone/backbone.jl")
include("residue.jl")
include("chain.jl")
include("protein.jl")
include("assign.jl")
include("io.jl")

end