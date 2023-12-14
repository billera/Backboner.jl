export Backbone, atom_coords

"""
    Backbone{A, T <: Real} <: AbstractBackbone{3, A, T}

A wrapper for a 3xNxL array of coordinates of atoms.
Backbone{3} is used to store 3-dimensional coordinates of the contiguous backbone atoms (N, Cα, C) of a protein chain.
"""
struct Backbone{A, T <: Real} <: AbstractBackbone{3, A, T}
    coords::AbstractArray{T, 3}
    atomnames::Vector{Symbol} # use StaticVector instead?

    function Backbone{A}(
        coords::AbstractArray{T, 3};
        atomnames::Vector{Symbol} = [Symbol("ATOM$i") for i in 1:A],
    ) where {A, T}
        @assert size(coords, 1) == 3 "coords must have 3 coordinates per atom"
        @assert size(coords, 2) == A "The second dimension of coords must be $A"
        @assert size(coords, 3) == length(atomnames) "atomnames must have the same length as the third dimension of coords (length)"
        return new{A, T}(coords, atomnames)
    end

    function Backbone{new_A}(backbone::Backbone{A}) where {A, new_A}
        @assert A*length(backbone) % new_A == 0 "Cannot convert Backbone{$A} to Backbone{$new_A}. The number of atoms ($(A*length(backbone))) must be a multiple of $new_A"
        return reshape(backbone.coords, 3, new_A, :)
    end

    Backbone(coords::AbstractArray{T, 3}) where T = Backbone{size(coords, 2)}(coords)
end

@inline Base.:(==)(bb1::Backbone, bb2::Backbone) = bb1.coords == bb2.coords
@inline Base.size(bb::Backbone) = size(bb.coords)
@inline Base.length(backbone::Backbone) = size(backbone, 3)
@inline Base.getindex(backbone::Backbone, i, j, k) = backbone.coords[i, j, k]
@inline Base.getindex(backbone::Backbone, i::Integer) = view(backbone.coords, :, :, i)
@inline Base.getindex(backbone::Backbone, r::UnitRange{Int}) = Backbone(view(backbone.coords, :, :, r))
@inline Base.reshape(backbone::Backbone, dims...) = Backbone(reshape(backbone.coords, dims...))

"""
    atom_coords(backbone, i)

Returns the coordinates of specific columns of atoms in a backbone.
"""
@inline atom_coords(backbone::Backbone, i) = view(backbone.coords, :, i, :)

@inline function Base.getproperty(backbone::Backbone, property::Symbol)
    if property in fieldnames
        return @inline getfield(backbone, property)
    elseif property in backbone.atomnames
        return atom_coords(backbone, findfirst(==(property), backbone.atomnames))
    else
        error("Backbone does not have a field or atom named $property")
    end
end

include("rotations.jl")
include("oxygen.jl")
include("bonds.jl")
include("dihedrals.jl")