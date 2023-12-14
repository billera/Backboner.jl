export Residue, backbone_atom_coords

struct Residue{T}
    index::Integer
    backbone::Backbone{3, T}
    oxygen::AbstractVector{T}
    aa::Char
    ss::Char
end

backbone_atom_coords(residue::Residue) = residue.backbone[residue.index]

function Base.getproperty(residue::Residue, property::Symbol)
    if property in fieldnames
        return getfield(residue, property)
    elseif property in residue.backbone.atomnames
        return getindex(residue.backbone, property)[residue.index]
    
    else
        error("Residue does not have a field or atom named $property")
    end
end

function Base.summary(residue::Residue)
    index = lpad(string(residue.index), length(string(length(residue.backbone))))
    aa3 = get(THREE_LETTER_AA_CODES, residue.aa, "XXX")
    ss = residue.ss
    return "Residue $index $aa3 $ss"
end

Base.show(io::IO, residue::Residue) = print(io, summary(residue))