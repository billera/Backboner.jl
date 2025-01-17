export atom_distances, carbonyl_nitrogen_distances

function backbone_bond_vectors(backbone::Backbone{N}) where N
    @assert N >= 3 "backbone needs at least the N, Cα and C atoms to calculate bond vectors"
    backbone3 = Backbone(atom_coord_matrix(backbone, 1:3))
    flattened_coords = reshape(backbone3.coords, 3, :)
    bond_vectors = @view(flattened_coords[:, 2:end]) .- @view(flattened_coords[:, 1:end-1])
    return bond_vectors
end

function atom_displacements(backbone::Backbone{N}, atom1::Integer, atom2::Integer, residue_offset::Integer=0) where N
    @assert 1 <= atom1 <= N && 1 <= atom2 <= N "Backbone{$N} does not have atoms $atom1 and $atom2"
    coords = backbone.coords
    displacements = @view(coords[:, atom1, 1:end-residue_offset]) .- @view(coords[:, atom2, 1+residue_offset:end])
    return displacements
end

"""
    atom_distances(backbone::Backbone, atom1::Integer, atom2::Integer, residue_offset::Integer)

Calculate the distances between all pairs of two types atoms in a backbone, e.g.
the distances between all pairs of contiguous carbonyl and nitrogen atoms.
atom1 and atom2 are the indices of the atoms in the backbone, and residue_offset
is the number of residues between the atoms (0 by default).

Returns a vector of distances of length (length(backbone) - residue_offset).
"""
function atom_distances(backbone::Backbone{N}, atom1::Integer, atom2::Integer, residue_offset::Integer=0) where N
    @assert 1 <= atom1 <= N && 1 <= atom2 <= N "Backbone{$N} does not have atoms $atom1 and $atom2"
    displacements = atom_displacements(backbone, atom1, atom2, residue_offset)
    distances = reshape(mapslices(norm, displacements, dims=1), :)
    return distances
end

"""
    carbonyl_nitrogen_distances(backbone::Backbone)

Calculate the distances between all pairs of contiguous carbonyl and nitrogen atoms in a backbone.
Returns a vector of distances of length (length(backbone) - 1).
"""
carbonyl_nitrogen_distances(backbone::Backbone) = atom_distances(backbone, 3, 1, 1)