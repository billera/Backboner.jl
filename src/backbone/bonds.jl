export atom_distances, carbonyl_nitrogen_distances

"""
    atom_displacements(backbone::Backbone, atom1::Integer, atom2::Integer, residue_offset::Integer)
"""
function atom_displacements(backbone::Backbone{A}, atom1::Integer, atom2::Integer, residue_offset::Integer=0) where A
    @assert 1 <= atom1 <= A && 1 <= atom2 <= A "Backbone{$A} does not have atoms $atom1 and $atom2"
    coords = backbone.coords
    displacements = @view(coords[:, atom2, 1+residue_offset:end]) .- @view(coords[:, atom1, 1:end-residue_offset])
    return displacements
end

"""
    atom_distances(backbone::Backbone, atom1::Integer, atom2::Integer, residue_offset::Integer)

Calculate the distances between all pairs of two types of atoms in a backbone, e.g.
the distances between all pairs of contiguous carbonyl and nitrogen atoms.
atom1 and atom2 are the indices of the atoms in the backbone, and residue_offset
is the number of residues between the atoms (0 by default).

To calculate the distances between all pairs of contiguous carbonyl and nitrogen atoms,
use atom_distances(backbone, 3, 1, 1), since the carbonyl oxygen is the third atom in the backbone,
the nitrogen is the first atom, and the nitrogen is on the next residue, so the offset is one.

Returns a vector of distances of length (length(backbone) - residue_offset).
"""
function atom_distances(backbone::Backbone{A}, atom1::Integer, atom2::Integer, residue_offset::Integer=0) where A
    @assert 1 <= atom1 <= A && 1 <= atom2 <= A "Backbone{$A} does not have atoms $atom1 and/or $atom2"
    displacements = atom_displacements(backbone, atom1, atom2, residue_offset)
    distances = reshape(mapslices(norm, displacements, dims=1), :)
    return distances
end

bond_vectors(backbone::Backbone) = atom_displacements(Backbone{1}(backbone), 1, 1, 1)
bond_lengths(backbone::Backbone) = atom_distances(Backbone{1}(backbone), 1, 1, 1)

"""
    carbonyl_nitrogen_distances(backbone::Backbone)

Calculate the distances between all pairs of contiguous carbonyl and nitrogen atoms in a backbone.
Returns a vector of distances of length (length(backbone) - 1).
"""
carbonyl_nitrogen_distances(backbone::Backbone) = atom_distances(backbone, 3, 1, 1)