export assign_secondary_structure!, assign_secondary_structure

import AssigningSecondaryStructure: assign_secondary_structure!, assign_secondary_structure

"""
    assign_secondary_structure!(protein)

Uses a simplified version of DSSP to fill the secondary structure vector of each chain with '-' (coil/loop), 'H' (helix), and 'E' (strand).
"""
function assign_secondary_structure!(protein::Protein)
    coords_vector = [cat(chain.backbone.coords, chain.oxygen, dims=2) for chain in protein]
    ss_vectors = assign_secondary_structure(coords_vector)
    for (chain, ssvector) in zip(protein, ss_vectors)
        @assert length(chain.ssvector) == length(ssvector)
        chain.ssvector .= ssvector
    end
    return protein
end

"""
    assign_secondary_structure(protein)

Returns a new protein with secondary structure assigned.
"""
function assign_secondary_structure(protein::Protein)
    new_protein = deepcopy(protein)
    assign_secondary_structure!(new_protein)
    return new_protein
end