export assign_secondary_structure!

import AssigningSecondaryStructure as ASS

"""
    assign_secondary_structure!(protein)
"""
function assign_secondary_structure!(protein::Protein)
    ss_num_vectors = ASS.assign_secondary_structure([chain.backbone.coords for chain in protein])
    for (chain, ss_num_vector) in zip(protein, ss_num_vectors)
        ssvector = SecondaryStructure.(ss_num_vector)
        @assert length(chain.ssvector) == length(ssvector)
        chain.ssvector .= ssvector
    end
end