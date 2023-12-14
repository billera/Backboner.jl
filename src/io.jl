export read_pdb, write_pdb

const THREE_LETTER_AA_CODES = Dict(
    'A' => "ALA", 'R' => "ARG", 'N' => "ASN", 'D' => "ASP",
    'C' => "CYS", 'Q' => "GLN", 'E' => "GLU", 'G' => "GLY",
    'H' => "HIS", 'I' => "ILE", 'L' => "LEU", 'K' => "LYS",
    'M' => "MET", 'F' => "PHE", 'P' => "PRO", 'S' => "SER",
    'T' => "THR", 'W' => "TRP", 'Y' => "TYR", 'V' => "VAL",
)

const ONE_LETTER_AA_CODES = Dict(v => k for (k, v) in THREE_LETTER_AA_CODES)

function collect_residues(atoms::Vector{PDBTools.Atom})
    residues = Vector{PDBTools.Atom}[]
    i = 1
    while i <= length(atoms) - 3 # Ensure there are at least four atoms left to process
        # Check if the next four atoms are N, CA, C, O in order
        if atoms[i].name == "N" && atoms[i+1].name == "CA" && atoms[i+2].name == "C" && atoms[i+3].name == "O" &&
                all(==(PDBTools.resnum(atoms[i])), PDBTools.resnum.(atoms[i+1:i+3]))
            push!(residues, atoms[i:i+3])
            i += 4
        else
            i += 1
        end
    end
    return residues
end

function Chain(atoms::Vector{PDBTools.Atom})
    id = PDBTools.chain(atoms[1])
    @assert allequal(PDBTools.chain.(atoms)) "atoms must be from the same chain"

    atoms_by_residue = collect_residues(atoms)
    backbone_coords = zeros(Float32, (3, 3, length(atoms_by_residue)))
    oxygens = zeros(Float32, (3, length(atoms_by_residue)))
    for (i, residue_atoms) in enumerate(atoms_by_residue)
        for (j, atom) in enumerate(residue_atoms)
            if j == 4
                oxygens[:, i] = [atom.x, atom.y, atom.z]
            else
                backbone_coords[:, j, i] = [atom.x, atom.y, atom.z]
            end
        end
    end
    backbone = Backbone(backbone_coords; atomnames=[:nitrogen, :alphacarbon, :carbon])

    aavector = [get(ONE_LETTER_AA_CODES, atom.resname, 'X') for atom in atoms if atom.name == "CA"]
    return Chain(id, backbone, aavector=aavector, oxygens=oxygens)
end

function Protein(atoms::Vector{PDBTools.Atom})
    filter!(a -> a.name in ["N", "CA", "C", "O"], atoms)
    ids = PDBTools.chain.(atoms)
    chains = [Chain(atoms[ids .== id]) for id in unique(ids)]
    return Protein(chains)
end

"""
    read_pdb(filename::String)

Assumes that each residue starts with four atoms: N, CA, C, O.
"""
read_pdb(filename::String) = Protein(PDBTools.readPDB(filename))

function write_pdb(protein::Protein, filename, header=:auto, footer=:auto)
    atoms = PDBTools.Atom[]
    index = 0
    residue_index = 0
    for chain in protein
        L = length(chain)
        for (resnum, (atom_coord_matrix, aa)) in enumerate(zip(eachslice(chain.backbone.coords, dims=3), chain.aavector))
            resname = get(THREE_LETTER_AA_CODES, aa, "XXX")
            residue_index += 1
            for (name, pos) in zip(["N", "CA", "C", "O"], [eachcol(atom_coord_matrix); chain.repeated[:oxygen][:, resnum]])
                index += 1
                atom = PDBTools.Atom(
                    index = index,
                    name = name,
                    resname = resname,
                    chain = chain.id,
                    resnum = resnum,
                    residue = residue_index,
                    x = pos[1], y = pos[2], z = pos[3],
                )
                push!(atoms, atom)
            end
        end

        # reuses the magic vector from src/backbone/oxygen.jl
        # note: the magic vector is meant to be used to calculate the O atom position,
        # but this is basically using it to get the next N position, 
        # so it's a hacky way to get an OXT atom position
        index += 1
        last_N, last_CA, last_C = eachcol(chain.backbone.coords[:, :, L])
        last_O = chain.repeated[:oxygen][:, L]
        rot_matrix = get_rotation_matrix(last_CA, last_C, last_O)
        OXT_pos = rot_matrix' \ magic_vector + last_C
        OXT_atom = PDBTools.Atom(
            index = index,
            name = "OXT",
            resname = get(THREE_LETTER_AA_CODES, chain.aavector[end], "XXX"),
            chain = chain.id,
            resnum = L,
            residue = residue_index,
            x = OXT_pos[1], y = OXT_pos[2], z = OXT_pos[3],
        )
        push!(atoms, OXT_atom)
    end
    PDBTools.writePDB(atoms, filename, header=header, footer=footer)
end