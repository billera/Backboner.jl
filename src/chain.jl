export Chain

# TODO: use symbols for IDs?

"""
    Chain{T} <: AbstractVector{Residue}

A chain has an identifier (usually a single letter) and holds the backbone atom coordinates, amino acid sequence, and secondary structures of a protein chain. 
"""
struct Chain{T} <: AbstractVector{Residue}
    len::Int
    id::AbstractString
    backbone::Backbone{3, T}
    aavector::Vector{Char}
    ssvector::Vector{Char}
    oxygens::AbstractMatrix{T}

    function Chain(
        id::AbstractString,
        backbone::Backbone{3, T};
        aavector::Vector{Char} = fill('G', length(backbone)),
        ssvector::Union{Vector{Char}, Vector{<:Integer}} = fill(' ', length(backbone)),
        oxygens::AbstractMatrix{T} = get_oxygens(backbone),
    ) where T
        len = length(backbone)
        @assert len == length(aavector) == length(ssvector) "backbone, aavector, and ssvector must have the same length"
        ssvector isa Vector{<:Integer} && (ssvector = get.("-HE", ssvector, ' '))
        return new{T}(len, id, backbone, aavector, ssvector, oxygens)
    end

    Chain(backbone::Backbone; kwargs...) = Chain("_", backbone; kwargs...) 
end

@inline Base.:(==)(chain1::Chain, chain2::Chain) = chain1.id == chain2.id && chain1.backbone == chain2.backbone && chain1.ssvector == chain2.ssvector
@inline Base.length(chain::Chain) = length(chain.backbone)
@inline Base.size(chain::Chain) = (length(chain),)
@inline Base.getindex(chain::Chain, i::Integer) = Residue(i, chain.backbone, chain.aavector[i], chain.ssvector[i])

function Base.getproperty(chain::Chain, property::Symbol)
    if property in fieldnames
        return getfield(chain, property)
    elseif property in chain.backbone.atomnames
        return getindex(chain.backbone, property)
    else
        error("Chain does not have a field or atom named $property")
    end
end

Base.summary(chain::Chain) = "Chain $(chain.id) with $(length(chain)) residue$(length(chain) == 1 ? "" : "s")"
Base.show(io::IO, chain::Chain) = print(io, summary(chain))

has_complete_ss(chain::Chain) = has_complete_ss(chain.ssvector)
