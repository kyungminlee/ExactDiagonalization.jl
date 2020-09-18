export AbstractHilbertSpaceRepresentation
export scalartype, bintype


"""
    AbstractHilbertSpaceRepresentation{S}
"""
abstract type AbstractHilbertSpaceRepresentation{S<:Number} end


scalartype(lhs::AbstractHilbertSpaceRepresentation{S}) where S = S
Base.valtype(lhs::AbstractHilbertSpaceRepresentation{S}) where S = S
bintype(lhs::AbstractHilbertSpaceRepresentation{S}) where S = bintype(typeof(lhs))


bitwidth(lhs::AbstractHilbertSpaceRepresentation) = bitwidth(basespace(lhs))
