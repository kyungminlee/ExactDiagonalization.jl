export AbstractHilbertSpaceRepresentation
export scalartype, bintype


"""
    AbstractHilbertSpaceRepresentation{S}
"""
abstract type AbstractHilbertSpaceRepresentation{S<:Number} end


scalartype(::AbstractHilbertSpaceRepresentation{S}) where S = S
Base.valtype(::AbstractHilbertSpaceRepresentation{S}) where S = S
bintype(lhs::AbstractHilbertSpaceRepresentation{S}) where S = bintype(typeof(lhs))


bitwidth(lhs::AbstractHilbertSpaceRepresentation) = bitwidth(basespace(lhs))
