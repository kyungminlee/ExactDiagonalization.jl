export AbstractHilbertSpaceRepresentation
export scalartype, bintype


"""
    AbstractHilbertSpaceRepresentation{S}
"""
abstract type AbstractHilbertSpaceRepresentation{S<:Number} end


scalartype(lhs::AbstractHilbertSpaceRepresentation{S}) where S = S
Base.valtype(lhs::AbstractHilbertSpaceRepresentation{S}) where S = S
bintype(lhs::AbstractHilbertSpaceRepresentation{S}) where S = bintype(typeof(lhs))


function bitwidth(lhs::AbstractHilbertSpaceRepresentation{S}) where S
    return bitwidth(basespace(lhs))
end
