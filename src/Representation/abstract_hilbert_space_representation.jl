export AbstractHilbertSpaceRepresentation
export scalartype, bintype

"""
    AbstractHilbertSpaceRepresentation{S}
"""
abstract type AbstractHilbertSpaceRepresentation{S<:Number} end

# need scalartype for the hilbert space representation also, since layer we will have
# symmetry reduction

import Base.valtype

scalartype(lhs ::AbstractHilbertSpaceRepresentation{S}) where S = S
valtype(lhs ::AbstractHilbertSpaceRepresentation{S}) where S = S
bintype(lhs ::AbstractHilbertSpaceRepresentation{S}) where S = bintype(typeof(lhs)) ::DataType
