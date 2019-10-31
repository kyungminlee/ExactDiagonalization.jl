export AbstractHilbertSpaceRepresentation
export bintype

abstract type AbstractHilbertSpaceRepresentation{S<:Number} end

# need scalartype for the hilbert space representation also, since layer we will have
# symmetry reduction

@inline scalartype(lhs ::AbstractHilbertSpaceRepresentation{S}) where S = S
@inline bintype(lhs ::AbstractHilbertSpaceRepresentation{S}) where S = bintype(typeof(lhs)) ::DataType
