export AbstractHilbertSpaceRepresentation
export bintype

abstract type AbstractHilbertSpaceRepresentation end

# need scalartype for the hilbert space representation also, since layer we will have
# symmetry reduction

@inline scalartype(lhs ::AbstractHilbertSpaceRepresentation) = scalartype(typeof(lhs))
@inline bintype(lhs ::AbstractHilbertSpaceRepresentation) = bintype(typeof(lhs))
