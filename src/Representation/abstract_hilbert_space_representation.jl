export AbstractHilbertSpaceRepresentation
export bintype

abstract type AbstractHilbertSpaceRepresentation end

# need eltype for the hilbert space representation also, since layer we will have
# symmetry reduction
import Base.eltype
@inline eltype(lhs ::AbstractHilbertSpaceRepresentation) = eltype(typeof(lhs))
@inline eltype(lhs ::Type{AbstractHilbertSpaceRepresentation}) = error("eltype not implemented")

@inline bintype(lhs ::AbstractHilbertSpaceRepresentation) = bintype(typeof(lhs))
@inline bintype(lhs ::Type{AbstractHilbertSpaceRepresentation}) = error("bintype not implemented")
