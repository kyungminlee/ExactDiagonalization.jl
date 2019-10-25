export ReducedHilbertSpaceRepresentation
export bintype

struct ReducedHilbertSpaceRepresentation{HSR <:HilbertSpaceRepresentation, BR, C<:Complex}
  parent ::HSR
  translation_group ::TranslationGroup
  basis_list ::Vector{BR}
  basis_mapping ::Vector{NamedTuple{(:index, :amplitude), Tuple{Int, C}}} # has size of parent dimension. each index item contains index at reduced basis, or -1 if not included
end

@inline scalartype(arg ::ReducedHilbertSpaceRepresentation{HSR, BR, C}) where {HSR, BR, C} = C
@inline scalartype(arg ::Type{ReducedHilbertSpaceRepresentation{HSR, BR, C}}) where {HSR, BR, C} = C

@inline bintype(arg ::ReducedHilbertSpaceRepresentation{HSR, BR, C}) where {HSR, BR, C} = BR
@inline bintype(arg ::Type{ReducedHilbertSpaceRepresentation{HSR, BR, C}}) where {HSR, BR, C} = BR

dimension(arg ::ReducedHilbertSpaceRepresentation{HSR, BR, C}) where {HSR, BR, C} = length(arg.basis_list)
