export ReducedHilbertSpaceRealization
export bintype

struct ReducedHilbertSpaceRealization{QN, BR, C<:Complex}
  parent ::HilbertSpaceRealization{QN, BR}
  translation_group ::TranslationGroup
  basis_list ::Vector{BR}
  basis_mapping ::Vector{NamedTuple{(:index, :amplitude), Tuple{Int, C}}} # has size of parent dimension. each index item contains index at reduced basis, or -1 if not included
end

import Base.eltype
eltype(arg ::ReducedHilbertSpaceRealization{QN, BR, C}) where {QN, BR, C} = C
eltype(arg ::Type{ReducedHilbertSpaceRealization{QN, BR, C}}) where {QN, BR, C} = C

bintype(arg ::ReducedHilbertSpaceRealization{QN, BR, C}) where {QN, BR, C} = BR
bintype(arg ::Type{ReducedHilbertSpaceRealization{QN, BR, C}}) where {QN, BR, C} = BR

dimension(arg ::ReducedHilbertSpaceRealization{QN, BR, C}) where {QN, BR, C} = length(arg.basis_list)
