export ReducedHilbertSpaceRepresentation
export bintype
export symmetry_reduce, symmetry_unreduce

import TightBindingLattice.TranslationGroup

struct ReducedHilbertSpaceRepresentation{HSR <:HilbertSpaceRepresentation, BR, C<:Complex} <:AbstractHilbertSpaceRepresentation{C}
  parent ::HSR
  translation_group ::TranslationGroup
  basis_list ::Vector{BR}
  #basis_mapping ::Vector{NamedTuple{(:index, :amplitude), Tuple{Int, C}}} # has size of parent dimension. each index item contains index at reduced basis, or -1 if not included
  basis_mapping_index ::Vector{Int}
  basis_mapping_amplitude ::Vector{C}
end

scalartype(arg ::Type{ReducedHilbertSpaceRepresentation{HSR, BR, C}}) where {HSR, BR, C} = C ::DataType
bintype(arg ::Type{ReducedHilbertSpaceRepresentation{HSR, BR, C}}) where {HSR, BR, C} = BR ::DataType

dimension(arg ::ReducedHilbertSpaceRepresentation{HSR, BR, C}) where {HSR, BR, C} = length(arg.basis_list)
