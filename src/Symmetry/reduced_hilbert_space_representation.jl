export ReducedHilbertSpaceRepresentation
export bintype
export symmetry_reduce, symmetry_unreduce

import TightBindingLattice.TranslationGroup


"""
    ReducedHilbertSpaceRepresentation{HSR, BR, C}

Representation of the symmetry-reduced hilbert space.
Currently only supports Translation group (i.e. Abelian group).

# Members
- `parent ::HSR`
- `translation_group::TranslationGroup`
- `basis_list::Vector{BR}`
- `basis_mapping_index::Vector{Int}`
- `basis_mapping_amplitude::Vector{C}`
"""
struct ReducedHilbertSpaceRepresentation{HSR <:HilbertSpaceRepresentation, BR, C<:Complex} <:AbstractHilbertSpaceRepresentation{C}
  parent ::HSR
  translation_group ::TranslationGroup
  basis_list ::Vector{BR}
  basis_mapping_index ::Vector{Int} # has size of parent dimension. each index item contains index at reduced basis, or -1 if not included
  basis_mapping_amplitude ::Vector{C}
end

scalartype(arg ::Type{ReducedHilbertSpaceRepresentation{HSR, BR, C}}) where {HSR, BR, C} = C ::DataType
bintype(arg ::Type{ReducedHilbertSpaceRepresentation{HSR, BR, C}}) where {HSR, BR, C} = BR ::DataType

"""
    dimension(arg ::ReducedHilbertSpaceRepresentation{HSR, BR, C}) -> Int

Dimension of the given reduced hilbert space representation, i.e. number of basis elements.
"""
dimension(arg ::ReducedHilbertSpaceRepresentation{HSR, BR, C}) where {HSR, BR, C} = length(arg.basis_list)
