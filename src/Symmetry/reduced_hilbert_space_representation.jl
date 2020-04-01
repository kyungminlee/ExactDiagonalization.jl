export ReducedHilbertSpaceRepresentation
export bintype
export symmetry_reduce, symmetry_unreduce

import TightBindingLattice.TranslationSymmetry


"""
    ReducedHilbertSpaceRepresentation{HSR, BR, C}

Representation of the symmetry-reduced hilbert space.
Currently only supports Translation group (i.e. Abelian group).
```
"""
struct ReducedHilbertSpaceRepresentation{HSR<:HilbertSpaceRepresentation, BR, C<:Complex} <:AbstractHilbertSpaceRepresentation{C}
  parent::HSR
  #translation_group::TranslationGroup
  translation_symmetry::TranslationSymmetry
  basis_list::Vector{BR}
  basis_mapping_index::Vector{Int} # has size of parent dimension. each index item contains index at reduced basis, or -1 if not included
  basis_mapping_amplitude::Vector{C}
end


import Base.valtype

valtype(arg::Type{ReducedHilbertSpaceRepresentation{HSR, BR, C}}) where {HSR, BR, C} = C
scalartype(arg::Type{ReducedHilbertSpaceRepresentation{HSR, BR, C}}) where {HSR, BR, C} = C
bintype(arg::Type{ReducedHilbertSpaceRepresentation{HSR, BR, C}}) where {HSR, BR, C} = BR


basespace(lhs::ReducedHilbertSpaceRepresentation{HSR, BR, C}) where {HSR, BR, C} = basespace(lhs.parent)


"""
    dimension(arg::ReducedHilbertSpaceRepresentation{HSR, BR, C}) -> Int

Dimension of the given reduced hilbert space representation, i.e. number of basis elements.
"""
dimension(arg::ReducedHilbertSpaceRepresentation{HSR, BR, C}) where {HSR, BR, C} = length(arg.basis_list)
