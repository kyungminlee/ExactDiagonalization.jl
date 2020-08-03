export ReducedHilbertSpaceRepresentation
export bintype
export symmetry_reduce, symmetry_unreduce


"""
    ReducedHilbertSpaceRepresentation{HSR, BR, C}

Representation of the symmetry-reduced hilbert space.
Currently only supports Translation group (i.e. Abelian group).
```
"""
struct ReducedHilbertSpaceRepresentation{
    HSR<:HilbertSpaceRepresentation,
    SIC<:AbstractSymmetryIrrepComponent,
    BR, C<:Complex
}<:AbstractHilbertSpaceRepresentation{C}
    parent::HSR
    symmetry_irrep_component::SIC
    basis_list::Vector{BR}
    basis_mapping_index::Vector{Int} # has size of parent dimension. each index item contains index at reduced basis, or -1 if not included
    basis_mapping_amplitude::Vector{C}
end


Base.valtype(arg::Type{ReducedHilbertSpaceRepresentation{HSR, SIC, BR, C}}) where {HSR, SIC, BR, C} = C
scalartype(arg::Type{ReducedHilbertSpaceRepresentation{HSR, SIC, BR, C}}) where {HSR, SIC, BR, C} = C
bintype(arg::Type{ReducedHilbertSpaceRepresentation{HSR, SIC, BR, C}}) where {HSR, SIC, BR, C} = BR


basespace(lhs::ReducedHilbertSpaceRepresentation) = basespace(lhs.parent)


"""
    dimension(arg::ReducedHilbertSpaceRepresentation{HSR, BR, C}) -> Int

Dimension of the given reduced hilbert space representation, i.e. number of basis elements.
"""
dimension(arg::ReducedHilbertSpaceRepresentation) = length(arg.basis_list)
