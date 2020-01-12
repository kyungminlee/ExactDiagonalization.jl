# Representation

## HilbertSpaceRepresentation

A [`HilbertSpaceRepresentation`](@ref) is a representation of the Hilbert space,
with the list of basis vectors in a ascending order of their binary representations, and a lookup table for them.
The [`HilbertSpaceRepresentation`](@ref) can be constructed using [`represent_array`](@ref)
which uses [`FrozenSortedArrayIndex`](@ref) for the lookup table,
or [`represent_dict`](@ref) which uses `Dict`.

```@docs
represent(::HilbertSpace)
represent_array(::HilbertSpace)
represent_dict(::HilbertSpace)
```
You can also explicitly supply a list of basis vectors
```@docs
represent(::HilbertSpace,::Vector{UInt})
represent_array(::HilbertSpace,::Vector{UInt})
represent_dict(::HilbertSpace,::Vector{UInt})
```

## OperatorRepresentation

An [`OperatorRepresentation`](@ref) is a representation of an operator in the given
Hilbert space representation.

```@docs
represent(::HilbertSpaceRepresentation,::AbstractOperator)
```
