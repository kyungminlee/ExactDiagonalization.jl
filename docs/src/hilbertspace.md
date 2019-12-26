# Hilbert space


## Site

The definition of a quantum many-body problem starts by defining the Hilbert space.
The [`Site`](@ref) serves as a unit Hilbert space, and the Hilbert space for whole
system can be constructed by taking the tensor product of them.

A [`Site`](@ref) can be constructed out of a set of [`State`](@ref). For example,
```julia
spinsite = Site{Int}([State{Int}("Up", 1), State{Int}("Dn", -1)])
```
constructs a two-state site with spin-half degrees of freedom.
The type parameter `Int` is the type of the Abelian quantum number, which, in this
case, is $2S_z$. Each basis vector is represented as a (0-based) binary number,
corresponding to their order in the constructor. For the example above, the up-state
is represented by a `0` and the down-state is represented by a `1`.

## HilbertSpace

We can combine multiple sites to form a [`HilbertSpace`](@ref). To construct a Hilbert
space from the spin-half sites as defined above,
```julia
hilbert_space = HilbertSpace([spinsite, spinsite, spinsite, spinsite])
```

Note that all the basis vectors of the Hilbert space will be represented as a binary
number, where each `Site` occupies a fixed location and width. e.g.
```
|↑↑↑↑⟩ = |0000⟩
|↑↑↑↓⟩ = |0001⟩
|↑↑↓↑⟩ = |0010⟩
       ⋮
|↓↓↓↓⟩ = |1111⟩
```

## HilbertSpaceSector

A sub Hilbert space, in terms of the Abelian quantum number, can be constructed
using [`HilbertSpaceSector`](@ref), by specifying the value of the quantum number
```julia
hilbert_space_sector = HilbertSpaceSector(hilbert_space, 0)
```
or a set of quantum number values, if you need to for whatever reason
```julia
hilbert_space_sector = HilbertSpaceSector(hilbert_space, [0,2,4])
```
