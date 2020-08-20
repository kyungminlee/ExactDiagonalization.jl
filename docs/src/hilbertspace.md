# Hilbert space


## Site

The definition of a quantum many-body problem starts by defining the Hilbert space.
The [`Site`](@ref) serves as a unit Hilbert space,
and the Hilbert space for whole system can be constructed by taking the tensor product of them.

A [`Site`](@ref) can be constructed out of a set of [`State`](@ref).
For example,
```julia-repl
julia> spinsite = Site([State("Up", 1), State("Dn", -1)])
Site{Tuple{Int64}}(State{Tuple{Int64}}[State{Tuple{Int64}}("Up", (1,)), State{Tuple{Int64}}("Dn", (-1,))])
```
constructs a two-state site with spin-half degrees of freedom.
The type parameter `Tuple{Int}` represents the type of Abelian quantum number.
which is is $2S_z$ in this case.
When there are more than one conserved quantum numbers, they can be combined:
e.g. `Tuple{Int, Int}`, to represent the charge and total $S_z$, for example.
Each basis vector is represented as a binary number,
corresponding to their order in the constructor (0-based).
For the example above, the up-state is represented as `0` and the down-state is represented as `1`.

## HilbertSpace

We can combine multiple sites to form a [`HilbertSpace`](@ref).
To construct a Hilbert space from the spin-half sites as defined above,
```julia
hilbert_space = HilbertSpace([spinsite, spinsite, spinsite, spinsite])
```

Note that all the basis vectors of the Hilbert space will be represented as a binary number,
where each `Site` occupies a fixed location and width. e.g.
```
|↑↑↑↑⟩ = |0000⟩
|↑↑↑↓⟩ = |0001⟩
|↑↑↓↑⟩ = |0010⟩
       ⋮
|↓↓↓↓⟩ = |1111⟩
```
The number of bits assigned for each site is determined by `Int(ceil(log2(length(site.states)))`, and can be accessed by [`bitwidth`](@ref).

## HilbertSpaceSector

A subspace of the whole Hilbert space, in terms of the Abelian quantum number, can be constructed using [`HilbertSpaceSector`](@ref), by specifying the value of the quantum number as an integer if the Hilbert space has a single integral quantum number,
```julia
hilbert_space_sector = HilbertSpaceSector(hilbert_space, 0)
```
or as a tuple
```julia
hilbert_space_sector = HilbertSpaceSector(hilbert_space, (0,))
```
You can also allow more than one quantum number values, if you need to for whatever reason
```julia
hilbert_space_sector = HilbertSpaceSector(hilbert_space, [(0,), (2,), (4,)])
```
or more shortly,
```julia
hilbert_space_sector = HilbertSpaceSector(hilbert_space, [0, 2, 4])
```
This example creates a subspace whose quantum numbers can be one of 0, 2, and 4.