export HilbertSpace
export quantum_number_sectors, get_quantum_number, extract, compress, update, get_state, get_state_index
export get_bitmask
export bitwidth
export scalartype
export qntype
export basespace

"""
    HilbertSpace{QN}

Abstract Hilbert space with quantum number type `QN`.

# Examples
```jldoctest
julia> using ExactDiagonalization

julia> spin_site = Site([State{Int64}("Up", +1), State{Int64}("Dn", -1)]);

julia> hs = HilbertSpace([spin_site, spin_site])
HilbertSpace{Int64}(Site{Int64}[Site{Int64}(State{Int64}[State{Int64}("Up", 1), State{Int64}("Dn", -1)], GenericSiteType), Site{Int64}(State{Int64}[State{Int64}("Up", 1), State{Int64}("Dn", -1)], GenericSiteType)], [1, 1], [0, 1, 2])
```
"""
struct HilbertSpace{QN} <: AbstractHilbertSpace
  sites ::Vector{Site{QN}}
  bitwidths ::Vector{Int}
  bitoffsets ::Vector{Int}

  HilbertSpace{QN}() where {QN} = new{QN}([], [], [0])

  function HilbertSpace(sites ::AbstractArray{Site{QN}, 1}) where QN
    bitwidths = map(bitwidth, sites)
    bitoffsets = Int[0, cumsum(bitwidths)...]
    new{QN}(sites, bitwidths, bitoffsets)
  end

  function HilbertSpace{QN}(sites ::AbstractArray{Site{QN}, 1}) where QN
    bitwidths = map(bitwidth, sites)
    bitoffsets = Int[0, cumsum(bitwidths)...]
    new{QN}(sites, bitwidths, bitoffsets)
  end
end


"""
    scalartype(arg ::Type{HilbertSpace{QN}})

Returns the scalar type of the given hilbert space type.
For HilbertSpace{QN}, it is always `Bool`.
"""
scalartype(arg ::Type{<:HilbertSpace}) = Bool
scalartype(arg ::HilbertSpace) = Bool


import Base.valtype
"""
    valtype(arg ::Type{HilbertSpace{QN}})

Returns the `valtype` (scalar type) of the given hilbert space type.
"""
valtype(arg ::Type{<:HilbertSpace}) = Bool
valtype(arg ::HilbertSpace) = Bool


"""
    qntype(arg ::Type{HilbertSpace{QN}})

Returns the quantum number type of the given hilbert space type.
"""
qntype(arg ::Type{HilbertSpace{QN}}) where QN = QN
qntype(arg ::HilbertSpace{QN}) where QN = QN


"""
    basespace(hs)

Get the base space of the HilbertSpace `hs`, which is itself.
"""
basespace(hs::HilbertSpace) = hs


"""
Total number of bits

```jldoctest
julia> using ExactDiagonalization

julia> spin_site = Site([State{Int64}("Up", +1), State{Int64}("Dn", -1)])
Site{Int64}(State{Int64}[State{Int64}("Up", 1), State{Int64}("Dn", -1)], GenericSiteType)

julia> hs = HilbertSpace{Int64}([spin_site, spin_site, spin_site,])
HilbertSpace{Int64}(Site{Int64}[Site{Int64}(State{Int64}[State{Int64}("Up", 1), State{Int64}("Dn", -1)], GenericSiteType), Site{Int64}(State{Int64}[State{Int64}("Up", 1), State{Int64}("Dn", -1)], GenericSiteType), Site{Int64}(State{Int64}[State{Int64}("Up", 1), State{Int64}("Dn", -1)], GenericSiteType)], [1, 1, 1], [0, 1, 2, 3])

julia> bitwidth(hs)
3
```
"""
bitwidth(hs::HilbertSpace) = hs.bitoffsets[end]


import Base.==
function (==)(lhs ::HilbertSpace{Q1}, rhs ::HilbertSpace{Q2}) where {Q1, Q2}
  return lhs.sites == rhs.sites
end


function get_bitmask(hs ::HilbertSpace, isite ::Integer, binary_type ::Type{BR}=UInt) ::BR where {BR<:Unsigned}
  return make_bitmask(hs.bitoffsets[isite+1], hs.bitoffsets[isite], BR)
end


"""
    quantum_number_sectors
"""
function quantum_number_sectors(hs ::HilbertSpace{QN})::Vector{QN} where QN
  qns = Set{QN}([zero(QN)])
  for site in hs.sites
    qns_next = Set{QN}()
    for state in site.states, q in qns
      push!(qns_next, q .+ state.quantum_number)
    end
    qns = qns_next
  end
  return sort(collect(qns))
end


"""
    get_quantum_number
"""
function get_quantum_number(hs ::HilbertSpace{QN}, binrep ::BR) where {QN, BR}
  sum(
    let i = get_state_index(hs, binrep, isite)
      site.states[i].quantum_number
    end
    for (isite, site) in enumerate(hs.sites)
  )
end


function get_quantum_number(hs ::HilbertSpace{QN}, indexarray ::AbstractArray{I, 1}) where {QN, I <:Integer}
    sum(
      site.states[indexarray[isite]].quantum_number
      for (isite, site) in enumerate(hs.sites)
    )
end


"""
Convert binary representation to an array of indices (of states)

# Examples
```jldoctest
julia> using ExactDiagonalization

julia> spin_site = Site([State{Int64}("Up", +1), State{Int64}("Dn", -1)]);

julia> hs = HilbertSpace([spin_site, spin_site]);

julia> extract(hs, 0x03)
CartesianIndex(2, 2)
```
"""
function extract(hs ::HilbertSpace{QN}, binrep ::BR) ::CartesianIndex where {QN, BR<:Unsigned}
  out = Int[]
  for (isite, site) in enumerate(hs.sites)
    @inbounds mask = make_bitmask(hs.bitwidths[isite], BR)
    index = Int(binrep & mask) + 1
    @boundscheck if !(1 <= index <= length(site.states))
      throw(BoundsError(1:length(site.states), index))
    end
    push!(out, index)
    binrep = binrep >> hs.bitwidths[isite]
  end
  return CartesianIndex(out...)
end


"""
    compress(hs, indexarray ::CartesianIndex, binary_type=UInt)

Convert a cartesian index (a of state) to its binary representation

# Examples
```jldoctest
julia> using ExactDiagonalization

julia> spin_site = Site([State{Int64}("Up", +1), State{Int64}("Dn", -1)]);

julia> hs = HilbertSpace([spin_site, spin_site]);

julia> compress(hs, CartesianIndex(2,2))
0x0000000000000003
```
"""
function compress(hs ::HilbertSpace{QN}, indexarray ::CartesianIndex, binary_type ::Type{BR}=UInt) where {QN, BR<:Unsigned}
  if length(indexarray) != length(hs.sites)
    throw(ArgumentError("length of indexarray should be the number of sites"))
  end
  binrep = zero(BR)
  for (isite, site) in enumerate(hs.sites)
    @boundscheck if !(1 <= indexarray[isite] <= dimension(site))
      throw(BoundsError(1:dimension(site), indexarray[isite]))
    end
    binrep |= (BR(indexarray[isite] - 1) << hs.bitoffsets[isite] )
  end
  return binrep
end


"""
    update(hs, binrep, isite, new_state_index)

Update the binary representation of a basis state
by changing the state at site `isite` to a new local state specified by
`new_state_index`.
"""
@inline function update(hs ::HilbertSpace, binrep ::BR, isite ::Integer, new_state_index ::Integer) where {BR<:Unsigned}
  @boundscheck if !(1 <= new_state_index <= dimension(hs.sites[isite]))
    throw(BoundsError(1:dimension(hs.sites[isite]), new_state_index))
  end
  mask = make_bitmask(hs.bitoffsets[isite+1], hs.bitoffsets[isite], BR)
  return (binrep & (~mask)) | (BR(new_state_index-1) << hs.bitoffsets[isite])
end


"""
    get_state_index(hs, binrep, isite)

Get the *index of* the local state at site `isite` for the basis state
represented by `binrep`.
"""
function get_state_index(hs ::HilbertSpace, binrep ::BR, isite ::Integer) where {BR<:Unsigned}
  return Int( ( binrep >> hs.bitoffsets[isite] ) & make_bitmask(hs.bitwidths[isite], BR) ) + 1
end


"""
    get_state(hs, binrep, isite)

Get the local state at site `isite` for the basis state
represented by `binrep`. Returns an object of type `State`
"""
function get_state(hs ::HilbertSpace, binrep ::BR, isite ::Integer) where {BR<:Unsigned}
  return hs.sites[isite].states[get_state_index(hs, binrep, isite)]
end


# import Base.iterate
# @inline function iterate(hs ::HilbertSpace{QN}) where {QN}
#   subiterator = Iterators.product((1:length(site.states) for site in hs.sites)...)
#   next = Base.iterate(subiterator)
#   next === nothing && return nothing
#   value, next_substate = next
#   return (Int[value...], (subiterator, next_substate))
# end
#
# import Base.iterate
# @inline function iterate(hs ::HilbertSpace{QN}, state) where {QN}
#   (subiterator, substate) = state
#   next = Base.iterate(subiterator, substate)
#   next === nothing && return nothing
#   value, next_substate = next
#   return (Int[value...], (subiterator, next_substate))
# end


import Base.keys
function keys(hs ::HilbertSpace{QN}) where QN
  return CartesianIndices( ((1:length(site.states) for site in hs.sites)..., ) )
end
