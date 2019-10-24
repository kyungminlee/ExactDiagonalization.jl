export HilbertSpace
export quantum_number_sectors, get_quantum_number, extract, compress, update, get_state, get_state_index
export get_bitmask
export bitwidth
export qntype
export basespace

"""
    HilbertSpace{QN}

Abstract Hilbert space with quantum number type `QN`.

# Examples
```jldoctest
using ExactDiagonalization
julia> spin_site = Site{Int64}([State{Int64}("Up", +1), State{Int64}("Dn", -1)])
Site{Int64}(State{Int64}[State{Int64}("Up", 1), State{Int64}("Dn", -1)])

julia> hs = HilbertSpace{Int64}([spin_site, spin_site])
HilbertSpace{Int64}(Site{Int64}[Site{Int64}(State{Int64}[State{Int64}("Up", 1), State{Int64}("Dn", -1)]), Site{Int64}(State{Int64}[State{Int64}("Up", 1), State{Int64}("Dn", -1)])], [1, 1], [0, 1, 2])
```
"""
struct HilbertSpace{QN} <: AbstractHilbertSpace
  sites ::Vector{Site{QN}}
  bitwidths ::Vector{Int}
  bitoffsets ::Vector{Int}

  HilbertSpace{QN}() where {QN} = new{QN}([], [], [0])
  function HilbertSpace(sites ::AbstractArray{Site{QN}, 1}) where QN
    bitwidths = Int[bitwidth(site) for site in sites]
    bitoffsets = Int[0, cumsum(bitwidths)...]
    new{QN}(sites, bitwidths, bitoffsets)
  end

  function HilbertSpace{QN}(sites ::AbstractArray{Site{QN}, 1}) where QN
    bitwidths = Int[bitwidth(site) for site in sites]
    bitoffsets = Int[0, cumsum(bitwidths)...]
    new{QN}(sites, bitwidths, bitoffsets)
  end
end

import Base.eltype
@inline eltype(arg ::HilbertSpace{QN}) where QN = Bool
@inline eltype(arg ::Type{HilbertSpace{QN}}) where QN = Bool

"""
    qntype


"""
@inline qntype(arg ::HilbertSpace{QN}) where QN = QN
@inline qntype(arg ::Type{HilbertSpace{QN}}) where QN = QN


"""
Total number of bits

```jldoctest
using ExactDiagonalization
julia> spin_site = Site{Int64}([State{Int64}("Up", +1), State{Int64}("Dn", -1)])
Site{Int64}(State{Int64}[State{Int64}("Up", 1), State{Int64}("Dn", -1)])

julia> hs = HilbertSpace{Int64}([spin_site, spin_site, spin_site,])
HilbertSpace{Int64}(Site{Int64}[Site{Int64}(State{Int64}[State{Int64}("Up", 1), State{Int64}("Dn", -1)]), Site{Int64}(State{Int64}[State{Int64}("Up", 1), State{Int64}("Dn", -1)]), Site{Int64}(State{Int64}[State{Int64}("Up", 1), State{Int64}("Dn", -1)])], [1, 1, 1], [0, 1, 2, 3])

julia> bitwidth(hs)
3
```
"""
@inline bitwidth(hs::HilbertSpace) = hs.bitoffsets[end]


@inline basespace(hs::HilbertSpace) = hs

import Base.==
function (==)(lhs ::HilbertSpace{Q1}, rhs ::HilbertSpace{Q2}) where {Q1, Q2}
  return (Q1 == Q2) && (lhs.sites == rhs.sites) #&& (lhs.bitwidths == rhs.bitwidths) && (lhs.bitoffsets == rhs.bitoffsets)
end

function get_bitmask(hs ::HilbertSpace, isite ::Integer; dtype ::DataType=UInt)
  return make_bitmask(hs.bitoffsets[isite+1], hs.bitoffsets[isite]; dtype=dtype)
end

"""
    quantum_number_sectors
"""
function quantum_number_sectors(hs ::HilbertSpace{QN})::Vector{QN} where QN
  qns = Set{QN}([zero(QN)])
  for site in hs.sites
    qns_next = Set{QN}()
    for state in site.states
      for q in qns
        push!(qns_next, q .+ state.quantum_number)
      end
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
    let
      i = get_state_index(hs, binrep, isite)
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

Examples
≡≡≡≡≡≡≡≡≡≡

```
```
"""
function extract(hs ::HilbertSpace{QN}, binrep ::U) where {QN, U <:Unsigned}
  out = Int[]
  for (isite, site) in enumerate(hs.sites)
    mask = make_bitmask(hs.bitwidths[isite]; dtype=U)
    index = Int(binrep & mask) + 1
    @boundscheck if !(1 <= index <= length(site.states))
      throw(BoundsError(1:length(site.states), index))
    end
    push!(out, index)
    binrep = binrep >> hs.bitwidths[isite]
  end
  return out
end


"""
Convert an array of indices (of states) to binary representation
"""
function compress(hs ::HilbertSpace{QN}, indexarray ::AbstractVector{I}; BR::DataType=UInt) where {QN, I<:Integer}
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


function update(hs ::HilbertSpace, binrep ::U, isite ::Integer, new_state_index ::Integer) where {U<:Unsigned}
  @boundscheck if !(1 <= new_state_index <= dimension(hs.sites[isite]))
    throw(BoundsError(1:dimension(hs.sites[isite]), new_state_index))
  end
  mask = make_bitmask(hs.bitoffsets[isite+1], hs.bitoffsets[isite]; dtype=U)
  return (binrep & (~mask)) | (U(new_state_index-1) << hs.bitoffsets[isite])
end

function get_state_index(hs ::HilbertSpace, binrep ::U, isite ::Integer) where {U<:Unsigned}
  return Int( ( binrep >> hs.bitoffsets[isite] ) & make_bitmask(hs.bitwidths[isite]; dtype=U) ) + 1
end

function get_state(hs ::HilbertSpace, binrep ::U, isite ::Integer) where {U<:Unsigned}
  return hs.sites[isite].states[get_state_index(hs, binrep, isite)]
end

import Base.iterate
@inline function iterate(hs ::HilbertSpace{QN}) where {QN}
  subiterator = Iterators.product((1:length(site.states) for site in hs.sites)...)
  next = Base.iterate(subiterator)
  next === nothing && return nothing
  value, next_substate = next
  return (Int[value...], (subiterator, next_substate))
end

import Base.iterate
@inline function iterate(hs ::HilbertSpace{QN}, state) where {QN}
  (subiterator, substate) = state
  next = Base.iterate(subiterator, substate)
  next === nothing && return nothing
  value, next_substate = next
  return (Int[value...], (subiterator, next_substate))
end
