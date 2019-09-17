export AbstractHilbertSpace
export add_site!, quantum_number_sectors, get_quantum_number, extract, compress, update, get_state, get_state_index
export get_bitmask


"""
    AbstractHilbertSpace{QN}

Abstract Hilbert space with quantum number type `QN`.

# Examples
```jldoctest
julia> spin_site = Site{Int64}([State{Int64}("Up", +1), State{Int64}("Dn", -1)])
Site{Int64}(State{Int64}[State{Int64}("Up", 1), State{Int64}("Dn", -1)])

julia> hs = AbstractHilbertSpace{Int64}([spin_site, spin_site])
AbstractHilbertSpace{Int64}(Site{Int64}[Site{Int64}(State{Int64}[State{Int64}("Up", 1), State{Int64}("Dn", -1)]), Site{Int64}(State{Int64}[State{Int64}("Up", 1), State{Int64}("Dn", -1)])], [1, 1], [0, 1, 2])
```
"""
struct AbstractHilbertSpace{QN}
  sites ::Vector{Site{QN}}
  bitwidths ::Vector{Int}
  bitoffsets ::Vector{Int}

  AbstractHilbertSpace{QN}() where {QN} = new{QN}([], [], [0])
  function AbstractHilbertSpace(sites ::AbstractArray{Site{QN}, 1}) where QN
    bitwidths = Int[bitwidth(site) for site in sites]
    bitoffsets = Int[0, cumsum(bitwidths)...]
    new{QN}(sites, bitwidths, bitoffsets)
  end

  function AbstractHilbertSpace{QN}(sites ::AbstractArray{Site{QN}, 1}) where QN
    bitwidths = Int[bitwidth(site) for site in sites]
    bitoffsets = Int[0, cumsum(bitwidths)...]
    new{QN}(sites, bitwidths, bitoffsets)
  end  
end

import Base.==

function ==(lhs ::AbstractHilbertSpace{Q1}, rhs ::AbstractHilbertSpace{Q2}) where {Q1, Q2}
  return (Q1 == Q2) && (lhs.sites == rhs.sites) #&& (lhs.bitwidths == rhs.bitwidths) && (lhs.bitoffsets == rhs.bitoffsets)
end

function get_bitmask(hs ::AbstractHilbertSpace, isite ::Integer; dtype ::DataType=UInt)
  return make_bitmask(hs.bitoffsets[isite+1], hs.bitoffsets[isite]; dtype=dtype)
end

# function add_site!(hs ::AbstractHilbertSpace{QN}, site ::Site{QN}) where QN
#   push!(hs.sites, site)
#   bw = bitwidth(site)
#   push!(hs.bitwidths, bw)
#   push!(hs.bitoffsets, hs.bitoffsets[end] + bw)
# end

function quantum_number_sectors(hs ::AbstractHilbertSpace{QN}) where QN
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

function get_quantum_number(hs ::AbstractHilbertSpace{QN}, binrep ::BR) where {QN, BR}
  sum(
    let
      i = get_state_index(hs, binrep, isite)
      site.states[i].quantum_number
    end
    for (isite, site) in enumerate(hs.sites)
  )
end

function get_quantum_number(hs ::AbstractHilbertSpace{QN}, indexarray ::AbstractArray{I, 1}) where {BinRep, QN, I <:Integer}
    sum(
        site.states[indexarray[isite]].quantum_number
        for (isite, site) in enumerate(hs.sites)
    )
end

# function extract_binrep(hilbert_space ::AbstractHilbertSpace{QN}, isite ::Integer, binrep ::U) where {QN, U <:Unsigned}
#   mask = make_bitmask(hilbert_space.bidwidths[isite]; dtype=U)
#   return (binrep >> hilbert_space.bitoffsets[isite]) & mask
# end

"""
Convert binary representation to an array of indices (of states)

Examples
≡≡≡≡≡≡≡≡≡≡

```
```
"""
function extract(hs ::AbstractHilbertSpace{QN}, binrep ::U) where {QN, U <:Unsigned}
  out = Int[]
  for (isite, site) in enumerate(hs.sites)
    mask = make_bitmask(hs.bitwidths[isite]; dtype=U)
    index = Int(binrep & mask) + 1
    @assert (1 <= index <= length(site.states))
    push!(out, index)
    binrep = binrep >> hs.bitwidths[isite]
  end
  return out
end


"""
Convert an array of indices (of states) to binary representation 
"""
function compress(hs ::AbstractHilbertSpace{QN}, indexarray ::AbstractVector{I}; BR::DataType=UInt) where {QN, I<:Integer}
  @assert length(indexarray) == length(hs.sites)

  binrep = zero(BR)
  for (isite, site) in enumerate(hs.sites)
    binrep |= (BR(indexarray[isite] - 1) << hs.bitoffsets[isite] )
  end
  return binrep
end


function update(hs ::AbstractHilbertSpace, binrep ::U, isite ::Integer, new_state_index ::Integer) where {U <:Unsigned}
  mask = make_bitmask(hs.bitoffsets[isite+1], hs.bitoffsets[isite]; dtype=U)
  return (binrep & (~mask)) | (U(new_state_index-1) << hs.bitoffsets[isite] )
end

function get_state_index(hs ::AbstractHilbertSpace, binrep ::U, isite ::Integer) where {U<:Unsigned}
  return Int( ( binrep >> hs.bitoffsets[isite] ) & make_bitmask(hs.bitwidths[isite]; dtype=U) ) + 1
end

function get_state(hs ::AbstractHilbertSpace, binrep ::U, isite ::Integer) where {U<:Unsigned}
  return hs.sites[get_state_index(hs, binrep, isite)]
end
