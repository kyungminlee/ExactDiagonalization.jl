mutable struct AbstractHilbertSpace{QN}
  sites ::Vector{Site{QN}}
  bitwidths ::Vector{UInt}
  bitoffsets ::Vector{UInt}

  AbstractHilbertSpace{QN}() where {QN} = new{QN}([], [], [0])
  function AbstractHilbertSpace(sites ::AbstractArray{Site{QN}, 1}) where QN
      bitwidths = UInt[bitwidth(site) for site in sites]
      bitoffsets = UInt[0, cumsum(bitwidths)...]
      new{QN}(sites, bitwidths, bitoffsets)
  end
end

function add_site!(hs ::AbstractHilbertSpace{QN}, site ::Site{QN}) where QN
  push!(hs.sites, site)
  bw = bitwidth(site)
  push!(hs.bitwidths, bw)
  push!(hs.bitoffsets, hs.bitoffsets[end] + bw)
end

function quantum_number_sectors(hs ::AbstractHilbertSpace{QN}) where QN
  qns = Set{QN}([zero(QN)])
  for site in hs.sites
    qns_next = Set{QN}()
    for state in site.states
      for q in qns
        push!(qns_next, q + state.quantum_number)
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
  mask = make_bitmask(hs.bitoffsets[isite+1]; lsb=hs.bitoffsets[isite], dtype=U)
  return (binrep & (~mask)) | (U(new_state_index-1) << hs.bitoffsets[isite] )
end

function get_state_index(hs ::AbstractHilbertSpace, binrep ::U, isite ::Integer) where {U<:Unsigned}
  return Int( ( binrep >> hs.bitoffsets[isite] ) & make_bitmask(hs.bitwidths[isite]; dtype=U) ) + 1
end

function get_state(hs ::AbstractHilbertSpace, binrep ::U, isite ::Integer) where {U<:Unsigned}
  return hs.sites[get_state_index(hs, binrep, isite)]
end
