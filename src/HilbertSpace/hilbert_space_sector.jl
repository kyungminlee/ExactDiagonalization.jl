export HilbertSpaceSector

struct HilbertSpaceSector{QN} <: AbstractHilbertSpace
  parent ::HilbertSpace{QN}
  allowed_quantum_numbers ::Set{QN}

  function HilbertSpaceSector(parent ::HilbertSpace{QN}) where QN
    sectors = quantum_number_sectors(parent)
    new{QN}(parent, Set(sectors))
  end

  function HilbertSpaceSector(parent ::HilbertSpace{QN}, allowed::QN) where QN
    sectors = Set{QN}(quantum_number_sectors(parent))
    new{QN}(parent, intersect(sectors, Set([allowed])))
  end

  function HilbertSpaceSector(parent ::HilbertSpace{QN},
                              allowed::Union{AbstractSet{QN}, AbstractVector{QN}}) where QN
    sectors = Set{QN}(quantum_number_sectors(parent))
    new{QN}(parent, intersect(sectors, Set(allowed)))
  end
end

import Base.eltype
@inline eltype(arg ::HilbertSpaceSector{QN}) where QN = Bool
@inline eltype(arg ::Type{HilbertSpaceSector{QN}}) where QN = Bool

export qntype
@inline qntype(arg ::HilbertSpaceSector{QN}) where QN = QN
@inline qntype(arg ::Type{HilbertSpaceSector{QN}}) where QN = QN

export basespace
@inline basespace(hs::HilbertSpaceSector) = basespace(hs.parent)

import Base.==
function ==(lhs ::HilbertSpaceSector{Q1}, rhs ::HilbertSpaceSector{Q2}) where {Q1, Q2}
  return (Q1 == Q2) && basespace(lhs) == basespace(rhs) && lhs.allowed_quantum_numbers == rhs.allowed_quantum_numbers
end


for fname in [:get_bitmask,
              :quantum_number_sectors,
              :get_quantum_number,
              :extract,
              :compress,
              :update,
              :get_state_index,
              :get_state]
  eval(:(@inline $fname(hss ::HilbertSpaceSector, args...;kwargs...) = $fname(hss.parent, args...; kwargs...)))
end

# function extract(hss ::HilbertSpaceSector{QN}, binrep ::U) where {QN, U <:Unsigned}
#   return extract(hss.parent, binrep)
# end
#
# function compress(hss ::HilbertSpaceSector{QN}, indexarray ::AbstractVector{I}; BR::DataType=UInt) where {QN, I<:Integer}
#   return compress(hss.parent, indexarray; BR=BR)
# end
#
# function update(hss ::HilbertSpaceSector, binrep ::U, isite ::Integer, new_state_index ::Integer) where {U<:Unsigned}
#   return update(hss.parent, binrep, isite, new_state_index)
# end
#
# function get_state_index(hss ::HilbertSpaceSector, binrep ::U, isite ::Integer) where {U<:Unsigned}
#   return get_state_index(hss.parent, binrep, isite)
# end
#
# function get_state(hss ::HilbertSpaceSector, binrep ::U, isite ::Integer) where {U<:Unsigned}
#   return get_state(hss.parent, binrep, isite)
# end
#
# @inline bitwidth(hss::HilbertSpaceSector) = bitwidth(hss.parent)

# import Base.iterate
# @inline function iterate(hss ::HilbertSpaceSector{QN}) where {QN}
#   @warn "Use of iterate for HilbertSpaceSector is deprecated"
#   error("Not implemented")
# end
