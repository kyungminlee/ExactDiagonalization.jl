export HilbertSpaceSector
export scalartype
export qntype
export basespace

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

scalartype(arg ::HilbertSpaceSector{QN}) where QN = Bool
scalartype(arg ::Type{HilbertSpaceSector{QN}}) where QN = Bool

import Base.valtype
valtype(arg ::HilbertSpaceSector{QN}) where QN = Bool
valtype(arg ::Type{HilbertSpaceSector{QN}}) where QN = Bool

qntype(arg ::HilbertSpaceSector{QN}) where QN = QN
qntype(arg ::Type{HilbertSpaceSector{QN}}) where QN = QN

basespace(hs::HilbertSpaceSector{QN}) where QN = basespace(hs.parent) ::HilbertSpace{QN}

import Base.==
function ==(lhs ::HilbertSpaceSector{Q1}, rhs ::HilbertSpaceSector{Q2}) where {Q1, Q2}
  return basespace(lhs) == basespace(rhs) && lhs.allowed_quantum_numbers == rhs.allowed_quantum_numbers
end


for fname in [:get_bitmask,
              :quantum_number_sectors,
              :get_quantum_number,
              :extract,
              :compress,
              :update,
              :get_state_index,
              :get_state]
  @eval begin
    @inline $fname(hss ::HilbertSpaceSector, args...;kwargs...) = $fname(hss.parent, args...; kwargs...)
  end
end
