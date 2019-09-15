export State, Site
export bitwidth, get_state, dimension

using StaticArrays

## TODO: Think about this
AbstractQuantumNumber = Union{Int, SVector{N, Int} where N}

import Base.==

struct State{QN<:AbstractQuantumNumber}
  name ::String
  quantum_number ::QN
  
  State(name ::AbstractString) = new{Int}(name, 0)
  State(name ::AbstractString, quantum_number ::QN) where {QN} = new{QN}(name, quantum_number)
  State{QN}(name ::AbstractString, quantum_number ::QN) where {QN} = new{QN}(name, quantum_number)
end

function ==(lhs ::State{Q1}, rhs ::State{Q2}) where {Q1, Q2}
  return (Q1 == Q2) && (lhs.name == rhs.name) && (lhs.quantum_number == rhs.quantum_number)
end

struct Site{QN<:AbstractQuantumNumber}
  states ::Vector{State{QN}}

  Site(states ::AbstractArray{State{QN}, 1}) where QN = new{QN}(states)
  Site{QN}(states ::AbstractArray{State{QN}, 1}) where QN = new{QN}(states)
end

function ==(lhs ::Site{Q1}, rhs ::Site{Q2}) where {Q1, Q2}
  return (Q1 == Q2) && (lhs.states == rhs.states)
end

bitwidth(site ::Site) ::Int = Int(ceil(log2(length(site.states))))
dimension(site ::Site) ::Int = length(site.states)

function get_state(site ::Site{QN}, binrep ::BR) where {QN, BR<:Unsigned}
  @assert let
    nd = bitwidth(site)
    mask = make_bitmask(nd; dtype=BR)
    mask & binrep == binrep
  end
  @assert binrep < length(site.states)
  return site.states[binrep + 1]
end

# function get_state(site ::Site{QN}, bitrep ::BitArray{1}) where QN
#   @assert let
#     nd = bitwidth(site)
#     length(bitrep == nd)
#   end
#   @assert length(bitrep.chunks) == 1
#   binrep = bitrep.chunks[1]
#   @assert binrep < length(site.states)
#   return site.states[binrep + 1]
# end
