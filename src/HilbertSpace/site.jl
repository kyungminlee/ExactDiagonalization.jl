export AbstractHilbertSpace
export State, Site
export bitwidth, get_state, dimension
export qntype

using StaticArrays

abstract type AbstractHilbertSpace end

## TODO: Think about this
AbstractQuantumNumber = Union{Int, SVector{N, Int} where N}

@inline qntype(arg ::T) where T = qntype(T)

"""
    State{QN}

State with quantum number type `QN`.

# Examples
```jldoctest
julia> using ExactDiagonalization, StaticArrays

julia> up = State{Int}("Up", 1)
State{Int64}("Up", 1)

julia> State("Dn", SVector{2, Int}([-1, 1]))
State{SArray{Tuple{2},Int64,1,2}}("Dn", [-1, 1])
```
"""
struct State{QN<:AbstractQuantumNumber}
  name ::String
  quantum_number ::QN

  State(name ::AbstractString) = new{Int}(name, 0)
  State(name ::AbstractString, quantum_number ::QN) where {QN} = new{QN}(name, quantum_number)
  State{QN}(name ::AbstractString, quantum_number ::QN) where {QN} = new{QN}(name, quantum_number)
end

import Base.==
function ==(lhs ::State{Q1}, rhs ::State{Q2}) where {Q1, Q2}
  return (Q1 == Q2) && (lhs.name == rhs.name) && (lhs.quantum_number == rhs.quantum_number)
end

@inline qntype(::Type{State{QN}}) where QN = QN

"""
    Site{QN}

A site with quantum number type `QN`.

# Examples
```jldoctest
julia> using ExactDiagonalization

julia> up = State{Int}("Up", 1); dn = State{Int}("Dn", -1);

julia> Site([up, dn])
Site{Int64}(State{Int64}[State{Int64}("Up", 1), State{Int64}("Dn", -1)])
```
"""
struct Site{QN<:AbstractQuantumNumber} <: AbstractHilbertSpace
  states ::Vector{State{QN}}

  Site(states ::AbstractArray{State{QN}, 1}) where QN = new{QN}(states)
  Site{QN}(states ::AbstractArray{State{QN}, 1}) where QN = new{QN}(states)
end

@inline qntype(::Type{Site{QN}}) where QN = QN

import Base.==
function ==(lhs ::Site{Q1}, rhs ::Site{Q2}) where {Q1, Q2}
  return (Q1 == Q2) && (lhs.states == rhs.states)
end

"""
    bitwidth(site ::Site)

Number of bits necessary to represent the states of the given site.
"""
bitwidth(site ::Site) ::Int = Int(ceil(log2(length(site.states))))

"""
    dimension(site ::Site)

Hilbert space dimension of a given site ( = number of states).
"""
dimension(site ::Site) ::Int = length(site.states)


"""
    get_state(site ::Site{QN}, binrep ::BR) where {QN, BR<:Unsigned}

Returns the state of `site` represented by the bits `binrep`.
"""
function get_state(site ::Site{QN}, binrep ::BR) where {QN, BR<:Unsigned}
  return site.states[Int(binrep) + 1]
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
