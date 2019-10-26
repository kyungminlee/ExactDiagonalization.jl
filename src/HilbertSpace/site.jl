export AbstractHilbertSpace
export State, Site
export bitwidth, get_state, dimension
export qntype
export quantum_number_sectors
export get_quantum_number
export compress
export get_state_index

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
@inline function get_state(site::Site, binrep::U) where {U<:Unsigned}
  return site.states[Int(binrep+1)]
end

@inline function compress(site ::Site, i ::Integer; BR::DataType=UInt)
  @boundscheck 1 <= i <= dimension(site) || throw(BoundsError("attempt to access a $(dimension(site))-state site at index $i"))
  return BR(i-1)
end

@inline function get_state_index(site::Site, binrep::U) where {U<:Unsigned}
  i = Int(binrep+1)
  @boundscheck 0 <= i <= dimension(site) || throw(BoundsError("attempt to access a $(dimension(site))-state site at index $i"))
  return i
end

@inline function quantum_number_sectors(site ::Site{QN})::Vector{QN} where QN
  return sort(collect(Set([state.quantum_number for state in site.states])))
end

@inline function get_quantum_number(site ::Site{QN}, i ::Integer)::QN where QN
  return site.states[i].quantum_number
end


import Base.keys
keys(site::Site{QN}) where QN = 1:dimension(site)
# import Base.iterate
# @inline iterate(site::Site{QN}) where {QN} = Base.iterate(1:length(site.states))
# @inline iterate(site::Site{QN}, i) where QN = Base.iterate(1:length(site.states), i)

# import Base.length
# @inline length(site::Site{QN}) where QN = length(site.states)
