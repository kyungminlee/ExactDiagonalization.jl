export AbstractParticle
export Fermion, HardcoreBoson, Boson

export getspeciessymb
export maxoccupancy
export exchangesign


abstract type AbstractParticle end

struct Fermion{Species} <: AbstractParticle end # Species is Symbol
struct HardcoreBoson{Species} <: AbstractParticle end
struct Boson{Species, Max} <: AbstractParticle end

getspeciessymb(p::Type{<:AbstractParticle}) = p.parameters[1]

# Default is 1 for different types
exchangesign(::Type{T1}, ::Type{T2}) where {T1 <:AbstractParticle, T2 <:AbstractParticle} = 1
exchangesign(::Type{T}, ::Type{T}) where T <: Fermion = -1
exchangesign(::T1, ::T2) where {T1 <:AbstractParticle, T2 <:AbstractParticle} = exchangesign(T1, T2)

maxoccupancy(::Type{<:Fermion}) = 1
maxoccupancy(::Type{<:HardcoreBoson}) = 1
maxoccupancy(::Type{Boson{S, M}}) where {S, M} = M::Int

export ParticleSiteType
struct ParticleSiteType{ParticleType<:AbstractParticle} <: AbstractSiteType end

maxoccupancy(::Type{ParticleSiteType{P}}) where P <: AbstractParticle = maxoccupancy(P)



#function splitparticle(hs::HilbertSpace, binrep::U) where {U<:Unsigned}
#end

#
# abstract type AbstractHilbertSpace end
# abstract type AbstractSite{QN} <: AbstractHilbertSpace end
#
#
# struct ParticleState{QN<:AbstractQuantumNumber}
#   name ::String
#   quantum_number::QN
#
#   ParticleState(name ::AbstractString) = new{Int}(name, 0)
#   ParticleState(name ::AbstractString, quantum_number ::QN) where {QN} = new{QN}(name, quantum_number)
#   ParticleState{QN}(name ::AbstractString, quantum_number ::QN) where {QN} = new{QN}(name, quantum_number)
# end
#
#
#
# struct ParticleSite{P<:AbstractParticle, QN<:AbstractQuantumNumber} <: AbstractSite{QN}
#   states ::Vector{ParticleState{QN}}
#   function ParticleSite(::Type{P}, states::AbstractVector{ParticleState{QN}}) where {P<:AbstractParticle, QN}
#     if length(states) != maxoccupancy(P)+1
#       throw(ArgumentError("number of states $(length(states)) does not match the maximum occupancy of the particle $P $(maxoccupancy(P))"))
#     end
#     new{P, QN}(states)
#   end
# end
#
#
#
# getparticlespecies(p::ParticleSite{P}) where P = P
# dimension(site::ParticleSite{P}) where P = maxoccupancy(P)+1
# bitwidth(site::ParticleSite) = Int(ceil(log2(dimension(site))))
# get_state(site::ParticleSite, binrep::U) where {U<:Unsigned} = site.states[Int(binrep)+1]
#
# @inline function compress(site ::ParticleSite, state_index ::Integer, binary_type::Type{BR}=UInt) where {BR<:Unsigned}
#   @boundscheck 1 <= state_index <= dimension(site) || throw(BoundsError("attempt to access a $(dimension(site))-state site at index $state_index"))
#   return BR(state_index-1)
# end
#
# @inline function get_state_index(site::ParticleSite, binrep::U) where {U<:Unsigned}
#   i = Int(binrep+1)
#   @boundscheck 1 <= i <= dimension(site) || throw(BoundsError("attempt to access a $(dimension(site))-state site at index $i"))
#   return i
# end
#
# function quantum_number_sectors(site ::ParticleSite{P, QN})::Vector{QN} where {P, QN}
#   return sort(collect(Set([state.quantum_number for state in site.states])))
# end
#
# function get_quantum_number(site ::ParticleSite{P, QN}, state_index ::Integer)::QN where {P, QN}
#   return site.states[state_index].quantum_number
# end
#
# import Base.keys
# keys(site::ParticleSite) = Base.OneTo(dimension(site))
#
#
#
# struct ParticleHilbertSpace{S<:Tuple{Vararg{ParticleSite{<:AbstractParticle, <:AbstractQuantumNumber}}}, QN} <: AbstractHilbertSpace
#   sites :: S
#   bitwidths :: Vector{Int}
#   bitoffsets :: Vector{Int}
#
#   function ParticleHilbertSpace(sites::Vararg{ParticleSite{<:AbstractParticle, QN}}) where QN<:AbstractQuantumNumber
#     S = typeof(sites)
#     bitwidths = map(bitwidth, collect(sites))
#     bitoffsets = Int[0, cumsum(bitwidths)...]
#     new{S, QN}(sites, bitwidths, bitoffsets)
#   end
# end
#
#
# scalartype(arg::Type{<:ParticleHilbertSpace}) = Bool
# scalartype(arg::ParticleHilbertSpace) = Bool
#
# valtype(arg ::Type{<:ParticleHilbertSpace}) = Bool
# valtype(arg ::ParticleHilbertSpace) = Bool
#
# qntype(arg ::Type{ParticleHilbertSpace{S, QN}}) where {S, QN} = QN
# qntype(arg ::ParticleHilbertSpace{S, QN}) where {S, QN} = QN
#
# basespace(hs::ParticleHilbertSpace) = hs
#
# bitwidth(hs::ParticleHilbertSpace) = hs.bitoffsets[end]
#
# function get_bitmask(hs ::ParticleHilbertSpace, isite ::Integer, binary_type ::Type{BR}=UInt) ::BR where {BR<:Unsigned}
#   return make_bitmask(hs.bitoffsets[isite+1], hs.bitoffsets[isite], BR)
# end
#
# function quantum_number_sectors(hs ::ParticleHilbertSpace{S, QN})::Vector{QN} where {S, QN}
#   qns = Set{QN}([zero(QN)])
#   for site in hs.sites
#     qns_next = Set{QN}()
#     for state in site.states, q in qns
#       push!(qns_next, q .+ state.quantum_number)
#     end
#     qns = qns_next
#   end
#   return sort(collect(qns))
# end
#
#
#
#
# function get_quantum_number(hs ::ParticleHilbertSpace{S, QN}, indexarray ::AbstractArray{I, 1}) where {S, QN, I <:Integer}
#     sum(
#       site.states[indexarray[isite]].quantum_number
#       for (isite, site) in enumerate(hs.sites)
#     )
# end
#
#
# function extract(hs ::HilbertSpace{QN}, binrep ::BR) ::CartesianIndex where {QN, BR<:Unsigned}
#   out = Int[]
#   for (isite, site) in enumerate(hs.sites)
#     @inbounds mask = make_bitmask(hs.bitwidths[isite], BR)
#     index = Int(binrep & mask) + 1
#     @boundscheck if !(1 <= index <= length(site.states))
#       throw(BoundsError(1:length(site.states), index))
#     end
#     push!(out, index)
#     binrep = binrep >> hs.bitwidths[isite]
#   end
#   return CartesianIndex(out...)
# end
#
# function compress(hs ::ParticleHilbertSpace, indexarray ::CartesianIndex, binary_type ::Type{BR}=UInt) where {BR<:Unsigned}
#   if length(indexarray) != length(hs.sites)
#     throw(ArgumentError("length of indexarray should be the number of sites"))
#   end
#   binrep = zero(BR)
#   for (isite, site) in enumerate(hs.sites)
#     @boundscheck if !(1 <= indexarray[isite] <= dimension(site))
#       throw(BoundsError(1:dimension(site), indexarray[isite]))
#     end
#     binrep |= (BR(indexarray[isite] - 1) << hs.bitoffsets[isite] )
#   end
#   return binrep
# end
#
#
#
# """
#     update(hs, binrep, isite, new_state_index)
#
# Update the binary representation of a basis state
# by changing the state at site `isite` to a new local state specified by
# `new_state_index`.
# """
# @inline function update(hs ::ParticleHilbertSpace, binrep ::BR, isite ::Integer, new_state_index ::Integer) where {BR<:Unsigned}
#   @boundscheck if !(1 <= new_state_index <= dimension(hs.sites[isite]))
#     throw(BoundsError(1:dimension(hs.sites[isite]), new_state_index))
#   end
#   mask = make_bitmask(hs.bitoffsets[isite+1], hs.bitoffsets[isite], BR)
#   return (binrep & (~mask)) | (BR(new_state_index-1) << hs.bitoffsets[isite])
# end
#
#
# function get_state_index(hs ::ParticleHilbertSpace, binrep ::BR, isite ::Integer) where {BR<:Unsigned}
#   return Int( ( binrep >> hs.bitoffsets[isite] ) & make_bitmask(hs.bitwidths[isite], BR) ) + 1
# end
#
#
# function get_state(hs ::ParticleHilbertSpace, binrep ::BR, isite ::Integer) where {BR<:Unsigned}
#   return hs.sites[isite].states[get_state_index(hs, binrep, isite)]
# end
#
#
# import Base.keys
# function keys(hs ::ParticleHilbertSpace)
#   return CartesianIndices(((1:length(site.states) for site in hs.sites)...,))
# end
#
#
#
#
#
# """
#     HilbertSpaceSector{QN}
#
# Hilbert space sector.
# """
# struct ParticleHilbertSpaceSector{QN} <: AbstractHilbertSpace
#   parent ::ParticleHilbertSpace{QN}
#   allowed_quantum_numbers ::Set{QN}
#
#   function ParticleHilbertSpaceSector(parent ::ParticleHilbertSpace{QN}) where QN
#     sectors = quantum_number_sectors(parent)
#     new{QN}(parent, Set(sectors))
#   end
#
#   function HilbertSpaceSector(parent ::HilbertSpace{QN}, allowed::QN) where QN
#     sectors = Set{QN}(quantum_number_sectors(parent))
#     new{QN}(parent, intersect(sectors, Set([allowed])))
#   end
#
#   function HilbertSpaceSector(parent ::HilbertSpace{QN},
#                               allowed::Union{AbstractSet{QN}, AbstractVector{QN}}) where QN
#     sectors = Set{QN}(quantum_number_sectors(parent))
#     new{QN}(parent, intersect(sectors, Set(allowed)))
#   end
# end
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# Electron = Fermion{:e}
# Alpha = HardcoreBoson{:α}
#
# electron_site = ParticleSite(Electron, [ParticleState("Empty", 0), ParticleState("Occupied", 1)])
#
# hs = ParticleHilbertSpace(electron_site, electron_site, electron_site, electron_site)
