export DenseState

"""
    DenseState

Implemented with dense vector, accessed through binary representation
"""
mutable struct DenseState{S<:Number, QN, BR<:Unsigned} #<:AbstractDict{BR, S}
  hilbert_space_realization ::HilbertSpaceRealization{QN, BR}
  components ::Vector{S}

  function DenseState{S, BR}(hsr ::HilbertSpaceRealization{QN, BR}) where {S, QN, BR}
    return new{S, QN, BR}(hsr, zeros(S,dimension(hsr)))
  end

  function DenseState{S, BR}(hsr ::HilbertSpaceRealization{QN, BR}, binrep ::BR) where {S, QN, BR}
    components = zeros(S, dimension(hsr))
    i = hsr.basis_lookup[binrep]
    components[i] = one(S)
    return new{S, QN, BR}(hs, components)
  end

  function DenseState{S, BR}(hsr ::HilbertSpaceRealization{QN, BR}, components::Pair{BR, <:Number}...) where {S, QN, BR}
    compos = zeros(S, dimension(hsr))
    for (cf, cs) in components
      icf = hsr.basis_lookup[cf]
      compos[cf] = cs
    end
    return new{S, QN, BR}(hs, compos)
  end

  function DenseState{S, BR}(hsr ::HilbertSpaceRealization{QN, BR}, components::AbstractDict{BR, S2}) where {S, QN, BR, S2}
    # TODO bound checking
    compos = zeros(S, dimension(hsr))
    for (cf, cs) in components
      icf = hsr.basis_lookup[cf]
      compos[icf] = cs
    end
    return new{S, QN, BR}(hs, compos)
  end

  function DenseState{S, BR}(hsr ::HilbertSpaceRealization{QN, BR},
                             components::AbstractVector{S2}) where {S, QN, BR, S2<:Number}
    @boundscheck dimension(hsr) == length(components) || throw(ArgumentError("components should match dimension"))
    return new{S, QN, BR}(hs, components)
  end

  function DenseState(hsr ::HilbertSpaceRealization{QN, BR},
                      components::AbstractVector{S}) where {S<:Number, QN, BR}
    @boundscheck dimension(hsr) == length(components) || throw(ArgumentError("components should match dimension"))
    return new{S, QN, BR}(hs, components)
  end
end

import Base.==
function (==)(lhs ::DenseState{S1, QN, BR}, rhs::DenseState{S2, QN, BR}) where {S1, S2, QN, BR}
  if lhs.hilbert_space_realization !== rhs.hilbert_space_realization
    return false
  end
  return lhs.components == rhs.components
end


import Base.isapprox
function isapprox(lhs ::DenseState{S1, QN, BR}, rhs::DenseState{S2, QN, BR};
                  atol=sqrt(eps(Float64)), rtol=sqrt(eps(Float64))) where {S1, S2, QN, BR}
  if lhs.hilbert_space_realization !== rhs.hilbert_space_realization
    return false
  end
  return isapprox(lhs.components, rhs.components; atol=atol, rtol=rtol)
end


import Base.copy
function copy(arg ::DenseState{S, QN, BR}) where {S, QN, BR}
  return DenseState{S, QN, BR}(arg.hilbert_space_realization, copy(arg.components))
end


import Base.eltype
eltype(arg ::DenseState{S, QN, BR}) where {S, QN, BR} = S


import Base.getindex
function Base.getindex(iter ::DenseState{S, QN, BR}, binrep ::BR) ::S where {S, QN, BR}
  i = iter.hilbert_space_realization.basis_lookup[binrep]
  return iter.components[i]
end

import Base.setindex!
function Base.setindex!(iter ::DenseState{S, QN, BR}, value ::Number, binrep ::BR) where {S, QN, BR}
  i = iter.hilbert_space_realization.basis_lookup[binrep]
  iter.components[i] = value
  return iter
end

import Base.iterate
function Base.iterate(iter ::DenseState{S, QN, BR}, i::Int=1) ::Union{Nothing, Tuple{Pair{BR, S}, Int}} where {S, QN, BR}
  i > length(components) && return nothing
  return (iter.hilbert_space_realization.basis_list[i] => components[i], i+1)
end


import Base.real, Base.imag, Base.conj

real(arg ::DenseState{R, QN, BR}) where {R<:Real, QN, BR} = copy(arg)
imag(arg ::DenseState{R, QN, BR}) where {R<:Real, QN, BR} = DenseState{R, QN, BR}(arg.hilbert_space_realiation)
conj(arg ::DenseState{R, QN, BR}) where {R<:Real, QN, BR} = copy(arg)

function real(arg ::DenseState{Complex{R}, QN, BR}) where {R<:Real, QN, BR}
  return DenseState{R, QN, BR}(arg.hilbert_space_realization, real.(arg.components))
end

function Base.imag(arg ::DenseState{Complex{R}, QN, BR}) where {R<:Real, QN, BR}
  return DenseState{R, QN, BR}(arg.hilbert_space_realization, imag.(arg.components))
end

function Base.conj(arg ::DenseState{Complex{R}, QN, BR}) where {R<:Real, QN, BR}
  return DenseState{Complex{R}, QN, BR}(arg.hilbert_space_realization, conj.(arg.components))
end


import Base.+, Base.-, Base.*, Base./, Base.\
(+)(arg ::DenseState{S, QN, BR}) where {S, QN, BR} = copy(arg)

function (-)(arg ::DenseState{S, QN, BR}) where {S, QN, BR}
  return DenseState{S, QN, BR}(arg.hilbert_space_realization, -arg.components)
end

function (+)(lhs ::DenseState{S1, QN, BR}, rhs ::DenseState{S2, QN, BR}) where {S1, S2, QN, BR}
  if lhs.hilbert_space_realization !== rhs.hilbert_space_realization
    throw(ArgumentError("HilbertSpaceRealization of lhs and rhs of + do not match"))
  end
  S = promote_type(S1, S2)
  return DenseState{S, QN, BR}(lhs.hilbert_space_realization, lhs.components .+ rhs.components)
end

function (-)(lhs ::DenseState{S1, QN, BR}, rhs ::DenseState{S2, QN, BR}) where {S1, S2, QN, BR}
  if lhs.hilbert_space_realization !== rhs.hilbert_space_realization
    throw(ArgumentError("HilbertSpaceRealization of lhs and rhs of + do not match"))
  end
  S = promote_type(S1, S2)
  return DenseState{S, QN, BR}(lhs.hilbert_space_realization, lhs.components .- rhs.components)
end

function (*)(arg::DenseState{S, QN, BR}, v::T) where {S, T <:Number, QN, BR}
  S2 = promote_type(S, T)
  return DenseState{S2, QN, BR}(arg.hilbert_space_realization, arg.components .* v)
end

function (*)(v::T, arg::DenseState{S, QN, BR}) where {S, T <:Number, QN, BR}
  S2 = promote_type(S, T)
  return DenseState{S2, QN, BR}(arg.hilbert_space_realization, v .* arg.components)
end

function (/)(arg::DenseState{S, QN, BR}, v::T) where {S, T <:Number, QN, BR}
  return DenseState(arg.hilbert_space_realization, arg.components ./ v)
end

function (\)(v::T, arg::DenseState{S, QN, BR}) where {S, T <:Number, QN, BR}
  return DenseState(arg.hilbert_space_realization, v .\ arg.components)
end

import LinearAlgebra.norm
norm(arg ::DenseState{S, QN, BR}) where {S, QN, BR} = norm(arg.components)

import LinearAlgebra.normalize
function normalize(arg ::DenseState{S, QN, BR}) where {S, QN, BR}
  nv = norm(arg.components)
  components = arg.components ./ nv
  return DenseState(arg.hilbert_space_realization, components)
end

import LinearAlgebra.normalize!
function normalize!(arg ::DenseState{S, QN, BR}) where {S, QN, BR}
  nv = norm(arg.components)
  arg.components ./= nv
  return arg
end
