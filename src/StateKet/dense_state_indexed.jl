export DenseState

mutable struct DenseState{Scalar<:Number, QN, BR}
  hilbert_space_realization ::HilbertSpaceRealization{QN, BR}
  components ::Vector{Scalar}

  function DenseState{Scalar, BR}(hsr ::HilbertSpaceRealization{QN, BR}) where {Scalar, QN, BR}
    return new{Scalar, QN, BR}(hsr, zeros(Scalar,dimension(hsr)))
  end

  function DenseState{Scalar, BR}(hsr ::HilbertSpaceRealization{QN, BR}, binrep ::BR) where {Scalar, QN, BR}
    components = zeros(Scalar, dimension(hsr))
    i = hsr.basis_lookup[binrep]
    components[i] = one(Scalar)
    return new{Scalar, BR}(hs, components)
  end

  function DenseState{Scalar, BR}(hsr ::HilbertSpaceRealization{QN, BR}, components::Pair{BR, <:Number}...) where {Scalar, QN, BR}
    compos = zeros(Scalar, dimension(hsr))
    for (cf, cs) in components
      icf = hsr.basis_lookup[cf]
      compos[cf] = cs
    end
    return new{Scalar, BR}(hs, compos)
  end

  function DenseState{Scalar, BR}(hsr ::HilbertSpaceRealization{QN, BR}, components::AbstractDict{BR, S2}) where {Scalar, QN, BR, S2}
    # TODO bound checking
    compos = zeros(Scalar, dimension(hsr))
    for (cf, cs) in components
      icf = hsr.basis_lookup[cf]
      compos[cf] = cs
    end
    return new{Scalar, BR}(hs, compos)
  end

  function DenseState{Scalar, BR}(hsr ::HilbertSpaceRealization{QN, BR}, components::AbstractVector{S2}) where {Scalar, QN, BR, S2}
    @boundscheck dimension(hsr) == length(components) || throw(ArgumentError("components should match dimension"))
    return new{Scalar, BR}(hs, compos)
  end

  function DenseState(hsr ::HilbertSpaceRealization{QN, BR}, components::AbstractVector{Scalar}) where {Scalar<:Number, QN, BR}
    @boundscheck dimension(hsr) == length(components) || throw(ArgumentError("components should match dimension"))
    return new{Scalar, BR}(hs, compos)
  end
end


import Base.==
function (==))(lhs ::DenseState{S1, QN, BR}, rhs::DenseState{S2, QN, BR}; atol=sqrt(eps(Float64)), rtol=sqrt(eps(Float64))) where {S1, S2, BR, QN}
  if lhs.hilbert_space_realization !== rhs.hilbert_space_realization
    return false
  end
  return lhs.components == rhs.components
end


import Base.isapprox
function isapprox(lhs ::DenseState{S1, QN, BR}, rhs::DenseState{S2, QN, BR}; atol=sqrt(eps(Float64)), rtol=sqrt(eps(Float64))) where {S1, S2, BR, QN}
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
function Base.getindex(iter::DenseState{S, QN, BR}, binrep::BR) ::S where {S, QN, BR}
  i = iter.hilbert_space_realization.basis_lookup[binrep]
  return iter.components[i]
end

import Base.setindex!
function Base.setindex!(iter::DenseState{S, QN, BR}, value, binrep::BR) where {S, QN, BR}
  i = iter.hilbert_space_realization.basis_lookup[binrep]
  iter.components[i] = value
  return iter
end

import Base.iterate
function Base.iterate(iter ::DenseState{S, QN, BR}, i::Int=1) ::Union{Nothing, Tuple{Tuple{BR, S}, Int}}
  i > length(components) && return nothing
  return ((iter.hilbert_space_realization.basis_list[i], components[i]), i+1)
end



import Base.real, Base.imag, Base.conj
function Base.real(arg ::DenseState{S, QN, BR}) where {S, QN, BR}
  return DenseState{real(S), QN, BR}(arg.hilbert_space_realization, real.(arg.components))
end
function Base.imag(arg ::DenseState{S, QN, BR}) where {S, QN, BR}
  return DenseState{real(S), QN, BR}(arg.hilbert_space_realization, imag.(arg.components))
end
function Base.conj(arg ::DenseState{S, QN, BR}) where {S, QN, BR}
  return DenseState{S, QN, BR}(arg.hilbert_space_realization, conj.(arg.components))
end


import Base.+, Base.-, Base.*, Base./
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
  S2 = promote_type(S, T)
  return DenseState{S2, QN, BR}(arg.hilbert_space_realization, arg.components ./ v)
end

function (\)(v::T, arg::DenseState{S, QN, BR}) where {S, T <:Number, QN, BR}
  S2 = promote_type(S, T)
  return DenseState{S2, QN, BR}(arg.hilbert_space_realization, v .\ arg.components)
end
