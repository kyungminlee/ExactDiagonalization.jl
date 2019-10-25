export SparseState
export choptol!
export bintype

using LinearAlgebra

"""
    struct SparseState{Scalar<:Number, BR}

Represents a row vector. Free.
"""
mutable struct SparseState{Scalar<:Number, BR}
  hilbert_space ::HilbertSpace
  components ::Dict{BR, Scalar}
  function SparseState{Scalar, BR}(hs ::HilbertSpace) where {Scalar, BR}
    return new{Scalar, BR}(hs, Dict{BR, Scalar}())
  end

  function SparseState{Scalar, BR}(hs ::HilbertSpace, components ::Dict{BR, Scalar}) where {Scalar, BR}
    return new{Scalar, BR}(hs, components)
  end

  function SparseState(hs ::HilbertSpace, components ::Dict{BR, Scalar}) where {Scalar, BR}
    return new{Scalar, BR}(hs, components)
  end

  function SparseState{Scalar, BR}(hs ::HilbertSpace, binrep ::BR) where {Scalar, BR}
    return new{Scalar, BR}(hs, Dict{BR, Scalar}(binrep => one(Scalar)))
  end

  function SparseState{Scalar, BR}(hs ::HilbertSpace, components::Pair{BR2, <:Number}...) where {Scalar, BR, BR2<:Unsigned}
    return new{Scalar, BR}(hs, Dict{BR, Scalar}(components))
  end

  function SparseState{Scalar, BR}(hs ::HilbertSpace, components ::AbstractDict{BR, S2}) where {Scalar, BR, S2}
    return new{Scalar, BR}(hs, components)
  end
end

import Base.getindex, Base.setindex!

function Base.getindex(state ::SparseState{Scalar, BR}, basis ::BR2) where {Scalar, BR, BR2 <:Unsigned}
  # TODO: check hilbert space
  return get(state.components, basis, zero(Scalar))
end

function Base.setindex!(state ::SparseState{Scalar, BR}, value ::S, basis ::BR2) where {Scalar, BR, S<:Number, BR2 <:Unsigned}
  # TODO: check hilbert space
  Base.setindex!(state.components, value, basis)
  return state
end


scalartype(::SparseState{Scalar, BR}) where {Scalar, BR} = Scalar
scalartype(::Type{SparseState{Scalar, BR}}) where {Scalar, BR} = Scalar
bintype(::SparseState{Scalar, BR}) where {Scalar, BR} = BR
bintype(::Type{SparseState{Scalar, BR}}) where {Scalar, BR} = BR


import Base.==
function (==)(lhs ::SparseState{S1, BR}, rhs::SparseState{S2, BR}) where {S1, S2, BR}
  return (lhs.hilbert_space === rhs.hilbert_space) && (lhs.components == rhs.components)
end

import Base.isapprox
function isapprox(lhs ::SparseState{S1, BR}, rhs::SparseState{S2, BR}; atol=sqrt(eps(Float64)), rtol=sqrt(eps(Float64))) where {S1, S2, BR}
  if lhs.hilbert_space !== rhs.hilbert_space
    return false
  end

  all_keys = union(keys(lhs.components), keys(rhs.components))
  for k in all_keys
    lv = get(lhs.components, k, zero(S1))
    rv = get(rhs.components, k, zero(S2))
    if ! isapprox(lv, rv; atol=atol, rtol=rtol)
      return false
    end
  end

  return true
end

import Base.copy
function copy(arg ::SparseState{S, BR}) where {S, BR}
  return SparseState{S, BR}(arg.hilbert_space, copy(arg.components))
end

import Base.real, Base.imag, Base.conj
real(arg ::SparseState{R, BR}) where {R<:Real, BR} = copy(arg)
imag(arg ::SparseState{R, BR}) where {R<:Real, BR} = SparseState{R, BR}(arg.hilbert_space)
conj(arg ::SparseState{R, BR}) where {R<:Real, BR} = copy(arg)

function real(arg ::SparseState{Complex{R}, BR}) where {R<:Real, BR}
  return SparseState{R, BR}(arg.hilbert_space,
                            Dict{BR, R}((k, real(v)) for (k, v) in arg.components))
end

function imag(arg ::SparseState{Complex{R}, BR}) where {R<:Real, BR}
  return SparseState{R, BR}(arg.hilbert_space,
                            Dict{BR, R}((k, imag(v)) for (k, v) in arg.components))
end

function conj(arg ::SparseState{Complex{R}, BR}) where {R<:Real, BR}
  return SparseState{Complex{R}, BR}(arg.hilbert_space,
                                     Dict{BR, Complex{R}}((k, conj(v)) for (k, v) in arg.components))
end

import Base.-, Base.+, Base.*, Base./, Base.\

(+)(arg ::SparseState{S, BR}) where {S, BR} = copy(arg)

function (-)(arg ::SparseState{S, BR}) where {S, BR}
  return SparseState{S, BR}(arg.hilbert_space,
                            Dict{BR, S}((k, -v) for (k, v) in arg.components))
end

function (+)(lhs ::SparseState{S1, BR}, rhs ::SparseState{S2, BR}) where {S1, S2, BR}
  lhs.hilbert_space === rhs.hilbert_space || throw(ArgumentError("Hilbert spaces of lhs and rhs of + should match"))

  S3 = promote_type(S1, S2)
  components = Dict{BR, S3}(lhs.components)
  for (b, v) in rhs.components
    components[b] = get(components, b, zero(S3)) + v
  end
  return SparseState{S3, BR}(lhs.hilbert_space, components)
end

function (-)(lhs ::SparseState{S1, BR}, rhs ::SparseState{S2, BR}) where {S1, S2, BR}
  lhs.hilbert_space === rhs.hilbert_space || throw(ArgumentError("Hilbert spaces of lhs and rhs of + should match"))

  S3 = promote_type(S1, S2)
  components = Dict{BR, S3}(lhs.components)
  for (b, v) in rhs.components
    components[b] = get(components, b, zero(S3)) - v
  end
  return SparseState{S3, BR}(lhs.hilbert_space, components)
end

function (*)(lhs ::SparseState{S1, BR}, rhs ::S2) where {S1, S2<:Number, BR}
  return SparseState(lhs.hilbert_space,
                     Dict(k => v * rhs for (k, v) in lhs.components))
end

function (*)(lhs ::S1, rhs ::SparseState{S2, BR}) where {S1<:Number, S2<:Number, BR}
  return SparseState(rhs.hilbert_space,
                     Dict(k => lhs * v for (k, v) in rhs.components))
end

function (/)(lhs ::SparseState{S1, BR}, rhs ::S2) where {S1, S2<:Number, BR}
  return SparseState(lhs.hilbert_space,
                     Dict(k => v / rhs for (k, v) in lhs.components))
end

function (\)(lhs ::S1, rhs ::SparseState{S2, BR}) where {S1<:Number, S2<:Number, BR}
  return SparseState(rhs.hilbert_space,
                     Dict(k => lhs \ v for (k, v) in rhs.components))
end

import Base.convert
function convert(type ::Type{SparseState{S1, BR}}, obj::SparseState{S2, BR}) where {S1, S2, BR}
  return SparseState{S1, BR}(obj.hilbert_space, Dict{BR, S1}(obj.components))
end

import Base.eltype
Base.eltype(::Type{SparseState{Scalar, BR}}) where {Scalar, BR} = Pair{BR, Scalar}

import Base.isempty
isempty(psi::SparseState{S, BR}) where {S, BR} = isempty(psi.components)

import Base.length
length(psi::SparseState{S, BR}) where {S, BR} = length(psi.components)

import Base.iterate
function Base.iterate(iter ::SparseState{S, BR}) where {S, BR}
  return Base.iterate(iter.components)
end

function Base.iterate(iter ::SparseState{S, BR}, i) where {S, BR}
  return Base.iterate(iter.components, i)
end

function choptol!(arg ::SparseState{S1, BR}, tol::Real) where {S1, BR}
  to_delete = [k for (k, v) in arg.components if isapprox(v, 0; atol=tol)]
  for k in to_delete
    delete!(arg.components, k)
  end
end

import LinearAlgebra.norm
function norm(arg ::SparseState{S1, BR}) where {S1, BR}
  if isempty(arg.components)
    return zero(real(S1))
  else
    return norm(values(arg.components))
  end
end

import LinearAlgebra.normalize
function normalize(arg ::SparseState{S1, BR}) where {S1, BR}
  norm_val = norm(arg)
  S2 = promote_type(typeof(norm_val), S1)
  components = Dict{BR, S2}(k => v/norm_val for (k, v) in arg.components)
  return SparseState{S2, BR}(arg.hilbert_space, components)
end

import LinearAlgebra.normalize!
function normalize!(arg ::SparseState{S1, BR}) where {S1, BR}
  norm_val = norm(arg)
  for (k, v) in arg.components
    arg[k] = v / norm_val
  end
  arg
end


# TODO Broadcasting
