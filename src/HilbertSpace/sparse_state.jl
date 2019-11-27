export SparseState
export choptol!
export bintype

using LinearAlgebra

"""
    struct SparseState{Scalar<:Number, BR}

Represents a row vector. Free.
"""
mutable struct SparseState{Scalar<:Number, BR}
  components ::Dict{BR, Scalar}
  function SparseState{Scalar, BR}() where {Scalar, BR}
    return new{Scalar, BR}(Dict{BR, Scalar}())
  end

  function SparseState{Scalar, BR}(components ::Dict{BR, Scalar}) where {Scalar, BR}
    return new{Scalar, BR}(components)
  end

  function SparseState(components ::Dict{BR, Scalar}) where {Scalar, BR}
    return new{Scalar, BR}(components)
  end

  function SparseState{Scalar, BR}(binrep ::BR) where {Scalar, BR}
    return new{Scalar, BR}(Dict{BR, Scalar}(binrep => one(Scalar)))
  end

  function SparseState{Scalar, BR}(components::Pair{BR2, <:Number}...) where {Scalar, BR, BR2<:Unsigned}
    return new{Scalar, BR}(Dict{BR, Scalar}(components))
  end

  function SparseState{Scalar, BR}(components ::AbstractDict{BR, S2}) where {Scalar, BR, S2}
    return new{Scalar, BR}(components)
  end
end

import Base.getindex, Base.setindex!

function Base.getindex(state ::SparseState{Scalar, BR}, basis ::BR2) where {Scalar, BR, BR2 <:Unsigned}
  return get(state.components, basis, zero(Scalar))
end

function Base.setindex!(state ::SparseState{Scalar, BR}, value ::S, basis ::BR2) where {Scalar, BR, S<:Number, BR2 <:Unsigned}
  Base.setindex!(state.components, value, basis)
  return state
end

scalartype(::SparseState{Scalar, BR}) where {Scalar, BR} = Scalar
scalartype(::Type{SparseState{Scalar, BR}}) where {Scalar, BR} = Scalar

import Base.valtype
valtype(::SparseState{Scalar, BR}) where {Scalar, BR} = Scalar
valtype(::Type{SparseState{Scalar, BR}}) where {Scalar, BR} = Scalar

bintype(::SparseState{Scalar, BR}) where {Scalar, BR} = BR
bintype(::Type{SparseState{Scalar, BR}}) where {Scalar, BR} = BR

import Base.==
function (==)(lhs ::SparseState{S1, BR}, rhs::SparseState{S2, BR}) where {S1, S2, BR}
  return lhs.components == rhs.components
end

import Base.isapprox
function isapprox(lhs ::SparseState{S1, BR}, rhs::SparseState{S2, BR}; atol=sqrt(eps(Float64)), rtol=sqrt(eps(Float64))) where {S1, S2, BR}
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
  return SparseState{S, BR}(copy(arg.components))
end

import Base.real, Base.imag, Base.conj
real(arg ::SparseState{R, BR}) where {R<:Real, BR} = copy(arg)
imag(arg ::SparseState{R, BR}) where {R<:Real, BR} = SparseState{R, BR}()
conj(arg ::SparseState{R, BR}) where {R<:Real, BR} = copy(arg)

function real(arg ::SparseState{Complex{R}, BR}) where {R<:Real, BR}
  return SparseState{R, BR}(Dict{BR, R}((k, real(v)) for (k, v) in arg.components))
end

function imag(arg ::SparseState{Complex{R}, BR}) where {R<:Real, BR}
  return SparseState{R, BR}(Dict{BR, R}((k, imag(v)) for (k, v) in arg.components))
end

function conj(arg ::SparseState{Complex{R}, BR}) where {R<:Real, BR}
  return SparseState{Complex{R}, BR}(Dict{BR, Complex{R}}((k, conj(v)) for (k, v) in arg.components))
end

import Base.-, Base.+, Base.*, Base./, Base.\

(+)(arg ::SparseState{S, BR}) where {S, BR} = copy(arg)

function (-)(arg ::SparseState{S, BR}) where {S, BR}
  return SparseState{S, BR}(Dict{BR, S}((k, -v) for (k, v) in arg.components))
end

function (+)(lhs ::SparseState{S1, BR}, rhs ::SparseState{S2, BR}) where {S1, S2, BR}
  S3 = promote_type(S1, S2)
  components = Dict{BR, S3}(lhs.components)
  for (b, v) in rhs.components
    components[b] = get(components, b, zero(S3)) + v
  end
  return SparseState{S3, BR}(components)
end

function (-)(lhs ::SparseState{S1, BR}, rhs ::SparseState{S2, BR}) where {S1, S2, BR}
  S3 = promote_type(S1, S2)
  components = Dict{BR, S3}(lhs.components)
  for (b, v) in rhs.components
    components[b] = get(components, b, zero(S3)) - v
  end
  return SparseState{S3, BR}(components)
end

function (*)(lhs ::SparseState{S1, BR}, rhs ::S2) where {S1, S2<:Number, BR}
  return SparseState(Dict(k => v * rhs for (k, v) in lhs.components))
end

function (*)(lhs ::S1, rhs ::SparseState{S2, BR}) where {S1<:Number, S2<:Number, BR}
  return SparseState(Dict(k => lhs * v for (k, v) in rhs.components))
end

function (/)(lhs ::SparseState{S1, BR}, rhs ::S2) where {S1, S2<:Number, BR}
  return SparseState(Dict(k => v / rhs for (k, v) in lhs.components))
end

function (\)(lhs ::S1, rhs ::SparseState{S2, BR}) where {S1<:Number, S2<:Number, BR}
  return SparseState(Dict(k => lhs \ v for (k, v) in rhs.components))
end

import Base.convert
function convert(type ::Type{SparseState{S1, BR}}, obj::SparseState{S2, BR}) where {S1, S2, BR}
  return SparseState{S1, BR}(Dict{BR, S1}(obj.components))
end

import Base.eltype
Base.eltype(::Type{SparseState{Scalar, BR}}) where {Scalar, BR} = Pair{BR, Scalar}

import Base.isempty
Base.isempty(psi::SparseState{S, BR}) where {S, BR} = isempty(psi.components)

import Base.length
Base.length(psi::SparseState{S, BR}) where {S, BR} = length(psi.components)

import Base.iterate
Base.iterate(iter ::SparseState{S, BR}) where {S, BR} = Base.iterate(iter.components)
Base.iterate(iter ::SparseState{S, BR}, i) where {S, BR} =  Base.iterate(iter.components, i)

function choptol!(arg ::SparseState{S1, BR}, tol::Real) where {S1, BR}
  to_delete = [k for (k, v) in arg.components if isapprox(v, 0; atol=tol)]
  for k in to_delete
    delete!(arg.components, k)
  end
  arg
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
  return SparseState{S2, BR}(components)
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
