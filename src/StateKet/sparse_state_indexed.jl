export SparseStateIndexed
export clean!

using LinearAlgebra

"""
    struct SparseStateIndexed{S<:Number, BR}

Represents a row vector. Not Free. Access by Integer (1-based index)
"""
mutable struct SparseStateIndexed{S<:Number, QN, BR}
  hilbert_space_realization ::HilbertSpaceRealization{QN, BR}
  components ::Dict{Int, S}

  function SparseStateIndexed{S, QN, BR}(
      hsr ::HilbertSpaceRealization{QN, BR},
      components ::Dict{Int, S}) where {S, QN, BR}
    return new{S, QN, BR}(hsr, components)
  end

  function SparseStateIndexed(
      hsr ::HilbertSpaceRealization{QN, BR},
      components ::Dict{Int, S}) where {S, QN, BR}
    return new{S, QN, BR}(hsr, components)
  end

  function SparseStateIndexed{S, QN, BR}(
      hsr ::HilbertSpaceRealization{QN, BR}) where {S, QN, BR}
    return new{S, QN, BR}(hsr, Dict{Int, S}())
  end

  function SparseStateIndexed{S, QN, BR}(
      hsr ::HilbertSpaceRealization{QN, BR},
      idx ::Integer) where {S, QN, BR}
    return new{S, QN, BR}(hsr, Dict{Int, S}(idx => one(S)))
  end

  function SparseStateIndexed{S, QN, BR}(
      hsr ::HilbertSpaceRealization{QN, BR},
      components::Pair{<:Integer, <:Number}...) where {S, QN, BR}
    return new{S, QN, BR}(hsr, Dict{Int, S}(components))
  end

  function SparseStateIndexed{S, QN, BR}(
      hsr ::HilbertSpaceRealization{QN, BR},
      components ::AbstractDict{<:Integer, S2}) where {S, QN, BR, S2 <:Number}
    # TODO bound checking
    return new{S, QN, BR}(hsr, Dict{Int, S}(components))
  end
end

import Base.getindex, Base.setindex!

@inline function Base.getindex(state ::SparseStateIndexed{S, QN, BR}, idx ::Integer) where {S, QN, BR}
  hsr = state.hilbert_space_realization
  @boundscheck (idx <= 0 || idx > dimension(hsr)) || throw(BoundsError("attempt to access SparseStateIndexed of dimension $(dimension(hsr)) at index [$idx]"))
  return get(state.components, idx, zero(S))
end

function Base.setindex!(state ::SparseStateIndexed{S, QN, BR}, value ::S, idx ::Integer) where {S, QN, BR}
  hsr = state.hilbert_space_realization
  @boundscheck (idx <= 0 || idx > dimension(hsr)) || throw(BoundsError("attempt to access SparseStateIndexed of dimension $(dimension(hsr)) at index [$idx]"))
  state.components[idx] = value
  return state
end

import Base.eltype
Base.eltype(state ::SparseStateIndexed{S, QN, BR}) where {S, QN, BR} = S


import Base.isempty
isempty(psi::SparseStateIndexed{S, QN, BR}) where {S, QN, BR} = isempty(psi.components)


import Base.==
function (==)(lhs ::SparseStateIndexed{S1, QN, BR}, rhs::SparseStateIndexed{S2, QN, BR}) where {S1, S2, QN, BR}
  return ((lhs.hilbert_space_realization === rhs.hilbert_space_realization)
          && (lhs.components == rhs.components))
end

import Base.isapprox
function isapprox(lhs ::SparseStateIndexed{S1, QN, BR}, rhs::SparseStateIndexed{S2, QN, BR}; atol=sqrt(eps(Float64)), rtol=sqrt(eps(Float64))) where {S1, S2, QN, BR}
  if lhs.hilbert_space_realization !== rhs.hilbert_space_realization
    return false
  end
  all_keys = union(keys(lhs.components), keys(rhs.components))
  for k in all_keys
    lv = get(lhs.components, k, zero(S1))
    rv = get(rhs.components, k, zero(S2))
    isapprox(lv, rv; atol=atol, rtol=rtol) || return false
  end
  return true
end

import Base.copy
function copy(arg ::SparseStateIndexed{S, QN, BR}) where {S, QN, BR}
  return SparseStateIndexed{S, BR}(arg.hilbert_space_realization, copy(arg.components))
end

import Base.real, Base.imag, Base.conj
real(arg ::SparseStateIndexed{R, QN, BR}) where {R<:Real, QN, BR} = copy(arg)
imag(arg ::SparseStateIndexed{R, QN, BR}) where {R<:Real, QN, BR} = SparseStateIndexed{R, QN, BR}(arg.hilbert_space_realization)
conj(arg ::SparseStateIndexed{R, QN, BR}) where {R<:Real, QN, BR} = copy(arg)

function real(arg ::SparseStateIndexed{Complex{R}, QN, BR}) where {R<:Real, QN, BR}
  hsr = arg.hilbert_space_realization
  components = Dict{Int, R}((k, real(v)) for (k, v) in arg.components)
  return SparseStateIndexed{R, QN, BR}(hsr, components)
end

function imag(arg ::SparseStateIndexed{Complex{R}, QN, BR}) where {R<:Real, QN, BR}
  hsr = arg.hilbert_space_realization
  components = Dict{Int, R}((k, imag(v)) for (k, v) in arg.components)
  return SparseStateIndexed{R, QN, BR}(hsr, components)
end

function conj(arg ::SparseStateIndexed{Complex{R}, QN, BR}) where {R<:Real, QN, BR}
  C = Complex{R}
  hsr = arg.hilbert_space_realization
  components = Dict{Int, C}((k, conj(v)) for (k, v) in arg.components)
  return SparseStateIndexed{C, QN, BR}(hsr, components)
end

import Base.-, Base.+, Base.*, Base./, Base.\

(+)(arg ::SparseStateIndexed{S, QN, BR}) where {S, QN, BR} = copy(arg)

function (-)(arg ::SparseStateIndexed{S, QN, BR}) where {S, QN, BR}
  return SparseStateIndexed{S, QN, BR}(arg.hilbert_space_realization,
                                       Dict{Int, S}((k, -v) for (k, v) in arg.components))
end

function (+)(lhs ::SparseStateIndexed{S1, QN, BR}, rhs ::SparseStateIndexed{S2, QN, BR}) where {S1, S2, QN, BR}
  if lhs.hilbert_space_realization !== rhs.hilbert_space_realization
    throw(ArgumentError("Hilbert spaces of lhs and rhs of + should match"))
  end
  S3 = promote_type(S1, S2)
  hsr = arg.hilbert_space_realization
  components = Dict{Int, S3}((k, S3(v)) for (k, v) in lhs.components)
  zero_value = zero(S3)
  for (i, v) in rhs.components
    components[i] = get(components, i, zero_value) + v
  end
  return SparseStateIndexed{S, QN, BR}(hsr, components)
end

function (-)(lhs ::SparseStateIndexed{S1, QN, BR}, rhs ::SparseStateIndexed{S2, QN, BR}) where {S1, S2, QN, BR}
  if lhs.hilbert_space_realization !== rhs.hilbert_space_realization
    throw(ArgumentError("Hilbert spaces of lhs and rhs of + should match"))
  end
  S3 = promote_type(S1, S2)
  hsr = arg.hilbert_space_realization
  components = Dict{Int, S3}((k, S3(v)) for (k, v) in lhs.components)
  zero_value = zero(S3)
  for (i, v) in rhs.components
    components[i] = get(components, i, zero_value) - v
  end
  return SparseStateIndexed{S, QN, BR}(hsr, components)
end

function (*)(lhs ::SparseStateIndexed{S1, QN, BR}, rhs ::S2) where {S1, S2<:Number, QN, BR}
  hsr = lhs.hilbert_space_realization
  return SparseStateIndexed(hsr, Dict((k, v * rhs) for (k, v) in lhs.components))
end

function (*)(lhs ::S1, rhs ::SparseStateIndexed{S2, QN, BR}) where {S1<:Number, S2<:Number, QN, BR}
  hsr = rhs.hilbert_space_realization
  return SparseStateIndexed(hsr, Dict((k, lhs * v) for (k, v) in rhs.components))
end

function (/)(lhs ::SparseStateIndexed{S1, QN, BR}, rhs ::S2) where {S1, S2<:Number, QN, BR}
  hsr = lhs.hilbert_space_realization
  return SparseStateIndexed(hsr, Dict((k, v / rhs) for (k, v) in lhs.components))
end

function (\)(lhs ::S1, rhs ::SparseStateIndexed{S2, QN, BR}) where {S1<:Number, S2<:Number, QN, BR}
  hsr = rhs.hilbert_space_realization
  return SparseStateIndexed(hsr, Dict((k, lhs \ v) for (k, v) in rhs.components))
end

import Base.convert
function convert(type ::Type{SparseStateIndexed{S1, QN, BR}},
                 obj::SparseStateIndexed{S2, QN, BR}) where {S1, S2, QN, BR}
  return SparseStateIndexed{S1, QN, BR}(obj.hilbert_space_realization, Dict{Int, S1}(obj.components))
end


import Base.iterate
Base.iterate(arg ::SparseStateIndexed{S, QN, BR}) where {S, QN, BR} = Base.iterate(arg.components)
Base.iterate(arg ::SparseStateIndexed{S, QN, BR}, i::Int) where {S, QN, BR} = Base.iterate(arg.components, i)

function clean!(arg ::SparseStateIndexed{S, QN, BR}; tol=sqrt(eps(Float64))) where {S, QN, BR}
  to_delete = [k for (k, v) in arg.components if isapprox(v, 0; atol=tol)]
  for k in to_delete
    delete!(arg.components, k)
  end
end

import LinearAlgebra.norm
function norm(arg ::SparseStateIndexed{S, QN, BR}) where {S, QN, BR}
  if isempty(arg.components)
    return zero(real(S))
  else
    return norm(values(arg.components))
  end
end

import LinearAlgebra.normalize
function normalize(arg ::SparseStateIndexed{S, BR}) where {S, QN, BR}
  norm_val = norm(arg)
  components = Dict((k, v/norm_val) for (k, v) in arg.components)
  return SparseStateIndexed{S2, BR}(arg.hilbert_space, components)
end

import LinearAlgebra.normalize!
function normalize!(arg ::SparseStateIndexed{S, BR}) where {S, QN, BR}
  norm_val = norm(arg)
  for (k, v) in arg.components
    arg[k] = v / norm_val
  end
  arg
end


# TODO Broadcasting
