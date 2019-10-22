export SparseState
export clean!

using LinearAlgebra

"""
    struct SparseState{Scalar<:Number, BR}

Represents a row vector?
"""
mutable struct SparseState{Scalar<:Number, BR}
  hilbert_space ::HilbertSpace
  components ::DefaultDict{BR, Scalar, Scalar}
  function SparseState{Scalar, BR}(hs ::HilbertSpace) where {Scalar, BR}
    return new{Scalar, BR}(hs, DefaultDict{BR, Scalar, Scalar}(zero(Scalar)))
  end

  function SparseState{Scalar, BR}(hs ::HilbertSpace, binrep ::BR) where {Scalar, BR}
    components = DefaultDict{BR, Scalar, Scalar}(zero(Scalar))
    components[binrep] = one(Scalar)
    return new{Scalar, BR}(hs, components)
  end

  function SparseState{Scalar, BR}(hs ::HilbertSpace, components::Pair{BR, <:Number}...) where {Scalar, BR}
    compos = DefaultDict{BR, Scalar, Scalar}(zero(Scalar))
    #compos[component.first] = component.second
    for (cf, cs) in components
      compos[cf] = cs
    end
    return new{Scalar, BR}(hs, compos)
  end

  function SparseState{Scalar, BR}(hs ::HilbertSpace, components ::AbstractDict{BR, S2}) where {Scalar, BR, S2}
    # TODO bound checking
    return new{Scalar, BR}(hs, components)
  end
end

import Base.getindex, Base.setindex!

function Base.getindex(state ::SparseState{Scalar, BR}, basis ::BR2) where {Scalar, BR, BR2}
  # TODO: check hilbert space
  return state.components[basis]
end

function Base.setindex!(state ::SparseState{Scalar, BR}, value ::S, basis ::BR2) where {Scalar, BR, S<:Number, BR2}
  # TODO: check hilbert space
  Base.setindex!(state.components, value, basis)
  state
end

import Base.eltype
Base.eltype(state ::SparseState{Scalar, BR}) where {Scalar, BR} = Scalar


import Base.isempty
isempty(psi::SparseState{S, BR}) where {S, BR} = isempty(psi.components)


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
    lv = haskey(lhs.components, k) ? lhs.components[k] : zero(S1)
    rv = haskey(rhs.components, k) ? rhs.components[k] : zero(S2)
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
real(arg ::SparseState{R, BR}) where {R<:Real, BR} = arg
imag(arg ::SparseState{R, BR}) where {R<:Real, BR} = SparseState{R, BR}(arg.hilbert_space)
conj(arg ::SparseState{R, BR}) where {R<:Real, BR} = arg

function real(arg ::SparseState{Complex{R}, BR}) where {R<:Real, BR}
  return SparseState{R, BR}(arg.hilbert_space,
                            DefaultDict{BR, R, R}(zero(R),
                                                  [(k, real(v)) for (k, v) in arg.components]
                                                 )
                           )
end

function imag(arg ::SparseState{Complex{R}, BR}) where {R<:Real, BR}
  return SparseState{R, BR}(arg.hilbert_space,
                            DefaultDict{BR, R, R}(zero(R),
                                                  [(k, imag(v)) for (k, v) in arg.components]
                                                 )
                           )
end

function conj(arg ::SparseState{Complex{R}, BR}) where {R<:Real, BR}
  return SparseState{Complex{R}, BR}(arg.hilbert_space,
                                     DefaultDict{BR, Complex{R}, Complex{R}}(zero(Complex{R}),
                                                                             [(k, conj(v)) for (k, v) in arg.components]
                                                                            )
                                    )
end

import Base.-, Base.+, Base.*, Base./

(+)(arg ::SparseState{S, BR}) where {S, BR} = copy(arg)

function (-)(arg ::SparseState{S, BR}) where {S, BR}
  out = SparseState{S, BR}(arg.hilbert_space)
  for (b, v) in arg.components
    out[b] = -v
  end
  return out
end

function (+)(lhs ::SparseState{S1, BR}, rhs ::SparseState{S2, BR}) where {S1, S2, BR}
  S3 = promote_type(S1, S2)
  if lhs.hilbert_space !== rhs.hilbert_space
    throw(ArgumentError("Hilbert spaces of lhs and rhs of + should match"))
  end
  out = SparseState{S3, BR}(lhs.hilbert_space)
  for (b, v) in lhs.components
    out[b] = v
  end
  for (b, v) in rhs.components
    out[b] += v
  end
  return out
end

function (-)(lhs ::SparseState{S1, BR}, rhs ::SparseState{S2, BR}) where {S1, S2, BR}
  S3 = promote_type(S1, S2)
  if lhs.hilbert_space !== rhs.hilbert_space
    throw(ArgumentError("Hilbert spaces of lhs and rhs of + should match"))
  end
  out = SparseState{S3, BR}(lhs.hilbert_space)
  for (b, v) in lhs.components
    out[b] = v
  end
  for (b, v) in rhs.components
    out[b] -= v
  end
  return out
end

function (*)(lhs ::SparseState{S1, BR}, rhs ::S2) where {S1, S2<:Number, BR}
  S3 = promote_type(S1, S2)
  out = SparseState{S3, BR}(lhs.hilbert_space)
  for (b, v) in lhs.components
    out[b] = v * rhs
  end
  return out
end

function (*)(lhs ::S1, rhs ::SparseState{S2, BR}) where {S1<:Number, S2<:Number, BR}
  S3 = promote_type(S1, S2)
  out = SparseState{S3, BR}(rhs.hilbert_space)
  for (b, v) in rhs.components
    out[b] = lhs * v
  end
  return out
end

function (/)(lhs ::SparseState{S1, BR}, rhs ::S2) where {S1, S2<:Number, BR}
  S3 = promote_type(S1, S2)
  out = SparseState{S3, BR}(lhs.hilbert_space)
  for (b, v) in lhs.components
    out[b] = v / rhs
  end
  return out
end

import Base.convert
function convert(type ::Type{SparseState{S1, BR}}, obj::SparseState{S2, BR}) where {S1, S2, BR}
  state = SparseState{S1, BR}(obj.hilbert_space)
  for (k, v) in obj.components
    state[k] = v
  end
  return state
end


import Base.iterate
function Base.iterate(iter ::SparseState{S, BR}, i::Int=iter.components.idxfloor) ::Tuple{Tuple{BR, S}, Int} where {S, BR}
  return Base.iterate(iter.components, i)
end


function clean!(arg ::SparseState{S1, BR}; tol=sqrt(eps(Float64))) where {S1, BR}
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
  components = DefaultDict{BR, S2, S2}(zero(S2), [(k, v/norm_val) for (k, v) in arg.components])
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
