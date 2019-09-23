export SparseState

"""
    struct SparseState{Scalar<:Number, BR}

Represents a row vector?
"""
mutable struct SparseState{Scalar<:Number, BR}
  hilbert_space ::AbstractHilbertSpace
  components ::DefaultDict{BR, Scalar, Scalar}
  function SparseState{Scalar, BR}(hs ::AbstractHilbertSpace) where {Scalar, BR}
    return new{Scalar, BR}(hs, DefaultDict{BR, Scalar, Scalar}(zero(Scalar)))
  end

  function SparseState{Scalar, BR}(hs ::AbstractHilbertSpace, binrep ::BR) where {Scalar, BR}
    components = DefaultDict{BR, Scalar, Scalar}(zero(Scalar))
    components[binrep] = one(Scalar)
    return new{Scalar, BR}(hs, components)
  end

  function SparseState{Scalar, BR}(hs ::AbstractHilbertSpace, component::Pair{BR, S2}, rest...) where {Scalar, BR, S2}
    components = DefaultDict{BR, Scalar, Scalar}(zero(Scalar))
    components[component.first] = component.second
    for (cf, cs) in rest
      components[cf] = cs
    end
    return new{Scalar, BR}(hs, components)
  end

  function SparseState{Scalar, BR}(hs ::AbstractHilbertSpace, components ::AbstractDict{BR, S2}) where {Scalar, BR, S2}
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


import Base.==
function (==)(lhs ::SparseState{S1, BR}, rhs::SparseState{S2, BR}) where {S1, S2, BR}
  return (lhs.hilbert_space == rhs.hilbert_space) && (lhs.components == rhs.components)
end


import Base.real, Base.imag, Base.conj
real(arg ::SparseState{R, BR}) where {R<:Real, BR} = arg
imag(arg ::SparseState{R, BR}) where {R<:Real, BR} = SparseState{R, BR}(arg.hilbert_space)
conj(arg ::SparseState{R, BR}) where {R<:Real, BR} = arg

function real(arg ::SparseState{Complex{R}, BR}) where {R<:Real, BR}
  components = DefaultDict{BR, R, R}(zero(R), [(k, real(v)) for (k, v) in arg.components])
  return SparseState{R, BR}(arg.hilbert_space, components)
end

function imag(arg ::SparseState{Complex{R}, BR}) where {R<:Real, BR}
  components = DefaultDict{BR, R, R}(zero(R), [(k, imag(v)) for (k, v) in arg.components])
  return SparseState{R, BR}(arg.hilbert_space, components)
end

function conj(arg ::SparseState{Complex{R}, BR}) where {R<:Real, BR}
  components = DefaultDict{BR, Complex{R}, Complex{R}}(zero(Complex{R}), [(k, conj(v)) for (k, v) in arg.components])
  return SparseState{Complex{R}, BR}(arg.hilbert_space, components)
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
