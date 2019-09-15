export SparseState

mutable struct SparseState{BR, Scalar<:Number}
  hilbert_space ::AbstractHilbertSpace
  components ::DefaultDict{BR, Scalar, Scalar}
  function SparseState{BR, Scalar}(hs ::AbstractHilbertSpace) where {BR, Scalar <: Number}
    return new{BR, Scalar}(hs, DefaultDict{BR, Scalar, Scalar}(zero(Scalar)))
  end

  function SparseState{BR, Scalar}(hs ::AbstractHilbertSpace, binrep ::BR) where {BR, Scalar<:Number}
    components = DefaultDict{BR, Scalar, Scalar}(zero(Scalar))
    components[binrep] = one(Scalar)
    return new{BR, Scalar}(hs, components)
  end
end

import Base.getindex, Base.setindex!

function Base.getindex(state ::SparseState{BR, Scalar}, basis ::BR) where {BR, Scalar <:Number}
  return state.components[basis]
end

function Base.setindex!(state ::SparseState{BR, Scalar}, value ::Scalar, basis ::BR) where {BR, Scalar <:Number}
  Base.setindex!(state.components, value, basis)
end

import Base.+

function +(lhs ::SparseState{BR, S1}, rhs ::SparseState{BR, S2}) where {BR, S1 <:Number, S2 <:Number}
  S3 = promote_type(S1, S2)
  if lhs.hilbert_space != rhs.hilbert_space
    throw(ArgumentError("Hilbert spaces of lhs and rhs of + should match"))
  end
  out = SparseState{BR, S3}(lhs.hilbert_space)
  for (b, v) in lhs.components
    out[b] += v
  end
  for (b, v) in rhs.components
    out[b] += v
  end
  return out 
end