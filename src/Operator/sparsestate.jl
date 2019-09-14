export SparseState

mutable struct SparseState{BinRep, Scalar <: Number}
  hilbert_space ::AbstractHilbertSpace
  components ::DefaultDict{BinRep, Scalar, Scalar}
  function SparseState{BinRep, Scalar}(hs ::AbstractHilbertSpace) where {BinRep, Scalar <: Number}
    return new{BinRep, Scalar}(hs, DefaultDict{BinRep, Scalar, Scalar}(zero(Scalar)))
  end

  function SparseState{BinRep, Scalar}(hs ::AbstractHilbertSpace, binrep ::BinRep) where {BinRep, Scalar<:Number}
    components = DefaultDict{BinRep, Scalar, Scalar}(zero(Scalar))
    components[binrep] = one(Scalar)
    return new{BinRep, Scalar}(hs, components)
  end
end

import Base.getindex, Base.setindex!

function Base.getindex(state ::SparseState{BinRep, Scalar}, basis ::BinRep) where {BinRep, Scalar <:Number}
  return state.components[basis]
end

function Base.setindex!(state ::SparseState{BinRep, Scalar}, value ::Scalar, basis ::BinRep) where {BinRep, Scalar <:Number}
  Base.setindex!(state.components, value, basis)
end

