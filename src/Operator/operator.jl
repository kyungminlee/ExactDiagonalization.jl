export GenericOperator

abstract type AbstractOperator{S} end

mutable struct GenericOperator{S} <: AbstractOperator{S}
  hilbert_space ::AbstractHilbertSpace
  terms ::Vector{AbstractOperator{S}}
  function GenericOperator{S}(hs ::AbstractHilbertSpace) where {S}
    return new{S}(hs, AbstractOperator{S}[])
  end
end

# import Base.+
# function +(op1 ::KPO, op2 ::KPO)
#   out = GenericOperator
#   push!(out.terms, op1)
#   push!(out.terms, op2)
# end
import Base.push!
function push!(op1 ::GenericOperator{S}, op2 ::AbstractOperator{S}) where {S}
  if (op1.hilbert_space != op2.hilbert_space)
    throw(ArgumentError("Hilbert spaces of lhs and rhs must match"))
  end
  push!(op1.terms, op2)
end


import Base.*

function *(lhs ::SparseState{BR, S1}, rhs ::GenericOperator{S2}) where {BR<:Unsigned, S1<:Number, S2<:Number}
  if (lhs.hilbert_space != rhs.hilbert_space)
    throw(ArgumentError("Hilbert spaces of lhs and rhs must match"))
  end
  S3 = promote_type(S1, S2)
  output = SparseState{BR, S3}(lhs.hilbert_space)
  for term in rhs.terms
    partial_output = lhs * term
    for (b, v) in partial_output.components
      output.components[b] += v
    end
  end
  return output
end
