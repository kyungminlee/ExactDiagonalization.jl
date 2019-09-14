export GenericOperator

mutable struct GenericOperator
  hilbert_space ::AbstractHilbertSpace
  terms ::Vector{KPO}
  function GenericOperator(hs ::AbstractHilbertSpace)
    return new(hs, KPO[])
  end
end

# import Base.+
# function +(op1 ::KPO, op2 ::KPO)
#   out = GenericOperator
#   push!(out.terms, op1)
#   push!(out.terms, op2)
# end
import Base.push!
function push!(op1 ::GenericOperator, op2 ::KPO)
  @assert op1.hilbert_space == op2.hilbert_space
  push!(op1.terms, op2)
end

