mutable struct GenericOperator
  terms ::Vector{KPO}
  function GenericOperator()
    return new(KPO[])
  end
end
