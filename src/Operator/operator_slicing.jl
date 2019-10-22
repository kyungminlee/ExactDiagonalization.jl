export get_slice

@inline function get_slice(br ::BR, nullop ::NullOperator) ::Vector{Tuple{BR, Bool}} where {BR}
  return Tuple{BR, Bool}[]
end

@inline function get_slice(nullop ::NullOperator, br ::BR) ::Vector{Tuple{BR, Bool}} where {BR}
  return Tuple{BR, Bool}[]
end

@inline function get_slice(pureop ::PureOperator{S, BR}, bcol ::BR) ::Vector{Tuple{BR, S}} where {S, BR}
  if (bcol & pureop.bitmask) == pureop.bitcol
    brow = (bcol & ~pureop.bitmask) | pureop.bitrow
    return Tuple{BR, S}[(brow, pureop.amplitude)]
  else
    return Tuple{BR, S}[]
  end
end

@inline function get_slice(brow ::BR, pureop ::PureOperator{S, BR}) ::Vector{Tuple{BR, S}} where {S, BR}
  if (brow & pureop.bitmask) == pureop.bitrow
    bcol = (brow & ~pureop.bitmask) | pureop.bitcol
    return Tuple{BR, S}[(bcol, pureop.amplitude)]
  else
    return Tuple{BR, S}[]
  end
end

@inline function get_slice(sumop ::SumOperator{S, BR}, bcol ::BR) ::Vector{Tuple{BR, S}} where {S, BR}
  return vcat((get_slice(t, bcol) for t in sumop.terms)...)
end

@inline function get_slice(brow ::BR, pureop ::SumOperator{S, BR}) ::Vector{Tuple{BR, S}} where {S, BR}
  return vcat((get_slice(brow, t) for t in sumop.terms)...)
end
