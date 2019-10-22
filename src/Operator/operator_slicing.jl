export get_slice

@inline function get_slice(br ::BR, nullop ::NullOperator) where {BR}
  return Pair{BR, Bool}[]
end

@inline function get_slice(nullop ::NullOperator, br ::BR) where {BR}
  return Pair{BR, Bool}[]
end

@inline function get_slice(pureop ::PureOperator{S, BR}, bcol ::BR) where {S, BR}
  match(b ::BR) ::Bool = (b & pureop.bitmask) == pureop.bitcol
  element(b ::BR) ::Pair{BR, S} = (((b & ~pureop.bitmask) | pureop.bitrow) => pureop.amplitude)
  return (element(b) for b in [bcol] if match(b))
  # if (bcol & pureop.bitmask) == pureop.bitcol
  #   brow = (bcol & ~pureop.bitmask) | pureop.bitrow
  #   return Tuple{BR, S}[(brow, pureop.amplitude)]
  # else
  #   return Tuple{BR, S}[]
  # end
end

@inline function get_slice(brow ::BR, pureop ::PureOperator{S, BR}) where {S, BR}
  match(b ::BR) ::Bool = (b & pureop.bitmask) == pureop.bitrow
  element(b ::BR) ::Pair{BR, S} = (((b & ~pureop.bitmask) | pureop.bitcol) => pureop.amplitude)
  return (element(b) for b in [brow] if match(b))
  #
  # if (brow & pureop.bitmask) == pureop.bitrow
  #   bcol = (brow & ~pureop.bitmask) | pureop.bitcol
  #   return Tuple{BR, S}[(bcol, pureop.amplitude)]
  # else
  #   return Tuple{BR, S}[]
  # end
end

@inline function get_slice(sumop ::SumOperator{S, BR}, bcol ::BR) where {S, BR}
  return Base.Iterators.flatten(get_slice(t, bcol) for t in sumop.terms)
end

@inline function get_slice(brow ::BR, pureop ::SumOperator{S, BR}) where {S, BR}
  return Base.Iterators.flatten(get_slice(brow, t) for t in sumop.terms)
end
