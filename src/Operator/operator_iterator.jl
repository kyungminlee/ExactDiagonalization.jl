export get_row_iterator
export get_column_iterator
export get_element

@inline function get_row_iterator(nullop ::NullOperator, br ::BR) where {BR<:Unsigned}
  return Pair{Bool, Bool}[]
end

@inline function get_column_iterator(nullop ::NullOperator, bc ::BR) where {BR<:Unsigned}
  return Pair{Bool, Bool}[]
end

@inline function get_column_iterator(pureop ::PureOperator{S, BR}, bcol ::BR2) where {S, BR<:Unsigned, BR2<:Unsigned}
  match(b ::BR2) ::Bool = (b & pureop.bitmask) == pureop.bitcol
  element(b ::BR2) ::Pair{BR, S} = (((b & ~pureop.bitmask) | pureop.bitrow) => pureop.amplitude)
  #range = match(bcol) ? (1:1) : (1:0)
  #return (element(bcol) for x in range)
  return match(bcol) ? Pair{BR, S}[element(bcol)] : Pair{BR, S}[]
end

@inline function get_row_iterator(pureop ::PureOperator{S, BR}, brow ::BR2) where {S, BR<:Unsigned, BR2<:Unsigned}
  match(b ::BR2) ::Bool = (b & pureop.bitmask) == pureop.bitrow
  element(b ::BR2) ::Pair{BR, S} = (((b & ~pureop.bitmask) | pureop.bitcol) => pureop.amplitude)
  return match(brow) ? Pair{BR, S}[element(brow)] : Pair{BR, S}[]
end

@inline function get_column_iterator(sumop ::SumOperator{S, BR}, bcol ::BR2) where {S, BR<:Unsigned, BR2<:Unsigned}
  match(pureop::PureOperator{S, BR}, b ::BR2) ::Bool = (b & pureop.bitmask) == pureop.bitcol
  element(pureop::PureOperator{S, BR}, b ::BR2) ::Pair{BR, S} = (((b & ~pureop.bitmask) | pureop.bitrow) => pureop.amplitude)
  return (element(t, bcol) for t in sumop.terms if match(t, bcol))
  #return Base.Iterators.flatten(get_column_iterator(t, bcol) for t in sumop.terms)
end

@inline function get_row_iterator(sumop::SumOperator{S, BR}, brow ::BR2) where {S, BR<:Unsigned, BR2<:Unsigned}
  match(pureop::PureOperator{S, BR}, b ::BR2) ::Bool = (b & pureop.bitmask) == pureop.bitrow
  element(pureop::PureOperator{S, BR}, b ::BR2) ::Pair{BR, S} = (((b & ~pureop.bitmask) | pureop.bitcol) => pureop.amplitude)
  return (element(t, brow) for t in sumop.terms if match(t, brow))
  #return Base.Iterators.flatten(get_row_iterator(t, brow) for t in sumop.terms)
end

@inline function get_element(nullop::NullOperator, br ::BR2, bc::BR3) ::Bool where {BR2 <:Unsigned, BR3<:Unsigned}
  return false
end

@inline function get_element(pureop::PureOperator{S, BR}, br ::BR2, bc::BR3) ::S where {S, BR, BR2 <:Unsigned, BR3<:Unsigned}
  match(b ::BR2) ::Bool = (b & pureop.bitmask) == pureop.bitcol
  element(b ::BR2) ::Pair{BR, S} = (((b & ~pureop.bitmask) | pureop.bitrow) => pureop.amplitude)
  if ((br & pureop.bitmask) == pureop.bitrow) && ((bc & pureop.bitmask) == pureop.bitcol)
    return pureop.amplitude
  else
    return zero(S)
  end
end

@inline function get_element(sumop::SumOperator{S, BR}, br ::BR2, bc::BR3) where {S, BR, BR2 <:Unsigned, BR3<:Unsigned}
  return sum(get_element(op, br, bc)::S for op in sumop.terms)
end
