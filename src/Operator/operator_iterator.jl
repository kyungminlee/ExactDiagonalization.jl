export get_row_iterator
export get_column_iterator
export get_element

@inline function get_row_iterator(nullop ::NullOperator, br ::BR) where {BR<:Unsigned}
  return (zero(BR) => false for i in 1:0)
end

@inline function get_column_iterator(nullop ::NullOperator, bc ::BR) where {BR<:Unsigned}
  return (zero(BR) => false for i in 1:0)
end

function get_column_iterator(pureop ::PureOperator{S, BR}, bcol ::BR2) where {S, BR<:Unsigned, BR2<:Unsigned}
  match ::Bool = (bcol & pureop.bitmask) == pureop.bitcol
  element ::Pair{BR, S} = (((bcol & ~pureop.bitmask) | pureop.bitrow) => pureop.amplitude)
  return (element for i in 1:(match ? 1 : 0))
end

function get_row_iterator(pureop ::PureOperator{S, BR}, brow ::BR2) where {S, BR<:Unsigned, BR2<:Unsigned}
  match ::Bool = (brow & pureop.bitmask) == pureop.bitrow
  element ::Pair{BR, S} = (((brow & ~pureop.bitmask) | pureop.bitcol) => pureop.amplitude)
  return (element for i in 1:(match ? 1 : 0))
end

function get_column_iterator(sumop ::SumOperator{S, BR}, bcol ::BR2) where {S, BR<:Unsigned, BR2<:Unsigned}
  match(pureop::PureOperator{S, BR}) ::Bool = (bcol & pureop.bitmask) == pureop.bitcol
  element(pureop::PureOperator{S, BR}) ::Pair{BR, S} = (((bcol & ~pureop.bitmask) | pureop.bitrow) => pureop.amplitude)
  return (element(t) for t in sumop.terms if match(t))
end

function get_row_iterator(sumop::SumOperator{S, BR}, brow ::BR2) where {S, BR<:Unsigned, BR2<:Unsigned}
  match(pureop::PureOperator{S, BR}) ::Bool = (brow & pureop.bitmask) == pureop.bitrow
  element(pureop::PureOperator{S, BR}) ::Pair{BR, S} = (((brow & ~pureop.bitmask) | pureop.bitcol) => pureop.amplitude)
  return (element(t) for t in sumop.terms if match(t))
end

@inline function get_element(nullop::NullOperator, br ::BR2, bc::BR3) ::Bool where {BR2 <:Unsigned, BR3<:Unsigned}
  return false
end

@inline function get_element(pureop::PureOperator{S, BR}, br ::BR2, bc::BR3) ::S where {S, BR, BR2 <:Unsigned, BR3<:Unsigned}
  if ((br & pureop.bitmask) == pureop.bitrow) && ((bc & pureop.bitmask) == pureop.bitcol)
    return pureop.amplitude
  else
    return zero(S)
  end
end

function get_element(sumop::SumOperator{S, BR}, br ::BR2, bc::BR3) where {S, BR, BR2 <:Unsigned, BR3<:Unsigned}
  return sum(get_element(op, br, bc)::S for op in sumop.terms)
end
