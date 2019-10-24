export get_row_iterator
export get_column_iterator

@inline function get_row_iterator(nullop ::NullOperator, br ::BR) where {BR<:Unsigned}
  return Pair{BR, Bool}[]
end

@inline function get_column_iterator(nullop ::NullOperator, bc ::BR) where {BR<:Unsigned}
  return Pair{BR, Bool}[]
end

@inline function get_iterator(nullop ::NullOperator)
  return Pair{BR, Bool}[]
end


@inline function get_column_iterator(pureop ::PureOperator{S, BR}, bcol ::BR) where {S, BR<:Unsigned}
  match(b ::BR) ::Bool = (b & pureop.bitmask) == pureop.bitcol
  element(b ::BR) ::Pair{BR, S} = (((b & ~pureop.bitmask) | pureop.bitrow) => pureop.amplitude)
  return (element(b)  for b in [bcol] if match(b))
end

@inline function get_row_iterator(pureop ::PureOperator{S, BR}, brow ::BR) where {S, BR<:Unsigned}
  match(b ::BR) ::Bool = (b & pureop.bitmask) == pureop.bitrow
  element(b ::BR) ::Pair{BR, S} = (((b & ~pureop.bitmask) | pureop.bitcol) => pureop.amplitude)
  return (element(b) for b in [brow] if match(b))
end

@inline function get_iterator(pureop::PureOperator{S, BR}) where {S, BR<:Unsigned}
  return Base.Iterators.flatten(
    (brow, bcol) => ampl for (brow, ampl) in
    (get_column_iterator(pureop, bcol) for bcol in pureop.hilbert_space.basis_list) )
end


@inline function get_column_iterator(sumop ::SumOperator{S, BR}, bcol ::BR) where {S, BR<:Unsigned}
  return Base.Iterators.flatten(get_column_iterator(t, bcol) for t in sumop.terms)
end

@inline function get_row_iterator(sumop::SumOperator{S, BR}, brow ::BR) where {S, BR<:Unsigned}
  return Base.Iterators.flatten(get_row_iterator(brow, t) for t in sumop.terms)
end

@inline function get_iterator(sumop::SumOperator{S, BR}) where {S, BR<:Unsigned}
  return Base.Iterators.flatten(
    (brow, bcol) => ampl for (brow, ampl) in
    (get_column_iterator(sumop, bcol) for bcol in pureop.hilbert_space.basis_list) )
end
