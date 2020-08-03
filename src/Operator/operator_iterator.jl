export get_row_iterator
export get_column_iterator
export get_element


"""
    get_row_iterator(op, br)

Returns an iterator over the elements of the row corresponding to bit representation `br`.
"""
function get_row_iterator(nullop::NullOperator, br::BR) where {BR<:Unsigned}
    return ((zero(UInt8) => false) for i in 1:0)
end


"""
    get_column_iterator(op, bcol)

Returns an iterator over the elements of the column corresponding to bit representation `bc`.
"""
function get_column_iterator(nullop::NullOperator, bcol::BR) where {BR<:Unsigned}
    return ((zero(UInt8) => false) for i in 1:0)
end


function get_row_iterator(pureop::PureOperator{S, BR}, brow::BR2) where {S, BR<:Unsigned, BR2<:Unsigned}
    match::Bool = (brow & pureop.bitmask) == pureop.bitrow
    element::Pair{BR, S} = (((brow & ~pureop.bitmask) | pureop.bitcol) => pureop.amplitude)
    return (element for i in 1:(match ? 1 : 0))
end


function get_column_iterator(pureop::PureOperator{S, BR}, bcol::BR2) where {S, BR<:Unsigned, BR2<:Unsigned}
    match::Bool = (bcol & pureop.bitmask) == pureop.bitcol
    element::Pair{BR, S} = (((bcol & ~pureop.bitmask) | pureop.bitrow) => pureop.amplitude)
    return (element for i in 1:(match ? 1 : 0))
end


function get_row_iterator(sumop::SumOperator{S, BR}, brow::BR2) where {S, BR<:Unsigned, BR2<:Unsigned}
    let brow::BR = BR(brow), terms::Vector{PureOperator{S, BR}} = sumop.terms
        match(pureop::PureOperator{S, BR})::Bool = ((brow & pureop.bitmask) == pureop.bitrow)
        element(pureop::PureOperator{S, BR})::Pair{BR, S} = (((brow & ~pureop.bitmask) | pureop.bitcol) => pureop.amplitude)
        return (element(t) for t in terms if match(t))
    end
end


function get_column_iterator(sumop::SumOperator{S, BR}, bcol::BR2) where {S, BR<:Unsigned, BR2<:Unsigned}
    let bcol::BR = BR(bcol), terms::Vector{PureOperator{S, BR}} = sumop.terms
        match(pureop::PureOperator{S, BR})::Bool = (bcol & pureop.bitmask) == pureop.bitcol
        element(pureop::PureOperator{S, BR})::Pair{BR, S} = (((bcol & ~pureop.bitmask) | pureop.bitrow) => pureop.amplitude)
        return (element(t) for t in terms if match(t))
    end
end


function get_element(nullop::NullOperator, br::BR2, bc::BR3)::Bool where {BR2 <:Unsigned, BR3<:Unsigned}
    return false
end


function get_element(pureop::PureOperator{S, BR}, br::BR2, bc::BR3)::S where {S, BR, BR2 <:Unsigned, BR3<:Unsigned}
    if ((br & pureop.bitmask) == pureop.bitrow) && (((br & ~pureop.bitmask) | pureop.bitcol) == bc)
        return pureop.amplitude
    else
        return zero(S)
    end
end


function get_element(sumop::SumOperator{S, BR}, br::BR2, bc::BR3) where {S, BR, BR2 <:Unsigned, BR3<:Unsigned}
    element(op::PureOperator{S, BR}) = get_element(op, br, bc)
    return mapreduce(element, +, sumop.terms; init=zero(S))
end
