export AbstractOperator
export get_row_iterator, get_column_iterator, get_iterator

abstract type AbstractOperator end

@inline get_row_iterator(op::AbstractOperator, br::BR) where {BR<:Unsigned} = error("get_row_iterator not implemented for $(typeof(op))")
@inline get_column_iterator(op::AbstractOperator, bc::BR) where {BR<:Unsigned} = error("get_column_iterator not implemented for $(typeof(op))")
#@inline get_iterator(op::AbstractOperator) = error("get_iterator not implemented for $(typeof(op))")

#=
 UNARY OPERATORS
 +    NO PO SO
 -    NO PO SO
 real NO PO SO
 imag NO PO SO
=#

#=
 Binary operators

 +/- NO PO SO
 NO  NO PO SO
 PO  PO SO SO
 SO  SO SO SO

  *  NO PO SO
 NO  NO PO SO
 PO  PO PO SO
 SO  SO SO SO
=#

import Base.-, Base.+

(-)(lhs ::AbstractOperator, rhs::AbstractOperator) = (lhs) + (-rhs)
(+)(op ::AbstractOperator) = op
