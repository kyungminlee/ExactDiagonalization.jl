export AbstractOperator
export get_row_iterator, get_column_iterator, get_iterator

abstract type AbstractOperator end

export scalartype
export bintype
@inline scalartype(lhs::AbstractOperator) = scalartype(typeof(lhs)) ::DataType
@inline bintype(lhs::AbstractOperator) = bintype(typeof(lhs)) ::DataType

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
@inline (-)(lhs ::AbstractOperator, rhs::AbstractOperator) = (lhs) + (-rhs)
@inline (+)(op ::AbstractOperator) = op
