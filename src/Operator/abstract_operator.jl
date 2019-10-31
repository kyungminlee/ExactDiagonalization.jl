export AbstractOperator
export get_row_iterator, get_column_iterator, get_iterator

abstract type AbstractOperator{S<:Number} end

export scalartype
export bintype
@inline scalartype(lhs::AbstractOperator{S}) where S = S
@inline bintype(lhs::AbstractOperator{S}) where S = bintype(typeof(lhs)) ::DataType

@inline scalartype(lhs::Type{<:AbstractOperator{S}}) where S = S

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
@inline (-)(lhs ::AbstractOperator{S1}, rhs::AbstractOperator{S2}) where {S1, S2} = (lhs) + (-rhs)
@inline (+)(op ::AbstractOperator{S}) where S = op

import LinearAlgebra.issymmetric
function issymmetric(arg::AbstractOperator{S}) where S
  return simplify(arg - transpose(arg)) == NullOperator()
end

import LinearAlgebra.ishermitian
function ishermitian(arg::AbstractOperator{S}) where S
  return simplify(arg - adjoint(arg)) == NullOperator()
end
