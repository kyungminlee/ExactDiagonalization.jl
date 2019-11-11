export AbstractOperator
export get_row_iterator, get_column_iterator, get_iterator

abstract type AbstractOperator{S<:Number} end

export scalartype
export bintype
scalartype(lhs::AbstractOperator{S}) where S = S
bintype(lhs::AbstractOperator{S}) where S = bintype(typeof(lhs)) ::DataType

scalartype(lhs::Type{<:AbstractOperator{S}}) where S = S

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
(-)(lhs ::AbstractOperator{S1}, rhs::AbstractOperator{S2}) where {S1, S2} = (lhs) + (-rhs)
(+)(op ::AbstractOperator{S}) where S = op

import LinearAlgebra.issymmetric
function issymmetric(arg::AbstractOperator{S}) where S
  return simplify(arg - transpose(arg)) == NullOperator()
end

import LinearAlgebra.ishermitian
function ishermitian(arg::AbstractOperator{S}) where S
  return simplify(arg - adjoint(arg)) == NullOperator()
end

import Base.^
function ^(lhs ::AbstractOperator{S}, p ::Integer) where S
  p <= 0 && error("Non-positive power for AbstractOperator not supported")

  # smallest nonzero power
  pow = simplify(lhs)
  while (p & 0x1) == 0
    pow = simplify(pow * pow)
    if isa(pow, NullOperator)
      return pow
    end
    p = p >> 1
  end

  out = pow
  p = p >> 1
  while p > 0
    pow = simplify(pow * pow)
    if (p & 0b1) != 0
      out = simplify(out * pow)
    end
    if isa(out, NullOperator)
      return out
    end
    p = p >> 1
  end
  return out
end
