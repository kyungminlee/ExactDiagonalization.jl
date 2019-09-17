export NullOperator

struct NullOperator <:AbstractOperator end

import Base.-, Base.+, Base.*

(-)(op ::NullOperator) = op

(*)(lhs ::NullOperator, rhs ::NullOperator) where {S<:Number, BR} = lhs
(*)(lhs ::AbstractOperator, rhs ::NullOperator) = rhs
(*)(lhs ::NullOperator, rhs ::AbstractOperator) = lhs
(*)(lhs ::Number, rhs ::NullOperator)::NullOperator = rhs
(*)(lhs ::NullOperator, rhs ::Number)::NullOperator = lhs

(+)(lhs ::NullOperator, rhs ::NullOperator) = lhs
(+)(lhs ::AbstractOperator, rhs ::NullOperator) = lhs
(+)(lhs ::NullOperator, rhs ::AbstractOperator) = rhs

import Base.real, Base.imag
real(arg::NullOperator) = arg
imag(arg::NullOperator) = arg

import Base.isless
# null operator is smaller than any other operators
isless(lhs ::NullOperator, rhs ::NullOperator) = false
isless(lhs ::NullOperator, rhs ::AbstractOperator) = true
isless(lhs ::AbstractOperator, rhs ::NullOperator) = false

export prettyprintln

function prettyprintln(op::NullOperator)
  println("NullOperator")
end