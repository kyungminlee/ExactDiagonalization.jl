export NullOperator

struct NullOperator <:AbstractOperator end

import Base.-, Base.+, Base.*, Base.==

(-)(op ::NullOperator) = op

(*)(lhs ::NullOperator, rhs ::NullOperator) = lhs

(*)(lhs ::AbstractOperator, rhs ::NullOperator) = rhs
(*)(lhs ::NullOperator, rhs ::AbstractOperator) = lhs

(*)(lhs ::Number, rhs ::NullOperator)::NullOperator = rhs
(*)(lhs ::NullOperator, rhs ::Number)::NullOperator = lhs

(+)(lhs ::NullOperator, rhs ::NullOperator) = lhs
(+)(lhs ::AbstractOperator, rhs ::NullOperator) = lhs
(+)(lhs ::NullOperator, rhs ::AbstractOperator) = rhs

(==)(lhs ::NullOperator, rhs::NullOperator) = true
# (==)(lhs ::NullOperator, rhs::AbstractOperator) = false
# (==)(lhs ::AbstractOperator, rhs::NullOperator) = false

import Base.real, Base.imag, Base.conj, Base.transpose
real(arg::NullOperator) = arg
imag(arg::NullOperator) = arg
conj(arg::NullOperator) = arg
transpose(arg::NullOperator) = arg


import Base.<
# null operator is smaller than any other operators
(<)(lhs ::NullOperator, rhs ::NullOperator) = false
(<)(lhs ::NullOperator, rhs ::AbstractOperator) = true
(<)(lhs ::AbstractOperator, rhs ::NullOperator) = false
