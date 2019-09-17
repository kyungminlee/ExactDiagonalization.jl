export AbstractOperator

abstract type AbstractOperator end

import Base.-
import Base.+
(-)(lhs ::AbstractOperator, rhs::AbstractOperator) = (lhs) + (-rhs)
(+)(op ::AbstractOperator) = op
