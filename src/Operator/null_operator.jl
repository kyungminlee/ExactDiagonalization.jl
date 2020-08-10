export NullOperator


"""
    NullOperator

A null operator, i.e. 0.
"""
struct NullOperator<:AbstractOperator{Bool} end


bintype(lhs::Type{NullOperator}) = UInt8


Base.:(-)(op::NullOperator) = op

Base.:(*)(lhs::NullOperator, rhs::NullOperator) = lhs

Base.:(*)(lhs::AbstractOperator, rhs::NullOperator) = rhs
Base.:(*)(lhs::NullOperator, rhs::AbstractOperator) = lhs

Base.:(*)(lhs::Number, rhs::NullOperator)::NullOperator = rhs
Base.:(*)(lhs::NullOperator, rhs::Number)::NullOperator = lhs
Base.:(\)(lhs::Number, rhs::NullOperator)::NullOperator = rhs
Base.:(/)(lhs::NullOperator, rhs::Number)::NullOperator = lhs
Base.:(//)(lhs::NullOperator, rhs::Number)::NullOperator = lhs

Base.:(+)(lhs::NullOperator, rhs::NullOperator) = lhs
Base.:(+)(lhs::AbstractOperator, rhs::NullOperator) = lhs
Base.:(+)(lhs::NullOperator, rhs::AbstractOperator) = rhs

Base.:(==)(lhs::NullOperator, rhs::NullOperator) = true


Base.real(arg::NullOperator) = arg
Base.imag(arg::NullOperator) = arg
Base.conj(arg::NullOperator) = arg
Base.transpose(arg::NullOperator) = arg
Base.adjoint(arg::NullOperator) = arg


# null operator is less than any other operators
Base.isless(lhs::NullOperator, rhs::NullOperator) = false
Base.isless(lhs::NullOperator, rhs::AbstractOperator) = true
Base.isless(lhs::AbstractOperator, rhs::NullOperator) = false
