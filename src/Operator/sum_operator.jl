
export SumOperator
export bintype

struct SumOperator{Scalar<:Number, BR <:Unsigned} <:AbstractOperator{Scalar}
  terms ::Vector{PureOperator{Scalar, BR}}

  function SumOperator{S, BR}(terms) where {S, BR}
    return new{S, BR}(terms)
  end
end

@inline scalartype(lhs ::Type{SumOperator{S, BR}}) where {S, BR} = S ::DataType
@inline bintype(lhs ::Type{SumOperator{S, BR}}) where {S, BR} = BR ::DataType

# === 1/6 Equality ===

import Base.==
@inline function (==)(lhs ::SumOperator{S1, BR}, rhs::SumOperator{S2, BR}) where {S1, S2, BR}
  return (lhs.terms == rhs.terms)
end


# === 2/6 Unary functions ===

import Base.real, Base.imag, Base.conj, Base.transpose, Base.adjoint

@inline real(arg ::SumOperator{S, BR}) where {S<:Real, BR} = arg
@inline imag(arg ::SumOperator{S, BR}) where {S<:Real, BR} = SumOperator{S, BR}([])
@inline conj(arg ::SumOperator{S, BR}) where {S<:Real, BR} = arg

@inline real(arg ::SumOperator{Complex{S}, BR}) where {S<:Real, BR} = SumOperator{S, BR}(real.(arg.terms))
@inline imag(arg ::SumOperator{Complex{S}, BR}) where {S<:Real, BR} = SumOperator{S, BR}(imag.(arg.terms))
@inline conj(arg ::SumOperator{Complex{S}, BR}) where {S<:Real, BR} = SumOperator{Complex{S}, BR}(conj.(arg.terms))

@inline transpose(arg ::SumOperator{S, BR}) where {S, BR} = SumOperator{S, BR}(transpose.(arg.terms))
@inline adjoint(arg ::SumOperator{S, BR}) where {S, BR} = SumOperator{S, BR}(adjoint.(arg.terms))


# === 3/6 Scalar Operators ===

import Base.+, Base.-, Base.*, Base./, Base.\

@inline (-)(arg ::SumOperator{S, BR}) where {S, BR} = SumOperator{S, BR}(-arg.terms)


@inline function (*)(lhs ::S1, rhs ::SumOperator{S2, BR}) where {S1<:Number, S2<:Number, BR}
  S = promote_type(S1, S2)
  SumOperator{S, BR}(lhs .* rhs.terms)
end

@inline function (*)(lhs ::SumOperator{S1, BR}, rhs ::S2) where {S1<:Number, S2<:Number, BR}
  S = promote_type(S1, S2)
  SumOperator{S, BR}(lhs.terms .* rhs)
end

@inline function (\)(lhs ::S1, rhs ::SumOperator{S2, BR}) where {S1<:Number, S2<:Number, BR}
  S = promote_type(S1, S2)
  SumOperator{S, BR}(lhs .\ rhs.terms)
end

@inline function (/)(lhs ::SumOperator{S1, BR}, rhs ::S2) where {S1<:Number, S2<:Number, BR}
  S = promote_type(S1, S2)
  SumOperator{S, BR}(lhs.terms ./ rhs)
end


# === 4/6 Operator Products ===

function (*)(lhs::SumOperator{S1, BR}, rhs::PureOperator{S2, BR}) where {S1, S2, BR}
  S3 = promote_type(S1, S2)
  terms ::Vector{PureOperator{S3, BR}} = filter((x) -> !isa(x, NullOperator), [t * rhs for t in lhs.terms])
  return SumOperator{S3, BR}(terms)
end

function (*)(lhs::PureOperator{S1, BR}, rhs::SumOperator{S2, BR}) where {S1, S2, BR}
  S3 = promote_type(S1, S2)
  terms ::Vector{PureOperator{S3, BR}} = filter((x) -> !isa(x, NullOperator), [lhs * t for t in rhs.terms])
  return SumOperator{S3, BR}(terms)
end

function (*)(lhs::SumOperator{S1, BR}, rhs::SumOperator{S2, BR}) where {S1, S2, BR}
  S3 = promote_type(S1, S2)
  terms ::Vector{PureOperator{S3, BR}} = filter((x) -> !isa(x, NullOperator), [tl * tr for tl in lhs.terms, tr in rhs.terms])
  return SumOperator{S3, BR}(vec(terms))
end


# === 5/6 Operator Sums ===

@inline function (+)(lhs::PureOperator{S1, BR}, rhs::PureOperator{S2, BR}) where {S1, S2, BR}
  S = promote_type(S1, S2)
  return SumOperator{S, BR}(PureOperator{S, BR}[lhs, rhs])
end

@inline function (+)(lhs::SumOperator{S1, BR}, rhs::PureOperator{S2, BR}) where {S1, S2, BR}
  S3 = promote_type(S1, S2)
  return SumOperator{S3, BR}(PureOperator{S3, BR}[lhs.terms..., rhs])
end

@inline function (+)(lhs::PureOperator{S1, BR}, rhs::SumOperator{S2, BR}) where {S1, S2, BR}
  S3 = promote_type(S1, S2)
  return SumOperator{S3, BR}(PureOperator{S3, BR}[lhs, rhs.terms...])
end

@inline function (+)(lhs::SumOperator{S1, BR}, rhs::SumOperator{S2, BR}) where {S1, S2, BR}
  S3 = promote_type(S1, S2)
  return SumOperator{S3, BR}(PureOperator{S3, BR}[lhs.terms..., rhs.terms...])
end


# === 6/6 Conversion ===

import Base.promote_rule
@inline function promote_rule(lhs::Type{SumOperator{S1, BR}}, rhs::Type{SumOperator{S2, BR}}) where {S1, S2, BR}
  S3 = promote_type(S1, S2)
  return SumOperator{S3, BR}
end

import Base.convert
@inline function convert(type ::Type{SumOperator{S1, BR}}, obj::SumOperator{S2, BR}) where {S1, S2, BR}
  return SumOperator{S1, BR}([convert(PureOperator{S1, BR}, t) for t in obj.terms])
end