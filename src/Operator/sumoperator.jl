
export SumOperator
export bintype

struct SumOperator{Scalar<:Number, BR <:Unsigned} <:AbstractOperator
  hilbert_space ::HilbertSpace
  terms ::Vector{PureOperator{Scalar, BR}}

  function SumOperator{S, BR}(hs ::HilbertSpace, terms) where {S, BR}
    @boundscheck if any(!isa(t,  NullOperator) && t.hilbert_space != hs for t in terms)
      throw(ArgumentError("Hilbert spaces don't match"))
    end
    return new{S, BR}(hs, terms)
  end
end


import Base.eltype
#eltype(lhs ::SumOperator{S, BR}) where {S, BR} = S
@inline eltype(lhs ::Type{SumOperator{S, BR}}) where {S, BR} = S

@inline bintype(lhs ::SumOperator{S, BR}) where {S, BR} = BR
@inline bintype(lhs ::Type{SumOperator{S, BR}}) where {S, BR} = BR

import Base.size
function size(arg::SumOperator{S, BR}) ::Tuple{Int, Int} where {S, BR}
  dim = dimension(arg.hilbert_space)
  return (dim, dim)
end

# === 1/6 Equality ===

import Base.==
function (==)(lhs ::SumOperator{S1, BR}, rhs::SumOperator{S2, BR}) where {S1, S2, BR}
  return (lhs.hilbert_space == rhs.hilbert_space) && (lhs.terms == rhs.terms)
end


# === 2/6 Unary functions ===

import Base.real, Base.imag, Base.conj, Base.transpose, Base.adjoint

real(arg ::SumOperator{S, BR}) where {S<:Real, BR} = arg
imag(arg ::SumOperator{S, BR}) where {S<:Real, BR} = SumOperator{S, BR}(arg.hilbert_space, [])
conj(arg ::SumOperator{S, BR}) where {S<:Real, BR} = arg

real(arg ::SumOperator{Complex{S}, BR}) where {S<:Real, BR} = SumOperator{S, BR}(arg.hilbert_space, real.(arg.terms))
imag(arg ::SumOperator{Complex{S}, BR}) where {S<:Real, BR} = SumOperator{S, BR}(arg.hilbert_space, imag.(arg.terms))
conj(arg ::SumOperator{Complex{S}, BR}) where {S<:Real, BR} = SumOperator{Complex{S}, BR}(arg.hilbert_space, conj.(arg.terms))

transpose(arg ::SumOperator{S, BR}) where {S, BR} = SumOperator{S, BR}(arg.hilbert_space, transpose.(arg.terms))
adjoint(arg ::SumOperator{S, BR}) where {S, BR} = SumOperator{S, BR}(arg.hilbert_space, adjoint.(arg.terms))


# === 3/6 Scalar Operators ===

import Base.+, Base.-, Base.*, Base./, Base.\

(-)(arg ::SumOperator{S, BR}) where {S, BR} = SumOperator{S, BR}(arg.hilbert_space, -arg.terms)


function (*)(lhs ::S1, rhs ::SumOperator{S2, BR}) where {S1<:Number, S2<:Number, BR}
  S = promote_type(S1, S2)
  SumOperator{S, BR}(rhs.hilbert_space, lhs .* rhs.terms)
end

function (*)(lhs ::SumOperator{S1, BR}, rhs ::S2) where {S1<:Number, S2<:Number, BR}
  S = promote_type(S1, S2)
  SumOperator{S, BR}(lhs.hilbert_space, lhs.terms .* rhs)
end

function (\)(lhs ::S1, rhs ::SumOperator{S2, BR}) where {S1<:Number, S2<:Number, BR}
  S = promote_type(S1, S2)
  SumOperator{S, BR}(rhs.hilbert_space, lhs .\ rhs.terms)
end

function (/)(lhs ::SumOperator{S1, BR}, rhs ::S2) where {S1<:Number, S2<:Number, BR}
  S = promote_type(S1, S2)
  SumOperator{S, BR}(lhs.hilbert_space, lhs.terms ./ rhs)
end


# === 4/6 Operator Products ===

function (*)(lhs::SumOperator{S1, BR}, rhs::PureOperator{S2, BR}) where {S1, S2, BR}
  @boundscheck if lhs.hilbert_space != rhs.hilbert_space
    throw(ArgumentError("Hilbert spaces don't match"))
  end
  S3 = promote_type(S1, S2)
  terms ::Vector{PureOperator{S3, BR}} = filter((x) -> !isa(x, NullOperator), [t * rhs for t in lhs.terms])
  return SumOperator{S3, BR}(lhs.hilbert_space, terms)
end

function (*)(lhs::PureOperator{S1, BR}, rhs::SumOperator{S2, BR}) where {S1, S2, BR}
  @boundscheck if lhs.hilbert_space != rhs.hilbert_space
    throw(ArgumentError("Hilbert spaces don't match"))
  end
  S3 = promote_type(S1, S2)
  terms ::Vector{PureOperator{S3, BR}} = filter((x) -> !isa(x, NullOperator), [lhs * t for t in rhs.terms])
  return SumOperator{S3, BR}(rhs.hilbert_space, terms)
end

function (*)(lhs::SumOperator{S1, BR}, rhs::SumOperator{S2, BR}) where {S1, S2, BR}
  @boundscheck if lhs.hilbert_space != rhs.hilbert_space
    throw(ArgumentError("Hilbert spaces don't match"))
  end
  S3 = promote_type(S1, S2)
  terms ::Vector{PureOperator{S3, BR}} = filter((x) -> !isa(x, NullOperator), [tl * tr for tl in lhs.terms, tr in rhs.terms])
  return SumOperator{S3, BR}(lhs.hilbert_space, vec(terms))
end


# === 5/6 Operator Sums ===

function (+)(lhs::PureOperator{S1, BR}, rhs::PureOperator{S2, BR}) where {S1, S2, BR}
  S = promote_type(S1, S2)
  return SumOperator{S, BR}(lhs.hilbert_space, PureOperator{S, BR}[lhs, rhs])
end

function (+)(lhs::SumOperator{S1, BR}, rhs::PureOperator{S2, BR}) where {S1, S2, BR}
  @boundscheck if lhs.hilbert_space != rhs.hilbert_space
    throw(ArgumentError("Hilbert spaces don't match"))
  end
  S3 = promote_type(S1, S2)
  return SumOperator{S3, BR}(lhs.hilbert_space, PureOperator{S3, BR}[lhs.terms..., rhs])
end

function (+)(lhs::PureOperator{S1, BR}, rhs::SumOperator{S2, BR}) where {S1, S2, BR}
  @boundscheck if lhs.hilbert_space != rhs.hilbert_space
    throw(ArgumentError("Hilbert spaces don't match"))
  end
  S3 = promote_type(S1, S2)
  return SumOperator{S3, BR}(lhs.hilbert_space, PureOperator{S3, BR}[lhs, rhs.terms...])
end

function (+)(lhs::SumOperator{S1, BR}, rhs::SumOperator{S2, BR}) where {S1, S2, BR}
  @boundscheck if lhs.hilbert_space != rhs.hilbert_space
    throw(ArgumentError("Hilbert spaces don't match"))
  end
  S3 = promote_type(S1, S2)
  return SumOperator{S3, BR}(lhs.hilbert_space, PureOperator{S3, BR}[lhs.terms..., rhs.terms...])
end


# === 6/6 Conversion ===

import Base.promote_rule
function promote_rule(lhs::Type{SumOperator{S1, BR}}, rhs::Type{SumOperator{S2, BR}}) where {S1, S2, BR}
  S3 = promote_type(S1, S2)
  return SumOperator{S3, BR}
end

import Base.convert
function convert(type ::Type{SumOperator{S1, BR}}, obj::SumOperator{S2, BR}) where {S1, S2, BR}
  return SumOperator{S1, BR}(obj.hilbert_space, [convert(PureOperator{S1, BR}, t) for t in obj.terms])
end
