export PureOperator
export pure_operator
export bintype

using LinearAlgebra

struct PureOperator{Scalar<:Number, BR<:Unsigned} <:AbstractOperator
  hilbert_space ::HilbertSpace
  bitmask ::BR
  bitrow ::BR  # Row
  bitcol ::BR  # Column
  amplitude ::Scalar

  function PureOperator{S, BR}(hilbert_space ::HilbertSpace,
                               bitmask, bitrow, bitcol, amplitude::S) where {S, BR}
    if (~bitmask) & bitrow != zero(BR)
      throw(ArgumentError("every bit of bitrow not in bitmask should be set to zero"))
    elseif (~bitmask) & bitcol != zero(BR)
      throw(ArgumentError("every bit of bitcol not in bitmask should be set to zero"))
    end
    return new{S, BR}(hilbert_space, bitmask, bitrow, bitcol, amplitude)
  end

  function PureOperator{S, BR}(hilbert_space ::HilbertSpace,
                               s::UniformScaling{S}) where {S, BR}
    return new{S, BR}(hilbert_space, zero(BR), zero(BR), zero(BR), s.λ)
  end
end


import Base.eltype
eltype(lhs ::PureOperator{S, BR}) where {S, BR} = S
eltype(lhs ::Type{PureOperator{S, BR}}) where {S, BR} = S

bintype(lhs ::PureOperator{S, BR}) where {S, BR} = BR
bintype(lhs ::Type{PureOperator{S, BR}}) where {S, BR} = BR


# === 1/6 (In)equality ===

import Base.==
function (==)(lhs ::PureOperator{S1, BR}, rhs::PureOperator{S2, BR}) where {S1, S2, BR}
  return ((lhs.hilbert_space == rhs.hilbert_space) &&
          (lhs.bitmask == rhs.bitmask) &&
          (lhs.bitrow == rhs.bitrow) &&
          (lhs.bitcol == rhs.bitcol) &&
          (lhs.amplitude == rhs.amplitude))
end

import Base.<
function (<)(lhs ::PureOperator{S1, BR}, rhs ::PureOperator{S2, BR}) where {S1, S2, BR}
  if lhs.bitmask < rhs.bitmask
    return true
  elseif lhs.bitmask > rhs.bitmask
    return false
  end

  if lhs.bitrow < rhs.bitrow
    return true
  elseif lhs.bitrow > rhs.bitrow
    return false
  end

  if lhs.bitcol < rhs.bitcol
    return true
  elseif lhs.bitcol > rhs.bitcol
    return false
  end

  return abs(lhs.amplitude) < abs(rhs.amplitude)
end

# === 2/6 Unary functions ===

import Base.real, Base.imag, Base.conj, Base.transpose, Base.adjoint

real(arg ::PureOperator{R, BR}) where {R<:Real, BR} = arg
imag(arg ::PureOperator{R, BR}) where {R<:Real, BR} = PureOperator{R, BR}(arg.hilbert_space, arg.bitmask, arg.bitrow, arg.bitcol, zero(R))
conj(arg ::PureOperator{R, BR}) where {R<:Real, BR} = arg

real(arg ::PureOperator{Complex{R}, BR}) where {R<:Real, BR} = PureOperator{R, BR}(arg.hilbert_space, arg.bitmask, arg.bitrow, arg.bitcol, real(arg.amplitude))
imag(arg ::PureOperator{Complex{R}, BR}) where {R<:Real, BR} = PureOperator{R, BR}(arg.hilbert_space, arg.bitmask, arg.bitrow, arg.bitcol, imag(arg.amplitude))
conj(arg ::PureOperator{Complex{R}, BR}) where {R<:Real, BR} = PureOperator{Complex{R}, BR}(arg.hilbert_space, arg.bitmask, arg.bitrow, arg.bitcol, conj(arg.amplitude))

transpose(arg ::PureOperator{S, BR}) where {S, BR} = PureOperator{S, BR}(arg.hilbert_space, arg.bitmask, arg.bitcol, arg.bitrow, arg.amplitude)
adjoint(arg ::PureOperator{S, BR}) where {S, BR} = PureOperator{S, BR}(arg.hilbert_space, arg.bitmask, arg.bitcol, arg.bitrow, conj(arg.amplitude))


# === 3/6 Scalar Operators ===

import Base.+, Base.-, Base.*, Base./, Base.\

(-)(op ::PureOperator{S, BR}) where {S, BR} = PureOperator{S, BR}(op.hilbert_space, op.bitmask, op.bitrow, op.bitcol, -op.amplitude)


function (*)(lhs ::S1, rhs ::PureOperator{S2, BR}) where {S1<:Number, S2<:Number, BR}
  S = promote_type(S1, S2)
  return PureOperator{S, BR}(rhs.hilbert_space, rhs.bitmask, rhs.bitrow, rhs.bitcol, lhs * rhs.amplitude)
end

function (*)(lhs ::PureOperator{S1, BR}, rhs ::S2) where {S1<:Number, S2<:Number, BR}
  S = promote_type(S1, S2)
  return PureOperator{S, BR}(lhs.hilbert_space, lhs.bitmask, lhs.bitrow, lhs.bitcol, lhs.amplitude * rhs)
end

function (\)(lhs ::S1, rhs ::PureOperator{S2, BR}) where {S1<:Number, S2<:Number, BR}
  S = promote_type(S1, S2)
  return PureOperator{S, BR}(rhs.hilbert_space, rhs.bitmask, rhs.bitrow, rhs.bitcol, lhs \ rhs.amplitude)
end

function (/)(lhs ::PureOperator{S1, BR}, rhs ::S2) where {S1<:Number, S2<:Number, BR}
  S = promote_type(S1, S2)
  return PureOperator{S, BR}(lhs.hilbert_space, lhs.bitmask, lhs.bitrow, lhs.bitcol, lhs.amplitude / rhs)
end


# === 4/6 Operator Products ===

function (*)(lhs ::PureOperator{S1, BR}, rhs ::PureOperator{S2, BR}) where {S1<:Number, S2<:Number, BR}
  if lhs.hilbert_space !== rhs.hilbert_space
    throw(ArgumentError("Hilbert spaces don't match"))
  end
  S3 = promote_type(S1, S2)

  onlylhs_bitmask   =   lhs.bitmask  & (~rhs.bitmask)
  onlyrhs_bitmask   = (~lhs.bitmask) &   rhs.bitmask
  intersect_bitmask =   lhs.bitmask  &   rhs.bitmask
  union_bitmask     =   lhs.bitmask  |   rhs.bitmask

  if (lhs.bitcol & intersect_bitmask) != (rhs.bitrow & intersect_bitmask)
    return NullOperator()
  else
    new_bitmask = union_bitmask
    new_bitrow = lhs.bitrow | (rhs.bitrow & onlyrhs_bitmask)
    new_bitcol = rhs.bitcol | (lhs.bitcol & onlylhs_bitmask)
    new_amplitude = lhs.amplitude * rhs.amplitude
    return PureOperator{S3, BR}(lhs.hilbert_space,
                                new_bitmask,
                                new_bitrow,
                                new_bitcol,
                                new_amplitude)
  end
end


# === 6/6 Conversion ===

import Base.promote_rule
function promote_rule(lhs::Type{PureOperator{S1, BR}}, rhs::Type{PureOperator{S2, BR}}) where {S1, S2, BR}
  S3 = promote_type(S1, S2)
  return PureOperator{S3, BR}
end

import Base.convert
function convert(type ::Type{PureOperator{S1, BR}}, obj::PureOperator{S2, BR}) where {S1, S2, BR}
  return PureOperator{S1, BR}(obj.hilbert_space, obj.bitmask, obj.bitrow, obj.bitcol, convert(S1, obj.amplitude))
end

function pure_operator(hilbert_space ::HilbertSpace,
                       isite ::Integer,
                       istate_row ::Integer,
                       istate_col ::Integer,
                       amplitude::Number=1;
                       dtype ::DataType=UInt)
  @boundscheck let
    site = hilbert_space.sites[isite]
    state_row = site.states[istate_row]
    state_col = site.states[istate_col]
  end
  bm = get_bitmask(hilbert_space, isite; dtype=dtype)
  br = dtype(istate_row - 1) << hilbert_space.bitoffsets[isite]
  bc = dtype(istate_col - 1) << hilbert_space.bitoffsets[isite]
  return PureOperator{typeof(amplitude), dtype}(hilbert_space, bm, br, bc, amplitude)
end

import Base.size
function size(arg::PureOperator{S, BR}) ::Tuple{Int, Int} where {S, BR}
  dim = dimension(arg.hilbert_space)
  return (dim, dim)
end
