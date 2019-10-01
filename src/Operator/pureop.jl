export PureOperator
export pure_operator

struct PureOperator{Scalar<:Number, BR<:Unsigned} <:AbstractOperator
  hilbert_space ::HilbertSpace
  bitmask ::BR
  bitrow ::BR  # Row
  bitcol ::BR  # Column
  amplitude ::Scalar
  # TODO: fermion sign

  function PureOperator{S, BR}(hilbert_space ::HilbertSpace,
                               bitmask, bitrow, bitcol, amplitude::S) where {S, BR}
    if (~bitmask) & bitrow != 0x0
      throw(ArgumentError("every bit of bitrow not in bitmask should be set to zero"))
    elseif (~bitmask) & bitcol != 0x0
      throw(ArgumentError("every bit of bitcol not in bitmask should be set to zero"))
    end
    return new{S, BR}(hilbert_space, bitmask, bitrow, bitcol, amplitude)
  end
end

#const OptionalPureOperator{Scalar, BR} = Union{PureOperator{Scalar, BR}, NullOperator} where {Scalar, BR}

import Base.==

function (==)(lhs ::PureOperator{S1, BR}, rhs::PureOperator{S2, BR}) where {S1, S2, BR}
  return ((lhs.hilbert_space == rhs.hilbert_space) &&
          (lhs.bitmask == rhs.bitmask) &&
          (lhs.bitrow == rhs.bitrow) &&
          (lhs.bitcol == rhs.bitcol) &&
          (lhs.amplitude == rhs.amplitude))
end

(-)(op ::PureOperator{S, BR}) where {S, BR} = PureOperator{S, BR}(op.hilbert_space, op.bitmask, op.bitrow, op.bitcol, -op.amplitude)

function (*)(lhs ::S1, rhs ::PureOperator{S2, BR}) where {S1<:Number, S2<:Number, BR}
  S = promote_type(S1, S2)
  return PureOperator{S, BR}(rhs.hilbert_space, rhs.bitmask, rhs.bitrow, rhs.bitcol, lhs * rhs.amplitude)
end

function (*)(lhs ::PureOperator{S1, BR}, rhs ::S2) where {S1<:Number, S2<:Number, BR}
  S = promote_type(S1, S2)
  return PureOperator{S, BR}(lhs.hilbert_space, lhs.bitmask, lhs.bitrow, lhs.bitcol, lhs.amplitude * rhs)
end

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


import Base.real, Base.imag, Base.conj, Base.transpose, Base.adjoint

real(arg ::PureOperator{R, BR}) where {R<:Real, BR} = arg
imag(arg ::PureOperator{R, BR}) where {R<:Real, BR} = PureOperator{R, BR}(arg.hilbert_space, arg.bitmask, arg.bitrow, arg.bitcol, zero(R))
conj(arg ::PureOperator{R, BR}) where {R<:Real, BR} = arg

real(arg ::PureOperator{Complex{R}, BR}) where {R<:Real, BR} = PureOperator{R, BR}(arg.hilbert_space, arg.bitmask, arg.bitrow, arg.bitcol, real(arg.amplitude))
imag(arg ::PureOperator{Complex{R}, BR}) where {R<:Real, BR} = PureOperator{R, BR}(arg.hilbert_space, arg.bitmask, arg.bitrow, arg.bitcol, imag(arg.amplitude))
conj(arg ::PureOperator{Complex{R}, BR}) where {R<:Real, BR} = PureOperator{Complex{R}, BR}(arg.hilbert_space, arg.bitmask, arg.bitrow, arg.bitcol, conj(arg.amplitude))

transpose(arg ::PureOperator{S, BR}) where {S, BR} = PureOperator{S, BR}(arg.hilbert_space, arg.bitmask, arg.bitcol, arg.bitrow, arg.amplitude)
adjoint(arg ::PureOperator{S, BR}) where {S, BR} = PureOperator{S, BR}(arg.hilbert_space, arg.bitmask, arg.bitcol, arg.bitrow, conj(arg.amplitude))


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
                       istate_source ::Integer,
                       istate_target ::Integer,
                       amplitude::Number=1;
                       dtype ::DataType=UInt)
  site = hilbert_space.sites[isite]
  state_source = site.states[istate_source]
  state_target = site.states[istate_target]

  bm = get_bitmask(hilbert_space, isite; dtype=dtype)
  bs = dtype(istate_source - 1) << hilbert_space.bitoffsets[isite]
  bt = dtype(istate_target - 1) << hilbert_space.bitoffsets[isite]
  return PureOperator{typeof(amplitude), dtype}(hilbert_space, bm, bs, bt, amplitude)
end