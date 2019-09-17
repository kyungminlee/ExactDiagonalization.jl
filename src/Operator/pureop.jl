#abstract type AbstractOperator end
#struct AbstractHilbertSpace end

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
export PureOperator
export pure_operator
export OptionalPureOperator


struct PureOperator{Scalar<:Number, BR<:Unsigned} <:AbstractOperator
  hilbert_space ::AbstractHilbertSpace
  bitmask ::BR
  bitsource ::BR  # Row
  bittarget ::BR  # Column
  amplitude ::Scalar
  # TODO: fermion sign

  function PureOperator{S, BR}(hilbert_space ::AbstractHilbertSpace,
                               bitmask, bitsource, bittarget, amplitude::S) where {S, BR}
    if (~bitmask) & bitsource != 0x0
      throw(ArgumentError("every bit of bitsource not in bitmask should be set to zero"))
    elseif (~bitmask) & bittarget != 0x0
      throw(ArgumentError("every bit of bittarget not in bitmask should be set to zero"))
    end
    return new{S, BR}(hilbert_space, bitmask, bitsource, bittarget, amplitude)
  end
end

function prettyprintln(op::PureOperator{S, BR}; prefix::AbstractString="") where {S, BR}
  println(prefix, "PureOperator")
  println(prefix, "M: ", string(op.bitmask, base=2, pad=op.hilbert_space.bitoffsets[end]))
  println(prefix, "S: ", string(op.bitsource, base=2, pad=op.hilbert_space.bitoffsets[end]))
  println(prefix, "T: ", string(op.bittarget, base=2, pad=op.hilbert_space.bitoffsets[end]))
end

const OptionalPureOperator{Scalar, BR} = Union{PureOperator{Scalar, BR}, NullOperator} where {Scalar, BR}

(-)(op ::PureOperator{S, BR}) where {S, BR} = PureOperator{S, BR}(op.hilbert_space, op.bitmask, op.bitsource, op.bittarget, -op.amplitude)

function (*)(lhs ::S1, rhs ::PureOperator{S2, BR}) where {S1<:Number, S2<:Number, BR}
  S = promote_type(S1, S2)
  PureOperator{S, BR}(rhs.hilbert_space, rhs.bitmask, rhs.bitsource, rhs.bittarget, lhs * rhs.amplitude)
end

function (*)(lhs ::PureOperator{S1, BR}, rhs ::S2) where {S1<:Number, S2<:Number, BR}
  S = promote_type(S1, S2)
  PureOperator{S, BR}(lhs.hilbert_space, lhs.bitmask, lhs.bitsource, lhs.bittarget, lhs.amplitude * rhs)
end

function (*)(lhs ::PureOperator{S1, BR}, rhs ::PureOperator{S2, BR}) where {S1<:Number, S2<:Number, BR}
  S3 = promote_type(S1, S2)

  onlylhs_bitmask   =   lhs.bitmask  & (~rhs.bitmask)
  onlyrhs_bitmask   = (~lhs.bitmask) &   rhs.bitmask
  intersect_bitmask =   lhs.bitmask  &   rhs.bitmask
  union_bitmask     =   lhs.bitmask  |   rhs.bitmask

  if (lhs.bittarget & intersect_bitmask) != (rhs.bitsource & intersect_bitmask)
    return NullOperator()
  else
    new_bitmask = union_bitmask
    new_bitsource = lhs.bitsource | (rhs.bitsource & onlyrhs_bitmask)
    new_bittarget = rhs.bittarget | (lhs.bittarget & onlylhs_bitmask)
    new_amplitude = lhs.amplitude * rhs.amplitude
    if isapprox(new_amplitude, 0)
      return NullOperator()
    else
      return PureOperator{S3, BR}(lhs.hilbert_space, 
                                  new_bitmask,
                                  new_bitsource,
                                  new_bittarget,
                                  new_amplitude)
    end
  end
end

real(arg ::PureOperator{S, BR}) where {S<:Real, BR} = arg
imag(arg ::PureOperator{S, BR}) where {S<:Real, BR} = NullOperator()

function real(arg ::PureOperator{Complex{R}, BR}) where {R<:Real, BR}
  return PureOperator{R, BR}(arg.hilbert_space,
                             arg.bitmask,
                             arg.bitsource,
                             arg.bittarget,
                             real(arg.amplitude))
end

function imag(arg ::PureOperator{Complex{R}, BR}) where {R<:Real, BR}
  return PureOperator{R, BR}(arg.hilbert_space,
                             arg.bitmask,
                             arg.bitsource,
                             arg.bittarget,
                             imag(arg.amplitude))
end


function isless(lhs ::PureOperator{S, BR}, rhs ::PureOperator{S, BR}) where {S, BR}
  if lhs.bitmask < rhs.bitmask
    return true
  elseif lhs.bitmask > rhs.bitmask
    return false
  end

  if lhs.bitsource < rhs.bitsource
    return true
  elseif lhs.bitsource > rhs.bitsource
    return false
  end

  if lhs.bittarget < rhs.bittarget
    return true
  elseif lhs.bittarget > rhs.bittarget
    return false
  end

  return false
end

import Base.promote_rule
function promote_rule(lhs::Type{PureOperator{S1, BR}}, rhs::Type{PureOperator{S2, BR}}) where {S1, S2, BR}
  S3 = promote_type(S1, S2)
  return PureOperator{S3, BR}
end

import Base.convert
function convert(type ::Type{PureOperator{S1, BR}}, obj::PureOperator{S2, BR}) where {S1, S2, BR}
  return PureOperator{S1, BR}(obj.hilbert_space, obj.bitmask, obj.bitsource, obj.bittarget, convert(S1, obj.amplitude))
end

function pure_operator(hilbert_space ::AbstractHilbertSpace,
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