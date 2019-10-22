#abstract type AbstractOperator end
#struct HilbertSpace end

export FermionPureOperator
export fermion_pure_operator

NOT IMPLEMENTED

#=
# TODO: Fermion operator should be compatible with
# Symmetry operation, (i.e. permutation between different orbitals)
# Fermion mask does not work very well.
# e.g.
# Sites :  1 -- 2 -- 3 -- 4 --
# Fermions:
# c†(1) :  +
# c†(2) :  z -- +
# c†(3) :  z -- z -- +
# c†(4) :  z -- z -- z -- +
c†(4)c(2): . -- - -- z -- +

e.g.)

c†(3) c(0)
  543210
r __1__0
c __0__1
m __1__1
f ___11_
F 111111 (all the fermions)

1. go to every bit of m,
2. check whether fmask turns on or off.
   0: left off, right 0
   3: left 0, right on

3. every on-bit of m is mapped to a different bit
   take the mapping with it.

e.g. add 4:
   0 -> 4
   3 -> 1

c†(1) c(4) = - c(4) c†(1)

  543210
r _0__1_
c _1__0_
m _1__1_
f 1____1
   (or)
f __11__
F 111111 (all the fermions)

   0->4: left off, right 0
   3->1: left 0, right on

# let "normal ordering" be such that
# the top bit of f starts with 0

Fp = F & (~m) :
Fm = F & m    :

Fp ^ f <=> f' : ordering possibilities (always choose the smaller fmask)
count_ones(Fm & r) == count_ones(Fm & c)  (mod 2) (fermion parity)

=#
struct FermionPureOperator{Scalar<:Number, BR<:Unsigned} <:AbstractOperator
  hilbert_space ::HilbertSpace
  bitmask ::BR
  bitfmask ::BR
  bitrow ::BR  # Row
  bitcol ::BR  # Column
  amplitude ::Scalar

  function FermionPureOperator{S, BR}(hilbert_space ::HilbertSpace,
                               bitmask, bitfmask, bitrow, bitcol, amplitude::S) where {S, BR}
    if (~bitmask) & bitrow != 0x0
      throw(ArgumentError("every bit of bitrow not in bitmask should be set to zero"))
    elseif (~bitmask) & bitcol != 0x0
      throw(ArgumentError("every bit of bitcol not in bitmask should be set to zero"))
    elseif bitmask & bitfmask != 0x0
      throw(ArgumentError("bitmask and bitfermionmask should not overlap"))
    end

    return new{S, BR}(hilbert_space, bitmask, bitfmask, bitrow, bitcol, amplitude)
  end
end

import Base.==

function (==)(lhs ::FermionPureOperator{S1, BR}, rhs::FermionPureOperator{S2, BR}) where {S1, S2, BR}
  return ((lhs.hilbert_space == rhs.hilbert_space) &&
          (lhs.bitmask == rhs.bitmask) &&
          (lhs.bitfermionmask == rhs.bitfermionmask) &&
          (lhs.bitrow == rhs.bitrow) &&
          (lhs.bitcol == rhs.bitcol) &&
          (lhs.amplitude == rhs.amplitude))
end

(-)(op ::FermionPureOperator{S, BR}) where {S, BR} = FermionPureOperator{S, BR}(op.hilbert_space, op.bitmask, op.bitfermionmask, op.bitrow, op.bitcol, -op.amplitude)

function (*)(lhs ::S1, rhs ::FermionPureOperator{S2, BR}) where {S1<:Number, S2<:Number, BR}
  S = promote_type(S1, S2)
  return FermionPureOperator{S, BR}(rhs.hilbert_space, rhs.bitmask, rhs.bitfermionmask, rhs.bitrow, rhs.bitcol, lhs * rhs.amplitude)
end

function (*)(lhs ::FermionPureOperator{S1, BR}, rhs ::S2) where {S1<:Number, S2<:Number, BR}
  S = promote_type(S1, S2)
  return FermionPureOperator{S, BR}(lhs.hilbert_space, lhs.bitmask, lhs.bitfermionmask, lhs.bitrow, lhs.bitcol, lhs.amplitude * rhs)
end

function (*)(lhs ::FermionPureOperator{S1, BR}, rhs ::FermionPureOperator{S2, BR}) where {S1<:Number, S2<:Number, BR}
  if lhs.hilbert_space !== rhs.hilbert_space
    throw(ArgumentError("Hilbert spaces don't match"))
  end
  S3 = promote_type(S1, S2)

  onlylhs_bitmask   =   lhs.bitmask  & (~rhs.bitmask)
  onlyrhs_bitmask   = (~lhs.bitmask) &   rhs.bitmask
  intersect_bitmask =   lhs.bitmask  &   rhs.bitmask
  union_bitmask     =   lhs.bitmask  |   rhs.bitmask


  fm1 = lhs.bitfermionmask ^ rhs.bitfermionmask
  fm2 = lhs.bitfermion

  fm_left1 = lhs.bitfermionmask & (~rhs.bitmask)
  fm_left2 = lhs.bitfermionmask & ( rhs.bitmask)  # compare this with lhs.bitcol

  fm_right1 = (~lhs.bitmask) & rhs.bitfermionmask
  fm_right2 = ( lhs.bitmask) & rhs.bitfermionmask # compare this with rhs.source

  new_bitfermionmask = fm_left1 ⊻ fm_left2

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

import Base.real, Base.imag, Base.conj, Base.transpose

real(arg ::PureOperator{S, BR}) where {S<:Real, BR} = arg
imag(arg ::PureOperator{S, BR}) where {S<:Real, BR} = PureOperator{S, BR}(arg.hilbert_space, arg.bitmask, arg.bitrow, arg.bitcol, zero(S))

function real(arg ::PureOperator{Complex{R}, BR}) where {R<:Real, BR}
  return PureOperator{R, BR}(arg.hilbert_space,
                             arg.bitmask,
                             arg.bitrow,
                             arg.bitcol,
                             real(arg.amplitude))
end

function imag(arg ::PureOperator{Complex{R}, BR}) where {R<:Real, BR}
  return PureOperator{R, BR}(arg.hilbert_space,
                             arg.bitmask,
                             arg.bitrow,
                             arg.bitcol,
                             imag(arg.amplitude))
end

function conj(arg ::PureOperator{S, BR}) where {S, BR}
  return PureOperator{S, BR}(arg.hilbert_space,
                             arg.bitmask,
                             arg.bitrow,
                             arg.bitcol,
                             conj(arg.amplitude))
end

function transpose(arg ::PureOperator{S, BR}) where {S, BR}
  return PureOperator{S, BR}(arg.hilbert_space,
                             arg.bitmask,
                             arg.bitcol, # switch
                             arg.bitrow, # order
                             arg.amplitude)
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
