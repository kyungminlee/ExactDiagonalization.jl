export symmetry_apply
export is_invariant

import TightBindingLattice.AbstractSymmetryOperation
import TightBindingLattice.AbstractSymmetryGroup


## AbstractSymmetryOperation

### HilbertSpaceSector
function symmetry_apply(hss ::HilbertSpaceSector{QN}, symop ::AbstractSymmetryOperation, args...; kwargs...) where {QN}
  return symmetry_apply(hss.parent, symop, args...; kwargs...)
end

function is_invariant(hss ::HilbertSpaceSector{QN}, symop ::AbstractSymmetryOperation, args...; kwargs...) where {QN}
  return is_invariant(hss.parent, symop, args...; kwargs...)
end

function is_invariant(hss::HilbertSpaceSector{QN}, symgroup ::AbstractSymmetryGroup, args...; kwargs...) where {QN}
  return is_invariant(hss.parent, symgroup, args...; kwargs...)
end

### generic symmetry operations for NullOperator and SumOperator
function symmetry_apply(hs ::HilbertSpace{QN}, symop ::AbstractSymmetryOperation, op::NullOperator) where {QN}
  return op
end

function symmetry_apply(hs::HilbertSpace{QN}, symop ::AbstractSymmetryOperation, op::SumOperator{S, BR}) where {QN, S<:Number, BR<:Unsigned}
  terms = collect(symmetry_apply(hs, symop, t) for t in op.terms)
  return SumOperator{S, BR}(terms)
end

## Permutation
### Binary Representation
function symmetry_apply(hs::HilbertSpace{QN}, permutation ::Permutation, bitrep ::BR) ::BR where {QN, BR<:Unsigned}
  out = zero(BR)
  for (i, j) in enumerate(permutation.map)
    out |= ( ( (bitrep >> hs.bitoffsets[i]) & make_bitmask(hs.bitwidths[i]) ) << hs.bitoffsets[j] )
  end
  return out
end

### Operator
function symmetry_apply(hs::HilbertSpace{QN}, permutation ::Permutation, op::PureOperator{S, BR}) where {QN, S<:Number, BR<:Unsigned}
  bm = symmetry_apply(hs, permutation, op.bitmask)
  br = symmetry_apply(hs, permutation, op.bitrow)
  bc = symmetry_apply(hs, permutation, op.bitcol)
  am = op.amplitude
  return PureOperator{S, BR}(bm, br, bc, am)
end


## isinvariant
function is_invariant(hs::HilbertSpace{QN}, symop ::AbstractSymmetryOperation, op::AbstractOperator) where {QN}
  return simplify(op - symmetry_apply(hs, symop, op)) == NullOperator()
end

function is_invariant(hs::HilbertSpace{QN}, symgroup ::AbstractSymmetryGroup, op::AbstractOperator) where {QN}
  return all(is_invariant(hs, g, op) for g in symgroup.generators)
end
