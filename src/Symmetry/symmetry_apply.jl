export symmetry_apply
export is_invariant

## Hilbert Space

function symmetry_apply(hs::HilbertSpace, permutation ::Permutation, bitrep ::BR) where {BR}
  out = zero(BR)
  for (i, j) in enumerate(permutation.map)
    out |= ( ( (bitrep >> hs.bitoffsets[i]) & make_bitmask(hs.bitwidths[i]) ) << hs.bitoffsets[j] )
  end
  return out
end

@inline function symmetry_apply(hss::HilbertSpaceSector, permutation ::Permutation, bitrep ::BR) where {BR}
  return symmetry_apply(hss.parent, permutation, bitrep)
end

## Operator

function symmetry_apply(permutation ::Permutation, op::NullOperator)
  return op
end

function symmetry_apply(permutation ::Permutation, op::PureOperator{S, BR}) where {S,BR}
  hs = op.hilbert_space
  bm = symmetry_apply(op.hilbert_space, permutation, op.bitmask)
  br = symmetry_apply(op.hilbert_space, permutation, op.bitrow)
  bc = symmetry_apply(op.hilbert_space, permutation, op.bitcol)
  am = op.amplitude
  return PureOperator{S, BR}(hs, bm, br, bc, am)
end

function symmetry_apply(permutation ::Permutation, op::SumOperator{S, BR}) where {S, BR}
  terms = collect(symmetry_apply(permutation, t) for t in op.terms)
  return SumOperator{S, BR}(op.hilbert_space, terms)
end


### isinvariant

function is_invariant(permutation ::Permutation, op::AbstractOperator)
  return simplify(op - symmetry_apply(permutation, op)) == NullOperator()
end

function is_invariant(trans_group ::TranslationGroup, op::AbstractOperator)
  return all(
    is_invariant(g, op) for g in trans_group.generators
  )
end
