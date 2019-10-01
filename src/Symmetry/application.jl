export apply_symmetry
export is_invariant

function apply_symmetry(hs::HilbertSpace, permutation ::Permutation, bitrep ::BR) where {BR}
  out = zero(BR)
  for (i, j) in enumerate(permutation.map)
    out |= ( ( (bitrep >> hs.bitoffsets[i]) & make_bitmask(hs.bitwidths[i]) ) << hs.bitoffsets[j] )
  end
  return out
end


function apply_symmetry(permutation ::Permutation, op::NullOperator)
  return op
end


function apply_symmetry(permutation ::Permutation, op::PureOperator{S, BR}) where {S,BR}
  hs = op.hilbert_space
  bm = apply_symmetry(op.hilbert_space, permutation, op.bitmask)
  br = apply_symmetry(op.hilbert_space, permutation, op.bitrow)
  bc = apply_symmetry(op.hilbert_space, permutation, op.bitcol)
  am = op.amplitude
  return PureOperator{S, BR}(hs, bm, br, bc, am)
end


function apply_symmetry(permutation ::Permutation, op::SumOperator{S, BR}) where {S, BR}
  terms = collect(apply_symmetry(permutation, t) for t in op.terms)
  return SumOperator{S, BR}(op.hilbert_space, terms)
end


function is_invariant(permutation ::Permutation, op::AbstractOperator)
  return simplify(op - apply_symmetry(permutation, op)) == NullOperator()
end

function is_invariant(trans_group ::TranslationGroup, op::AbstractOperator)
  return all(
    is_invariant(g, op) for g in trans_group.generators
  )
end