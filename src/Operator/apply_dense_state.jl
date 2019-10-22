export apply, apply!, apply_unsafe!


function apply_unsafe!(out::DenseState{S1, QN, BR}, nullop ::NullOperator, psi::DenseState{S2, QN, BR}) where {S1, S2, QN, BR}
  return out
end

function apply_unsafe!(out::DenseState{S1, QN, BR}, psi::DenseState{S2, QN, BR}, nullop ::NullOperator) where {S1, S2, QN, BR}
  return out
end

function apply_unsafe!(out::DenseState{S1, QN, BR}, pureop ::PureOperator{S2, BR}, psi::DenseState{S3, QN, BR}) where {S1, S2, S3, QN, BR}
  for (b, v) in zip(psi.hilbert_space_realization.basis_list, psi.components)
    if (b & pureop.bitmask) == pureop.bitcol
      b2 = (b & ~pureop.bitmask) | pureop.bitrow
      out[b2] += pureop.amplitude * v
    end
  end
  out
end

function apply_unsafe!(out::DenseState{S1, QN, BR}, psi::DenseState{S3, QN, BR}, pureop ::PureOperator{S2, BR}) where {S1, S2, S3, QN, BR}
  for (b, v) in zip(psi.hilbert_space_realization.basis_list, psi.components)
    if (b & pureop.bitmask) == pureop.bitrow
      b2 = (b & ~pureop.bitmask) | pureop.bitcol
      out[b2] += v * pureop.amplitude
    end
  end
  out
end

function apply_unsafe!(out::DenseState{S1, QN, BR}, sumop ::SumOperator{S2, BR}, psi::DenseState{S3, QN, BR}) where {S1, S2, S3, QN, BR}
  for t in sumop.terms
    apply_unsafe!(out, t, psi)
  end
  out
end

function apply_unsafe!(out::DenseState{S1, QN, BR}, psi::DenseState{S3, QN, BR}, sumop ::SumOperator{S2, BR}) where {S1, S2, S3, QN, BR}
  for t in sumop.terms
    apply_unsafe!(out, psi, t)
  end
  out
end



"""
    apply!

Apply operator to `psi` and add it to `out`.
"""
function apply!(out::DenseState{S1, QN, BR}, nullop ::NullOperator, psi::DenseState{S2, QN, BR}) where {S1, S2, QN, BR}
  return out
end

function apply!(out::DenseState{S1, QN, BR}, psi::DenseState{S2, QN, BR}, nullop ::NullOperator) where {S1, S2, QN, BR}
  return out
end

function apply!(out::DenseState{S1, QN, BR}, pureop ::PureOperator{S2, BR}, psi::DenseState{S3, QN, BR}) where {S1, S2, S3, QN, BR}
  if out.hilbert_space_realization !== psi.hilbert_space_realization || pureop.hilbert_space !== psi.hilbert_space_realization.hilbert_space
    throw(ArgumentError("Hilbert spaces of lhs and rhs of + should match"))
  end
  apply_unsafe!(out, pureop, psi)
end

function apply!(out::DenseState{S1, QN, BR}, psi::DenseState{S3, QN, BR}, pureop ::PureOperator{S2, BR}) where {S1, S2, S3, QN, BR}
  if out.hilbert_space_realization !== psi.hilbert_space_realization || pureop.hilbert_space !== psi.hilbert_space_realization.hilbert_space
    throw(ArgumentError("Hilbert spaces of lhs and rhs of + should match"))
  end
  apply_unsafe!(out, psi, pureop)
end


function apply!(out::DenseState{S1, QN, BR}, sumop ::SumOperator{S2, BR}, psi::DenseState{S3, QN, BR}) where {S1, S2, S3, QN, BR}
  if out.hilbert_space_realization !== psi.hilbert_space_realization || sumop.hilbert_space !== psi.hilbert_space_realization.hilbert_space
    throw(ArgumentError("Hilbert spaces should match"))
  end
  apply_unsafe!(out, sumop, psi)
end

function apply!(out::DenseState{S1, QN, BR}, psi::DenseState{S3, QN, BR}, sumop ::SumOperator{S2, BR}) where {S1, S2, S3, QN, BR}
  if out.hilbert_space_realization !== psi.hilbert_space_realization || sumop.hilbert_space !== psi.hilbert_space_realization.hilbert_space
    throw(ArgumentError("Hilbert spaces should match"))
  end
  apply_unsafe!(out, psi, sumop)
end


function apply(pureop ::NullOperator, psi::DenseState{S2, QN, BR}) where {S2, QN, BR}
  return DenseState{S2, QN, BR}(psi.hilbert_space)
end

function apply(psi::DenseState{S2, QN, BR}, pureop ::NullOperator) where {S2, QN, BR}
  return DenseState{S2, QN, BR}(psi.hilbert_space)
end


function apply(pureop ::PureOperator{S1, BR}, psi::DenseState{S2, QN, BR}) where {S1, S2, QN, BR}
  if pureop.hilbert_space !== psi.hilbert_space_realiation.hilbert_space
    throw(ArgumentError("Hilbert spaces of lhs and rhs of + should match"))
  end
  S3 = promote_type(S1, S2)
  out = DenseState{S3, QN, BR}(psi.hilbert_space_realization)
  apply_unsafe!(out, pureop, psi)
  return out
end

function apply(psi::DenseState{S2, QN, BR}, pureop ::PureOperator{S1, BR}) where {S1, S2, QN, BR}
  if pureop.hilbert_space !== psi.hilbert_space_realiation.hilbert_space
    throw(ArgumentError("Hilbert spaces of lhs and rhs of + should match"))
  end
  S3 = promote_type(S1, S2)
  out = DenseState{S3, QN, BR}(psi.hilbert_space_realization)
  apply_unsafe!(out, psi, pureop)
  return out
end


function apply(sumop ::SumOperator{S1, BR}, psi::DenseState{S2, QN, BR}) where {S1, S2, QN, BR}
  if sumop.hilbert_space !== psi.hilbert_space_realiation.hilbert_space
    throw(ArgumentError("Hilbert spaces of lhs and rhs of + should match"))
  end
  S3 = promote_type(S1, S2)
  out = DenseState{S3, QN, BR}(psi.hilbert_space_realization)
  for t in sumop.terms
    apply_unsafe!(out, t, psi)
  end
  return out
end


function apply(psi::DenseState{S2, QN, BR}, sumop ::SumOperator{S1, BR}) where {S1, S2, QN, BR}
  S3 = promote_type(S1, S2)
  if sumop.hilbert_space !== psi.hilbert_space_realiation.hilbert_space
    throw(ArgumentError("Hilbert spaces of lhs and rhs of + should match"))
  end
  out = DenseState{S3, QN, BR}(psi.hilbert_space_realization)
  for t in sumop.terms
    apply_unsafe!(out, psi, t)
  end
  return out
end
