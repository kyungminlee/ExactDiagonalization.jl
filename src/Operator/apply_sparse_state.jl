export apply, apply!, apply_unsafe!

function apply_unsafe!(out::SparseState{S1, BR}, nullop ::NullOperator, psi::SparseState{S2, BR}) where {S1, S2, BR}
  return out
end

function apply_unsafe!(out::SparseState{S1, BR}, psi::SparseState{S2, BR}, nullop ::NullOperator) where {S1, S2, BR}
  return out
end

function apply_unsafe!(out::SparseState{S1, BR}, pureop ::PureOperator{S2, BR}, psi::SparseState{S3, BR}) where {S1, S2, S3, BR}
  for (b, v) in psi.components
    if (b & pureop.bitmask) == pureop.bitcol
      b2 = (b & ~pureop.bitmask) | pureop.bitrow
      out[b2] += pureop.amplitude * v
    end
  end
  out
end

function apply_unsafe!(out::SparseState{S1, BR}, psi::SparseState{S3, BR}, pureop ::PureOperator{S2, BR}) where {S1, S2, S3, BR}
  for (b, v) in psi.components
    if (b & pureop.bitmask) == pureop.bitrow
      b2 = (b & ~pureop.bitmask) | pureop.bitcol
      out[b2] += v * pureop.amplitude
    end
  end
  out
end

function apply_unsafe!(out::SparseState{S1, BR}, sumop ::SumOperator{S2, BR}, psi::SparseState{S3, BR}) where {S1, S2, S3, BR}
  for t in sumop.terms
    apply_unsafe!(out, t, psi)
  end
  out
end

function apply_unsafe!(out::SparseState{S1, BR}, psi::SparseState{S3, BR}, sumop ::SumOperator{S2, BR}) where {S1, S2, S3, BR}
  for t in sumop.terms
    apply_unsafe!(out, psi, t)
  end
  out
end


"""
    apply!

Apply operator to `psi` and add it to `out`.
"""
function apply!(out::SparseState{S1, BR}, nullop ::NullOperator, psi::SparseState{S2, BR}) where {S1, S2, BR}
  return out
end

function apply!(out::SparseState{S1, BR}, psi::SparseState{S2, BR}, nullop ::NullOperator) where {S1, S2, BR}
  return out
end

function apply!(out::SparseState{S1, BR}, pureop ::PureOperator{S2, BR}, psi::SparseState{S3, BR}) where {S1, S2, S3, BR}
  if pureop.hilbert_space !== psi.hilbert_space || out.hilbert_space !== psi.hilbert_space
    throw(ArgumentError("Hilbert spaces of lhs and rhs of + should match"))
  end
  apply_unsafe!(out, pureop, psi)
end

function apply!(out::SparseState{S1, BR}, psi::SparseState{S3, BR}, pureop ::PureOperator{S2, BR}) where {S1, S2, S3, BR}
  if pureop.hilbert_space !== psi.hilbert_space || out.hilbert_space !== psi.hilbert_space
    throw(ArgumentError("Hilbert spaces of lhs and rhs of + should match"))
  end
  apply_unsafe!(out, psi, pureop)
end


function apply!(out::SparseState{S1, BR}, sumop ::SumOperator{S2, BR}, psi::SparseState{S3, BR}) where {S1, S2, S3, BR}
  if sumop.hilbert_space !== psi.hilbert_space || out.hilbert_space !== psi.hilbert_space
    throw(ArgumentError("Hilbert spaces should match"))
  end
  apply_unsafe!(out, sumop, psi)
end

function apply!(out::SparseState{S1, BR}, psi::SparseState{S3, BR}, sumop ::SumOperator{S2, BR}) where {S1, S2, S3, BR}
  if sumop.hilbert_space !== psi.hilbert_space || out.hilbert_space !== psi.hilbert_space
    throw(ArgumentError("Hilbert spaces should match"))
  end
  apply_unsafe!(out, psi, sumop)
end


function apply(pureop ::NullOperator, psi::SparseState{S2, BR}) where {S2, BR}
  return SparseState{S2, BR}(psi.hilbert_space)
end

function apply(psi::SparseState{S2, BR}, pureop ::NullOperator) where {S2, BR}
  return SparseState{S2, BR}(psi.hilbert_space)
end


function apply(pureop ::PureOperator{S1, BR}, psi::SparseState{S2, BR}) where {S1, S2, BR}
  S3 = promote_type(S1, S2)
  if pureop.hilbert_space !== psi.hilbert_space
    throw(ArgumentError("Hilbert spaces of lhs and rhs of + should match"))
  end
  out = SparseState{S3, BR}(psi.hilbert_space)
  apply_unsafe!(out, pureop, psi)
  return out
end

function apply(psi::SparseState{S2, BR}, pureop ::PureOperator{S1, BR}) where {S1, S2, BR}
  S3 = promote_type(S1, S2)
  if pureop.hilbert_space !== psi.hilbert_space
    throw(ArgumentError("Hilbert spaces of lhs and rhs of + should match"))
  end
  out = SparseState{S3, BR}(psi.hilbert_space)
  apply_unsafe!(out, psi, pureop)
  return out
end


function apply(sumop ::SumOperator{S1, BR}, psi::SparseState{S2, BR}) where {S1, S2, BR}
  S3 = promote_type(S1, S2)
  if sumop.hilbert_space !== psi.hilbert_space
    throw(ArgumentError("Hilbert spaces of lhs and rhs of + should match"))
  end
  out = SparseState{S3, BR}(psi.hilbert_space)
  for t in sumop.terms
    apply_unsafe!(out, t, psi)
  end
  return out
end


function apply(psi::SparseState{S2, BR}, sumop ::SumOperator{S1, BR}) where {S1, S2, BR}
  S3 = promote_type(S1, S2)
  if sumop.hilbert_space !== psi.hilbert_space
    throw(ArgumentError("Hilbert spaces of lhs and rhs of + should match"))
  end
  out = SparseState{S3, BR}(psi.hilbert_space)
  for t in sumop.terms
    apply_unsafe!(out, psi, t)
  end
  return out
end
