export apply, apply!, apply_unsafe!

"""
    apply!

Apply operator to `psi` and add it to `out`.
"""
function apply!(out::SparseState{S1, BR},
                nullop ::NullOperator,
                psi::SparseState{S2, BR}) where {S1, S2, BR}
  return out
end


function apply!(out::SparseState{S1, BR},
                psi::SparseState{S2, BR},
                nullop ::NullOperator) where {S1, S2, BR}
  return out
end


function apply!(out::SparseState{S1, BR},
                pureop ::PureOperator{S2, BR},
                psi::SparseState{S3, BR}) where {S1, S2, S3, BR}
  for (b, v) in psi.components
    if (b & pureop.bitmask) == pureop.bitcol
      b2 = (b & ~pureop.bitmask) | pureop.bitrow
      out[b2] += pureop.amplitude * v
    end
  end
  out
end


function apply!(out::SparseState{S1, BR},
                psi::SparseState{S3, BR},
                pureop ::PureOperator{S2, BR}) where {S1, S2, S3, BR}
  for (b, v) in psi.components
    if (b & pureop.bitmask) == pureop.bitrow
      b2 = (b & ~pureop.bitmask) | pureop.bitcol
      out[b2] += v * pureop.amplitude
    end
  end
  out
end


function apply!(out::SparseState{S1, BR},
                sumop ::SumOperator{S2, BR},
                psi::SparseState{S3, BR}) where {S1, S2, S3, BR}
  for t in sumop.terms
    apply!(out, t, psi)
  end
  out
end


function apply!(out::SparseState{S1, BR},
                psi::SparseState{S3, BR},
                sumop ::SumOperator{S2, BR}) where {S1, S2, S3, BR}
  for t in sumop.terms
    apply!(out, psi, t)
  end
  out
end


function apply(pureop ::NullOperator, psi::SparseState{S2, BR}) where {S2, BR}
  return SparseState{S2, BR}()
end


function apply(psi::SparseState{S2, BR}, pureop ::NullOperator) where {S2, BR}
  return SparseState{S2, BR}()
end


function apply(pureop ::PureOperator{S1, BR}, psi::SparseState{S2, BR}) where {S1, S2, BR}
  S3 = promote_type(S1, S2)
  out = SparseState{S3, BR}()
  apply!(out, pureop, psi)
  return out
end


function apply(psi::SparseState{S2, BR}, pureop ::PureOperator{S1, BR}) where {S1, S2, BR}
  S3 = promote_type(S1, S2)
  out = SparseState{S3, BR}()
  apply!(out, psi, pureop)
  return out
end


function apply(sumop ::SumOperator{S1, BR}, psi::SparseState{S2, BR}) where {S1, S2, BR}
  S3 = promote_type(S1, S2)
  out = SparseState{S3, BR}()
  for t in sumop.terms
    apply!(out, t, psi)
  end
  return out
end


function apply(psi::SparseState{S2, BR}, sumop ::SumOperator{S1, BR}) where {S1, S2, BR}
  S3 = promote_type(S1, S2)
  out = SparseState{S3, BR}()
  for t in sumop.terms
    apply!(out, psi, t)
  end
  return out
end
