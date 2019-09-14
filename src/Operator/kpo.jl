export KroneckerProductOperator
export clean!


struct KroneckerProductOperator{Scalar <:Number}
  hilbert_space ::AbstractHilbertSpace
  amplitude ::Scalar
  operators ::Dict{Int, Matrix{Scalar}}
end

KPO = KroneckerProductOperator

function clean!(op ::KPO; tol=sqrt(eps(Float64)))
  keys_to_delete = [k for (k, v) in op.operators if isapprox(v, I)]
  for k in keys_to_delete
    delete!(k, op)
  end
end

import Base.*

"""
O3 = O1 * O2
"""
function *(lhs ::KPO{S1}, rhs ::KPO{S2}) where {S1<:Number, S2<:Number}
  @assert(lhs.hilbert_space == rhs.hilbert_space)
  S3 = promote_type(S1, S2)

  kl = keys(lhs.operators)
  kr = keys(rhs.operators)
  common = intersect(kl, kr)
  complete = union(kl, kr)
  only_lhs = setdiff(kl, kr)
  only_rhs = setdiff(kr, kl)

  output_operators = Dict{Int, Matrix{S3}}()

  for k in common
    @assert size(lhs.operators[k]) == size(rhs.operators[k])
    L = lhs.operators[k]
    R = rhs.operators[k]
    output_operators[k] = L * R
  end
  for k in only_lhs
    output_operators[k] = lhs.operators[k]
  end
  for k in only_rhs
    output_operators[k] = rhs.operators[k]
  end
  return KPO{S3}(lhs.hilbert_space, lhs.amplitude * rhs.amplitude, output_operators)
end

"""
O2 = O1 * 0.1
"""
function *(lhs ::KPO{S1}, rhs::S2) where {S1<:Number, S2<:Number}
  S3 = promote_type(S1, S2)
  return KPO{S3}(lhs.hilbert_space, lhs.amplitude * rhs, lhs.operators)
end

"""
O2 = 0.1 * O1
"""
function *(lhs ::S1, rhs::KPO{S2}) where {S1<:Number, S2<:Number}
  S3 = promote_type(S1, S2)
  return KPO{S3}(rhs.hilbert_space, lhs * rhs.amplitude, rhs.operators)
end


"""
<ψ'| = <ψ| O 
"""
function *(lhs::SparseState{BinRep, SS1}, rhs ::KPO{OS}) where {OS<:Number, BinRep, SS1 <:Number}
  OutScalar = promote_type(OS, SS)
  SS = SparseState{BinRep, OutScalar}

  @assert lhs.hilbert_space == rhs.hilbert_space
  hs = lhs.hilbert_space

  ψ = lhs
  for (isite, site_op) in rhs.operators
    @assert 1 <= isite <= length(hs.sites)
    ψp = SS(lhs.hilbert_space)
    site_dim = dimension(hs.sites[isite])
    for (r ::BinRep, amplitude) in ψ.components
      sri = get_state_index(hs, r, isite)
      for sci in 1:site_dim
        value = site_op[sri, sci]
        if ! isapprox(value, 0)
          c = update(hs, r, isite, sci)
          ψp[c] += value * amplitude
        end
      end
    end
    ψ = ψp
  end
  for (k,v) in ψ.components
    ψ.components[k] = v * rhs.amplitude
  end
  return ψ
end


