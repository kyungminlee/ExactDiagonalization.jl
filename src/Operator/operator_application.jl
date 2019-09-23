export apply, apply!


import Base.isempty
isempty(psi::SparseState{S, BR}) where {S, BR} = isempty(psi.components)

function apply!(out::SparseState{S1, BR}, nullop ::NullOperator, psi::SparseState{S2, BR}) where {S1, S2, BR}
end


function apply!(out::SparseState{S1, BR}, pureop ::PureOperator{S2, BR}, psi::SparseState{S3, BR}) where {S1, S2, S3, BR}
  if pureop.hilbert_space !== psi.hilbert_space || out.hilbert_space !== psi.hilbert_space
    throw(ArgumentError("Hilbert spaces of lhs and rhs of + should match"))
  end
  apply_unsafe!(out, pureop, psi)
end


function apply!(out::SparseState{S1, BR}, sumop ::SumOperator{S2, BR}, psi::SparseState{S3, BR}) where {S1, S2, S3, BR}
  if sumop.hilbert_space !== psi.hilbert_space || out.hilbert_space !== psi.hilbert_space 
    throw(ArgumentError("Hilbert spaces should match"))
  end
  return apply_unsafe!(out, sumop, psi)
end


function apply_unsafe!(out::SparseState{S1, BR}, pureop ::PureOperator{S2, BR}, psi::SparseState{S3, BR}) where {S1, S2, S3, BR}
  for (b, v) in psi.components
    if (b & pureop.bitmask) == pureop.bitsource
      b2 = (b & ~pureop.bitmask) | pureop.bittarget
      out[b2] += pureop.amplitude * v
    end
  end
end

function apply_unsafe!(out::SparseState{S1, BR}, sumop ::SumOperator{S2, BR}, psi::SparseState{S3, BR}) where {S1, S2, S3, BR}
  for t in sumop.terms
    apply!(out, t, psi)
  end
end


function apply(pureop ::NullOperator, psi::SparseState{S2, BR}) where {S2, BR}
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





# # TODO: Replace with more efficient
# function apply_naive(hs ::AbstractHilbertSpace,
#                op ::KroneckerProductOperator{OS},
#                binrep ::U) where {OS <:Number, U <: Unsigned}
#   intermediate = apply(hs, op, extract(hs, binrep))
#   out = Dict{U, OS}()
#   for (k, v) in intermediate
#     out[compress(hs, k; BR=U)] = v
#   end
#   return out
# end

# # TODO: Replace with more efficient
# function apply_naive(hs ::AbstractHilbertSpace,
#                ops ::AbstractArray{KroneckerProductOperator{OS}},
#                binrep ::U) where {OS<:Number, U <:Unsigned}
#   intermediate = apply(hs, ops, extract(hs, binrep))
#   out = Dict{U, OS}()
#   for (k, v) in intermediate
#     out[compress(hs, k; BR=U)] = v
#   end
#   return out
# end



# """
# Returns dict
# """
# function apply(hs ::AbstractHilbertSpace,
#                op ::KroneckerProductOperator{OS},
#                indexarray ::AbstractArray{I, 1}) where {OS <:Number, I<:Integer}
#   zero_scalar = zero(OS)

#   @assert length(indexarray) == length(hs.sites)
#   output = Dict(indexarray => 1.0)

#   for (isite, siteop) in op.operators
#     #@show isite
#     #@show length(hs.sites)
#     @assert 1 <= isite <= length(hs.sites)
#     nextoutput = DefaultDict{Vector{Int}, OS}(zero_scalar)
#     sitedim = length(hs.sites[isite].states)
#     for (indexarray, amplitude) in output
#       rowindex = indexarray[isite]
#       for colindex in 1:sitedim
#         value = siteop[rowindex, colindex]
#         if ! isapprox(value, 0)
#           newindexarray = copy(indexarray)
#           newindexarray[isite] = colindex
#           nextoutput[newindexarray] += value * amplitude
#         end
#       end
#     end
#     output = nextoutput
#   end

#   for (k, v) in output
#     output[k] = v * op.amplitude
#   end

#   return output
# end

# function apply(hs ::AbstractHilbertSpace,
#                ops ::AbstractArray{KroneckerProductOperator{S}},
#                indexarray ::AbstractArray{I, 1}) where {S<:Number, I<:Integer}
#   zero_scalar = zero(S)
  
#   @assert length(indexarray) == length(hs.sites)
#   results = DefaultDict{Vector{Int}, S}(zero_scalar)
#   for op in ops
#     output = apply(hs, op, indexarray)
#     for (indexarray, amplitude) in output
#       results[indexarray] += amplitude
#     end
#   end
#   return results
# end

# function apply(hs ::AbstractHilbertSpace,
#                op ::KroneckerProductOperator{OS},
#                psi ::AbstractDict{Vector{Int}, SS}) where {OS<:Number, SS<:Number}
#   OutScalar = promote_type(OS, SS)
#   zero_scalar = zero(OutScalar)
#   results = DefaultDict{Vector{Int}, OutScalar}(zero_scalar)
#   for (indexarray, amplitude) in psi
#     output = apply(hs, op, indexarray)
#     for (ϕ, α) in output
#       results[ϕ] += amplitude * α
#     end
#   end
#   return results
# end

# function apply(hs ::AbstractHilbertSpace,
#                ops ::AbstractArray{KroneckerProductOperator{OS}},
#                psi ::AbstractDict{Vector{Int}, SS}) where {OS<:Number, SS<:Number}
#   OutScalar = promote_type(OS, SS)
#   zero_scalar = zero(OutScalar)
#   results = DefaultDict{Vector{Int}, OutScalar}(zero_scalar)

#   for (indexarray, amplitude) in psi
#     output = apply(hs, ops, indexarray)
#     for (ϕ, α) in output
#       results[indexarray] += amplitude * α 
#     end
#   end
#   return results
# end
