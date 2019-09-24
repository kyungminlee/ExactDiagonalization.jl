
# """
# Returns dict
# """
# function apply(hs ::HilbertSpace,
#                op ::KroneckerProductOperator{OS},
#                binrep ::BR) where {OS<:Number, BR<:Unsigned}
#   zero_scalar = zero(OS)
#   output = Dict{BR, OS}(binrep => 1.0)
#   for (isite, siteop) in op.operators
#     @assert 1 <= isite <= length(hs.sites)
#     nextoutput = DefaultDict{BR, OS}(zero_scalar)
#     sitedim = length(hs.sites[isite].states)
#     for (r ::BR, amplitude ::OS) in output
#       sri = get_state_index(hs, r, isite) # site row index
#       for sci in 1:sitedim # site col index
#         value = siteop[sri, sci]
#         if ! isapprox(value, 0)
#           c = update(hs, r, isite, sci)
#           nextoutput[c] += value * amplitude
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


# function apply(hs ::HilbertSpace,
#                ops ::AbstractArray{KroneckerProductOperator{OS}},
#                binrep ::BR) where {OS<:Number, BR <:Unsigned}
#   zero_scalar = zero(OS)
#   results = DefaultDict{BR, OS}(zero_scalar)
#   for op in ops
#     output = apply(hs, op, binrep)
#     for (row, amplitude) in output
#       results[row] += amplitude
#     end
#   end
#   return results
# end


# function apply(hs ::HilbertSpace,
#                ops ::AbstractArray{KroneckerProductOperator{OS}},
#                psi ::AbstractDict{BR, SS}) where {OS<:Number, SS<:Number, BR<:Unsigned}
#   OutScalar = promote_type(OS, SS)
#   zero_scalar = zero(OutScalar)
#   results = DefaultDict{BR, OutScalar}(zero_scalar)

#   for (binrep, amplitude) in psi
#     output = apply(hs, ops, binrep)
#     for (ϕ, α) in output
#       results[binrep] += amplitude * α 
#     end
#   end
#   return results
# end