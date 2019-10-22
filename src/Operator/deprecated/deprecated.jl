# export KroneckerProductOperator
# export clean!

# struct KroneckerProductOperator{Scalar<:Number} <:AbstractOperator
#   hilbert_space ::HilbertSpace
#   amplitude ::Scalar
#   operators ::Dict{Int, Matrix{Scalar}}

#   function KroneckerProductOperator(
#       hs ::HilbertSpace,
#       am ::S1,
#       ops::AbstractDict{I, Matrix{S2}}) where {S1<:Number, I<:Integer, S2<:Number}

#       @warn "KroneckerProductOperator is deprecated"

#       S3 = promote_type(S1, S2)
#       n_sites = length(hs.sites)

#       for (i_site, matrix) in ops
#         if !( 1 <= i_site <= n_sites)
#           throw(ArgumentError("site index $(i_site) is not within range"))
#         end
#       end

#       return new{S3}(hs, S3(am), Dict{Int, Matrix{S3}}(i=>m for (i,m) in ops))
#   end

#   function KroneckerProductOperator{S3}(
#     hs ::HilbertSpace,
#     am ::S1,
#     ops::AbstractDict{I, Matrix{S2}}) where {S1<:Number, I<:Integer, S2<:Number, S3<:Number}

#     @warn "KroneckerProductOperator is deprecated"

#     n_sites = length(hs.sites)

#     for (i_site, matrix) in ops
#       if !( 1 <= i_site <= n_sites)
#         throw(ArgumentError("site index $(i_site) is not within range"))
#       end
#     end

#     return new{S3}(hs, S3(am), Dict{Int, Matrix{S3}}(i=>m for (i,m) in ops))
# end
# end

# KPO = KroneckerProductOperator

# # function clean!(op ::KPO; tol=sqrt(eps(Float64)))
# #   keys_to_delete = [k for (k, v) in op.operators if isapprox(v, I)]
# #   for k in keys_to_delete
# #     delete!(k, op)
# #   end
# # end

# import Base.*

# """
# O3 = O1 * O2
# """
# function *(lhs ::KPO{S1}, rhs ::KPO{S2}) where {S1<:Number, S2<:Number}
#   @assert(lhs.hilbert_space == rhs.hilbert_space)
#   S3 = promote_type(S1, S2)

#   kl = keys(lhs.operators)
#   kr = keys(rhs.operators)
#   common = intersect(kl, kr)
#   complete = union(kl, kr)
#   only_lhs = setdiff(kl, kr)
#   only_rhs = setdiff(kr, kl)

#   output_operators = Dict{Int, Matrix{S3}}()

#   for k in common
#     @assert size(lhs.operators[k]) == size(rhs.operators[k])
#     L = lhs.operators[k]
#     R = rhs.operators[k]
#     output_operators[k] = L * R
#   end
#   for k in only_lhs
#     output_operators[k] = lhs.operators[k]
#   end
#   for k in only_rhs
#     output_operators[k] = rhs.operators[k]
#   end
#   return KPO{S3}(lhs.hilbert_space, lhs.amplitude * rhs.amplitude, output_operators)
# end

# """
# O2 = O1 * 0.1
# """
# function *(lhs ::KPO{S1}, rhs::S2) where {S1<:Number, S2<:Number}
#   S3 = promote_type(S1, S2)
#   return KPO{S3}(lhs.hilbert_space, lhs.amplitude * rhs, lhs.operators)
# end

# """
# O2 = 0.1 * O1
# """
# function *(lhs ::S1, rhs::KPO{S2}) where {S1<:Number, S2<:Number}
#   S3 = promote_type(S1, S2)
#   return KPO{S3}(rhs.hilbert_space, lhs * rhs.amplitude, rhs.operators)
# end


# """
# <ψ'| = <ψ| O
# """
# function *(lhs::SparseState{BR, SS1}, rhs ::KPO{OS}) where {OS<:Number, BR, SS1 <:Number}
#   OutScalar = promote_type(OS, SS1)
#   SS = SparseState{OutScalar, BR}

#   @assert lhs.hilbert_space == rhs.hilbert_space
#   hs = lhs.hilbert_space

#   ψ = lhs
#   for (isite, site_op) in rhs.operators
#     @assert 1 <= isite <= length(hs.sites)
#     ψp = SS(lhs.hilbert_space)
#     site_dim = dimension(hs.sites[isite])
#     for (r ::BR, amplitude) in ψ.components
#       sri = get_state_index(hs, r, isite)
#       for sci in 1:site_dim
#         value = site_op[sri, sci]
#         if ! isapprox(value, 0)
#           c = update(hs, r, isite, sci)
#           ψp[c] += value * amplitude
#         end
#       end
#     end
#     ψ = ψp
#   end
#   for (k,v) in ψ.components
#     ψ.components[k] = v * rhs.amplitude
#   end
#   return ψ
# end


# import Base.convert

# function convert(type::Type{SumOperator{S, BR}}, obj::KroneckerProductOperator{S}) where {S, BR}
#   hs = obj.hilbert_space
#   bm = zero(BR)
#   for (i, op) in obj.operators
#     bm = bm | make_bitmask(hs.bitoffsets[i+1], hs.bitoffsets[i]; dtype=BR)
#   end
#   isites = sort(collect(keys(obj.operators)))

#   nonzeros = Vector{Tuple{Int, Int, S}}[ [] for i in isites]
#   for (i, isite) in enumerate(isites)
#     mat = obj.operators[isite]
#     n = size(mat)[1]
#     for irow in 1:n, icol in 1:n
#       if ! isapprox(mat[irow, icol], 0)
#         push!(nonzeros[i], (irow, icol, mat[irow, icol]))
#       end
#     end
#   end

#   terms = PureOperator{S, BR}[]
#   for pairs in Iterators.product(nonzeros...)
#     #@show p
#     bm = zero(BR)
#     bs = zero(BR)
#     bt = zero(BR)
#     am = obj.amplitude
#     for (isite, p) in zip(isites, pairs)
#       bm |= make_bitmask(hs.bitoffsets[isite+1], hs.bitoffsets[isite]; dtype=BR)
#       bs |= BR(p[1]-1) << hs.bitoffsets[isite]
#       bt |= BR(p[2]-1) << hs.bitoffsets[isite]
#       am *= p[3]
#     end
#     push!(terms, PureOperator{S, BR}(hs, bm, bs, bt, am))
#   end

#   return SumOperator{S, BR}(hs, terms)

# end




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


# export apply, apply!, apply_unsafe!
#
#
# function apply_unsafe!(out::DenseState{S1, QN, BR}, nullop ::NullOperator, psi::DenseState{S2, QN, BR}) where {S1, S2, QN, BR}
#   return out
# end
#
# function apply_unsafe!(out::DenseState{S1, QN, BR}, psi::DenseState{S2, QN, BR}, nullop ::NullOperator) where {S1, S2, QN, BR}
#   return out
# end
#
# function apply_unsafe!(out::DenseState{S1, QN, BR}, pureop ::PureOperator{S2, BR}, psi::DenseState{S3, QN, BR}) where {S1, S2, S3, QN, BR}
#   for (b, v) in zip(psi.hilbert_space_realization.basis_list, psi.components)
#     if (b & pureop.bitmask) == pureop.bitcol
#       b2 = (b & ~pureop.bitmask) | pureop.bitrow
#       out[b2] += pureop.amplitude * v
#     end
#   end
#   out
# end
#
# function apply_unsafe!(out::DenseState{S1, QN, BR}, psi::DenseState{S3, QN, BR}, pureop ::PureOperator{S2, BR}) where {S1, S2, S3, QN, BR}
#   for (b, v) in zip(psi.hilbert_space_realization.basis_list, psi.components)
#     if (b & pureop.bitmask) == pureop.bitrow
#       b2 = (b & ~pureop.bitmask) | pureop.bitcol
#       out[b2] += v * pureop.amplitude
#     end
#   end
#   out
# end
#
# function apply_unsafe!(out::DenseState{S1, QN, BR}, sumop ::SumOperator{S2, BR}, psi::DenseState{S3, QN, BR}) where {S1, S2, S3, QN, BR}
#   for t in sumop.terms
#     apply_unsafe!(out, t, psi)
#   end
#   out
# end
#
# function apply_unsafe!(out::DenseState{S1, QN, BR}, psi::DenseState{S3, QN, BR}, sumop ::SumOperator{S2, BR}) where {S1, S2, S3, QN, BR}
#   for t in sumop.terms
#     apply_unsafe!(out, psi, t)
#   end
#   out
# end
#
#
#
# """
#     apply!
#
# Apply operator to `psi` and add it to `out`.
# """
# function apply!(out::DenseState{S1, QN, BR}, nullop ::NullOperator, psi::DenseState{S2, QN, BR}) where {S1, S2, QN, BR}
#   return out
# end
#
# function apply!(out::DenseState{S1, QN, BR}, psi::DenseState{S2, QN, BR}, nullop ::NullOperator) where {S1, S2, QN, BR}
#   return out
# end
#
# function apply!(out::DenseState{S1, QN, BR}, pureop ::PureOperator{S2, BR}, psi::DenseState{S3, QN, BR}) where {S1, S2, S3, QN, BR}
#   if out.hilbert_space_realization !== psi.hilbert_space_realization || pureop.hilbert_space !== psi.hilbert_space_realization.hilbert_space
#     throw(ArgumentError("Hilbert spaces of lhs and rhs of + should match"))
#   end
#   apply_unsafe!(out, pureop, psi)
# end
#
# function apply!(out::DenseState{S1, QN, BR}, psi::DenseState{S3, QN, BR}, pureop ::PureOperator{S2, BR}) where {S1, S2, S3, QN, BR}
#   if out.hilbert_space_realization !== psi.hilbert_space_realization || pureop.hilbert_space !== psi.hilbert_space_realization.hilbert_space
#     throw(ArgumentError("Hilbert spaces of lhs and rhs of + should match"))
#   end
#   apply_unsafe!(out, psi, pureop)
# end
#
#
# function apply!(out::DenseState{S1, QN, BR}, sumop ::SumOperator{S2, BR}, psi::DenseState{S3, QN, BR}) where {S1, S2, S3, QN, BR}
#   if out.hilbert_space_realization !== psi.hilbert_space_realization || sumop.hilbert_space !== psi.hilbert_space_realization.hilbert_space
#     throw(ArgumentError("Hilbert spaces should match"))
#   end
#   apply_unsafe!(out, sumop, psi)
# end
#
# function apply!(out::DenseState{S1, QN, BR}, psi::DenseState{S3, QN, BR}, sumop ::SumOperator{S2, BR}) where {S1, S2, S3, QN, BR}
#   if out.hilbert_space_realization !== psi.hilbert_space_realization || sumop.hilbert_space !== psi.hilbert_space_realization.hilbert_space
#     throw(ArgumentError("Hilbert spaces should match"))
#   end
#   apply_unsafe!(out, psi, sumop)
# end
#
#
# function apply(pureop ::NullOperator, psi::DenseState{S2, QN, BR}) where {S2, QN, BR}
#   return DenseState{S2, QN, BR}(psi.hilbert_space)
# end
#
# function apply(psi::DenseState{S2, QN, BR}, pureop ::NullOperator) where {S2, QN, BR}
#   return DenseState{S2, QN, BR}(psi.hilbert_space)
# end
#
#
# function apply(pureop ::PureOperator{S1, BR}, psi::DenseState{S2, QN, BR}) where {S1, S2, QN, BR}
#   if pureop.hilbert_space !== psi.hilbert_space_realiation.hilbert_space
#     throw(ArgumentError("Hilbert spaces of lhs and rhs of + should match"))
#   end
#   S3 = promote_type(S1, S2)
#   out = DenseState{S3, QN, BR}(psi.hilbert_space_realization)
#   apply_unsafe!(out, pureop, psi)
#   return out
# end
#
# function apply(psi::DenseState{S2, QN, BR}, pureop ::PureOperator{S1, BR}) where {S1, S2, QN, BR}
#   if pureop.hilbert_space !== psi.hilbert_space_realiation.hilbert_space
#     throw(ArgumentError("Hilbert spaces of lhs and rhs of + should match"))
#   end
#   S3 = promote_type(S1, S2)
#   out = DenseState{S3, QN, BR}(psi.hilbert_space_realization)
#   apply_unsafe!(out, psi, pureop)
#   return out
# end
#
#
# function apply(sumop ::SumOperator{S1, BR}, psi::DenseState{S2, QN, BR}) where {S1, S2, QN, BR}
#   if sumop.hilbert_space !== psi.hilbert_space_realiation.hilbert_space
#     throw(ArgumentError("Hilbert spaces of lhs and rhs of + should match"))
#   end
#   S3 = promote_type(S1, S2)
#   out = DenseState{S3, QN, BR}(psi.hilbert_space_realization)
#   for t in sumop.terms
#     apply_unsafe!(out, t, psi)
#   end
#   return out
# end
#
#
# function apply(psi::DenseState{S2, QN, BR}, sumop ::SumOperator{S1, BR}) where {S1, S2, QN, BR}
#   S3 = promote_type(S1, S2)
#   if sumop.hilbert_space !== psi.hilbert_space_realiation.hilbert_space
#     throw(ArgumentError("Hilbert spaces of lhs and rhs of + should match"))
#   end
#   out = DenseState{S3, QN, BR}(psi.hilbert_space_realization)
#   for t in sumop.terms
#     apply_unsafe!(out, psi, t)
#   end
#   return out
# end




# # TODO: Replace with more efficient
# function apply_naive(hs ::HilbertSpace,
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
# function apply_naive(hs ::HilbertSpace,
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
# function apply(hs ::HilbertSpace,
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

# function apply(hs ::HilbertSpace,
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

# function apply(hs ::HilbertSpace,
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

# function apply(hs ::HilbertSpace,
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
