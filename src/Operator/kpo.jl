# export KroneckerProductOperator
# export clean!

# struct KroneckerProductOperator{Scalar<:Number} <:AbstractOperator
#   hilbert_space ::AbstractHilbertSpace
#   amplitude ::Scalar
#   operators ::Dict{Int, Matrix{Scalar}}

#   function KroneckerProductOperator(
#       hs ::AbstractHilbertSpace,
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
#     hs ::AbstractHilbertSpace,
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
