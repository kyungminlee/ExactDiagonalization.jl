export apply
export materialize, materialize_parallel

# """
# Returns dict
# """
# function apply(hs ::AbstractHilbertSpace,
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


# function apply(hs ::AbstractHilbertSpace,
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


# function apply(hs ::AbstractHilbertSpace,
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

import Base.isempty
isempty(psi::SparseState{S2, BR}) where {S2, BR} = isempty(psi.components)

function apply(pureop ::NullOperator, psi::SparseState{S2, BR}) where {S2, BR}
  return SparseState{S2, BR}(psi.hilbert_space)
end


function apply(pureop ::PureOperator{S1, BR}, psi::SparseState{S2, BR}) where {S1, S2, BR}
  S3 = promote_type(S1, S2)
  if pureop.hilbert_space != psi.hilbert_space
    throw(ArgumentError("Hilbert spaces of lhs and rhs of + should match"))
  end
  out = SparseState{S3, BR}(psi.hilbert_space)
  for (b, v) in psi.components
    if (b & pureop.bitmask) == pureop.bitsource
      b2 = (b & ~pureop.bitmask) | pureop.bittarget
      out[b2] += pureop.amplitude * v
    end
  end
  return out
end

function apply(sumop ::SumOperator{S1, BR}, psi::SparseState{S2, BR}) where {S1, S2, BR}
  S3 = promote_type(S1, S2)
  if sumop.hilbert_space != psi.hilbert_space
    throw(ArgumentError("Hilbert spaces of lhs and rhs of + should match"))
  end
  out = SparseState{S3, BR}(psi.hilbert_space)
  for t in sumop.terms
    out += apply(t, psi)
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


# function materialize(
#         hsb ::ConcreteHilbertSpaceBlock{BinRep, QN},
#         ops ::AbstractArray{KroneckerProductOperator{OS}};
#         Scalar::DataType=OS) where {BinRep, QN, OS<:Number}
#   hs = hsb.hilbert_space
#   rows = Int[]
#   cols = Int[]
#   vals = Scalar[]
#   err = 0.0
#   for (irow, row) in enumerate(hsb.basis_list)
#     out = apply(hs, ops, row)
#     for (col, amplitude) in out
#       if ! isapprox(amplitude, 0)
#         if !haskey(hsb.basis_lookup, col)
#           err += abs(amplitude)^2
#         else 
#           icol = hsb.basis_lookup[col]
#           push!(rows, irow)
#           push!(cols, icol)
#           push!(vals, amplitude)
#         end        
#       end
#     end
#   end
  
#   if Scalar <:Complex
#     if isapprox( norm(imag.(vals)), 0)
#       vals = real.(vals)
#     end
#   end
#   n = dimension(hsb)
#   return (sparse(rows, cols, vals, n, n), err)
# end



# function materialize(
#         hsb ::ConcreteHilbertSpace{QN, BR},
#         ops ::AbstractArray{KroneckerProductOperator{OS}};
#         Scalar::DataType=OS) where {QN, BR<:Unsigned, OS<:Number}
#   hs = hsb.hilbert_space
#   rows = Int[]
#   cols = Int[]
#   vals = Scalar[]
#   err = 0.0
#   for (irow, row) in enumerate(hsb.basis_list)
#     out = apply(hs, ops, row)
#     for (col, amplitude) in out
#       if ! isapprox(amplitude, 0)
#         if !haskey(hsb.basis_lookup, col)
#           err += abs(amplitude)^2
#         else 
#           icol = hsb.basis_lookup[col]
#           push!(rows, irow)
#           push!(cols, icol)
#           push!(vals, amplitude)
#         end
#       end
#     end
#   end
  
#   if Scalar <:Complex
#     if isapprox( maximum(abs.(imag.(vals))), 0)
#       vals = real.(vals)
#     end
#   end
#   n = dimension(hsb)
#   return (sparse(rows, cols, vals, n, n), err)
# end



# function materialize_parallel(
#         hsb ::ConcreteHilbertSpace{QN, BR},
#         ops ::AbstractArray{KroneckerProductOperator{OS}};
#         Scalar::DataType=OS) where {QN, BR<:Unsigned, OS<:Number}
#   hs = hsb.hilbert_space
#   err = 0.0

#   nthreads = Threads.nthreads()
#   local_rows = [ Int[] for i in 1:nthreads]
#   local_cols = [ Int[] for i in 1:nthreads]
#   local_vals = [ Scalar[] for i in 1:nthreads]
#   local_err = [0.0 for i in 1:nthreads]

#   n_basis = length(hsb.basis_list)

#   Threads.@threads for irow in 1:n_basis
#     id = Threads.threadid()
#     row = hsb.basis_list[irow]
#     out = apply(hs, ops, row)
#     for (col, amplitude) in out
#       if ! isapprox(amplitude, 0)
#         if !haskey(hsb.basis_lookup, col)
#           local_err[id] += abs(amplitude)^2
#         else 
#           icol = hsb.basis_lookup[col]
#           push!(local_rows[id], irow)
#           push!(local_cols[id], icol)
#           push!(local_vals[id], amplitude)
#         end
#       end
#     end
#   end
#   rows = vcat(local_rows...)
#   cols = vcat(local_cols...)
#   vals = vcat(local_vals...)
#   err = sum(local_err)
  
#   if Scalar <:Complex
#     if isapprox( norm(imag.(vals)), 0)
#       vals = real.(vals)
#     end
#   end
#   n = dimension(hsb)
#   return (sparse(rows, cols, vals, n, n), err)
# end


function materialize(
  chs ::ConcreteHilbertSpace{QN, BR},
  sumop ::SumOperator{S, BR}) where {QN, BR<:Unsigned, S<:Number}
  hs = chs.hilbert_space
  rows = Int[]
  cols = Int[]
  vals = S[]
  
  err = 0.0
  
  for (irow, row) in enumerate(chs.basis_list)
    ψrow = SparseState{S, BR}(hs, row)
    ψcol = apply(sumop, ψrow)
    for (col, amplitude) in ψcol.components
      if ! isapprox(amplitude, 0)
        if !haskey(chs.basis_lookup, col)
          err += abs(amplitude)^2
        else 
          icol = chs.basis_lookup[col]
          push!(rows, irow)
          push!(cols, icol)
          push!(vals, amplitude)
        end
      end
    end
  end

  if S <:Complex
    if isapprox( maximum(abs.(imag.(vals))), 0)
      vals = real.(vals)
    end
  end
  n = dimension(chs)
  return (sparse(rows, cols, vals, n, n), err)
end




function materialize_parallel(
  chs ::ConcreteHilbertSpace{QN, BR},
  sumop ::SumOperator{S, BR}) where {QN, BR<:Unsigned, S<:Number}

  hs = chs.hilbert_space

  nthreads = Threads.nthreads()
  local_rows = [ Int[] for i in 1:nthreads]
  local_cols = [ Int[] for i in 1:nthreads]
  local_vals = [ S[] for i in 1:nthreads]
  local_err = [0.0 for i in 1:nthreads]

  n_basis = length(chs.basis_list)

  Threads.@threads for irow in 1:n_basis
    id = Threads.threadid()
    row = chs.basis_list[irow]

    ψrow = SparseState{S, BR}(hs, row)
    ψcol = apply(sumop, ψrow)
    for (col, amplitude) in ψcol.components
      if ! isapprox(amplitude, 0)
        if !haskey(chs.basis_lookup, col)
          local_err[id] += abs(amplitude)^2
        else 
          icol = chs.basis_lookup[col]
          push!(local_rows[id], irow)
          push!(local_cols[id], icol)
          push!(local_vals[id], amplitude)
        end
      end
    end
  end
  
  rows = vcat(local_rows...)
  cols = vcat(local_cols...)
  vals = vcat(local_vals...)
  err = sum(local_err)

  if isempty(vals)
    vals = Float64[]
  elseif S <:Complex && isapprox(maximum(abs.(imag.(vals))), 0)
    vals = real.(vals)
  end
  n = dimension(chs)
  return (sparse(rows, cols, vals, n, n), err)
end
