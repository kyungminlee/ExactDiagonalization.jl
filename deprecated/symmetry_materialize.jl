
#
# function materialize(
#     rhsr ::ReducedHilbertSpaceRepresentation{QN, BR, C},
#     operator ::AbstractOperator;
#     tol::Real=sqrt(eps(Float64))) where {QN, BR, C}
#
#   @assert is_invariant(rhsr.translation_group, operator)
#   hs = rhsr.parent.hilbert_space
#
#   rows = Int[]
#   cols = Int[]
#   vals = C[]
#   err  = 0.0
#
#   n_basis = length(rhsr.basis_list)
#   for irow_r in 1:n_basis
#     brow = rhsr.basis_list[irow_r]
#     irow_p = rhsr.parent.basis_lookup[brow]
#     irow_r2, ampl_row = rhsr.basis_mapping[irow_p]
#     @assert irow_r == irow_r2 "$irow_r != $irow_r2"
#     ψrow = SparseState{C, BR}(hs, brow=>1/ampl_row)
#     ψcol = SparseState{C, BR}(hs)
#     apply_unsafe!(ψcol, ψrow, operator)
#     choptol!(ψcol, tol)
#
#     for (bcol, ampl) in ψcol.components
#       if ! haskey(rhsr.parent.basis_lookup, bcol)
#         err += abs(ampl)^2
#         continue
#       end
#
#       icol_p = rhsr.parent.basis_lookup[bcol]
#       (icol_r, ampl_col) = rhsr.basis_mapping[icol_p]
#       if !(icol_r > 0)
#         err += abs(ampl)^2
#         continue
#       end
#       push!(rows, irow_r)
#       push!(cols, icol_r)
#       push!(vals, ampl * ampl_col)
#     end
#   end
#
#   if isempty(vals)
#     vals = Float64[]
#   elseif isapprox(maximum(abs.(imag.(vals))), 0; atol=tol)
#     vals = real.(vals)
#   end
#   spmat = sparse(rows, cols, vals, n_basis, n_basis)
#   droptol!(spmat, tol)
#   return (spmat, err)
# end
#
#
#
# function materialize_parallel(
#     rhsr :: ReducedHilbertSpaceRepresentation{QN, BR, C},
#     operator ::AbstractOperator;
#     tol::Real=sqrt(eps(Float64))) where {QN, BR, C}
#   debug(LOGGER, "BEGIN materialize_parallel for ReducedHilbertSpaceRealiation")
#   debug(LOGGER, "Checking whether the operator is translationally invariant")
#   @assert is_invariant(rhsr.translation_group, operator)
#   hs = rhsr.parent.hilbert_space
#
#   nthreads = Threads.nthreads()
#   debug(LOGGER, "Number of threads: $nthreads")
#   local_rows = [ Int[] for i in 1:nthreads]
#   local_cols = [ Int[] for i in 1:nthreads]
#   local_vals = [ C[] for i in 1:nthreads]
#   local_err =  Float64[0.0 for i in 1:nthreads]
#
#   n_basis = length(rhsr.basis_list)
#   debug(LOGGER, "Starting materialization (parallel)")
#   Threads.@threads for irow_r in 1:n_basis
#     id = Threads.threadid()
#     brow = rhsr.basis_list[irow_r]
#     irow_p = rhsr.parent.basis_lookup[brow]
#     irow_r2, ampl_row = rhsr.basis_mapping[irow_p]
#     @assert irow_r == irow_r2 "$irow_r != $irow_r2"
#     ψrow = SparseState{C, BR}(hs, brow=>1/ampl_row)
#     ψcol = SparseState{C, BR}(hs)
#     apply_unsafe!(ψcol, ψrow, operator)
#     choptol!(ψcol, tol)
#
#     for (bcol, ampl) in ψcol.components
#       if ! haskey(rhsr.parent.basis_lookup, bcol)
#         local_err[id] += abs(ampl^2)
#         continue
#       end
#       icol_p = rhsr.parent.basis_lookup[bcol]
#       (icol_r, ampl_col) = rhsr.basis_mapping[icol_p]
#       if !(icol_r > 0)
#         local_err[id] += abs(ampl)^2
#         continue
#       end
#       push!(local_rows[id], irow_r)
#       push!(local_cols[id], icol_r)
#       push!(local_vals[id], ampl * ampl_col)
#     end
#   end
#   debug(LOGGER, "Finished materialization (parallel)")
#
#   rows ::Vector{Int} = vcat(local_rows...)
#   cols ::Vector{Int} = vcat(local_cols...)
#   vals ::Vector{C} = vcat(local_vals...)
#   err ::Float64 = sum(local_err)
#
#   if isempty(vals)
#     debug(LOGGER, "Matrix empty")
#     vals = Float64[]
#   elseif isapprox(maximum(abs.(imag.(vals))), 0; atol=tol)
#     debug(LOGGER, "Matrix purely real")
#     vals = real.(vals)
#   end
#   spmat = sparse(rows, cols, vals, n_basis, n_basis)
#   debug(LOGGER, "Number of nonzero elements: $(length(spmat.nzval))")
#   choptol!(spmat, tol)
#   debug(LOGGER, "Number of nonzero elements after choptol!: $(length(spmat.nzval))")
#   debug(LOGGER, "END materialize_parallel for ReducedHilbertSpaceRealiation")
#   return (spmat, err)
# end
