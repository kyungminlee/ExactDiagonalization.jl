export materialize, materialize_parallel


function materialize(
  hsr ::HilbertSpaceRepresentation{QN, BR},
  nop ::NullOperator) where {QN, BR<:Unsigned}
  n = dimension(hsr)
  return (sparse([], [], Float64[], n, n), 0.0)
end


function materialize_parallel(
    hsr ::HilbertSpaceRepresentation{QN, BR},
    nop ::NullOperator) where {QN, S<:Number, BR<:Unsigned}
  n = dimension(hsr)
  return (sparse([], [], Float64[], n, n), 0.0)
end


function materialize(
    hsr ::HilbertSpaceRepresentation{QN, BR},
    sumop ::SumOperator{S, BR};
    tol::Real=sqrt(eps(Float64))) where {QN, S<:Number, BR<:Unsigned}
  debug(LOGGER, "BEGIN materialize for HilbertSpaceRepresentation")
  hs = hsr.hilbert_space
  rows = Int[]
  cols = Int[]
  vals = S[]
  err = 0.0

  for (irow, row) in enumerate(hsr.basis_list)
    ψrow = SparseState{S, BR}(hs, row)
    ψcol = SparseState{S, BR}(hs)
    apply_unsafe!(ψcol, ψrow, sumop)
    for (col, amplitude) in ψcol.components
      isapprox(amplitude, 0) && continue

      if !haskey(hsr.basis_lookup, col)
        err += abs(amplitude)^2
      else
        icol = hsr.basis_lookup[col]
        push!(rows, irow)
        push!(cols, icol)
        push!(vals, amplitude)
      end
    end
  end

  if isempty(vals)
    debug(LOGGER, "Matrix empty")
    vals = Float64[]
  elseif S <:Complex && isapprox( maximum(abs.(imag.(vals))), 0)
    debug(LOGGER, "Matrix purely real")
    vals = real.(vals)
  end
  n = dimension(hsr)
  spmat = sparse(rows, cols, vals, n, n)
  debug(LOGGER, "Number of nonzero elements: $(length(spmat.nzval))")
  SparseArrays.droptol!(spmat, tol)
  debug(LOGGER, "Number of nonzero elements after choptol!: $(length(spmat.nzval))")
  debug(LOGGER, "END materialize for HilbertSpaceRepresentation")
  return (sparse(rows, cols, vals, n, n), err)
end


function materialize_parallel(
    hsr ::HilbertSpaceRepresentation{QN, BR},
    sumop ::SumOperator{S, BR};
    tol::Real=sqrt(eps(Float64))) where {QN, S<:Number, BR<:Unsigned}
  debug(LOGGER, "BEGIN materialize_parallel for HilbertSpaceRepresentation")
  hs = hsr.hilbert_space

  nthreads = Threads.nthreads()
  local_rows = [ Int[] for i in 1:nthreads]
  local_cols = [ Int[] for i in 1:nthreads]
  local_vals = [ S[] for i in 1:nthreads]
  local_err =  Float64[0.0 for i in 1:nthreads]

  n_basis = length(hsr.basis_list)

  Threads.@threads for irow in 1:n_basis
    id = Threads.threadid()
    row = hsr.basis_list[irow]

    ψrow = SparseState{S, BR}(hs, row)
    ψcol = SparseState{S, BR}(hs)
    apply_unsafe!(ψcol, ψrow, sumop)
    for (col, amplitude) in ψcol.components
      isapprox(amplitude, 0) && continue

      if !haskey(hsr.basis_lookup, col)
        local_err[id] += abs(amplitude)^2
      else
        icol = hsr.basis_lookup[col]
        push!(local_rows[id], irow)
        push!(local_cols[id], icol)
        push!(local_vals[id], amplitude)
      end
    end
  end

  rows ::Vector{Int} = vcat(local_rows...)
  cols ::Vector{Int} = vcat(local_cols...)
  vals ::Vector{S} = vcat(local_vals...)
  err ::Float64 = sum(local_err)

  if isempty(vals)
    debug(LOGGER, "Matrix empty")
    vals = Float64[]
  elseif S <:Complex && isapprox( maximum(abs.(imag.(vals))), 0)
    debug(LOGGER, "Matrix purely real")
    vals = real.(vals)
  end
  n = dimension(hsr)
  spmat = sparse(rows, cols, vals, n, n)
  debug(LOGGER, "Number of nonzero elements: $(length(spmat.nzval))")
  droptol!(spmat, tol)
  debug(LOGGER, "Number of nonzero elements after choptol!: $(length(spmat.nzval))")
  debug(LOGGER, "END materialize_parallel for HilbertSpaceRepresentation")
  return (spmat, err)
end
