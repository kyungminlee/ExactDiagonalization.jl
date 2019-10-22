export represent
export get_slice
export ReducedOperatorRepresentation

struct ReducedOperatorRepresentation{RH <:ReducedHilbertSpaceRealization, O <:AbstractOperator}
    reduced_hilbert_space_realization ::RH
    operator ::O
end

bintype(lhs ::ReducedOperatorRepresentation{RH, O}) where {RH, O} = bintype(RH)
bintype(lhs ::Type{ReducedOperatorRepresentation{RH, O}}) where {RH, O} = bintype(RH)

import Base.eltype
eltype(lhs ::ReducedOperatorRepresentation{RH, O}) where {RH, O} = promote_type(eltype(RH), eltype(O))
eltype(lhs ::Type{ReducedOperatorRepresentation{RH, O}}) where {RH, O} = promote_type(eltype(RH), eltype(O))


import Base.size
function size(arg::ReducedOperatorRepresentation{RH, O}) ::Tuple{Int, Int} where {RH, O}
  dim = dimension(arg.reduced_hilbert_space_realization)
  return (dim, dim)
end

function represent(rhsr ::RH, op ::O) where {RH <:ReducedHilbertSpaceRealization, O <:AbstractOperator}
    return ReducedOperatorRepresentation{RH, O}(rhsr, op)
end

function get_slice(irow_r ::Integer, opr ::ReducedOperatorRepresentation{RH, O}) where {RH, O}
  rhsr = opr.reduced_hilbert_space_realization
  hsr = rhsr.parent
  S = eltype(opr.operator)

  brow = rhsr.basis_list[irow_r]
  irow_p = hsr.basis_lookup[brow]
  irow_r2, ampl_row = rhsr.basis_mapping[irow_p]
  @assert irow_r == irow_r2 "$irow_r != $irow_r2"

  zero_val = 0.0 #zero(S)
  missing_val = (-1 => zero_val)

  function f(bcol, ampl)
    haskey(rhsr.parent.basis_lookup, bcol) || return missing_val
    icol_p = rhsr.parent.basis_lookup[bcol]
    (icol_r, ampl_col) = rhsr.basis_mapping[icol_p]
    (icol_r > 0) || return missing_val
    return icol_r => (ampl / ampl_row) * ampl_col
  end
  full_iter = (f(bcol, ampl) for (bcol, ampl) in get_slice(brow, opr.operator))

  return ( (icol_r, ampl) for (icol_r, ampl) in full_iter if icol_r >= 0)
end


function get_slice(opr ::ReducedOperatorRepresentation{RH, O}, icol_r ::Integer) where {RH, O}
  rhsr = opr.reduced_hilbert_space_realization
  hsr = rhsr.parent
  S = eltype(opr.operator)

  bcol = rhsr.basis_list[icol_r]
  icol_p = hsr.basis_lookup[bcol]
  icol_r2, ampl_col = rhsr.basis_mapping[icol_p]
  @assert icol_r == icol_r2 "$icol_r != $icol_r2"

  zero_val = 0.0 #zero(S)
  missing_val = (-1 => zero_val)

  function f(brow, ampl)
    haskey(rhsr.parent.basis_lookup, brow) || return missing_val
    irow_p = rhsr.parent.basis_lookup[brow]
    (irow_r, ampl_row) = rhsr.basis_mapping[irow_p]
    (irow_r > 0) || return missing_val
    return irow_r => (ampl / ampl_col) * ampl_row
  end
  full_iter = (f(brow, ampl) for (brow, ampl) in get_slice(opr.operator, bcol))

  return ( (irow_r, ampl) for (irow_r, ampl) in full_iter if irow_r >= 0)
end


import SparseArrays.sparse
function sparse(opr::ReducedOperatorRepresentation{H, O}; tol ::Real=sqrt(eps(Float64))) where {H, O}
  S = promote_type(eltype(opr), eltype(H))
  m, n = size(opr)
  colptr = zeros(Int, n+1)
  rowval = Int[]
  nzval = S[]

  colptr[1] = 1
  for icol in 1:n

    colvec = Dict{Int, S}()
    for (irow, ampl) in get_slice(opr, icol)
      colvec[irow] = get(colvec, irow, zero(S)) + ampl
    end
    to_delete = Int[irow for (irow, ampl) in colvec if abs(ampl) < tol]
    for irow in to_delete
      delete!(colvec, irow)
    end

    colptr[icol+1] = colptr[icol] + length(colvec)
    sorted_items = sort(collect(colvec), by = item -> item[1])
    append!(rowval, irow for (irow, ampl) in sorted_items)
    append!(nzval, ampl for (irow, ampl) in sorted_items)
  end
  return SparseMatrixCSC{S, Int}(m, n, colptr, rowval, nzval)
end
