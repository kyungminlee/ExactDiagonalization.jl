export represent
export get_slice
export ReducedOperatorRepresentation
export bintype

struct ReducedOperatorRepresentation{RHSR <:ReducedHilbertSpaceRepresentation, O <:AbstractOperator}
  reduced_hilbert_space_realization ::RHSR
  operator ::O
end

bintype(lhs ::ReducedOperatorRepresentation{RHSR, O}) where {RHSR, O} = bintype(RHSR)
bintype(lhs ::Type{ReducedOperatorRepresentation{RHSR, O}}) where {RHSR, O} = bintype(RHSR)

import Base.eltype
@inline eltype(lhs ::ReducedOperatorRepresentation{RHSR, O}) where {RHSR, O} = promote_type(eltype(RHSR), eltype(O))
@inline eltype(lhs ::Type{ReducedOperatorRepresentation{RHSR, O}}) where {RHSR, O} = promote_type(eltype(RHSR), eltype(O))

import Base.size
function size(arg::ReducedOperatorRepresentation{RHSR, O}) ::Tuple{Int, Int} where {RHSR, O}
  dim = dimension(arg.reduced_hilbert_space_realization)
  return (dim, dim)
end

function represent(rhsr ::RHSR, op ::O) where {RHSR <:ReducedHilbertSpaceRepresentation, O <:AbstractOperator}
    return ReducedOperatorRepresentation{RHSR, O}(rhsr, op)
end


function get_slice(irow_r ::Integer, opr ::ReducedOperatorRepresentation{RHSR, O}) where {RHSR, O}
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


function get_slice(opr ::ReducedOperatorRepresentation{RHSR, O}, icol_r ::Integer) where {RHSR, O}
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
function sparse(opr::ReducedOperatorRepresentation{RHSR, O}; tol ::Real=sqrt(eps(Float64))) where {RHSR, O}
  S = promote_type(eltype(opr), eltype(RHSR))
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
    choptol!(colvec, tol)

    colptr[icol+1] = colptr[icol] + length(colvec)
    sorted_items = sort(collect(colvec), by = item -> item[1])
    append!(rowval, irow for (irow, ampl) in sorted_items)
    append!(nzval, ampl for (irow, ampl) in sorted_items)
  end
  return SparseMatrixCSC{S, Int}(m, n, colptr, rowval, nzval)
end
