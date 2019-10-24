export ReducedOperatorRepresentation
export bintype
export represent
export get_row_iterator
export get_column_iterator
export get_row
export get_column

struct ReducedOperatorRepresentation{RHSR <:ReducedHilbertSpaceRepresentation, O <:AbstractOperator} <:AbstractOperatorRepresentation
  reduced_hilbert_space_realization ::RHSR
  operator ::O
end

@inline spacetype(lhs::Type{ReducedOperatorRepresentation{RHSR, O}}) where {RHSR, O} = RHSR
@inline operatortype(lhs ::Type{ReducedOperatorRepresentation{RHSR, O}}) where {RHSR, O} = O
@inline get_space(lhs ::ReducedOperatorRepresentation{RHSR, O}) where {RHSR, O} = lhs.reduced_hilbert_space_realization

function represent(rhsr ::RHSR, op ::O) where {RHSR <:ReducedHilbertSpaceRepresentation, O <:AbstractOperator}
    return ReducedOperatorRepresentation{RHSR, O}(rhsr, op)
end


function get_row_iterator(opr ::ReducedOperatorRepresentation{RHSR, O},
                          irow_r ::Integer;
                          include_all::Bool=false) where {RHSR, O}
  rhsr = opr.reduced_hilbert_space_realization
  hsr = rhsr.parent
  S = eltype(opr.operator)
  dim = dimension(opr.reduced_hilbert_space_realization)

  brow = rhsr.basis_list[irow_r]
  irow_p = hsr.basis_lookup[brow]
  irow_r2, ampl_row = rhsr.basis_mapping[irow_p]
  @assert irow_r == irow_r2 "$irow_r != $irow_r2"

  zero_val = 0.0 #zero(S)
  missing_val = (-1 => zero_val)

  function element(bcol, ampl)
    haskey(rhsr.parent.basis_lookup, bcol) || return missing_val
    icol_p = rhsr.parent.basis_lookup[bcol]
    (icol_r, ampl_col) = rhsr.basis_mapping[icol_p]
    (icol_r > 0) || return missing_val
    return icol_r => (ampl / ampl_row) * ampl_col
  end
  full_iter = (element(bcol, ampl) for (bcol, ampl) in get_row_iterator(opr.operator, brow))

  if include_all
    return full_iter
  else
    return (icol_r => ampl for (icol_r, ampl) in full_iter if 1 <= icol_r <= dim)
  end
end


function get_column_iterator(opr ::ReducedOperatorRepresentation{RHSR, O}, icol_r ::Integer;
                             include_all::Bool=false) where {RHSR, O}
  rhsr = opr.reduced_hilbert_space_realization
  hsr = rhsr.parent
  S = eltype(opr.operator)
  dim = dimension(opr.reduced_hilbert_space_realization)

  bcol = rhsr.basis_list[icol_r]
  icol_p = hsr.basis_lookup[bcol]
  icol_r2, ampl_col = rhsr.basis_mapping[icol_p]
  @assert icol_r == icol_r2 "$icol_r != $icol_r2"

  zero_val = 0.0 #zero(S)
  missing_val = (-1 => zero_val)

  function element(brow, ampl)
    haskey(rhsr.parent.basis_lookup, brow) || return missing_val
    irow_p = rhsr.parent.basis_lookup[brow]
    (irow_r, ampl_row) = rhsr.basis_mapping[irow_p]
    (irow_r > 0) || return missing_val
    return irow_r => (ampl / ampl_col) * ampl_row
  end
  full_iter = (element(brow, ampl) for (brow, ampl) in get_column_iterator(opr.operator, bcol))

  if include_all
    return full_iter
  else
    return (irow_r => ampl for (irow_r, ampl) in full_iter if 1 <= irow_r <= dim)
  end
end
