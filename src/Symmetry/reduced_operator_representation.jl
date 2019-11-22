export ReducedOperatorRepresentation
export bintype
export represent
export get_row_iterator
export get_column_iterator
export get_row
export get_column

struct ReducedOperatorRepresentation{RHSR <:ReducedHilbertSpaceRepresentation, O <:AbstractOperator, S<:Number, BR<:Unsigned} <:AbstractOperatorRepresentation{S}
  reduced_hilbert_space_representation ::RHSR
  operator ::O
  function ReducedOperatorRepresentation(rhsr ::RHSR, op::O) where {RHSR <:ReducedHilbertSpaceRepresentation, O <:AbstractOperator}
    S = promote_type(scalartype(RHSR), scalartype(O))
    BR = bintype(RHSR)
    return new{RHSR, O, S, BR}(rhsr, op)
  end
end

spacetype(lhs::Type{ReducedOperatorRepresentation{RHSR, O, S, BR}}) where {RHSR, O, S, BR} = RHSR
operatortype(lhs ::Type{ReducedOperatorRepresentation{RHSR, O, S, BR}}) where {RHSR, O, S, BR} = O
get_space(lhs ::ReducedOperatorRepresentation{RHSR, O, S, BR}) where {RHSR, O, S, BR} = lhs.reduced_hilbert_space_representation ::RHSR

function represent(rhsr ::RHSR, op ::O) where {RHSR <:ReducedHilbertSpaceRepresentation, O <:AbstractOperator}
  return ReducedOperatorRepresentation(rhsr, op)
end

import Base.show
function show(io::IO, ::MIME"text/plain", arg::ReducedOperatorRepresentation{RHSR, O, S, BR})  where {RHSR, O, S, BR}
  print(io, string(typeof(arg)), "(", arg.reduced_hilbert_space_representation, ", ", arg.operator, ")")
end

function get_row_iterator(opr ::ReducedOperatorRepresentation{RHSR, O, S, BR},
                          irow_r ::Integer) where {RHSR, O, S, BR}
  rhsr = opr.reduced_hilbert_space_representation
  hsr = rhsr.parent
  dim = dimension(opr.reduced_hilbert_space_representation)

  brow = rhsr.basis_list[irow_r]
  irow_p = hsr.basis_lookup[brow]
  irow_r2, ampl_row = rhsr.basis_mapping[irow_p]
  inv_ampl_row = one(S) / ampl_row

  zero_val = zero(S)

  full_iter = let
    basis_lookup = rhsr.parent.basis_lookup
    basis_mapping ::Vector{NamedTuple{(:index, :amplitude), Tuple{Int, S}}} = rhsr.basis_mapping
    operator = opr.operator
    function element(bcol::BR, ampl::S) ::Pair{Int, S}
      icol_p = get(basis_lookup, bcol, -1)
      (icol_p > 0) || return (-1 => ampl)
      (icol_r, ampl_col) = basis_mapping[icol_p]
      (icol_r > 0) || return (-1 => ampl)
      return icol_r => ampl * ampl_col * inv_ampl_row
    end
    (element(bcol, ampl) for (bcol::BR, ampl::S) in get_row_iterator(operator, brow))
  end
  return full_iter
end

function get_column_iterator(opr ::ReducedOperatorRepresentation{RHSR, O, S, BR},
                             icol_r ::Integer) where {RHSR, O, S, BR}
  rhsr = opr.reduced_hilbert_space_representation
  hsr = rhsr.parent
  dim = dimension(opr.reduced_hilbert_space_representation)
  bcol = rhsr.basis_list[icol_r]
  icol_p = hsr.basis_lookup[bcol]
  icol_r2, ampl_col = rhsr.basis_mapping[icol_p]
  #@assert icol_r == icol_r2 "$icol_r != $icol_r2"

  inv_ampl_col = one(S) / ampl_col
  zero_val = zero(S)

  full_iter = let
    basis_lookup = rhsr.parent.basis_lookup
    basis_mapping ::Vector{NamedTuple{(:index, :amplitude), Tuple{Int, S}}} = rhsr.basis_mapping
    operator = opr.operator
    function element(brow::BR, ampl::S) ::Pair{Int, S}
      irow_p = get(basis_lookup, brow, -1)
      (irow_p > 0) || return (-1 => ampl)
      (irow_r, ampl_row) = basis_mapping[irow_p]
      (irow_r > 0) || return (-1 => ampl)
      return irow_r => (ampl * ampl_row) * inv_ampl_col
    end
    (element(brow, ampl) for (brow::BR, ampl::S) in get_column_iterator(operator, bcol))
  end
  return full_iter
end

# TODO: better implementation
function get_element(opr ::ReducedOperatorRepresentation{RHSR, O, S, BR},
                     irow_r ::Integer, icol_r ::Integer) where {RHSR, O, S, BR}
  rhsr = opr.reduced_hilbert_space_representation
  @boundscheck let
   dim = length(rhsr.basis_list)
   if irow_r <= 0 || irow_r > dim || icol_r <= 0 || icol_r > dim
     throw(BoundsError(opr, [irow_r, icol_r]))
   end
  end
  return @inbounds sum(ampl for (irow_r2, ampl::S) in get_column_iterator(opr, icol_r) if irow_r2 == irow_r)
end
