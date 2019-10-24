export OperatorRepresentation
export represent

struct OperatorRepresentation{HSR <:HilbertSpaceRepresentation, O <:AbstractOperator} <: AbstractOperatorRepresentation
  hilbert_space_representation ::HSR
  operator ::O

  function OperatorRepresentation(hsr ::HSR, op ::O) where {HSR<:HilbertSpaceRepresentation, O<:AbstractOperator}
    new{HSR, O}(hsr, op)
  end
  function OperatorRepresentation{HSR, O}(hsr ::HSR, op ::O) where {HSR<:HilbertSpaceRepresentation, O<:AbstractOperator}
    new{HSR, O}(hsr, op)
  end
end

function represent(hsr ::HSR, op ::O) where {HSR<:HilbertSpaceRepresentation, O<:AbstractOperator}
  return OperatorRepresentation{HSR, O}(hsr, op)
end


@inline spacetype(lhs::Type{OperatorRepresentation{HSR, O}}) where {HSR, O}= HSR
@inline operatortype(lhs ::Type{OperatorRepresentation{HSR, O}}) where {HSR, O} = O
@inline get_space(lhs ::OperatorRepresentation{HSR, O}) where {HSR, O} = lhs.hilbert_space_representation


import Base.==
function (==)(opr1 ::HilbertSpaceRepresentation{HSR1, O1},
              opr2 ::HilbertSpaceRepresentation{HSR2, O2}) where {HSR1, O1, HSR2, O2}
  return ((opr1.hilbert_space_representation == opr2.hilbert_space_representation) &&
          (opr1.operator == opr2.operator))
end


## iterators

"""
May contain duplicates
"""
function get_row_iterator(opr ::OperatorRepresentation{HSR, O},
                          irow ::Integer;
                          include_all::Bool=false) where {HSR, O}
  hsr = opr.hilbert_space_representation
  brow = hsr.basis_list[irow]
  iter = (get(hsr.basis_lookup, bcol, -1) => amplitude
            for (bcol, amplitude) in get_row_iterator(opr.operator, brow))
  if include_all
    return iter
  else
    dim = size(opr, 2)
    return (icol => amplitude for (icol, amplitude) in iter if 1 <= icol <= dim)
  end
end


function get_column_iterator(opr ::OperatorRepresentation{HSR, O}, icol ::Integer; include_all::Bool=false) where {HSR, O}
  hsr = opr.hilbert_space_representation
  bcol = hsr.basis_list[icol]
  iter = (get(hsr.basis_lookup, brow, -1) => amplitude
            for (brow, amplitude) in get_column_iterator(opr.operator, bcol))
  if include_all
    return iter
  else
    dim = size(opr, 1)
    return (irow => amplitude for (irow, amplitude) in iter if 1 <= irow <= dim)
  end
end


function apply_unsafe!(out ::Vector{S1},
                       opr ::OperatorRepresentation{HSR, O},
                       state ::AbstractVector{S2};
                       range ::AbstractVector{<:Integer}=1:dimension(opr.hilbert_space_representation)
                       ) where {HSR, O, S1<:Number, S2<:Number}
  hsr = opr.hilbert_space_representation
  dim = length(hsr.basis_list)
  err = zero(Float64)
  for icol in range
    bcol = hsr.basis_list[icol]
    v = state[icol]
    for ((brow, _), amplitude) in get_column_iterator(opr.operator, bcol)
      irow = get(hsr.basis_lookup, brow, -1)
      if 1 <= irow <= dim
        out[irow] += amplitude * v
      else
        err += abs(amplitude*v)^2
      end
    end
  end
  sqrt(err)
end

function apply_unsafe!(out ::Vector{S1},
                       state ::AbstractVector{S2},
                       opr ::OperatorRepresentation{HSR, O};
                       range ::AbstractVector{<:Integer}=1:dimension(opr.hilbert_space_representation)
                       ) where {HSR, O, S1<:Number, S2<:Number}
  hsr = opr.hilbert_space_representation
  dim = length(hsr.basis_list)
  err = zero(Float64)
  for irow in range
    brow = hsr.basis_list[irow]
    v = state[irow]
    for (bcol, amplitude) in get_row_iterator(brow, opr.operator)
      icol = get(hsr.basis_lookup, bcol, -1)
      if 1 <= icol <= dim
        out[icol] += v * amplitude
      else
        err += abs(v * amplitude)^2
      end
    end
  end
  sqrt(err)
end


import Base.*
function (*)(opr ::OperatorRepresentation{HSR, O}, state ::AbstractVector{S}) where {HSR, O, S<:Number}
  hsr = opr.hilbert_space_representation
  n = length(hsr.basis_list)
  T = promote_type(S, eltype(O))
  out = zeros(T, n)
  err = apply_unsafe!(out, opr, state)
  return out
end


import Base.*
function (*)(state ::AbstractVector{S}, opr ::OperatorRepresentation{HSR, O}) where {HSR, O, S<:Number}
  hsr = opr.hilbert_space_representation
  n = len(hsr.basis_list)
  T = promote_type(S, eltype(O))
  out = zeros(T, n)
  err = apply_unsafe!(out, state, opr)
  out
end



# function splitblock(n ::Integer, b ::Integer) ::Vector{Int}
#   (n < 0) && throw(ArgumentError("n cannot be negative"))
#   (b <= 0) && throw(ArgumentError("b must be positive"))
#   #b >= n && return    # too many blocks
#   blocksize = n ÷ b
#   blocks = blocksize * ones(Int, b)
#   r = n - (blocksize * b)
#   blocks[1:r] .+= 1
#   return blocks
# end
#
# function apply_thread_unsafe!(out ::Vector{S1},
#                               opr ::OperatorRepresentation{HS, O},
#                               state ::AbstractVector{S2}
#                               ) where {HS, O, S1<:Number, S2<:Number}
#
#   nthreads = Threads.nthreads()
#   nthreads == 1 && return apply_unsafe!(out, opr, stage, range)
#
#   hsr = opr.hilbert_space_representation
#   dim = length(hsr.basis_list)
#   counts = splitblock(dim, nthreads)
#   offsets = cumsum([1, counts...])
#   locks = [Threads.SpinLock() for i in 1:nthreads]
#   offdiagonals = [Dict(Int, S1) for i in 1:nthreads]
#
#   err = zero(Float64)
#   Threads.@threads for b in 1:nthreads
#     lower = offsets[b]
#     upper = offsets[b+1]-1
#     for icol in lower:upper
#       bcol = hsr.basis_list[icol]
#       v = state[icol]
#       for (brow, amplitude) in get_column_iterator(opr.operator, bcol)
#         irow = get(hsr.basis_lookup, brow, -1)
#
#         if lower <= irow <= upper
#           out[irow] += amplitude * v
#         elseif 1 <= irow <= dim
#           off diagonal
#         else
#           err += abs(amplitude*v)^2
#         end
#       end # for brow
#     end # for icol
#   end # for b
#   sqrt(err)
# end
