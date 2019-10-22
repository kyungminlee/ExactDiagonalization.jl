
export AbstractOperatorRepresentation
export OperatorRepresentation
export represent
export bintype

abstract type AbstractOperatorRepresentation end

struct OperatorRepresentation{HSR <:HilbertSpaceRepresentation, O <:AbstractOperator} <: AbstractOperatorRepresentation
  hilbert_space_realization ::HSR
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


import Base.size
function size(arg::OperatorRepresentation{HSR, O}) ::Tuple{Int, Int} where {HSR, O}
  dim = dimension(arg.hilbert_space_realization)
  return (dim, dim)
end

@inline bintype(lhs ::OperatorRepresentation{HSR, O}) where {HSR, O} = bintype(HSR)
@inline bintype(lhs ::Type{OperatorRepresentation{HSR, O}}) where {HSR, O} = bintype(HSR)

import Base.eltype
@inline eltype(lhs ::OperatorRepresentation{HSR, O}) where {HSR, O} = eltype(O)
@inline eltype(lhs ::Type{OperatorRepresentation{HSR, O}}) where {HSR, O} = eltype(O)

function get_slice(opr ::OperatorRepresentation{HSR, O}, icol ::Integer) where {HSR, O}
  hsr = opr.hilbert_space_realization
  dim = length(hsr.basis_list)
  bcol = hsr.basis_list[icol]
  iter = ((get(hsr.basis_lookup, brow, -1), amplitude) for (brow, amplitude) in get_slice(opr.operator, bcol))
  return ((irow => amplitude) for (irow, amplitude) in iter if 1 <= irow <= dim)
end

"""
May contain duplicates
"""
function get_slice(irow ::Integer, opr ::OperatorRepresentation{HSR, O}) where {HSR, O}
  hsr = opr.hilbert_space_realization
  dim = length(hsr.basis_list)
  brow = hsr.basis_list[irow]
  iter = ((get(hsr.basis_lookup, bcol, -1), amplitude) for (bcol, amplitude) in get_slice(brow, opr.operator))
  return ((irow => amplitude) for (icol, amplitude) in iter if 1 <= icol <= dim)
end

function apply_unsafe!(out ::Vector{S1},
                       opr ::OperatorRepresentation{HSR, O},
                       state ::AbstractVector{S2};
                       range ::AbstractVector{<:Integer}=1:dimension(opr.hilbert_space_realization)
                       ) where {HSR, O, S1<:Number, S2<:Number}
  hsr = opr.hilbert_space_realization
  dim = length(hsr.basis_list)
  err = zero(Float64)
  for icol in range
    bcol = hsr.basis_list[icol]
    v = state[icol]
    for ((brow, _), amplitude) in get_slice(opr.operator, bcol)
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
                       range ::AbstractVector{<:Integer}=1:dimension(opr.hilbert_space_realization)
                       ) where {HSR, O, S1<:Number, S2<:Number}
  hsr = opr.hilbert_space_realization
  dim = length(hsr.basis_list)
  err = zero(Float64)
  for irow in range
    brow = hsr.basis_list[irow]
    v = state[irow]
    for (bcol, amplitude) in get_slice(brow, opr.operator)
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
#   hsr = opr.hilbert_space_realization
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
#       for (brow, amplitude) in get_slice(opr.operator, bcol)
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


import Base.*
function (*)(opr ::OperatorRepresentation{HSR, O}, state ::AbstractVector{S}) where {HSR, O, S<:Number}
  hsr = opr.hilbert_space_realization
  n = length(hsr.basis_list)
  T = promote_type(S, eltype(O))
  out = zeros(T, n)
  err = apply_unsafe!(out, opr, state)
  return out
end


import Base.*
function (*)(state ::AbstractVector{S}, opr ::OperatorRepresentation{HSR, O}) where {HSR, O, S<:Number}
  hsr = opr.hilbert_space_realization
  n = len(hsr.basis_list)
  T = promote_type(S, eltype(O))
  out = zeros(T, n)
  err = apply_unsafe!(out, state, opr)
  out
end


import SparseArrays.sparse
function sparse(opr::OperatorRepresentation{HSR, O}; tol ::Real=sqrt(eps(Float64))) where {HSR, O}
  S = eltype(opr)
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
