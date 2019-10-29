
export AbstractOperatorRepresentation
export spacetype, operatortype
export bintype
export get_space
export get_row, get_column
export sparse_serial, sparse_parallel

# AbstractOperatorRepresentation

abstract type AbstractOperatorRepresentation end

## typetraits

# a subclass of AbstractOperatorRepresentation should implement
# spacetype, operatortype, and get_space.

@inline spacetype(lhs::AbstractOperatorRepresentation) = spacetype(typeof(lhs)) ::Type{<:AbstractHilbertSpaceRepresentation}
@inline operatortype(lhs ::AbstractOperatorRepresentation) = operatortype(typeof(lhs)) ::Type{<:AbstractOperator}

#if not specialized
# @inline spacetype(lhs::Type{AbstractOperatorRepresentation}) = error("spacetype not implemented")
# @inline operatortype(lhs ::Type{AbstractOperatorRepresentation}) = error("operatortype not implemented")
#@inline get_space(lhs::AbstractOperatorRepresentation) = error("get_space not implemented")

@inline bintype(lhs ::AbstractOperatorRepresentation) = bintype(typeof(lhs)) ::DataType
@inline bintype(lhs ::Type{<:AbstractOperatorRepresentation}) = bintype(spacetype(lhs)) ::DataType

#@inline scalartype(lhs::AbstractOperatorRepresentation) ::DataType = scalartype(typeof(lhs)) ::DataType
@inline scalartype(lhs::AbstractOperatorRepresentation) ::DataType = promote_type(scalartype(spacetype(lhs)), scalartype(operatortype(lhs))) ::DataType
@inline scalartype(lhs::Type{<:AbstractOperatorRepresentation}) ::DataType = promote_type(scalartype(spacetype(lhs)), scalartype(operatortype(lhs))) ::DataType

import Base.size
@inline function size(arg::AbstractOperatorRepresentation) ::Tuple{Int, Int}
  dim = dimension(get_space(arg))
  return (dim, dim)
end

@inline size(arg::AbstractOperatorRepresentation, i::Integer) = size(arg)[i]


import Base.==
@inline function (==)(lhs ::AbstractOperatorRepresentation,
              rhs ::AbstractOperatorRepresentation)
  return ((get_space(lhs) == get_space(rhs)) && (lhs.operator == rhs.operator))
end

import Base.+, Base.-, Base.*

for uniop in [:+, :-]
  expr = :(
  @inline function ($uniop)(lhs ::AbstractOperatorRepresentation)
    return represent(lhs.hilbert_space_representation, ($uniop)(lhs.operator))
  end
  )
  eval(expr)
end

for binop in [:+, :-, :*]
  expr = :(
  @inline function ($binop)(lhs ::AbstractOperatorRepresentation,
                            rhs ::AbstractOperatorRepresentation)
    @boundscheck if (get_space(lhs) != get_space(rhs))
      throw(ArgumentError("The two OperatorRepresentation s do not have the same HilbertSpaceRepresentation"))
    end
    return represent(lhs.hilbert_space_representation, ($binop)(lhs.operator, rhs.operator))
  end
  )
  eval(expr)
end



import Base.Matrix
function Matrix(opr::AbstractOperatorRepresentation)
  S = scalartype(opr)
  m, n = size(opr)
  out = zeros(S, (m, n))
  Threads.@threads for icol in 1:n
    for (irow, ampl) in get_column_iterator(opr, icol; include_all=false)
      out[irow, icol] += ampl
    end
  end
  return out
end


import SparseArrays.sparse

function sparse(opr::AbstractOperatorRepresentation; tol ::Real=sqrt(eps(Float64)))
  sp = Threads.nthreads() == 1 ? sparse_serial : sparse_parallel
  return sp(opr; tol=tol)
end

function sparse_serial(opr::AbstractOperatorRepresentation; tol ::Real=sqrt(eps(Float64)))
  S = scalartype(opr)
  m, n = size(opr)
  colptr = zeros(Int, n+1)
  rowval = Int[]
  nzval = S[]

  colptr[1] = 1
  for icol in 1:n
    colvec = Dict{Int, S}()
    for (irow, ampl) in get_column_iterator(opr, icol; include_all=false)
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

function sparse_parallel(opr::AbstractOperatorRepresentation; tol ::Real=sqrt(eps(Float64)))
  S = scalartype(opr)
  m, n = size(opr)

  colsize = zeros(Int, n)

  nblocks = Threads.nthreads()
  block_ranges = splitrange(1:n, nblocks) # guaranteed to be ordered
  spinlock = Threads.SpinLock()

  local_rowval = Vector{Int}[Int[] for i in 1:nblocks]
  local_nzval = Vector{S}[S[] for i in 1:nblocks]

  Threads.@threads for iblock in 1:nblocks
    subrange = block_ranges[iblock]
    for icol in subrange
      colvec = Dict{Int, S}()
      for (irow, ampl) in get_column_iterator(opr, icol; include_all=false)
        colvec[irow] = get(colvec, irow, zero(S)) + ampl
      end
      choptol!(colvec, tol)
      sorted_items = sort(collect(colvec), by=item->item[1])
      colsize[icol] = length(sorted_items)
      append!(local_rowval[iblock], irow for (irow, ampl) in sorted_items)
      append!(local_nzval[iblock], ampl for (irow, ampl) in sorted_items)
    end
  end

  colptr = cumsum(Int[1, colsize...])
  rowval = vcat(local_rowval...)
  nzval = vcat(local_nzval...)
  return SparseMatrixCSC{S, Int}(m, n, colptr, rowval, nzval)
end



function get_row(opr ::AbstractOperatorRepresentation, irow::Integer)
  S = scalartype(opr)
  Z = zero(S)
  dim = size(opr, 2)
  items = Dict{Int, S}()
  for (icol, val) in get_row_iterator(opr, irow; include_all=false)
    items[icol] = get(items, icol, Z) + val
  end
  choptol!(items, sqrt(eps(Float64)))
  return sparsevec(items, dim)
end

function get_column(opr ::AbstractOperatorRepresentation, icol::Integer)
  S = scalartype(opr)
  Z = zero(S)
  dim = size(opr, 1)
  items = Dict{Int, S}()
  for (irow, val) in get_column_iterator(opr, icol; include_all=false)
    items[irow] = get(items, irow, Z) + val
  end
  choptol!(items, sqrt(eps(Float64)))
  return sparsevec(items, dim)
end



import Base.getindex
@inline function getindex(oprep ::AbstractOperatorRepresentation, ::Colon, icol::Integer)
  return get_column(oprep, icol)
end

@inline function getindex(oprep ::AbstractOperatorRepresentation, irow::Integer, ::Colon)
  return get_row(oprep, irow)
end

@inline function getindex(oprep ::AbstractOperatorRepresentation, ::Colon, ::Colon)
  return sparse(oprep)
end
