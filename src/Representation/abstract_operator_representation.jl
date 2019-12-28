
export AbstractOperatorRepresentation
export spacetype, operatortype
export bintype
export get_space
export dimension, bitwidth
export get_row, get_column
export sparse_serial, sparse_parallel


"""
    AbstractOperatorRepresentation{S}
"""
abstract type AbstractOperatorRepresentation{S} <: AbstractMatrix{S} end


## typetraits

import Base.valtype
valtype(lhs::Type{<:AbstractOperatorRepresentation{T}}) where T = T
valtype(lhs::AbstractOperatorRepresentation{T}) where T = T

scalartype(lhs::Type{<:AbstractOperatorRepresentation{T}}) where T = T
scalartype(lhs::AbstractOperatorRepresentation{T}) where T = T

bintype(lhs ::Type{<:AbstractOperatorRepresentation{T}}) where T = bintype(spacetype(lhs)) ::DataType
bintype(lhs ::AbstractOperatorRepresentation{T}) where T = bintype(typeof(lhs)) ::DataType


# a subclass of AbstractOperatorRepresentation should implement
# spacetype, operatortype, and get_space.
spacetype(lhs::AbstractOperatorRepresentation{T}) where T = spacetype(typeof(lhs)) ::Type{<:AbstractHilbertSpaceRepresentation}
operatortype(lhs ::AbstractOperatorRepresentation{T}) where T = operatortype(typeof(lhs)) ::Type{<:AbstractOperator}


dimension(lhs::AbstractOperatorRepresentation{S}) where S = dimension(get_space(lhs))
bitwidth(lhs::AbstractOperatorRepresentation{S}) where S = bitwidth(get_space(lhs))


import Base.size

function size(arg::AbstractOperatorRepresentation{T}) ::Tuple{Int, Int} where T
  dim = dimension(get_space(arg))
  return (dim, dim)
end

size(arg::AbstractOperatorRepresentation{T}, i::Integer) where T = size(arg)[i]


import Base.==
function (==)(lhs ::AbstractOperatorRepresentation{T1},
                      rhs ::AbstractOperatorRepresentation{T2}) where {T1, T2}
  return ((get_space(lhs) == get_space(rhs)) && (lhs.operator == rhs.operator))
end


import Base.+, Base.-, Base.*, Base./, Base.\

for uniop in [:+, :-]
  @eval begin
    function ($uniop)(lhs ::AbstractOperatorRepresentation{T}) where T
      return represent(get_space(lhs), ($uniop)(lhs.operator))
    end
  end
end

for binop in [:+, :-, :*]
  @eval begin
    function ($binop)(lhs ::AbstractOperatorRepresentation{T1},
                      rhs ::AbstractOperatorRepresentation{T2}) where {T1, T2}
      @boundscheck if (get_space(lhs) != get_space(rhs))
        throw(ArgumentError("The two OperatorRepresentation's do not have the same HilbertSpaceRepresentation"))
      end
      return represent(get_space(lhs), simplify(($binop)(lhs.operator, rhs.operator)))
    end
  end
end

function (*)(lhs ::AbstractOperatorRepresentation{T}, rhs ::Number) where {T}
  return represent(get_space(lhs), simplify(lhs.operator * rhs))
end

function (*)(lhs ::Number, rhs ::AbstractOperatorRepresentation{T}) where {T}
  return represent(get_space(rhs), simplify(lhs * rhs.operator))
end

function (/)(lhs ::AbstractOperatorRepresentation{T}, rhs ::Number) where {T}
  return represent(get_space(lhs), simplify(lhs.operator / rhs))
end

function (\)(lhs ::Number, rhs ::AbstractOperatorRepresentation{T}) where {T}
  return represent(get_space(rhs), simplify(lhs \ rhs.operator))
end

function simplify(arg::AbstractOperatorRepresentation{T}) where {T}
  return represent(get_space(arg), simplify(arg.operator))
end


import LinearAlgebra.ishermitian
function ishermitian(arg::AbstractOperatorRepresentation{S}) where S
  return ishermitian(arg.operator)
end


import LinearAlgebra.mul!
function LinearAlgebra.mul!(out ::AbstractVector{S1},
                            opr ::AbstractOperatorRepresentation{S2},
                            state ::AbstractVector{S3}) where {S1<:Number, S2<:Number, S3<:Number}
  fill!(out, zero(S1))
  apply!(out, opr, state)
  out
end


import Base.Matrix
function Matrix(opr::AbstractOperatorRepresentation{S}) where S
  m, n = size(opr)
  out = zeros(S, (m, n))
  Threads.@threads for icol in 1:n
    for (irow, ampl) in get_column_iterator(opr, icol)
      if 1 <= irow <= m
        out[irow, icol] += ampl
      end
    end
  end
  return out
end


import SparseArrays.sparse
function SparseArrays.sparse(opr::AbstractOperatorRepresentation{S}; tol ::Real=sqrt(eps(Float64))) where S
  sp = Threads.nthreads() == 1 ? sparse_serial : sparse_parallel
  return sp(opr; tol=tol)
end


function sparse_serial(opr::AbstractOperatorRepresentation{S}; tol ::Real=sqrt(eps(Float64))) where S
  m, n = size(opr)
  colptr = zeros(Int, n+1)
  rowval = Int[]
  nzval = S[]

  colptr[1] = 1
  for icol in 1:n
    colvec = Dict{Int, S}()
    for (irow, ampl) in get_column_iterator(opr, icol)
      if 1 <= irow <= m
        colvec[irow] = get(colvec, irow, zero(S)) + ampl
      end
    end
    choptol!(colvec, tol)
    colptr[icol+1] = colptr[icol] + length(colvec)
    sorted_items = sort(collect(colvec), by = item -> item[1])
    append!(rowval, irow for (irow, ampl) in sorted_items)
    append!(nzval, ampl for (irow, ampl) in sorted_items)
  end
  return SparseMatrixCSC{S, Int}(m, n, colptr, rowval, nzval)
end


function sparse_parallel(opr::AbstractOperatorRepresentation{S}; tol ::Real=sqrt(eps(Float64))) where S
  m, n = size(opr)

  colsize = zeros(Int, n)

  nblocks = Threads.nthreads()
  block_ranges = splitrange(1:n, nblocks) # guaranteed to be ordered

  local_rowval = Vector{Int}[Int[] for i in 1:nblocks]
  local_nzval = Vector{S}[S[] for i in 1:nblocks]

  Threads.@threads for iblock in 1:nblocks
    subrange = block_ranges[iblock]
    for icol in subrange
      colvec = Dict{Int, S}()
      for (irow, ampl) in get_column_iterator(opr, icol)
        if 1 <= irow <= m
          colvec[irow] = get(colvec, irow, zero(S)) + ampl
        end
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


function get_row(opr ::AbstractOperatorRepresentation{S}, irow::Integer) where S
  Z = zero(S)
  dim = size(opr, 2)
  items = Dict{Int, S}()
  for (icol, val) in get_row_iterator(opr, irow)
    if 1 <= icol <= dim
      items[icol] = get(items, icol, Z) + val
    end
  end
  choptol!(items, sqrt(eps(Float64)))
  return sparsevec(items, dim)
end


function get_column(opr ::AbstractOperatorRepresentation{S}, icol::Integer) where S
  Z = zero(S)
  dim = size(opr, 1)
  items = Dict{Int, S}()
  for (irow, val) in get_column_iterator(opr, icol)
    if 1 <= irow <= dim
      items[irow] = get(items, irow, Z) + val
    end
  end
  choptol!(items, sqrt(eps(Float64)))
  return sparsevec(items, dim)
end


import Base.getindex
@inline function getindex(oprep ::AbstractOperatorRepresentation{S}, ::Colon, icol::Integer) where S
  return get_column(oprep, icol)
end

@inline function getindex(oprep ::AbstractOperatorRepresentation{S}, irow::Integer, ::Colon) where S
  return get_row(oprep, irow)
end

@inline function getindex(oprep ::AbstractOperatorRepresentation{S}, ::Colon, ::Colon) where S
  return sparse(oprep)
end

@inline function getindex(oprep ::AbstractOperatorRepresentation{S}, irow::Integer, icol::Integer) where S
  return get_element(oprep, irow, icol)
end


"""
    apply!(out, opr, state)

Perform `out += opr * state`. Apply the operator representation `opr` to the
column vector `state` and *add* it to the column vector `out`.
Return sum of errors and sum of error-squared.
Call [`apply_serial!`](@ref) if `Threads.nthreads() == 1`, and [`apply_parallel!`](@ref) otherwise.
"""
function apply!(out ::AbstractVector{S1},
                opr ::AbstractOperatorRepresentation{S},
                state ::AbstractVector{S2}
               ) where {S, S1<:Number, S2<:Number}
  return Threads.nthreads() == 1 ?
             apply_serial!(out, opr, state) :
             apply_parallel!(out, opr, state)
end


"""
    apply!(out, opr, state)

Perform `out += opr * state`. Apply the operator representation `opr` to the
row vector `state` and *add* it to the row vector `out`.
Return sum of errors and sum of error-squared.
Call [`apply_serial!`](@ref) if `Threads.nthreads() == 1`, and [`apply_parallel!`](@ref) otherwise.
"""
function apply!(out ::AbstractVector{S1},
                state ::AbstractVector{S2},
                opr ::AbstractOperatorRepresentation{S}
               ) where {HSR, S, O, S1<:Number, S2<:Number}
  return Threads.nthreads() == 1 ?
             apply_serial!(out, state, opr) :
             apply_parallel!(out, state, opr)
end


"""
    apply_serial!(out, opr, state)

Perform `out += opr * state`. Apply the operator representation `opr` to the
column vector `state` and *add* it to the column vector `out`.
Return sum of errors and sum of error-squared.
Single-threaded version.
"""
function apply_serial!(out ::AbstractVector{S1},
                       opr ::AbstractOperatorRepresentation{S},
                       state ::AbstractVector{S2}
                      ) where {S, S1<:Number, S2<:Number}
  # w_r += sum_r  A_rc v_c
  nrows, ncols = size(opr)
  length(out) != nrows && throw(DimensionMismatch("out has length $(length(out)) != dimension $(nrows)"))
  length(state) != ncols && throw(DimensionMismatch("state has length $(length(state)) != dimension $(ncols)"))
  for irow in 1:nrows
    for (icol::Int, amplitude::S) in get_row_iterator(opr, irow)
      if 1 <= icol <= ncols
        @inbounds out[irow] += amplitude * state[icol]
      end
    end
  end
  out
end


"""
    apply_serial!(out, state, opr)

Perform `out += state * opr`. Apply the operator representation `opr` to the
row vector `state` and *add* it to the row vector `out`.
Return sum of errors and sum of error-squared.
Single-threaded version.
"""
function apply_serial!(out ::AbstractVector{S1},
                       state ::AbstractVector{S2},
                       opr ::AbstractOperatorRepresentation{S}
                      ) where {S, S1<:Number, S2<:Number}
  # w_c += sum_r v_r A_rc
  nrows, ncols = size(opr)
  length(out) != ncols && throw(DimensionMismatch("out has length $(length(out)) != dimension $(ncols)"))
  length(state) != nrows && throw(DimensionMismatch("state has length $(length(state)) != dimension $(nrows)"))
  for icol in 1:ncols
    for (irow::Int, amplitude::S) in get_column_iterator(opr, icol)
      if 1 <= irow <= nrows
        @inbounds out[icol] += state[irow] * amplitude
      end
    end
  end
  out
end


"""
    apply_parallel!(out, opr, state)

Perform `out += opr * state`. Apply the operator representation `opr` to the
column vector `state` and *add* it to the column vector `out`.
Return sum of errors and sum of error-squared.
Multi-threaded version.
"""
function apply_parallel!(out ::AbstractVector{S1},
                         opr ::AbstractOperatorRepresentation{S},
                         state ::AbstractVector{S2}
                         ) where {S, S1<:Number, S2<:Number}
  # w_r += sum_r  A_rc v_c
  nrows, ncols = size(opr)
  length(out) != nrows && throw(DimensionMismatch("out has length $(length(out)) != dimension $(nrows)"))
  length(state) != ncols && throw(DimensionMismatch("state has length $(length(state)) != dimension $(ncols)"))

  R = real(S)
  Threads.@threads for irow in 1:nrows
    for (icol::Int, amplitude::S) in get_row_iterator(opr, irow)
      if 1 <= icol <= ncols
        @inbounds out[irow] += amplitude * state[icol]
      end
    end
  end
  out
end


"""
    apply_parallel!(out, state, opr)

Perform `out += state * opr`. Apply the operator representation `opr` to the
row vector `state` and *add* it to the row vector `out`.
Return sum of errors and sum of error-squared.
Multi-threaded version.
"""
function apply_parallel!(out ::AbstractVector{S1},
                         state ::AbstractVector{S2},
                         opr ::AbstractOperatorRepresentation{S}
                         ) where {S, S1<:Number, S2<:Number}
  # w(c) += sum_(r) v(r) A(r,c)
  nrows, ncols = size(opr)
  length(out) != ncols && throw(DimensionMismatch("out has length $(length(out)) != dimension $(ncols)"))
  length(state) != nrows && throw(DimensionMismatch("state has length $(length(state)) != dimension $(nrows)"))

  R = real(S)
  Threads.@threads for icol in 1:ncols
    for (irow::Int, amplitude::S) in get_column_iterator(opr, icol)
      if 1 <= irow <= nrows
        @inbounds out[icol] += state[irow] * amplitude
      end
    end
  end
  out
end
