
export AbstractOperatorRepresentation
export spacetype, operatortype
export bintype
export get_space
export get_row, get_column
export sparse_serial, sparse_parallel

# AbstractOperatorRepresentation

abstract type AbstractOperatorRepresentation{T} <: AbstractMatrix{T} end

## typetraits

# a subclass of AbstractOperatorRepresentation should implement
# spacetype, operatortype, and get_space.

@inline spacetype(lhs::AbstractOperatorRepresentation{T}) where T = spacetype(typeof(lhs)) ::Type{<:AbstractHilbertSpaceRepresentation}
@inline operatortype(lhs ::AbstractOperatorRepresentation{T}) where T = operatortype(typeof(lhs)) ::Type{<:AbstractOperator}

#if not specialized
# @inline spacetype(lhs::Type{AbstractOperatorRepresentation}) = error("spacetype not implemented")
# @inline operatortype(lhs ::Type{AbstractOperatorRepresentation}) = error("operatortype not implemented")
#@inline get_space(lhs::AbstractOperatorRepresentation) = error("get_space not implemented")

@inline bintype(lhs ::AbstractOperatorRepresentation{T}) where T = bintype(typeof(lhs)) ::DataType
@inline bintype(lhs ::Type{<:AbstractOperatorRepresentation{T}}) where T = bintype(spacetype(lhs)) ::DataType

#@inline scalartype(lhs::AbstractOperatorRepresentation) ::DataType = scalartype(typeof(lhs)) ::DataType
@inline scalartype(lhs::AbstractOperatorRepresentation{T}) where T = T
@inline scalartype(lhs::Type{<:AbstractOperatorRepresentation{T}}) where T = T

import Base.size
@inline function size(arg::AbstractOperatorRepresentation{T}) ::Tuple{Int, Int} where T
  dim = dimension(get_space(arg))
  return (dim, dim)
end

@inline size(arg::AbstractOperatorRepresentation{T}, i::Integer) where T = size(arg)[i]


import Base.==
@inline function (==)(lhs ::AbstractOperatorRepresentation{T1},
                      rhs ::AbstractOperatorRepresentation{T2}) where {T1, T2}
  return ((get_space(lhs) == get_space(rhs)) && (lhs.operator == rhs.operator))
end

import Base.+, Base.-, Base.*

for uniop in [:+, :-]
  expr = :(
  @inline function ($uniop)(lhs ::AbstractOperatorRepresentation{T}) where T
    return represent(lhs.hilbert_space_representation, ($uniop)(lhs.operator))
  end
  )
  eval(expr)
end

for binop in [:+, :-, :*]
  expr = :(
  @inline function ($binop)(lhs ::AbstractOperatorRepresentation{T1},
                            rhs ::AbstractOperatorRepresentation{T2}) where {T1, T2}
    @boundscheck if (get_space(lhs) != get_space(rhs))
      throw(ArgumentError("The two OperatorRepresentation s do not have the same HilbertSpaceRepresentation"))
    end
    return represent(lhs.hilbert_space_representation, ($binop)(lhs.operator, rhs.operator))
  end
  )
  eval(expr)
end


import LinearAlgebra.ishermitian
function ishermitian(arg::AbstractOperatorRepresentation{S}) where S
  return ishermitian(arg.operator)
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


import LinearAlgebra.mul!
function LinearAlgebra.mul!(out ::AbstractVector{S1},
                            opr ::AbstractOperatorRepresentation{S2},
                            state ::AbstractVector{S3}) where {S1<:Number, S2<:Number, S3<:Number}
  fill!(out, zero(S1))
  apply!(out, opr, state)
  out
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
  spinlock = Threads.SpinLock()

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


"""
    apply!(out, opr, state)

Perform `out += opr * state`. Apply the operator representation `opr` to the
column vector `state` and *add* it to the column vector `out`.
Return sum of errors and sum of error-squared.
Call `apply_serial!` if `Threads.nthreads() == 1`, and `apply_parallel!` if greater.

# Arguments
- `out ::Vector{S1}`
- `opr ::OperatorRepresentation{HSR, O}`
- `state ::AbstractVector{S2}`
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
Call `apply_serial!` if `Threads.nthreads() == 1`, and `apply_parallel!` if greater.

# Arguments
- `out ::Vector{S1}`
- `state ::AbstractVector{S2}`
- `opr ::AbstractOperatorRepresentation{S}`
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
    apply_serial!(out, opr, state; range=1:size(opr, 2))

Perform `out += opr * state`. Apply the operator representation `opr` to the
column vector `state` and *add* it to the column vector `out`.
Return sum of errors and sum of error-squared.
Single-threaded version.

# Arguments
- `out ::Vector{S1}`
- `opr ::AbstractOperatorRepresentation{S}`
- `state ::AbstractVector{S2}`
- `range ::AbstractVector{<:Integer}=1:size(opr, 1)`
"""
function apply_serial!(out ::AbstractVector{S1},
                       opr ::AbstractOperatorRepresentation{S},
                       state ::AbstractVector{S2}
                      ) where {S, S1<:Number, S2<:Number}
  # w_r += sum_r  A_rc v_c
  nrows, ncols = size(opr)
  length(out) != nrows && throw(DimensionMismatch("out has length $(length(out)) != dimension $(nrows)"))
  length(state) != ncols && throw(DimensionMismatch("state has length $(length(state)) != dimension $(ncols)"))
  err, err_sq = zero(S1), zero(real(S1))
  for irow in 1:nrows
    for (icol, amplitude) in get_row_iterator(opr, irow)
      if 1 <= icol <= ncols
        @inbounds out[irow] += amplitude * state[icol]
      else
        @fastmath err += amplitude
        @fastmath err_sq += abs2(amplitude)
      end
    end
  end
  return (err, err_sq)
end

"""
    apply_serial!(out, state, opr; range=1:size(opr, 1))

Perform `out += state * opr`. Apply the operator representation `opr` to the
row vector `state` and *add* it to the row vector `out`.
Return sum of errors and sum of error-squared.
Single-threaded version.

# Arguments
- `out ::Vector{S1}`
- `state ::AbstractVector{S2}`
- `opr ::AbstractOperatorRepresentation{S}`
- `range ::AbstractVector{<:Integer}=1:size(opr, 1)`
"""
function apply_serial!(out ::AbstractVector{S1},
                       state ::AbstractVector{S2},
                       opr ::AbstractOperatorRepresentation{S}
                      ) where {S, S1<:Number, S2<:Number}
  # w_c += sum_r v_r A_rc
  nrows, ncols = size(opr)
  length(out) != ncols && throw(DimensionMismatch("out has length $(length(out)) != dimension $(ncols)"))
  length(state) != nrows && throw(DimensionMismatch("state has length $(length(state)) != dimension $(nrows)"))
  err, err_sq = zero(S1), zero(real(S1))
  for icol in 1:ncols
    for (irow, amplitude) in get_column_iterator(opr, icol)
      if 1 <= irow <= nrows
        @inbounds out[icol] += state[irow] * amplitude
      else
        @fastmath err += amplitude
        @fastmath err_sq += abs2(amplitude)
      end
    end
  end
  return (err, err_sq)
end


"""
    apply_parallel!(out, opr, state; range=1:size(opr, 2))

Perform `out += opr * state`. Apply the operator representation `opr` to the
column vector `state` and *add* it to the column vector `out`.
Return sum of errors and sum of error-squared.
Multi-threaded version.

# Arguments
- `out ::Vector{S1}`
- `opr ::AbstractOperatorRepresentation{S}`
- `state ::AbstractVector{S2}`
- `range ::AbstractVector{<:Integer}=1:size(opr, 2)`
"""
function apply_parallel!(out ::AbstractVector{S1},
                         opr ::AbstractOperatorRepresentation{S},
                         state ::AbstractVector{S2}
                         ) where {S, S1<:Number, S2<:Number}
  # w_r += sum_r  A_rc v_c
  nrows, ncols = size(opr)
  length(out) != nrows && throw(DimensionMismatch("out has length $(length(out)) != dimension $(nrows)"))
  length(state) != ncols && throw(DimensionMismatch("state has length $(length(state)) != dimension $(ncols)"))

  R1 = real(S1)

  nthreads = Threads.nthreads()
  local_err = zeros(S1, nthreads)
  local_err_sq = zeros(R1, nthreads)

  Threads.@threads for irow in 1:nrows
    id = Threads.threadid()
    l_err, l_err_sq = zero(S1), zero(R1)
    for (icol, amplitude) in get_row_iterator(opr, irow)
      if 1 <= icol <= ncols
        @inbounds out[irow] += amplitude * state[icol]
      else
        @fastmath l_err += amplitude
        @fastmath l_err_sq += abs2(amplitude)
      end
    end
    @fastmath local_err[id] += l_err
    @fastmath local_err_sq[id] += l_err_sq
  end
  @fastmath err = sum(local_err)
  @fastmath err_sq = sum(local_err_sq)
  return (err, err_sq)
end


"""
    apply_parallel!(out, state, opr; range=1:size(opr, 1))

Perform `out += state * opr`. Apply the operator representation `opr` to the
row vector `state` and *add* it to the row vector `out`.
Return sum of errors and sum of error-squared.
Multi-threaded version.

# Arguments
- `out ::Vector{S1}`
- `state ::AbstractVector{S2}`
- `opr ::OperatorRepresentation{HSR, O}`
- `range ::AbstractVector{<:Integer}=1:dimension(opr.hilbert_space_representation)`
"""
function apply_parallel!(out ::AbstractVector{S1},
                         state ::AbstractVector{S2},
                         opr ::AbstractOperatorRepresentation{S}
                         ) where {S, S1<:Number, S2<:Number}
  # w(c) += sum_(r) v(r) A(r,c)
  nrows, ncols = size(opr)
  length(out) != ncols && throw(DimensionMismatch("out has length $(length(out)) != dimension $(ncols)"))
  length(state) != nrows && throw(DimensionMismatch("state has length $(length(state)) != dimension $(nrows)"))

  R1 = real(S1)
  nthreads = Threads.nthreads()
  local_err = zeros(S1, nthreads)
  local_err_sq = zeros(R1, nthreads)

  Threads.@threads for icol in 1:ncols
    id = Threads.threadid()
    l_err, l_err_sq = zero(S1), zero(R1)
    for (irow, amplitude) in get_column_iterator(opr, icol)
      if 1 <= irow <= nrows
        @inbounds out[icol] += state[irow] * amplitude
      else
        @fastmath l_err += amplitude
        @fastmath l_err_sq += abs2(amplitude)
      end
    end
    @fastmath local_err[id] += l_err
    @fastmath local_err_sq[id] += l_err_sq
  end
  @fastmath err = sum(local_err)
  @fastmath err_sq = sum(local_err_sq)
  return (err, err_sq)
end


#
# """
#     apply_parallel!(out, opr, state; range=1:size(opr, 2))
#
# Perform `out += opr * state`. Apply the operator representation `opr` to the
# column vector `state` and *add* it to the column vector `out`.
# Return sum of errors and sum of error-squared.
# Multi-threaded version.
#
# # Arguments
# - `out ::Vector{S1}`
# - `opr ::AbstractOperatorRepresentation{S}`
# - `state ::AbstractVector{S2}`
# - `range ::AbstractVector{<:Integer}=1:size(opr, 2)`
# """
# function apply_parallel!(out ::AbstractVector{S1},
#                          opr ::AbstractOperatorRepresentation{S},
#                          state ::AbstractVector{S2};
#                          range ::AbstractVector{<:Integer}=1:size(opr, 2)
#                          ) where {S, S1<:Number, S2<:Number}
#   nrows, ncols = size(opr)
#   length(out) != nrows && throw(DimensionMismatch("out has length $(length(out)) != dimension $(nrows)"))
#   length(state) != ncols && throw(DimensionMismatch("state has length $(length(state)) != dimension $(ncols)"))
#   for i in range # check bounds first before multithreading
#     (1<=i<=ncols) || throw(BoundsError("attempt to access $ncols-element $(typeof(state)) at index [$i]"))
#   end
#
#   err ::S1 = zero(S1)
#   err_sq ::real(S1) = zero(real(S1))
#
#   nblocks = Threads.nthreads()
#   block_ranges = splitrange(range, nblocks)
#   spinlock = Threads.SpinLock()
#   Threads.@threads for iblock in 1:nblocks
#     subrange = block_ranges[iblock]
#     local_err, local_err_sq, local_offdiag = zero(S1), zero(real(S1)), Dict{Int, S1}()
#     for icol in subrange
#       @inbounds v = state[icol]
#       for (irow, amplitude) in get_column_iterator(opr, icol)
#         if 1 <= irow <= nrows
#           if irow in subrange
#             @inbounds out[irow] += amplitude * v
#           else
#             local_offdiag[irow] = get(local_offdiag, irow, zero(S1)) + amplitude * v
#           end
#         else
#           local_err += amplitude * v
#           local_err_sq += abs(amplitude * v)^2
#         end
#       end
#     end
#     lock(spinlock)
#     err += local_err
#     err_sq += local_err_sq
#     for (irow, ampl) in local_offdiag
#       @inbounds out[irow] += ampl
#     end
#     unlock(spinlock)
#   end
#   return (err ::S1, err_sq ::real(S1))
# end


# """
#     apply_parallel!(out, state, opr; range=1:size(opr, 1))
#
# Perform `out += state * opr`. Apply the operator representation `opr` to the
# row vector `state` and *add* it to the row vector `out`.
# Return sum of errors and sum of error-squared.
# Multi-threaded version.
#
# # Arguments
# - `out ::Vector{S1}`
# - `state ::AbstractVector{S2}`
# - `opr ::OperatorRepresentation{HSR, O}`
# - `range ::AbstractVector{<:Integer}=1:dimension(opr.hilbert_space_representation)`
# """
# function apply_parallel!(out ::AbstractVector{S1},
#                          state ::AbstractVector{S2},
#                          opr ::AbstractOperatorRepresentation{S};
#                          range ::AbstractVector{<:Integer}=1:size(opr, 1)
#                          ) where {S, S1<:Number, S2<:Number}
#   nrows, ncols = size(opr)
#   length(out) != ncols && throw(ArgumentError("out has length $(length(out)) != dimension $(ncols)"))
#   length(state) != nrows && throw(ArgumentError("state has length $(length(state)) != dimension $(nrows)"))
#   for i in range  # check bounds first before multithreading
#     (1<=i<=nrows) || throw(BoundsError("attempt to access $nrows-element $(typeof(state)) at index [$i]"))
#   end
#   err, err_sq = zero(S1), zero(real(S1))
#   nblocks = Threads.nthreads()
#   block_ranges = splitrange(range, nblocks)
#   spinlock = Threads.SpinLock()
#   Threads.@threads for iblock in 1:nblocks
#     subrange = block_ranges[iblock]
#     local_err, local_err_sq, local_offdiag = zero(S1), zero(real(S1)), Dict{Int, S1}()
#     for irow in subrange
#       @inbounds v = state[irow]
#       for (icol, amplitude) in get_row_iterator(opr, irow)
#         if 1 <= icol <= ncols
#           if icol in subrange
#             @inbounds out[icol] += amplitude * v
#           else
#             local_offdiag[icol] = get(local_offdiag, icol, zero(S1)) + v * amplitude
#           end
#         else
#           local_err += v * amplitude
#           local_err_sq += abs(v * amplitude)^2
#         end
#       end
#     end
#     lock(spinlock)
#     err += local_err
#     err_sq += local_err_sq
#     for (icol, ampl) in local_offdiag
#       @inbounds out[icol] += ampl
#     end
#     unlock(spinlock)
#   end
#   return (err, err_sq)
# end
