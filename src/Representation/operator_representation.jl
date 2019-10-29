export OperatorRepresentation
export represent
export apply!, apply_serial!, apply_parallel!

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


@inline spacetype(lhs::Type{OperatorRepresentation{HSR, O}}) where {HSR<:HilbertSpaceRepresentation, O<:AbstractOperator} = HSR
@inline operatortype(lhs ::Type{OperatorRepresentation{HSR, O}}) where {HSR<:HilbertSpaceRepresentation, O<:AbstractOperator} = O
@inline get_space(lhs ::OperatorRepresentation{HSR, O}) where {HSR<:HilbertSpaceRepresentation, O<:AbstractOperator} = lhs.hilbert_space_representation ::HSR


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

"""
    apply!(out, opr, state; range=1:size(opr, 2))

Perform `out += opr * state`. Apply the operator representation `opr` to the
column vector `state` and *add* it to the column vector `out`.
Return sum of errors and sum of error-squared.
Call `apply_serial!` if `Threads.nthreads() == 1`, and `apply_parallel!` if greater.

# Arguments
- `out ::Vector{S1}`
- `opr ::OperatorRepresentation{HSR, O}`
- `state ::AbstractVector{S2}`
- `range ::AbstractVector{<:Integer}=1:size(opr, 2)`
"""
function apply!(out ::Vector{S1},
                opr ::OperatorRepresentation{HSR, O},
                state ::AbstractVector{S2};
                range ::AbstractVector{<:Integer}=1:size(opr, 2)
               ) where {HSR, O, S1<:Number, S2<:Number}
  return Threads.nthreads() == 1 ?
             apply_serial!(out, opr, state; range=range) :
             apply_parallel!(out, opr, state; range=range)
end

"""
    apply!(out, opr, state; range=1:size(opr, 1))

Perform `out += opr * state`. Apply the operator representation `opr` to the
row vector `state` and *add* it to the row vector `out`.
Return sum of errors and sum of error-squared.
Call `apply_serial!` if `Threads.nthreads() == 1`, and `apply_parallel!` if greater.

# Arguments
- `out ::Vector{S1}`
- `state ::AbstractVector{S2}`
- `opr ::OperatorRepresentation{HSR, O}`
- `range ::AbstractVector{<:Integer}=1:size(opr, 1)`
"""
function apply!(out ::Vector{S1},
                state ::AbstractVector{S2},
                opr ::OperatorRepresentation{HSR, O};
                range ::AbstractVector{<:Integer}=1:size(opr, 1)
               ) where {HSR, O, S1<:Number, S2<:Number}
  return Threads.nthreads() == 1 ?
             apply_serial!(out, state, opr; range=range) :
             apply_parallel!(out, state, opr; range=range)
end

"""
    apply_serial!(out, opr, state; range=1:size(opr, 2))

Perform `out += opr * state`. Apply the operator representation `opr` to the
column vector `state` and *add* it to the column vector `out`.
Return sum of errors and sum of error-squared.
Single-threaded version.

# Arguments
- `out ::Vector{S1}`
- `opr ::OperatorRepresentation{HSR, O}`
- `state ::AbstractVector{S2}`
- `range ::AbstractVector{<:Integer}=1:size(opr, 2)`
"""
function apply_serial!(out ::Vector{S1},
                       opr ::OperatorRepresentation{HSR, O},
                       state ::AbstractVector{S2};
                       range ::AbstractVector{<:Integer}=1:size(opr, 2)
                      ) where {HSR, O, S1<:Number, S2<:Number}
  hsr = opr.hilbert_space_representation
  nrows, ncols = size(opr)
  length(out) != nrows && throw(ArgumentError("out has length $(length(out)) != dimension $(nrows)"))
  length(state) != ncols && throw(ArgumentError("state has length $(length(state)) != dimension $(ncols)"))
  err, err_sq = zero(S1), zero(real(S1))
  for icol in range
    v = state[icol]
    for (irow, amplitude) in get_column_iterator(opr, icol; include_all=true)
      if 1 <= irow <= nrows
        @inbounds out[irow] += amplitude * v
      else
        err += amplitude * v
        err_sq += abs(amplitude * v)^2
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
- `opr ::OperatorRepresentation{HSR, O}`
- `range ::AbstractVector{<:Integer}=1:size(opr, 1)`
"""
function apply_serial!(out ::Vector{S1},
                       state ::AbstractVector{S2},
                       opr ::OperatorRepresentation{HSR, O};
                       range ::AbstractVector{<:Integer}=1:size(opr, 1)
                      ) where {HSR, O, S1<:Number, S2<:Number}
  hsr = opr.hilbert_space_representation
  nrows, ncols = size(opr)
  length(out) != ncols && throw(ArgumentError("out has length $(length(out)) != dimension $(ncols)"))
  length(state) != nrows && throw(ArgumentError("state has length $(length(state)) != dimension $(nrows)"))
  err, err_sq = zero(S1), zero(real(S1))
  for irow in range
    v = state[irow]
    for (icol, amplitude) in get_row_iterator(opr, irow; include_all=true)
      if 1 <= icol <= ncols
        @inbounds out[icol] += v * amplitude
      else
        err += amplitude * v
        err_sq += abs(amplitude * v)^2
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
- `opr ::OperatorRepresentation{HSR, O}`
- `state ::AbstractVector{S2}`
- `range ::AbstractVector{<:Integer}=1:size(opr, 2)`
"""
function apply_parallel!(out ::Vector{S1},
                         opr ::OperatorRepresentation{HSR, O},
                         state ::AbstractVector{S2};
                         range ::AbstractVector{<:Integer}=1:size(opr, 2)
                         ) where {HSR, O, S1<:Number, S2<:Number}
  hsr = opr.hilbert_space_representation
  nrows, ncols = size(opr)
  length(out) != nrows && throw(ArgumentError("out has length $(length(out)) != dimension $(nrows)"))
  length(state) != ncols && throw(ArgumentError("state has length $(length(state)) != dimension $(ncols)"))
  for i in range # check bounds first before multithreading
    (1<=i<=ncols) || throw(BoundsError("attempt to access $ncols-element $(typeof(state)) at index [$i]"))
  end
  err, err_sq = zero(S1), zero(real(S1))
  nblocks = Threads.nthreads()
  block_ranges = splitrange(range, nblocks)
  spinlock = Threads.SpinLock()
  Threads.@threads for iblock in 1:nblocks
    subrange = block_ranges[iblock]
    local_err, local_err_sq, local_offdiag = zero(S1), zero(real(S1)), Dict{Int, S1}()
    for icol in subrange
      @inbounds v = state[icol]
      for (irow, amplitude) in get_column_iterator(opr, icol; include_all=true)
        if 1 <= irow <= nrows
          if irow in subrange
            @inbounds out[irow] += amplitude * v
          else
            local_offdiag[irow] = get(local_offdiag, irow, zero(S1)) + amplitude * v
          end
        else
          local_err += amplitude * v
          local_err_sq += abs(amplitude * v)^2
        end
      end
    end
    lock(spinlock)
    err += local_err
    err_sq += local_err_sq
    for (irow, ampl) in local_offdiag
      @inbounds out[irow] += ampl
    end
    unlock(spinlock)
  end
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
function apply_parallel!(out ::Vector{S1},
                         state ::AbstractVector{S2},
                         opr ::OperatorRepresentation{HSR, O};
                         range ::AbstractVector{<:Integer}=1:size(opr, 1)
                         ) where {HSR, O, S1<:Number, S2<:Number}
  hsr = opr.hilbert_space_representation
  nrows, ncols = size(opr)
  length(out) != ncols && throw(ArgumentError("out has length $(length(out)) != dimension $(ncols)"))
  length(state) != nrows && throw(ArgumentError("state has length $(length(state)) != dimension $(nrows)"))
  for i in range  # check bounds first before multithreading
    (1<=i<=nrows) || throw(BoundsError("attempt to access $nrows-element $(typeof(state)) at index [$i]"))
  end
  err, err_sq = zero(S1), zero(real(S1))
  nblocks = Threads.nthreads()
  block_ranges = splitrange(range, nblocks)
  spinlock = Threads.SpinLock()
  Threads.@threads for iblock in 1:nblocks
    subrange = block_ranges[iblock]
    local_err, local_err_sq, local_offdiag = zero(S1), zero(real(S1)), Dict{Int, S1}()
    for irow in subrange
      @inbounds v = state[irow]
      for (icol, amplitude) in get_row_iterator(opr, irow; include_all=true)
        if 1 <= icol <= ncols
          if icol in subrange
            @inbounds out[icol] += amplitude * v
          else
            local_offdiag[icol] = get(local_offdiag, icol, zero(S1)) + v * amplitude
          end
        else
          local_err += v * amplitude
          local_err_sq += abs(v * amplitude)^2
        end
      end
    end
    lock(spinlock)
    err += local_err
    err_sq += local_err_sq
    for (icol, ampl) in local_offdiag
      @inbounds out[icol] += ampl
    end
    unlock(spinlock)
  end
  return (err, err_sq)
end


import Base.*
function (*)(opr ::OperatorRepresentation{HSR, O}, state ::AbstractVector{S}) where {HSR, O, S<:Number}
  hsr = opr.hilbert_space_representation
  n = dimension(hsr)
  T = promote_type(S, scalartype(O))
  out = zeros(T, n)
  err = apply!(out, opr, state)
  return out
end


import Base.*
function (*)(state ::AbstractVector{S}, opr ::OperatorRepresentation{HSR, O}) where {HSR, O, S<:Number}
  hsr = opr.hilbert_space_representation
  n = dimension(hsr)
  T = promote_type(S, scalartype(O))
  out = zeros(T, n)
  err = apply!(out, state, opr)
  out
end
