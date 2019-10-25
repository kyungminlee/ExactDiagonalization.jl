export OperatorRepresentation
export represent
export apply!, apply_threaded!

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


function apply!(out ::Vector{S1},
                opr ::OperatorRepresentation{HSR, O},
                state ::AbstractVector{S2};
                range ::AbstractVector{<:Integer}=1:dimension(opr.hilbert_space_representation)
               ) where {HSR, O, S1<:Number, S2<:Number}
  hsr = opr.hilbert_space_representation
  nrows, ncols = size(opr)
  length(out) != nrows && throw(ArgumentError("out has length $(length(out)) != dimension $(nrows)"))
  length(state) != ncols && throw(ArgumentError("state has length $(length(state)) != dimension $(ncols)"))
  dim = dimension(hsr)
  @assert nrows == dim && ncols == dim

  err = zero(S1)
  err_sq = zero(real(S1))
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
  return err, err_sq
end

function apply!(out ::Vector{S1},
                state ::AbstractVector{S2},
                opr ::OperatorRepresentation{HSR, O};
                range ::AbstractVector{<:Integer}=1:dimension(opr.hilbert_space_representation)
               ) where {HSR, O, S1<:Number, S2<:Number}
  hsr = opr.hilbert_space_representation
  nrows, ncols = size(opr)
  length(out) != ncols && throw(ArgumentError("out has length $(length(out)) != dimension $(ncols)"))
  length(state) != nrows && throw(ArgumentError("state has length $(length(state)) != dimension $(nrows)"))
  dim = dimension(hsr)
  @assert nrows == dim && ncols == dim
            
  err = zero(S1)
  err_sq = zero(real(S1))
  for irow in range
    v = state[irow]
    for (icol, amplitude) in get_row_iterator(opr, irow; include_all=true)
      if 1 <= icol <= dim
        @inbounds out[icol] += v * amplitude
      else
        err += amplitude * v
        err_sq += abs(amplitude * v)^2
      end
    end
  end
  return err, err_sq
end


"""
    splitblock

Split n into b blocks

"""
function splitblock(n ::Integer, b ::Integer) ::Vector{Int}
  (n < 0) && throw(ArgumentError("n cannot be negative"))
  (b <= 0) && throw(ArgumentError("b must be positive"))
  blocksize = n ÷ b
  blocks = blocksize * ones(Int, b)
  r = n - (blocksize * b)
  blocks[1:r] .+= 1
  return blocks
end

function splitrange(range::AbstractVector{<:Integer}, b::Integer)
  #nblocks = min(b, length(range))
  nblocks = b
  block_sizes = splitblock(length(range), nblocks)
  block_offsets = cumsum([1, block_sizes...])
  return [range[block_offsets[iblock]:block_offsets[iblock+1]-1] for iblock in 1:nblocks]
end

function apply_threaded!(out ::Vector{S1},
                         opr ::OperatorRepresentation{HSR, O},
                         state ::AbstractVector{S2}
                         ;range ::AbstractVector{<:Integer}=1:size(opr, 2)
                         ) where {HSR, O, S1<:Number, S2<:Number}
  hsr = opr.hilbert_space_representation
  nrows, ncols = size(opr)
  length(out) != nrows && throw(ArgumentError("out has length $(length(out)) != dimension $(nrows)"))
  length(state) != ncols && throw(ArgumentError("state has length $(length(state)) != dimension $(ncols)"))
  dim = dimension(hsr)
  @assert nrows == dim && ncols == dim

  # check bounds first before multithreading
  for i in range
    (1<=i<=ncols) || throw(BoundsError("attempt to access $ncols-element $(typeof(state)) at index [$i]"))
  end
  
  nblocks = Threads.nthreads()
  block_ranges = splitrange(range, nblocks)

  err = zero(S1)
  err_sq = zero(real(S1))
  spinlock = Threads.SpinLock()
  Threads.@threads for iblock in 1:nblocks
    local_err = zero(S1)
    local_err_sq = zero(S1)
    local_offdiag = Dict{Int, S1}()

    subrange = block_ranges[iblock]
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



function apply_threaded!(out ::Vector{S1},
                         state ::AbstractVector{S2},
                         opr ::OperatorRepresentation{HSR, O};
                         range ::AbstractVector{<:Integer}=1:size(opr, 1)
                         ) where {HSR, O, S1<:Number, S2<:Number}
  hsr = opr.hilbert_space_representation
  nrows, ncols = size(opr)
  length(out) != ncols && throw(ArgumentError("out has length $(length(out)) != dimension $(ncols)"))
  length(state) != nrows && throw(ArgumentError("state has length $(length(state)) != dimension $(nrows)"))

  dim = dimension(hsr)
  @assert nrows == dim && ncols == dim

  # check bounds first before multithreading
  for i in range
    (1<=i<=nrows) || throw(BoundsError("attempt to access $nrows-element $(typeof(state)) at index [$i]"))
  end
  
  nblocks = Threads.nthreads()
  block_ranges = splitrange(range, nblocks)

  err = zero(S1)
  err_sq = zero(real(S1))
  
  spinlock = Threads.SpinLock()

  Threads.@threads for iblock in 1:nblocks
    local_err = zero(S1)
    local_err_sq = zero(S1)
    local_offdiag = Dict{Int, S1}()

    subrange = block_ranges[iblock]
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
  T = promote_type(S, eltype(O))
  out = zeros(T, n)
  err = apply!(out, opr, state)
  return out
end


import Base.*
function (*)(state ::AbstractVector{S}, opr ::OperatorRepresentation{HSR, O}) where {HSR, O, S<:Number}
  hsr = opr.hilbert_space_representation
  n = dimension(hsr)
  T = promote_type(S, eltype(O))
  out = zeros(T, n)
  err = apply!(out, state, opr)
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
