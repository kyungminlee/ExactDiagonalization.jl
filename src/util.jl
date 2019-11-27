export make_bitmask
export choptol!
export merge_vec

function make_bitmask(msb ::Integer,
                      binary_type::Type{BR}=UInt) where {BR <:Unsigned}
  mask = BR(0x1) << msb - BR(0x1)
  return mask
end

function make_bitmask(msb ::Integer,
                      lsb ::Integer,
                      binary_type::Type{BR}=UInt) where {BR <:Unsigned}
  mask = BR(0x1) << msb - BR(0x1)
  submask = BR(0x1) << lsb - BR(0x1)
  return mask ⊻ submask
end

function merge_vec(x::Vector{T}, y::Vector{T})::Vector{T} where {T}
  (nx, ny) = (1, 1)
  z = T[]
  sizehint!(z, length(x) + length(y))
  while nx <= length(x) && ny <= length(y)
    if x[nx] < y[ny]
      push!(z, x[nx])
      nx += 1
    elseif y[ny] < x[nx]
      push!(z, y[ny])
      ny += 1
    else
      push!(z, x[nx])
      push!(z, y[ny])
      nx += 1
      ny += 1
    end
  end
  if nx <= length(x)
    [ push!(z, x[n]) for n = nx:length(x) ]
  elseif ny <= length(y)
    [ push!(z, y[n]) for n = ny:length(y) ]
  end
  return z
end

function choptol!(d ::Dict{K, V}, tol::Real) where {K, V<:Number}
  to_delete = K[k for (k, v) in d if abs(v) < tol]
  for k in to_delete
    delete!(d, k)
  end
  d
end



"""
    splitblock

Split n into b blocks.

# Arguments
- `n ::Integer`: the number of elements to split.
- `b ::Integer`: the number of blocks.
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
  nblocks = b
  block_sizes = splitblock(length(range), nblocks)
  block_offsets = cumsum([1, block_sizes...])
  return [range[block_offsets[iblock]:block_offsets[iblock+1]-1] for iblock in 1:nblocks]
end
