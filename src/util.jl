export IntegerModulo

export make_bitmask
export choptol!
export merge_vec

struct IntegerModulo{N} <: Integer
  value ::Int
  IntegerModulo{N}(value::Integer) where N = new{N}(mod(value, N))
end

import Base.+, Base.-, Base.*
(+)(lhs::IntegerModulo{N}) where N = lhs
(-)(lhs::IntegerModulo{N}) where N = IntegerModulo{N}(-lhs.value)

(+)(lhs::IntegerModulo{N}, rhs::IntegerModulo{N}) where N = IntegerModulo{N}(lhs.value + rhs.value)
(+)(lhs::IntegerModulo{N}, rhs::Integer) where N = IntegerModulo{N}(lhs.value + rhs)
(+)(lhs::Integer, rhs::IntegerModulo{N}) where N = IntegerModulo{N}(lhs + rhs.value)

(*)(lhs::IntegerModulo{N}, rhs::IntegerModulo{N}) where N = IntegerModulo{N}(lhs.value * rhs.value)
(*)(lhs::IntegerModulo{N}, rhs::Integer) where N = IntegerModulo{N}(lhs.value * rhs)
(*)(lhs::Integer, rhs::IntegerModulo{N}) where N = IntegerModulo{N}(lhs * rhs.value)




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


export elmax, elmin
export elmaximum, elminimum

elmax(x::S, y::S) where {S<:Number} = max(x,y)
elmax(x::T, y::T) where {T<:Tuple{<:Number}} = (max(first(x), first(y)),)

function elmax(x ::T, y::T) where {T<:Tuple{<:Number, <:Number, Vararg{<:Number}}}
  return (max(first(x), first(y)), elmax(x[2:end], y[2:end])...)
end

elmin(x::S, y::S) where {S<:Number} = max(x,y)
elmin(x::T, y::T) where {T<:Tuple{<:Number}} = (max(first(x), first(y)),)

function elmin(x ::T, y::T) where {T<:Tuple{<:Number, <:Number, Vararg{<:Number}}}
  return (max(first(x), first(y)), elmax(x[2:end], y[2:end])...)
end

elmaximum(arr::AbstractArray) = mapreduce(identity, elmax, arr)
elminimum(arr::AbstractArray) = mapreduce(identity, elmin, arr)
