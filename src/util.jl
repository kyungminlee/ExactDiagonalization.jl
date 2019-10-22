export make_bitmask

function make_bitmask(msb ::Integer,
                      lsb ::Integer=0;
                      dtype ::DataType=UInt)
  mask = dtype(0x1) << msb - dtype(0x1)
  submask = dtype(0x1) << lsb - dtype(0x1)
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
