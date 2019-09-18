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
  while nx <= length(x) && ny <= length(y)
    if x[nx] < y[ny]
      push!(z, x[nx])
      nx += 1
    elseif y[ny] < x[nx]
      push!(z, y[ny])
      ny += 1
    else
      push!(z, x[nx])
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

function bitcount(i ::UInt64) ::Int
	i = i - ((i >> 1) & 0x5555555555555555);
	i = (i & 0x3333333333333333) + ((i >> 2) & 0x3333333333333333);
	i = (i + (i >> 4)) & 0x0f0f0f0f0f0f0f0f;
	i = i + (i >> 8);
	i = i + (i >> 16);
	i = i + (i >> 32);
  return Int(i & 0x7f);
end
