abstract type AbstractFrozenSortedArray{K, V} <: AbstractDict{K, V} end

struct FrozenSortedArrayIndex{K} <:AbstractFrozenSortedArray{K, Int}
  keys ::Vector{K}
  function FrozenSortedArrayIndex{K}(keys::Vector{K}) where K
    if issorted(keys)
      return new{K}(keys)
    else
      throw(ArgumentError("vals must be sorted"))
    end
  end

  function FrozenSortedArrayIndex(keys::Vector{K}) where K
    return FrozenSortedArrayIndex{K}(keys)
  end
end

function fs_index(arr::FrozenSortedArrayIndex{K}, key) ::Int where K
  idxrange ::UnitRange{Int} = searchsorted(arr.keys, key) ::UnitRange{Int}
  isempty(idxrange) && return 0
  return idxrange[1]
end

function Base.getindex(arr::FrozenSortedArrayIndex{K}, key) ::Int where K
  idx = fs_index(arr, key)
  @boundscheck idx <= 0 && throw(ArgumentError("key $key not found"))
  return idx
end

function Base.haskey(arr::FrozenSortedArrayIndex{K}, key) ::Bool where K
  idx = fs_index(arr, key)
  return idx > 0
end

import Base.iterate
function iterate(iter::FrozenSortedArrayIndex{K}, state::Int=1) ::Union{Nothing, Tuple{Pair{K, Int}, Int}} where {K}
  state > length(iter.keys) && return nothing
  return ((iter.keys[state] => state), state+1)
end
