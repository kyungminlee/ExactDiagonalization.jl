export AbstractFrozenSortedArray
export FrozenSortedArrayIndex

abstract type AbstractFrozenSortedArray{K, V} <: AbstractDict{K, V} end

struct FrozenSortedArrayIndex{K} <:AbstractFrozenSortedArray{K, Int}
  keys ::Vector{K}
  function FrozenSortedArrayIndex{K}(keys::Vector{K}) where K
    issorted(keys) || throw(ArgumentError("vals must be sorted"))

    if length(keys) > 1
      for i in 2:length(keys)
        keys[i-1] == keys[i] && throw(ArgumentError("vals contains duplicates $(keys[i])"))
      end
    end
    return new{K}(keys)
  end

  function FrozenSortedArrayIndex(keys::Vector{K}) where K
    return FrozenSortedArrayIndex{K}(keys)
  end
end

@inline function fs_index(arr::FrozenSortedArrayIndex{K}, key) ::Int where K
  idxrange ::UnitRange{Int} = searchsorted(arr.keys, key) ::UnitRange{Int}
  isempty(idxrange) && return 0
  return idxrange[1]
end

import Base.getindex
@inline function Base.getindex(arr::FrozenSortedArrayIndex{K}, key) ::Int where K
  idx = fs_index(arr, key)
  @boundscheck idx <= 0 && throw(KeyError("key $key not found"))
  return idx
end

import Base.haskey
@inline function Base.haskey(arr::FrozenSortedArrayIndex{K}, key) ::Bool where K
  idx = fs_index(arr, key)
  return idx > 0
end

import Base.get
@inline function Base.get(arr::FrozenSortedArrayIndex{K}, key, default::Int) ::Int where K
  idx = fs_index(arr, key)
  return (idx > 0) ? idx : default
end

import Base.iterate
function iterate(iter::FrozenSortedArrayIndex{K}, state::Int=1) ::Union{Nothing, Tuple{Pair{K, Int}, Int}} where {K}
  state > length(iter.keys) && return nothing
  return ((iter.keys[state] => state), state+1)
end

import Base.eltype
@inline eltype(iter::Type{FrozenSortedArrayIndex{K}}) where K = Pair{K, Int}

import Base.length
@inline length(iter::FrozenSortedArrayIndex{K}) where K = length(iter.keys)

import Base.keys
@inline keys(iter::FrozenSortedArrayIndex{K}) where K = iter.keys
