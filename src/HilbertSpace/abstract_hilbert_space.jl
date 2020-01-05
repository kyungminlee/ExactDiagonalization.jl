export AbstractHilbertSpace
abstract type AbstractHilbertSpace end

## TODO: Think about this
AbstractQuantumNumber = Union{Int, SVector{N, Int} where N}

export AbstractSiteType

abstract type AbstractSiteType end
