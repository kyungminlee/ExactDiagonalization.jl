export AbstractQuantumNumber
export AbstractHilbertSpace
export AbstractSiteType


## TODO: Think about this
AbstractQuantumNumber = Integer

abstract type AbstractHilbertSpace{QN<:Tuple{Vararg{<:AbstractQuantumNumber}}} end

abstract type AbstractSiteType end
