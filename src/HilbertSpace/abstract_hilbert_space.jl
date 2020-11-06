export AbstractQuantumNumber
export AbstractHilbertSpace

## TODO: Think about this
AbstractQuantumNumber = Integer

abstract type AbstractHilbertSpace{QN<:Tuple{Vararg{<:AbstractQuantumNumber}}} end
