export AbstractHilbertSpaceRepresentation
export scalartype, bintype


"""
    AbstractHilbertSpaceRepresentation{S}
"""
abstract type AbstractHilbertSpaceRepresentation{S<:Number} end


import Base.valtype

scalartype(lhs ::AbstractHilbertSpaceRepresentation{S}) where S = S
valtype(lhs ::AbstractHilbertSpaceRepresentation{S}) where S = S
bintype(lhs ::AbstractHilbertSpaceRepresentation{S}) where S = bintype(typeof(lhs)) ::DataType


bitwidth(lhs ::AbstractHilbertSpaceRepresentation{S}) where S = bitwidth(basespace(lhs)) ::Int
