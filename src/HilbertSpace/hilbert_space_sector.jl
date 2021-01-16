export HilbertSpaceSector
export scalartype
export qntype
export basespace
export bitwidth


"""
    HilbertSpaceSector{QN}

Hilbert space sector.
"""
struct HilbertSpaceSector{HS<:AbstractHilbertSpace, QN<:Tuple{Vararg{<:AbstractQuantumNumber}}}<:AbstractHilbertSpace{QN}
    parent::HS
    allowed_quantum_numbers::Set{QN}

    """
        HilbertSpaceSector(parent::HilbertSpace{QN}) where QN
    """
    function HilbertSpaceSector(parent::HS) where {HS<:AbstractHilbertSpace}
        QN = qntype(HS)
        sectors = quantum_number_sectors(parent)
        return new{HS, QN}(parent, Set(sectors))
    end

    """
        HilbertSpaceSector(parent::HilbertSpace{QN}, allowed::Integer) where {QN<:Tuple{<:Integer}}
    """
    function HilbertSpaceSector(parent::HS, allowed::Integer) where {HS<:AbstractHilbertSpace{<:Tuple{<:Integer}}}
        QN = qntype(HS)
        sectors = Set{QN}(quantum_number_sectors(parent))
        return new{HS, QN}(parent, intersect(sectors, Set([(allowed,)])))
    end

    """
        HilbertSpaceSector(parent::HilbertSpace{QN}, allowed::QN) where QN
    """
    function HilbertSpaceSector(parent::AbstractHilbertSpace{QN}, allowed::QN) where {QN}
        HS = typeof(parent)
        sectors = Set{QN}(quantum_number_sectors(parent))
        return new{HS, QN}(parent, intersect(sectors, Set([allowed])))
    end

    """
        HilbertSpaceSector(parent::HilbertSpace{QN}, allowed::Union{AbstractSet{<:Integer}, AbstractVector{<:Integer}}) where QN
    """
    function HilbertSpaceSector(
        parent::AbstractHilbertSpace{QN},
        allowed::Union{AbstractSet{<:Integer}, AbstractVector{<:Integer}}
    ) where {QN}
        HS = typeof(parent)
        sectors = Set{QN}(quantum_number_sectors(parent))
        return new{HS, QN}(parent, intersect(sectors, Set((x,) for x in allowed)))
    end

    function HilbertSpaceSector(
        parent::AbstractHilbertSpace{QN},
        allowed::Union{AbstractSet{QN}, AbstractVector{QN}}
    ) where QN
        HS = typeof(parent)
        sectors = Set{QN}(quantum_number_sectors(parent))
        return new{HS, QN}(parent, intersect(sectors, Set(allowed)))
    end
end


"""
    scalartype(arg::Type{HilbertSpaceSector{HS, QN}})

Returns the scalar type of the given hilbert space sector type.
For HilbertSpaceSector{QN}, it is always `Bool`.
"""
scalartype(::Type{HilbertSpaceSector{HS, QN}}) where {HS, QN} = Bool
scalartype(::HilbertSpaceSector{HS, QN}) where {HS, QN} = Bool


"""
    valtype(arg::Type{HilbertSpaceSector{HS, QN}})

Returns the `valtype` (scalar type) of the given hilbert space sector type.
For HilbertSpaceSector{QN}, it is always `Bool`.
"""
Base.valtype(::Type{HilbertSpaceSector{HS, QN}}) where {HS, QN} = Bool
Base.valtype(::HilbertSpaceSector{HS, QN}) where {HS, QN} = Bool


"""
    qntype(arg::Type{HilbertSpaceSector{QN}})

Returns the quantum number type of the given hilbert space sector type.
"""
qntype(::Type{HilbertSpaceSector{HS, QN}}) where {HS, QN} = QN
qntype(::HilbertSpaceSector{HS, QN}) where {HS, QN} = QN


"""
    basespace(hss)

Get the base space of the `HilbertSpaceSector`, which is
its parent `HilbertSpace` (with no quantum number restriction).
"""
basespace(hss::HilbertSpaceSector{HS, QN}) where {HS, QN} = basespace(hss.parent)::HS

#bitwidth(hss::HilbertSpaceSector) = bitwidth(basespace(hss))

function Base.:(==)(lhs::HilbertSpaceSector{HS, Q1}, rhs::HilbertSpaceSector{HS, Q2}) where {HS, Q1, Q2}
    return (
        basespace(lhs) == basespace(rhs) &&
        lhs.allowed_quantum_numbers == rhs.allowed_quantum_numbers
    )
end


for fname in [
    :bitwidth,
    :get_bitmask,
    :quantum_number_sectors,
    :get_quantum_number,
    :extract,
    :compress,
    :update,
    :get_state_index,
    :get_state,
]
    @eval begin
        """
            $($fname)(hss::HilbertSpaceSector, args...;kwargs...)

        Call `$($fname)` with basespace of `hss`.
        """
        @inline $fname(hss::HilbertSpaceSector, args...;kwargs...) = $fname(hss.parent, args...; kwargs...)
    end
end
