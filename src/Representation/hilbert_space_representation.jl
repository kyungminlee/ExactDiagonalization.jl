export HilbertSpaceRepresentation
export dimension
export represent, represent_array, represent_dict


"""
    HilbertSpaceRepresentation{HS, BR, DictType}

# Fields
```
hilbert_space :: HS
basis_list    :: Vector{BR}
basis_lookup  :: DictType
```
"""
struct HilbertSpaceRepresentation{HS <:AbstractHilbertSpace,
                                  BR <:Unsigned,
                                  DictType <:AbstractDict{BR, <:Int}
                                  } <:AbstractHilbertSpaceRepresentation{Bool}
  hilbert_space ::HS
  basis_list ::Vector{BR}
  basis_lookup ::DictType

  function HilbertSpaceRepresentation(
      hilbert_space ::HilbertSpace{QN},
      basis_list::AbstractVector{BR},
      basis_lookup::DictType) where {QN, BR <:Unsigned, DictType<:AbstractDict{BR, <:Int}}
    if sizeof(BR)*8 <= bitwidth(hilbert_space)
      # equality added such that the MSB checks overflow
      throw(ArgumentError("type $(BR) not enough to represent the hilbert space (need $(bitwidth(hilbert_space)) bits)"))
    end
    return new{HilbertSpace{QN}, BR, DictType}(hilbert_space, basis_list, basis_lookup)
  end

  function HilbertSpaceRepresentation(
      hilbert_space_sector ::HilbertSpaceSector{QN},
      basis_list::AbstractVector{BR},
      basis_lookup::DictType) where {QN, BR <:Unsigned, DictType<:AbstractDict{BR, <:Int}}
    if sizeof(BR)*8 <= bitwidth(hilbert_space_sector)
      # equality added such that the MSB checks overflow
      throw(ArgumentError("type $(BR) not enough to represent the hilbert space (need $(bitwidth(hilbert_space_sector)) bits)"))
    end
    return new{HilbertSpace{QN}, BR, DictType}(hilbert_space_sector.parent, basis_list, basis_lookup)
  end
end


import Base.valtype
valtype(lhs ::Type{HilbertSpaceRepresentation{HS, BR, DT}}) where {HS, BR, DT} = Bool
scalartype(lhs ::Type{HilbertSpaceRepresentation{HS, BR, DT}}) where {HS, BR, DT} = Bool
bintype(lhs ::Type{HilbertSpaceRepresentation{HS, BR, DT}}) where {HS, BR, DT} = BR


basespace(lhs::HilbertSpaceRepresentation{HS, BR, DT}) where {HS, BR, DT} = lhs.hilbert_space ::HS


"""
    dimension

Dimension of the Concrete Hilbert space, i.e. number of basis vectors.
"""
dimension(hsr ::HilbertSpaceRepresentation) = length(hsr.basis_list)


import Base.==
function (==)(lhs ::HilbertSpaceRepresentation{H1, B1, D1},
              rhs ::HilbertSpaceRepresentation{H2, B2, D2}) where {H1, B1, D1, H2, B2, D2}
  return basespace(lhs) == basespace(rhs) && (lhs.basis_list == rhs.basis_list)
end


function checkvalidbasis(hsr::HilbertSpaceRepresentation{HS, BR, DT}) where {HS, BR, DT}
  for (ivec, bvec) in enumerate(hsr.basis_list)
    ivec2 = hsr.basis_lookup[bvec]
    @assert ivec == ivec2
  end
end


function hs_get_basis_list(hs ::HilbertSpace{QN}, binary_type::Type{BR}=UInt)::Vector{BR} where {QN, BR<:Unsigned}
  if sizeof(BR) * 8 <= bitwidth(hs)
    throw(ArgumentError("type $(BR) not enough to represent the hilbert space (need $(bitwidth(hs)) bits)"))
  end
  basis_list = BR[]
  for indexarray in keys(hs)
    push!(basis_list, compress(hs, indexarray, BR))
  end
  sort!(basis_list)
  return basis_list
end


"""
    hs_get_basis_list(hss, binary_type=UInt)

Generate a basis for the `HilbertSpaceSector`.
"""
function hs_get_basis_list(hss::HilbertSpaceSector{QN}, binary_type::Type{BR}=UInt) ::Vector{BR} where {QN, BR<:Unsigned}
  hs = hss.parent
  if sizeof(BR) * 8 <= bitwidth(hs)
    throw(ArgumentError("type $(BR) not enough to represent the hilbert space (need $(bitwidth(hs)) bits)"))
  end
  if isempty(intersect(hss.allowed_quantum_numbers, quantum_number_sectors(hs)))
    return BR[]
  end

  quantum_numbers = Vector{QN}[QN[state.quantum_number for state in site.states] for site in hs.sites]

  n_sites = length(hs.sites)

  # `qn_possible[i]` contains all possible quantum numbers in the subspace
  # spanned by sites 1 ⊗ 2 ⊗ ... ⊗ (i-1).
  # `qn_requested[i]` contains quantum numbers that are requested in the subspace
  # spanned by sites 1 ⊗ 2 ⊗ ... ⊗ (i-1).
  # `qn_schedule` is the intersection of these two, and are the ones we need to
  # generate the basis.
  qn_possible ::Vector{Vector{QN}} = let qn_possible = Vector{Vector{QN}}(undef, n_sites+1)
    qn_possible[1] = [tuplezero(QN)]
    for i in eachindex(hs.sites)
      a = Set{QN}()
      for qi in quantum_numbers[i], qa in qn_possible[i]
        push!(a, qa .+ qi)
      end
      qn_possible[i+1] = sort(collect(a))
    end
    qn_possible
  end

  qn_requested ::Vector{Vector{QN}} = let qn_requested = Vector{Vector{QN}}(undef, n_sites+1)
    qn_requested[n_sites+1] = sort(collect(hss.allowed_quantum_numbers))
    for i in n_sites:-1:1
      a = Set{QN}()
      for qi in quantum_numbers[i], qa in qn_requested[i+1]
        push!(a, qa .- qi)
      end
      qn_requested[i] = sort(collect(a))
    end
    qn_requested
  end

  qn_schedule = [intersect(x,y) for (x, y) in zip(qn_requested, qn_possible)]

  sector_basis_list = Dict{QN, Vector{BR}}(tuplezero(QN) => BR[BR(0x0)])
  new_sector_basis_list = Dict{QN, Vector{BR}}()

  #sl = Threads.SpinLock()
  for i in 1:n_sites
    empty!(new_sector_basis_list)
    # @threads
    for q in qn_schedule[i+1]
      new_sector_basis_list_q = BR[]
      for (i_state, q_curr) in enumerate(quantum_numbers[i])
        q_prev ::QN = q .- q_curr
        if haskey(sector_basis_list, q_prev)
          append!(new_sector_basis_list_q,
                  (s | (BR(i_state-1) << hs.bitoffsets[i])) for s in sector_basis_list[q_prev])
        end
      end
      #lock(sl)
      new_sector_basis_list[q] = new_sector_basis_list_q
      #unlock(sl)
    end
    sector_basis_list, new_sector_basis_list = new_sector_basis_list, sector_basis_list
  end

  basis_list ::Vector{BR} = let
    basis_list = BR[]
    sector_basis_list ::Dict{QN, Vector{BR}} = sector_basis_list

    qs = collect(keys(sector_basis_list))
    for q in qs
      states = pop!(sector_basis_list, q)
      basis_list = merge_vec(basis_list, states)
      GC.gc()
    end
    basis_list
  end
  @assert issorted(basis_list)
  return basis_list
end


"""
    represent(hs, binary_type=UInt)

Make a HilbertSpaceRepresentation with all the basis vectors of the specified HilbertSpace.
This function defaults to [`represent_array`](@ref).
"""
represent(hs::AbstractHilbertSpace, binary_type::Type{BR}=UInt) where {BR <:Unsigned} = represent_array(hs, binary_type)


"""
    represent_array(hs, binary_type=UInt)

Make a HilbertSpaceRepresentation with all the basis vectors of the specified HilbertSpace
using [`FrozenSortedArrayIndex`](@ref).
"""
function represent_array(hs::AbstractHilbertSpace, binary_type::Type{BR}=UInt) where {BR <:Unsigned}
  basis_list = hs_get_basis_list(hs, BR)
  basis_lookup = FrozenSortedArrayIndex{BR}(basis_list)
  return HilbertSpaceRepresentation(basespace(hs), basis_list, basis_lookup)
end


"""
    represent(hs, basis_list)

Make a HilbertSpaceRepresentation with the provided list of basis vectors.
This defaults to [`represent_array`](@ref).
"""
represent(hs ::AbstractHilbertSpace, basis_list ::AbstractVector{BR}) where {BR<:Unsigned} = represent_array(hs, basis_list)


"""
    represent_array(hs, basis_list)

Make a HilbertSpaceRepresentation with the provided list of basis vectors
using [`FrozenSortedArrayIndex`](@ref).
"""
function represent_array(hs ::AbstractHilbertSpace,
                         basis_list ::AbstractVector{BR}) where {BR<:Unsigned}
  if !issorted(basis_list)
    basis_list = sort(basis_list)
  end
  basis_lookup = FrozenSortedArrayIndex{BR}(basis_list)
  return HilbertSpaceRepresentation(basespace(hs), basis_list, basis_lookup)
end


"""
    represent_dict(hs, binary_type=UInt)

Make a HilbertSpaceRepresentation with the provided list of basis vectors
using `Dict{binary_type, Int}`.
"""
function represent_dict(hs::AbstractHilbertSpace, binary_type::Type{BR}=UInt) where {BR<:Unsigned}
  basis_list = hs_get_basis_list(hs, BR)
  basis_lookup = Dict{BR, Int}()
  sizehint!(basis_lookup, length(basis_list) + length(basis_list) ÷ 4)
  for (ibasis, basis) in enumerate(basis_list)
    basis_lookup[basis] = ibasis
  end
  return HilbertSpaceRepresentation(basespace(hs), basis_list, basis_lookup)
end


"""
    represent_dict(hs, basis_list)

Make a HilbertSpaceRepresentation with the provided list of basis vectors using `Dict`.
"""
function represent_dict(hs ::AbstractHilbertSpace,
                        basis_list ::AbstractVector{BR}) where {BR<:Unsigned}
  if !issorted(basis_list)
    basis_list = sort(basis_list)
  end
  basis_lookup = Dict{BR, Int}()
  sizehint!(basis_lookup, length(basis_list) + length(basis_list) ÷ 4)
  for (ibasis, basis) in enumerate(basis_list)
    basis_lookup[basis] = ibasis
  end
  return HilbertSpaceRepresentation(basespace(hs), basis_list, basis_lookup)
end
