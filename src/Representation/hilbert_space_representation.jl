export HilbertSpaceRepresentation
export dimension
export represent, represent_dict

struct HilbertSpaceRepresentation{HS <:AbstractHilbertSpace,
                                  BR <:Unsigned,
                                  DictType <:AbstractDict{BR, <:Int}
                                  } <:AbstractHilbertSpaceRepresentation
  hilbert_space ::HS
  basis_list ::Vector{BR}
  basis_lookup ::DictType

  function HilbertSpaceRepresentation(
      hilbert_space ::HilbertSpace{QN},
      basis_list::AbstractVector{BR},
      basis_lookup::DictType) where {QN, BR <:Unsigned, DictType<:AbstractDict{BR, <:Int}}
    #if count_ones(typemax(BR)) <= bitwidth(hilbert_space)
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
    #if count_ones(typemax(BR)) <= bitwidth(hilbert_space_sector)
    if sizeof(BR)*8 <= bitwidth(hilbert_space_sector)
      # equality added such that the MSB checks overflow
      throw(ArgumentError("type $(BR) not enough to represent the hilbert space (need $(bitwidth(hilbert_space_sector)) bits)"))
    end
    return new{HilbertSpace{QN}, BR, DictType}(hilbert_space_sector.parent, basis_list, basis_lookup)
  end
end

@inline scalartype(lhs ::Type{HilbertSpaceRepresentation{HS, BR, DT}}) where {HS, BR, DT} = Bool
@inline bintype(lhs ::Type{HilbertSpaceRepresentation{HS, BR, DT}}) where {HS, BR, DT} = BR ::DataType

@inline basespace(lhs::HilbertSpaceRepresentation{HS, BR, DT}) where {HS, BR, DT} = lhs.hilbert_space ::HS

@inline bitwidth(hss::HilbertSpaceSector{QN}) where QN = bitwidth(hss.parent) ::Int


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

"""
    dimension

Dimension of the Concrete Hilbert space, i.e. number of basis vectors.
"""
@inline dimension(hsr ::HilbertSpaceRepresentation) = length(hsr.basis_list)



function hs_get_basis_list(hs ::HilbertSpace{QN}; BR ::DataType=UInt) ::Vector{BR} where {QN}
  if sizeof(BR) * 8 <= bitwidth(hs)
    throw(ArgumentError("type $(BR) not enough to represent the hilbert space (need $(bitwidth(hs)) bits)"))
  end
  basis_list = BR[]
  for indexarray in keys(hs)
    push!(basis_list, compress(hs, indexarray))
  end
  sort!(basis_list)
  return basis_list
end


function hs_get_basis_list(hss::HilbertSpaceSector{QN}; BR::DataType=UInt) ::Vector{BR} where {QN}
  hs = hss.parent
  if sizeof(BR) * 8 <= bitwidth(hs)
    throw(ArgumentError("type $(BR) not enough to represent the hilbert space (need $(bitwidth(hs)) bits)"))
  end
  allowed = hss.allowed_quantum_numbers
  sectors = Set(quantum_number_sectors(hs))
  if isempty(intersect(allowed, sectors))
    return BR[]
  end

  quantum_numbers = [[state.quantum_number for state in site.states] for site in hs.sites]
  possible_quantum_numbers = [Set([zero(QN)])]  # PQN[i]: possible QN left of i

  n_sites = length(hs.sites)
  for i in 1:n_sites
    pq = Set{QN}(q1 + q2 for q1 in possible_quantum_numbers[i], q2 in quantum_numbers[i])
    push!(possible_quantum_numbers, pq)
  end

  function generate(i ::Int, allowed ::AbstractSet{QN})
    (i == 0) && return (zero(QN) in allowed) ? Dict(zero(QN) => [BR(0x0)]) : Dict()
    allowed_prev = Set{QN}()
    for q1 in quantum_numbers[i], q2 in allowed
      q = q2 - q1
      (q in possible_quantum_numbers[i]) && push!(allowed_prev, q)
    end
    result_prev = generate(i-1, allowed_prev)
    result = Dict{QN, Vector{BR}}()
    for (i_state, q_curr) in enumerate(quantum_numbers[i])
      for (q_prev, states_prev) in result_prev
        q = q_prev + q_curr
        if q in allowed
          if !haskey(result, q)
            result[q] = BR[]
          end
          append!(result[q], (s | (BR(i_state-1) << hs.bitoffsets[i])) for s in states_prev)
        end
      end
    end
    return result
  end
  basis_list ::Vector{BR} = let
    basis_list = BR[]
    result = generate(n_sites, allowed)
    for (q, states) in result
      basis_list = merge_vec(basis_list, states)
    end
    basis_list
  end
  sort!(basis_list)
  return basis_list
end


"""
    represent(hs; BR=UInt)

Make a HilbertSpaceRepresentation with all the basis vectors of the specified HilbertSpaceSector.

# Arguments
- `hs ::AbstractHilbertSpace`
- `BR ::DataType=UInt`: Binary representation type
"""
function represent(hs::AbstractHilbertSpace; BR::DataType=UInt)
  basis_list = hs_get_basis_list(hs; BR=BR)
  basis_lookup = FrozenSortedArrayIndex{BR}(basis_list)
  return HilbertSpaceRepresentation(basespace(hs), basis_list, basis_lookup)
end

"""
    represent(hs, basis_list)

Make a HilbertSpaceRepresentation with the provided list of basis vectors

# Arguments
- `hs ::AbstractHilbertSpace`
- `basis_list ::AbstractVector{BR}`
"""
function represent(hs ::AbstractHilbertSpace,
                   basis_list ::AbstractVector{BR}) where {BR<:Unsigned}
  if !issorted(basis_list)
    basis_list = sort(basis_list)
  end
  basis_lookup = FrozenSortedArrayIndex{BR}(basis_list)
  return HilbertSpaceRepresentation(basespace(hs), basis_list, basis_lookup)
end


"""
    represent_dict(hs; BR=UInt)

Make a HilbertSpaceRepresentation with all the basis vectors of the specified HilbertSpace.

# Arguments
- `hs ::AbstractHilbertSpace`
- `BR ::DataType=UInt`: Binary representation type
"""
function represent_dict(hs::AbstractHilbertSpace; BR::DataType=UInt)
  basis_list = hs_get_basis_list(hs; BR=BR)
  basis_lookup = Dict{BR, Int}(basis => ibasis for (ibasis, basis) in enumerate(basis_list))
  return HilbertSpaceRepresentation(basespace(hs), basis_list, basis_lookup)
end


"""
    represent_dict(hs, basis_list)

Make a HilbertSpaceRepresentation with the provided list of basis vectors using Dict

# Arguments
- `hs ::HilbertSpace{QN}`: Abstract Hilbert space
- `basis_list ::AbstractVector{BR}`
"""
function represent_dict(hs ::AbstractHilbertSpace,
                        basis_list ::AbstractVector{BR}) where {BR<:Unsigned}
  if !issorted(basis_list)
    basis_list = sort(basis_list)
  end
  basis_lookup = Dict{BR, Int}(basis => ibasis for (ibasis, basis) in enumerate(basis_list))
  return HilbertSpaceRepresentation(basespace(hs), basis_list, basis_lookup)
end
