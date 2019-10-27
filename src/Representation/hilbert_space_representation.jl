export HilbertSpaceRepresentation
export dimension, represent, materialize

struct HilbertSpaceRepresentation{HS <:AbstractHilbertSpace, BR <:Unsigned} <:AbstractHilbertSpaceRepresentation
  hilbert_space ::HS
  basis_list ::Vector{BR}
  basis_lookup ::FrozenSortedArrayIndex{BR}

  function HilbertSpaceRepresentation(
      hilbert_space ::HilbertSpace{QN},
      basis_list::AbstractVector{BR},
      basis_lookup::FrozenSortedArrayIndex{BR}) where {QN, BR <:Unsigned}
    if count_ones(typemax(BR)) <= bitwidth(hilbert_space)
      # equality added such that the MSB checks overflow
      throw(ArgumentError("type $(BR) not enough to represent the hilbert space (need $(bitwidth(hilbert_space)) bits)"))
    end
    return new{HilbertSpace{QN}, BR}(hilbert_space, basis_list, basis_lookup)
  end

  function HilbertSpaceRepresentation(
      hilbert_space_sector ::HilbertSpaceSector{QN},
      basis_list::AbstractVector{BR},
      basis_lookup::FrozenSortedArrayIndex{BR}) where {QN, BR <:Unsigned}
    if count_ones(typemax(BR)) <= bitwidth(hilbert_space_sector)
      # equality added such that the MSB checks overflow
      throw(ArgumentError("type $(BR) not enough to represent the hilbert space (need $(bitwidth(hilbert_space_sector)) bits)"))
    end
    return new{HilbertSpace{QN}, BR}(hilbert_space_sector.parent, basis_list, basis_lookup)
  end
end

@inline scalartype(lhs ::Type{HilbertSpaceRepresentation{HS, BR}}) where {HS, BR} = Bool
@inline bintype(lhs ::Type{HilbertSpaceRepresentation{HS, BR}}) where {HS, BR} = BR ::DataType

@inline basespace(lhs::HilbertSpaceRepresentation{HS, BR}) where {HS<:AbstractHilbertSpace, BR<:UInt} = lhs.hilbert_space ::HS

@inline bitwidth(hss::HilbertSpaceSector{QN}) where QN = bitwidth(hss.parent) ::Int


import Base.==
function (==)(lhs ::HilbertSpaceRepresentation{H1, B1},
              rhs ::HilbertSpaceRepresentation{H2, B2}) where {H1, B1, H2, B2}
  return basespace(lhs) == basespace(rhs) && (lhs.basis_list == rhs.basis_list)
end

function checkvalidbasis(hsr::HilbertSpaceRepresentation{HS, BR}) where {HS <:AbstractHilbertSpace, BR <:Unsigned}
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

"""
    represent(hs; BR ::DataType=UInt)

Make a HilbertSpaceRepresentation with all the basis vectors of the specified HilbertSpace.

# Arguments
- `hs ::HilbertSpace{QN}`: Abstract Hilbert space
- `BR ::DataType=UInt`: Binary representation type
"""
function represent(hs ::HilbertSpace{QN}; BR ::DataType=UInt) where {QN}
  if count_ones(typemax(BR)) <= bitwidth(hs)
    throw(ArgumentError("type $(BR) not enough to represent the hilbert space (need $(bitwidth(hs)) bits)"))
  end
  HS = HilbertSpace{QN}
  basis_list = BR[]
  for indexarray in keys(hs)
    push!(basis_list, compress(hs, indexarray))
  end
  sort!(basis_list)
  basis_lookup = FrozenSortedArrayIndex{BR}(basis_list)
  return HilbertSpaceRepresentation(hs, basis_list, basis_lookup)
end

"""
    represent(hs; BR ::DataType=UInt)

Make a HilbertSpaceRepresentation with all the basis vectors of the specified HilbertSpace.

# Arguments
- `hs ::HilbertSpace{QN}`: Abstract Hilbert space
- `allowed`: Allowed quantum numbers
- `BR ::DataType=UInt`: Binary representation type
"""
function represent(hss::HilbertSpaceSector{QN}; BR::DataType=UInt) where {QN}
  hs = hss.parent
  allowed = hss.allowed_quantum_numbers
  sectors = Set(quantum_number_sectors(hs))
  if isempty(intersect(allowed, sectors))
    return HilbertSpaceRepresentation(hs, BR[], FrozenSortedArrayIndex{BR}(BR[]))
  end

  quantum_numbers = [[state.quantum_number for state in site.states] for site in hs.sites]
  possible_quantum_numbers = [Set([zero(QN)])]  # PQN[i]: possible QN left of i

  n_sites = length(hs.sites)
  for i in 1:n_sites
    pq = Set{QN}(q1 + q2 for q1 in possible_quantum_numbers[i], q2 in quantum_numbers[i])
    push!(possible_quantum_numbers, pq)
  end

  function generate(i ::Int, allowed ::AbstractSet{QN})
    if i == 0
      return (zero(QN) in allowed) ? Dict(zero(QN) => [BR(0x0)]) : Dict()
    end
    allowed_prev = Set{QN}()
    for q1 in quantum_numbers[i], q2 in allowed
      q = q2 - q1
      if q in possible_quantum_numbers[i]
        push!(allowed_prev, q)
      end
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
  basis_lookup = FrozenSortedArrayIndex{BR}(basis_list)
  return HilbertSpaceRepresentation(hs, basis_list, basis_lookup)
end


function represent(hs ::AbstractHilbertSpace,
                   basis_list ::AbstractVector{BR}) where {BR<:Unsigned}
  #HS = HilbertSpace{QN}
  basis_list = sort(basis_list)
  #basis_lookup = Dict{BR, Int}(basis => ibasis for (ibasis, basis) in enumerate(basis_list))
  basis_lookup = FrozenSortedArrayIndex{BR}(basis_list)
  return HilbertSpaceRepresentation(basespace(hs), basis_list, basis_lookup)
end
