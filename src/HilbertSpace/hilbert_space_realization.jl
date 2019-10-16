export HilbertSpaceRealization
export dimension, realize, materialize

using MinimalPerfectHash

struct HilbertSpaceRealization{QN, BR}
  hilbert_space :: HilbertSpace{QN}
  basis_list ::Vector{BR}
  basis_lookup ::MinimalPerfectHash.CHD{BR, Int}
end

import Base.==
function (==)(lhs ::HilbertSpaceRealization{H1, B1}, rhs ::HilbertSpaceRealization{H2, B2}) where {H1, B1, H2, B2}
  return (H1 == H2) && (B1 == B2) && (lhs.hilbert_space == rhs.hilbert_space) && (lhs.basis_list == rhs.basis_list)
end

"""
    dimension

Dimension of the Concrete Hilbert space, i.e. number of basis vectors.
"""
dimension(hsr ::HilbertSpaceRealization) = length(hsr.basis_list)

"""
    realize(hs; BR ::DataType=UInt)

Make a HilbertSpaceRealization with all the basis vectors of the specified HilbertSpace.

# Arguments
- `hs ::HilbertSpace{QN}`: Abstract Hilbert space
- `BR ::DataType=UInt`: Binary representation type
"""
function realize(hs ::HilbertSpace{QN}; BR ::DataType=UInt) where {QN}
  basis_list = BR[]
  for indexarray in Iterators.product((1:length(site.states) for site in hs.sites)...)
    indexarray = Int[indexarray...]
    push!(basis_list, compress(hs, indexarray))
  end
  #basis_lookup = Dict{BR, Int}(basis => ibasis for (ibasis, basis) in enumerate(basis_list))
  basis_lookup = MinimalPerfectHash.CHD(basis_list, collect(1:length(basis_list)))
  return HilbertSpaceRealization{QN, BR}(hs, basis_list, basis_lookup)
end

# function concretize_naive(
#     hs ::HilbertSpace{QN},
#     qn ::QN;
#     BR ::DataType=UInt) where {QN}
#   sectors = quantum_number_sectors(hs)
#   if ! (qn in sectors)
#     return HilbertSpaceRealization{QN, BR}(hs, [], Dict())
#   end
#   basis_list = BR[]
#   for indexarray in Iterators.product((1:length(site.states) for site in hs.sites)...)
#     indexarray = Int[indexarray...]
#     q = get_quantum_number(hs, indexarray)
#     if q == qn
#       push!(basis_list, compress(hs, indexarray))
#     end
#   end
#   basis_lookup = Dict{BR, Int}()
#   for (ibasis, basis) in enumerate(basis_list)
#     basis_lookup[basis] = ibasis
#   end
#   return HilbertSpaceRealization{QN, BR}(hs, basis_list, basis_lookup)
# end

function realize(
    hs::HilbertSpace{QN},
    qn::QN;
    BR::DataType=UInt) where {QN}
  return realize(hs, [qn]; BR=BR)
end

"""
    realize(hs; BR ::DataType=UInt)

Make a HilbertSpaceRealization with all the basis vectors of the specified HilbertSpace.

# Arguments
- `hs ::HilbertSpace{QN}`: Abstract Hilbert space
- `allowed`: Allowed quantum numbers
- `BR ::DataType=UInt`: Binary representation type
"""
function realize(
    hs::HilbertSpace{QN},
    allowed::Union{AbstractSet{QN}, AbstractVector{QN}};
    BR::DataType=UInt) where {QN}
  allowed = Set(allowed)
  sectors = Set(quantum_number_sectors(hs))
  if isempty(intersect(allowed, sectors))
    return HilbertSpaceRealization{QN, BR}(hs, [], MinimalPerfectHash.CHD(BR[], Int[]))
  end

  quantum_numbers = [[state.quantum_number for state in site.states] for site in hs.sites]
  possible_quantum_numbers = [Set([zero(QN)])]  # PQN[i]: possible QN left of i

  n_sites = length(hs.sites)
  for i in 1:n_sites
    # pq = Set{QN}()
    # for q1 in possible_quantum_numbers[i], q2 in quantum_numbers[i]
    #   push!(pq, q1 .+ q2)
    # end
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

    result = DefaultDict{QN, Vector{BR}}(Vector{BR})
    for (i_state, q_curr) in enumerate(quantum_numbers[i])
      for (q_prev, states_prev) in result_prev
        q = q_prev + q_curr
        if q in allowed
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

  #basis_lookup = Dict{BR, Int}(basis => ibasis for (ibasis, basis) in enumerate(basis_list))
  basis_lookup = MinimalPerfectHash.CHD(basis_list, collect(1:length(basis_list)))
  return HilbertSpaceRealization{QN, BR}(hs, basis_list, basis_lookup)
end


function realize(hs ::HilbertSpace{QN},
                    basis_list ::AbstractArray{BR}) where {QN, BR<:Unsigned}
  #basis_lookup = Dict{BR, Int}(basis => ibasis for (ibasis, basis) in enumerate(basis_list))
  basis_lookup = MinimalPerfectHash.CHD(basis_list, collect(1:length(basis_list)))
  return HilbertSpaceRealization{QN, BR}(hs, basis_list, basis_lookup)
end
