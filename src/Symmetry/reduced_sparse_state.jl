export ReducedSparseState
export get_component
export set_component!

struct ReducedSparseState
  reduced_hilbert_space_realization ::ReducedHilbertSpaceRealization
  components ::DefaultDict{Int, ComplexF64, ComplexF64}
  function ReducedSparseState(rhsr ::ReducedHilbertSpaceRealization)
    return new(rhsr, DefaultDict{Int, ComplexF64}(zero(ComplexF64)))
  end
end

"""
|R(001)> = bᵢ(001) |001> + bᵢ(010) |010> + bᵢ(100) |100> + ...
(symmetry related)

|ψᵣ> = a(001) |R(001)> + a(011) |R(011)> + ...

# Return
  (a(R(bvec)) * b_i(bvec)
"""
function get_component(rss ::ReducedSparseState, bvec::UInt) ::ComplexF64
  lookup = rss.reduced_hilbert_space_realization.basis_lookup
  ivec, amplitude = lookup[bvec]
  #(_, amplitude_parent) = lookup[ivec]
  return rss.components[ivec] * amplitude
end


"""
|R(001)> = bᵢ(001) |001> + bᵢ(010) |010> + bᵢ(100) |100> + ...
|ψᵣ> = a(001) |R(001)> + a(011) |R(011)> + ...

Set a(R(bvec)) = value / b(bvec)
"""
function set_component!(rss ::ReducedSparseState, value ::ComplexF64, bvec ::UInt)
  list = rss.reduced_hilbert_space_realization.basis_list
  lookup = rss.reduced_hilbert_space_realization.basis_lookup
  ivec, amplitude = lookup[bvec]
  (_, amplitude_parent) = lookup[list[ivec]]
  rss.components[ivec] = value / amplitude
  return rss
end


# import Base.setindex!
# function setindex!(rss ::ReducedSparseState, value ::ComplexF64, bvec ::UInt)
#   list = rss.reduced_hilbert_space_realization.basis_list
#   lookup = rss.reduced_hilbert_space_realization.basis_lookup
#   ivec, amplitude = lookup[bvec]
#   (_, amplitude_parent) = lookup[list[ivec]]
#   rss.components[ivec] = value / amplitude * amplitude_parent
#   return rss
# end

# import Base.getindex
# function getindex(rss::ReducedSparseState, bvec::UInt)
#   lookup = rss.reduced_hilbert_space_realization.basis_lookup
#   ivec, amplitude = lookup[bvec]
#   (_, amplitude_parent) = lookup[ivec]
#   return rss.components[ivec] * amplitude / amplitude_parent
# end

import LinearAlgebra.norm
import LinearAlgebra.normalize!

function norm(rss ::ReducedSparseState)
  list = rss.reduced_hilbert_space_realization.basis_list
  lookup = rss.reduced_hilbert_space_realization.basis_lookup

  Ng = length(rss.reduced_hilbert_space_realization.translation_group.elements)
  norm_sq = sum(abs(v)^2 for v in values(rss.components))
  return sqrt(norm_sq)
end

function normalize!(rss ::ReducedSparseState)
  norm_value = norm(rss)
  for (ivec, ampl) in rss.components
    rss.components[ivec] = ampl / norm_value
  end
  return rss
end
