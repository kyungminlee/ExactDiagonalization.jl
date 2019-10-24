export ReducedSparseStateIndexed
# export get_component
# export set_component!

struct ReducedSparseStateIndexed
  reduced_hilbert_space_realization ::ReducedHilbertSpaceRepresentation
  components ::DefaultDict{Int, ComplexF64, ComplexF64}
  function ReducedSparseStateIndexed(rhsr ::ReducedHilbertSpaceRepresentation)
    return new(rhsr, DefaultDict{Int, ComplexF64}(zero(ComplexF64)))
  end
end

#=
|R(001)> = bᵢ(001) |001> + bᵢ(010) |010> + bᵢ(100) |100> + ...
(symmetry related)

|ψᵣ> = a(001) |R(001)> + a(011) |R(011)> + ...

# Return
  (a(R(bvec)) * b_i(bvec)
=#
import Base.getindex
function getindex(rss ::ReducedSparseStateIndexed, bvec::UInt) ::ComplexF64
  parent_lookup = rss.reduced_hilbert_space_realization.hilbert_space.basis_lookup
  reduced_lookup = rss.reduced_hilbert_space_realization.basis_lookup
  ivec_parent = parent_lookup[bvec]
  #lookup = rss.reduced_hilbert_space_realization.basis_lookup
  ivec, amplitude = reduced_lookup[ivec_parent]
  #(_, amplitude_parent) = reduced_lookup[ivec]
  return rss.components[ivec] * amplitude
end


#=
|R(001)> = bᵢ(001) |001> + bᵢ(010) |010> + bᵢ(100) |100> + ...
|ψᵣ> = a(001) |R(001)> + a(011) |R(011)> + ...

Set a(R(bvec)) = value / b(bvec)
=#
import Base.setindex!
function setindex!(rss ::ReducedSparseStateIndexed, value ::Number, bvec ::UInt)
  #list = rss.reduced_hilbert_space_realization.basis_list
  parent_lookup = rss.reduced_hilbert_space_realization.hilbert_space.basis_lookup
  reduced_lookup = rss.reduced_hilbert_space_realization.basis_lookup
  ivec_parent = parent_lookup[bvec]
  ivec, amplitude = reduced_lookup[ivec_parent]
  #(_, amplitude_parent) = reduced_lookup[ivec]
  rss.components[ivec] = value / amplitude
  return rss
end


import LinearAlgebra.norm
import LinearAlgebra.normalize!

function norm(rss ::ReducedSparseStateIndexed)
  list = rss.reduced_hilbert_space_realization.basis_list
  lookup = rss.reduced_hilbert_space_realization.basis_lookup

  Ng = length(rss.reduced_hilbert_space_realization.translation_group.elements)
  norm_sq = sum(abs(v)^2 for v in values(rss.components))
  return sqrt(norm_sq)
end

function normalize!(rss ::ReducedSparseStateIndexed)
  norm_value = norm(rss)
  for (ivec, ampl) in rss.components
    rss.components[ivec] = ampl / norm_value
  end
  return rss
end
