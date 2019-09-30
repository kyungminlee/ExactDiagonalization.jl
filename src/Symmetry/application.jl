
function apply_symmetry(hs::HilbertSpace, permutation ::Permutation, bitrep ::BR) where {BR}
  out = zero(BR)
  for (i, j) in enumerate(permutation.map)
    out |= ( ( (bitrep >> hs.bitoffsets[i]) & make_bitmask(hs.bitwidths[i]) ) << hs.bitoffsets[j] )
  end
  return out
end

