export symmetry_reduce, symmetry_unreduce


raw"""
    symmetry_unreduce

```math
\begin{pmatrix} l_1 \\ l_2 \\ l_3 \\ \vdots \\ l_n \end{pmatrix}
=
\begin{pmatrix}
. & \cdots & . \\
. & \cdots & . \\
. & \cdots & . \\
  & \dots & \\
. & \cdots &
\end{pmatrix}
\begin{pmatrix} s_1 \\ \vdots \\ s_m \end{pmatrix}
```
"""
function symmetry_unreduce(
    rhsr::ReducedHilbertSpaceRepresentation{HSR, BR, C},
    small_vector::AbstractVector{Si}
  ) where {HSR, BR, C, Si<:Number}
  if length(small_vector) != dimension(rhsr)
    throw(DimensionMismatch("Dimension of the input vector should match the reduced representation"))
  end
  So = promote_type(C, Si)
  large_vector = zeros(So, dimension(rhsr.parent))
  for (i_p, i_r) in enumerate(rhsr.basis_mapping_index)
    if i_r > 0
      ampl = rhsr.basis_mapping_amplitude[i_p]
      # H_r = U† H U
      large_vector[i_p] += ampl * small_vector[i_r]
    end
  end
  return large_vector
end


"""
    symmetry_reduce(rhsr, large_vector)

Reduce a large vector into the reduced hilbert space representation.
Simply throw away components that don't fit.
"""
function symmetry_reduce(
    rhsr::ReducedHilbertSpaceRepresentation{HSR, BR, C},
    large_vector::AbstractVector{Si}
  ) where {HSR, BR, C, Si<:Number}
  if length(large_vector) != dimension(rhsr.parent)
    throw(DimensionMismatch("Dimension of the input vector should match the larger representation"))
  end
  So = promote_type(C, Si)
  small_vector = zeros(So, dimension(rhsr))

  # basis mapping
  # (i_p | i_r | ampl) indicates : U_(p, r) = ampl
  for (i_p, i_r) in enumerate(rhsr.basis_mapping_index)
    if i_r > 0
      ampl = rhsr.basis_mapping_amplitude[i_p]
      # H_r = U† H U
      small_vector[i_r] += conj(ampl) * large_vector[i_p]
    end
  end
  return small_vector
end
