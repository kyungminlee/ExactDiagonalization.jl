export symmetry_reduce, symmetry_reduce_serial, symmetry_reduce_parallel
export symmetry_reduce, symmetry_unreduce

import TightBindingLattice.TranslationGroup


"""
    symmetry_reduce(hsr, trans_group, frac_momentum, complex_type=ComplexF64, tol=sqrt(eps(Float64)))

Symmetry-reduce the HilbertSpaceRepresentation using translation group.

"""
function symmetry_reduce(
    hsr::HilbertSpaceRepresentation{QN, BR, DT},
    trans_group::TranslationGroup,
    fractional_momentum::AbstractVector{<:Rational},
    complex_type::Type{ComplexType}=ComplexF64;
    tol::Real=sqrt(eps(Float64))) where {QN, BR, DT, ComplexType<:Complex}
  symred = Threads.nthreads() == 1 ? symmetry_reduce_serial : symmetry_reduce_parallel
  return symred(hsr, trans_group, fractional_momentum, ComplexType; tol=tol)
end


"""
    symmetry_reduce_serial(hsr, trans_group, frac_momentum, complex_type=ComplexF64, tol=sqrt(eps(Float64)))

Symmetry-reduce the HilbertSpaceRepresentation using translation group (single threaded).

"""
function symmetry_reduce_serial(
    hsr::HilbertSpaceRepresentation{QN, BR, DT},
    trans_group::TranslationGroup,
    fractional_momentum::AbstractVector{<:Rational},
    complex_type::Type{ComplexType}=ComplexF64;
    tol::FloatType=sqrt(eps(Float64))
    ) where {QN, BR, DT, ComplexType<:Complex, FloatType<:Real}

  HSR = HilbertSpaceRepresentation{QN, BR, DT}

  ik = let
    match(k::Vector{Rational{Int}})::Bool = k == fractional_momentum
    findfirst(match, trans_group.fractional_momenta)
  end

  if ik === nothing
    throw(ArgumentError("fractional momentum $(fractional_momentum)" *
                        " not an irrep of the translation group"))
  end

  n_basis = length(hsr.basis_list)

  basis_mapping_representative = Vector{Int}(undef, n_basis)
  fill!(basis_mapping_representative, -1)
  basis_mapping_amplitude = zeros(ComplexType, n_basis)

  phases = conj.(trans_group.character_table[ik, :])
  group_size = length(trans_group.elements)

  size_estimate = let
    denom = max(1, length(trans_group.fractional_momenta) - 1)
    n_basis ÷ denom
  end

  reduced_basis_list = BR[]
  sizehint!(reduced_basis_list, size_estimate)

  visited = falses(n_basis)

  basis_states = Vector{BR}(undef, group_size)
  basis_amplitudes = Dict{BR, ComplexType}()
  sizehint!(basis_amplitudes, group_size + group_size ÷ 2)

  is_identity = [is_compatible(fractional_momentum, t) for t in trans_group.translations]

  for ivec_p in 1:n_basis
    visited[ivec_p] && continue
    bvec = hsr.basis_list[ivec_p]

    compatible = true
    for i in 2:group_size
      g = trans_group.elements[i]
      bvec_prime = symmetry_apply(hsr.hilbert_space, g, bvec)
      if bvec_prime < bvec
        compatible = false
        break
      elseif bvec_prime == bvec && !is_identity[i]
        compatible = false
        break
      end
      basis_states[i] = bvec_prime
    end
    (!compatible) && continue
    basis_states[1] = bvec

    push!(reduced_basis_list, bvec)

    empty!(basis_amplitudes)
    for i in 1:group_size
      p = phases[i]
      bvec_prime = basis_states[i]
      basis_amplitudes[bvec_prime] = p # They're all the same.
    end
    inv_norm = 1.0 / sqrt(float(length(basis_amplitudes)))

    for (bvec_prime, amplitude) in basis_amplitudes
      ivec_p_prime = hsr.basis_lookup[bvec_prime]
      visited[ivec_p_prime] = true
      basis_mapping_representative[ivec_p_prime] = ivec_p
      basis_mapping_amplitude[ivec_p_prime] = amplitude * inv_norm
    end
  end

  basis_mapping_index = Vector{Int}(undef, n_basis)
  fill!(basis_mapping_index, -1)

  for (ivec_r, bvec) in enumerate(reduced_basis_list)
    ivec_p = hsr.basis_lookup[bvec]
    basis_mapping_index[ivec_p] = ivec_r
  end

  for (ivec_p_prime, ivec_p) in enumerate(basis_mapping_representative)
    (ivec_p <= 0) && continue  # not in this irrep
    (ivec_p_prime == ivec_p) && continue  # already in the lookup
    ivec_r = basis_mapping_index[ivec_p]
    basis_mapping_index[ivec_p_prime] = ivec_r
  end

  return ReducedHilbertSpaceRepresentation{HSR, BR, ComplexType}(hsr, trans_group, reduced_basis_list,
                                                                 basis_mapping_index, basis_mapping_amplitude)
end


"""
    symmetry_reduce_parallel(hsr, trans_group, frac_momentum, complex_type=ComplexF64, tol=sqrt(eps(Float64)))

Symmetry-reduce the HilbertSpaceRepresentation using translation group (multi-threaded).

"""
function symmetry_reduce_parallel(
    hsr::HilbertSpaceRepresentation{QN, BR, DT},
    trans_group::TranslationGroup,
    fractional_momentum::AbstractVector{<:Rational},
    complex_type::Type{ComplexType}=ComplexF64;
    tol::Real=sqrt(eps(Float64))
    ) where {QN, BR, DT, ComplexType<:Complex}

  HSR = HilbertSpaceRepresentation{QN, BR, DT}

  @debug "BEGIN symmetry_reduce_parallel"
  ik = let
    match(k::Vector{Rational{Int}})::Bool  = k == fractional_momentum
    findfirst(match, trans_group.fractional_momenta)
  end

  if ik === nothing
    throw(ArgumentError("fractional momentum $(fractional_momentum)" *
                        " not an irrep of the translation group"))
  end

  n_basis = length(hsr.basis_list)
  @debug "Original Hilbert space dimension: $n_basis"

  # basis_mapping_index and basis_mapping_amplitude contain information about
  # which basis vector of the larger Hilbert space is included
  # in the basis of the smaller Hilbert space with what amplitude.
  basis_mapping_representative = Vector{Int}(undef, n_basis)
  fill!(basis_mapping_representative, -1)
  basis_mapping_amplitude = zeros(ComplexType, n_basis)

  nthreads = Threads.nthreads()
  size_estimate = n_basis ÷ max(1, length(trans_group.fractional_momenta) - 1)
  @debug "Estimate for the reduced Hilbert space dimension: $size_estimate"

  local_reduced_basis_list = Vector{Vector{BR}}(undef, nthreads)
  for i in 1:nthreads
    local_reduced_basis_list[i] = BR[]
    sizehint!(local_reduced_basis_list[i], size_estimate ÷ nthreads + 1)
  end

  # |σᵨ⟩ = (1/Norm) ∑ χᵨ(g)* Rᵨ(g) |σ⟩
  phases = conj.(trans_group.character_table[ik, :])
  group_size = length(trans_group.elements)

  visited = zeros(UInt8, n_basis) # use UInt8 rather than Bool for thread safety

  local_basis_states = Matrix{BR}(undef, (nthreads, group_size))
  local_basis_amplitudes = Vector{Dict{BR, ComplexType}}(undef, nthreads)
  for id in 1:nthreads
    local_basis_amplitudes[id] = Dict{BR, ComplexType}()
    sizehint!(local_basis_amplitudes[id], group_size)
  end

  # Load balancing
  #   The representatives are the smaller binary numbers.
  #   Distribute them equally among threads.
  reorder = Int[]
  sizehint!(reorder, n_basis)
  nblocks = (n_basis + nthreads - 1) ÷ nthreads
  for i in 1:nthreads, j in 1:nblocks
    k = i + nthreads * (j-1)
    if 1 <= k <= n_basis
      push!(reorder, k)
    end
  end

  # whether the element is in the kernel of the representation
  is_identity = [is_compatible(fractional_momentum, t) for t in trans_group.translations]

  @debug "Starting reduction (parallel)"
  Threads.@threads for itemp in 1:n_basis
    ivec_p = reorder[itemp]
    (visited[ivec_p] != 0x0) && continue

    id = Threads.threadid()
    bvec = hsr.basis_list[ivec_p]

    # A basis binary representation is incompatible with the reduced Hilbert space if
    # (1) it is not the smallest among the its star, or
    # (2) its star is smaller than the representation
    compatible = true
    for i in 2:group_size
      g = trans_group.elements[i]
      bvec_prime = symmetry_apply(hsr.hilbert_space, g, bvec)
      if bvec_prime < bvec
        compatible = false
        break
      elseif bvec_prime == bvec && !is_identity[i]
        compatible = false
        break
      end
      local_basis_states[id, i] = bvec_prime
    end
    (!compatible) && continue
    local_basis_states[id, 1] = bvec

    push!(local_reduced_basis_list[id], bvec)

    empty!(local_basis_amplitudes[id])
    for i in 1:group_size
      p = phases[i]
      bvec_prime = local_basis_states[id, i]
      local_basis_amplitudes[id][bvec_prime] = p # Same bvec_prime, same p.
    end
    inv_norm = 1.0 / sqrt(float(length(local_basis_amplitudes[id])))

    for (bvec_prime, amplitude) in local_basis_amplitudes[id]
      ivec_p_prime = hsr.basis_lookup[bvec_prime]
      visited[ivec_p_prime] = 0x1
      basis_mapping_representative[ivec_p_prime] = ivec_p
      basis_mapping_amplitude[ivec_p_prime] = amplitude*inv_norm
    end
  end
  @debug "Finished reduction (parallel)"

  reduced_basis_list = BR[]
  sizehint!(reduced_basis_list, sum(length(x) for x in local_reduced_basis_list))
  while !isempty(local_reduced_basis_list)
    lbl = pop!(local_reduced_basis_list)
    append!(reduced_basis_list, lbl)
  end
  @debug "Collected basis list"

  sort!(reduced_basis_list)
  @debug "Sorted basis list"

  # Basis vectors of the unreduced Hilbert space that are
  # not included in the reduced Hilbert space are marked as -1
  basis_mapping_index = Vector{Int}(undef, n_basis)
  fill!(basis_mapping_index, -1)

  Threads.@threads for ivec_r in eachindex(reduced_basis_list)
    bvec = reduced_basis_list[ivec_r]
    ivec_p = hsr.basis_lookup[bvec]
    basis_mapping_index[ivec_p] = ivec_r
  end
  @debug "Collected basis lookup (diagonal)"

  Threads.@threads for ivec_p_prime in eachindex(basis_mapping_representative)
    ivec_p = basis_mapping_representative[ivec_p_prime]
    (ivec_p <= 0) && continue  # not in this irrep
    (ivec_p_prime == ivec_p) && continue  # already in the lookup
    ivec_r = basis_mapping_index[ivec_p]
    basis_mapping_index[ivec_p_prime] = ivec_r
  end
  @debug "Collected basis lookup (offdiagonal)"

  @debug "END symmetry_reduce_parallel"
  return ReducedHilbertSpaceRepresentation{HSR, BR, ComplexType}(hsr, trans_group, reduced_basis_list,
                                                                 basis_mapping_index, basis_mapping_amplitude)
end


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
