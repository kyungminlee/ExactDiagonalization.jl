export ReducedHilbertSpaceRealization
export symmetry_reduce, symmetry_reduce_parallel
export materialize, materialize_parallel

struct ReducedHilbertSpaceRealization{QN, BR, C<:Complex}
  parent_hilbert_space_realization ::HilbertSpaceRealization{QN, BR}
  translation_group ::TranslationGroup
  basis_list ::Vector{BR}
  basis_lookup ::MinimalPerfectHash.CHD{BR, NamedTuple{(:index, :amplitude), Tuple{Int, C}}}
  #basis_lookup ::Dict{BR, NamedTuple{(:index, :amplitude), Tuple{Int, C}}}
end

function symmetry_reduce(hsr ::HilbertSpaceRealization{QN, BR},
                trans_group ::TranslationGroup,
                fractional_momentum ::AbstractVector{Rational};
                ComplexType::DataType=ComplexF64) where {QN, BR}
  ik = findfirst(collect(
    trans_group.fractional_momenta[ik] == fractional_momentum
    for ik in 1:length(trans_group.fractional_momenta) ))

  isnothing(ik) && throw(ArgumentError("fractional momentum $(fractional_momentum) not an irrep of the translation group"))

  phases = trans_group.character_table[ik, :]
  n_basis = length(hsr.basis_list)

  reduced_basis_list = BR[]
  ParentAmplitudeType = NamedTuple{(:parent,:amplitude), Tuple{Int,ComplexType}}
  parent_amplitude_list = ParentAmplitudeType[(parent=-1,amplitude=zero(ComplexType)) for i in 1:n_basis]

  size_estimate = let
    denom = max(1, length(trans_group.fractional_momenta) - 1)
    n_basis ÷ denom
  end
  sizehint!(reduced_basis_list, size_estimate)

  visited = falses(n_basis)
  for ivec in 1:n_basis
    visited[ivec] && continue
    bvec = hsr.basis_list[ivec]

    compatible = true
    ψ = SparseState{ComplexType, BR}(hsr.hilbert_space)
    for i in 1:length(trans_group.elements)
      t = trans_group.translations[i]
      g = trans_group.elements[i]

      bvec_prime = apply_symmetry(hsr.hilbert_space, g, bvec)

      if bvec_prime < bvec
        compatible = false
        break
      elseif bvec_prime == bvec && !is_compatible(fractional_momentum, t)
        compatible = false
        break
      end

      p = phases[i]
      ψ[bvec_prime] += p
    end
    (!compatible) && continue

    clean!(ψ)
    @assert !isempty(ψ)

    ivec_primes = Int[hsr.basis_lookup[bvec_prime] for bvec_prime in keys(ψ.components)]

    if any(visited[ivec_primes])
      continue
    else
      visited[ivec_primes] .= true
    end

    normalize!(ψ)
    push!(reduced_basis_list, bvec)
    for (bvec_prime, amplitude) in ψ.components
      ivec_prime = hsr.basis_lookup[bvec_prime]
      parent_amplitude_list[ivec_prime] = (parent=ivec, amplitude=amplitude)
    end
  end

  sort!(reduced_basis_list)

  ItemType = NamedTuple{(:index, :amplitude), Tuple{Int, ComplexType}}
  reduced_basis_lookup = Dict{BR, ItemType}(
                              let
                                ivec = hsr.basis_lookup[bvec]
                                amplitude = parent_amplitude_list[ivec].amplitude
                                bvec => (index=ivec_r, amplitude=amplitude)
                              end for (ivec_r, bvec) in enumerate(reduced_basis_list))
  sizehint!(reduced_basis_lookup, length(reduced_basis_list))

  for (ivec_prime, (ivec, amplitude)) in enumerate(parent_amplitude_list)
    (ivec == -1) && continue  # not in this irrep
    (ivec_prime == ivec) && continue # already in the lookup

    bvec_prime = hsr.basis_list[ivec_prime]
    bvec = hsr.basis_list[ivec]
    ivec_r = reduced_basis_lookup[bvec].index
    reduced_basis_lookup[bvec_prime] = (index=ivec_r, amplitude=amplitude)
  end
  reduced_basis_lookup = MinimalPerfectHash.CHD{BR, NamedTuple{(:index, :amplitude), Tuple{Int, C}}}(reduced_basis_lookup)

  return ReducedHilbertSpaceRealization{QN, BR, ComplexType}(hsr, trans_group, reduced_basis_list, reduced_basis_lookup)
end


function symmetry_reduce_parallel(hsr ::HilbertSpaceRealization{QN, BR},
                trans_group ::TranslationGroup,
                fractional_momentum ::AbstractVector{Rational};
                ComplexType::DataType=ComplexF64) where {QN, BR}
  debug(LOGGER, "BEGIN symmetry_reduce_parallel")
  ik = findfirst(collect(
    trans_group.fractional_momenta[ik] == fractional_momentum
    for ik in 1:length(trans_group.fractional_momenta) ))

  isnothing(ik) && throw(ArgumentError("fractional momentum $(fractional_momentum) not an irrep of the translation group"))
  phases = trans_group.character_table[ik, :]
  n_basis = length(hsr.basis_list)
  debug(LOGGER, "Original Hilbert space dimension: $n_basis")

  visit_lock = Threads.SpinLock()
  nthreads = Threads.nthreads()
  local_reduced_basis_list = [BR[] for i in 1:nthreads]
  ParentAmplitudeType = NamedTuple{(:parent,:amplitude), Tuple{Int,ComplexType}}
  parent_amplitude_list = ParentAmplitudeType[(parent=-1,amplitude=zero(ComplexType)) for i in 1:n_basis]

  size_estimate = let
    denom = max(1, length(trans_group.fractional_momenta) - 1)
    n_basis ÷ denom
  end
  debug(LOGGER, "Estimate for the reduced Hilbert space dimension (local): $size_estimate")
  for i in eachindex(local_reduced_basis_list)
    sizehint!(local_reduced_basis_list[i], size_estimate ÷ nthreads + 1)
  end

  # Load balancing (the representatives are the smaller binary numbers)
  reorder = Int[]
  sizehint!(reorder, n_basis)
  nblocks = (n_basis + nthreads - 1) ÷ nthreads
  for i in 1:nthreads, j in 1:nblocks
    k = i + nthreads * (j-1)
    if 1 <= k <= n_basis
      push!(reorder, k)
    end
  end
  #@assert length(reorder) == n_basis
  @assert sort(reorder) == 1:n_basis

  visited = falses(n_basis)
  debug(LOGGER, "Starting reduction (parallel)")
  Threads.@threads for itemp in 1:n_basis
    ivec = reorder[itemp]
    visited[ivec] && continue
    id = Threads.threadid()
    bvec = hsr.basis_list[ivec]

    compatible = true
    ψ = SparseState{ComplexType, BR}(hsr.hilbert_space)
    for i in 1:length(trans_group.elements)
      t = trans_group.translations[i]
      g = trans_group.elements[i]

      bvec_prime = apply_symmetry(hsr.hilbert_space, g, bvec)

      if bvec_prime < bvec
        compatible = false
        break
      elseif bvec_prime == bvec && !is_compatible(fractional_momentum, t)
        compatible = false
        break
      end

      p = phases[i]
      ψ[bvec_prime] += p
    end
    (!compatible) && continue

    clean!(ψ)
    @assert !isempty(ψ)

    ivec_primes = Int[hsr.basis_lookup[bvec_prime] for bvec_prime in keys(ψ.components)]

    lock(visit_lock)
    if any(visited[ivec_primes])
      unlock(visit_lock)
      continue
    else
      visited[ivec_primes] .= true
      unlock(visit_lock)
    end

    normalize!(ψ)
    push!(local_reduced_basis_list[id], bvec)
    for (bvec_prime, amplitude) in ψ.components
      ivec_prime = hsr.basis_lookup[bvec_prime]
      parent_amplitude_list[ivec_prime] = (parent=ivec, amplitude=amplitude)
    end
  end
  debug(LOGGER, "Finished reduction (parallel)")

  debug(LOGGER, "Collecting basis list")
  reduced_basis_list = BR[]
  while !isempty(local_reduced_basis_list)
    lbl = pop!(local_reduced_basis_list)
    append!(reduced_basis_list, lbl)
  end
  sort!(reduced_basis_list)

  debug(LOGGER, "Collecting basis lookup (diagonal)")
  ItemType = NamedTuple{(:index, :amplitude), Tuple{Int, ComplexType}}
  reduced_basis_lookup = Dict{BR, ItemType}(
                              let
                                ivec = hsr.basis_lookup[bvec]
                                amplitude = parent_amplitude_list[ivec].amplitude
                                bvec => (index=ivec_r, amplitude=amplitude)
                              end for (ivec_r, bvec) in enumerate(reduced_basis_list))
  sizehint!(reduced_basis_lookup, length(reduced_basis_list))

  debug(LOGGER, "Collecting basis lookup (offdiagonal)")
  for (ivec_prime, (ivec, amplitude)) in enumerate(parent_amplitude_list)
    (ivec == -1) && continue  # not in this irrep
    (ivec_prime == ivec) && continue # already in the lookup

    bvec_prime = hsr.basis_list[ivec_prime]
    bvec = hsr.basis_list[ivec]
    ivec_r = reduced_basis_lookup[bvec].index
    reduced_basis_lookup[bvec_prime] = (index=ivec_r, amplitude=amplitude)
  end
  debug(LOGGER, "END symmetry_reduce_parallel")
  reduced_basis_lookup = MinimalPerfectHash.CHD{BR, NamedTuple{(:index, :amplitude), Tuple{Int, C}}}(reduced_basis_lookup)
  return ReducedHilbertSpaceRealization{QN, BR, ComplexType}(hsr, trans_group, reduced_basis_list, reduced_basis_lookup)
end


function materialize(rhsr :: ReducedHilbertSpaceRealization{QN, BR, C},
                     operator ::AbstractOperator;
                     tol::Real=sqrt(eps(Float64))) where {QN, BR, C}

  @assert is_invariant(rhsr.translation_group, operator)
  hs = rhsr.parent_hilbert_space_realization.hilbert_space

  rows = Int[]
  cols = Int[]
  vals = C[]
  err  = 0.0

  n_basis = length(rhsr.basis_list)
  for irow in 1:n_basis
    brow = rhsr.basis_list[irow]
    ampl_row = rhsr.basis_lookup[brow].amplitude
    ψrow = SparseState{C, BR}(hs, brow=>1/ampl_row)
    ψcol = SparseState{C, BR}(hs)
    apply_unsafe!(ψcol, ψrow, operator)
    clean!(ψcol)

    for (bcol, ampl) in ψcol.components
      if ! haskey(rhsr.basis_lookup, bcol)
        err += abs(ampl.^2)
        continue
      end

      (icol, ampl_col) = rhsr.basis_lookup[bcol]
      push!(rows, irow)
      push!(cols, icol)
      push!(vals, ampl * ampl_col)
    end
  end

  if isempty(vals)
    vals = Float64[]
  elseif isapprox(maximum(abs.(imag.(vals))), 0; atol=tol)
    vals = real.(vals)
  end
  spmat = sparse(rows, cols, vals, n_basis, n_basis)
  droptol!(spmat, tol)
  return (spmat, err)
end



function materialize_parallel(rhsr :: ReducedHilbertSpaceRealization{QN, BR, C},
                     operator ::AbstractOperator;
                     tol::Real=sqrt(eps(Float64))) where {QN, BR, C}
  debug(LOGGER, "BEGIN materialize_parallel for ReducedHilbertSpaceRealiation")
  debug(LOGGER, "Checking whether the operator is translationally invariant")
  @assert is_invariant(rhsr.translation_group, operator)
  hs = rhsr.parent_hilbert_space_realization.hilbert_space

  nthreads = Threads.nthreads()
  debug(LOGGER, "Number of threads: $nthreads")
  local_rows = [ Int[] for i in 1:nthreads]
  local_cols = [ Int[] for i in 1:nthreads]
  local_vals = [ C[] for i in 1:nthreads]
  local_err =  Float64[0.0 for i in 1:nthreads]

  n_basis = length(rhsr.basis_list)
  debug(LOGGER, "Starting materialization (parallel)")
  Threads.@threads for irow in 1:n_basis
    id = Threads.threadid()
    brow = rhsr.basis_list[irow]

    ampl_row = rhsr.basis_lookup[brow].amplitude
    ψrow = SparseState{C, BR}(hs, brow=>1/ampl_row)
    ψcol = SparseState{C, BR}(hs)
    apply_unsafe!(ψcol, ψrow, operator)
    clean!(ψcol)

    for (bcol, ampl) in ψcol.components
      if ! haskey(rhsr.basis_lookup, bcol)
        local_err[id] += abs(ampl.^2)
        continue
      end

      (icol, ampl_col) = rhsr.basis_lookup[bcol]
      push!(local_rows[id], irow)
      push!(local_cols[id], icol)
      push!(local_vals[id], ampl * ampl_col)
    end
  end
  debug(LOGGER, "Finished materialization (parallel)")

  rows ::Vector{Int} = vcat(local_rows...)
  cols ::Vector{Int} = vcat(local_cols...)
  vals ::Vector{C} = vcat(local_vals...)
  err ::Float64 = sum(local_err)

  if isempty(vals)
    debug(LOGGER, "Matrix empty")
    vals = Float64[]
  elseif isapprox(maximum(abs.(imag.(vals))), 0; atol=tol)
    debug(LOGGER, "Matrix purely real")
    vals = real.(vals)
  end
  spmat = sparse(rows, cols, vals, n_basis, n_basis)
  debug(LOGGER, "Number of nonzero elements: $(length(spmat.nzval))")
  droptol!(spmat, tol)
  debug(LOGGER, "Number of nonzero elements after droptol!: $(length(spmat.nzval))")
  debug(LOGGER, "END materialize_parallel for ReducedHilbertSpaceRealiation")
  return (spmat, err)
end
