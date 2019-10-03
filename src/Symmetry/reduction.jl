export ReducedHilbertSpaceRealization
export symmetry_reduce, symmetry_reduce_parallel
export materialize, materialize_parallel

struct ReducedHilbertSpaceRealization{QN, BR, C<:Complex}
  parent_hilbert_space_realization ::HilbertSpaceRealization{QN, BR}
  translation_group ::TranslationGroup
  basis_list ::Vector{UInt}
  basis_lookup ::Dict{UInt, NamedTuple{(:index, :amplitude), Tuple{Int, C}}}
end

function symmetry_reduce(hsr ::HilbertSpaceRealization{QN, BR},
                trans_group ::TranslationGroup,
                fractional_momentum ::AbstractVector{Rational};
                ComplexType::DataType=ComplexF64) where {QN, BR}
  ik = findfirst(collect(
    trans_group.fractional_momenta[ik] == fractional_momentum
    for ik in 1:length(trans_group.fractional_momenta) ))
  
  isnothing(ik) && throw(ArgumentError("fractional momentum $(fractional_momentum) not an irrep of the translation group"))

  # check if fractional momentum is compatible with translation group ?
  #k = float.(fractional_momentum) .* 2π
  phases = trans_group.character_table[ik, :]
  #[ cis(dot(k, t)) for t in trans_group.translations]
  reduced_basis_list = Set{UInt}()
  parent_amplitude = Dict()

  for bvec in hsr.basis_list
    if haskey(parent_amplitude, bvec)
      continue
    end

    ψ = SparseState{ComplexType, UInt}(hsr.hilbert_space)
    identity_translations = Vector{Int}[]
    for i in 1:length(trans_group.elements)
      t = trans_group.translations[i]
      g = trans_group.elements[i]
      p = phases[i]

      bvec_prime = apply_symmetry(hsr.hilbert_space, g, bvec)
      ψ[bvec_prime] += p
      if bvec_prime == bvec
        push!(identity_translations, t)
      end
    end
    if !is_compatible(fractional_momentum, identity_translations)
      continue
    end
    clean!(ψ)
    @assert !isempty(ψ)

    normalize!(ψ)
    push!(reduced_basis_list, bvec)

    for (bvec_prime, amplitude) in ψ.components
      parent_amplitude[bvec_prime] = (parent=bvec, amplitude=amplitude)
    end
  end
  reduced_basis_list = sort(collect(reduced_basis_list))
  reduced_basis_lookup = Dict(bvec => (index=ivec, amplitude=parent_amplitude[bvec].amplitude)
                              for (ivec, bvec) in enumerate(reduced_basis_list))

  for (bvec_prime, (bvec, amplitude)) in parent_amplitude
    bvec_prime == bvec && continue
    reduced_basis_lookup[bvec_prime] = (index=reduced_basis_lookup[bvec].index, amplitude=amplitude)
  end
  return ReducedHilbertSpaceRealization{QN, BR, ComplexType}(hsr, trans_group, reduced_basis_list, reduced_basis_lookup)
end


function symmetry_reduce_parallel(hsr ::HilbertSpaceRealization{QN, BR},
                trans_group ::TranslationGroup,
                fractional_momentum ::AbstractVector{Rational};
                ComplexType::DataType=ComplexF64) where {QN, BR}
  ik = findfirst(collect(
    trans_group.fractional_momenta[ik] == fractional_momentum
    for ik in 1:length(trans_group.fractional_momenta) ))
  
  isnothing(ik) && throw(ArgumentError("fractional momentum $(fractional_momentum) not an irrep of the translation group"))

  # check if fractional momentum is compatible with translation group ?
  #k = float.(fractional_momentum) .* 2π
  phases = trans_group.character_table[ik, :]
  #[ cis(dot(k, t)) for t in trans_group.translations]

  n_basis = length(hsr.basis_list)

  mutex = Threads.Mutex()

  nthreads = Threads.nthreads()
  local_reduced_basis_list = [UInt[] for i in 1:nthreads]
  ParentAmplitudeDictType = Dict{BR, NamedTuple{(:parent,:amplitude), Tuple{BR,ComplexType}}}
  local_parent_amplitude = [ParentAmplitudeDictType() for i in 1:nthreads]

  visited = BitArray(undef, n_basis)
  visited .= false

  Threads.@threads for ivec in 1:n_basis
    visited[ivec] && continue
    id = Threads.threadid()
    bvec = hsr.basis_list[ivec]

    ψ = SparseState{ComplexType, UInt}(hsr.hilbert_space)
    identity_translations = Vector{Int}[]
    for i in 1:length(trans_group.elements)
      t = trans_group.translations[i]
      g = trans_group.elements[i]
      p = phases[i]

      bvec_prime = apply_symmetry(hsr.hilbert_space, g, bvec)
      ψ[bvec_prime] += p
      if bvec_prime == bvec
        push!(identity_translations, t)
      end
    end
    (!is_compatible(fractional_momentum, identity_translations)) && continue
    (bvec != minimum(keys(ψ.components))) && continue

    clean!(ψ)
    @assert !isempty(ψ)

    ivec_primes = [hsr.basis_lookup[bvec_prime] for bvec_prime in keys(ψ.components)]

    lock(mutex)
    if any(visited[ivec_primes]) 
      unlock(mutex)
      continue
    else
      visited[ivec_primes] .= true
      unlock(mutex)
    end

    normalize!(ψ)
    push!(local_reduced_basis_list[id], bvec)
    for (bvec_prime, amplitude) in ψ.components
      local_parent_amplitude[id][bvec_prime] = (parent=bvec, amplitude=amplitude)
    end
  end

  reduced_basis_list ::Vector{BR} = vcat(local_reduced_basis_list...)
  sort!(reduced_basis_list)

  parent_amplitude = ParentAmplitudeDictType()
  merge!(parent_amplitude, local_parent_amplitude...)
  reduced_basis_lookup = Dict(bvec => (index=ivec, amplitude=parent_amplitude[bvec].amplitude)
                              for (ivec, bvec) in enumerate(reduced_basis_list))

  for (bvec_prime, (bvec, amplitude)) in parent_amplitude
    bvec_prime == bvec && continue
    reduced_basis_lookup[bvec_prime] = (index=reduced_basis_lookup[bvec].index, amplitude=amplitude)
  end
  return ReducedHilbertSpaceRealization{QN, BR, ComplexType}(hsr, trans_group, reduced_basis_list, reduced_basis_lookup)
end


function materialize(rhsr :: ReducedHilbertSpaceRealization{QN, BR, C},
                     operator ::AbstractOperator;
                     tol::Real=sqrt(eps(Float64))) where {QN, BR, C}
  # TODO CHECK IF THe OPERATOR HAS TRANSLATION SYMMETRY
  rows = Int[]
  cols = Int[]
  vals = ComplexF64[]
  
  hs = rhsr.parent_hilbert_space_realization.hilbert_space
  err = 0.0

  for (irow, brow) in enumerate(rhsr.basis_list)
    ampl_row = rhsr.basis_lookup[brow].amplitude
    ψrow = SparseState{C, UInt}(hs, brow=>1/ampl_row)
    ψcol = SparseState{C, UInt}(hs)
    apply!(ψcol, ψrow, operator)
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

  if !isempty(vals) && maximum(imag.(vals)) < tol
    vals = real.(vals)
  end

  n = length(rhsr.basis_list)
  return (sparse(rows, cols, vals, n, n), err)
end



function materialize_parallel(rhsr :: ReducedHilbertSpaceRealization{QN, BR, C},
                     operator ::AbstractOperator;
                     tol::Real=sqrt(eps(Float64))) where {QN, BR, C}
  # TODO CHECK IF THe OPERATOR HAS TRANSLATION SYMMETRY
  hs = rhsr.parent_hilbert_space_realization.hilbert_space

  nthreads = Threads.nthreads()
  local_rows = [ Int[] for i in 1:nthreads]
  local_cols = [ Int[] for i in 1:nthreads]
  local_vals = [ C[] for i in 1:nthreads]
  local_err =  Float64[0.0 for i in 1:nthreads]

  n_basis = length(rhsr.basis_list)
  
  Threads.@threads for irow in 1:n_basis
    id = Threads.threadid()
    brow = rhsr.basis_list[irow]

    ampl_row = rhsr.basis_lookup[brow].amplitude
    ψrow = SparseState{C, UInt}(hs, brow=>1/ampl_row)
    ψcol = SparseState{C, UInt}(hs)
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

  rows ::Vector{Int} = vcat(local_rows...) 
  cols ::Vector{Int} = vcat(local_cols...) 
  vals ::Vector{C} = vcat(local_vals...) 
  err ::Float64 = sum(local_err) 

  n = length(rhsr.basis_list)
  if isempty(vals)
    vals = Float64[]
  elseif isapprox( maximum(abs.(imag.(vals))), 0)
    vals = real.(vals)
  end
  return (sparse(rows, cols, vals, n, n), err)
end

