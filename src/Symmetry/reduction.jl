export ReducedHilbertSpaceRealization
export symmetry_reduce
export materialize

struct ReducedHilbertSpaceRealization
  parent_hilbert_space_realization ::HilbertSpaceRealization
  translation_group ::TranslationGroup
  basis_list ::Vector{UInt}
  basis_lookup ::Dict{UInt, NamedTuple{(:index, :amplitude), Tuple{Int, ComplexF64}}}
end

function symmetry_reduce(hsr ::HilbertSpaceRealization,
                trans_group ::TranslationGroup,
                fractional_momentum ::AbstractVector{Rational})
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

    ψ = SparseState{ComplexF64, UInt}(hsr.hilbert_space)
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
  return ReducedHilbertSpaceRealization(hsr, trans_group, reduced_basis_list, reduced_basis_lookup)
end


function materialize(rhsr :: ReducedHilbertSpaceRealization,
                     operator ::AbstractOperator;
                     tol::Real=sqrt(eps(Float64)))
  # TODO CHECK IF THe OPERATOR HAS TRANSLATION SYMMETRY
  rows = Int[]
  cols = Int[]
  vals = ComplexF64[]
  
  hs = rhsr.parent_hilbert_space_realization.hilbert_space
  err = 0.0

  for (irow, brow) in enumerate(rhsr.basis_list)
    ampl_row = rhsr.basis_lookup[brow].amplitude
    ψrow = SparseState{ComplexF64, UInt}(hs, brow=>1/ampl_row)
    ψcol = SparseState{ComplexF64, UInt}(hs)
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

  if maximum(imag.(vals)) < tol
    vals = real.(vals)
  end

  n = length(rhsr.basis_list)
  return (sparse(rows, cols, vals, n, n), err)
end
