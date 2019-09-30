
struct TranslationGroup <: AbstractSymmetryGroup
  generators ::Vector{Permutation}
  
  translations ::Vector{Vector{Int}}
  elements ::Vector{Permutation}
  fractional_momenta ::Vector{Vector{Rational}}
  #representations ::Vector{ Vector{ComplexF64} }
  
  function TranslationGroup(generators::AbstractArray{Permutation})
    for g1 in generators, g2 in generators
      @assert g1 * g2 == g2 * g1
    end
    shape = [g.cycle_length for g in generators]
    translations = vcat( collect( Iterators.product([0:g.cycle_length-1 for g in generators]...) )...)
    translations = [ [x...] for x in translations]
    elements = [prod(gen^d for (gen, d) in zip(generators, dist)) for (ig, dist) in enumerate(translations)]
    momentum(sub) = [x//d for (x, d) in zip(sub, shape)]
    fractional_momenta = [momentum(sub) for sub in translations]
    #momentum(sub) = [x//d for (x, d) in zip(sub, shape)]
    #momenta = [momentum(sub) for sub in translations]
    #representations = collect(Iterators.product([0:g.cycle_length-1 for g in generators]...))
    return new(generators, translations, elements, fractional_momenta)
  end
end

function is_compatible(
    fractional_momentum ::AbstractVector{Rational},
    identity_translation ::AbstractVector{<:Integer}
    )
  value = sum( i * j for (i,j) in zip(fractional_momentum, identity_translation))
  return mod(value, 1) == 0
end


function is_compatible(
    fractional_momentum ::AbstractVector{Rational},
    identity_translations ::AbstractVector{<:AbstractVector{<:Integer}}
  )
  return all(is_compatible(fractional_momentum, t) for t in identity_translations)
end

