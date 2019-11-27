using SparseArrays
using LinearAlgebra
using ExactDiagonalization
using TightBindingLattice
using MinimalPerfectHash

QN = Int;
up = State{QN}("Up", 1);
dn = State{QN}("Dn",-1);
spin_site = Site{QN}([up, dn]);

n_sites = 16;
hs = HilbertSpace([spin_site for i in 1:n_sites]);
hss = HilbertSpaceSector(hs, 0)

hsr = represent_dict(hss);

function pauli_matrix(hs::HilbertSpace, isite ::Integer, j ::Symbol)
  if j == :x
    return pure_operator(hs, isite, 1, 2, 1, UInt) + pure_operator(hs, isite, 2, 1, 1, UInt)
  elseif j == :y
    return pure_operator(hs, isite, 1, 2, -im, UInt) + pure_operator(hs, isite, 2, 1, im, UInt)
  elseif j == :z
    return pure_operator(hs, isite, 1, 1, 1, UInt) + pure_operator(hs, isite, 2, 2, -1, UInt)
  elseif j == :+
    return pure_operator(hs, isite, 1, 2, 1, UInt)
  elseif j == :-
    return pure_operator(hs, isite, 2, 1, 1, UInt)
  else
    throw(ArgumentError("pauli matrix of type $(j) not supported"))
  end
end;

σ = Dict( (isite, j) => pauli_matrix(hs, isite, j) for isite in 1:n_sites, j in [:x, :y, :z, :+, :-]);
Sx = sum(σ[i,:x] for i in 1:n_sites)
Sy = sum(σ[i,:y] for i in 1:n_sites)
Sz = sum(σ[i,:z] for i in 1:n_sites)

spin_squared = simplify( Sx^2 + Sy^2 + Sz^2 )
j1 = sum(σ[(i, j)] * σ[( mod(i, n_sites) + 1 , j)] for i in 1:n_sites, j in [:x, :y, :z]);

translation_group = TranslationGroup([Permutation([ mod(i, n_sites)+1 for i in 1:n_sites])])
ks = translation_group.fractional_momenta
rhsr = symmetry_reduce(hsr, translation_group, ks[1])

j1_redrep = represent(rhsr, j1)
j1_redrep_sparse = sparse(j1_redrep)

using BenchmarkTools
using Arpack




function simplify_combine(so ::SumOperator{S, BR}; tol::Real=sqrt(eps(Float64))) where {S, BR}
  bw = sizeof(BR)*8
  diagonal_terms = [t for t in so.terms if t.bitrow == t.bitcol]
  offdiagonal_terms = [t for t in so.terms if t.bitrow != t.bitcol]

  @label simplify_combine_loop_start

  for ib in 0:(bw-1)
    mask = make_bitmask(ib+1, ib)

    for (i1, t1) in enumerate(diagonal_terms)
      (t1.bitmask & mask) == 0x0 && continue
      t1.bitrow != t1.bitcol && continue

      for i2 in (i1+1):length(diagonal_terms)
        t2 = diagonal_terms[i2]
        (t2.bitmask & mask) == 0x0 && continue
        t2.bitrow != t2.bitcol && continue

        if ( ((t1.bitmask & ~mask) == (t2.bitmask & ~mask)) &&
             ((t1.bitrow  & ~mask) == (t2.bitrow  & ~mask)) &&
             ((t1.bitrow  &  mask) != (t2.bitrow  &  mask)) )

          if isapprox(t1.amplitude, t2.amplitude; atol=tol)
            i_small, i_large = i1, i2
            new_bitmask = t1.bitmask & (~mask)
            new_bitrow  = t1.bitrow  & (~mask)
            new_bitcol  = t1.bitcol  & (~mask)
            new_amplitude = t1.amplitude

            diagonal_terms[i_small] = PureOperator{S, BR}(new_bitmask, new_bitrow, new_bitcol, new_amplitude)
            deleteat!(diagonal_terms, i_large)

            @goto simplify_combine_loop_start
            #dirty = true
            #break
          else
            # different amplitude. do nothing
          end # if amplitude
        end # if bits

      end # for i2
    end # for i1
  end # for ib

  #end # while dirty
  return simplify(SumOperator{S, BR}([diagonal_terms..., offdiagonal_terms...]))
end


spin_squared_2 = simplify_combine(spin_squared)

# eigs(j1_redrep)
# @btime eigs(j1_redrep)
# @timev eigs(j1_redrep)
