using SparseArrays
using LinearAlgebra
using ExactDiagonalization

QN = Int;
up = State{QN}("Up", 1);
dn = State{QN}("Dn",-1);
spin_site = Site{QN}([up, dn]);

n_sites = 5;
hs = HilbertSpace([spin_site for i in 1:n_sites]);
hss = HilbertSpaceSector(hs, 1)
hsr = represent(hss);

function pauli_matrix(hs::HilbertSpace, isite ::Integer, j ::Symbol)
  if j == :x
    return pure_operator(hs, isite, 1, 2, 1; dtype=UInt) + pure_operator(hs, isite, 2, 1, 1; dtype=UInt)
  elseif j == :y
    return pure_operator(hs, isite, 1, 2, -im; dtype=UInt) + pure_operator(hs, isite, 2, 1, im; dtype=UInt)
  elseif j == :z
    return pure_operator(hs, isite, 1, 1, 1; dtype=UInt) + pure_operator(hs, isite, 2, 2, -1; dtype=UInt)
  elseif j == :+
    return pure_operator(hs, isite, 1, 2, 1; dtype=UInt)
  elseif j == :-
    return pure_operator(hs, isite, 2, 1, 1; dtype=UInt)
  else
    throw(ArgumentError("pauli matrix of type $(j) not supported"))
  end
end;

σ = Dict( (isite, j) => pauli_matrix(hs, isite, j) for isite in 1:n_sites, j in [:x, :y, :z, :+, :-]);

j1 = sum(σ[(i, j)] * σ[( mod(i, n_sites) + 1 , j)] for i in 1:n_sites, j in [:x, :y, :z]);

j1_r = represent(hsr, j1)

n = dimension(hsr)
for icol in 1:n
  bcol = hsr.basis_list[icol]
  @show bcol
  for (brow, v) in collect(get_slice(j1, bcol))
    println(get(hsr.basis_lookup, brow, -1) => v)
  end
  @show icol
  for (irow, v) in collect(get_slice(j1_r, icol))
    println(irow => v)
  end
  println()
  # ket = zeros(ComplexF64, n)
  # bra = zeros(ComplexF64, n)
  # ket[irow] = 1.0
  # apply_unsafe!(bra, j1_r, ket)
  # println(real.(bra))
end


translation_group = TranslationGroup([Permutation([ mod(i, n_sites)+1 for i in 1:n_sites])])
ks = translation_group.fractional_momenta

j1_rS = Hermitian(Matrix(sparse(j1_r)))
@show eigvals(j1_rS)

for k in ks
  @show k
  rhsr = symmetry_reduce(hsr, translation_group, k)
  @show dimension(rhsr)
  j1_s = represent(rhsr, j1)
  j1_sS = Hermitian(Matrix(sparse(j1_s)))
  @show eigvals(j1_sS)
end
