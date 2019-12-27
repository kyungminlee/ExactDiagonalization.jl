# Examples

## S=1/2 Heisenberg Chain

```julia
using SparseArrays
using LinearAlgebra
using Arpack
using Plots
using ExactDiagonalization
using TightBindingLattice

println("# S=1/2 Heisenberg Chain")
n_sites = 8;
(hs, σ) = ExactDiagonalization.Toolkit.spin_half_system(n_sites)
println("Quantum number sectors (2Sz) : ", quantum_number_sectors(hs))

S = Dict(μ => sum(σ(i, μ) for i in 1:n_sites) for μ in [:x, :y, :z])
spin_squared = simplify(S[:x]^2 + S[:y]^2 + S[:z]^2)
j1 = simplify(sum(σ(i, j) * σ(mod(i, n_sites)+1 , j) for i in 1:n_sites for j in [:x, :y, :z]))
hs_rep = represent(hs);     # Use FrozenSortedArrayIndex{UInt} for basis lookup

println("## All sectors at once")
plt = plot(size=(300, 400))
begin
  j1_rep = represent(hs_rep, j1)
  eigenvalues, eigenvectors = eigs(j1_rep; nev=32, which=:SR, )
  eigenvalues = real.(eigenvalues)
  println("E : ", eigenvalues[1:5])
  scatter!(plt,
           zeros(size(eigenvalues)), eigenvalues,
           markershape=:hline,
           markersize=10,
           markerstrokecolor=:red,
           legend=:none)
end

println("## Sz=0, each momentum sectors")
hs_sector = HilbertSpaceSector(hs, 0)
hs_rep = represent_dict(hs_sector) # Use Dict{UInt, Int} for basis lookup
translation_group = TranslationGroup([Permutation([ mod(i, n_sites)+1 for i in 1:n_sites])])
ks = translation_group.fractional_momenta
for (ik, k) in enumerate(ks)
  hs_redrep = symmetry_reduce(hs_rep, translation_group, k)
  j1_redrep = represent(hs_redrep, j1)
  #j1_redrep_sparse = sparse(j1_redrep)   # Make a sparse matrix
  j1_redrep_dense = Matrix(j1_redrep)   # Make a dense matrix
  eigenvalues = eigvals(Hermitian(j1_redrep_dense))
  println("E(k=", join(string.(k), ","), ") : ", eigenvalues[1:5])
  scatter!(plt,
           ones(size(eigenvalues)).*ik, eigenvalues,
           markershape=:hline,
           markersize=10,
           markerstrokecolor=:blue,
           legend=:none)
end
xticks!(plt, collect(0:n_sites), ["All", string.(1:n_sites)...])
xlims!(plt, -1, n_sites+1)
xlabel!(plt, "Momentum index")
ylabel!(plt, "Energy")
```