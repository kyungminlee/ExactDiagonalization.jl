using Test
using ExactDiagonalization

using TightBindingLattice
using ExactDiagonalization.Toolkit: pauli_matrix

@testset "symmetry_apply" begin
  QN = Int

  # Test State and Site
  up = State("Up", 1)
  dn = State("Dn",-1)
  spin_site = Site([up, dn])

  hs = HilbertSpace(repeat([spin_site], 4))
  hss = HilbertSpaceSector(hs, 0)
  p = Permutation([2,3,4,1])
  @test symmetry_apply(hs, p, 0b0001) == 0b0010
  @test symmetry_apply(hss, p, 0b0101) == 0b1010

  nop = NullOperator()
  pop1 = PureOperator{Float64, UInt}(0b1101, 0b0101, 0b1100, 2.0)
  pop2 = PureOperator{Float64, UInt}(0b0010, 0b0000, 0b0010, 3.0)
  sop = pop1 + pop2

  @test symmetry_apply(hs, p, nop) == nop
  @test symmetry_apply(hs, p, pop1) == PureOperator{Float64, UInt}(0b1011, 0b1010, 0b1001, 2.0)
  @test symmetry_apply(hs, p, sop) == SumOperator{Float64, UInt}([
      symmetry_apply(hs, p, pop1), symmetry_apply(hs, p, pop2)])

  σ = Dict( (isite, j) => pauli_matrix(hs, isite, j) for isite in 1:4, j in [:x, :y, :z, :+, :-])
  j1 = simplify(sum(σ[i, j] * σ[mod(i, 4) + 1 , j] for i in 1:4, j in [:x, :y, :z]))

  @test is_invariant(hs, p, nop)
  @test !is_invariant(hs, p, pop1)
  @test !is_invariant(HilbertSpaceSector(hs, 0), p, pop1)
  @test !is_invariant(hs, p, sop)
  @test !is_invariant(HilbertSpaceSector(hs, 0), p, sop)
  @test is_invariant(hs, p, j1)
  @test is_invariant(HilbertSpaceSector(hs, 0), p, j1)

  tg = TranslationGroup(p)
  @test is_invariant(hs, tg, nop)
  @test !is_invariant(hs, tg, pop1)
  @test !is_invariant(HilbertSpaceSector(hs, 0), tg, pop1)
  @test !is_invariant(hs, tg, sop)
  @test !is_invariant(HilbertSpaceSector(hs, 0), tg, sop)
  @test is_invariant(hs, tg, j1)
  @test is_invariant(HilbertSpaceSector(hs, 0), tg, j1)
end
