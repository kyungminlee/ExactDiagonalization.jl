using Test
using ExactDiagonalization

using TightBindingLattice
using ExactDiagonalization.Toolkit: pauli_matrix

@testset "symmetry_apply" begin

  unitcell = make_unitcell(1.0; OrbitalType=String)
  addorbital!(unitcell, "Spin", FractCoord([0], [0.0]))
  lattice = make_lattice(unitcell, 4)

  # Test State and Site
  up = State("Up", 1)
  dn = State("Dn",-1)
  spin_site = Site([up, dn])

  hs = HilbertSpace(repeat([spin_site], 4))
  hss = HilbertSpaceSector(hs, 0)
  transop = SitePermutation([2,3,4,1])
  invop = SitePermutation([1,4,3,2])

  @test symmetry_apply(hs, transop, 0b0001) == 0b0010
  @test symmetry_apply(hss, transop, 0b0101) == 0b1010

  nop = NullOperator()
  pop1 = PureOperator{Float64, UInt}(0b1101, 0b0101, 0b1100, 2.0)
  pop2 = PureOperator{Float64, UInt}(0b0010, 0b0000, 0b0010, 3.0)
  sop = pop1 + pop2

  @test symmetry_apply(hs, transop, nop) == nop
  @test symmetry_apply(hs, transop, pop1) == PureOperator{Float64, UInt}(0b1011, 0b1010, 0b1001, 2.0)
  @test symmetry_apply(hs, transop, sop) == SumOperator{Float64, UInt}([
      symmetry_apply(hs, transop, pop1), symmetry_apply(hs, transop, pop2)])

  @test symmetry_apply(hs, invop, nop) == nop
  @test symmetry_apply(hs, invop, pop1) == PureOperator{Float64, UInt}(0b0111, 0b0101, 0b0110, 2.0)
  @test symmetry_apply(hs, invop, sop) == SumOperator{Float64, UInt}([
      symmetry_apply(hs, invop, pop1), symmetry_apply(hs, invop, pop2)])

  σ = Dict( (isite, j) => pauli_matrix(hs, isite, j) for isite in 1:4, j in [:x, :y, :z, :+, :-])
  j1 = simplify(sum(σ[i, j] * σ[mod(i, 4) + 1 , j] for i in 1:4, j in [:x, :y, :z]))

  @test is_invariant(hs, transop, nop)
  @test !is_invariant(hs, transop, pop1)
  @test !is_invariant(HilbertSpaceSector(hs, 0), transop, pop1)
  @test !is_invariant(hs, transop, sop)
  @test !is_invariant(HilbertSpaceSector(hs, 0), transop, sop)
  @test is_invariant(hs, transop, j1)
  @test is_invariant(HilbertSpaceSector(hs, 0), transop, j1)

  @test is_invariant(hs, invop, nop)
  @test !is_invariant(hs, invop, pop1)
  @test !is_invariant(HilbertSpaceSector(hs, 0), invop, pop1)
  @test !is_invariant(hs, invop, sop)
  @test !is_invariant(HilbertSpaceSector(hs, 0), invop, sop)
  @test is_invariant(hs, invop, j1)
  @test is_invariant(HilbertSpaceSector(hs, 0), invop, j1)

  tsymbed = translation_symmetry_embedding(lattice)

  @test is_invariant(hs, tsymbed, nop)
  @test !is_invariant(hs, tsymbed, pop1)
  @test !is_invariant(HilbertSpaceSector(hs, 0), tsymbed, pop1)
  @test !is_invariant(hs, tsymbed, sop)
  @test !is_invariant(HilbertSpaceSector(hs, 0), tsymbed, sop)
  @test is_invariant(hs, tsymbed, j1)
  @test is_invariant(HilbertSpaceSector(hs, 0), tsymbed, j1)

  psymbed = embed(lattice, project(PointSymmetryDatabase.get(2), [1 0 0;]))

  @test is_invariant(hs, psymbed, nop)
  @test !is_invariant(hs, psymbed, pop1)
  @test !is_invariant(HilbertSpaceSector(hs, 0), psymbed, pop1)
  @test !is_invariant(hs, psymbed, sop)
  @test !is_invariant(HilbertSpaceSector(hs, 0), psymbed, sop)
  @test is_invariant(hs, psymbed, j1)
  @test is_invariant(HilbertSpaceSector(hs, 0), psymbed, j1)

end
