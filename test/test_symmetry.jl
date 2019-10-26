using Test
using ExactDiagonalization

using LinearAlgebra

include("pauli_matrix.jl")

@testset "Permutation" begin
  @test_throws ArgumentError Permutation([1,2,4])
  @test_throws OverflowError Permutation([mod(x, 4096)+1 for x in 1:4096])
  p0 = Permutation([1,2,3,4])
  p1 = Permutation([2,3,4,1])
  p2 = Permutation([3,4,1,2])
  p3 = Permutation([4,1,2,3])

  @test p1 * p2 == p3
  @test p1 != p3
  @test p1^0 == p0
  @test p1^1 == p1
  @test p1^2 == p2
  @test p1^3 == p3

  @test p0.cycle_length == 1
  @test p1.cycle_length == 4
  @test p2.cycle_length == 2
  @test p3.cycle_length == 4

  @test_throws ArgumentError Permutation([1,2,3,4]) * Permutation([1,2,3,4,5])
  @test hash(Permutation(Int[1,2,3,4])) == hash(Int[1,2,3,4])
end

@testset "translation" begin
  @testset "constructor exceptions" begin
    @test_throws ArgumentError TranslationGroup([Permutation([2,3,4,1]), Permutation([3,4,1,2])])
    @test_throws ArgumentError TranslationGroup([Permutation([2,1,3,4]), Permutation([1,3,2,4])])
  end


  t1 = Permutation([2,3,1, 5,6,4])
  t2 = Permutation([4,5,6, 1,2,3])
  g = TranslationGroup([t1, t2])

  @test g.generators == [t1, t2]
  @test g.translations == [[0,0], [1,0], [2,0], [0,1], [1,1], [2,1]]
  @test Set(g.elements) == Set([t1^d1*t2^d2 for d1 in 0:2 for d2 in 0:1])
  @test length(Set(g.elements)) == 2*3
  @test g.fractional_momenta == [[0//3, 0//2], [1//3, 0//2], [2//3, 0//2],
                                 [0//3, 1//2], [1//3, 1//2], [2//3, 1//2]]
  χ = [cis(2π * (k ⋅ t)) for k in g.fractional_momenta, t in g.translations]
  @test isapprox(g.character_table, χ; atol=sqrt(eps(Float64)))

  @test is_compatible([0//1, 0//1], [0,0])
  @test !is_compatible([0//1, 1//2], [0,1])

  @test is_compatible([0//1, 0//1], [[0,0]])
  @test !is_compatible([0//1, 1//2], [[0,0], [0,1]])
end

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
