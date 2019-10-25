using Test
using ExactDiagonalization

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
end

@testset "translation" begin
  @testset "constructor exceptions" begin
    t1 = Permutation([2,3,4,1])
    t2 = Permutation([3,4,1,2])
    @test_throws ArgumentError TranslationGroup([t1, t2])

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
  #@show g.character_table 
end
