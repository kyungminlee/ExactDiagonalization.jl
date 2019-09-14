using Test
using ExactDiagonalization

@testset "SpinHalfHeisenbergChain" begin
  QN = Int
  up = State{QN}("Up", 1)
  dn = State{QN}("Dn",-1)
  spin_site = Site{QN}([up, dn])
  n = 4
  hs1 = AbstractHilbertSpace{QN}(repeat([spin_site], n))
  hs2 = AbstractHilbertSpace{QN}()
  for i in 1:n
      add_site!(hs2, spin_site)
  end
  @test length(hs1.sites) == length(hs2.sites)
  @test hs1.bitwidths == hs2.bitwidths
  @test hs1.bitoffsets == hs2.bitoffsets
end
