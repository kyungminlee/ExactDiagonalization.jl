using Test
using ExactDiagonalization

@testset "SpinHalfHeisenbergChain" begin
  # Set up
  QN = Int

  # Test State and Site
  up = State{QN}("Up", 1)
  dn = State{QN}("Dn",-1)
  spin_site = Site{QN}([up, dn])

  @test up.quantum_number == 1
  @test dn.quantum_number == -1
  @test bitwidth(spin_site) == 1
  @test get_state(spin_site, 0b0) == up
  @test get_state(spin_site, 0b1) == dn

  let
    ze = State{QN}("Ze", 0)
    spin_one_site = Site{QN}([up, ze, dn])
    @test bitwidth(spin_one_site) == 2
  end

  # Test AbstractHilbertSpace
  n = 4
  hs = AbstractHilbertSpace(repeat([spin_site], n))
  hs2 = AbstractHilbertSpace{QN}()
  for i in 1:n
    add_site!(hs2, spin_site)
  end

  @test length(hs.sites) == length(hs2.sites)
  @test hs.bitwidths == hs2.bitwidths
  @test hs.bitoffsets == hs2.bitoffsets

  # Test ConcreteHilbertSpace
  chs = concretize(hs; BR=UInt64)
  @test dimension(chs) == 2^n
  @test chs.basis_list == collect(UInt64(0):UInt64(2^n-1))
  @test all(chs.basis_lookup[basis] == ibasis for (ibasis, basis) in enumerate(chs.basis_list))

  sector_bases = Dict{Int, Vector{UInt}}(
    +4 => [0b0000],
    +2 => [0b0001, 0b0010, 0b0100, 0b1000], 
     0 => [0b0011, 0b0101, 0b0110, 0b1001, 0b1010, 0b1100],
    -2 => [0b0111, 0b1011, 0b1101, 0b1110],
    -4 => [0b1111],
  )
  sectors = quantum_number_sectors(hs)
  @test sectors == sort(collect(keys(sector_bases)))
  for (q, target_basis) in sector_bases
    chs1 = concretize(hs, Set([q]); BR=UInt64)
    chs2 = concretize(hs, target_basis)
    chs3 = concretize(hs, UInt64[0x1])
    @test chs1 == chs1
    @test chs1 == chs2
    @test chs1 != chs3
    @test !(chs1 == chs3)
  end



end
