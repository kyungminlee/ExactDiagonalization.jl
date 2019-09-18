using Test
using ExactDiagonalization
using StaticArrays

@testset "ConcreteHilbertSpace" begin
  QN = SVector{2, Int}
  em = State("Em", QN( 0, 0))  # charge and spin
  up = State("Up", QN( 1, 1))
  dn = State("Dn", QN( 1,-1))
  spin_site = Site([up, dn])
  site = Site([em, up, dn])
  hs = AbstractHilbertSpace([site, spin_site, spin_site]) # f s s
  sectors = quantum_number_sectors(hs)
  @test sectors == QN[[2, -2], [2, 0], [2, 2], [3, -3], [3, -1], [3, 1], [3, 3]]
  @testset "concretize" begin
    chs_all = concretize(hs)
    @test chs_all.basis_list == UInt[
      0b0000, 0b0001, 0b0010,
      0b0100, 0b0101, 0b0110,
      0b1000, 0b1001, 0b1010,
      0b1100, 0b1101, 0b1110,
    ]
    @test all(ibasis == chs_all.basis_lookup[basis] for (ibasis, basis) in enumerate(chs_all.basis_list))
    @test chs_all == concretize(hs, sectors)

    @show sectors
    @test concretize(hs, QN([ 2,-2])).basis_list == [0b1100]
    @test concretize(hs, QN([ 2, 0])).basis_list == [0b0100, 0b1000]
    @test concretize(hs, QN([ 2, 2])).basis_list == [0b0000]
    @test concretize(hs, QN([ 3,-3])).basis_list == [0b1110]
    @test concretize(hs, QN([ 3,-1])).basis_list == [0b0110, 0b1010, 0b1101] # UDD DUD DDU
    @test concretize(hs, QN([ 3, 1])).basis_list == [0b0010, 0b0101, 0b1001] # UUD UDU DUU
    @test concretize(hs, QN([ 3, 3])).basis_list == [0b0001] # UUU

    for qn in sectors
      chs = concretize(hs, qn)
      @test all(ibasis == chs.basis_lookup[basis] for (ibasis, basis) in enumerate(chs.basis_list))
    end

    @test concretize(hs, QN[[ 3, 3], [ 3, 1]]).basis_list == [0b0001, 0b0010, 0b0101, 0b1001]
    @test concretize(hs, Set{QN}([[ 3, 3], [ 3, 1]])).basis_list == [0b0001, 0b0010, 0b0101, 0b1001]
  end
end