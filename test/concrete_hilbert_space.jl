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
    @test dimension(chs_all) == length(chs_all.basis_list)
    @test chs_all == concretize(hs, sectors)

    # empty QN
    @test concretize(hs, QN[]).basis_list == []

    qn_basis = Dict{QN, Vector{UInt}}(
                  QN([ 2,-2]) => [0b1100],
                  QN([ 2, 0]) => [0b0100, 0b1000],
                  QN([ 2, 2]) => [0b0000],
                  QN([ 3,-3]) => [0b1110],
                  QN([ 3,-1]) => [0b0110, 0b1010, 0b1101], # UDD DUD DDU
                  QN([ 3, 1]) => [0b0010, 0b0101, 0b1001], # UUD UDU DUU
                  QN([ 3, 3]) => [0b0001], # UUU
              )

    # Test each sector explicitly
    for (qn, basis_list) in qn_basis
      chs = concretize(hs, qn)
      @test chs.basis_list == basis_list
      @test all(ibasis == chs.basis_lookup[basis] for (ibasis, basis) in enumerate(chs.basis_list))
      @test dimension(chs) == length(chs.basis_list)
      chs2 = concretize(hs, basis_list)
      @test chs == chs2
    end

    for qn in sectors
      chs = concretize(hs, qn)
      @test chs_all != chs
      @test all(ibasis == chs.basis_lookup[basis] for (ibasis, basis) in enumerate(chs.basis_list))
      @test dimension(chs) == length(chs.basis_list)
      @test chs == concretize(hs, chs.basis_list)
    end

    @test concretize(hs, QN[[ 3, 3], [ 3, 1]]).basis_list == [0b0001, 0b0010, 0b0101, 0b1001]
    @test concretize(hs, Set{QN}([[ 3, 3], [ 3, 1]])).basis_list == [0b0001, 0b0010, 0b0101, 0b1001]
  end
end