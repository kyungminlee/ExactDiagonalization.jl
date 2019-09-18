using Test
using ExactDiagonalization

using StaticArrays

@testset "State" begin
  @test isimmutable(State)
  let
    s = State("MyState")
    @test typeof(s) == State{Int}
    @test s.name == "MyState"
    @test s.quantum_number == 0
  end

  let
    s = State("Up", 1)
    @test typeof(s) == State{Int}
    @test s.name == "Up"
    @test s.quantum_number == 1
  end

  let
    s = State{SVector{2, Int}}("Up", SVector{2, Int}([1, 1]))
    @test s.name == "Up"
    @test s.quantum_number == [1, 1]
  end
end

@testset "Site" begin
  @test isimmutable(Site)
  @testset "spin-half" begin
    up = State("Up", 1)
    dn = State("Dn",-1)
    site1 = Site([up, dn])
    site2 = Site{Int}([up, dn])
    @test site1 == site2
    @test dimension(site1) == 2
    @test get_state(site1, 0x0000000) == up
    @test get_state(site1, 0x0000001) == dn
  end
  
  @testset "spin-charge" begin
    QN = SVector{2, Int}
    em = State("Em", QN( 0, 0))
    up = State("Up", QN( 1, 1))
    dn = State("Dn", QN(-1, 1))
    ud = State("UpDn", QN( 0, 2))
    @test_throws MethodError State{QN}("X", 1)

    site = Site([em, up, dn])
    @test_throws MethodError Site([em, up, dn, State("X", 1)])
    @test bitwidth(site) == 2
    @test dimension(site) == 3
    @test get_state(site, 0b0000000) == em
    @test get_state(site, 0b0000001) == up
    @test get_state(site, 0b0000010) == dn
    @test_throws BoundsError get_state(site, 0b0000011) #TODO best exception type?
  end
end

@testset "HilbertSpace" begin
  @test isimmutable(AbstractHilbertSpace)
  @testset "spinhalf" begin
    QN = Int
    up = State("Up", 1)
    dn = State{QN}("Dn",-1)
    spin_site = Site([up, dn])
    @test AbstractHilbertSpace{QN}().sites == []
    @test AbstractHilbertSpace{QN}().bitwidths == []
    @test AbstractHilbertSpace{QN}().bitoffsets == [0]
    hs = AbstractHilbertSpace([spin_site, spin_site, spin_site, spin_site])
    hs2 = AbstractHilbertSpace{QN}([spin_site, spin_site, spin_site, spin_site])
    @test hs == hs2

    @test get_bitmask(hs, 1) == 0b0001
    @test get_bitmask(hs, 2) == 0b0010
    @test get_bitmask(hs, 3) == 0b0100
    @test get_bitmask(hs, 4) == 0b1000
    @test_throws BoundsError get_bitmask(hs, 5)
    @test_throws BoundsError get_bitmask(hs,-1)

    @test quantum_number_sectors(hs) == [-4, -2, 0, 2, 4]
    @test get_quantum_number(hs, 0b0000) == 4
    @test get_quantum_number(hs, 0b0001) == 2
    @test get_quantum_number(hs, 0b0010) == 2
    @test get_quantum_number(hs, 0b0011) == 0
    @test get_quantum_number(hs, 0b0100) == 2
    @test get_quantum_number(hs, 0b0101) == 0
    @test get_quantum_number(hs, 0b0110) == 0
    @test get_quantum_number(hs, 0b0111) == -2

    @test get_quantum_number(hs, [1,1,1,1]) == 4
    @test get_quantum_number(hs, [2,1,1,1]) == 2
    @test get_quantum_number(hs, [1,2,1,1]) == 2
    @test get_quantum_number(hs, [2,2,1,1]) == 0
    @test get_quantum_number(hs, [1,1,2,1]) == 2
    @test get_quantum_number(hs, [2,1,2,1]) == 0
    @test get_quantum_number(hs, [1,2,2,1]) == 0
    @test get_quantum_number(hs, [2,2,2,1]) == -2
  end

  @testset "charge-spin" begin
    QN = SVector{2, Int}
    em = State("Em", QN( 0, 0))  # charge and spin
    up = State("Up", QN( 1, 1))
    dn = State("Dn", QN( 1,-1))
    spin_site = Site([up, dn])
    site = Site([em, up, dn])
    hs = AbstractHilbertSpace([site, site, spin_site, site])
    @test hs.bitoffsets[end] == 2 + 2 + 1 + 2
    @test get_bitmask(hs, 1) == 0b0000011
    @test get_bitmask(hs, 2) == 0b0001100
    @test get_bitmask(hs, 3) == 0b0010000
    @test get_bitmask(hs, 4) == 0b1100000
    @test_throws BoundsError get_bitmask(hs, 5)
    @test_throws BoundsError get_bitmask(hs,-1)

    sectors = [[1, -1], [1, 1],
               [2, -2], [2, 0], [2, 2],
               [3, -3], [3, -1], [3, 1], [3, 3],
               [4, -4], [4, -2], [4, 0], [4, 2], [4, 4]]
    @test quantum_number_sectors(hs) == sectors

    @test get_state_index(hs, 0b0000000, 1) == 1
    @test get_state_index(hs, 0b0000001, 1) == 2
    @test get_state_index(hs, 0b0000010, 1) == 3
    @test get_state_index(hs, 0b0000011, 1) == 4

    @test get_state_index(hs, 0b0000011, 2) == 1
    @test get_state_index(hs, 0b0000111, 2) == 2
    @test get_state_index(hs, 0b0001011, 2) == 3
    @test get_state_index(hs, 0b0001111, 2) == 4

    @test get_state_index(hs, 0b0000000, 3) == 1
    @test get_state_index(hs, 0b0010000, 3) == 2
    @test get_state_index(hs, 0b1101111, 3) == 1
    @test get_state_index(hs, 0b1111111, 3) == 2

    @test get_quantum_number(hs, 0b0000000) == [1, 1]
    @test get_quantum_number(hs, 0b0000001) == [2, 2]
    @test get_quantum_number(hs, 0b0000010) == [2, 0]
    @test_throws BoundsError get_quantum_number(hs, 0b0000011)
    @test get_quantum_number(hs, 0b0000100) == [2, 2] # em - up - up - em
    @test get_quantum_number(hs, 0b0000101) == [3, 3] # up - up - up - em
    @test get_quantum_number(hs, 0b0000110) == [3, 1] # dn - up - up - em
    @test_throws BoundsError get_quantum_number(hs, 0b0000111)

    @test update(hs, 0b1111110, 1, 1) == 0b1111100
    @test update(hs, 0b1111110, 1, 2) == 0b1111101
    @test update(hs, 0b1111110, 1, 3) == 0b1111110
    @test_throws BoundsError update(hs, 0b1111110, 1, 4)
  end


end