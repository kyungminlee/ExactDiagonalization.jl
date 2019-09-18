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

    @test extract(hs, 0b0000000) == [1,1,1,1]
    @test extract(hs, 0b0000001) == [2,1,1,1]
    @test extract(hs, 0b0000010) == [3,1,1,1]
    @test_throws BoundsError extract(hs, 0b0000011)

    @test extract(hs, 0b0010000) == [1,1,2,1]
    @test extract(hs, 0b0010001) == [2,1,2,1]
    @test extract(hs, 0b0010010) == [3,1,2,1]
    @test_throws BoundsError extract(hs, 0b0010011)

    @test compress(hs, [1,1,2,1]) == 0b0010000
    @test compress(hs, [2,1,2,1]) == 0b0010001
    @test compress(hs, [3,1,2,1]) == 0b0010010
    @test_throws BoundsError compress(hs, [4, 1, 2, 1])
    @test_throws ArgumentError compress(hs, [1, 2])
  end
end


@testset "SparseState" begin
  QN = SVector{2, Int}
  em = State("Em", QN( 0, 0))  # charge and spin
  up = State("Up", QN( 1, 1))
  dn = State("Dn", QN( 1,-1))
  spin_site = Site([up, dn])
  site = Site([em, up, dn])
  hs = AbstractHilbertSpace([site, site, spin_site, site])
  hs2 = AbstractHilbertSpace([site, spin_site, site, site])
  
  @testset "constructor" begin
    ψ1 = SparseState{ComplexF64, UInt}(hs, UInt(0b0010001))
    @test ψ1.hilbert_space == hs

    @test ψ1.components == Dict(UInt(0x11) => 1.0 + 0.0im)
    ψ1[UInt(0x0)] # DefaultDict
    @test ψ1[UInt(0x11)] == 1.0 + 0.0im
    @test ψ1.components != Dict(UInt(0x11) => 1.0 + 0.0im)
    
    ψ2 = SparseState{ComplexF64, UInt}(hs)
    ψ2[UInt(0b0010001)] = 1.0
    @test ψ2.components == Dict(UInt(0x11) => 1.0 + 0.0im)

    ψ3 = SparseState{ComplexF64, UInt}(hs, UInt(0b0010001) => 2.0)
    @test ψ3.components == Dict(UInt(0x11) => 2.0 + 0.0im)

    ψ4 = SparseState{ComplexF64, UInt}(hs, UInt(0b0000001) => 2.0, UInt(0b0010001) => 3.0 )
    @test ψ4.components == Dict(UInt(0x11) => 3.0 + 0.0im, UInt(0x1) => 2.0 + 0.0im)
  end

  @testset "hilbert" begin
    ψ1 = SparseState{ComplexF64, UInt}(hs, UInt(0b0010001))
    ψ2 = SparseState{ComplexF64, UInt}(hs2, UInt(0b0010001))
    @test ψ1 != ψ2
    @test_throws ArgumentError ψ1 + ψ2 
    @test_throws ArgumentError ψ1 - ψ2 
  end

  @testset "equality" begin
    ψ1 = SparseState{ComplexF64, UInt}(hs, UInt(0b0010001))
    ψ2 = SparseState{ComplexF64, UInt}(hs)
    ψ2[UInt(0b0010001)] = 1.0
    ψ3 = SparseState{ComplexF64, UInt}(hs, UInt(0b0010000))
    ψ4 = SparseState{ComplexF64, UInt}(hs)
    ψ4[UInt(0b0010001)] = 2.0

    @test ψ1 == ψ2
    @test ψ1 != ψ3
    @test ψ1 != ψ4
    
    ψ2[UInt(0x0)]
    @test ψ1 != ψ2
  end  

  @testset "unary" begin
    @testset "real" begin
      ψ1 = SparseState{Float64, UInt}(hs)
      ψ1[UInt(0b0000001)] = 3.0
      
      @test +ψ1 == ψ1
      @test -ψ1 != ψ1
      @test (-ψ1).components == Dict(UInt(0b0000001) => -3.0)
      @test real(ψ1) == ψ1
      @test imag(ψ1) != ψ1
      @test conj(ψ1) == ψ1
      @test real(ψ1).components == Dict(UInt(0b0000001) => 3.0)
      @test imag(ψ1).components == Dict()
      @test conj(ψ1).components == Dict(UInt(0b0000001) => 3.0)
    end

    @testset "complex" begin
      ψ1 = SparseState{ComplexF64, UInt}(hs)
      ψ1[UInt(0b0000001)] = 3.0 + im
      
      @test +ψ1 == ψ1
      @test -ψ1 != ψ1
      @test (-ψ1).components == Dict(UInt(0b0000001) => -3.0 - im)
      @test real(ψ1) != ψ1
      @test imag(ψ1) != ψ1
      @test conj(ψ1) != ψ1
      @test real(ψ1).components == Dict(UInt(0b0000001) => 3.0)
      @test imag(ψ1).components == Dict(UInt(0b0000001) => 1.0)
      @test conj(ψ1).components == Dict(UInt(0b0000001) => 3.0 - im)
    end
  end # unary

  @testset "binary" begin
    @testset "scalar" begin
      ψ1 = SparseState{Float64, UInt}(hs, UInt(0b0000001) => 2.0)
      @test (ψ1 * 2).hilbert_space == ψ1.hilbert_space
      @test (2 * ψ1).hilbert_space == ψ1.hilbert_space
      @test (ψ1 / 2).hilbert_space == ψ1.hilbert_space

      @test (ψ1 * 3).components == Dict(UInt(0x001) => 6.0)
      @test (3 * ψ1).components == Dict(UInt(0x001) => 6.0)
      @test (ψ1 * 3.0).components == Dict(UInt(0x001) => 6.0)
      @test (3.0 * ψ1).components == Dict(UInt(0x001) => 6.0)
      @test (ψ1 * (3.0+im)).components == Dict(UInt(0x001) => 6.0 + 2.0im)
      @test ((3.0+im) * ψ1).components == Dict(UInt(0x001) => 6.0 + 2.0im)

      @test (ψ1 / 4).components == Dict(UInt(0x001) => 0.5)
      @test (ψ1 / 4.0).components == Dict(UInt(0x001) => 0.5)
      @test (ψ1 / (4.0+0.0im)).components == Dict(UInt(0x001) => 0.5 + 0.0im)
    end
    @testset "sum" begin
      ψ1 = SparseState{Float64, UInt}(hs, UInt(0b0000001) => 2.0)
      ψ2 = SparseState{Float64, UInt}(hs, UInt(0b0000001) => 0.25, UInt(0b0001001) => 4.0)
      ψ3 = ψ1 + ψ2
      ψ4 = ψ2 + ψ1
      @test typeof(ψ3) == SparseState{Float64, UInt}
      @test typeof(ψ4) == SparseState{Float64, UInt}
      @test ψ3 == ψ4
      @test ψ3.components == Dict(UInt(0b0000001) => 2.25, UInt(0b0001001) => 4.0)

      ψ3 = ψ1 - ψ2
      ψ4 = ψ2 - ψ1
      @test typeof(ψ3) == SparseState{Float64, UInt}
      @test typeof(ψ4) == SparseState{Float64, UInt}

      @test ψ3 != ψ4
      @test ψ3 == -ψ4
      @test ψ3.components == Dict(UInt(0b0000001) => 1.75, UInt(0b0001001) => -4.0)
      @test ψ4.components == Dict(UInt(0b0000001) =>-1.75, UInt(0b0001001) => 4.0)
    end # sametype

    @testset "intertype" begin
      ψ1 = SparseState{Float64, UInt}(hs, UInt(0b0000001) => 2.0)
      ψ2 = SparseState{ComplexF64, UInt}(hs, UInt(0b0000001) => 4.0 + 0.25im, UInt(0b1001) => 3.0)
      ψ3 = ψ1 + ψ2
      ψ4 = ψ2 + ψ1
      @test ψ3 == ψ4
      @test ψ3.components == Dict(UInt(0b0000001) => 6.00 + 0.25im, UInt(0b1001) => 3.0)
      @test typeof(ψ3) == SparseState{ComplexF64, UInt}
      @test typeof(ψ4) == SparseState{ComplexF64, UInt}

      ψ3 = ψ1 - ψ2
      ψ4 = ψ2 - ψ1
      @test ψ3 == -ψ4
      @test ψ3.components == Dict(UInt(0b0000001) => -2.00 - 0.25im, UInt(0b1001) => -3.0)
      @test typeof(ψ3) == SparseState{ComplexF64, UInt}
      @test typeof(ψ4) == SparseState{ComplexF64, UInt}
    end # intertype
  end # binary
end # SparseState