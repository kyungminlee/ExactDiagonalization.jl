using Test
using ExactDiagonalization
using StaticArrays

@testset "HSR" begin

  @testset "spin half" begin
    up = State("Up", +1)
    dn = State("Dn", -1)
    spinsite = Site([up, dn])
    @testset "constructor" begin
      hilbert_space = HilbertSpace([spinsite for i in 1:4])
      basis_list = UInt[0b0011, 0b0101, 0b0110, 0b1001, 0b1010, 0b1100] # for Sz=0 sector
      basis_lookup = FrozenSortedArrayIndex(basis_list)

      hsr1 = represent(hilbert_space, basis_list)
      hsr2 = represent(HilbertSpaceSector(hilbert_space, 0))
      hsr3 = HilbertSpaceRepresentation(hilbert_space, basis_list, basis_lookup)
      hsr4 = HilbertSpaceRepresentation(HilbertSpaceSector(hilbert_space, 0), basis_list, basis_lookup)
      @test hsr1 == hsr2
      @test hsr1 == hsr3
      @test hsr1 == hsr4
      @test basespace(hsr1) === hilbert_space
    end

    @testset "constructor-exceptions" begin
      hilbert_space = HilbertSpace([spinsite for i in 1:9])
      basis_list = UInt8[0b0101]
      @test_throws ArgumentError HilbertSpaceRepresentation(hilbert_space, basis_list, FrozenSortedArrayIndex(basis_list))
      @test_throws ArgumentError represent(hilbert_space, UInt8)
      @test_throws ArgumentError HilbertSpaceRepresentation(HilbertSpaceSector(hilbert_space, 0), basis_list, FrozenSortedArrayIndex(basis_list))
      @test_throws ArgumentError represent(HilbertSpaceSector(hilbert_space, 0), UInt8)
    end

    @testset "typetraits" begin
      hilbert_space = HilbertSpace([spinsite for i in 1:4])
      hsr = represent(hilbert_space, UInt32)
      @test scalartype(hsr) === Bool
      @test valtype(hsr) === Bool
      @test bintype(hsr) === UInt32
      @test scalartype(typeof(hsr)) === Bool
      @test valtype(typeof(hsr)) === Bool
      @test bintype(typeof(hsr)) === UInt32
    end

    @testset "properties" begin
      hilbert_space = HilbertSpace([spinsite for i in 1:4])
      hsr = represent(hilbert_space, UInt32)
      @test dimension(hsr) == 16
      @test bitwidth(hsr) == 4
    end

    @testset "validity" begin
      up = State("Up", 0)
      dn = State("Dn", 0)
      spinsite = Site([up, dn])
      hilbert_space = HilbertSpace([spinsite for i in 1:9])
      let # missing key
        basis_list = UInt[0x0]
        basis_lookup = FrozenSortedArrayIndex(UInt[0x1])
        hsr = HilbertSpaceRepresentation(hilbert_space, basis_list, basis_lookup)
        @test_throws KeyError ExactDiagonalization.checkvalidbasis(hsr)
      end
      let # wrong index
        basis_list = UInt[0x1]
        basis_lookup = FrozenSortedArrayIndex(UInt[0x0, 0x1])
        hsr = HilbertSpaceRepresentation(hilbert_space, basis_list, basis_lookup)
        @test_throws AssertionError ExactDiagonalization.checkvalidbasis(hsr)
      end
    end
  end # testset spin half

  @testset "tJ" begin
    QN = Tuple{Int, Int}
    em = State("Em", ( 0, 0))  # charge and spin
    up = State("Up", ( 1, 1))
    dn = State("Dn", ( 1,-1))
    spin_site = Site([up, dn])
    site = Site([em, up, dn])
    hs = HilbertSpace([site, spin_site, spin_site]) # f s s
    sectors = quantum_number_sectors(hs)
    @test sectors == QN[(2, -2), (2, 0), (2, 2), (3, -3), (3, -1), (3, 1), (3, 3)]

    @testset "represent" begin
      for represent in [represent_array, represent_dict]
        hsr_all = represent(hs)
        @test hsr_all.basis_list == UInt[
          0b0000, 0b0001, 0b0010,
          0b0100, 0b0101, 0b0110,
          0b1000, 0b1001, 0b1010,
          0b1100, 0b1101, 0b1110,
        ]
        @test all(ibasis == hsr_all.basis_lookup[basis] for (ibasis, basis) in enumerate(hsr_all.basis_list))
        @test dimension(hsr_all) == length(hsr_all.basis_list)
        @test hsr_all == represent(HilbertSpaceSector(hs, sectors))

        # empty QN
        @test represent(HilbertSpaceSector(hs, QN[])).basis_list == []

        qn_basis = Dict{QN, Vector{UInt}}(
                      ( 2,-2) => [0b1100],
                      ( 2, 0) => [0b0100, 0b1000],
                      ( 2, 2) => [0b0000],
                      ( 3,-3) => [0b1110],
                      ( 3,-1) => [0b0110, 0b1010, 0b1101], # UDD DUD DDU
                      ( 3, 1) => [0b0010, 0b0101, 0b1001], # UUD UDU DUU
                      ( 3, 3) => [0b0001], # UUU
                  )

        # Test each sector explicitly
        for (qn, basis_list) in qn_basis
          hsr = represent(HilbertSpaceSector(hs, qn))
          @test hsr.basis_list == basis_list
          @test all(ibasis == hsr.basis_lookup[basis] for (ibasis, basis) in enumerate(hsr.basis_list))
          @test dimension(hsr) == length(hsr.basis_list)
          hsr2 = represent(hs, basis_list)
          @test hsr == hsr2
        end

        for qn in sectors
          hsr = represent(HilbertSpaceSector(hs, qn))
          @test hsr_all != hsr
          @test all(ibasis == hsr.basis_lookup[basis] for (ibasis, basis) in enumerate(hsr.basis_list))
          @test dimension(hsr) == length(hsr.basis_list)
          @test hsr == represent(hs, hsr.basis_list)
        end

        @test represent(HilbertSpaceSector(hs, QN[( 3, 3), ( 3, 1)])).basis_list == [0b0001, 0b0010, 0b0101, 0b1001]
        @test represent(HilbertSpaceSector(hs, Set{QN}([( 3, 3), ( 3, 1)]))).basis_list == [0b0001, 0b0010, 0b0101, 0b1001]
        @test represent(hs, UInt[0b0001, 0b0010, 0b0101, 0b1001]).basis_list == UInt[0b0001, 0b0010, 0b0101, 0b1001]
        # order
        @test represent(hs, UInt[0b0101, 0b0001, 0b0010,0b1001]).basis_list == UInt[0b0001, 0b0010, 0b0101, 0b1001]
      end
    end
  end # tJ
end
