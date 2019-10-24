using Test
using ExactDiagonalization

@testset "util" begin
  @testset "make_bitmask" begin
    for dtype in [UInt8, UInt16, UInt32, UInt64, UInt128]
      @test isa(make_bitmask(4; dtype=dtype), dtype)
      @test isa(make_bitmask(4, 2; dtype=dtype), dtype)
    end
    @test make_bitmask(5) == 0b11111
    @test make_bitmask(5, 3) == 0b11000
  end

  #=
  @testset "bitcount" begin
    @test ExactDiagonalization.bitcount(UInt(0b0110101)) == 4
    @test ExactDiagonalization.bitcount(UInt(0b0100001)) == 2
    for dtype in [UInt8, UInt16, UInt32, UInt64, UInt128]
      values = rand(dtype, 64)
      for v in values
        c1 = count(x -> x == '1', string(v, base=2))
        c2 = ExactDiagonalization.bitcount(v)
        @test c1 == c2
      end
    end
  end
  =#

  @testset "merge_vec" begin
    let
      a = Int[1,5,7,8]
      b = Int[2,3,4,6,9]
      c = ExactDiagonalization.merge_vec(a, b)
      @test c == [1,2,3,4,5,6,7,8,9]
    end

    let
      a = Int[1,2,3]
      b = Int[4,5,6,7]
      c = ExactDiagonalization.merge_vec(a, b)
      @test c == [1,2,3,4,5,6,7]
    end

    let
      a = Int[4,5,6,7]
      b = Int[1,2,3]
      c = ExactDiagonalization.merge_vec(a, b)
      @test c == [1,2,3,4,5,6,7]
    end

    let
      a = Int[1,2,3]
      b = Int[4,5,6,7]
      c = ExactDiagonalization.merge_vec(a, b)
      @test c == [1,2,3,4,5,6,7]
    end

    let
      a = Int[1,2,3,4]
      b = Int[5,6,7]
      c = ExactDiagonalization.merge_vec(a, b)
      @test c == [1,2,3,4,5,6,7]
    end

    let
      a = Int[1,2,3]
      b = Int[4,5,6,7]
      c = ExactDiagonalization.merge_vec(a, b)
      @test c == [1,2,3,4,5,6,7]
    end

    let
      a = Int[1,1,1,1]
      b = Int[1,1,1,1,1]
      c = ExactDiagonalization.merge_vec(a, b)
      @test c == [1,1,1,1,1,1,1,1,1]
    end

  end

  @testset "choptol!" begin
    let
      d = Dict("A" => 0.0, "B" => 1E-9)
      choptol!(d, 1E-12)
      @test d == Dict("B"=>1E-9)
    end
    let
      d = Dict("A" => 0.0, "B" => 1E-9)
      choptol!(d, 1E-6)
      @test d == Dict()
    end
  end

end



@testset "FrozenSortedArray" begin
  @test_throws ArgumentError FrozenSortedArrayIndex{UInt}(UInt[0x2, 0x1, 0x3])
  @test_throws ArgumentError FrozenSortedArrayIndex{UInt}(UInt[0x1, 0x1, 0x3])
  fsa = FrozenSortedArrayIndex{UInt}(UInt[0x2, 0x4, 0x6])
  @test fsa[0x2] == 1
  @test fsa[0x4] == 2
  @test fsa[0x6] == 3
  @test collect(fsa) == [0x2 =>1, 0x4=>2, 0x6=>3]
  @test length(fsa) == 3
  @test_throws KeyError fsa[0x1]
  @test haskey(fsa, 0x2)
  @test !haskey(fsa, 0x3)
  @test get(fsa, 0x1, -1) == -1
  @test keys(fsa) == UInt[0x2, 0x4, 0x6]
end
