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

end
