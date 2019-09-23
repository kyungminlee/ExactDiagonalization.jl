using Test
using ExactDiagonalization

using StaticArrays

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

  @testset "convert" begin
    ψr = SparseState{Float64, UInt}(hs, UInt(0b00100001) => 2.0)
    ψc = SparseState{ComplexF64, UInt}(hs, UInt(0b00100001) => 2.0 + 0.0im)
    ψc_r = convert(SparseState{Float64, UInt}, ψc)
    ψr_c = convert(SparseState{ComplexF64, UInt}, ψr)
    @test typeof(ψr) == typeof(ψc_r)
    @test ψr == ψc_r
    @test typeof(ψc) == typeof(ψr_c)
    @test ψc == ψr_c
    
    ψc2 = SparseState{ComplexF64, UInt}(hs, UInt(0b00100001) => 1.0 + 2.0im)
    @test_throws InexactError convert(SparseState{Float64, UInt}, ψc2)
  end

  @testset "unary" begin
    @testset "real" begin
      ψ1 = SparseState{Float64, UInt}(hs)
      @test isempty(ψ1)
      ψ1[UInt(0b0000001)] = 3.0
      @test !isempty(ψ1)
      
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
      @test isempty(ψ1)
      ψ1[UInt(0b0000001)] = 3.0 + im
      @test !isempty(ψ1)

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