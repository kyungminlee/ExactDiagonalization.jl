using Test
using LinearAlgebra
using ExactDiagonalization
using StaticArrays

@testset "SparseState" begin
  QN = SVector{2, Int}
  em = State("Em", QN( 0, 0))  # charge and spin
  up = State("Up", QN( 1, 1))
  dn = State("Dn", QN( 1,-1))
  spin_site = Site([up, dn])
  site = Site([em, up, dn])
  hs = HilbertSpace([site, site, spin_site, site])
  hs2 = HilbertSpace([site, spin_site, site, site])

  @testset "constructor" begin
    ψ1 = SparseState{ComplexF64, UInt}(hs, UInt(0b0010001))
    @test ψ1.hilbert_space == hs

    @test ψ1.components == Dict(UInt(0x11) => 1.0 + 0.0im)
    @test ψ1[UInt(0x11)] == 1.0 + 0.0im
    @test ψ1.components == Dict(UInt(0x11) => 1.0 + 0.0im)

    ψ2 = SparseState{ComplexF64, UInt}(hs)
    ψ2[UInt(0b0010001)] = 1.0
    @test ψ2.components == Dict(UInt(0x11) => 1.0 + 0.0im)

    ψ3 = SparseState{ComplexF64, UInt}(hs, UInt(0b0010001) => 2.0)
    @test ψ3.components == Dict(UInt(0x11) => 2.0 + 0.0im)

    ψ4 = SparseState{ComplexF64, UInt}(hs, UInt(0b0000001) => 2.0, UInt(0b0010001) => 3.0 )
    @test ψ4.components == Dict(UInt(0x11) => 3.0 + 0.0im, UInt(0x1) => 2.0 + 0.0im)

    ψ5 = SparseState{ComplexF64, UInt}(hs, Dict(UInt(0b0000001) => 2.0, UInt(0b0010001) => 3.0 ))
    @test ψ5.components == Dict(UInt(0x11) => 3.0 + 0.0im, UInt(0x1) => 2.0 + 0.0im)
  end

  @testset "type" begin
    ψ1 = SparseState{ComplexF64, UInt32}(hs, UInt32(0b0010001))
    @test scalartype(ψ1) === ComplexF64
    @test scalartype(typeof(ψ1)) === ComplexF64
    @test bintype(ψ1) === UInt32
    @test bintype(typeof(ψ1)) === UInt32
    @test eltype(ψ1) === Pair{UInt32, ComplexF64}
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
  end

  @testset "isapprox" begin
    ψ1 = SparseState{ComplexF64, UInt}(hs, UInt(0b0010001) => 1.0 + 0.0im)
    @test isapprox(ψ1, ψ1)
    let
      ψ2 = SparseState{ComplexF64, UInt}(hs2, UInt(0b0010001) => 1.0 + 0.0im)
      @test !isapprox(ψ1, ψ2)
    end
    let
      ψ2 = copy(ψ1)
      ψ2[UInt(0b0010001)] += 1E-13
      @test isapprox(ψ1, ψ2)
      @test ψ1 != ψ2
    end
    let
      ψ2 = copy(ψ1)
      ψ2[UInt(0b0000000)] = 0
      @test isapprox(ψ1, ψ2)
      @test ψ1 != ψ2
    end
  end

  @testset "iterate" begin
    ψ = SparseState{ComplexF64, UInt}(hs, Dict(UInt(0b0000001) => 2.0, UInt(0b0010001) => 3.0 ))
    @test Dict(collect(ψ)) == Dict(UInt(0b0000001) => 2.0+0.0im, UInt(0b0010001) => 3.0+0.0im)
    @test Set([k for (k, v) in ψ]) == Set(UInt[0b0000001, 0b0010001])
    @test Set([v for (k, v) in ψ]) == Set([2.0+0.0im, 3.0+0.0im])
    @test length(ψ) == 2
  end

  @testset "choptol!" begin
    let
      d = SparseState{ComplexF64, UInt}(hs, UInt(0b0000001) => 0.0, UInt(0b0010001) => 1E-9)
      choptol!(d, 1E-12)
      @test d.components == Dict(UInt(0b0010001)=>1E-9 + 0.0im)
    end
    let
      d = SparseState{ComplexF64, UInt}(hs, UInt(0b0000001) => 0.0, UInt(0b0010001) => 1E-9)
      choptol!(d, 1E-6)
      @test d.components == Dict()
    end
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

    @testset "normalize" begin
      ψ0 = SparseState{ComplexF64, UInt}(hs)
      @test norm(ψ0) == 0.0

      ψ = SparseState{ComplexF64, UInt}(hs, UInt(0b0000) => 3.0 + 4.0im)
      @test isapprox(norm(ψ), 5.0)
      ψ2 = normalize(ψ)
      @test isapprox(ψ2, SparseState{ComplexF64, UInt}(hs, UInt(0b0000) => 0.6 + 0.8im))

      @test isapprox(normalize!(ψ), ψ2)
      @test isapprox(ψ, ψ2)
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
      @test (4 \ ψ1).components == Dict(UInt(0x001) => 0.5)
      @test (4.0 \ ψ1).components == Dict(UInt(0x001) => 0.5)
      @test ((4.0+0.0im) \ ψ1 ).components == Dict(UInt(0x001) => 0.5 + 0.0im)
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
