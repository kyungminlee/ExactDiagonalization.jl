using Test
using ExactDiagonalization

using StaticArrays

@testset "NullOperator" begin
  nop = NullOperator()
  nop2 = NullOperator()

  @test isa(nop, NullOperator)
  @test nop == nop2

  @testset "unary" begin
    @test -nop == nop
    @test +nop == nop
    @test real(nop) == nop
    @test imag(nop) == nop
    @test conj(nop) == nop
    @test transpose(nop) == nop
  end

  @testset "binary" begin
    @test nop + nop == nop
    @test nop - nop == nop
    @test nop * nop == nop

    @test nop * 2 == nop
    @test nop * 2.0 == nop
    @test nop * (2.0 + 1.0im) == nop

    @test 2 * nop == nop
    @test 2.0 * nop == nop
    @test (2.0 + 1.0im) * nop == nop
  end
end

@testset "PureOperator" begin
  QN = Int
  up = State("Up", 1)
  dn = State{QN}("Dn",-1)
  spin_site = Site{QN}([up, dn])

  hs = AbstractHilbertSpace{QN}([spin_site for i in 1:4])
  hs2 = AbstractHilbertSpace{QN}([spin_site, spin_site])
  nop = NullOperator()

  @testset "type" begin
    types = [Int32, Int64, Float32, Float64, ComplexF32, ComplexF64]
    for t1 in types, t2 in types
      pop = PureOperator{t1, UInt}(hs, 0b0010, 0b0000, 0b0000, t1(2))
      pop.amplitude == t1(2)
    end
  end

  @testset "equality" begin
    pop1 = PureOperator{ComplexF64, UInt}(hs, 0b0010, 0b0000, 0b0000, 2.0 + 3.0im)
    pop2 = PureOperator{ComplexF64, UInt}(hs, 0b0010, 0b0000, 0b0000, 2.0 + 3.0im)
    pop3 = PureOperator{ComplexF64, UInt}(hs, 0b0010, 0b0000, 0b0000, 2.0 + 1.0im)
    @test pop1 == pop2
    @test pop1 != pop3
  end

  @testset "equality-intertype" begin
    pop1 = PureOperator{Int, UInt}(hs, 0b0010, 0b0000, 0b0000, 2)
    pop2 = PureOperator{Float64, UInt}(hs, 0b0010, 0b0000, 0b0000, 2.0)
    pop3 = PureOperator{ComplexF64, UInt}(hs, 0b0010, 0b0000, 0b0000, 2.0 + 0.0im)

    # same amplitude in different types
    @test pop1 == pop2
    @test pop1 == pop3
    @test pop2 == pop3
    
    # different bitfields
    @test pop2 != PureOperator{Float64, UInt}(hs, 0b0001, 0b0000, 0b0000, 2.0)
    @test pop2 != PureOperator{Float64, UInt}(hs, 0b0010, 0b0010, 0b0000, 2.0)
    @test pop2 != PureOperator{Float64, UInt}(hs, 0b0010, 0b0000, 0b0010, 2.0)

    # different hilbert space
    @test pop2 != PureOperator{Float64, UInt}(hs2, 0b0010, 0b0000, 0b0000, 2.0)
  end

  @testset "unary" begin
    @testset "real" begin
      pop1 = PureOperator{Float64, UInt}(hs, 0b0010, 0b0010, 0b0000, 2.0)
      pop2 = PureOperator{Float64, UInt}(hs, 0b0010, 0b0010, 0b0000, 2.0)
      @test pop1 == pop2
      @test (+pop1).amplitude == 2.0
      @test (-pop1).amplitude == -2.0
      @test (real(pop1)) == pop1
      @test (imag(pop1)) != NullOperator()
      @test (imag(pop1)) == PureOperator{Float64, UInt}(hs, 0b0010, 0b0010, 0b0000, 0.0)
      @test conj(pop1) == pop1
      @test transpose(pop1) == PureOperator{Float64, UInt}(hs, 0b0010, 0b0000, 0b0010, 2.0)
    end

    @testset "complex" begin
      pop1 = PureOperator{ComplexF64, UInt}(hs, 0b0010, 0b0010, 0b0000, 2.0 + 3.0im)
      pop2 = PureOperator{ComplexF64, UInt}(hs, 0b0010, 0b0010, 0b0000, 2.0 + 3.0im)
      @test pop1 == pop2
      @test (+pop1).amplitude == 2.0 + 3.0im
      @test (-pop1).amplitude == -2.0 - 3.0im
      @test (real(pop1)).amplitude == 2.0
      @test (imag(pop1)).amplitude == 3.0

      @test real(pop1) == PureOperator{Float64, UInt}(hs, 0b0010, 0b0010, 0b0000, 2.0)
      @test imag(pop1) == PureOperator{Float64, UInt}(hs, 0b0010, 0b0010, 0b0000, 3.0)
      @test conj(pop1) == PureOperator{ComplexF64, UInt}(hs, 0b0010, 0b0010, 0b0000, 2.0 - 3.0im)
      @test transpose(pop1) == PureOperator{ComplexF64, UInt}(hs, 0b0010, 0b0000, 0b0010, 2.0 + 3.0im)
    end
  end

  @testset "binary" begin
    @testset "product" begin
      @testset "types" begin
        types = [Int, Float64, ComplexF64]
        for t1 in types, t2 in types
          pop1 = PureOperator{t1, UInt}(hs, 0b0010, 0b0000, 0b0000, t1(2))
          pop2 = PureOperator{t2, UInt}(hs, 0b0010, 0b0000, 0b0000, t2(2))
          t3 = promote_type(t1, t2)
          pop3 = pop1 * pop2
          pop4 = pop2 * pop1
          @test pop1 != NullOperator()
          @test pop2 != NullOperator()
          @test pop3 != NullOperator()
          @test pop4 != NullOperator()
          @test isa(pop3, PureOperator{t3, UInt})
          @test isa(pop4, PureOperator{t3, UInt})
        end
      end

      @testset "nullop" begin
        pop = PureOperator{Float64, UInt}(hs, 0b0010, 0b0000, 0b0010, 2.0)
        nop = NullOperator()

        @test pop * nop == nop
        @test nop * pop == nop
        @test pop + nop == pop
        @test nop + pop == pop
      end

      @testset "disjoint" begin
        pop1 = PureOperator{Float64, UInt}(hs, 0b0010, 0b0000, 0b0010, 2.0)
        pop2 = PureOperator{Float64, UInt}(hs, 0b1000, 0b1000, 0b0000, 3.0)        
        @test pop1 * pop2 == PureOperator{Float64, UInt}(hs, 0b1010, 0b1000, 0b0010, 6.0)
      end

      @testset "compatible" begin
        pop1 = PureOperator{Float64, UInt}(hs, 0b0010, 0b0000, 0b0010, 2.0)
        pop2 = PureOperator{Float64, UInt}(hs, 0b1010, 0b1010, 0b1000, 3.0)
        @test pop1 * pop2 == PureOperator{Float64, UInt}(hs, 0b1010, 0b1000, 0b1000, 6.0)
      end

      @testset "conflict" begin
        pop1 = PureOperator{Float64, UInt}(hs, 0b0010, 0b0000, 0b0010, 2.0)
        pop2 = PureOperator{Float64, UInt}(hs, 0b0010, 0b0000, 0b0000, 3.0)
        @test pop1 * pop2 == NullOperator()
      end

    end
  end
end

@testset "SumOperator" begin
  QN = Int
  up = State("Up", 1)
  dn = State{QN}("Dn",-1)
  spin_site = Site{QN}([up, dn])

  hs = AbstractHilbertSpace{QN}([spin_site for i in 1:4])
  hs2 = AbstractHilbertSpace{QN}([spin_site, spin_site])
  nop = NullOperator()

  @testset "constructor" begin
    pop1 = PureOperator{Float64, UInt}(hs, 0b0010, 0b0000, 0b0000, 2.0)
    pop2 = PureOperator{Int, UInt}(hs, 0b0001, 0b0001, 0b0001, 3)
    sop = SumOperator{Float64, UInt}(hs, [pop1, pop2])
    @test sop.hilbert_space == pop1.hilbert_space
    @test sop.terms[1] == pop1
    @test sop.terms[2] == pop2

    # empty
    SumOperator{ComplexF64, UInt}(hs, [])
    
    pop3 = PureOperator{ComplexF64, UInt}(hs, 0b0010, 0b0000, 0b0010, 3.0 + 4.0im)
    @test_throws InexactError SumOperator{Float64, UInt}(hs, [pop1, pop3])
    
    # hilbert space mismatch
    pop4 = PureOperator{Float64, UInt}(hs2, 0b0001, 0b0001, 0b0001, 3.0)
    @test_throws ArgumentError SumOperator{Float64, UInt}(hs, [pop1, pop4])
  end
  
  @testset "equality" begin
    pop1 = PureOperator{Float64, UInt}(hs, 0b0010, 0b0000, 0b0000, 2.0)
    pop2 = PureOperator{Float64, UInt}(hs, 0b0010, 0b0000, 0b0010, 3.0)
    pop3 = PureOperator{Float64, UInt}(hs, 0b0010, 0b0000, 0b0000, 4.0)
    sop1 = SumOperator{Float64, UInt}(hs, [pop1, pop2])
    sop2 = SumOperator{Float64, UInt}(hs, [pop1, pop3])
    @test sop1 == sop1
    @test sop1 != sop2

    # compare across types
    pop4 = PureOperator{ComplexF64, UInt}(hs, 0b0010, 0b0000, 0b0010, 3.0+0.0im)
    sop3 = SumOperator{ComplexF64, UInt}(hs, [pop1, pop4])
    @test sop1 == sop3
  end
  
  @testset "unary" begin
    pop1 = PureOperator{Float64, UInt}(hs, 0b0010, 0b0000, 0b0000, 2.0)
    pop2 = PureOperator{ComplexF64, UInt}(hs, 0b0010, 0b0000, 0b0010, 3.0 + 4.0im)
    @test_throws InexactError SumOperator{Float64, UInt}(hs, [pop1, pop2])
    sop = SumOperator{ComplexF64, UInt}(hs, [pop1, pop2])
    
    @test (+sop).hilbert_space == sop.hilbert_space
    @test (-sop).hilbert_space == sop.hilbert_space
    @test real(sop).hilbert_space == sop.hilbert_space
    @test imag(sop).hilbert_space == sop.hilbert_space

    @test (+sop).terms == [pop1, pop2]
    @test (-sop).terms == [-pop1, -pop2]

    @test real(sop).terms == [ real(pop1), real(pop2)]
    @test imag(sop).terms == [ imag(pop1), imag(pop2)]

    @test conj(sop).terms == [ conj(pop1), conj(pop2)]
    @test transpose(sop).terms == [ transpose(pop1), transpose(pop2)]
  end

  @testset "binary" begin
    @testset "sum" begin
      pop1 = PureOperator{Float64, UInt}(hs, 0b0010, 0b0000, 0b0000, 2.0)
      pop2 = PureOperator{ComplexF64, UInt}(hs, 0b0010, 0b0000, 0b0010, 3.0 + 4.0im)
      @test pop1 + pop2 == SumOperator{ComplexF64, UInt}(hs, [pop1, pop2])
      @test pop1 - pop2 == SumOperator{ComplexF64, UInt}(hs, [pop1,-pop2])

      sop = SumOperator{ComplexF64, UInt}(hs, [pop1, pop2])
      pop3 = PureOperator{Float64, UInt}(hs, 0b0010, 0b0010, 0b0010, 5.0)
      @test sop + pop3 == SumOperator{ComplexF64, UInt}(hs, [pop1, pop2, pop3])
      @test sop - pop3 == SumOperator{ComplexF64, UInt}(hs, [pop1, pop2,-pop3])
      @test pop3 + sop == SumOperator{ComplexF64, UInt}(hs, [pop3, pop1, pop2])
      @test pop3 - sop == SumOperator{ComplexF64, UInt}(hs, [pop3,-pop1,-pop2])

      @test_throws ArgumentError sop + SumOperator{ComplexF64, UInt}(hs2, [])
      @test_throws ArgumentError sop - SumOperator{ComplexF64, UInt}(hs2, [])

      pop4 = PureOperator{ComplexF64, UInt}(hs, 0b0100, 0b0100, 0b0000, 6.0+1.0im)
      sop2 = SumOperator{ComplexF64, UInt}(hs, [pop3, pop4])
      @test sop + sop2 == SumOperator{ComplexF64, UInt}(hs, [pop1, pop2, pop3, pop4])
      @test sop - sop2 == SumOperator{ComplexF64, UInt}(hs, [pop1, pop2,-pop3,-pop4])
    end

    @testset "product" begin
      pop1 = PureOperator{Float64, UInt}(hs, 0b0010, 0b0000, 0b0000, 2.0)
      pop2 = PureOperator{Float64, UInt}(hs, 0b0010, 0b0000, 0b0010, 3.0)
      pop3 = PureOperator{ComplexF64, UInt}(hs, 0b0010, 0b0000, 0b0010, 3.0 + 4.0im)
      pop4 = PureOperator{ComplexF64, UInt}(hs, 0b0010, 0b0010, 0b0000, 5.0 + 6.0im)
      sop = pop1 + pop2

      pop5 = pop1 * pop3
      pop6 = pop2 * pop3
      @test pop6 == NullOperator()
      @test sop * pop3 == SumOperator{ComplexF64, UInt}(hs, [pop5])

      pop7 = pop4 * pop1
      pop8 = pop4 * pop2
      @test pop7 != NullOperator()
      @test pop8 != NullOperator()
      @test pop4 * sop == SumOperator{ComplexF64, UInt}(hs, [pop7, pop8])

      sop2 = pop1 + pop2 + pop4
      nonzeroterms = filter((x) -> !isa(x, NullOperator), 
                            [pop1*pop1,
                             pop2*pop1,
                             pop4*pop1,
                             pop1*pop2,
                             pop2*pop2,
                             pop4*pop2,
                             pop1*pop4,
                             pop2*pop4,
                             pop4*pop4,
                            ])
      @test length(nonzeroterms) == 5
      @test sop2 * sop2 == SumOperator{ComplexF64, UInt}(hs, nonzeroterms)
    end
  end
end

  # popA00 = PureOperator{Float64, UInt}(hs, 0b0010, 0b0000, 0b0000, 2.0)
  # popA01 = PureOperator{ComplexF64, UInt}(hs, 0b0010, 0b0000, 0b0010, 3.0+im)
  # popA10 = PureOperator{ComplexF64, UInt}(hs, 0b0010, 0b0010, 0b0000, 4.0+2im)
  # popA11 = PureOperator{ComplexF64, UInt}(hs, 0b0010, 0b0010, 0b0010, 5.0+0im)

  # popB00 = PureOperator{ComplexF64, UInt}(hs, 0b1000, 0b0000, 0b0000, 2.0+im)
  # popB01 = PureOperator{ComplexF64, UInt}(hs, 0b1000, 0b0000, 0b1000, 3.0+im)
  # popB10 = PureOperator{ComplexF64, UInt}(hs, 0b1000, 0b1000, 0b0000, 4.0+2im)
  # popB11 = PureOperator{ComplexF64, UInt}(hs, 0b1000, 0b1000, 0b1000, 5.0+0im)

  # -popA00
  # +popA00
  # real(popA00)
  # imag(popA00)

  # pop = popA00
  # sop = pop + pop

  # +sop
  # -sop
  # real(sop)
  # imag(sop)

  # nop + nop
  # nop + pop
  # nop + sop

  # pop + nop
  # pop + pop
  # pop + sop

  # sop + nop
  # sop + pop
  # sop + sop

