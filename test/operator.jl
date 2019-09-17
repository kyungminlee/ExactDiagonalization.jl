using Test
using ExactDiagonalization

using StaticArrays

@testset "PureOperator" begin
  QN = Int
  up = State("Up", 1)
  dn = State{QN}("Dn",-1)
  spin_site = Site{QN}([up, dn])

  hs = AbstractHilbertSpace{QN}([spin_site for i in 1:4])

  nop = NullOperator()
  pop0   = PureOperator{Float64, UInt}(hs, 0b0010, 0b0000, 0b0000, 0.0)
  popA00 = PureOperator{Float64, UInt}(hs, 0b0010, 0b0000, 0b0000, 2.0)
  popA01 = PureOperator{ComplexF64, UInt}(hs, 0b0010, 0b0000, 0b0010, 3.0+im)
  popA10 = PureOperator{ComplexF64, UInt}(hs, 0b0010, 0b0010, 0b0000, 4.0+2im)
  popA11 = PureOperator{ComplexF64, UInt}(hs, 0b0010, 0b0010, 0b0010, 5.0+0im)

  popB00 = PureOperator{ComplexF64, UInt}(hs, 0b1000, 0b0000, 0b0000, 2.0+im)
  popB01 = PureOperator{ComplexF64, UInt}(hs, 0b1000, 0b0000, 0b1000, 3.0+im)
  popB10 = PureOperator{ComplexF64, UInt}(hs, 0b1000, 0b1000, 0b0000, 4.0+2im)
  popB11 = PureOperator{ComplexF64, UInt}(hs, 0b1000, 0b1000, 0b1000, 5.0+0im)

  -popA00
  +popA00
  real(popA00)
  imag(popA00)

  pop = popA00
  sop = pop + pop

  +sop
  -sop
  real(sop)
  imag(sop)

  nop + nop
  nop + pop
  nop + sop

  pop + nop
  pop + pop
  pop + sop

  sop + nop
  sop + pop
  sop + sop


end