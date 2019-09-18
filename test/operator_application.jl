using Test
using ExactDiagonalization

function pauli_matrix(hs::AbstractHilbertSpace, isite ::Integer, j ::Symbol)
  if j == :x
    return pure_operator(hs, isite, 1, 2, 1; dtype=UInt) + pure_operator(hs, isite, 2, 1, 1; dtype=UInt)
  elseif j == :y
    return pure_operator(hs, isite, 1, 2, -im; dtype=UInt) + pure_operator(hs, isite, 2, 1, im; dtype=UInt)
  elseif j == :z
    return pure_operator(hs, isite, 1, 1, 1; dtype=UInt) + pure_operator(hs, isite, 2, 2, -1; dtype=UInt)
  elseif j == :+
    return pure_operator(hs, isite, 1, 2, 1; dtype=UInt)
  elseif j == :-
    return pure_operator(hs, isite, 2, 1, 1; dtype=UInt)
  else
    throw(ArgumentError("pauli matrix of type $(j) not supported"))
  end
end

@testset "apply" begin
  QN = Int
  up = State("Up", QN( 1))
  dn = State("Dn", QN(-1))
  spin_site = Site([up, dn])
  hs = AbstractHilbertSpace([spin_site, spin_site, spin_site, spin_site])
  hs2 = AbstractHilbertSpace([spin_site, spin_site, spin_site,])
  chs = concretize(hs, 0)
  psi = SparseState{Float64, UInt}(hs, UInt(0b0011) => 2.0, UInt(0b0101) => 10.0)
  
  σ(i::Integer, j::Symbol) = pauli_matrix(hs, i, j)
  
  @testset "hilbert" begin
    psi2 = SparseState{Float64, UInt}(hs2, UInt(0b0011) => 2.0, UInt(0b0101) => 10.0)
    @test_throws ArgumentError apply(σ(1, :+), psi2)
    @test_throws ArgumentError apply(σ(1, :x), psi2)
  end
  
  @test apply(σ(1, :+), psi) == SparseState{Float64, UInt}(hs)
  @test apply(σ(1, :-), psi) == SparseState{Float64, UInt}(hs, UInt(0b0010) => 2.0, UInt(0b0100) => 10.0)

  @test apply(σ(2, :+), psi) == SparseState{Float64, UInt}(hs, UInt(0b0111) => 10.0)
  @test apply(σ(2, :-), psi) == SparseState{Float64, UInt}(hs, UInt(0b0001) => 2.0)

  for i1 in 1:4, j1 in [:x, :y, :z, :+, :-]
    for i2 in 1:4, j2 in [:x, :y, :z, :+, :-]
      ϕ1 = apply(σ(i2, j2), apply(σ(i1, j1), psi))
      ϕ2 = apply(σ(i1, j1) * σ(i2,j2), psi)
      @test ϕ1 == ϕ2
    end
  end
end
# @testset "apply" begin
#   QN = SVector{2, Int}
#   em = State("Em", QN( 0, 0))  # charge and spin
#   up = State("Up", QN( 1, 1))
#   dn = State("Dn", QN( 1,-1))
#   spin_site = Site([up, dn])
#   site = Site([em, up, dn])
#   hs = AbstractHilbertSpace([site, spin_site, spin_site]) # f s s

#   chs = concretize(hs, QN([3, 1]))
# end