using Test
using ExactDiagonalization

@testset "simplify" begin
  function pauli_matrix(hs::AbstractHilbertSpace, isite ::Integer, j ::Symbol)
    if j == :x
      return pure_operator(hs, isite, 1, 2, 1.0; dtype=UInt) + pure_operator(hs, isite, 2, 1, 1.0; dtype=UInt)
    elseif j == :y
      return pure_operator(hs, isite, 1, 2, -1.0im; dtype=UInt) + pure_operator(hs, isite, 2, 1, 1.0im; dtype=UInt)
    elseif j == :z
      return pure_operator(hs, isite, 1, 1, 1.0; dtype=UInt) + pure_operator(hs, isite, 2, 2, -1.0; dtype=UInt)
    elseif j == :+
      return pure_operator(hs, isite, 1, 2, 1.0; dtype=UInt)
    elseif j == :-
      return pure_operator(hs, isite, 2, 1, 1.0; dtype=UInt)
    else
      throw(ArgumentError("pauli matrix of type $(j) not supported"))
    end
  end
  
  @testset "spinhalf" begin
    QN = Int
    up = State("Up", QN( 1))
    dn = State("Dn", QN(-1))
    spin_site = Site([up, dn])
    n = 4
    hs = AbstractHilbertSpace([spin_site, spin_site, spin_site, spin_site])

    σ(i::Integer, j::Symbol) = pauli_matrix(hs, i, j)

    heisenberg = sum(2 * σ(i, :+) * σ( mod(i, n) + 1, :-) for i in 1:n)
    heisenberg += sum(2 * σ(i, :-) * σ( mod(i, n) + 1, :+) for i in 1:n)
    heisenberg += sum(σ(i, :z) * σ( mod(i, n) + 1, :z) for i in 1:n)
    heisenberg2 = sum(σ(i, j) * σ(mod(i,n)+1, j) for i in 1:n, j in [:x, :y, :z])
    heisenberg3 = simplify(heisenberg)
    heisenberg4 = simplify(heisenberg2)
    
    @test heisenberg != heisenberg2
    @test length(heisenberg.terms) < length(heisenberg2.terms)
    @test Set(heisenberg.terms) != Set(heisenberg2.terms)
    @test Set(heisenberg.terms) == Set(heisenberg3.terms)
    @test Set(heisenberg.terms) == Set(heisenberg4.terms)
    @test heisenberg3 == heisenberg4

    @test simplify(heisenberg3 + heisenberg4) == simplify(2 * heisenberg3)
    @test simplify(heisenberg3 - heisenberg3) == NullOperator()
  end
end