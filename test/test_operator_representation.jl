using Test
using ExactDiagonalization

using StaticArrays
using LinearAlgebra


function pauli_matrix(hs::HilbertSpace, isite ::Integer, j ::Symbol)
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
end;



@testset "OperatorRepresentation" begin
  @testset "spin half" begin
    up = State("Up", +1)
    dn = State("Dn", -1)
    spinsite = Site([up, dn])
    n_sites = 4
    hs = HilbertSpace([spinsite for i in 1:n_sites])
    hsr = represent(hs)

    σ = Dict((isite, j) => pauli_matrix(hs, isite, j)
               for isite in 1:n_sites, j in [:x, :y, :z, :+, :-]);

    @testset "constructor" begin
      opr1 = OperatorRepresentation(hsr, σ[1, :x])
      opr2 = OperatorRepresentation{typeof(hsr), typeof(σ[1,:x])}(hsr, σ[1, :x])
      @show opr1.hilbert_space_representation == opr2.hilbert_space_representation
      @show opr1.operator == opr2.operator
      #@test opr1 == opr2
    end

  end
end
