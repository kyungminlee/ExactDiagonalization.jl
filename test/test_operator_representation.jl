using Test
using ExactDiagonalization

using StaticArrays
using LinearAlgebra
using SparseArrays

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
      opr3 = represent(hsr, σ[1, :x])
      @test opr1.hilbert_space_representation == hsr
      @test opr1.operator == σ[1, :x]
      @test opr2.hilbert_space_representation == hsr
      @test opr2.operator == σ[1, :x]
      @test opr1 == opr2
      @test opr1 == opr3
      @test get_space(opr1) === hsr
    end

    @testset "typetraits" begin
      opr1 = OperatorRepresentation(hsr, σ[1, :x])
      opr2 = OperatorRepresentation(hsr, σ[1, :y])
      OR = typeof(opr1)
      @test spacetype(OR) == typeof(hsr)
      @test operatortype(OR) == typeof(σ[1, :x])
      @test spacetype(opr1) == typeof(hsr)
      @test operatortype(opr1) == typeof(σ[1, :x])

      @test eltype(opr1) === Int
      @test eltype(opr2) === Complex{Int}
    end

    @testset "size" begin
      opr = OperatorRepresentation(hsr, σ[1, :x])
      dim = dimension(hsr)
      @test size(opr) == (dim, dim)
      @test size(opr, 1) == dim
      @test size(opr, 2) == dim
      @test_throws BoundsError size(opr, 3)
    end


    @testset "unary operator" begin
      op = σ[1, :x]
      opr = OperatorRepresentation(hsr, op)
      @test +opr == opr
      @test -opr == represent(hsr, -op)
    end

    @testset "binary operator" begin
      op1 = σ[1, :x]
      op2 = σ[3, :z]
      opr1 = OperatorRepresentation(hsr, op1)
      opr2 = OperatorRepresentation(hsr, op2)
      @test opr1 + opr2 == represent(hsr, op1 + op2)
      @test opr1 - opr2 == represent(hsr, op1 - op2)
      @test opr1 * opr2 == represent(hsr, op1 * op2)
    end

    @testset "iterator" begin
      opr = OperatorRepresentation(hsr, σ[2, :+]) # non-Hermitian
      dim = dimension(hsr)
      σ₊ = [0 1; 0 0]
      σ₀ = [1 0; 0 1]
      H0 = kron(σ₀, σ₀, σ₊, σ₀)
      opr_s = sparse(opr)
      @test opr_s == H0
      @test opr[:,:] == opr_s
      @test isa(opr_s, SparseMatrixCSC)
      for i in 1:dim
        @test get_row(opr, i) == H0[i, :]
        @test get_column(opr, i) == H0[:, i]
        @test get_row(opr, i) != get_column(opr, i) # all rows and columns are different

        @test get_row(opr, i) == opr[i, :]
        @test get_column(opr, i) == opr[:, i]
      end
    end # testset iterator

  end
end
