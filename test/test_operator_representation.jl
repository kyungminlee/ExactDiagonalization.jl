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
      @test bintype(OR) == UInt
      @test bintype(opr1) == UInt

      @test eltype(opr1) === Int
      @test eltype(opr2) === Complex{Int}
      @test eltype(typeof(opr1)) === Int
      @test eltype(typeof(opr2)) === Complex{Int}
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

      hs3 = HilbertSpace([spinsite, spinsite])
      op3 = pure_operator(hs3, 1, 1, 1)
      hsr3 = represent(hs3)
      opr3 = represent(hsr3, op3)
      @test_throws ArgumentError opr1 + opr3
      @test_throws ArgumentError opr1 - opr3
      @test_throws ArgumentError opr1 * opr3
    end


    @testset "iterator" begin
      hs2 = HilbertSpace([spinsite, spinsite])
      op2 = pure_operator(hs2, 2, 1, 2)*2.0 + pure_operator(hs2, 2, 2, 1)*3.0 # 2 σ⁻₂ + 3 σ⁺₂
      hsr2 = represent(hs, UInt[0b00, 0b01])
      opr2 = represent(hsr2, op2)
      @test collect(get_row_iterator(opr2, 1; include_all=true)) == [-1 => 2.0]
      @test collect(get_row_iterator(opr2, 2; include_all=true)) == [-1 => 2.0]
      @test collect(get_column_iterator(opr2, 1; include_all=true)) == [-1 => 3.0]
      @test collect(get_column_iterator(opr2, 2; include_all=true)) == [-1 => 3.0]
      @test collect(get_row_iterator(opr2, 1)) == []
      @test collect(get_row_iterator(opr2, 2)) == []
      @test collect(get_column_iterator(opr2, 1)) == []
      @test collect(get_column_iterator(opr2, 2)) == []
    end

    @testset "get" begin
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

    @testset "apply" begin # complete space
      hsr_0 = represent(hs)
      dim = dimension(hsr_0)

      @testset "type" begin
        opr_x = represent(hsr_0, σ[2, :x])

        @test typeof(opr_x * rand(Int, dim)) === Vector{Int}
        @test typeof(opr_x * rand(Float64, dim)) === Vector{Float64}
        @test typeof(opr_x * rand(ComplexF64, dim)) === Vector{ComplexF64}

        opr_y = represent(hsr_0, σ[2, :y]) # operator has eltype Complex{Int}
        @test typeof(opr_y * rand(Int, dim)) === Vector{Complex{Int}}
        @test typeof(opr_y * rand(Float64, dim)) === Vector{ComplexF64}
        @test typeof(opr_y * rand(ComplexF64, dim)) === Vector{ComplexF64}
      end

      @testset "value" begin
        opr = represent(hsr_0, σ[2, :+])
        
        σ₊ = [0 1; 0 0]
        σ₀ = [1 0; 0 1]
        
        op_dense = kron(σ₀, σ₀, σ₊, σ₀)
        state = rand(ComplexF64, dim)
        
        out1 = zeros(ComplexF64, dim)

        # matrix * columnvector
        out0 = op_dense * state
        apply_unsafe!(out1, opr, state)
        out2 = opr * state
        @test isapprox(out0, out1, atol=1E-6)
        @test isapprox(out0, out2, atol=1E-6)
        
        # add to the previous (do not overwrite)
        apply_unsafe!(out1, opr, state)
        @test !isapprox(out0, out1, atol=1E-6)
        
        # rowvector * matrix
        out1[:] .= zero(ComplexF64)
        
        out0 = transpose( transpose(state) * op_dense )
        apply_unsafe!(out1, state, opr)
        out2 = state * opr
        @test isapprox(out0, out1, atol=1E-6)
        @test isapprox(out0, out2, atol=1E-6)

        # add to the previous (do not overwrite)
        apply_unsafe!(out1, state, opr)
        @test !isapprox(out0, out1, atol=1E-6)
      end
    end

  end # testset spin half
end # testset OperatorRepresentation
