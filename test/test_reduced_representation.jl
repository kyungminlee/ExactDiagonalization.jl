using Test
using ExactDiagonalization

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
end

@testset "RedRep" begin
  QN = Int

  # Test State and Site
  up = State("Up", 1)
  dn = State("Dn",-1)
  spin_site = Site([up, dn])
  n = 4
  hs = HilbertSpace(repeat([spin_site], n))
  σ = Dict( (isite, j) => pauli_matrix(hs, isite, j) for isite in 1:n, j in [:x, :y, :z, :+, :-])
  j1 = sum(σ[i, j] * σ[mod(i, n) + 1 , j] for i in 1:n, j in [:x, :y, :z])

  hsr = represent(HilbertSpaceSector(hs, 0))
  translation_group = TranslationGroup([Permutation([2,3,4,1])])
  
  @testset "RHSR" begin
    rhsr = symmetry_reduce(hsr, translation_group, [0//1])
    @test scalartype(rhsr) === ComplexF64
    @test scalartype(typeof(rhsr)) === ComplexF64
    @test bintype(rhsr) === UInt
    @test bintype(typeof(rhsr)) === UInt
  end

  @testset "ROR" begin
    rhsr = symmetry_reduce(hsr, translation_group, [0//1])
    j1_r = represent(rhsr, j1)
    @testset "typetraits" begin
      @test scalartype(j1) === Complex{Int}
      @test scalartype(j1_r) === ComplexF64
      @test scalartype(typeof(j1)) === Complex{Int}
      @test scalartype(typeof(j1_r)) === ComplexF64

      @test bintype(j1_r) === UInt
      @test bintype(typeof(j1_r)) === UInt
      
      @test spacetype(j1_r) === typeof(rhsr)
      @test operatortype(j1_r) === typeof(j1)
      @test spacetype(typeof(j1_r)) === typeof(rhsr)
      @test operatortype(typeof(j1_r)) === typeof(j1)
      @test get_space(j1_r) === rhsr
    end
  end # testset ROR


  @testset "ROR-ALL" begin
    tol = sqrt(eps(Float64))
    j1_mat = Matrix(sparse(represent(hsr, j1)))
    @test isapprox(j1_mat, adjoint(j1_mat); atol=tol)
    eigenvalues1 = eigvals(Hermitian(j1_mat))
    eigenvalues2 = Float64[]
    for k in translation_group.fractional_momenta
      rhsr = symmetry_reduce(hsr, translation_group, k)
      j1_redrep = represent(rhsr, j1)
      j1_redmat = Matrix(sparse(j1_redrep))
      @test isapprox(j1_redmat, adjoint(j1_redmat); atol=tol)
      append!(eigenvalues2, eigvals(Hermitian(j1_redmat)))
    end
    sort!(eigenvalues2)
    @test length(eigenvalues1) == length(eigenvalues2)
    @test isapprox(eigenvalues1, eigenvalues2; atol=tol)
  end # testset iterator

end # testset RedOp

@testset "RedOp-nontriv" begin
# TODO: at an angle
end

