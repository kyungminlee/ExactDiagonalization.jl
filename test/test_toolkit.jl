using Test
using ExactDiagonalization

@testset "Toolkit" begin
  @testset "SpinHalf" begin
    n_sites = 4

    QN = Int
    up = State("Up", 1)
    dn = State("Dn",-1)
    spin_site = Site([up, dn])

    hs1 = HilbertSpace([spin_site for i in 1:n_sites])
    (hs2, pauli) = ExactDiagonalization.Toolkit.spin_half_system(n_sites)

    @test hs1 == hs2
    for i_site in 1:n_sites
      @test pauli(i_site, :x) == pure_operator(hs1, i_site, 1, 2, 1, UInt) + pure_operator(hs1, i_site, 2, 1, 1, UInt)
      @test pauli(i_site, :y) == pure_operator(hs1, i_site, 1, 2, -im, UInt) + pure_operator(hs1, i_site, 2, 1, +im, UInt)
      @test pauli(i_site, :z) == pure_operator(hs1, i_site, 1, 1, 1, UInt) + pure_operator(hs1, i_site, 2, 2, -1, UInt)
      @test pauli(i_site, :+) == pure_operator(hs1, i_site, 1, 2, 1, UInt)
      @test pauli(i_site, :-) == pure_operator(hs1, i_site, 2, 1, 1, UInt)
    end
    @test_throws ArgumentError pauli(1, :unknown)
  end

  @testset "SpinSystem" begin
    r2, r3, r6 = sqrt(2), sqrt(3), sqrt(6)

    spin_matrices = Dict(
      1//2 => Dict(
          :x => [0 1; 1 0] ./2,
          :y => [0 1; -1 0] ./(2im),
          :z => [1 0; 0 -1] ./2,
          :+ => [0 1; 0 0],
          :- => [0 0; 1 0],
        ),
      1//1 => Dict(
          :x => [0 1 0; 1 0 1; 0 1 0] ./ r2,
          :y => [0 1 0; -1 0 1; 0 -1 0] ./ (r2*1im),
          :z => [1 0 0; 0 0 0; 0 0 -1],
          :+ => [0 1 0; 0 0 1; 0 0 0] .* r2,
          :- => [0 0 0; 1 0 0; 0 1 0] .* r2,
        ),
      3//2 => Dict(
          :x => [0 r3 0 0; r3 0 2 0; 0 2 0 r3; 0 0 r3 0] ./ 2,
          :y => [0 r3 0 0;-r3 0 2 0; 0 -2 0 r3; 0 0 -r3 0] ./ (2im),
          :z => [3/2 0 0 0; 0 1/2 0 0; 0 0 -1/2 0; 0 0 0 -3/2],
          :+ => [0 r3 0 0; 0 0 2 0; 0 0 0 r3; 0 0 0 0],
          :- => [0 0 0 0; r3 0 0 0; 0 2 0 0; 0 0 r3 0],
        ),
      2//1 => Dict(
          :x => [0 2 0 0 0; 2 0 r6 0 0; 0 r6 0 r6 0; 0 0 r6 0 2; 0 0 0 2 0] ./ 2,
          :y => [0 2 0 0 0; -2 0 r6 0 0; 0 -r6 0 r6 0; 0 0 -r6 0 2; 0 0 0 -2 0] ./ (2im),
          :z => [2 0 0 0 0; 0 1 0 0 0; 0 0 0 0 0; 0 0 0 -1 0; 0 0 0 0 -2],
          :+ => [0 2 0 0 0; 0 0 r6 0 0; 0 0 0 r6 0; 0 0 0 0 2; 0 0 0 0 0],
          :- => [0 0 0 0 0; 2 0 0 0 0; 0 r6 0 0 0; 0 0 r6 0 0; 0 0 0 2 0],
        ),
    )

    @testset "SingleSite" begin
      for S in [1//2, 1//1, 3//2, 2//1]
        (hs, spin) = ExactDiagonalization.Toolkit.spin_system(1, S)
        hsr = represent(hs)
        for μ in [:x, :y, :z, :+, :-]
          @test isapprox(spin_matrices[S][μ], Matrix(represent(hsr, spin(1, μ))))
        end
      end # for S
    end # SingleSite
  end

  @testset "product_state" begin
    (hs, pauli) = ExactDiagonalization.Toolkit.spin_half_system(4)
    @test_throws ArgumentError ExactDiagonalization.Toolkit.product_state(hs, [[1.0, 0.0], [0.0, 0.0]]) # too few
    @test_throws ArgumentError ExactDiagonalization.Toolkit.product_state(hs, [[1.0, 0.0], [0.0, 1.0], [0.0, 1.0], [0.0, 1.0], [0.0, 1.0]]) # too many
    @test_throws ArgumentError ExactDiagonalization.Toolkit.product_state(hs, [[1.0, 0.0], [0.0, 1.0], [0.0, 1.0, 0.0], [0.0, 1.0]])

    local_states = [[1.0, 0.0], [1.0 + 10.0im, 2.0], [1.0, 0.0], [0.0, 1.0]]
    @testset "HilbertSpace" begin
      psi1 = ExactDiagonalization.Toolkit.product_state(hs, local_states)
      psi2 = SparseState{ComplexF64, UInt}(0b1000 => 1.0 + 10.0im, 0b1010 => 2.0)
      @test isapprox(psi1, psi2; atol=1E-8)
    end

    @testset "HilbertSpaceRepresentation" begin
      hsr = represent(hs)
      psi1 = ExactDiagonalization.Toolkit.product_state(hsr, local_states)
      psi2 = zeros(ComplexF64, 16)
      psi2[9] = 1.0 + 10.0im
      psi2[11] = 2.0
      @test isapprox(psi1, psi2; atol=1E-8)
    end
    @testset "sector-rep" begin
      hssr = represent(HilbertSpaceSector(hs, 0))
      # 0011, 0101, 0110, 1001, 1010, 1100
      psi1 = ExactDiagonalization.Toolkit.product_state(hssr, local_states)
      psi2 = zeros(ComplexF64, 6)
      psi2[5] = 2.0
      @test isapprox(psi1, psi2; atol=1E-8)
    end
  end
end
