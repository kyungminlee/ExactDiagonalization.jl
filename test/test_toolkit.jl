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

  @testset "product_state" begin
    (hs, pauli) = ExactDiagonalization.Toolkit.spin_half_system(4)
    local_states = [[1.0, 0.0], [1.0 + 10.0im, 2.0], [1.0, 0.0], [0.0, 1.0]]
    psi1 = ExactDiagonalization.Toolkit.product_state(hs, local_states)
    psi2 = SparseState{ComplexF64, UInt}(0b1000 => 1.0 + 10.0im, 0b1010 => 2.0)
    @test isapprox(psi1, psi2; atol=1E-8)
  end
end
