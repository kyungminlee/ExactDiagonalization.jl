using Test
using ExactDiagonalization

using LinearAlgebra
using TightBindingLattice
using ExactDiagonalization.Toolkit: pauli_matrix

@testset "symmetry_reduce" begin
  QN = Tuple{Int}
  up = State("Up", 1);
  dn = State("Dn",-1);
  spin_site = Site{QN}([up, dn]);

  n_sites = 4;
  hs = HilbertSpace([spin_site for i in 1:n_sites]);
  hss = HilbertSpaceSector(hs, 0)
  hsr = represent(hss);

  translation_group = TranslationGroup([Permutation([2,3,4,1])])

  for symred in [symmetry_reduce, symmetry_reduce_serial, symmetry_reduce_parallel]
    rhsr = symred(hsr, translation_group, [0//1])
    @test rhsr.basis_list == UInt[0b0011, 0b0101]
    @test rhsr.parent === hsr

    rhsr = symred(hsr, translation_group, [1//4])
    @test rhsr.basis_list == UInt[0b0011]
    @test rhsr.parent === hsr

    rhsr = symred(hsr, translation_group, [2//4])
    @test rhsr.basis_list == UInt[0b0011, 0b0101]
    @test rhsr.parent === hsr

    rhsr = symred(hsr, translation_group, [3//4])
    @test rhsr.basis_list == UInt[0b0011]
    @test rhsr.parent === hsr

    @test_throws ArgumentError symred(hsr, translation_group, [1//5])

    tol = sqrt(eps(Float64))
    for k in translation_group.fractional_momenta
      rhsr = symred(hsr, translation_group, k)
      for (i_p, b) in enumerate(hsr.basis_list)
        if b in rhsr.basis_list
          @test 1 <= rhsr.basis_mapping_index[i_p] <= dimension(rhsr)
          @test rhsr.basis_list[rhsr.basis_mapping_index[i_p]] == b
          @test isapprox(imag(rhsr.basis_mapping_amplitude[i_p]), 0; atol=tol)
        end
      end
      for (i_p, i_r) in enumerate(rhsr.basis_mapping_index)
        if i_r == -1
          @test ! (rhsr.parent.basis_list[i_p] in rhsr.basis_list)
        end
      end

      @test_throws DimensionMismatch symmetry_unreduce(rhsr, zeros(ComplexF64, dimension(rhsr)+1))
      @test_throws DimensionMismatch symmetry_reduce(rhsr, zeros(ComplexF64, dimension(hsr)+1))

      sv = rand(ComplexF64, dimension(rhsr))
      lv = symmetry_unreduce(rhsr, sv)
      sv2 = symmetry_reduce(rhsr, lv)
      @test isapprox(sv, sv2; atol=tol)
    end
  end

  @testset "convention" begin
    # want  |ψ(k)⟩ = ∑ exp(+ikx) |ψ(x)⟩
    QN = Tuple{Int}
    up = State("Up", 1);
    dn = State("Dn",-1);
    spin_site = Site{QN}([up, dn]);
    n_sites = 7;
    hs = HilbertSpace([spin_site for i in 1:n_sites]);

    let
      hss = HilbertSpaceSector(hs, 5)
      hsr = represent(hss)
      p = Permutation([2,3,4,5,6,7,1])
      translation_group = TranslationGroup(p)
      @test symmetry_apply(hs, p, 0b0000001) == 0b0000010
      rhsr = symmetry_reduce(hsr, translation_group, [1//7])
      @test hsr.basis_list == UInt[0b0000001, 0b0000010, 0b0000100, 0b0001000, 0b0010000, 0b0100000, 0b1000000]
      @test rhsr.basis_list == UInt[0b0000001]
      ψk = symmetry_unreduce(rhsr, [1.0])
      @test isapprox(ψk, [cis(2π * i/n_sites) / sqrt(n_sites) for i in 0:(n_sites-1)]; atol=sqrt(eps(Float64)))
    end

    let
      hss = HilbertSpaceSector(hs, -5)
      hsr = represent(hss)
      p = Permutation([2,3,4,5,6,7,1])
      translation_group = TranslationGroup(p)
      @test symmetry_apply(hs, p, 0b0000001) == 0b0000010
      rhsr = symmetry_reduce(hsr, translation_group, [1//7])
      # opposite order
      @test hsr.basis_list == UInt[0b0111111, 0b1011111, 0b1101111, 0b1110111, 0b1111011, 0b1111101, 0b1111110]
      @test rhsr.basis_list == UInt[0b0111111]
      ψk = symmetry_unreduce(rhsr, [1.0])
      @test isapprox(ψk, [cis(-2π * i/n_sites) / sqrt(n_sites) for i in 0:(n_sites-1)]; atol=sqrt(eps(Float64)))
    end

  end

  # σ = Dict( (isite, j) => pauli_matrix(hs, isite, j) for isite in 1:n_sites, j in [:x, :y, :z, :+, :-]);
  # j1 = simplify( sum(σ[(i, j)] * σ[( mod(i, n_sites) + 1 , j)] for i in 1:n_sites, j in [:x, :y, :z]) );
end
