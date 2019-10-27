using Test
using ExactDiagonalization

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

@testset "symmetry_reduce" begin
  QN = Int;
  up = State{QN}("Up", 1);
  dn = State{QN}("Dn",-1);
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
          @test 1 <= rhsr.basis_mapping[i_p].index <= dimension(rhsr)
          @test rhsr.basis_list[rhsr.basis_mapping[i_p].index] == b
          @test isapprox(imag(rhsr.basis_mapping[i_p].amplitude), 0; atol=tol)
        end
      end
      for (i_p, (i_r, amplitude)) in enumerate(rhsr.basis_mapping)
        if i_r == -1
          @test ! (rhsr.parent.basis_list[i_p] in rhsr.basis_list)
        end
      end

      sv = rand(ComplexF64, dimension(rhsr))
      lv = symmetry_unreduce(rhsr, sv)
      sv2 = symmetry_reduce(rhsr, lv)
      @test isapprox(sv, sv2; atol=tol)
    end
  end
  # σ = Dict( (isite, j) => pauli_matrix(hs, isite, j) for isite in 1:n_sites, j in [:x, :y, :z, :+, :-]);
  # j1 = simplify( sum(σ[(i, j)] * σ[( mod(i, n_sites) + 1 , j)] for i in 1:n_sites, j in [:x, :y, :z]) );
end
