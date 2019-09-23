using Test
using ExactDiagonalization

@testset "operatorapplication" begin
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
    hs = AbstractHilbertSpace([spin_site, spin_site, spin_site, spin_site])
    σ(i::Integer, j::Symbol) = pauli_matrix(hs, i, j)

    chs = concretize(hs, 0)
    psi = SparseState{Float64, UInt}(hs, UInt(0b0011) => 2.0, UInt(0b0101) => 10.0)
    
    @testset "hilbert" begin
      hs2 = AbstractHilbertSpace([spin_site, spin_site, spin_site,])
      psi2 = SparseState{Float64, UInt}(hs2, UInt(0b0011) => 2.0, UInt(0b0101) => 10.0)
      @test_throws ArgumentError apply(σ(1, :+), psi2)
      @test_throws ArgumentError apply(σ(1, :x), psi2)
    end

    @testset "apply" begin
      @test apply(NullOperator(), psi) == SparseState{Float64, UInt}(hs)
    
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

    @testset "apply!" begin
      let
        psi2 = SparseState{Float64, UInt}(hs, UInt(0b0010) => 0.25)
        apply!(psi2, NullOperator(), psi)
        @test psi2 == SparseState{Float64, UInt}(hs, UInt(0b0010) => 0.25)
      end

      let
        psi2 = SparseState{Float64, UInt}(hs2, UInt(0b0011) => 2.0, UInt(0b0101) => 10.0)
        @test_throws ArgumentError apply!(psi2, σ(1, :+), psi)
        @test_throws ArgumentError apply!(psi2, σ(1, :x), psi)
      end

      let
        psi2 = SparseState{Float64, UInt}(hs, UInt(0b0010) => 0.25)
        apply!(psi2, σ(1, :+), psi)
        @test psi2 == SparseState{Float64, UInt}(hs, UInt(0b0010) => 0.25)
      end

      let
        psi2 = SparseState{Float64, UInt}(hs, UInt(0b0010) => 0.25)
        apply!(psi2, σ(1, :-), psi)
        @test psi2 == SparseState{Float64, UInt}(hs, UInt(0b0010) => 2.25, UInt(0b0100) => 10.0)
      end

      let
        psi2 = SparseState{Float64, UInt}(hs, UInt(0b0010) => 0.25)
        apply!(psi2, σ(2, :+), psi)
        @test psi2 == SparseState{Float64, UInt}(hs, UInt(0b0010) => 0.25, UInt(0b0111) => 10.0)
      end

      let
        psi2 = SparseState{Float64, UInt}(hs, UInt(0b0010) => 0.25)
        apply!(psi2, σ(2, :-), psi)
        @test psi2 == SparseState{Float64, UInt}(hs, UInt(0b0010) => 0.25, UInt(0b0001) => 2.0)
      end

      for i1 in 1:4, j1 in [:x, :y, :z, :+, :-]
        ϕ1 = apply(σ(i1, j1), psi)
        ϕ2 = SparseState{ComplexF64, UInt}(hs)
        apply!(ϕ2, σ(i1, j1), psi)
        @test ϕ1 == ϕ2
        ϕ3 = SparseState{Float64, UInt}(hs)
      end

      let
        psi2 = SparseState{ComplexF64, UInt}(hs, UInt(0b0010) => 0.25)
        psi3 = apply(σ(2, :y), psi2)
        psi4 = SparseState{Float64, UInt}(hs)
        @test_throws InexactError apply!(psi4, σ(2, :y), psi2)
        # TODO more
      end

    end


    @testset "materialize" begin
      n = 4
      heisenberg = sum(2 * σ(i, :+) * σ( mod(i, n) + 1, :-) for i in 1:n)
      heisenberg += sum(2 * σ(i, :-) * σ( mod(i, n) + 1, :+) for i in 1:n)
      heisenberg += sum(σ(i, :z) * σ( mod(i, n) + 1, :z) for i in 1:n)
      field_x = sum(σ(i, :x) for i in 1:n)
      field_y = sum(σ(i, :y) for i in 1:n)
      field_z = sum(σ(i, :z) for i in 1:n)
      @testset "zerosector" begin
        chs = concretize(hs, 0)
        H, ϵ = materialize(chs, heisenberg)
        @test H ≈ [ 0  2  0  0  2  0;
                    2 -4  2  2  0  2;
                    0  2  0  0  2  0;
                    0  2  0  0  2  0;
                    2  0  2  2 -4  2;
                    0  2  0  0  2  0]
        Hp, ϵp = materialize_parallel(chs, heisenberg)
        @test Hp ≈ H

        heisenberg2 = sum(σ(i, j) * σ(mod(i,n)+1, j) for i in 1:n, j in [:x, :y, :z])
        H2, ϵ2 = materialize(chs, heisenberg2)
        @test H2 ≈ H
        H3, ϵ3 = materialize_parallel(chs, heisenberg2)
        @test H3 ≈ H

        let
          Bx, ϵx = materialize(chs, field_x)
          By, ϵy = materialize(chs, field_y)
          Bz, ϵz = materialize(chs, field_z)
          Bxp, ϵxp = materialize_parallel(chs, field_x)
          Byp, ϵyp = materialize_parallel(chs, field_y)
          Bzp, ϵzp = materialize_parallel(chs, field_z)
          @test ϵx > 0
          @test ϵy > 0
          @test ϵz ≈ 0
          @test ϵx ≈ ϵxp
          @test ϵy ≈ ϵyp
          @test ϵz ≈ ϵzp
        end
      end

      @testset "allsector" begin
        PAULI_MATRIX = [
          [0.0 1.0; 1.0 0.0], [0.0 -1.0im; 1.0im 0.0], [1.0 0.0; 0.0 -1.0], [1.0 0.0; 0.0 1.0]
        ]
        chs = concretize(hs)
        H, ϵ = materialize(chs, heisenberg)
        H2  = sum( kron(PAULI_MATRIX[i], PAULI_MATRIX[i], PAULI_MATRIX[4], PAULI_MATRIX[4]) for i in 1:3)
        H2 += sum( kron(PAULI_MATRIX[4], PAULI_MATRIX[i], PAULI_MATRIX[i], PAULI_MATRIX[4]) for i in 1:3)
        H2 += sum( kron(PAULI_MATRIX[4], PAULI_MATRIX[4], PAULI_MATRIX[i], PAULI_MATRIX[i]) for i in 1:3)
        H2 += sum( kron(PAULI_MATRIX[i], PAULI_MATRIX[4], PAULI_MATRIX[4], PAULI_MATRIX[i]) for i in 1:3)
        @test H ≈ H2

        let
          Bx, _ = materialize(chs, field_x)
          By, _ = materialize(chs, field_y)
          Bz, _ = materialize(chs, field_z)
          @test isa(Bx.nzval, Array{Float64, 1})
          @test isa(By.nzval, Array{ComplexF64, 1})
          @test isa(Bz.nzval, Array{Float64, 1})
        end
        let
          Bx, _ = materialize_parallel(chs, field_x)
          By, _ = materialize_parallel(chs, field_y)
          Bz, _ = materialize_parallel(chs, field_z)
          @test isa(Bx.nzval, Array{Float64, 1})
          @test isa(By.nzval, Array{ComplexF64, 1})
          @test isa(Bz.nzval, Array{Float64, 1})
        end
      end
    end
  end
end
