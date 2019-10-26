using Test
using Suppressor
using ExactDiagonalization

@testset "prettyprintln" begin
  function pauli_matrix(hs::HilbertSpace, isite ::Integer, j ::Symbol)
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
    hs = HilbertSpace([spin_site, spin_site, spin_site, spin_site])
    σ(i::Integer, j::Symbol) = pauli_matrix(hs, i, j)

    buf = IOBuffer()

    val = NullOperator()
    result1 = @capture_out prettyprintln(val)
    prettyprintln(buf, val)
    result2 = String(take!(buf))
    @test result1 == result2
    @test result1 == "NullOperator\n"

    val = σ(1, :+)
    result1 = @capture_out prettyprintln(val)
    prettyprintln(buf, val)
    result2 = String(take!(buf))
    @test result1 == result2
    @test result1 == join(["PureOperator",
                           "| M: 0000000000000000000000000000000000000000000000000000000000000001",
                           "| R: 0000000000000000000000000000000000000000000000000000000000000000",
                           "| C: 0000000000000000000000000000000000000000000000000000000000000001",
                           "| A: 1.0", ""], "\n")
    val = σ(2, :x)
    result1 = @capture_out prettyprintln(val)
    prettyprintln(buf, val)
    result2 = String(take!(buf))
    @test result1 == result2
    @test result1 == join(["SumOperator",
                           "| PureOperator",
                           "| | M: 0000000000000000000000000000000000000000000000000000000000000010",
                           "| | R: 0000000000000000000000000000000000000000000000000000000000000000",
                           "| | C: 0000000000000000000000000000000000000000000000000000000000000010",
                           "| | A: 1.0",
                           "| PureOperator",
                           "| | M: 0000000000000000000000000000000000000000000000000000000000000010",
                           "| | R: 0000000000000000000000000000000000000000000000000000000000000010",
                           "| | C: 0000000000000000000000000000000000000000000000000000000000000000",
                           "| | A: 1.0",
                           ""], "\n")

    # val = SparseState{Float64, UInt}(hs, UInt(0b0010)=>0.2, UInt(0b100) => 0.3 )
    # result1 = @capture_out prettyprintln(val)
    # prettyprintln(buf, val)
    # result2 = String(take!(buf))
    # @test result1 == result2
    # @test result1 == join(["SparseState",
    #                        "| 0010 : 0.2",
    #                        "| 0100 : 0.3", ""], "\n")
  end
end
