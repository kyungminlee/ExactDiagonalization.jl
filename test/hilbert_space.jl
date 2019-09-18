using Test
using ExactDiagonalization

using StaticArrays

@testset "HilbertSpace" begin

  @testset "Spin" begin
    QN = Int
    up = State("Up", 1)
    dn = State{QN}("Dn",-1)
  end

  @testset "SpinCharge" begin
    QN = SVector{2, Int}
    em = State("Em", QN( 0, 0))
    up = State("Up", QN( 1, 1))
    dn = State("Dn", QN(-1, 1))
    ud = State("UpDn", QN( 0, 2))
    @test_throws MethodError State{QN}("X", 1)

    site = Site([em, up, dn])
    @test_throws MethodError Site([em, up, dn, State("X", 1)])
    
    #@show State("X", 1)
    #@show State("UpDn", QN( 0, 2))
  end


end