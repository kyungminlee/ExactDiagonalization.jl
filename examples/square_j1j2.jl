using LinearAlgebra
using SparseArrays
using Arpack
using JLD2
using ExactDiagonalization
using ArgParse
using Memento

main_logger = Memento.config!("info"; fmt="[{date} | {level} | {name}]: {msg}")

function make_square_lattice(n1 ::Integer, n2 ::Integer)
  n_sites = n1 * n2
  
  nearest_neighbor_bond_types = [(1,0), (0,1)]
  nearest_neighbor_pairs = Tuple{Int, Int}[]
  for i1 in 1:n1
      for i2 in 1:n2
          i = (i1-1) * n2 + (i2-1) + 1
          for (d1, d2) in nearest_neighbor_bond_types
              j1 = mod(i1 + d1 - 1, n1) + 1
              j2 = mod(i2 + d2 - 1, n2) + 1
              j = (j1-1) * n2 + (j2-1) + 1

              push!(nearest_neighbor_pairs, (i,j))
          end
      end
  end

  second_nearest_neighbor_bond_types = [(1,1), (-1,1)]
  second_nearest_neighbor_pairs = Tuple{Int, Int}[]
  for i1 in 1:n1
      for i2 in 1:n2
          i = (i1-1) * n2 + (i2-1) + 1
          for (d1, d2) in second_nearest_neighbor_bond_types
              j1 = mod(i1 + d1 - 1, n1) + 1
              j2 = mod(i2 + d2 - 1, n2) + 1
              j = (j1-1) * n2 + (j2-1) + 1

              push!(second_nearest_neighbor_pairs, (i,j))
          end
      end
  end

  #=
        o
      / | \
    o - o - o 
      \ | /
        o
  =#
  chiral_triplet_types = [
      (( 1, 0), ( 0, 1)),
      (( 0, 1), (-1, 0)),
      ((-1, 0), ( 0,-1)),
      (( 0,-1), ( 1, 0)),
  ]
  chiral_triplets = Tuple{Int64, Int64, Int64}[]
  for i1 in 1:n1, i2 in 1:n2
    i = (i1-1) * n2 + (i2-1) + 1
    for ((d1, d2), (e1, e2)) in chiral_triplet_types
      j1 = mod(i1 + d1 - 1, n1) + 1
      j2 = mod(i2 + d2 - 1, n2) + 1
      j = (j1-1) * n2 + (j2-1) + 1

      k1 = mod(i1 + e1 - 1, n1) + 1
      k2 = mod(i2 + e2 - 1, n2) + 1
      k = (k1-1) * n2 + (k2-1) + 1
      push!(chiral_triplets, (i,j,k))
    end
  end
  return (nearest_neighbor_pairs, second_nearest_neighbor_pairs, chiral_triplets)
end

function make_J1J2J3_hamiltonian(n1::Integer, n2::Integer)
  QN = Int
  up = State{QN}("Up", 1)
  dn = State{QN}("Dn",-1)
  spin_site = Site{QN}([up, dn])

  n_sites = n1 * n2
  info(main_logger, "Building AbstractHilbertSpace")
  hs = AbstractHilbertSpace([spin_site for i in 1:n_sites])


  PAULI_MATRICES = [ Float64[0 1.0; 1.0 0.0],
                     ComplexF64[0.0 -1.0*im; 1.0*im 0.0],
                     Float64[1.0 0.0; 0.0 -1.0] ]

  function σ(i ::Integer, j::Symbol) ::KroneckerProductOperator{ComplexF64}
    if j == :x
      return KroneckerProductOperator{ComplexF64}(hs, 1.0, Dict(i=>PAULI_MATRICES[1]))
    elseif j == :y
      return KroneckerProductOperator{ComplexF64}(hs, 1.0, Dict(i=>PAULI_MATRICES[2]))
    elseif j == :z
      return KroneckerProductOperator{ComplexF64}(hs, 1.0, Dict(i=>PAULI_MATRICES[3]))
    elseif j == :+
      return KroneckerProductOperator{ComplexF64}(hs, 1.0, Dict(i=>[0.0 1.0; 0.0 0.0]))
    elseif j == :-
      return KroneckerProductOperator{ComplexF64}(hs, 1.0, Dict(i=>[0.0 0.0; 1.0 0.0]))
    else
      throw(ArgumentError("pauli matrix of type $(j) not supported"))
    end
  end

  (nearest_neighbor_pairs, second_nearest_neighbor_pairs, chiral_triplets) = make_square_lattice(n1, n2)

  j1_terms = KroneckerProductOperator{ComplexF64}[]
  for (i,j) in nearest_neighbor_pairs
    push!(j1_terms, 2 * σ(i, :+) * σ(j, :-))
    push!(j1_terms, 2 * σ(i, :-) * σ(j, :+))
    push!(j1_terms, σ(i, :z) * σ(j, :z))
  end

  j2_terms = KroneckerProductOperator{ComplexF64}[]
  for (i,j) in second_nearest_neighbor_pairs
    push!(j2_terms, 2 * σ(i, :+) * σ(j, :-))
    push!(j2_terms, 2 * σ(i, :-) * σ(j, :+))
    push!(j2_terms, σ(i, :z) * σ(j, :z))
  end

  j3_terms = KroneckerProductOperator{ComplexF64}[]
  for (i,j,k) in chiral_triplets
    push!(j3_terms, 2*im * σ(i, :z) * σ(j, :+) * σ(k, :-))
    push!(j3_terms,-2*im * σ(i, :z) * σ(j, :-) * σ(k, :+))

    push!(j3_terms, 2*im * σ(j, :z) * σ(k, :+) * σ(i, :-))
    push!(j3_terms,-2*im * σ(j, :z) * σ(k, :-) * σ(i, :+))
    
    push!(j3_terms, 2*im * σ(k, :z) * σ(i, :+) * σ(j, :-))
    push!(j3_terms,-2*im * σ(k, :z) * σ(i, :-) * σ(j, :+))
  end
  return (hs, j1_terms, j2_terms, j3_terms)
end

function make_J1J2J3_hamiltonian2(n1::Integer, n2::Integer)

  QN = Int
  up = State{QN}("Up", 1)
  dn = State{QN}("Dn",-1)
  spin_site = Site{QN}([up, dn])

  n_sites = n1 * n2
  info(main_logger, "Building AbstractHilbertSpace")
  hs = AbstractHilbertSpace([spin_site for i in 1:n_sites])
  
  PAULI_MATRICES = [ Float64[0 1.0; 1.0 0.0],
                     ComplexF64[0.0 -1.0*im; 1.0*im 0.0],
                     Float64[1.0 0.0; 0.0 -1.0] ]

  function σraw(i ::Integer, j::Symbol) ::KroneckerProductOperator{ComplexF64}
    if j == :x
      return KroneckerProductOperator{ComplexF64}(hs, 1.0, Dict(i=>PAULI_MATRICES[1]))
    elseif j == :y
      return KroneckerProductOperator{ComplexF64}(hs, 1.0, Dict(i=>PAULI_MATRICES[2]))
    elseif j == :z
      return KroneckerProductOperator{ComplexF64}(hs, 1.0, Dict(i=>PAULI_MATRICES[3]))
    elseif j == :+
      return KroneckerProductOperator{ComplexF64}(hs, 1.0, Dict(i=>[0.0 1.0; 0.0 0.0]))
    elseif j == :-
      return KroneckerProductOperator{ComplexF64}(hs, 1.0, Dict(i=>[0.0 0.0; 1.0 0.0]))
    else
      throw(ArgumentError("pauli matrix of type $(j) not supported"))
    end
  end

  function σ(i ::Integer, j ::Symbol)
    return convert(SumOperator{ComplexF64, UInt}, σraw(i,j))
  end
  (nearest_neighbor_pairs, second_nearest_neighbor_pairs, chiral_triplets) = make_square_lattice(n1, n2)

  j1 = NullOperator()
  j1_terms = []
  for (i,j) in nearest_neighbor_pairs
    # j1 += 2 * σ(i, :+) * σ(j, :-)
    # j1 += 2 * σ(i, :-) * σ(j, :+)
    # j1 += σ(i, :z) * σ(j, :z)
    push!(j1_terms, 2 * σ(i, :+) * σ(j, :-))
    push!(j1_terms, 2 * σ(i, :-) * σ(j, :+))
    push!(j1_terms, σ(i, :z) * σ(j, :z))
  end

  j2 = NullOperator()
  j2_terms = []
  for (i,j) in second_nearest_neighbor_pairs
    #j2 += 2 * σ(i, :+) * σ(j, :-)
    #j2 += 2 * σ(i, :-) * σ(j, :+)
    #j2 += σ(i, :z) * σ(j, :z)
    push!(j2_terms, 2 * σ(i, :+) * σ(j, :-))
    push!(j2_terms, 2 * σ(i, :-) * σ(j, :+))
    push!(j2_terms, σ(i, :z) * σ(j, :z))
  end

  j3 = NullOperator()
  j3_terms  =[]
  for (i,j,k) in chiral_triplets
    # j3 += 2*im * σ(i, :z) * σ(j, :+) * σ(k, :-)
    # j3 += -2*im * σ(i, :z) * σ(j, :-) * σ(k, :+)

    # j3 += 2*im * σ(j, :z) * σ(k, :+) * σ(i, :-)
    # j3 += -2*im * σ(j, :z) * σ(k, :-) * σ(i, :+)
    
    # j3 += 2*im * σ(k, :z) * σ(i, :+) * σ(j, :-)
    # j3 += -2*im * σ(k, :z) * σ(i, :-) * σ(j, :+)

    push!(j3_terms, 2*im * σ(i, :z) * σ(j, :+) * σ(k, :-))
    push!(j3_terms,-2*im * σ(i, :z) * σ(j, :-) * σ(k, :+))

    push!(j3_terms, 2*im * σ(j, :z) * σ(k, :+) * σ(i, :-))
    push!(j3_terms,-2*im * σ(j, :z) * σ(k, :-) * σ(i, :+))
    
    push!(j3_terms, 2*im * σ(k, :z) * σ(i, :+) * σ(j, :-))
    push!(j3_terms,-2*im * σ(k, :z) * σ(i, :-) * σ(j, :+))
  end
  j1 = sum(j1_terms)
  j2 = sum(j2_terms)
  j3 = sum(j3_terms)

  return (hs, j1, j2, j3)
end


function parse_commandline()
  s = ArgParseSettings()
  @add_arg_table s begin
    "--n1"
      arg_type = Int
      required = true
    "--n2"
      arg_type = Int
      required = true
    "--sz"
      arg_type = Int
      required = true
  end
  return parse_args(s)
end


function main()
  args = parse_commandline()

  #(n1, n2) = (4, 4)
  n1 = args["n1"]
  n2 = args["n2"]
  qn = args["sz"]

  info(main_logger, "Parameteters : n1 = $n1, n2 = $n2, sz = $qn")
  
  info(main_logger, "Making J1, J2, J3 terms")
  #(hs, j1_terms, j2_terms, j3_terms) = make_J1J2J3_hamiltonian(n1, n2)
  (hs, j1, j2, j3) = make_J1J2J3_hamiltonian2(n1, n2)





  #sectors = quantum_number_sectors(hs)
  #sectors = [x for x in sectors if x >= 0]
  #for qn in sectors
  info(main_logger, "Sector = $qn")
  info(main_logger, "Concretizing Hilbert Space"); flush(stdout)
  chs = concretize(hs, Set([qn]))
  info(main_logger, "Hilbert space dimension = $(length(chs.basis_list))")

  info(main_logger, "Materializing J1"); flush(stdout)
  j1_sparse, ϵ = materialize_parallel(chs, j1_terms)
  info(main_logger, "Number of nonzero values = $(length(j1_sparse.nzval))")
  info(main_logger, "Bound for error = $ϵ")

  info(main_logger, "Materializing J2"); flush(stdout)
  j2_sparse, ϵ = materialize_parallel(chs, j2_terms)
  info(main_logger, "Number of nonzero values = $(length(j2_sparse.nzval))")
  info(main_logger, "Bound for error = $ϵ")

  info(main_logger, "Materializing J3"); flush(stdout)
  j3_sparse, ϵ = materialize_parallel(chs, j3_terms)
  info(main_logger, "Number of nonzero values = $(length(j3_sparse.nzval))")
  info(main_logger, "Bound for error = $ϵ")

  J1 = 1.0
  J2s = 0:0.2:2
  J3s = 0:0.2:2
  spectrum = Dict()
  for J2 in J2s, J3 in J3s
    info(main_logger, "Starting (J1, J2, J3) = ($J1, $J2, $J3)")
    info(main_logger, "Constructing the total Hamiltonian")
    H = J1 * j1_sparse + J2 * j2_sparse + J3 * j3_sparse

    if size(H)[1] <= 20
      info(main_logger, "Diagonalizing dense Hamiltonian (size = $(size(H)[1]))")
      spectrum[(J1, J2, J3)] = sort(eigvals(Hermitian(Matrix(H))))
    else
      info(main_logger, "Diagonalizing sparse Hamiltonian (size = $(size(H)[1]))")
      (eigenvalues, eigenvectors) = eigs(H; which=:SR)
      spectrum[(J1, J2, J3)] = sort(real.(eigenvalues))
    end
  end # for J2, J3
  filename = "spectrum_$(n1)_$(n2)_$(qn).jld2"
  info(main_logger, "Saving $filename")
  @save filename n1, n2, qn, spectrum
  #end # for qn
end

if PROGRAM_FILE == @__FILE__ 
  main()
end