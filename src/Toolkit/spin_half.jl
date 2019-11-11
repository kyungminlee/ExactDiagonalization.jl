using ExactDiagonalization

"""
    spin_half_system(n_sites)

Create a Hilbert space of spin 1/2 system of `n_sites`
"""
function spin_half_system(n_sites ::Integer)
  QN = Int
  up = State{QN}("Up", 1)
  dn = State{QN}("Dn", -1)
  spin_site = Site{QN}([up, dn])
  hilbert_space = HilbertSpace([spin_site for i in 1:n_sites])
  pauli(isite ::Integer, j::Symbol) = pauli_matrix(hilbert_space, isite, j)
  return (hilbert_space, pauli)
end

"""
    pauli_matrix(hs::HilbertSpace, isite ::Integer, j ::Symbol)

Return an Operator for Pauli Matrix at site `isite`.

# Arguments
- `hs ::HilbertSpace`
- `isite ::Integer` : site index
- `j ::Symbol`: one of `:x`, `:y`, `:z`, `:+`, `:-`
"""
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
