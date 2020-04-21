"""
    spin_half_system(n_sites)

Create a Hilbert space of spin 1/2 system of `n_sites`

# Returns
- `(hilbert_space, pauli)`
"""
function spin_half_system(n_sites::Integer)
  up = State("Up", 1)
  dn = State("Dn", -1)
  spin_site = Site([up, dn])
  hilbert_space = HilbertSpace([spin_site for i in 1:n_sites])
  pauli(isite::Integer, j::Symbol) = pauli_matrix(hilbert_space, isite, j)
  return (hilbert_space, pauli)
end


"""
    pauli_matrix(hs, isite, j)

Return an Operator for Pauli Matrix at site `isite`.
`j` is one of `:x`, `:y`, `:z`, `:+`, `:-`.
"""
function pauli_matrix(hs::HilbertSpace, isite::Integer, j::Symbol)
  if j == :x
    return pure_operator(hs, isite, 1, 2, 1, UInt) + pure_operator(hs, isite, 2, 1, 1, UInt)
  elseif j == :y
    return pure_operator(hs, isite, 1, 2, -im, UInt) + pure_operator(hs, isite, 2, 1, im, UInt)
  elseif j == :z
    return pure_operator(hs, isite, 1, 1, 1, UInt) + pure_operator(hs, isite, 2, 2, -1, UInt)
  elseif j == :+
    return pure_operator(hs, isite, 1, 2, 1, UInt)
  elseif j == :-
    return pure_operator(hs, isite, 2, 1, 1, UInt)
  else
    throw(ArgumentError("pauli matrix of type $(j) not supported"))
  end
end
