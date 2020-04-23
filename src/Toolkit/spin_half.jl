"""
    spin_half_system(n_sites)

Create a Hilbert space of spin 1/2 system of `n_sites`

# Returns
- `(hilbert_space, pauli)`
"""
function spin_half_system(n_sites::Integer, ::Type{BR}=UInt) where {BR<:Unsigned}
  if n_sites > sizeof(BR) * 8
    throw(ArgumentError("spin half system of $n_sites sites cannot be expressed using $BR"))
  end
  up = State("Up", 1)
  dn = State("Dn", -1)
  spin_site = Site([up, dn])
  hilbert_space = HilbertSpace([spin_site for i in 1:n_sites])
  pauli(isite::Integer, j::Symbol) = pauli_matrix(hilbert_space, isite, j, BR)
  return (hilbert_space, pauli)
end


"""
    pauli_matrix(hs, isite, j)

Return an Operator for Pauli Matrix at site `isite`.
`j` is one of `:x`, `:y`, `:z`, `:+`, `:-`.
"""
function pauli_matrix(hs::HilbertSpace, isite::Integer, j::Symbol, ::Type{BR}=UInt) where {BR<:Unsigned}
  if j == :x
    return pure_operator(hs, isite, 1, 2, 1, BR) + pure_operator(hs, isite, 2, 1, 1, BR)
  elseif j == :y
    return pure_operator(hs, isite, 1, 2, -im, BR) + pure_operator(hs, isite, 2, 1, im, BR)
  elseif j == :z
    return pure_operator(hs, isite, 1, 1, 1, BR) + pure_operator(hs, isite, 2, 2, -1, BR)
  elseif j == :+
    return pure_operator(hs, isite, 1, 2, 1, BR)
  elseif j == :-
    return pure_operator(hs, isite, 2, 1, 1, BR)
  else
    throw(ArgumentError("pauli matrix of type $(j) not supported"))
  end
end
