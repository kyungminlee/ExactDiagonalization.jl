function spin_system(n_sites::Integer, S::Integer, ::Type{BR}=UInt) where {BR<:Unsigned}
  return spin_system(n_sites, S//1, BR)
end

"""
    spin_system(n_sites, S=1//2)

Create a Hilbert space of spin 1/2 system of `n_sites`

# Returns
- `(hilbert_space, pauli)`
"""
function spin_system(n_sites::Integer, S::Rational=1//2, ::Type{BR}=UInt) where {BR<:Unsigned}
  if (S <= 0)
    throw(ArgumentError("S needs to positive"))
  elseif !(denominator(S) == 1 || denominator(S) == 2)
    throw(ArgumentError("S needs to be either integer or half-integer"))
  end
  site = Site([State("Sz=$Sz", Int(2*Sz)) for Sz in -S:1:S])
  hilbert_space = HilbertSpace([site for i in 1:n_sites])
  twoS = Int(S*2)

  if ceil(Int, log2(twoS + 1 - eps(Float64)))*n_sites > sizeof(BR) * 8
    throw(ArgumentError("spin $S system of $n_sites sites cannot be expressed using $BR"))
  end

  Ssq = 0.25 * twoS * (twoS+2)
  Sr = float(S)
  function spin(isite::Integer, j::Symbol)
    A(k::Integer) = sqrt(Ssq - (Sr-k)*(Sr-k+1)) # offdiagonal amplitude
    # m = S - k
    if j == :+
      return sum(pure_operator(hilbert_space, isite, k, k+1, A(k), BR) for k in 1:twoS)
    elseif j == :-
      return sum(pure_operator(hilbert_space, isite, k+1, k, A(k), BR) for k in 1:twoS)
    elseif j == :x
      Sp = sum(pure_operator(hilbert_space, isite, k, k+1, A(k), BR) for k in 1:twoS)
      Sm = sum(pure_operator(hilbert_space, isite, k+1, k, A(k), BR) for k in 1:twoS)
      return (Sp + Sm) * 0.5
    elseif j == :y
      Sp = sum(pure_operator(hilbert_space, isite, k, k+1, A(k), BR) for k in 1:twoS)
      Sm = sum(pure_operator(hilbert_space, isite, k+1, k, A(k), BR) for k in 1:twoS)
      return (Sp - Sm) * (-0.5im)
    elseif j == :z
      return sum(pure_operator(hilbert_space, isite, k+1, k+1, Sr - k, BR) for k in 0:twoS)
    else
      throw(ArgumentError("spin component of type $(j) not supported"))
    end
  end
  return (hilbert_space, spin)
end
