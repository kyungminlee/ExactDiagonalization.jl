export simplify

simplify(op ::NullOperator; tol::AbstractFloat=0.0) = op

function simplify(op ::PureOperator{S, BR}; tol ::AbstractFloat=sqrt(eps(Float64))) where {S<:Real, BR}
  if isapprox(op.amplitude, 0.0; rtol=tol, atol=tol)
    return NullOperator()
  else
    return op
  end
end

function simplify(op ::PureOperator{Complex{S}, BR}; tol ::AbstractFloat=sqrt(eps(Float64))) where {S<:Real, BR}
  if isapprox(abs(op.amplitude), 0.0; rtol=tol, atol=tol)
    return NullOperator()
  end
  if isapprox(imag(op.amplitude), 0; rtol=tol, atol=tol)
    return real(op)
  end
  return op
end


function simplify(so ::SumOperator{S, BR}; tol::AbstractFloat=sqrt(eps(Float64))) where {S, BR}

  terms ::Vector{PureOperator{S, BR}} = filter((x) -> !isa(x, NullOperator), simplify.(so.terms))
  
  if isempty(terms)
    return NullOperator()
  end
  
  sort!(terms; lt=(<))
  new_terms = PureOperator{S, BR}[]

  hs ::AbstractHilbertSpace = terms[1].hilbert_space
  bm ::BR = terms[1].bitmask
  bs ::BR = terms[1].bitsource
  bt ::BR = terms[1].bittarget
  am ::S  = terms[1].amplitude

  for term in terms[2:end]
    if (bm == term.bitmask) && (bs == term.bitsource) && (bt == term.bittarget)
      am += term.amplitude
    else
      if ! isapprox(am, 0; rtol=tol, atol=tol)
        push!(new_terms, PureOperator{S, BR}(hs, bm, bs, bt, am))
      end
      bm = term.bitmask
      bs = term.bitsource
      bt = term.bittarget
      am = term.amplitude
    end
  end

  if ! isapprox(am, 0; rtol=tol, atol=tol)
    push!(new_terms, PureOperator{S, BR}(hs, bm, bs, bt, am))
  end

  if isempty(new_terms)
    return NullOperator()
  end

  if S <: Complex && isapprox(maximum(abs(imag(t.amplitude)) for t in new_terms), 0; atol=tol, rtol=tol)
    R = real(S)
    if length(new_terms) == 1
      return real(new_terms[1])
    else
      return SumOperator{R, BR}(so.hilbert_space, real.(new_terms))
    end
  else
    if length(new_terms) == 1
      return new_terms[1]
    else
      return SumOperator{S, BR}(so.hilbert_space, new_terms)
    end
  end
  
end