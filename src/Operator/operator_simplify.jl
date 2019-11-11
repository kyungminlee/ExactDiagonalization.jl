export simplify

simplify(op ::NullOperator; tol::Real=0.0) = op

function simplify(op ::PureOperator{S, BR}; tol ::Real=Base.rtoldefault(S)) where {S<:Real, BR}
  if isapprox(op.amplitude, zero(S); atol=tol)
    return NullOperator()
  else
    return op
  end
end

function simplify(op ::PureOperator{Complex{S}, BR}; tol ::Real=Base.rtoldefault(S)) where {S<:Real, BR}
  if isapprox(abs(op.amplitude), zero(S); atol=tol)
    return NullOperator()
  end
  if isapprox(imag(op.amplitude), zero(S); atol=tol)
    return real(op)
  end
  return op
end


function simplify(so ::SumOperator{S, BR}; tol ::Real=Base.rtoldefault(real(S))) where {S, BR}
  terms ::Vector{PureOperator{S, BR}} = filter((x) -> !isa(x, NullOperator), simplify.(so.terms))

  isempty(terms) && return NullOperator()

  sort!(terms; lt=(<))
  new_terms = PureOperator{S, BR}[]

  bm ::BR = terms[1].bitmask
  br ::BR = terms[1].bitrow
  bc ::BR = terms[1].bitcol
  am ::S  = terms[1].amplitude

  for term in terms[2:end]
    if (bm == term.bitmask) && (br == term.bitrow) && (bc == term.bitcol)
      am += term.amplitude
    else
      if ! isapprox(am, zero(S); atol=tol)
        push!(new_terms, PureOperator{S, BR}(bm, br, bc, am))
      end
      bm = term.bitmask
      br = term.bitrow
      bc = term.bitcol
      am = term.amplitude
    end
  end

  if ! isapprox(am, zero(S); atol=tol)
    push!(new_terms, PureOperator{S, BR}(bm, br, bc, am))
  end

  if isempty(new_terms)
    return NullOperator()
  end

  R = real(S)
  if S <: Complex && isapprox(maximum(abs(imag(t.amplitude)) for t in new_terms), zero(R); atol=tol)
    if length(new_terms) == 1
      return real(new_terms[1])
    else
      return SumOperator{R, BR}(real.(new_terms))
    end
  else
    if length(new_terms) == 1
      return new_terms[1]
    else
      return SumOperator{S, BR}(new_terms)
    end
  end

end
