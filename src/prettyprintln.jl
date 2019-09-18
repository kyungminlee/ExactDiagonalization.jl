export prettyprintln

function prettyprintln(op::NullOperator)
  println("NullOperator")
end

function prettyprintln(op::PureOperator{S, BR}; prefix::AbstractString="") where {S, BR}
  println(prefix, "PureOperator")
  println(prefix, "| M: ", string(op.bitmask, base=2, pad=op.hilbert_space.bitoffsets[end]))
  println(prefix, "| S: ", string(op.bitsource, base=2, pad=op.hilbert_space.bitoffsets[end]))
  println(prefix, "| T: ", string(op.bittarget, base=2, pad=op.hilbert_space.bitoffsets[end]))
  println(prefix, "| A: ", op.amplitude) 
end

function prettyprintln(op::SumOperator{S, BR}; prefix::AbstractString="") where {S, BR}
  println(prefix, "SumOperator")
  for t in op.terms
    prettyprintln(t; prefix=string(prefix, "| "))
  end
end

function prettyprintln(psi::SparseState{S, BR}; prefix::AbstractString="") where {S, BR}
  println(prefix, "SparseState")
  bs = sort(collect(keys(psi.components)))
  for b in bs
    println(prefix, "  ", string(b, base=2, pad=psi.hilbert_space.bitoffsets[end]), " : ", psi.components[b])
  end
end