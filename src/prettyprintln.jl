export prettyprintln

prettyprintln(x...) = prettyprintln(stdout, x...)

function prettyprintln(io::IO, op::NullOperator)
  println(io, "NullOperator")
end

function prettyprintln(io::IO, op::PureOperator{S, BR}; prefix::AbstractString="") where {S, BR}
  println(io, prefix, "PureOperator")
  println(io, prefix, "| M: ", string(op.bitmask, base=2, pad=op.hilbert_space.bitoffsets[end]))
  println(io, prefix, "| S: ", string(op.bitsource, base=2, pad=op.hilbert_space.bitoffsets[end]))
  println(io, prefix, "| T: ", string(op.bittarget, base=2, pad=op.hilbert_space.bitoffsets[end]))
  println(io, prefix, "| A: ", op.amplitude) 
end

function prettyprintln(io::IO, op::SumOperator{S, BR}; prefix::AbstractString="") where {S, BR}
  println(io, prefix, "SumOperator")
  for t in op.terms
    prettyprintln(io, t; prefix=string(prefix, "| "))
  end
end

function prettyprintln(io ::IO, psi::SparseState{S, BR}; prefix::AbstractString="") where {S, BR}
  println(io, prefix, "SparseState")
  bs = sort(collect(keys(psi.components)))
  for b in bs
    println(io, prefix, "| ", string(b, base=2, pad=psi.hilbert_space.bitoffsets[end]), " : ", psi.components[b])
  end
end