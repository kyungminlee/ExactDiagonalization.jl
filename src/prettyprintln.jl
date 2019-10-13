export prettyprintln

prettyprintln(x; kwargs...) = prettyprintln(stdout, x; kwargs...)

function prettyprintln(io::IO, op::NullOperator; prefix::AbstractString="")
  println(io, prefix, "NullOperator")
end

function prettyprintln(io::IO, op::PureOperator{S, BR}; prefix::AbstractString="") where {S, BR}
  println(io, prefix, "PureOperator")
  println(io, prefix, "| M: ", string(op.bitmask, base=2, pad=op.hilbert_space.bitoffsets[end]))
  println(io, prefix, "| R: ", string(op.bitrow, base=2, pad=op.hilbert_space.bitoffsets[end]))
  println(io, prefix, "| C: ", string(op.bitcol, base=2, pad=op.hilbert_space.bitoffsets[end]))
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
  n = psi.hilbert_space.bitoffsets[end]
  for b in bs
    println(io, prefix, "| ", string(b, base=2, pad=n), " : ", psi.components[b])
  end
end

function prettyprintln(io ::IO, hsr::HilbertSpaceRealization; prefix::AbstractString="")
  n = hsr.hilbert_space.bitoffsets[end]
  println(io, prefix, "HilbertSpaceRealization")
  for bvec in hsr.basis_list
    println(io, prefix, "| ", string(bvec, base=2, pad=n))
  end
end