module ExactDiagonalization

using StaticArrays
using SparseArrays
using Memento

const LOGGER = getlogger(@__MODULE__)

function __init__()
  Memento.register(LOGGER)
end

include("util.jl")
include("frozensortedarray.jl")

include("HilbertSpace.jl")

include("Operator.jl")

include("Representation.jl")

include("Symmetry.jl")

include("prettyprintln.jl")

end
