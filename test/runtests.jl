using Test
using ExactDiagonalization

include("test_util.jl")

include("test_hilbert_space.jl")
include("test_hilbert_space_sector.jl")
include("test_operator.jl")
#include("test_sparse_state.jl")
include("test_operator_simplify.jl")
include("test_operator_application.jl")

include("test_hilbert_space_representation.jl")
include("test_operator_representation.jl")

include("test_symmetry.jl")

include("test_prettyprintln.jl")
