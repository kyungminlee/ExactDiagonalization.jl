var documenterSearchIndex = {"docs":
[{"location":"hilbertspace/#Hilbert-space-1","page":"Hilbert space","title":"Hilbert space","text":"","category":"section"},{"location":"representation/#Representation-1","page":"Representation","title":"Representation","text":"","category":"section"},{"location":"symmetry/#Symmetry-1","page":"Symmetry","title":"Symmetry","text":"","category":"section"},{"location":"operator/#Operator-1","page":"Operator","title":"Operator","text":"","category":"section"},{"location":"#ExactDiagonalization-1","page":"Home","title":"ExactDiagonalization","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"Implements exact diagonalization","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Modules = [ExactDiagonalization]","category":"page"},{"location":"#ExactDiagonalization.HilbertSpace","page":"Home","title":"ExactDiagonalization.HilbertSpace","text":"HilbertSpace{QN}\n\nAbstract Hilbert space with quantum number type QN.\n\nExamples\n\njulia> using ExactDiagonalization\n\njulia> spin_site = Site{Int64}([State{Int64}(\"Up\", +1), State{Int64}(\"Dn\", -1)])\nSite{Int64}(State{Int64}[State{Int64}(\"Up\", 1), State{Int64}(\"Dn\", -1)])\n\njulia> hs = HilbertSpace{Int64}([spin_site, spin_site])\nHilbertSpace{Int64}(Site{Int64}[Site{Int64}(State{Int64}[State{Int64}(\"Up\", 1), State{Int64}(\"Dn\", -1)]), Site{Int64}(State{Int64}[State{Int64}(\"Up\", 1), State{Int64}(\"Dn\", -1)])], [1, 1], [0, 1, 2])\n\n\n\n\n\n","category":"type"},{"location":"#ExactDiagonalization.Permutation","page":"Home","title":"ExactDiagonalization.Permutation","text":"Permutation(perms ::AbstractVector{Int}; max_cycle=2048)\n\nCreate a permutation of integers from 1 to n. perms should be a permutation of 1:n.\n\nArguments\n\nperms: an integer vector containing a permutation of integers from 1 to n\nmax_cycle: maximum length of cycle\n\n\n\n\n\n","category":"type"},{"location":"#ExactDiagonalization.Site","page":"Home","title":"ExactDiagonalization.Site","text":"Site{QN}\n\nA site with quantum number type QN.\n\nExamples\n\njulia> using ExactDiagonalization\n\njulia> up = State{Int}(\"Up\", 1); dn = State{Int}(\"Dn\", -1);\n\njulia> Site([up, dn])\nSite{Int64}(State{Int64}[State{Int64}(\"Up\", 1), State{Int64}(\"Dn\", -1)])\n\n\n\n\n\n","category":"type"},{"location":"#ExactDiagonalization.SparseState","page":"Home","title":"ExactDiagonalization.SparseState","text":"struct SparseState{Scalar<:Number, BR}\n\nRepresents a row vector. Free.\n\n\n\n\n\n","category":"type"},{"location":"#ExactDiagonalization.State","page":"Home","title":"ExactDiagonalization.State","text":"State{QN}\n\nState with quantum number type QN.\n\nExamples\n\njulia> using ExactDiagonalization, StaticArrays\n\njulia> up = State{Int}(\"Up\", 1)\nState{Int64}(\"Up\", 1)\n\njulia> State(\"Dn\", SVector{2, Int}([-1, 1]))\nState{SArray{Tuple{2},Int64,1,2}}(\"Dn\", [-1, 1])\n\n\n\n\n\n","category":"type"},{"location":"#ExactDiagonalization.apply!-Union{Tuple{BR}, Tuple{S2}, Tuple{S1}, Tuple{SparseState{S1,BR},NullOperator,SparseState{S2,BR}}} where BR where S2 where S1","page":"Home","title":"ExactDiagonalization.apply!","text":"apply!\n\nApply operator to psi and add it to out.\n\n\n\n\n\n","category":"method"},{"location":"#ExactDiagonalization.apply!-Union{Tuple{S2}, Tuple{S1}, Tuple{O}, Tuple{HSR}, Tuple{Array{S1,1},AbstractArray{S2,1},OperatorRepresentation{HSR,O}}} where S2<:Number where S1<:Number where O where HSR","page":"Home","title":"ExactDiagonalization.apply!","text":"apply!(out, opr, state; range=1:size(opr, 1))\n\nPerform out += opr * state. Apply the operator representation opr to the row vector state and add it to the row vector out. Return sum of errors and sum of error-squared. Call apply_serial! if Threads.nthreads() == 1, and apply_parallel! if greater.\n\nArguments\n\nout ::Vector{S1}\nstate ::AbstractVector{S2}\nopr ::OperatorRepresentation{HSR, O}\nrange ::AbstractVector{<:Integer}=1:size(opr, 1)\n\n\n\n\n\n","category":"method"},{"location":"#ExactDiagonalization.apply!-Union{Tuple{S2}, Tuple{S1}, Tuple{O}, Tuple{HSR}, Tuple{Array{S1,1},OperatorRepresentation{HSR,O},AbstractArray{S2,1}}} where S2<:Number where S1<:Number where O where HSR","page":"Home","title":"ExactDiagonalization.apply!","text":"apply!(out, opr, state; range=1:size(opr, 2))\n\nPerform out += opr * state. Apply the operator representation opr to the column vector state and add it to the column vector out. Return sum of errors and sum of error-squared. Call apply_serial! if Threads.nthreads() == 1, and apply_parallel! if greater.\n\nArguments\n\nout ::Vector{S1}\nopr ::OperatorRepresentation{HSR, O}\nstate ::AbstractVector{S2}\nrange ::AbstractVector{<:Integer}=1:size(opr, 2)\n\n\n\n\n\n","category":"method"},{"location":"#ExactDiagonalization.apply_parallel!-Union{Tuple{S2}, Tuple{S1}, Tuple{O}, Tuple{HSR}, Tuple{Array{S1,1},AbstractArray{S2,1},OperatorRepresentation{HSR,O}}} where S2<:Number where S1<:Number where O where HSR","page":"Home","title":"ExactDiagonalization.apply_parallel!","text":"apply_parallel!(out, state, opr; range=1:size(opr, 1))\n\nPerform out += state * opr. Apply the operator representation opr to the row vector state and add it to the row vector out. Return sum of errors and sum of error-squared. Multi-threaded version.\n\nArguments\n\nout ::Vector{S1}\nstate ::AbstractVector{S2}\nopr ::OperatorRepresentation{HSR, O}\nrange ::AbstractVector{<:Integer}=1:dimension(opr.hilbert_space_representation)\n\n\n\n\n\n","category":"method"},{"location":"#ExactDiagonalization.apply_parallel!-Union{Tuple{S2}, Tuple{S1}, Tuple{O}, Tuple{HSR}, Tuple{Array{S1,1},OperatorRepresentation{HSR,O},AbstractArray{S2,1}}} where S2<:Number where S1<:Number where O where HSR","page":"Home","title":"ExactDiagonalization.apply_parallel!","text":"apply_parallel!(out, opr, state; range=1:size(opr, 2))\n\nPerform out += opr * state. Apply the operator representation opr to the column vector state and add it to the column vector out. Return sum of errors and sum of error-squared. Multi-threaded version.\n\nArguments\n\nout ::Vector{S1}\nopr ::OperatorRepresentation{HSR, O}\nstate ::AbstractVector{S2}\nrange ::AbstractVector{<:Integer}=1:size(opr, 2)\n\n\n\n\n\n","category":"method"},{"location":"#ExactDiagonalization.apply_serial!-Union{Tuple{S2}, Tuple{S1}, Tuple{O}, Tuple{HSR}, Tuple{Array{S1,1},AbstractArray{S2,1},OperatorRepresentation{HSR,O}}} where S2<:Number where S1<:Number where O where HSR","page":"Home","title":"ExactDiagonalization.apply_serial!","text":"apply_serial!(out, state, opr; range=1:size(opr, 1))\n\nPerform out += state * opr. Apply the operator representation opr to the row vector state and add it to the row vector out. Return sum of errors and sum of error-squared. Single-threaded version.\n\nArguments\n\nout ::Vector{S1}\nstate ::AbstractVector{S2}\nopr ::OperatorRepresentation{HSR, O}\nrange ::AbstractVector{<:Integer}=1:size(opr, 1)\n\n\n\n\n\n","category":"method"},{"location":"#ExactDiagonalization.apply_serial!-Union{Tuple{S2}, Tuple{S1}, Tuple{O}, Tuple{HSR}, Tuple{Array{S1,1},OperatorRepresentation{HSR,O},AbstractArray{S2,1}}} where S2<:Number where S1<:Number where O where HSR","page":"Home","title":"ExactDiagonalization.apply_serial!","text":"apply_serial!(out, opr, state; range=1:size(opr, 2))\n\nPerform out += opr * state. Apply the operator representation opr to the column vector state and add it to the column vector out. Return sum of errors and sum of error-squared. Single-threaded version.\n\nArguments\n\nout ::Vector{S1}\nopr ::OperatorRepresentation{HSR, O}\nstate ::AbstractVector{S2}\nrange ::AbstractVector{<:Integer}=1:size(opr, 2)\n\n\n\n\n\n","category":"method"},{"location":"#ExactDiagonalization.bitwidth-Tuple{HilbertSpace}","page":"Home","title":"ExactDiagonalization.bitwidth","text":"Total number of bits\n\njulia> using ExactDiagonalization\n\njulia> spin_site = Site{Int64}([State{Int64}(\"Up\", +1), State{Int64}(\"Dn\", -1)])\nSite{Int64}(State{Int64}[State{Int64}(\"Up\", 1), State{Int64}(\"Dn\", -1)])\n\njulia> hs = HilbertSpace{Int64}([spin_site, spin_site, spin_site,])\nHilbertSpace{Int64}(Site{Int64}[Site{Int64}(State{Int64}[State{Int64}(\"Up\", 1), State{Int64}(\"Dn\", -1)]), Site{Int64}(State{Int64}[State{Int64}(\"Up\", 1), State{Int64}(\"Dn\", -1)]), Site{Int64}(State{Int64}[State{Int64}(\"Up\", 1), State{Int64}(\"Dn\", -1)])], [1, 1, 1], [0, 1, 2, 3])\n\njulia> bitwidth(hs)\n3\n\n\n\n\n\n","category":"method"},{"location":"#ExactDiagonalization.bitwidth-Tuple{Site}","page":"Home","title":"ExactDiagonalization.bitwidth","text":"bitwidth(site ::Site)\n\nNumber of bits necessary to represent the states of the given site.\n\n\n\n\n\n","category":"method"},{"location":"#ExactDiagonalization.dimension-Tuple{HilbertSpaceRepresentation}","page":"Home","title":"ExactDiagonalization.dimension","text":"dimension\n\nDimension of the Concrete Hilbert space, i.e. number of basis vectors.\n\n\n\n\n\n","category":"method"},{"location":"#ExactDiagonalization.dimension-Tuple{Site}","page":"Home","title":"ExactDiagonalization.dimension","text":"dimension(site ::Site)\n\nHilbert space dimension of a given site ( = number of states).\n\n\n\n\n\n","category":"method"},{"location":"#ExactDiagonalization.extract-Union{Tuple{U}, Tuple{QN}, Tuple{HilbertSpace{QN},U}} where U<:Unsigned where QN","page":"Home","title":"ExactDiagonalization.extract","text":"Convert binary representation to an array of indices (of states)\n\nExamples ≡≡≡≡≡≡≡≡≡≡\n\n\n\n\n\n\n\n","category":"method"},{"location":"#ExactDiagonalization.get_quantum_number-Union{Tuple{BR}, Tuple{QN}, Tuple{HilbertSpace{QN},BR}} where BR where QN","page":"Home","title":"ExactDiagonalization.get_quantum_number","text":"get_quantum_number\n\n\n\n\n\n","category":"method"},{"location":"#ExactDiagonalization.get_row_iterator-Union{Tuple{O}, Tuple{HSR}, Tuple{OperatorRepresentation{HSR,O},Integer}} where O where HSR","page":"Home","title":"ExactDiagonalization.get_row_iterator","text":"May contain duplicates\n\n\n\n\n\n","category":"method"},{"location":"#ExactDiagonalization.get_state-Union{Tuple{U}, Tuple{Site,U}} where U<:Unsigned","page":"Home","title":"ExactDiagonalization.get_state","text":"get_state(site ::Site{QN}, binrep ::BR) where {QN, BR<:Unsigned}\n\nReturns the state of site represented by the bits binrep.\n\n\n\n\n\n","category":"method"},{"location":"#ExactDiagonalization.is_compatible-Tuple{AbstractArray{#s30,1} where #s30<:Rational,AbstractArray{#s24,1} where #s24<:Integer}","page":"Home","title":"ExactDiagonalization.is_compatible","text":"is_compatible\n\nCheck whether the fractional momentum ([0, 1)ᴺ) compatible with the identity translation. i.e. k¹ R¹ + k² R² + ... + kᴺ Rᴺ = 0 (mod 1)\n\nArguments\n\nfractional_momentum ::AbstractVector{Rational} : k\nidentity_translation ::AbstractVector{<:Integer} : R\n\n\n\n\n\n","category":"method"},{"location":"#ExactDiagonalization.qntype-Union{Tuple{HilbertSpace{QN}}, Tuple{QN}} where QN","page":"Home","title":"ExactDiagonalization.qntype","text":"qntype\n\n\n\n\n\n","category":"method"},{"location":"#ExactDiagonalization.quantum_number_sectors-Union{Tuple{HilbertSpace{QN}}, Tuple{QN}} where QN","page":"Home","title":"ExactDiagonalization.quantum_number_sectors","text":"quantum_number_sectors\n\n\n\n\n\n","category":"method"},{"location":"#ExactDiagonalization.represent-Union{Tuple{HilbertSpaceSector{QN}}, Tuple{QN}} where QN","page":"Home","title":"ExactDiagonalization.represent","text":"represent(hs; BR ::DataType=UInt)\n\nMake a HilbertSpaceRepresentation with all the basis vectors of the specified HilbertSpace.\n\nArguments\n\nhs ::HilbertSpace{QN}: Abstract Hilbert space\nallowed: Allowed quantum numbers\nBR ::DataType=UInt: Binary representation type\n\n\n\n\n\n","category":"method"},{"location":"#ExactDiagonalization.represent-Union{Tuple{HilbertSpace{QN}}, Tuple{QN}} where QN","page":"Home","title":"ExactDiagonalization.represent","text":"represent(hs; BR ::DataType=UInt)\n\nMake a HilbertSpaceRepresentation with all the basis vectors of the specified HilbertSpace.\n\nArguments\n\nhs ::HilbertSpace{QN}: Abstract Hilbert space\nBR ::DataType=UInt: Binary representation type\n\n\n\n\n\n","category":"method"},{"location":"#ExactDiagonalization.symmetry_reduce-Union{Tuple{Si}, Tuple{C}, Tuple{BR}, Tuple{HSR}, Tuple{ReducedHilbertSpaceRepresentation{HSR,BR,C},AbstractArray{Si,1}}} where Si<:Number where C where BR where HSR","page":"Home","title":"ExactDiagonalization.symmetry_reduce","text":"\n\n\n\n","category":"method"},{"location":"#ExactDiagonalization.symmetry_reduce_serial-Union{Tuple{BR}, Tuple{QN}, Tuple{HilbertSpaceRepresentation{QN,BR},TranslationGroup,AbstractArray{#s88,1} where #s88<:Rational}} where BR where QN","page":"Home","title":"ExactDiagonalization.symmetry_reduce_serial","text":"symmetry_reduce_serial(hsr, trans_group, frac_momentum; ComplexType=ComplexF64, tol=sqrt(eps(Float64)))\n\nSymmetry-reduce the HilbertSpaceRepresentation using translation group.\n\n\n\n\n\n","category":"method"},{"location":"#Base.:*-Tuple{Permutation,Permutation}","page":"Home","title":"Base.:*","text":"*(p1 ::Permutation, p2 ::Permutation)\n\nMultiply the two permutation. Return [p2.map[ p1.map[x] ] for x in 1:length(p1.map)].\n\nExamples\n\njulia> using ExactDiagonalization\n\njulia> Permutation([1,3,2]) * Permutation([2,1,3])\nPermutation([2, 3, 1], 3)\n\njulia> Permutation([2,1,3]) * Permutation([1,3,2])\nPermutation([3, 1, 2], 3)\n\n\n\n\n\n","category":"method"},{"location":"#Base.:^-Tuple{Permutation,Integer}","page":"Home","title":"Base.:^","text":"^(perm ::Permutation, pow ::Integer)\n\nExponentiate the permutation.\n\nExamples\n\n```jldoctest julia> using ExactDiagonalization\n\njulia> Permutation([2,3,4,1])^2 Permutation([3, 4, 1, 2], 2)\n\n\n\n\n\n","category":"method"},{"location":"#ExactDiagonalization.splitblock-Tuple{Integer,Integer}","page":"Home","title":"ExactDiagonalization.splitblock","text":"splitblock\n\nSplit n into b blocks.\n\nArguments\n\nn ::Integer: the number of elements to split.\nb ::Integer: the number of blocks.\n\n\n\n\n\n","category":"method"}]
}
