var documenterSearchIndex = {"docs":
[{"location":"hilbertspace/#Hilbert-space-1","page":"Hilbert space","title":"Hilbert space","text":"","category":"section"},{"location":"representation/#Representation-1","page":"Representation","title":"Representation","text":"","category":"section"},{"location":"symmetry/#Symmetry-1","page":"Symmetry","title":"Symmetry","text":"","category":"section"},{"location":"operator/#Operator-1","page":"Operator","title":"Operator","text":"","category":"section"},{"location":"#ExactDiagonalization-1","page":"Home","title":"ExactDiagonalization","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"Implements exact diagonalization","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Modules = [ExactDiagonalization]","category":"page"},{"location":"#ExactDiagonalization.HilbertSpace","page":"Home","title":"ExactDiagonalization.HilbertSpace","text":"HilbertSpace{QN}\n\nAbstract Hilbert space with quantum number type QN.\n\nExamples\n\njulia> using ExactDiagonalization\n\njulia> spin_site = Site{Int64}([State{Int64}(\"Up\", +1), State{Int64}(\"Dn\", -1)])\nSite{Int64}(State{Int64}[State{Int64}(\"Up\", 1), State{Int64}(\"Dn\", -1)])\n\njulia> hs = HilbertSpace{Int64}([spin_site, spin_site])\nHilbertSpace{Int64}(Site{Int64}[Site{Int64}(State{Int64}[State{Int64}(\"Up\", 1), State{Int64}(\"Dn\", -1)]), Site{Int64}(State{Int64}[State{Int64}(\"Up\", 1), State{Int64}(\"Dn\", -1)])], [1, 1], [0, 1, 2])\n\n\n\n\n\n","category":"type"},{"location":"#ExactDiagonalization.Site","page":"Home","title":"ExactDiagonalization.Site","text":"Site{QN}\n\nA site with quantum number type QN.\n\nExamples\n\njulia> using ExactDiagonalization\n\njulia> up = State{Int}(\"Up\", 1); dn = State{Int}(\"Dn\", -1);\n\njulia> Site([up, dn])\nSite{Int64}(State{Int64}[State{Int64}(\"Up\", 1), State{Int64}(\"Dn\", -1)])\n\n\n\n\n\n","category":"type"},{"location":"#ExactDiagonalization.SparseState","page":"Home","title":"ExactDiagonalization.SparseState","text":"struct SparseState{Scalar<:Number, BR}\n\nRepresents a row vector. Free.\n\n\n\n\n\n","category":"type"},{"location":"#ExactDiagonalization.State","page":"Home","title":"ExactDiagonalization.State","text":"State{QN}\n\nState with quantum number type QN.\n\nExamples\n\njulia> using ExactDiagonalization, StaticArrays\n\njulia> up = State{Int}(\"Up\", 1)\nState{Int64}(\"Up\", 1)\n\njulia> State(\"Dn\", SVector{2, Int}([-1, 1]))\nState{SArray{Tuple{2},Int64,1,2}}(\"Dn\", [-1, 1])\n\n\n\n\n\n","category":"type"},{"location":"#ExactDiagonalization.apply!-Union{Tuple{BR}, Tuple{S2}, Tuple{S1}, Tuple{SparseState{S1,BR},NullOperator,SparseState{S2,BR}}} where BR where S2 where S1","page":"Home","title":"ExactDiagonalization.apply!","text":"apply!\n\nApply operator to psi and add it to out.\n\n\n\n\n\n","category":"method"},{"location":"#ExactDiagonalization.bitwidth-Tuple{HilbertSpace}","page":"Home","title":"ExactDiagonalization.bitwidth","text":"Total number of bits\n\njulia> using ExactDiagonalization\n\njulia> spin_site = Site{Int64}([State{Int64}(\"Up\", +1), State{Int64}(\"Dn\", -1)])\nSite{Int64}(State{Int64}[State{Int64}(\"Up\", 1), State{Int64}(\"Dn\", -1)])\n\njulia> hs = HilbertSpace{Int64}([spin_site, spin_site, spin_site,])\nHilbertSpace{Int64}(Site{Int64}[Site{Int64}(State{Int64}[State{Int64}(\"Up\", 1), State{Int64}(\"Dn\", -1)]), Site{Int64}(State{Int64}[State{Int64}(\"Up\", 1), State{Int64}(\"Dn\", -1)]), Site{Int64}(State{Int64}[State{Int64}(\"Up\", 1), State{Int64}(\"Dn\", -1)])], [1, 1, 1], [0, 1, 2, 3])\n\njulia> bitwidth(hs)\n3\n\n\n\n\n\n","category":"method"},{"location":"#ExactDiagonalization.bitwidth-Tuple{Site}","page":"Home","title":"ExactDiagonalization.bitwidth","text":"bitwidth(site ::Site)\n\nNumber of bits necessary to represent the states of the given site.\n\n\n\n\n\n","category":"method"},{"location":"#ExactDiagonalization.dimension-Tuple{HilbertSpaceRepresentation}","page":"Home","title":"ExactDiagonalization.dimension","text":"dimension\n\nDimension of the Concrete Hilbert space, i.e. number of basis vectors.\n\n\n\n\n\n","category":"method"},{"location":"#ExactDiagonalization.dimension-Tuple{Site}","page":"Home","title":"ExactDiagonalization.dimension","text":"dimension(site ::Site)\n\nHilbert space dimension of a given site ( = number of states).\n\n\n\n\n\n","category":"method"},{"location":"#ExactDiagonalization.extract-Union{Tuple{U}, Tuple{QN}, Tuple{HilbertSpace{QN},U}} where U<:Unsigned where QN","page":"Home","title":"ExactDiagonalization.extract","text":"Convert binary representation to an array of indices (of states)\n\nExamples ≡≡≡≡≡≡≡≡≡≡\n\n\n\n\n\n\n\n","category":"method"},{"location":"#ExactDiagonalization.get_quantum_number-Union{Tuple{BR}, Tuple{QN}, Tuple{HilbertSpace{QN},BR}} where BR where QN","page":"Home","title":"ExactDiagonalization.get_quantum_number","text":"get_quantum_number\n\n\n\n\n\n","category":"method"},{"location":"#ExactDiagonalization.get_row_iterator-Union{Tuple{O}, Tuple{HSR}, Tuple{OperatorRepresentation{HSR,O},Integer}} where O where HSR","page":"Home","title":"ExactDiagonalization.get_row_iterator","text":"May contain duplicates\n\n\n\n\n\n","category":"method"},{"location":"#ExactDiagonalization.get_state-Union{Tuple{U}, Tuple{Site,U}} where U<:Unsigned","page":"Home","title":"ExactDiagonalization.get_state","text":"get_state(site ::Site{QN}, binrep ::BR) where {QN, BR<:Unsigned}\n\nReturns the state of site represented by the bits binrep.\n\n\n\n\n\n","category":"method"},{"location":"#ExactDiagonalization.is_compatible-Tuple{AbstractArray{#s56,1} where #s56<:Rational,AbstractArray{#s30,1} where #s30<:Integer}","page":"Home","title":"ExactDiagonalization.is_compatible","text":"is_compatible\n\nCheck whether the fractional momentum ([0, 1)ᴺ) compatible with the identity translation. i.e. k¹ R¹ + k² R² + ... + kᴺ Rᴺ = 0 (mod 1)\n\nArguments\n\nfractional_momentum ::AbstractVector{Rational} : k\nidentity_translation ::AbstractVector{<:Integer} : R\n\n\n\n\n\n","category":"method"},{"location":"#ExactDiagonalization.qntype-Union{Tuple{HilbertSpace{QN}}, Tuple{QN}} where QN","page":"Home","title":"ExactDiagonalization.qntype","text":"qntype\n\n\n\n\n\n","category":"method"},{"location":"#ExactDiagonalization.quantum_number_sectors-Union{Tuple{HilbertSpace{QN}}, Tuple{QN}} where QN","page":"Home","title":"ExactDiagonalization.quantum_number_sectors","text":"quantum_number_sectors\n\n\n\n\n\n","category":"method"},{"location":"#ExactDiagonalization.represent-Union{Tuple{HilbertSpaceSector{QN}}, Tuple{QN}} where QN","page":"Home","title":"ExactDiagonalization.represent","text":"represent(hs; BR ::DataType=UInt)\n\nMake a HilbertSpaceRepresentation with all the basis vectors of the specified HilbertSpace.\n\nArguments\n\nhs ::HilbertSpace{QN}: Abstract Hilbert space\nallowed: Allowed quantum numbers\nBR ::DataType=UInt: Binary representation type\n\n\n\n\n\n","category":"method"},{"location":"#ExactDiagonalization.represent-Union{Tuple{HilbertSpace{QN}}, Tuple{QN}} where QN","page":"Home","title":"ExactDiagonalization.represent","text":"represent(hs; BR ::DataType=UInt)\n\nMake a HilbertSpaceRepresentation with all the basis vectors of the specified HilbertSpace.\n\nArguments\n\nhs ::HilbertSpace{QN}: Abstract Hilbert space\nBR ::DataType=UInt: Binary representation type\n\n\n\n\n\n","category":"method"},{"location":"#ExactDiagonalization.splitblock-Tuple{Integer,Integer}","page":"Home","title":"ExactDiagonalization.splitblock","text":"splitblock\n\nSplit n into b blocks\n\n\n\n\n\n","category":"method"}]
}
