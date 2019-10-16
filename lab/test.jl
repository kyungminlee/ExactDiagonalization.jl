include("../src/util.jl")

L=zero(UInt64)
R=one(UInt64) << 30 - 1
basislist = UInt64[]
for b in L:R
    if count_ones(b) == 15
        push!(basislist, b)
    end
    if b % 10000000 == 0
        println(b, "\t", length(basislist))
    end
end

