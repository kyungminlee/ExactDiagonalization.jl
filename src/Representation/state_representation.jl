function represent(
    hsrep::HilbertSpaceRepresentation{HS, BR, DT},
    state::SparseState{S, BR2}
) where {HS, BR, DT, S, BR2}
    err = zero(S)
    err_sq = zero(real(S))
    out = zeros(S, dimension(hsrep))
    for (bvec, amplitude) in state
        i = get(hsrep.basis_lookup, BR(bvec), -1)
        if i > 0
            out[i] = amplitude
        else
            err += amplitude
            err_sq += abs2(amplitude)
        end
    end
    return (out, err, err_sq)
end


"""
    SparseState(hsrep, state_rep, tol=√eps(Float64))

Make a `SparseState` from a representation
"""
function SparseState(
    hsrep::HilbertSpaceRepresentation{HS, BR, DT},
    state_rep::AbstractVector{S},
    tol::Real=Base.rtoldefault(Float64)
) where {HS, BR, DT, S<:Number}
    if dimension(hsrep) != length(state_rep)
        throw(DimensionMismatch(
            "dimension of the Hilbert space representation ($(dimension(hsrep))) "*
            "does not match the length of the vector ($(length(state_rep)))"
        ))
    end
    out = SparseState{S, BR}()
    for (i, amplitude) in enumerate(state_rep)
        if !isapprox(amplitude, zero(S); atol=tol)
            out[hsrep.basis_list[i]] = amplitude
        end
    end
    return out
end
