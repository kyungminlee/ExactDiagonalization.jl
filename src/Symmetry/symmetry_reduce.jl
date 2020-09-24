export symmetry_reduce, symmetry_unreduce


"""
    symmetry_reduce(hsr, lattice, symmetry_irrep_component, complex_type=ComplexF64, tol=√ϵ)

Symmetry-reduce the HilbertSpaceRepresentation using translation group.

"""
function symmetry_reduce(
    hsr::HilbertSpaceRepresentation{QN, BR, DT},
    ssic::AbstractSymmetryIrrepComponent,
    complex_type::Type{ComplexType}=ComplexF64;
    tol::Real=Base.rtoldefault(Float64)
) where {QN, BR, DT, ComplexType<:Complex}
    symred = Threads.nthreads() == 1 ? symmetry_reduce_serial : symmetry_reduce_parallel
    return symred(hsr, ssic, ComplexType; tol=tol)
end

"""
    symmetry_reduce(rhsr, large_vector)

Reduce a large vector into the reduced hilbert space representation.
Simply throw away components that don't fit.
"""
function symmetry_reduce(
    rhsr::ReducedHilbertSpaceRepresentation{HSR, SIC, BR, C},
    large_vector::AbstractVector{Si}
) where {HSR, SIC, BR, C, Si<:Number}
    if length(large_vector) != dimension(rhsr.parent)
        throw(DimensionMismatch("Dimension of the input vector should match the larger representation"))
    end
    So = promote_type(C, Si)

    # basis mapping
    # (i_p | i_r | ampl) indicates : U_(p, r) = ampl
    n_p = length(rhsr.basis_mapping_index)
    local_small_vectors = [zeros(So, dimension(rhsr)) for tid in 1:Threads.nthreads()]
    #for (i_p, i_r) in enumerate(rhsr.basis_mapping_index)
    Threads.@threads for i_p in 1:n_p
        i_r = rhsr.basis_mapping_index[i_p]
        if i_r > 0
            tid = Threads.threadid()
            ampl = rhsr.basis_mapping_amplitude[i_p]
            # H_r = U† H U
            local_small_vectors[tid][i_r] += conj(ampl) * large_vector[i_p]
        end
    end
    for tid in 2:Threads.nthreads()
        local_small_vectors[1] += local_small_vectors[tid]
    end
    return local_small_vectors[1]
end




raw"""
    symmetry_unreduce

```math
\begin{pmatrix} l_1 \\ l_2 \\ l_3 \\ \vdots \\ l_n \end{pmatrix}
=
\begin{pmatrix}
. & \cdots & . \\
. & \cdots & . \\
. & \cdots & . \\
  & \dots & \\
. & \cdots &
\end{pmatrix}
\begin{pmatrix} s_1 \\ \vdots \\ s_m \end{pmatrix}
```
"""
function symmetry_unreduce(
    rhsr::ReducedHilbertSpaceRepresentation{HSR, SIC, BR, C},
    small_vector::AbstractVector{Si}
) where {HSR, SIC, BR, C, Si<:Number}
    if length(small_vector) != dimension(rhsr)
        throw(DimensionMismatch("Dimension of the input vector should match the reduced representation"))
    end
    So = promote_type(C, Si)
    large_vector = zeros(So, dimension(rhsr.parent))
    for (i_p, i_r) in enumerate(rhsr.basis_mapping_index)
        if i_r > 0
            ampl = rhsr.basis_mapping_amplitude[i_p]
            # H_r = U† H U
            large_vector[i_p] += ampl * small_vector[i_r]
        end
    end
    return large_vector
end
