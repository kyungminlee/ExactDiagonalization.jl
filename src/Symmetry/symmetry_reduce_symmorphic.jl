export symmetry_reduce_serial, symmetry_reduce_parallel




"""
    symmetry_reduce_serial(hsr, trans_group, frac_momentum, complex_type=ComplexF64, tol=sqrt(eps(Float64)))

Symmetry-reduce the HilbertSpaceRepresentation using translation group (single threaded).

"""
function symmetry_reduce_serial(
        hsr::HilbertSpaceRepresentation{QN, BR, DT},
        ssic::SymmorphicIrrepComponent{SymmetryEmbedding{TranslationSymmetry}, SymmetryEmbedding{PointSymmetry}},
        complex_type::Type{ComplexType}=ComplexF64;
        tol::Real=Base.rtoldefault(Float64)
        ) where {QN, BR, DT, ComplexType<:Complex}

    HSR = HilbertSpaceRepresentation{QN, BR, DT}

    n_basis = length(hsr.basis_list)

    basis_mapping_representative = Vector{Int}(undef, n_basis)
    fill!(basis_mapping_representative, -1)
    basis_mapping_amplitude = zeros(ComplexType, n_basis)

    tsym_group_size = group_order(ssic.component1.symmetry)
    psym_group_size = group_order(ssic.component2.symmetry)
    group_size = tsym_group_size * psym_group_size

    size_estimate = n_basis ÷ max(1, group_size - 1)

    reduced_basis_list = BR[]
    sizehint!(reduced_basis_list, size_estimate)

    visited = falses(n_basis)

    basis_states = Matrix{BR}(undef, tsym_group_size, psym_group_size)  # recomputed for every bvec
    basis_amplitudes = Dict{BR, ComplexType}()
    sizehint!(basis_amplitudes, group_size + group_size ÷ 2)

    tsym_symops_and_amplitudes = [(x, conj(y)) for (x, y) in get_irrep_iterator(ssic.component1)]
    psym_symops_and_amplitudes = [(x, conj(y)) for (x, y) in get_irrep_iterator(ssic.component2)]
    @assert length(tsym_symops_and_amplitudes) == tsym_group_size
    @assert length(psym_symops_and_amplitudes) == psym_group_size

    zeroC = zero(ComplexF64)
    for ivec_p in 1:n_basis
        visited[ivec_p] && continue
        bvec = hsr.basis_list[ivec_p]

        # basis_states[1,1] = bvec
        #  TODO: think about order tsym first or psym first
        #  I think translation needs to be done first, since we are considering
        #  the Bloch states and applying point operations to the Bloch states.
        for it in 1:tsym_group_size
            (tsym_symop, tsym_ampl) = tsym_symops_and_amplitudes[it]
            # isapprox(tsym_ampl, 0; atol=tol) && continue
            bvec_prime_t = symmetry_apply(hsr.hilbert_space, tsym_symop, bvec)
            for ip in 1:psym_group_size
                # (it, ip) == (1,1) && continue
                (psym_symop, psym_ampl) = psym_symops_and_amplitudes[ip]
                # isapprox(psym_ampl, 0; atol=tol) && continue
                bvec_prime = symmetry_apply(hsr.hilbert_space, psym_symop, bvec_prime_t)
                basis_states[it, ip] = bvec_prime
            end # for ip
        end # for it

        empty!(basis_amplitudes)
        for it in 1:tsym_group_size
            (_, tsym_ampl) = tsym_symops_and_amplitudes[it]
            isapprox(tsym_ampl, 0; atol=tol) && continue
            for ip in 1:psym_group_size
                (_, psym_ampl) = psym_symops_and_amplitudes[ip]
                isapprox(psym_ampl, 0; atol=tol) && continue
                bvec_prime = basis_states[it, ip]
                ampl = tsym_ampl * psym_ampl
                basis_amplitudes[bvec_prime] = get(basis_amplitudes, bvec_prime, zeroC) + ampl
            end
        end
        choptol!(basis_amplitudes, tol)
        !haskey(basis_amplitudes, bvec) && continue
        minimum(keys(basis_amplitudes)) != bvec && continue

        push!(reduced_basis_list, bvec)
        inv_norm = one(Float64) / sqrt(sum(abs2.(values(basis_amplitudes))))
        for (bvec_prime, amplitude) in basis_amplitudes
            ivec_p_prime = hsr.basis_lookup[bvec_prime]
            visited[ivec_p_prime] = true
            basis_mapping_representative[ivec_p_prime] = ivec_p
            basis_mapping_amplitude[ivec_p_prime] = amplitude * inv_norm
        end
    end # for ivec_p

    basis_mapping_index = Vector{Int}(undef, n_basis)
    fill!(basis_mapping_index, -1)

    for (ivec_r, bvec) in enumerate(reduced_basis_list)
        ivec_p = hsr.basis_lookup[bvec]
        basis_mapping_index[ivec_p] = ivec_r
    end

    for (ivec_p_prime, ivec_p) in enumerate(basis_mapping_representative)
        (ivec_p <= 0) && continue  # not in this irrep
        (ivec_p_prime == ivec_p) && continue  # already in the lookup
        ivec_r = basis_mapping_index[ivec_p]
        basis_mapping_index[ivec_p_prime] = ivec_r
    end

    return ReducedHilbertSpaceRepresentation{HSR, SymmorphicIrrepComponent{SymmetryEmbedding{TranslationSymmetry}, SymmetryEmbedding{PointSymmetry}}, BR, ComplexType}(
                hsr, ssic, reduced_basis_list,
                basis_mapping_index, basis_mapping_amplitude)
end




"""
    symmetry_reduce_parallel(hsr, trans_group, frac_momentum, complex_type=ComplexF64, tol=sqrt(eps(Float64)))

Symmetry-reduce the HilbertSpaceRepresentation using translation group (multi-threaded).

"""
function symmetry_reduce_parallel(
        hsr::HilbertSpaceRepresentation{QN, BR, DT},
        ssic::SymmorphicIrrepComponent{SymmetryEmbedding{TranslationSymmetry}, SymmetryEmbedding{PointSymmetry}},
        complex_type::Type{ComplexType}=ComplexF64;
        tol::Real=Base.rtoldefault(Float64)
        ) where {QN, BR, DT, ComplexType<:Complex}

    HSR = HilbertSpaceRepresentation{QN, BR, DT}
    @debug "BEGIN symmetry_reduce_parallel"

    n_basis = length(hsr.basis_list)
    @debug "Original Hilbert space dimension: $n_basis"

    # basis_mapping_index and basis_mapping_amplitude contain information about
    # which basis vector of the larger Hilbert space is included
    # in the basis of the smaller Hilbert space with what amplitude.
    basis_mapping_representative = Vector{Int}(undef, n_basis)
    fill!(basis_mapping_representative, -1)
    basis_mapping_amplitude = zeros(ComplexType, n_basis)

    tsym_group_size = group_order(ssic.component1.symmetry)
    psym_group_size = group_order(ssic.component2.symmetry)
    group_size = tsym_group_size * psym_group_size

    size_estimate = n_basis ÷ max(1, group_size - 1)
    @debug "Estimate for the reduced Hilbert space dimension: $size_estimate"

    nthreads = Threads.nthreads()
    local_reduced_basis_list = Vector{Vector{BR}}(undef, nthreads)
    for i in 1:nthreads
        local_reduced_basis_list[i] = BR[]
        sizehint!(local_reduced_basis_list[i], size_estimate ÷ nthreads + 1)
    end

    visited = zeros(UInt8, n_basis) # use UInt8 rather than Bool for thread safety

    local_basis_states = Array{BR, 3}(undef, (nthreads, tsym_group_size, psym_group_size))
    local_basis_amplitudes = Vector{Dict{BR, ComplexType}}(undef, nthreads)
    for id in 1:nthreads
        local_basis_amplitudes[id] = Dict{BR, ComplexType}()
        sizehint!(local_basis_amplitudes[id], group_size)
    end

    # Load balancing
    #   The representatives are the smaller binary numbers.
    #   Distribute them equally among threads.
    reorder = Int[]
    sizehint!(reorder, n_basis)
    nblocks = (n_basis + nthreads - 1) ÷ nthreads
    for i in 1:nthreads, j in 1:nblocks
        k = i + nthreads * (j-1)
        if 1 <= k <= n_basis
            push!(reorder, k)
        end
    end

    # tsic = ssic.translation
    # tsym = tsic.symmetry
    # tsym_irrep_index = tsic.irrep_index
    # # whether the element is in the kernel of the representation
    # orthogonal_momentum = tsym.orthogonal_coordinates[tsym_irrep_index]
    # orthogonal_shape = tsym.orthogonal_shape
    # is_identity = [iscompatible(orthogonal_momentum, orthogonal_shape, t)
    #                for t in tsym.orthogonal_coordinates]
    tsym_symops_and_amplitudes = [(x, conj(y)) for (x, y) in get_irrep_iterator(ssic.component1)]
    psym_symops_and_amplitudes = [(x, conj(y)) for (x, y) in get_irrep_iterator(ssic.component2)]

    @assert length(tsym_symops_and_amplitudes) == tsym_group_size
    @assert length(psym_symops_and_amplitudes) == psym_group_size

    @debug "Starting reduction (parallel)"
    zeroC = zero(ComplexF64)
    Threads.@threads for itemp in 1:n_basis
        ivec_p = reorder[itemp]
        (visited[ivec_p] != 0x0) && continue

        id = Threads.threadid()
        bvec = hsr.basis_list[ivec_p]

        # A basis binary representation is incompatible with the reduced Hilbert space if
        # (1) it is not the smallest among the its star, or
        # (2) its star is smaller than the representation
        local_basis_states[id, 1, 1] = bvec
        for it in 1:tsym_group_size
            (tsym_symop, tsym_ampl) = tsym_symops_and_amplitudes[it]
            # isapprox(tsym_ampl, 0; atol=tol) && continue
            bvec_prime_t = symmetry_apply(hsr.hilbert_space, tsym_symop, bvec)
            for ip in 1:psym_group_size
                # (it, ip) == (1, 1) && continue
                (psym_symop, psym_ampl) = psym_symops_and_amplitudes[ip]
                # isapprox(psym_ampl, 0; atol=tol) && continue
                bvec_prime = symmetry_apply(hsr.hilbert_space, psym_symop, bvec_prime_t)
                local_basis_states[id, it, ip] = bvec_prime
            end # for ip
        end # for it

        empty!(local_basis_amplitudes[id])
        for it in 1:tsym_group_size
            (_, tsym_ampl) = tsym_symops_and_amplitudes[it]
            isapprox(tsym_ampl, 0; atol=tol) && continue
            for ip in 1:psym_group_size
                (_, psym_ampl) = psym_symops_and_amplitudes[ip]
                isapprox(psym_ampl, 0; atol=tol) && continue
                bvec_prime = local_basis_states[id, it, ip]
                ampl = tsym_ampl * psym_ampl
                local_basis_amplitudes[id][bvec_prime] = get(local_basis_amplitudes[id], bvec_prime, zeroC) + ampl
            end
        end
        choptol!(local_basis_amplitudes[id], tol)
        !haskey(local_basis_amplitudes[id], bvec) && continue
        minimum(keys(local_basis_amplitudes[id])) != bvec && continue

        push!(local_reduced_basis_list[id], bvec)
        inv_norm = one(Float64) / sqrt(sum(abs2.(values(local_basis_amplitudes[id]))))
        for (bvec_prime, amplitude) in local_basis_amplitudes[id]
            ivec_p_prime = hsr.basis_lookup[bvec_prime]
            visited[ivec_p_prime] = 0x1
            basis_mapping_representative[ivec_p_prime] = ivec_p
            basis_mapping_amplitude[ivec_p_prime] = amplitude * inv_norm
        end
    end
    @debug "Finished reduction (parallel)"

    reduced_basis_list = BR[]
    sizehint!(reduced_basis_list, sum(length(x) for x in local_reduced_basis_list))
    while !isempty(local_reduced_basis_list)
        lbl = pop!(local_reduced_basis_list)
        append!(reduced_basis_list, lbl)
    end
    @debug "Collected basis list"

    sort!(reduced_basis_list)
    @debug "Sorted basis list"

    # Basis vectors of the unreduced Hilbert space that are
    # not included in the reduced Hilbert space are marked as -1
    basis_mapping_index = Vector{Int}(undef, n_basis)
    fill!(basis_mapping_index, -1)

    Threads.@threads for ivec_r in eachindex(reduced_basis_list)
        bvec = reduced_basis_list[ivec_r]
        ivec_p = hsr.basis_lookup[bvec]
        basis_mapping_index[ivec_p] = ivec_r
    end
    @debug "Collected basis lookup (diagonal)"

    Threads.@threads for ivec_p_prime in eachindex(basis_mapping_representative)
        ivec_p = basis_mapping_representative[ivec_p_prime]
        (ivec_p <= 0) && continue  # not in this irrep
        (ivec_p_prime == ivec_p) && continue  # already in the lookup
        ivec_r = basis_mapping_index[ivec_p]
        basis_mapping_index[ivec_p_prime] = ivec_r
    end
    @debug "Collected basis lookup (offdiagonal)"

    @debug "END symmetry_reduce_parallel"
    return ReducedHilbertSpaceRepresentation{HSR, SymmorphicIrrepComponent{SymmetryEmbedding{TranslationSymmetry}, SymmetryEmbedding{PointSymmetry}}, BR, ComplexType}(
                hsr, ssic, reduced_basis_list,
                basis_mapping_index, basis_mapping_amplitude)
end
