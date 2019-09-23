export materialize, materialize_parallel

# function materialize(
#         hsb ::ConcreteHilbertSpaceBlock{BinRep, QN},
#         ops ::AbstractArray{KroneckerProductOperator{OS}};
#         Scalar::DataType=OS) where {BinRep, QN, OS<:Number}
#   hs = hsb.hilbert_space
#   rows = Int[]
#   cols = Int[]
#   vals = Scalar[]
#   err = 0.0
#   for (irow, row) in enumerate(hsb.basis_list)
#     out = apply(hs, ops, row)
#     for (col, amplitude) in out
#       if ! isapprox(amplitude, 0)
#         if !haskey(hsb.basis_lookup, col)
#           err += abs(amplitude)^2
#         else 
#           icol = hsb.basis_lookup[col]
#           push!(rows, irow)
#           push!(cols, icol)
#           push!(vals, amplitude)
#         end        
#       end
#     end
#   end
  
#   if Scalar <:Complex
#     if isapprox( norm(imag.(vals)), 0)
#       vals = real.(vals)
#     end
#   end
#   n = dimension(hsb)
#   return (sparse(rows, cols, vals, n, n), err)
# end



# function materialize(
#         hsb ::ConcreteHilbertSpace{QN, BR},
#         ops ::AbstractArray{KroneckerProductOperator{OS}};
#         Scalar::DataType=OS) where {QN, BR<:Unsigned, OS<:Number}
#   hs = hsb.hilbert_space
#   rows = Int[]
#   cols = Int[]
#   vals = Scalar[]
#   err = 0.0
#   for (irow, row) in enumerate(hsb.basis_list)
#     out = apply(hs, ops, row)
#     for (col, amplitude) in out
#       if ! isapprox(amplitude, 0)
#         if !haskey(hsb.basis_lookup, col)
#           err += abs(amplitude)^2
#         else 
#           icol = hsb.basis_lookup[col]
#           push!(rows, irow)
#           push!(cols, icol)
#           push!(vals, amplitude)
#         end
#       end
#     end
#   end
  
#   if Scalar <:Complex
#     if isapprox( maximum(abs.(imag.(vals))), 0)
#       vals = real.(vals)
#     end
#   end
#   n = dimension(hsb)
#   return (sparse(rows, cols, vals, n, n), err)
# end



# function materialize_parallel(
#         hsb ::ConcreteHilbertSpace{QN, BR},
#         ops ::AbstractArray{KroneckerProductOperator{OS}};
#         Scalar::DataType=OS) where {QN, BR<:Unsigned, OS<:Number}
#   hs = hsb.hilbert_space
#   err = 0.0

#   nthreads = Threads.nthreads()
#   local_rows = [ Int[] for i in 1:nthreads]
#   local_cols = [ Int[] for i in 1:nthreads]
#   local_vals = [ Scalar[] for i in 1:nthreads]
#   local_err = [0.0 for i in 1:nthreads]

#   n_basis = length(hsb.basis_list)

#   Threads.@threads for irow in 1:n_basis
#     id = Threads.threadid()
#     row = hsb.basis_list[irow]
#     out = apply(hs, ops, row)
#     for (col, amplitude) in out
#       if ! isapprox(amplitude, 0)
#         if !haskey(hsb.basis_lookup, col)
#           local_err[id] += abs(amplitude)^2
#         else 
#           icol = hsb.basis_lookup[col]
#           push!(local_rows[id], irow)
#           push!(local_cols[id], icol)
#           push!(local_vals[id], amplitude)
#         end
#       end
#     end
#   end
#   rows = vcat(local_rows...)
#   cols = vcat(local_cols...)
#   vals = vcat(local_vals...)
#   err = sum(local_err)
  
#   if Scalar <:Complex
#     if isapprox( norm(imag.(vals)), 0)
#       vals = real.(vals)
#     end
#   end
#   n = dimension(hsb)
#   return (sparse(rows, cols, vals, n, n), err)
# end

function materialize(
  chs ::ConcreteHilbertSpace{QN, BR},
  nop ::NullOperator) where {QN, BR<:Unsigned}
  n = dimension(chs)
  return (sparse([], [], Float64[], n, n), 0.0)
end


function materialize_parallel(
  chs ::ConcreteHilbertSpace{QN, BR},
  nop ::NullOperator) where {QN, S<:Number, BR<:Unsigned}

  n = dimension(chs)
  return (sparse([], [], Float64[], n, n), 0.0)
end


function materialize(
  chs ::ConcreteHilbertSpace{QN, BR},
  sumop ::SumOperator{S, BR}) where {QN, S<:Number, BR<:Unsigned}
  hs = chs.hilbert_space
  rows = Int[]
  cols = Int[]
  vals = S[]
  err = 0.0
  
  for (irow, row) in enumerate(chs.basis_list)
    ψrow = SparseState{S, BR}(hs, row)
    ψcol = SparseState{S, BR}(hs)
    apply_unsafe!(ψcol, sumop, ψrow)    
    for (col, amplitude) in ψcol.components
      isapprox(amplitude, 0) && continue
      
      if !haskey(chs.basis_lookup, col)
        err += abs(amplitude)^2
      else 
        icol = chs.basis_lookup[col]
        push!(rows, irow)
        push!(cols, icol)
        push!(vals, amplitude)
      end
    end
  end

  if isempty(vals)
    vals = Float64[]
  elseif S <:Complex && isapprox( maximum(abs.(imag.(vals))), 0)
    vals = real.(vals)
  end
  n = dimension(chs)
  return (sparse(rows, cols, vals, n, n), err)
end


function materialize_parallel(
  chs ::ConcreteHilbertSpace{QN, BR},
  sumop ::SumOperator{S, BR}) where {QN, S<:Number, BR<:Unsigned}
  hs = chs.hilbert_space

  nthreads = Threads.nthreads()
  local_rows = [ Int[] for i in 1:nthreads]
  local_cols = [ Int[] for i in 1:nthreads]
  local_vals = [ S[] for i in 1:nthreads]
  local_err =  Float64[0.0 for i in 1:nthreads]

  n_basis = length(chs.basis_list)

  Threads.@threads for irow in 1:n_basis
    id = Threads.threadid()
    row = chs.basis_list[irow]

    ψrow = SparseState{S, BR}(hs, row)
    ψcol = SparseState{S, BR}(hs)
    apply_unsafe!(ψcol, sumop, ψrow)
    for (col, amplitude) in ψcol.components
      isapprox(amplitude, 0) && continue
      
      if !haskey(chs.basis_lookup, col)
        local_err[id] += abs(amplitude)^2
      else 
        icol = chs.basis_lookup[col]
        push!(local_rows[id], irow)
        push!(local_cols[id], icol)
        push!(local_vals[id], amplitude)
      end
    end
  end

  rows ::Vector{Int} = vcat(local_rows...) 
  cols ::Vector{Int} = vcat(local_cols...) 
  vals ::Vector{S} = vcat(local_vals...) 
  err ::Float64 = sum(local_err) 

  if isempty(vals)
    vals = Float64[]
  elseif S <:Complex && isapprox( maximum(abs.(imag.(vals))), 0)
    vals = real.(vals)
  end
  n = dimension(chs)
  return (sparse(rows, cols, vals, n, n), err)
end

