struct PMA{T} <: Factorization{T}
	U::Matrix{T}
	S::Vector{T}
	V::Matrix{T}
	VV::Matrix{T}
end

simplexkernel(n::Integer, w::Real=1.0) = w^2/(n*(n+1))*(ones(n,n)+I)



function simplexgraph2kernelmatrix(G::AbstractMatrix{Bool})
	N = size(G,1)
	# convert graph to simplex matrix
	K = zeros(N,N) # matrix describing all simplexes
	for j=1:size(G,2)
		ind = findall(G[:,j])
		K[ind,ind] .+= simplexkernel(length(ind))
	end
	K
end


# sparse input and output
function simplexgraph2kernelmatrix(G::AbstractSparseMatrix{Bool})
	N = size(G,1)
	# convert graph to simplex matrix

	I,J = Int[],Int[]
	V = Float64[]
	for j=1:size(G,2)
		ind = findall(G[:,j])

		# TODO: Could be built directly without constructing temporary simplexkernel matrix
		S = simplexkernel(length(ind))
		for c in 1:length(ind)
			append!(I, ind)
			append!(J, Iterators.repeated(ind[c],length(ind)))
			append!(V, S[:,c])
		end
	end
	sparse(I,J,V) # sums over repeated indices
end



function simplexgraph2kernelmatrixroot(G::AbstractMatrix{Bool})
	K = simplexgraph2kernelmatrix(G)
	# Factor K into AᵀA where A is symmetric too.
	# F = eigen(Symmetric(K))
	F = eigen(Symmetric(convert(Matrix,K))) # convert to matrix handles the case when K is sparse and is a no-op otherwise
	F.vectors*Diagonal(sqrt.(F.values))*F.vectors'
end




function _pma(X::AbstractMatrix, S::AbstractMatrix; nsv=typemax(Int))
	P,N = size(X)
	Y = X*S

	K = Symmetric(Y'Y)
	nsv = min(nsv, N)
	F = eigen(K, N-nsv+1:N)
	Σ = sqrt.(max.(0.,reverse(F.values)))
	VV = F.vectors[:,end:-1:1] # Coordinates of simplex equivalents, not interesting in practice.

	U = Y*VV ./ Σ'
	V = X'U ./ Σ' # Coordinates of original sample points in low dimensional space.

	PMA(U,Σ,V,VV)
end

pma(X::AbstractMatrix, G::AbstractMatrix{Bool}; kwargs...) = _pma(X, simplexgraph2kernelmatrixroot(G); kwargs...)




function _pma2(X::AbstractMatrix, dims::Integer, variableKernel, sampleKernelRoot)
	P,N = size(X)

	R2 = Symmetric(X'*(variableKernel*X))
	K = Symmetric(sampleKernelRoot'R2*sampleKernelRoot)
	dims = min(dims, N)
	F = eigen(K, N-dims+1:N)
	Σ = sqrt.(max.(0.,reverse(F.values)))
	W = F.vectors[:,end:-1:1] # Basis for sample simplices

	U = X*(sampleKernelRoot*W) ./ Σ' # Coordinates of original variables in low dimensional space.

	# R = cholesky(R2).U # doesn't work well since the matrix is only positively semidefinite and round errors will cause failure
	FR = eigen(R2)
	R = (FR.vectors * Diagonal(sqrt.(max.(0.,FR.values))))'

	Z = R*(sampleKernelRoot*W) ./ Σ' # Z columns are *a* base for the variable simplices, but is not interesting in itself.
	V = R'Z ./ Σ'

	U,Σ,V,Z,W
end


function pma2(X; dims=typemax(Int),
	             variableGraph=nothing, variableKernel=nothing,
	             sampleGraph=nothing, sampleKernelRoot=nothing)
	@assert any(!isequal(nothing),(variableGraph,variableKernel,sampleGraph,sampleKernelRoot)) "Specify at least one of the keyword arguments variableGraph and sampleGraph"

	@assert !(variableGraph!=nothing && variableKernel!=nothing) "Only one of the keyword arguments variableGraph and variableKernel can be used."
	@assert !(sampleGraph!=nothing && sampleKernelRoot!=nothing) "Only one of the keyword arguments sampleGraph and sampleKernelRoot can be used."

	variableGraph  != nothing && (variableKernel = simplexgraph2kernelmatrix(variableGraph))
	variableKernel == nothing && (variableKernel = I)

	sampleGraph      != nothing && (sampleKernelRoot = simplexgraph2kernelmatrixroot(sampleGraph))
	sampleKernelRoot == nothing && (sampleKernelRoot = I)

	_pma2(X, dims, variableKernel, sampleKernelRoot)
end


# function _pma2(X::AbstractMatrix, SVariables::AbstractMatrix, SSamples::AbstractMatrix; dim=typemax(Int))
# 	P,N = size(X)
# 	Y = SVariables*X*SSamples

# 	K = Symmetric(Y'Y)
# 	dim = min(dim, N)
# 	F = eigen(K, N-dim+1:N)
# 	Σ = sqrt.(max.(0.,reverse(F.values)))
# 	VV = F.vectors[:,end:-1:1] # Coordinates of sample   side simplex equivalents
# 	UU = Y*VV ./ Σ'            # Coordinates of variable side simplex equivalents


# 	U = X*SSamples*VV   ./ Σ' # Coordinates of original variables     in low dimensional space.
# 	# V = X'SVariables'UU ./ Σ' # Coordinates of original sample points in low dimensional space.
# 	V = X'*(SVariables'UU) ./ Σ' # Coordinates of original sample points in low dimensional space.

# 	# NxP PxP Pxd

# 	# U = Y*VV ./ Σ'
# 	# V = X'U ./ Σ' # Coordinates of original sample points in low dimensional space.

# 	U,Σ,V,UU,VV
# end

# pma2(X::AbstractMatrix, GVariables::AbstractMatrix, GSamples::AbstractMatrix; kwargs...) = _pma2(X, sparsegraph2simplices(GVariables), graph2simplices(GSamples); kwargs...)
