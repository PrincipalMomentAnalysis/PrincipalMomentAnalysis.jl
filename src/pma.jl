simplexkernel(n::Integer, w::Real=1.0) = w^2/(n*(n+1))*(ones(n,n)+I)


# TODO: improve implementation! Might be possible to avoid the matrix factorization step.
function graph2simplices(G)
	N = size(G,1)
	# convert graph to simplex matrix
	S2 = zeros(N,N) # matrix describing all simplexes
	for j=1:size(G,2)
		ind = findall(G[:,j])
		S2[ind,ind] .+= simplexkernel(length(ind))
	end

	# Factor S^2 into SᵀS where S is symmetric too.
	F = eigen(Symmetric(S2))
	F.vectors*Diagonal(sqrt.(F.values))*F.vectors'
end




function _pma(X::AbstractMatrix, S::AbstractMatrix; dim=typemax(Int))
	P,N = size(X)
	Y = X*S

	K = Symmetric(Y'Y)
	dim = min(dim, N)
	F = eigen(K, N-dim+1:N)
	Σ = sqrt.(max.(0.,reverse(F.values)))
	VV = F.vectors[:,end:-1:1] # Coordinates of simplex equivalents, not interesting in practice.

	U = Y*VV ./ Σ'
	V = X'U ./ Σ' # Coordinates of original sample points in low dimensional space.

	U,Σ,V,VV
end

pma(X::AbstractMatrix, G::AbstractMatrix; kwargs...) = _pma(X, graph2simplices(G); kwargs...)





# # --- test with graph for variable side too ---
# function sparsegraph2simplices(G)
# 	N = size(G,1)
# 	# convert graph to simplex matrix
# 	S = spzeros(N,N) # matrix describing all simplexes
# 	for j=1:N
# 		ind = findall(G[:,j])
# 		S[ind,ind] .+= simplexmatrix(length(ind))
# 	end
# 	S
# end


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
