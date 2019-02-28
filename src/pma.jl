
function simplexmatrix(n::Integer, w::Real=1.0)
	@assert n>=1
	n==1 && return w*ones(1,1)

	u = 1/n + (n-1)/(n*√(n+1))
	v = (1-u)/(n-1)

	# reweight
	u = w*u/n
	v = w*v/n

	u*I + v*(ones(n,n)-I)
end



function graph2simplices(G)
	N = size(G,1)
	# convert graph to simplex matrix
	S = zeros(N,N) # matrix describing all simplexes
	for j=1:N
		ind = findall(G[:,j])
		S[ind,ind] .+= simplexmatrix(length(ind))
	end
	S
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

	U,Σ,V
end

pma(X::AbstractMatrix, G::AbstractMatrix; kwargs...) = _pma(X, graph2simplices(G); kwargs...)
