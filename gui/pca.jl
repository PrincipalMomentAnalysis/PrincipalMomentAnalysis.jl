function pca(X::AbstractMatrix; dim=typemax(Int))
	P,N = size(X)

	K = Symmetric(X'X)
	dim = min(dim, N)
	F = eigen(K, N-dim+1:N)
	Σ = sqrt.(max.(0.,reverse(F.values)))
	V = F.vectors[:,end:-1:1]

	U = X*V ./ Σ'
	U,Σ,V
end
