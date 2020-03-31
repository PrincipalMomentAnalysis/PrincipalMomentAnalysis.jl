function svdbyeigen(A; nsv::Integer=3)
	P,N = size(A)
	K = Symmetric(N<=P ? A'A : A*A')
	M = size(K,1)
	F = eigen(K, M-nsv+1:M)
	S = sqrt.(max.(0.,reverse(F.values)))

	V = F.vectors[:,end:-1:1]
	N<=P ? SVD(A*V./S',S,V') : SVD(V,S,V'A./S)
end
