"""
	PMA <: Factorization

The output of `pma` representing the Principal Moment Analysis factorization of a matrix `A`.
If `F::PMA` is the factorization object, `U`, `S`, `V` and `Vt` can be obtained via `F.U`, `F.S`, `F.V` and `F.Vt` such as A ≈ U * Diagonal(S) * Vt'.

See also `pma`.
"""
struct PMA{T} <: Factorization{T}
	U::Matrix{T}
	S::Vector{T}
	Vt::Matrix{T}
	VV::Matrix{T}
end

function Base.getproperty(F::PMA,name::Symbol)
	name==:V && return F.Vt'
	getfield(F,name)
end

function simplices2kernelmatrix(sg::SimplexGraph)
	A = simplices2kernelmatrixroot(sg, simplify=false)
	A*A'
end

function simplices2kernelmatrixroot(sg::SimplexGraph; simplify=true)::SparseMatrixCSC
	s = vec(sum(sg.G,dims=1))

	L = sg.G * Diagonal(sqrt.(sg.w./(s.*(s.+1))))
	D = Diagonal(sqrt.(vec(sum(x->x^2,L,dims=2))))
	simplify || return [L D]

	QR = qr(vcat(sparse(L'), D))
	sparse(QR.R[:,invperm(QR.pcol)]')
end




function _pma(A::AbstractMatrix, S::AbstractMatrix; nsv::Integer=6)
	P,N = size(A)
	L = size(S,2)
	Y = A*S

	K = Symmetric(Y'Y)
	nsv = min(nsv, P, N)
	F = eigen(K, L-nsv+1:L)
	Σ = sqrt.(max.(0.,reverse(F.values)))
	VV = F.vectors[:,end:-1:1] # Coordinates of simplex equivalents, not interesting in practice.

	U = Y*VV ./ Σ'
	Vt = U'A ./ Σ # Coordinates of original sample points in low dimensional space.

	PMA(U,Σ,Vt,VV)
end

"""
	pma(A, G::SimplexGraph; nsv=6)

Computes the Principal Moment Analysis of the matrix `A` (variables × samples) using the sample SimplexGraph `G`.
Each column in `G` is a boolean vector representing one simplex. True means that a vertex is part of the simplex and false that it is not.
Set `nsv` to control the number of singular values and vectors returned.
Returns a `PMA` struct.

See also `PMA`, `SimplexGraph`.
"""
pma(A::AbstractMatrix, sg::SimplexGraph; kwargs...) = _pma(A, simplices2kernelmatrixroot(sg); kwargs...)


"""
	pma(A, G::AbstractMatrix{Bool}; nsv=6)

Computes the Principal Moment Analysis of the matrix `A` (variables × samples).
Each column in `G` is a boolean vector representing one simplex. True means that a vertex is part of the simplex and false that it is not.
Set `nsv` to control the number of singular values and vectors returned.
Returns a `PMA` struct.

See also `PMA`, `SimplexGraph`.
"""
pma(A::AbstractMatrix, G::AbstractMatrix{Bool}; kwargs...) = pma(A, SimplexGraph(G); kwargs...)
