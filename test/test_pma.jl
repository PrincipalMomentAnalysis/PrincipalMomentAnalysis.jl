@testset "PMA" begin

function factorizationcmp(F1,F2)
	@test size(F1.U)==size(F2.U)
	@test size(F1.S)==size(F2.S)
	@test size(F1.V)==size(F2.V)
	@test size(F1.Vt)==size(F2.Vt)
	@test F1.S ≈ F2.S
	sgn = sign.(diag(F1.U'F2.U))
	@test F1.U.*sgn' ≈ F2.U
	@test F1.V.*sgn' ≈ F2.V
	@test F1.Vt.*sgn ≈ F2.Vt
	nothing
end

extractdims(F::SVD,dims) = SVD(F.U[:,dims], F.S[dims], F.Vt[dims,:])

A = [-8 0 6 0 0 0; -1 9 6 -2 0 -1; -1 0 -2 0 0 0; 0 1 0 1 0 -2; -9 -3 0 0 0 -5; -6 3 -8 -2 -4 8; 0 -1 -1 0 -3 -1; 4 0 0 -6 0 2]

@testset "PCA" begin
	FAns = svd(A)
	FPMA = pma(A,I(6); nsv=6)
	factorizationcmp(FAns,FPMA)
	FPMA4 = pma(A,I(6); nsv=4)
	factorizationcmp(extractdims(FAns,1:4),FPMA4)
end

@testset "TotalSimplexMass" begin
	A2 = hcat(repeat(A[:,1],1,5), repeat(A[:,2],1,2), A[:,3:end]) # five copies of first sample and two copies of second

	G = Matrix(I(11))
	G[1:5,1:5] .= true # Connect (identical) samples 1-5
	G[6:7,6:7] .= true # Connect (identical) samples 6-7

	FAns = svd(A2)
	# FPMA = pma(A2,G; nsv=11) # TODO: Add this test case if we add support for computing U,V columns corresponding to zero-valued singular values.
	# factorizationcmp(FAns,FPMA)
	FPMA6 = pma(A2,G; nsv=6) # non-zero singular values
	factorizationcmp(extractdims(FAns,1:6),FPMA6)
	FPMA4 = pma(A2,G; nsv=4)
	factorizationcmp(extractdims(FAns,1:4),FPMA4)
end

@testset "TotalSimplexMass2" begin
	A2 = hcat(repeat(A[:,1],1,4), A[:,2:end]) # four copies of first sample

	G = Matrix(I(9))
	G[1:4,1:4] .= LowerTriangular(ones(4,4)) # Simplices of dimensions 1-4 connecting (some of) the 4 identical samples

	FAns = svd(A2)
	# FPMA = pma(A2,G; nsv=9) # TODO: Add this test case if we add support for computing U,V columns corresponding to zero-valued singular values.
	# factorizationcmp(FAns,FPMA)
	FPMA6 = pma(A2,G; nsv=6) # non-zero singular values
	factorizationcmp(extractdims(FAns,1:6),FPMA6)
	FPMA4 = pma(A2,G; nsv=4)
	factorizationcmp(extractdims(FAns,1:4),FPMA4)
end


@testset "DisjointCliques" begin
	G = Bool[1 0 0 0 0 0; 0 1 1 0 0 0; 0 1 1 0 0 0; 0 0 0 1 1 1; 0 0 0 1 1 1; 0 0 0 1 1 1]
	# S = [1 0 0 0 0 0; 0 2/3 1/3 0 0 0; 0 1/3 2/3 0 0 0; 0 0 0 1/2 1/4 1/4; 0 0 0 1/4 1/2 1/4; 0 0 0 1/4 1/4 1/2] # corresponding S matrix
	a,b,c,d = 1/2+sqrt(1/12), 1/2-sqrt(1/12), 2/3, 1/6
	Ssq = [1 0 0 0 0 0; 0 a b 0 0 0; 0 b a 0 0 0; 0 0 0 c d d; 0 0 0 d c d; 0 0 0 d d c ]
	F = svd(A*Ssq)
	FVt = (F.U'A)./F.S
	FAns = SVD(F.U, F.S, FVt)

	FPMA = pma(A,G; nsv=6)
	factorizationcmp(FAns,FPMA)

	FPMA4 = pma(A,G; nsv=4)
	factorizationcmp(extractdims(FAns,1:4),FPMA4)
end

@testset "SymmetricGraph" begin
	G = Bool[1 0 0 0 1 0; 0 1 1 0 1 1; 0 1 1 1 1 0; 0 0 1 1 1 0; 1 1 1 1 1 1; 0 1 0 0 1 1]

	S = [8/21 1/42 1/42 1/42 4/21 1/42; 1/42 29/70 13/105 31/420 29/140 11/70; 1/42 13/105 29/70 11/70 29/140 31/420; 1/42 31/420 11/70 11/35 11/70 1/42; 4/21 29/140 29/140 11/70 32/35 11/70; 1/42 11/70 31/420 1/42 11/70 11/35]
	Ssq = sqrt(S);
	F = svd(A*Ssq)
	FVt = (F.U'A)./F.S
	FAns = SVD(F.U, F.S, FVt)

	FPMA = pma(A,G; nsv=6)
	factorizationcmp(FAns,FPMA)

	FPMA4 = pma(A,G; nsv=4)
	factorizationcmp(extractdims(FAns,1:4),FPMA4)
end

@testset "AsymmetricGraph" begin
	G = Bool[1 0 0 0 0 0; 0 1 1 1 1 0; 0 1 1 0 1 1; 0 1 1 1 0 1; 0 0 0 0 1 1; 1 0 0 0 0 1]

	S = [1/3 0 0 0 0 1/6; 0 5/6 1/4 1/3 1/12 0; 0 1/4 3/5 13/60 2/15 1/20; 0 1/3 13/60 23/30 1/20 1/20; 0 1/12 2/15 1/20 4/15 1/20; 1/6 0 1/20 1/20 1/20 13/30]
	Ssq = sqrt(S);
	F = svd(A*Ssq)
	FVt = (F.U'A)./F.S
	FAns = SVD(F.U, F.S, FVt)

	FPMA = pma(A,G; nsv=6)
	factorizationcmp(FAns,FPMA)

	FPMA4 = pma(A,G; nsv=4)
	factorizationcmp(extractdims(FAns,1:4),FPMA4)
end

@testset "AsymmetricGraph2" begin
	G = Bool[1 0 0 0 0 0; 0 0 1 1 1 0; 0 1 1 0 1 1; 0 1 1 0 0 1; 0 0 0 0 1 1; 1 0 0 0 0 0]

	S = [1/3 0 0 0 0 1/6; 0 4/3 1/6 1/12 1/12 0; 0 1/6 5/6 1/3 1/6 0; 0 1/12 1/3 2/3 1/12 0; 0 1/12 1/6 1/12 1/3 0; 1/6 0 0 0 0 1/3]
	Ssq = sqrt(S);
	F = svd(A*Ssq)
	FVt = (F.U'A)./F.S
	FAns = SVD(F.U, F.S, FVt)

	FPMA = pma(A,G; nsv=6)
	factorizationcmp(FAns,FPMA)

	FPMA4 = pma(A,G; nsv=4)
	factorizationcmp(extractdims(FAns,1:4),FPMA4)
end



end
