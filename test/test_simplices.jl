@testset "Simplices" begin

@testset "ThreeGroups" begin
	GAnswer = Bool[1 1 0 0; 1 1 0 0; 0 0 1 0; 0 0 0 1]
	SAnswer = [2/3 1/3 0 0; 1/3 2/3 0 0; 0 0 1 0; 0 0 0 1]

	annot = ["A","A","B","C"]
	G = groupsimplices(annot)
	@test G==GAnswer
	S = simplices2kernelmatrix(G)
	@test S ≈ SAnswer

	# try permuted version
	perm = [2, 3, 4, 1]
	G = groupsimplices(annot[perm])
	@test G==GAnswer[perm,perm]
	S = simplices2kernelmatrix(G)
	@test S ≈ SAnswer[perm,perm]
end

@testset "OneSample" begin
	annot = ["A"]
	G = groupsimplices(annot)
	@test G==trues(1,1)
	S = simplices2kernelmatrix(G)
	@test S==ones(1,1)
end

@testset "OneGroup" begin
	annot = ["A", "A", "A"]
	G = groupsimplices(annot)
	@test G==trues(3,3)
	S = simplices2kernelmatrix(G)
	@test S ≈ (Diagonal(ones(3)) + ones(3,3))/4
end

@testset "UniqueGroups" begin
	annot = ["A", "B", "C", "D"]
	G = groupsimplices(annot)
	@test G==Diagonal(trues(4))
	S = simplices2kernelmatrix(G)
	@test S == Diagonal(ones(4))
end



@testset "OneTimeSeries" begin
	N = 6
	annot = repeat(["A"], N)
	t = 1:N
	G = timeseriessimplices(t, groupby=annot)
	@test G == SymTridiagonal(trues(N),trues(N-1))
	S = simplices2kernelmatrix(G)
	@test S ≈ [1/2 1/4 1/12 0 0 0; 1/4 2/3 1/6 1/12 0 0; 1/12 1/6 1/2 1/6 1/12 0; 0 1/12 1/6 1/2 1/6 1/12; 0 0 1/12 1/6 2/3 1/4; 0 0 0 1/12 1/4 1/2]
end

@testset "TwoTimeSeries" begin
	N = 6
	annot = repeat(["A","B"], inner=N)
	t = repeat(1:N, 2)
	G = timeseriessimplices(t, groupby=annot)
	T = SymTridiagonal(trues(N),trues(N-1))
	Z = falses(N,N)
	@test G == [T Z; Z T]
	S = simplices2kernelmatrix(G)
	S1 = [1/2 1/4 1/12 0 0 0; 1/4 2/3 1/6 1/12 0 0; 1/12 1/6 1/2 1/6 1/12 0; 0 1/12 1/6 1/2 1/6 1/12; 0 0 1/12 1/6 2/3 1/4; 0 0 0 1/12 1/4 1/2]
	@test S ≈ [S1 Z; Z S1]
end

@testset "OneTimeSeries2" begin
	N = 6
	annot = repeat(["A"], 2N)
	t = div.(0:2N-1,2)
	G = timeseriessimplices(t, groupby=annot)
	@test G == repeat(SymTridiagonal(trues(N),trues(N-1)), inner=(2,2))
	S = simplices2kernelmatrix(G)
	@test S ≈ [31/105 31/210 31/210 31/210 1/21 1/21 0 0 0 0 0 0; 31/210 31/105 31/210 31/210 1/21 1/21 0 0 0 0 0 0; 31/210 31/210 41/105 41/210 2/21 2/21 1/21 1/21 0 0 0 0; 31/210 31/210 41/210 41/105 2/21 2/21 1/21 1/21 0 0 0 0; 1/21 1/21 2/21 2/21 2/7 1/7 2/21 2/21 1/21 1/21 0 0; 1/21 1/21 2/21 2/21 1/7 2/7 2/21 2/21 1/21 1/21 0 0; 0 0 1/21 1/21 2/21 2/21 2/7 1/7 2/21 2/21 1/21 1/21; 0 0 1/21 1/21 2/21 2/21 1/7 2/7 2/21 2/21 1/21 1/21; 0 0 0 0 1/21 1/21 2/21 2/21 41/105 41/210 31/210 31/210; 0 0 0 0 1/21 1/21 2/21 2/21 41/210 41/105 31/210 31/210; 0 0 0 0 0 0 1/21 1/21 31/210 31/210 31/105 31/210; 0 0 0 0 0 0 1/21 1/21 31/210 31/210 31/210 31/105]
end

@testset "MixedTimeSeries" begin
	GAnswer = Bool[1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0; 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0; 0 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0; 0 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0; 0 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1; 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1; 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1]
	SAnswer = [11/42 11/84 11/84 1/21 1/21 1/21 0 0 0 0 0 0 0 0 0 0 0 0; 11/84 97/210 97/420 31/210 31/210 31/210 0 0 0 0 0 0 0 0 0 0 0 0; 11/84 97/420 97/210 31/210 31/210 31/210 0 0 0 0 0 0 0 0 0 0 0 0; 1/21 31/210 31/210 31/105 31/210 31/210 0 0 0 0 0 0 0 0 0 0 0 0; 1/21 31/210 31/210 31/210 31/105 31/210 0 0 0 0 0 0 0 0 0 0 0 0; 1/21 31/210 31/210 31/210 31/210 31/105 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 1/1 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 1/2 1/4 1/4 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 1/4 1/2 1/4 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 1/4 1/4 1/2 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 13/30 13/60 13/60 1/20 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 13/60 13/30 13/60 1/20 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 13/60 13/60 1/2 1/12 1/30 1/30 1/30 0; 0 0 0 0 0 0 0 0 0 0 1/20 1/20 1/12 11/30 2/15 2/15 2/15 1/10; 0 0 0 0 0 0 0 0 0 0 0 0 1/30 2/15 11/30 11/60 11/60 3/20; 0 0 0 0 0 0 0 0 0 0 0 0 1/30 2/15 11/60 11/30 11/60 3/20; 0 0 0 0 0 0 0 0 0 0 0 0 1/30 2/15 11/60 11/60 11/30 3/20; 0 0 0 0 0 0 0 0 0 0 0 0 0 1/10 3/20 3/20 3/20 3/10]

	annot = ["A", "A", "A", "A", "A", "A", "B", "C", "C", "C", "D", "D", "D", "D", "D", "D", "D", "D" ]
	t     = [ 0,   5,   5,   7,   7,   7,   8,   4,   4,   8,   1,   1,   2,   3,   4,   4,   4,   5  ]
	G = timeseriessimplices(t, groupby=annot)
	@test G == GAnswer
	S = simplices2kernelmatrix(G)
	@test S ≈ SAnswer

	# try permuted version
	perm = [3, 9, 13, 12, 2, 4, 17, 10, 7, 14, 11, 1, 6, 8, 5, 15, 18, 16]
	G = timeseriessimplices(t[perm], groupby=annot[perm])
	@test G==GAnswer[perm,perm]
	S = simplices2kernelmatrix(G)
	@test S ≈ SAnswer[perm,perm]
end


@testset "kNN" begin
	A = Float64[9 1 3 9 3 3 7 3 5 5; 0 2 4 7 1 1 8 4 7 1; 0 6 3 7 0 9 0 2 0 6; 8 8 7 6 5 9 1 4 2 5];

	G = neighborsimplices(A)
	@test G == Diagonal(trues(10))
	S = simplices2kernelmatrix(G)
	@test S == Diagonal(ones(10))

	G = neighborsimplices(A; k=1)
	@test G == Bool[1 0 0 0 0 0 0 0 0 0; 0 1 0 0 0 1 0 0 0 1; 0 0 1 0 0 0 0 1 0 0; 0 0 0 1 0 0 0 0 0 0; 1 0 0 0 1 0 0 0 0 0; 0 1 0 0 0 1 0 0 0 0; 0 0 0 0 0 0 1 0 1 0; 0 0 1 0 1 0 0 1 0 0; 0 0 0 0 0 0 1 0 1 0; 0 0 0 1 0 0 0 0 0 1]
	S = simplices2kernelmatrix(G)
	@test S ≈ [1/3 0 0 0 1/6 0 0 0 0 0; 0 1/1 0 0 0 1/3 0 0 0 1/6; 0 0 2/3 0 0 0 0 1/3 0 0; 0 0 0 1/3 0 0 0 0 0 1/6; 1/6 0 0 0 2/3 0 0 1/6 0 0; 0 1/3 0 0 0 2/3 0 0 0 0; 0 0 0 0 0 0 2/3 0 1/3 0; 0 0 1/3 0 1/6 0 0 1/1 0 0; 0 0 0 0 0 0 1/3 0 2/3 0; 0 1/6 0 1/6 0 0 0 0 0 2/3]

	G = neighborsimplices(A; k=2)
	@test G == Bool[1 0 0 0 0 0 0 0 0 0; 0 1 1 0 0 1 0 0 0 1; 1 1 1 1 1 0 0 1 0 1; 0 0 0 1 0 0 0 0 0 0; 1 0 0 0 1 0 0 1 0 0; 0 1 0 0 0 1 0 0 0 0; 0 0 0 0 0 0 1 0 1 0; 0 0 1 0 1 0 1 1 1 0; 0 0 0 0 0 0 1 0 1 0; 0 0 0 1 0 1 0 0 0 1]
	S = simplices2kernelmatrix(G)
	@test S ≈ [1/6 0 1/12 0 1/12 0 0 0 0 0; 0 2/3 1/4 0 0 1/6 0 1/12 0 1/6; 1/12 1/4 7/6 1/12 1/4 1/12 0 1/4 0 1/6; 0 0 1/12 1/6 0 0 0 0 0 1/12; 1/12 0 1/4 0 1/2 0 0 1/6 0 0; 0 1/6 1/12 0 0 1/3 0 0 0 1/12; 0 0 0 0 0 0 1/3 1/6 1/6 0; 0 1/12 1/4 0 1/6 0 1/6 5/6 1/6 0; 0 0 0 0 0 0 1/6 1/6 1/3 0; 0 1/6 1/6 1/12 0 1/12 0 0 0 1/2]

	G = neighborsimplices(A; k=6)
	@test G == Bool[1 0 0 0 1 0 0 0 0 0; 0 1 1 0 1 1 0 1 0 1; 1 1 1 1 1 1 1 1 1 1; 1 1 0 1 0 1 1 0 1 1; 1 1 1 0 1 1 1 1 1 1; 0 1 1 1 0 1 0 0 0 1; 0 0 0 1 0 0 1 1 1 0; 1 1 1 1 1 1 1 1 1 1; 1 0 1 1 1 0 1 1 1 0; 1 1 1 1 1 1 1 1 1 1]
	S = simplices2kernelmatrix(G)
	@test S ≈ [1/14 1/56 1/28 1/56 1/28 0 0 1/28 1/28 1/28; 1/56 3/14 3/28 3/56 3/28 1/14 1/56 3/28 3/56 3/28; 1/28 3/28 5/14 1/8 9/56 5/56 1/14 5/28 1/8 5/28; 1/56 3/56 1/8 1/4 3/28 1/14 3/56 1/8 1/14 1/8 ; 1/28 3/28 9/56 3/28 9/28 1/14 3/56 9/56 3/28 9/56; 0 1/14 5/56 1/14 1/14 5/28 1/56 5/56 1/28 5/56; 0 1/56 1/14 3/56 3/56 1/56 1/7 1/14 1/14 1/14; 1/28 3/28 5/28 1/8 9/56 5/56 1/14 5/14 1/8 5/28; 1/28 3/56 1/8 1/14 3/28 1/28 1/14 1/8 1/4 1/8 ; 1/28 3/28 5/28 1/8 9/56 5/56 1/14 5/28 1/8 5/14]

	G = neighborsimplices(A; k=9)
	@test G == trues(10,10)
	S = simplices2kernelmatrix(G)
	@test S ≈ (ones(10,10)+Diagonal(ones(10)))/11
end
@testset "kNNSparse" begin
	A = Float64[9 1 3 9 3 3 7 3 5 5; 0 2 4 7 1 1 8 4 7 1; 0 6 3 7 0 9 0 2 0 6; 8 8 7 6 5 9 1 4 2 5];

	G = sparseneighborsimplices(A)
	@test G == Diagonal(trues(10))
	S = simplices2kernelmatrix(G)
	@test S == Diagonal(ones(10))

	G = sparseneighborsimplices(A; k=1)
	@test G == Bool[1 0 0 0 0 0 0 0 0 0; 0 1 0 0 0 1 0 0 0 1; 0 0 1 0 0 0 0 1 0 0; 0 0 0 1 0 0 0 0 0 0; 1 0 0 0 1 0 0 0 0 0; 0 1 0 0 0 1 0 0 0 0; 0 0 0 0 0 0 1 0 1 0; 0 0 1 0 1 0 0 1 0 0; 0 0 0 0 0 0 1 0 1 0; 0 0 0 1 0 0 0 0 0 1]
	S = simplices2kernelmatrix(G)
	@test S ≈ [1/3 0 0 0 1/6 0 0 0 0 0; 0 1/1 0 0 0 1/3 0 0 0 1/6; 0 0 2/3 0 0 0 0 1/3 0 0; 0 0 0 1/3 0 0 0 0 0 1/6; 1/6 0 0 0 2/3 0 0 1/6 0 0; 0 1/3 0 0 0 2/3 0 0 0 0; 0 0 0 0 0 0 2/3 0 1/3 0; 0 0 1/3 0 1/6 0 0 1/1 0 0; 0 0 0 0 0 0 1/3 0 2/3 0; 0 1/6 0 1/6 0 0 0 0 0 2/3]

	G = sparseneighborsimplices(A; k=2)
	@test G == Bool[1 0 0 0 0 0 0 0 0 0; 0 1 1 0 0 1 0 0 0 1; 1 1 1 1 1 0 0 1 0 1; 0 0 0 1 0 0 0 0 0 0; 1 0 0 0 1 0 0 1 0 0; 0 1 0 0 0 1 0 0 0 0; 0 0 0 0 0 0 1 0 1 0; 0 0 1 0 1 0 1 1 1 0; 0 0 0 0 0 0 1 0 1 0; 0 0 0 1 0 1 0 0 0 1]
	S = simplices2kernelmatrix(G)
	@test S ≈ [1/6 0 1/12 0 1/12 0 0 0 0 0; 0 2/3 1/4 0 0 1/6 0 1/12 0 1/6; 1/12 1/4 7/6 1/12 1/4 1/12 0 1/4 0 1/6; 0 0 1/12 1/6 0 0 0 0 0 1/12; 1/12 0 1/4 0 1/2 0 0 1/6 0 0; 0 1/6 1/12 0 0 1/3 0 0 0 1/12; 0 0 0 0 0 0 1/3 1/6 1/6 0; 0 1/12 1/4 0 1/6 0 1/6 5/6 1/6 0; 0 0 0 0 0 0 1/6 1/6 1/3 0; 0 1/6 1/6 1/12 0 1/12 0 0 0 1/2]

	G = sparseneighborsimplices(A; k=6)
	@test G == Bool[1 0 0 0 1 0 0 0 0 0; 0 1 1 0 1 1 0 1 0 1; 1 1 1 1 1 1 1 1 1 1; 1 1 0 1 0 1 1 0 1 1; 1 1 1 0 1 1 1 1 1 1; 0 1 1 1 0 1 0 0 0 1; 0 0 0 1 0 0 1 1 1 0; 1 1 1 1 1 1 1 1 1 1; 1 0 1 1 1 0 1 1 1 0; 1 1 1 1 1 1 1 1 1 1]
	S = simplices2kernelmatrix(G)
	@test S ≈ [1/14 1/56 1/28 1/56 1/28 0 0 1/28 1/28 1/28; 1/56 3/14 3/28 3/56 3/28 1/14 1/56 3/28 3/56 3/28; 1/28 3/28 5/14 1/8 9/56 5/56 1/14 5/28 1/8 5/28; 1/56 3/56 1/8 1/4 3/28 1/14 3/56 1/8 1/14 1/8 ; 1/28 3/28 9/56 3/28 9/28 1/14 3/56 9/56 3/28 9/56; 0 1/14 5/56 1/14 1/14 5/28 1/56 5/56 1/28 5/56; 0 1/56 1/14 3/56 3/56 1/56 1/7 1/14 1/14 1/14; 1/28 3/28 5/28 1/8 9/56 5/56 1/14 5/14 1/8 5/28; 1/28 3/56 1/8 1/14 3/28 1/28 1/14 1/8 1/4 1/8 ; 1/28 3/28 5/28 1/8 9/56 5/56 1/14 5/28 1/8 5/14]

	G = sparseneighborsimplices(A; k=9)
	@test G == trues(10,10)
	S = simplices2kernelmatrix(G)
	@test S ≈ (ones(10,10)+Diagonal(ones(10)))/11
end


@testset "NNDist" begin
	A = Float64[9 1 0 4 0 9 0 8 7 2; 3 6 4 4 0 6 4 2 5 4; 0 7 7 7 8 7 6 5 0 2; 5 6 8 9 5 5 5 1 7 9]

	G = neighborsimplices(A; r=0.2)
	@test G == Diagonal(trues(10))
	S = simplices2kernelmatrix(G)
	@test S == Diagonal(ones(10))

	G = neighborsimplices(A; r=0.5)
	@test G == Bool[1 0 0 0 0 0 0 0 1 0; 0 1 1 1 0 0 1 0 0 0; 0 1 1 1 1 0 1 0 0 1; 0 1 1 1 0 0 1 0 0 1; 0 0 1 0 1 0 1 0 0 0; 0 0 0 0 0 1 0 1 0 0; 0 1 1 1 1 0 1 0 0 1; 0 0 0 0 0 1 0 1 0 0; 1 0 0 0 0 0 0 0 1 1; 0 0 1 1 0 0 1 0 1 1]
	S = simplices2kernelmatrix(G)
	@test S ≈ [1/2 0 0 0 0 0 0 0 1/4 1/12; 0 11/42 11/84 11/84 1/21 0 11/84 0 0 17/210; 0 11/84 52/105 23/140 11/84 0 26/105 0 1/30 4/35; 0 11/84 23/140 23/70 1/21 0 23/140 0 1/30 4/35; 0 1/21 11/84 1/21 11/42 0 11/84 0 0 1/21; 0 0 0 0 0 2/3 0 1/3 0 0; 0 11/84 26/105 23/140 11/84 0 52/105 0 1/30 4/35; 0 0 0 0 0 1/3 0 2/3 0 0; 1/4 0 1/30 1/30 0 0 1/30 0 17/30 7/60; 1/12 17/210 4/35 4/35 1/21 0 4/35 0 7/60 83/210]

	G = neighborsimplices(A; r=0.8)
	@test G == Bool[1 0 0 1 0 1 0 1 1 1; 0 1 1 1 1 1 1 1 1 1; 0 1 1 1 1 1 1 0 0 1; 1 1 1 1 1 1 1 1 1 1; 0 1 1 1 1 0 1 1 0 1; 1 1 1 1 0 1 1 1 1 1; 0 1 1 1 1 1 1 1 1 1; 1 1 0 1 1 1 1 1 1 0; 1 1 0 1 0 1 1 1 1 1; 1 1 1 1 1 1 1 0 1 1]
	S = simplices2kernelmatrix(G)
	@test S ≈ [383/2310 13/220 31/990 383/4620 3/88 383/4620 13/220 199/2772 383/4620 1913/27720; 13/220 811/3465 1237/13860 811/6930 2551/27720 2749/27720 811/6930 2441/27720 161/1980 953/9240; 31/990 1237/13860 1237/6930 1237/13860 361/4620 1979/27720 1237/13860 557/9240 53/990 1237/13860; 383/4620 811/6930 1237/13860 976/3465 2551/27720 487/3960 811/6930 443/3960 1457/13860 391/3080; 3/88 2551/27720 361/4620 2551/27720 2551/13860 257/3465 2551/27720 437/6930 223/3960 361/4620; 383/4620 2749/27720 1979/27720 487/3960 257/3465 487/1980 2749/27720 1303/13860 1457/13860 6/55; 13/220 811/6930 1237/13860 811/6930 2551/27720 2749/27720 811/3465 2441/27720 161/1980 953/9240; 199/2772 2441/27720 557/9240 443/3960 437/6930 1303/13860 2441/27720 443/1980 1303/13860 97/990; 383/4620 161/1980 53/990 1457/13860 223/3960 1457/13860 161/1980 1303/13860 1457/6930 281/3080; 1913/27720 953/9240 1237/13860 391/3080 361/4620 6/55 953/9240 97/990 281/3080 391/1540]

	G = neighborsimplices(A; r=1.0)
	@test G == trues(10,10)
	S = simplices2kernelmatrix(G)
	@test S ≈ (Diagonal(ones(10)) + ones(10,10))/11
end
@testset "NNDistSparse" begin
	A = Float64[9 1 0 4 0 9 0 8 7 2; 3 6 4 4 0 6 4 2 5 4; 0 7 7 7 8 7 6 5 0 2; 5 6 8 9 5 5 5 1 7 9]

	G = sparseneighborsimplices(A; r=0.2)
	@test G == Diagonal(trues(10))
	S = simplices2kernelmatrix(G)
	@test S == Diagonal(ones(10))

	G = sparseneighborsimplices(A; r=0.5)
	@test G == Bool[1 0 0 0 0 0 0 0 1 0; 0 1 1 1 0 0 1 0 0 0; 0 1 1 1 1 0 1 0 0 1; 0 1 1 1 0 0 1 0 0 1; 0 0 1 0 1 0 1 0 0 0; 0 0 0 0 0 1 0 1 0 0; 0 1 1 1 1 0 1 0 0 1; 0 0 0 0 0 1 0 1 0 0; 1 0 0 0 0 0 0 0 1 1; 0 0 1 1 0 0 1 0 1 1]
	S = simplices2kernelmatrix(G)
	@test S ≈ [1/2 0 0 0 0 0 0 0 1/4 1/12; 0 11/42 11/84 11/84 1/21 0 11/84 0 0 17/210; 0 11/84 52/105 23/140 11/84 0 26/105 0 1/30 4/35; 0 11/84 23/140 23/70 1/21 0 23/140 0 1/30 4/35; 0 1/21 11/84 1/21 11/42 0 11/84 0 0 1/21; 0 0 0 0 0 2/3 0 1/3 0 0; 0 11/84 26/105 23/140 11/84 0 52/105 0 1/30 4/35; 0 0 0 0 0 1/3 0 2/3 0 0; 1/4 0 1/30 1/30 0 0 1/30 0 17/30 7/60; 1/12 17/210 4/35 4/35 1/21 0 4/35 0 7/60 83/210]

	G = sparseneighborsimplices(A; r=0.8)
	@test G == Bool[1 0 0 1 0 1 0 1 1 1; 0 1 1 1 1 1 1 1 1 1; 0 1 1 1 1 1 1 0 0 1; 1 1 1 1 1 1 1 1 1 1; 0 1 1 1 1 0 1 1 0 1; 1 1 1 1 0 1 1 1 1 1; 0 1 1 1 1 1 1 1 1 1; 1 1 0 1 1 1 1 1 1 0; 1 1 0 1 0 1 1 1 1 1; 1 1 1 1 1 1 1 0 1 1]
	S = simplices2kernelmatrix(G)
	@test S ≈ [383/2310 13/220 31/990 383/4620 3/88 383/4620 13/220 199/2772 383/4620 1913/27720; 13/220 811/3465 1237/13860 811/6930 2551/27720 2749/27720 811/6930 2441/27720 161/1980 953/9240; 31/990 1237/13860 1237/6930 1237/13860 361/4620 1979/27720 1237/13860 557/9240 53/990 1237/13860; 383/4620 811/6930 1237/13860 976/3465 2551/27720 487/3960 811/6930 443/3960 1457/13860 391/3080; 3/88 2551/27720 361/4620 2551/27720 2551/13860 257/3465 2551/27720 437/6930 223/3960 361/4620; 383/4620 2749/27720 1979/27720 487/3960 257/3465 487/1980 2749/27720 1303/13860 1457/13860 6/55; 13/220 811/6930 1237/13860 811/6930 2551/27720 2749/27720 811/3465 2441/27720 161/1980 953/9240; 199/2772 2441/27720 557/9240 443/3960 437/6930 1303/13860 2441/27720 443/1980 1303/13860 97/990; 383/4620 161/1980 53/990 1457/13860 223/3960 1457/13860 161/1980 1303/13860 1457/6930 281/3080; 1913/27720 953/9240 1237/13860 391/3080 361/4620 6/55 953/9240 97/990 281/3080 391/1540]

	G = sparseneighborsimplices(A; r=1.0)
	@test G == trues(10,10)
	S = simplices2kernelmatrix(G)
	@test S ≈ (Diagonal(ones(10)) + ones(10,10))/11
end


@testset "kNNDist" begin
	A = Float64[1 0 -9 0 0 -6 3 -3 -1 7; -1 0 -5 1 1 0 -6 0 0 -1; 0 -5 2 2 0 -1 -2 0 0 4; -4 -2 0 -7 -5 -4 0 0 -5 0]

	G = neighborsimplices(A; k=2, r=0.35)
	@test G == Bool[1 1 0 1 1 0 1 1 1 1; 1 1 0 0 0 0 1 0 0 0; 0 0 1 0 0 0 0 0 0 0; 1 0 0 1 1 0 0 0 1 0; 1 1 0 1 1 0 0 0 1 0; 0 0 1 0 0 1 0 1 1 0; 0 0 0 0 0 0 1 0 0 1; 1 0 1 0 0 1 0 1 1 0; 1 0 0 1 1 1 0 1 1 0; 0 0 0 0 0 0 0 0 0 1]
	S = simplices2kernelmatrix(G)
	@test S ≈ [94/105 4/21 0 31/210 97/420 31/420 1/6 41/420 83/420 1/12; 4/21 8/21 0 1/42 3/28 0 1/12 1/42 1/42 0; 0 0 1/6 0 0 1/12 0 1/12 0 0; 31/210 1/42 0 31/105 31/210 1/42 0 1/21 31/210 0; 97/420 3/28 0 31/210 97/210 1/42 0 1/21 31/210 0; 31/420 0 1/12 1/42 1/42 101/210 0 101/420 11/70 0; 1/6 1/12 0 0 0 0 1/3 0 0 1/12; 41/420 1/42 1/12 1/21 1/21 101/420 0 37/70 19/105 0; 83/420 1/42 0 31/210 31/210 11/70 0 19/105 59/105 0; 1/12 0 0 0 0 0 1/12 0 0 1/6]
end
@testset "kNNDistSparse" begin
	A = Float64[1 0 -9 0 0 -6 3 -3 -1 7; -1 0 -5 1 1 0 -6 0 0 -1; 0 -5 2 2 0 -1 -2 0 0 4; -4 -2 0 -7 -5 -4 0 0 -5 0]

	G = sparseneighborsimplices(A; k=2, r=0.35)
	@test G == Bool[1 1 0 1 1 0 1 1 1 1; 1 1 0 0 0 0 1 0 0 0; 0 0 1 0 0 0 0 0 0 0; 1 0 0 1 1 0 0 0 1 0; 1 1 0 1 1 0 0 0 1 0; 0 0 1 0 0 1 0 1 1 0; 0 0 0 0 0 0 1 0 0 1; 1 0 1 0 0 1 0 1 1 0; 1 0 0 1 1 1 0 1 1 0; 0 0 0 0 0 0 0 0 0 1]
	S = simplices2kernelmatrix(G)
	@test S ≈ [94/105 4/21 0 31/210 97/420 31/420 1/6 41/420 83/420 1/12; 4/21 8/21 0 1/42 3/28 0 1/12 1/42 1/42 0; 0 0 1/6 0 0 1/12 0 1/12 0 0; 31/210 1/42 0 31/105 31/210 1/42 0 1/21 31/210 0; 97/420 3/28 0 31/210 97/210 1/42 0 1/21 31/210 0; 31/420 0 1/12 1/42 1/42 101/210 0 101/420 11/70 0; 1/6 1/12 0 0 0 0 1/3 0 0 1/12; 41/420 1/42 1/12 1/21 1/21 101/420 0 37/70 19/105 0; 83/420 1/42 0 31/210 31/210 11/70 0 19/105 59/105 0; 1/12 0 0 0 0 0 1/12 0 0 1/6]
end


@testset "kNNDistReduced" begin
    A = Float64[76 48 40 14 -15 7 -42 -57 -37 -34; 84 23 59 76 22 -8 -123 -44 -77 -12; -56 -28 36 -59 37 50 -21 8 23 10; -19 -59 -84 42 50 22 -5 38 -19 34; 33 56 -9 -97 56 2 -2 -44 -19 24; 36 98 -13 32 19 25 -81 53 -50 -119]

    G = neighborsimplices(A; k=2, r=0.35, dim=1)
    @test G==Bool[1 1 1 1 0 0 0 0 0 0; 1 1 1 1 0 0 0 0 0 0; 1 1 1 1 1 1 0 1 0 0; 1 1 1 1 1 1 0 0 0 0; 0 0 1 1 1 1 0 1 0 0; 0 0 1 1 1 1 0 1 0 0; 0 0 0 0 0 0 1 0 1 1; 0 0 1 0 1 1 0 1 1 1; 0 0 0 0 0 0 1 1 1 1; 0 0 0 0 0 0 1 1 1 1]
    S = simplices2kernelmatrix(G)
    @test S ≈ [17/60 17/120 17/120 17/120 1/24 1/24 0 1/56 0 0; 17/120 17/60 17/120 17/120 1/24 1/24 0 1/56 0 0; 17/120 17/120 13/28 5/24 37/280 37/280 0 13/120 1/42 1/42; 17/120 17/120 5/24 5/12 13/120 13/120 0 71/840 0 0; 1/24 1/24 37/280 13/120 37/140 37/280 0 13/120 1/42 1/42; 1/24 1/24 37/280 13/120 37/280 37/140 0 13/120 1/42 1/42; 0 0 0 0 0 0 11/30 1/10 11/60 11/60; 1/56 1/56 13/120 71/840 13/120 13/120 1/10 5/12 13/105 13/105; 0 0 1/42 0 1/42 1/42 11/60 13/105 29/70 29/140; 0 0 1/42 0 1/42 1/42 11/60 13/105 29/140 29/70]

    G = neighborsimplices(A; k=2, r=0.45, dim=3)
    @test G==Bool[1 1 1 1 0 0 0 0 0 0; 1 1 1 0 1 0 0 0 0 0; 1 1 1 0 1 1 0 0 0 0; 0 0 0 1 0 0 0 0 0 0; 0 1 1 0 1 1 0 1 1 0; 0 0 1 1 1 1 0 1 1 0; 0 0 0 0 0 0 1 0 1 1; 0 0 0 0 1 1 0 1 0 0; 0 0 0 0 1 1 1 0 1 1; 0 0 0 0 0 0 1 0 1 1]
    S = simplices2kernelmatrix(G)
    @test S ≈ [1/2 1/6 1/6 1/12 1/12 7/60 0 0 0 0; 1/6 8/21 4/21 0 3/28 2/35 0 1/42 1/42 0; 1/6 4/21 47/105 0 59/420 19/210 0 2/35 2/35 0; 1/12 0 0 1/6 0 1/12 0 0 0 0; 1/12 3/28 59/420 0 18/35 29/140 1/30 59/420 19/210 1/30; 7/60 2/35 19/210 1/12 29/140 61/105 1/30 59/420 19/210 1/30; 0 0 0 0 1/30 1/30 2/5 0 1/5 1/5; 0 1/42 2/35 0 59/420 59/420 0 59/210 2/35 0; 0 1/42 2/35 0 19/210 19/210 1/5 2/35 18/35 1/5; 0 0 0 0 1/30 1/30 1/5 0 1/5 2/5]

    G = neighborsimplices(A; k=2, r=0.55, dim=5)
    @test G==Bool[1 1 1 1 0 0 0 0 0 0; 1 1 0 0 0 0 0 0 0 0; 1 0 1 0 0 1 0 0 0 0; 0 0 0 1 0 0 0 0 0 0; 0 0 0 0 1 1 0 1 0 0; 0 1 1 0 1 1 0 1 1 0; 0 0 0 0 0 0 1 0 1 1; 0 0 0 1 1 1 0 1 1 0; 0 0 0 0 0 1 1 1 1 1; 0 0 0 0 0 0 1 0 1 1]
    S = simplices2kernelmatrix(G)
    @test S ≈ [2/3 1/6 1/6 1/12 0 1/6 0 1/12 0 0; 1/6 1/3 1/12 0 0 1/12 0 0 0 0; 1/6 1/12 2/5 0 1/30 7/60 0 1/30 1/30 0; 1/12 0 0 1/6 0 0 0 1/12 0 0; 0 0 1/30 0 1/3 1/6 0 1/6 1/12 0; 1/6 1/12 7/60 0 1/6 11/15 1/30 1/5 7/60 1/30; 0 0 0 0 0 1/30 2/5 1/30 1/5 1/5; 1/12 0 1/30 1/12 1/6 1/5 1/30 17/30 7/60 1/30; 0 0 1/30 0 1/12 7/60 1/5 7/60 17/30 1/5; 0 0 0 0 0 1/30 1/5 1/30 1/5 2/5]

    G = neighborsimplices(A'; k=1, r=0.65, dim=1)
    @test G==Bool[1 1 0 0 1 1; 1 1 0 0 0 1; 0 0 1 1 1 0; 0 0 1 1 1 0; 1 0 1 1 1 1; 1 1 0 0 1 1]
    S = simplices2kernelmatrix(G)
    @test S ≈ [13/30 11/60 1/30 1/30 2/15 13/60; 11/60 11/30 0 0 1/10 11/60; 1/30 0 2/5 1/5 1/5 1/30; 1/30 0 1/5 2/5 1/5 1/30; 2/15 1/10 1/5 1/5 3/5 2/15; 13/60 11/60 1/30 1/30 2/15 13/30]

    G = neighborsimplices(A'; k=1, r=0.75, dim=3)
    @test G==Bool[1 1 1 0 1 1; 1 1 0 0 0 0; 1 0 1 1 1 0; 0 0 1 1 0 0; 1 0 1 0 1 0; 1 0 0 0 0 1]
    S = simplices2kernelmatrix(G)
    @test S ≈ [1/1 1/5 1/6 1/20 1/6 1/5; 1/5 2/5 1/30 0 1/30 1/30; 1/6 1/30 2/3 13/60 1/6 1/30; 1/20 0 13/60 13/30 1/20 0; 1/6 1/30 1/6 1/20 1/3 1/30; 1/5 1/30 1/30 0 1/30 2/5]

    G = neighborsimplices(A'; k=1, r=0.85, dim=5)
    @test G==Bool[1 1 1 0 1 1; 1 1 0 0 0 1; 1 0 1 1 1 0; 0 0 1 1 1 0; 1 0 1 1 1 0; 1 1 0 0 0 1]
    S = simplices2kernelmatrix(G)
    @test S ≈ [3/5 1/5 2/15 1/10 2/15 1/5; 1/5 2/5 1/30 0 1/30 1/5; 2/15 1/30 13/30 11/60 13/60 1/30; 1/10 0 11/60 11/30 11/60 0; 2/15 1/30 13/60 11/60 13/30 1/30; 1/5 1/5 1/30 0 1/30 2/5]
end
@testset "kNNDistReducedSparse" begin
    A = Float64[76 48 40 14 -15 7 -42 -57 -37 -34; 84 23 59 76 22 -8 -123 -44 -77 -12; -56 -28 36 -59 37 50 -21 8 23 10; -19 -59 -84 42 50 22 -5 38 -19 34; 33 56 -9 -97 56 2 -2 -44 -19 24; 36 98 -13 32 19 25 -81 53 -50 -119]

    G = sparseneighborsimplices(A; k=2, r=0.35, dim=1)
    @test G==Bool[1 1 1 1 0 0 0 0 0 0; 1 1 1 1 0 0 0 0 0 0; 1 1 1 1 1 1 0 1 0 0; 1 1 1 1 1 1 0 0 0 0; 0 0 1 1 1 1 0 1 0 0; 0 0 1 1 1 1 0 1 0 0; 0 0 0 0 0 0 1 0 1 1; 0 0 1 0 1 1 0 1 1 1; 0 0 0 0 0 0 1 1 1 1; 0 0 0 0 0 0 1 1 1 1]
    S = simplices2kernelmatrix(G)
    @test S ≈ [17/60 17/120 17/120 17/120 1/24 1/24 0 1/56 0 0; 17/120 17/60 17/120 17/120 1/24 1/24 0 1/56 0 0; 17/120 17/120 13/28 5/24 37/280 37/280 0 13/120 1/42 1/42; 17/120 17/120 5/24 5/12 13/120 13/120 0 71/840 0 0; 1/24 1/24 37/280 13/120 37/140 37/280 0 13/120 1/42 1/42; 1/24 1/24 37/280 13/120 37/280 37/140 0 13/120 1/42 1/42; 0 0 0 0 0 0 11/30 1/10 11/60 11/60; 1/56 1/56 13/120 71/840 13/120 13/120 1/10 5/12 13/105 13/105; 0 0 1/42 0 1/42 1/42 11/60 13/105 29/70 29/140; 0 0 1/42 0 1/42 1/42 11/60 13/105 29/140 29/70]

    G = sparseneighborsimplices(A; k=2, r=0.45, dim=3)
    @test G==Bool[1 1 1 1 0 0 0 0 0 0; 1 1 1 0 1 0 0 0 0 0; 1 1 1 0 1 1 0 0 0 0; 0 0 0 1 0 0 0 0 0 0; 0 1 1 0 1 1 0 1 1 0; 0 0 1 1 1 1 0 1 1 0; 0 0 0 0 0 0 1 0 1 1; 0 0 0 0 1 1 0 1 0 0; 0 0 0 0 1 1 1 0 1 1; 0 0 0 0 0 0 1 0 1 1]
    S = simplices2kernelmatrix(G)
    @test S ≈ [1/2 1/6 1/6 1/12 1/12 7/60 0 0 0 0; 1/6 8/21 4/21 0 3/28 2/35 0 1/42 1/42 0; 1/6 4/21 47/105 0 59/420 19/210 0 2/35 2/35 0; 1/12 0 0 1/6 0 1/12 0 0 0 0; 1/12 3/28 59/420 0 18/35 29/140 1/30 59/420 19/210 1/30; 7/60 2/35 19/210 1/12 29/140 61/105 1/30 59/420 19/210 1/30; 0 0 0 0 1/30 1/30 2/5 0 1/5 1/5; 0 1/42 2/35 0 59/420 59/420 0 59/210 2/35 0; 0 1/42 2/35 0 19/210 19/210 1/5 2/35 18/35 1/5; 0 0 0 0 1/30 1/30 1/5 0 1/5 2/5]

    G = sparseneighborsimplices(A; k=2, r=0.55, dim=5)
    @test G==Bool[1 1 1 1 0 0 0 0 0 0; 1 1 0 0 0 0 0 0 0 0; 1 0 1 0 0 1 0 0 0 0; 0 0 0 1 0 0 0 0 0 0; 0 0 0 0 1 1 0 1 0 0; 0 1 1 0 1 1 0 1 1 0; 0 0 0 0 0 0 1 0 1 1; 0 0 0 1 1 1 0 1 1 0; 0 0 0 0 0 1 1 1 1 1; 0 0 0 0 0 0 1 0 1 1]
    S = simplices2kernelmatrix(G)
    @test S ≈ [2/3 1/6 1/6 1/12 0 1/6 0 1/12 0 0; 1/6 1/3 1/12 0 0 1/12 0 0 0 0; 1/6 1/12 2/5 0 1/30 7/60 0 1/30 1/30 0; 1/12 0 0 1/6 0 0 0 1/12 0 0; 0 0 1/30 0 1/3 1/6 0 1/6 1/12 0; 1/6 1/12 7/60 0 1/6 11/15 1/30 1/5 7/60 1/30; 0 0 0 0 0 1/30 2/5 1/30 1/5 1/5; 1/12 0 1/30 1/12 1/6 1/5 1/30 17/30 7/60 1/30; 0 0 1/30 0 1/12 7/60 1/5 7/60 17/30 1/5; 0 0 0 0 0 1/30 1/5 1/30 1/5 2/5]

    G = sparseneighborsimplices(A'; k=1, r=0.65, dim=1)
    @test G==Bool[1 1 0 0 1 1; 1 1 0 0 0 1; 0 0 1 1 1 0; 0 0 1 1 1 0; 1 0 1 1 1 1; 1 1 0 0 1 1]
    S = simplices2kernelmatrix(G)
    @test S ≈ [13/30 11/60 1/30 1/30 2/15 13/60; 11/60 11/30 0 0 1/10 11/60; 1/30 0 2/5 1/5 1/5 1/30; 1/30 0 1/5 2/5 1/5 1/30; 2/15 1/10 1/5 1/5 3/5 2/15; 13/60 11/60 1/30 1/30 2/15 13/30]

    G = sparseneighborsimplices(A'; k=1, r=0.75, dim=3)
    @test G==Bool[1 1 1 0 1 1; 1 1 0 0 0 0; 1 0 1 1 1 0; 0 0 1 1 0 0; 1 0 1 0 1 0; 1 0 0 0 0 1]
    S = simplices2kernelmatrix(G)
    @test S ≈ [1/1 1/5 1/6 1/20 1/6 1/5; 1/5 2/5 1/30 0 1/30 1/30; 1/6 1/30 2/3 13/60 1/6 1/30; 1/20 0 13/60 13/30 1/20 0; 1/6 1/30 1/6 1/20 1/3 1/30; 1/5 1/30 1/30 0 1/30 2/5]

    G = sparseneighborsimplices(A'; k=1, r=0.85, dim=5)
    @test G==Bool[1 1 1 0 1 1; 1 1 0 0 0 1; 1 0 1 1 1 0; 0 0 1 1 1 0; 1 0 1 1 1 0; 1 1 0 0 0 1]
    S = simplices2kernelmatrix(G)
    @test S ≈ [3/5 1/5 2/15 1/10 2/15 1/5; 1/5 2/5 1/30 0 1/30 1/5; 2/15 1/30 13/30 11/60 13/60 1/30; 1/10 0 11/60 11/30 11/60 0; 2/15 1/30 13/60 11/60 13/30 1/30; 1/5 1/5 1/30 0 1/30 2/5]
end

end
