@testset "Graph" begin

@testset "ThreeGroups" begin
	GAnswer = Bool[1 1 0 0; 1 1 0 0; 0 0 1 0; 0 0 0 1]
	SAnswer = [2/3 1/3 0 0; 1/3 2/3 0 0; 0 0 1 0; 0 0 0 1]

	annot = ["A","A","B","C"]
	G = buildgraph(annot)
	@test G==GAnswer
	S = simplexgraph2kernelmatrix(G)
	@test S ≈ SAnswer

	# try permuted version
	perm = [2, 3, 4, 1]
	G = buildgraph(annot[perm])
	@test G==GAnswer[perm,perm]
	S = simplexgraph2kernelmatrix(G)
	@test S ≈ SAnswer[perm,perm]
end

@testset "OneSample" begin
	annot = ["A"]
	G = buildgraph(annot)
	@test G==trues(1,1)
	S = simplexgraph2kernelmatrix(G)
	@test S==ones(1,1)
end

@testset "OneGroup" begin
	annot = ["A", "A", "A"]
	G = buildgraph(annot)
	@test G==trues(3,3)
	S = simplexgraph2kernelmatrix(G)
	@test S ≈ (Diagonal(ones(3)) + ones(3,3))/4
end

@testset "UniqueGroups" begin
	annot = ["A", "B", "C", "D"]
	G = buildgraph(annot)
	@test G==Diagonal(trues(4))
	S = simplexgraph2kernelmatrix(G)
	@test S == Diagonal(ones(4))
end



@testset "OneTimeSeries" begin
	N = 6
	annot = repeat(["A"], N)
	t = 1:N
	G = buildgraph(annot, t)
	@test G == SymTridiagonal(trues(N),trues(N-1))
	S = simplexgraph2kernelmatrix(G)
	@test S ≈ [1/2 1/4 1/12 0 0 0; 1/4 2/3 1/6 1/12 0 0; 1/12 1/6 1/2 1/6 1/12 0; 0 1/12 1/6 1/2 1/6 1/12; 0 0 1/12 1/6 2/3 1/4; 0 0 0 1/12 1/4 1/2]
end

@testset "TwoTimeSeries" begin
	N = 6
	annot = repeat(["A","B"], inner=N)
	t = repeat(1:N, 2)
	G = buildgraph(annot, t)
	T = SymTridiagonal(trues(N),trues(N-1))
	Z = falses(N,N)
	@test G == [T Z; Z T]
	S = simplexgraph2kernelmatrix(G)
	S1 = [1/2 1/4 1/12 0 0 0; 1/4 2/3 1/6 1/12 0 0; 1/12 1/6 1/2 1/6 1/12 0; 0 1/12 1/6 1/2 1/6 1/12; 0 0 1/12 1/6 2/3 1/4; 0 0 0 1/12 1/4 1/2]
	@test S ≈ [S1 Z; Z S1]
end

@testset "OneTimeSeries2" begin
	N = 6
	annot = repeat(["A"], 2N)
	t = div.(0:2N-1,2)
	G = buildgraph(annot, t)
	@test G == repeat(SymTridiagonal(trues(N),trues(N-1)), inner=(2,2))
	S = simplexgraph2kernelmatrix(G)
	@test S ≈ [31/105 31/210 31/210 31/210 1/21 1/21 0 0 0 0 0 0; 31/210 31/105 31/210 31/210 1/21 1/21 0 0 0 0 0 0; 31/210 31/210 41/105 41/210 2/21 2/21 1/21 1/21 0 0 0 0; 31/210 31/210 41/210 41/105 2/21 2/21 1/21 1/21 0 0 0 0; 1/21 1/21 2/21 2/21 2/7 1/7 2/21 2/21 1/21 1/21 0 0; 1/21 1/21 2/21 2/21 1/7 2/7 2/21 2/21 1/21 1/21 0 0; 0 0 1/21 1/21 2/21 2/21 2/7 1/7 2/21 2/21 1/21 1/21; 0 0 1/21 1/21 2/21 2/21 1/7 2/7 2/21 2/21 1/21 1/21; 0 0 0 0 1/21 1/21 2/21 2/21 41/105 41/210 31/210 31/210; 0 0 0 0 1/21 1/21 2/21 2/21 41/210 41/105 31/210 31/210; 0 0 0 0 0 0 1/21 1/21 31/210 31/210 31/105 31/210; 0 0 0 0 0 0 1/21 1/21 31/210 31/210 31/210 31/105]
end

@testset "MixedTimeSeries" begin
	GAnswer = Bool[1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0; 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0; 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0; 0 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0; 0 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0; 0 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 0; 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1; 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1; 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1; 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1]
	SAnswer = [11/42 11/84 11/84 1/21 1/21 1/21 0 0 0 0 0 0 0 0 0 0 0 0; 11/84 97/210 97/420 31/210 31/210 31/210 0 0 0 0 0 0 0 0 0 0 0 0; 11/84 97/420 97/210 31/210 31/210 31/210 0 0 0 0 0 0 0 0 0 0 0 0; 1/21 31/210 31/210 31/105 31/210 31/210 0 0 0 0 0 0 0 0 0 0 0 0; 1/21 31/210 31/210 31/210 31/105 31/210 0 0 0 0 0 0 0 0 0 0 0 0; 1/21 31/210 31/210 31/210 31/210 31/105 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 1/1 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 1/2 1/4 1/4 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 1/4 1/2 1/4 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 1/4 1/4 1/2 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 13/30 13/60 13/60 1/20 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 13/60 13/30 13/60 1/20 0 0 0 0; 0 0 0 0 0 0 0 0 0 0 13/60 13/60 1/2 1/12 1/30 1/30 1/30 0; 0 0 0 0 0 0 0 0 0 0 1/20 1/20 1/12 11/30 2/15 2/15 2/15 1/10; 0 0 0 0 0 0 0 0 0 0 0 0 1/30 2/15 11/30 11/60 11/60 3/20; 0 0 0 0 0 0 0 0 0 0 0 0 1/30 2/15 11/60 11/30 11/60 3/20; 0 0 0 0 0 0 0 0 0 0 0 0 1/30 2/15 11/60 11/60 11/30 3/20; 0 0 0 0 0 0 0 0 0 0 0 0 0 1/10 3/20 3/20 3/20 3/10]

	annot = ["A", "A", "A", "A", "A", "A", "B", "C", "C", "C", "D", "D", "D", "D", "D", "D", "D", "D" ]
	t     = [ 0,   5,   5,   7,   7,   7,   8,   4,   4,   8,   1,   1,   2,   3,   4,   4,   4,   5  ]
	G = buildgraph(annot,t)
	@test G == GAnswer
	S = simplexgraph2kernelmatrix(G)
	@test S ≈ SAnswer

	# try permuted version
	perm = [3, 9, 13, 12, 2, 4, 17, 10, 7, 14, 11, 1, 6, 8, 5, 15, 18, 16]
	G = buildgraph(annot[perm], t[perm])
	@test G==GAnswer[perm,perm]
	S = simplexgraph2kernelmatrix(G)
	@test S ≈ SAnswer[perm,perm]
end


@testset "kNN" begin
	GAnswer1 = Bool[1 0 0 0 0 0 0 0 0 0; 0 1 0 0 0 1 0 0 0 1; 0 0 1 0 0 0 0 1 0 0; 0 0 0 1 0 0 0 0 0 0; 1 0 0 0 1 0 0 0 0 0; 0 1 0 0 0 1 0 0 0 0; 0 0 0 0 0 0 1 0 1 0; 0 0 1 0 1 0 0 1 0 0; 0 0 0 0 0 0 1 0 1 0; 0 0 0 1 0 0 0 0 0 1]
	SAnswer1 = [1/3 0 0 0 1/6 0 0 0 0 0; 0 1/1 0 0 0 1/3 0 0 0 1/6; 0 0 2/3 0 0 0 0 1/3 0 0; 0 0 0 1/3 0 0 0 0 0 1/6; 1/6 0 0 0 2/3 0 0 1/6 0 0; 0 1/3 0 0 0 2/3 0 0 0 0; 0 0 0 0 0 0 2/3 0 1/3 0; 0 0 1/3 0 1/6 0 0 1/1 0 0; 0 0 0 0 0 0 1/3 0 2/3 0; 0 1/6 0 1/6 0 0 0 0 0 2/3]

	GAnswer2 = Bool[1 0 0 0 0 0 0 0 0 0; 0 1 1 0 0 1 0 0 0 1; 1 1 1 1 1 0 0 1 0 1; 0 0 0 1 0 0 0 0 0 0; 1 0 0 0 1 0 0 1 0 0; 0 1 0 0 0 1 0 0 0 0; 0 0 0 0 0 0 1 0 1 0; 0 0 1 0 1 0 1 1 1 0; 0 0 0 0 0 0 1 0 1 0; 0 0 0 1 0 1 0 0 0 1]
	SAnswer2 = [1/6 0 1/12 0 1/12 0 0 0 0 0; 0 2/3 1/4 0 0 1/6 0 1/12 0 1/6; 1/12 1/4 7/6 1/12 1/4 1/12 0 1/4 0 1/6; 0 0 1/12 1/6 0 0 0 0 0 1/12; 1/12 0 1/4 0 1/2 0 0 1/6 0 0; 0 1/6 1/12 0 0 1/3 0 0 0 1/12; 0 0 0 0 0 0 1/3 1/6 1/6 0; 0 1/12 1/4 0 1/6 0 1/6 5/6 1/6 0; 0 0 0 0 0 0 1/6 1/6 1/3 0; 0 1/6 1/6 1/12 0 1/12 0 0 0 1/2]

	GAnswer6 = Bool[1 0 0 0 1 0 0 0 0 0; 0 1 1 0 1 1 0 1 0 1; 1 1 1 1 1 1 1 1 1 1; 1 1 0 1 0 1 1 0 1 1; 1 1 1 0 1 1 1 1 1 1; 0 1 1 1 0 1 0 0 0 1; 0 0 0 1 0 0 1 1 1 0; 1 1 1 1 1 1 1 1 1 1; 1 0 1 1 1 0 1 1 1 0; 1 1 1 1 1 1 1 1 1 1]
	SAnswer6 = [1/14 1/56 1/28 1/56 1/28 0 0 1/28 1/28 1/28; 1/56 3/14 3/28 3/56 3/28 1/14 1/56 3/28 3/56 3/28; 1/28 3/28 5/14 1/8 9/56 5/56 1/14 5/28 1/8 5/28; 1/56 3/56 1/8 1/4 3/28 1/14 3/56 1/8 1/14 1/8 ; 1/28 3/28 9/56 3/28 9/28 1/14 3/56 9/56 3/28 9/56; 0 1/14 5/56 1/14 1/14 5/28 1/56 5/56 1/28 5/56; 0 1/56 1/14 3/56 3/56 1/56 1/7 1/14 1/14 1/14; 1/28 3/28 5/28 1/8 9/56 5/56 1/14 5/14 1/8 5/28; 1/28 3/56 1/8 1/14 3/28 1/28 1/14 1/8 1/4 1/8 ; 1/28 3/28 5/28 1/8 9/56 5/56 1/14 5/28 1/8 5/14]

	A = Float64[9 1 3 9 3 3 7 3 5 5; 0 2 4 7 1 1 8 4 7 1; 0 6 3 7 0 9 0 2 0 6; 8 8 7 6 5 9 1 4 2 5];

	G = neighborhoodgraph(A, 0, 0.0)
	@test G == Diagonal(trues(10))
	S = simplexgraph2kernelmatrix(G)
	@test S == Diagonal(ones(10))

	G = neighborhoodgraph(A, 1, 0.0)
	@test G == GAnswer1
	S = simplexgraph2kernelmatrix(G)
	@test S ≈ SAnswer1

	G = neighborhoodgraph(A, 2, 0.0)
	@test G == GAnswer2
	S = simplexgraph2kernelmatrix(G)
	@test S ≈ SAnswer2

	G = neighborhoodgraph(A, 6, 0.0)
	@test G == GAnswer6
	S = simplexgraph2kernelmatrix(G)
	@test S ≈ SAnswer6

	G = neighborhoodgraph(A, 9, 0.0)
	@test G == trues(10,10)
	S = simplexgraph2kernelmatrix(G)
	@test S ≈ (ones(10,10)+Diagonal(ones(10)))/11
end



end
