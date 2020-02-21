@testset "Graph" begin

@testset "Case1" begin
    annot = ["A","A","B","C"]
    G = buildgraph(annot)
    @test G==[true true 0 0; true true 0 0; 0 0 true 0; 0 0 0 true]
    S = simplexgraph2kernelmatrix(G)
    @test S ≈ [2/3 1/3 0 0; 1/3 2/3 0 0; 0 0 1 0; 0 0 0 1]
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

end