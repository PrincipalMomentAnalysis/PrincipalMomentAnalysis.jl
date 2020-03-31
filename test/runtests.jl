using PrincipalMomentAnalysis
using Test
using LinearAlgebra
using Statistics

# change back to this simple oneliner if we can just ignore problems in StaticArrays that is imported by NearestNeighbors
# @test [] == detect_ambiguities(Base, Core, PrincipalMomentAnalysis)
@testset "Method ambiguities" begin
	ambiguities = detect_ambiguities(Base, Core, PrincipalMomentAnalysis)
	filter!(x->x[1].module==PrincipalMomentAnalysis || x[2].module==PrincipalMomentAnalysis, ambiguities)
	@test ambiguities == []
end


const simplices2kernelmatrix     = PrincipalMomentAnalysis.simplices2kernelmatrix
const simplices2kernelmatrixroot = PrincipalMomentAnalysis.simplices2kernelmatrixroot

include("test_simplices.jl")
include("test_pma.jl")
include("test_util.jl")
