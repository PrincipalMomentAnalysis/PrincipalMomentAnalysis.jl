using PrincipalMomentAnalysis
using Test
using LinearAlgebra
using Statistics

@test [] == detect_ambiguities(Base, Core, PrincipalMomentAnalysis)

const simplices2kernelmatrix     = PrincipalMomentAnalysis.simplices2kernelmatrix
const simplices2kernelmatrixroot = PrincipalMomentAnalysis.simplices2kernelmatrixroot

include("test_simplices.jl")
include("test_pma.jl")
include("test_util.jl")
