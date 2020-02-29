using PrincipalMomentAnalysis
using Test
using LinearAlgebra

@test [] == detect_ambiguities(Base, Core, PrincipalMomentAnalysis)

const simplices2kernelmatrix     = PrincipalMomentAnalysis.simplices2kernelmatrix
const simplices2kernelmatrixroot = PrincipalMomentAnalysis.simplices2kernelmatrixroot

include("test_simplices.jl")
include("test_pma.jl")
