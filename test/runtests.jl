using PrincipalMomentAnalysis
using Test
using LinearAlgebra

@test [] == detect_ambiguities(Base, Core, PrincipalMomentAnalysis)

include("test_graph.jl")
include("test_pma.jl")
