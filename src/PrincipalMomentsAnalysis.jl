module PrincipalMomentsAnalysis

using LinearAlgebra
using Statistics

include("graph.jl")
include("pma.jl")

export
	pma,
	buildgraph,
	neighborhoodgraph,
	graph2simplices

end
