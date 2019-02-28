module PrincipalMomentsAnalysis

using LinearAlgebra

include("graph.jl")
include("pma.jl")

export
	pma,
	buildgraph,
	neighborhoodgraph,
	graph2simplices

end
