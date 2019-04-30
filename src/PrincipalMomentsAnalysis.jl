module PrincipalMomentAnalysis

using LinearAlgebra
using Statistics

include("graph.jl")
include("pma.jl")
include("util.jl")

export
	pma,
	buildgraph,
	neighborhoodgraph,
	graph2simplices,
	normalizemean!,
	normalizemeanstd!

end
