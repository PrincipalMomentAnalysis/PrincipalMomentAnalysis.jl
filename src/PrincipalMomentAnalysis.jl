module PrincipalMomentAnalysis

using LinearAlgebra
using SparseArrays
using Statistics
using Distributed

include("graph.jl")
include("pma.jl")
include("util.jl")
include("projectionscore.jl")

export
	pma,
	_pma, # WIP
	pma2, # WIP
	_pma2, # WIP
	buildgraph,
	neighborhoodgraph,
	sparseneighborhoodgraph,
	simplexgraph2kernelmatrix,
	simplexgraph2kernelmatrixroot,
	normalizemean!,
	normalizemeanstd!,
	projectionscore,
	projectionscorefiltered

end
