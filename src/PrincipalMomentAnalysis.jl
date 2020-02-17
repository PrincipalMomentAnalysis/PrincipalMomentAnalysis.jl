module PrincipalMomentAnalysis

using LinearAlgebra
using SparseArrays
using Statistics

include("graph.jl")
include("pma.jl")
include("util.jl")

export
	PMA,
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
	normalizemean,
	normalizemeanstd!,
	normalizemeanstd

end
