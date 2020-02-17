module PrincipalMomentAnalysis

using LinearAlgebra
using SparseArrays
using Statistics

include("graph.jl")
include("pma.jl")
include("util.jl")

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
	normalizemeanstd!

end
