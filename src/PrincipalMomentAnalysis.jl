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
