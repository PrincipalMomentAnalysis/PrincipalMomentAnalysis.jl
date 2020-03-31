module PrincipalMomentAnalysis

using LinearAlgebra
using SparseArrays
using Statistics
using NearestNeighbors

include("svd.jl")
include("simplices.jl")
include("pma.jl")
include("util.jl")

export
	PMA,
	pma,
	SimplexGraph,
	groupsimplices,
	timeseriessimplices,
	neighborsimplices,
	normalizemean!,
	normalizemean,
	normalizemeanstd!,
	normalizemeanstd

end
