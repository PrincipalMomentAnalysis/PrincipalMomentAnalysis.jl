module PrincipalMomentAnalysis

using LinearAlgebra
using SparseArrays
using Statistics

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
	neighborsimplices2,
	sparseneighborsimplices,
	sparseneighborsimplices2,
	normalizemean!,
	normalizemean,
	normalizemeanstd!,
	normalizemeanstd

end
