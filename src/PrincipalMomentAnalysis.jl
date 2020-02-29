module PrincipalMomentAnalysis

using LinearAlgebra
using SparseArrays
using Statistics

include("simplices.jl")
include("pma.jl")
include("util.jl")

export
	PMA,
	pma,
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
