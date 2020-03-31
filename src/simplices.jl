"""
	SimplexGraph{T<:AbstractMatrix{Bool}, S<:Union{Number,AbstractVector{<:Number}}}

Struct describing a set of simplices with vertices numbered from 1:N.
Optionally, each simplex can be assigned a weight.

Usage
--------------

	SimplexGraph(G, w=true)

Construct a SimplexGraph. `G` is a `MxN` matrix of booleans, where `M` is the number of vertices and `N` is the number of simplices.
Each column in `G` represents one simplex. True means that a vertex is part of the simplex and false that it is not.
Optionally, a vector `w` of length N can be used to specify a weight (total mass) for each simplex. The weight defaults to 1 (true).

See also `pma`, `groupsimplices`, `timeseriessimplices`, `neighborsimplices`.
"""
struct SimplexGraph{T<:AbstractMatrix{Bool}, S<:Union{Number,AbstractVector{<:Number}}}
	G::T
	w::S

	function SimplexGraph(G::T, w::S=true) where {T<:AbstractMatrix{Bool}, S<:Union{Number,AbstractVector{<:Number}}}
		@assert ndims(w)==0 || (size(w,1)==size(G,2))
		new{T,S}(G,w)
	end
end


"""
	groupsimplices(groupby::AbstractVector)

Create SimplexGraph connecting elements with identical values in the groupby vector.

# Examples
```
julia> sg = groupsimplices(["A","A","B","C","B"]);

julia> sg.G
5×3 BitArray{2}:
 1  0  0
 1  0  0
 0  1  0
 0  0  1
 0  1  0

julia> sg.w
3-element Array{Int64,1}:
 2
 2
 1
```
"""
function groupsimplices(groupby::AbstractVector)
	u = unique(groupby)
	ind = indexin(groupby,unique(groupby))
	G = falses(length(groupby), length(u)) # TODO: make sparse?
	w = zeros(Int,length(u))
	for (i,k) in enumerate(ind)
		G[i,k] = true
		w[k] += 1
	end
	SimplexGraph(G,w)
end


"""
	timeseriessimplices(time::AbstractVector; groupby::AbstractVector)

Create SimplexGraph connecting elements adjacent in time.
In case of ties, **all** elements at a unique timepoint will be connected to the **all** elements at the previous, current and next timepoints.
If `groupby` is specified, the elements are first divided by group, and then connected by time.

# Examples
```
julia> sg = timeseriessimplices([0.5, 1.0, 4.0]);

julia> sg.G
3×3 BitArray{2}:
 1  1  0
 1  1  1
 0  1  1

julia> sg = timeseriessimplices([0.5, 1.0, 1.0, 4.0]);

julia> sg.G
4×3 BitArray{2}:
 1  1  0
 1  1  1
 1  1  1
 0  1  1

julia> sg.w
3-element Array{Int64,1}:
 1
 2
 1

julia> sg = timeseriessimplices([2, 4, 6, 2, 4, 8]; groupby=["A","A","A","B","B","B"]);

julia> sg.G
6×6 BitArray{2}:
 1  1  0  0  0  0
 1  1  1  0  0  0
 0  1  1  0  0  0
 0  0  0  1  1  0
 0  0  0  1  1  1
 0  0  0  0  1  1
```
"""
function timeseriessimplices(time::AbstractVector; groupby::AbstractVector=falses(length(time)))
	N = length(groupby)

	G = BitVector() # make sparse?
	w = Int[]

	# TODO: Rewrite. Can be simplified and optimized.
	for g in unique(groupby)
		ind = findall( groupby.==g )
		time2 = time[ind]
		uniqueTimes = sort(unique(time2))

		currInd = Int[]
		nextInd = ind[time2.==uniqueTimes[1]]
		s = falses(length(time))
		s[nextInd] .= true

		for t in Iterators.drop(uniqueTimes,1)
			prevInd = currInd
			currInd = nextInd
			nextInd = ind[time2.==t]
			s[nextInd] .= true
			append!(G,s)
			push!(w, length(currInd))
			s[prevInd] .= false
		end
		if length(uniqueTimes)==2 # Two timepoints means the same simplex is used for both.
			w[end] += length(nextInd)
		else
			append!(G,s)
			push!(w, length(nextInd))
		end
	end

	G = reshape(G, length(time), :)
	SimplexGraph(G,w)
end



"""
	neighborsimplices(A::AbstractMatrix; k, r, dim, symmetric, normalizedist, groupby)

Create simplex graph connecting nearest neighbor samples.

# Inputs
* `A`: Data matrix (variables × samples).
* `k`: Number of nearest neighbors to connect. Default: `0`.
* `r`: Connected all neighbors with disctance `≤r`. Default: `0.0`.
* `dim`: Reduce the dimension to `dim` before computing distances. Useful to reduce noise. Default: Disabled.
* `symmetric`: Make the simplex graph symmetric. Default: `false`.
* `normalizedist`: Normalize distances to the scale [0.0,1.0] such that the maximal distance from a point to the origin is 0.5. Affects the `r` parameter. Default: `true`.
* `groupby`: Only connected samples within the specified groups. Default: Disabled.

# Examples
```
julia> sg = neighborsimplices([0 0 2 2; 0 1 1 0]; k=1); sg.G
4×4 BitArray{2}:
 1  1  0  0
 1  1  0  0
 0  0  1  1
 0  0  1  1

julia> sg = neighborsimplices([0 0 2 2; 0 1 1 0]; r=0.45); sg.G
4×4 BitArray{2}:
 1  1  0  1
 1  1  1  0
 0  1  1  1
 1  0  1  1
```
"""
function neighborsimplices(A::AbstractMatrix; dim::Integer=typemax(Int), kwargs...)
	if dim<minimum(size(A))
		F = svdbyeigen(A; nsv=dim)
		A = Diagonal(F.S)*F.Vt
	end
	_neighborsimplices(A; kwargs...)
end

function _neighborsimplices(A::AbstractMatrix; k::Integer=0, r::Real=0.0, symmetric=false, normalizedist=true, groupby=nothing)#falses(size(A,2)))
	k==0 && r==0.0 && return SimplexGraph(sparse(Diagonal(trues(size(A,2))))) # trivial case, avoid some computations.
	r>0.0 && normalizedist && (r *= 2*sqrt(maximum(sum(x->x^2, A, dims=1))))

	N = size(A,2)
	G = falses(N,N)

	if groupby===nothing
		_neighborsimplices!(G, A, identity, k, r)
	else
		for g in unique(groupby)
			ind = findall( groupby.==g )
			_neighborsimplices!(G, A[:,ind], x->ind[x], k, r)
		end
	end

	SimplexGraph(G)
end

function _neighborsimplices!(G, A, indfun, k, r)
	tree = BallTree(A)

	if k>0
		indices,_ = knn(tree, A, min(k+1,size(A,2)))
		for (j,I) in enumerate(indices)
			G[indfun.(I),indfun(j)] .= true
		end
	end

	if r>0.0
		indices = inrange(tree, A, r)
		for (j,I) in enumerate(indices)
			G[indfun.(I),indfun(j)] .= true
		end
	end
end
