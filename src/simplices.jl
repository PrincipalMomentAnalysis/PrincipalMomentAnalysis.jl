"""
	groupsimplices(groupby::AbstractVector)

Create simplex graph connecting elements with identical values in the groupby vector.

# Examples
```
julia> G = groupsimplices(["A","A","B","C","B"])
5×5 BitArray{2}:
 1  1  0  0  0
 1  1  0  0  0
 0  0  1  0  1
 0  0  0  1  0
 0  0  1  0  1
```
"""
function groupsimplices(groupby::AbstractVector)
	N = length(groupby)
	G = falses(N,N) # make sparse?
	for g in unique(groupby)
		ind = findall( groupby.==g )
		G[ind,ind] .= true
	end
	G
end


"""
	timeseriessimplices(time::AbstractVector; groupby::AbstractVector)

Create simplex graph connecting elements adjacent in time.
In case of ties, **All** elements at a unique timepoint will be connected to the **all** elements at the previous, current and next timepoints.
If `groupby` is specified, the elements are first divided by group, and then connected by time.

# Examples
```
julia> G = timeseriessimplices([0.5, 1.0, 4.0])
3×3 BitArray{2}:
 1  1  0
 1  1  1
 0  1  1

julia> G = timeseriessimplices([0.5, 1.0, 1.0, 4.0])
4×4 BitArray{2}:
 1  1  1  0
 1  1  1  1
 1  1  1  1
 0  1  1  1

julia> G = timeseriessimplices([2, 4, 6, 2, 4, 8]; groupby=["A","A","A","B","B","B"])
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

	G = falses(N,N) # make sparse?

	for g in unique(groupby)
		ind = findall( groupby.==g )
		time2 = time[ind]
		
		prevInd = Int[]
		for t in unique(sort(time2))
			currInd = ind[time2.==t]

			G[currInd,currInd] .= true
			G[currInd,prevInd] .= true
			G[prevInd,currInd] .= true
			prevInd = currInd
		end
	end

	G
end


"""
	neighborsimplices2(D2; k, r, symmetric, normalizedist, groupby)

Create simplex graph connecting nearest neighbors in given symmetric matrix where element `i,j` equals the squared distance between samples `i` and `j`.

# Inputs
* `D2`: Matrix of squared distances.
* `k`: Number of nearest neighbors to connect. Default: `0`.
* `r`: Connected all neighbors with disctance `≤r`. Default: `0.0`.
* `symmetric`: Make the simplex graph symmetric. Default: `false`.
* `normalizedist`: Normalize distances to the scale [0.0,1.0]. Affects the `r` parameter. Default: `true`.
* `groupby`: Only connected samples within the specified groups. Default: Disabled.
"""
function neighborsimplices2(D2::AbstractMatrix, k::Integer, r::Real; symmetric=false, normalizedist=true, groupby=falses(size(D2,1)))
	@assert issymmetric(D2)
	@assert all(>=(0.0), D2)
	N = size(D2,1)
	r2 = r*r * (normalizedist ? maximum(D2) : 1.0)

	uniqueGroups = unique(groupby)
	groupInds = Dict( g=>findall(groupby.==g) for g in uniqueGroups )

	G = falses(N,N)
	for j=1:N
		gInds = groupInds[groupby[j]]
		ind = gInds[sortperm(D2[gInds,j])]
		kk = max(k+1, searchsortedlast(D2[ind,j], r2)) # k+1 to include current node
		G[ind[1:min(kk,length(ind))], j] .= true
	end

	symmetric && (G .|= G')
	G
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
* `normalizedist`: Normalize distances to the scale [0.0,1.0]. Affects the `r` parameter. Default: `true`.
* `groupby`: Only connected samples within the specified groups. Default: Disabled.

# Examples
```
julia> G = neighborsimplices([0 0 2 2; 0 1 1 0]; k=1)
4×4 BitArray{2}:
 1  1  0  0
 1  1  0  0
 0  0  1  1
 0  0  1  1

julia> G = neighborsimplices([0 0 2 2; 0 1 1 0]; r=0.9)
4×4 BitArray{2}:
 1  1  0  1
 1  1  1  0
 0  1  1  1
 1  0  1  1
```
"""
function neighborsimplices(A::AbstractMatrix; k::Integer=0, r::Real=0.0, dim::Integer=typemax(Int), kwargs...)
	P,N = size(A)
	K = A'A

	if dim<size(K,1)
		F = eigen(Symmetric(K), N-dim+1:N)
		K = F.vectors*Diagonal(F.values)*F.vectors'
	end

	d = diag(K)
	D2 = Symmetric(max.(0., d .+ d' .- 2K)) # matrix of squared distances
	neighborsimplices2(D2,k,r;kwargs...)
end




function sparseneighborsimplices2(D2::Symmetric, k::Integer, r::Float64; symmetric=false, normalizedist=true)
	@assert all(>=(0.0), D2)
	N = size(D2,1)
	r2 = r*r * (normalizedist ? maximum(D2) : 1.0)

	I,J = Int[],Int[]
	for j=1:N
		ind = sortperm(D2[:,j])
		kk = max(k+1, searchsortedlast(D2[ind,j], r2)) # k+1 to include current node
		rows = ind[1:min(kk,length(ind))]
		append!(I,rows)
		append!(J,Iterators.repeated(j,length(rows)))
	end
	G = sparse(I,J,trues(length(I)))

	symmetric && (G .|= G')
	G
end
function sparseneighborsimplices(A::AbstractMatrix, k::Integer, r::Float64, dim::Integer=typemax(Int); kwargs...)
	P,N = size(A)


	K = if dim>=min(N,P)
		A'A # no dimension reduction needed
	elseif N<=P
		F = eigen(Symmetric(A'A), N-dim+1:N) # A'A is small
		F.vectors*Diagonal(F.values)*F.vectors'
	else # if P<N
		F = eigen(Symmetric(A*A'), P-dim+1:P) # AA' is small
		U = F.vectors
		VΣ = A'U
		VΣ*VΣ'
	end

	d = diag(K)
	D2 = Symmetric(max.(0., d .+ d' .- 2K)) # matrix of squared distances
	sparseneighborsimplices2(D2,k,r;kwargs...)
end
