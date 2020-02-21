function buildgraph(groupBy::AbstractVector)
	N = length(groupBy)
	G = falses(N,N) # make sparse?
	for g in unique(groupBy)
		ind = findall( groupBy.==g )
		G[ind,ind] .= true
	end
	G
end

function buildgraph(groupBy::AbstractVector, time::AbstractVector)
	N = length(groupBy)

	G = falses(N,N) # make sparse?

	for g in unique(groupBy)
		ind = findall( groupBy.==g )
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


function neighborhoodgraph(D2::Symmetric, k::Integer, r::Float64; symmetric=false, groupBy=ones(size(D2,1)))
	N = size(D2,1)
	r2 = r*r

	uniqueGroups = unique(groupBy)
	groupInds = Dict( g=>findall(groupBy.==g) for g in uniqueGroups )


	G = falses(N,N)
	for j=1:N
		gInds = groupInds[groupBy[j]]
		ind = gInds[sortperm(D2[gInds,j])]
		kk = max(k+1, searchsortedlast(D2[ind,j], r2)) # k+1 to include current node
		G[ind[1:min(kk,length(ind))], j] .= true
	end

	symmetric && (G .|= G')
	G
end
function neighborhoodgraph(X::AbstractMatrix, k::Integer, r::Float64, dim::Integer=typemax(Int); kwargs...)
	P,N = size(X)
	K = X'X

	if dim<size(K,1)
		F = eigen(Symmetric(K), N-dim+1:N)
		K = F.vectors*Diagonal(F.values)*F.vectors'
	end

	d = diag(K)
	D2 = Symmetric(max.(0., d .+ d' .- 2K)) # matrix of squared distances
	neighborhoodgraph(D2,k,r;kwargs...)
end




function sparseneighborhoodgraph(D2::Symmetric, k::Integer, r::Float64; symmetric=false)
	N = size(D2,1)
	r2 = r*r

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
function sparseneighborhoodgraph(X::AbstractMatrix, k::Integer, r::Float64, dim::Integer=typemax(Int); kwargs...)
	P,N = size(X)


	K = if dim>=min(N,P)
		X'X # no dimension reduction needed
	elseif N<=P
		F = eigen(Symmetric(X'X), N-dim+1:N) # X'X is small
		F.vectors*Diagonal(F.values)*F.vectors'
	else # if P<N
		F = eigen(Symmetric(X*X'), P-dim+1:P) # XX' is small
		U = F.vectors
		VΣ = X'U
		VΣ*VΣ'
	end

	d = diag(K)
	D2 = Symmetric(max.(0., d .+ d' .- 2K)) # matrix of squared distances
	sparseneighborhoodgraph(D2,k,r;kwargs...)
end
