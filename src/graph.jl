function buildgraph(groupBy::AbstractVector, time=ones(length(groupBy)))
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


function neighborhoodgraph(X::Matrix, k::Integer, r::Float64; symmetric=false)
	P,N = size(X)
	K = X'X
	d = diag(K)
	D2 = max(0., d .+ d' .- 2K) # matrix of squared distances
	r2 = r*r

	G = falses(N,N)
	for j=1:N

		ind = sortperm(D2[:,j])
		d2 = D2[ind,j]
		kk = max(k+1, searchsortedlast(d2, r2)) # k+1 to include current node
		G[ind[1:min(kk,end)], j] = true
	end

	symmetric && (G .|= G')
	G
end