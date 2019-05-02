function α(s::Integer,K::Symmetric{Float64})
    N = size(K,2)
    s = min(s,N) # any value s>N will get an α score of 1
    eig = reverse(eigvals(K, N-s+1:N)) # get the largest eigenvalues (as many as we need) and sort them in decreasing order
    num = sum(eig) # ∑ᵢ₌₁ˢ λᵢ²
    denom = tr(K) # ∑_ᵢ₌₁ⁿ λᵢ²
    sqrt(num/denom)
end
α(s::Integer, X::AbstractMatrix) = α(s,Symmetric(X'*X))


function _samplematrix(P::Int,N::Int,S::AbstractMatrix{Float64})
    X = randn(P,N)
    X .-= mean(X,dims=2)
    X ./= std(X,dims=2)
    X *= S
    X
end


function _αbootstrap(P::Int, N::Int, S::AbstractMatrix{Float64}, s::Integer, nbrIter::Integer)
    sum(i->α(s,_samplematrix(P,N,S)), 1:nbrIter) # bootstrap and sum
end


function _projectionscore(X::AbstractMatrix, S::AbstractMatrix{Float64}, s::Integer; nbrIter::Integer=10)
    @assert size(X,2)==size(S,1)==size(S,2)
    P = size(X,1)
    N = size(X,2)

    αB = NaN
    if nprocs()==1 || nbrIter==1
        αB=_αbootstrap(P,N,S,s,nbrIter)/nbrIter # no threading needed
    else
        # try to minimize data movement by spawning once per process
        W = workers()
        nbrWorkers = min( length(W), nbrIter ) # never use more than nbrIter workers
        nbrIterPerWorker = [ ((d,r)=divrem(nbrIter,nbrWorkers); d+(r>=i)) for i=1:nbrWorkers ] # divide iterations evenly among workers

        # TODO: do with pmap instead?
        refs = Array{Any}(nbrWorkers)
        for i=1:nbrWorkers
            refs[i] = @spawnat W[i] _αbootstrap(P,N,S,s,nbrIterPerWorker[i])
        end
        αB = sum(fetch, refs)/nbrIter # take mean over all bootstrapped α's
    end
    α(s,X*S) - αB
end

function projectionscore(X::AbstractMatrix, G::AbstractMatrix, s::Integer; kwargs...)
    @assert maximum(abs.(sum(X,dims=2))) < 1e-9 "X must be centered"
    @assert all(std(X,dims=2).-1.0 .< 1e-9) "All variables in X must have standard deviation 1"
    _projectionscore(X, graph2simplices(G), s; kwargs...)
end



# Assumes X is sorted by decreasing σ!
# Also assumes X has already been multiplied with Simplex matrix S.
function _αfiltered(X, s, σ, σThresholds)
    αs = zeros(length(s), length(σThresholds))

    # Build Kernel Matrix gradually, starting with the variables with highest σ
    K = zeros(size(X,2),size(X,2))
    prevInd = 1
    for i = length(σThresholds):-1:1
        # get range of variables to add
        currInd = prevInd
        while currInd<length(σ) && σ[currInd] >= σThresholds[i]
            currInd += 1
        end

        varRange = prevInd:currInd-1

        if !isempty(varRange)
            X2 = X[varRange,:]
            K += X2'X2
            αs[:,i] .= α(s, Symmetric(K))
        elseif i<length(σThresholds)
            αs[:,i] .= αs[:,i+1]
        end
        prevInd = currInd
    end

    αs
end

function _αfilteredsum(P::Int,N::Int,S::AbstractMatrix{Float64},s,σ,σThresholds,nbrIter)
    sum(i->_αfiltered(_samplematrix(P,N,S),s,σ,σThresholds), 1:nbrIter)
end


function _projectionscorefiltered(X::AbstractMatrix, S::AbstractMatrix, variableStds::AbstractVector{Float64}, s::Integer, σThresholds::AbstractVector;
                                 nbrIter::Integer=10)
    @assert all(diff(σThresholds).>=0) "σThresholds must be increasing"
    @assert length(variableStds) == size(X,1)

    σ = variableStds ./ maximum(variableStds) # scale to 0-1

    nbrSamples = size(X,2)

    # get rid of those variables that will not be used at any threshold
    baseFilter = σ.>=σThresholds[1]
    X = X[baseFilter,:]
    σ = σ[baseFilter]

    if length(σ)==0
        warn("No variables remain after filtering at lowest given threshold.")
        return zeros(length(s), length(σThresholds))
    end

    # sort X so that the variables are ordered by decreasing σ
    σPerm = sortperm(σ, rev=true)
    σ = σ[σPerm]
    X = X[σPerm,:]

    α = _αfiltered(X*S,s,σ,σThresholds)

    P,N = size(X)

    αB = NaN

    if nprocs()==1 || nbrIter==1  # no threading needed
        αB = _αfilteredsum(P,N,S,s,σ,σThresholds,nbrIter)
    else
        # try to minimize data movement by spawning once per process
        W = workers()
        nbrWorkers = min( length(W), nbrIter ) # never use more than nbrIter workers
        nbrIterPerWorker = [ ((d,r)=divrem(nbrIter,nbrWorkers); d+(r>=i)) for i=1:nbrWorkers ] # divide iterations evenly among workers

        refs = Array{Any}(nbrWorkers)
        for i=1:nbrWorkers
            refs[i] = @spawnat W[i] _αfilteredsum(P,N,S,s,σ,σThresholds,nbrIterPerWorker[i])
        end
        αB = sum(fetch, refs) # take mean over all bootstrapped α's
    end

    α - αB/nbrIter
end


function projectionscorefiltered(X::AbstractMatrix, G::AbstractMatrix, variableStds::AbstractVector{Float64},
                                 s::Integer, σThresholds::AbstractVector;
                                 kwargs...)
    nbrSamples = size(X,2)
    @assert maximum(abs.(sum(X,dims=2))) < 1e-9 "X must be centered"
    @assert all(std(X,dims=2).-1.0 .< 1e-9) "All variables in X must have standard deviation 1"
    _projectionscorefiltered(X,graph2simplices(G),variableStds,s,σThresholds;kwargs...)
end
