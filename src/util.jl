normalizemean!(X) = X .-= mean(X,dims=2)
normalizemean(X) = normalizemean!(copy(X))

function normalizemeanstd!(X)
	normalizemean!(X)
	X ./= max.(1e-12, std(X,dims=2)) # max to avoid div by zero
end
normalizemeanstd(X) = normalizemeanstd!(copy(X))
