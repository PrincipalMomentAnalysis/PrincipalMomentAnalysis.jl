"""
	normalizemean!(X)

In place removal of the variable mean from the P×N matrix `X` where P is the number of variables and N is the number of samples.
"""
normalizemean!(X::AbstractMatrix) = X .-= mean(X,dims=2)

"""
	normalizemean(X)

Remove the variable mean from the P×N matrix `X` where P is the number of variables and N is the number of samples.
"""
normalizemean(X::AbstractMatrix) = normalizemean!(copy(X))

"""
	normalizemeanstd!(X)

In place normalization of variables to be mean zero and standard deviation one of the P×N matrix `X` where P is the number of variables and N is the number of samples.
"""
function normalizemeanstd!(X::AbstractMatrix)
	normalizemean!(X::AbstractMatrix)
	X ./= max.(1e-12, std(X,dims=2)) # max to avoid div by zero
end
"""
	normalizemeanstd!(X)

Normalize variables to be mean zero and standard deviation one in the P×N matrix `X` where P is the number of variables and N is the number of samples.
"""
normalizemeanstd(X::AbstractMatrix) = normalizemeanstd!(copy(X))
