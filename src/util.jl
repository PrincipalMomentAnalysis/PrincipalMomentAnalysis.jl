"""
	normalizemean!(A)

In place removal of the variable mean from the P×N matrix `A` where P is the number of variables and N is the number of samples.
"""
normalizemean!(A::AbstractMatrix) = A .-= mean(A,dims=2)

"""
	normalizemean(A)

Remove the variable mean from the P×N matrix `A` where P is the number of variables and N is the number of samples.
"""
normalizemean(A::AbstractMatrix) = A .- mean(A,dims=2)

"""
	normalizemeanstd!(A)

In place normalization of variables to be mean zero and standard deviation one of the P×N matrix `A` where P is the number of variables and N is the number of samples.
"""
function normalizemeanstd!(A::AbstractMatrix)
	normalizemean!(A::AbstractMatrix)
	A ./= max.(1e-12, std(A,dims=2)) # max to avoid div by zero
end
"""
	normalizemeanstd!(A)

Normalize variables to be mean zero and standard deviation one in the P×N matrix `A` where P is the number of variables and N is the number of samples.
"""
function normalizemeanstd(A::AbstractMatrix)
	A = normalizemean(A)
	A ./= max.(1e-12, std(A,dims=2)) # max to avoid div by zero
end
