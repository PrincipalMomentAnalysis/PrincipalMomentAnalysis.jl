@testset "PMA" begin

function factorizationcmp(F1,F2)
	@test size(F1.U)==size(F2.U)
	@test size(F1.S)==size(F2.S)
	@test size(F1.V)==size(F2.V)
	@test F1.S ≈ F2.S
	sgn = sign.(diag(F1.U'F2.U))
	@test F1.U.*sgn' ≈ F2.U
	@test F1.V.*sgn' ≈ F2.V
	nothing
end

@testset "PCA" begin
    A = [-8 0 6 0 0 0; -1 9 6 -2 0 -1; -1 0 -2 0 0 0; 0 1 0 1 0 -2; -9 -3 0 0 0 -5; -6 3 -8 -2 -4 8; 0 -1 -1 0 -3 -1; 4 0 0 -6 0 2]
    SVD = svd(A)
    PMA = pma(A,I(6); nsv=6)
    factorizationcmp(SVD,PMA)
end

end
