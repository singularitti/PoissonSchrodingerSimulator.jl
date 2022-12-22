using SparseArrays

using Plots

using LastHomework
using LastHomework.ConjugateGradient
using LastHomework.Electrostatics

N = 128
𝛟₀ = SolutionVector(zeros(N^2));
set!(𝛟₀, Boundary(0));
set!(𝛟₀, InternalSquare(5));
𝛒₀ = ResidualVector(zeros(N^2));
set!(𝛒₀, PointCharges(-20));
A = sparse(DiscreteLaplacian(N));

𝛟, history = solve(A, -𝛒₀, 𝛟₀; maxiter=1000)
phimat = collect(reshape(𝛟, N, N))
regionheatmap(phimat)
surfaceplot(phimat)
surfaceplot!(collect(reshape(𝛟₀, N, N)))
