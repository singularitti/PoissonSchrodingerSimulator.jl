using Plots

using LastHomework
using LastHomework.ConjugateGradient
using LastHomework.Electrostatics

N = 16
ϕ₀ = SolutionMatrix(zeros(N, N));
set!(ϕ₀, Boundary(0));
set!(ϕ₀, InternalSquare(5));
surfaceplot(ϕ₀)
ρ₀ = ResidualMatrix(zeros(N, N));
set!(ρ₀, PointCharges(-20));
surfaceplot(ρ₀)
A = DiscreteLaplacian(N);

ϕ = solve(A, -vec(ρ₀), vec(ϕ₀); maxiter=1000)
regionheatmap(ϕ)
surfaceplot(ϕ)
