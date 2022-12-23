using Plots

using LastHomework
using LastHomework.ConjugateGradient
using LastHomework.Electrostatics

logger = Logger(1000)
L = 32
N = L + 1  # Grid size
ϕ₀ = zeros(N^2);
boundary = Boundary((N, N), 0)
square = InternalSquare((N, N), 5)
set!(ϕ₀, boundary);
set!(ϕ₀, square);
surfaceplot(ϕ₀)
ρ = zeros(N^2);
set!(ρ, PointCharges((N, N), -20));
surfaceplot(ρ)
A = DiscreteLaplacian(N);

ϕ = solve!(logger, A, ρ, ϕ₀; maxiter=500)
regionheatmap(ϕ)
surfaceplot(ϕ)
conjugacyplot(A, logger)
