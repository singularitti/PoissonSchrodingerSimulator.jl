using Plots

using LastHomework
using LastHomework.ConjugateGradient
using LastHomework.Electrostatics: PointCharge

logger = Logger(1000)
L = 32
N = L + 1  # Grid size
ϕ₀ = zeros(N^2);
boundary = Boundary((N, N), 0)
ϕ₀ = set(ϕ₀, boundary);
surfaceplot(ϕ₀)
ρ = zeros(N^2);
ρ = set(ρ, PointCharge((N, N), -20));
surfaceplot(ρ)
A = DiscreteLaplacian(N);

ϕ = solve!(logger, A, ρ, ϕ₀; maxiter=500)
regionheatmap(ϕ)
surfaceplot(ϕ)
