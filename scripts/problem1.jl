using Plots

using LastHomework
using LastHomework.ConjugateGradient
using LastHomework.Electrostatics

logger = Logger(1000)
L = 32
N = L + 1  # Grid size
ϕ₀ = zeros(N^2);
ρ = zeros(N^2);
A = DiscreteLaplacian(N);

ϕ = solve!(logger, A, ρ, ϕ₀; maxiter=500)
regionheatmap(reshape(ϕ, N, N))
surfaceplot(reshape(ϕ, N, N))
conjugacyplot(A, logger)
