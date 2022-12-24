using Plots

using LastHomework
using LastHomework.ConjugateGradient
using LastHomework.Electrostatics

logger = Logger(1000)
L = 128
N = L + 1  # Grid size

myreshape(𝐯) = reshape(𝐯, N, N)

𝛟₀ = zeros(N^2);
𝛒 = zeros(N^2);
A = DiscreteLaplacian(N);

𝛟 = solve!(logger, A, 𝛒, 𝛟₀; maxiter=500)
regionheatmap(myreshape(𝛟))
savefig("tex/plots/phi_heatmap.pdf")
surfaceplot(myreshape(𝛟))
savefig("tex/plots/phi_surface.pdf")
residualplot(logger)
savefig("tex/plots/residual.pdf")
# conjugacyplot(A, logger)
