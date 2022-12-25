using Plots

using LastHomework
using LastHomework.ConjugateGradient
using LastHomework.Electrostatics

step_cg = 500
logger = Logger(step_cg)
L = 128
N = L + 1  # Grid size

myreshape(𝐯) = reshape(𝐯, N, N)

𝛟₀ = zeros(N^2);
𝛒 = zeros(N^2);
A = DiscreteLaplacian(N);

𝛟 = solve!(logger, A, 𝛒, 𝛟₀; maxiter=step_cg)
regionheatmap(myreshape(𝛟₀))
savefig("tex/plots/phi0_heatmap.pdf")
regionheatmap(myreshape(𝛒))
savefig("tex/plots/rho_heatmap.pdf")
regionheatmap(myreshape(𝛟))
savefig("tex/plots/phi_heatmap.pdf")
surfaceplot(myreshape(𝛟); xmirror=true, ymirror=true, view_angle=(225, 30))
savefig("tex/plots/phi_surface.pdf")
residualplot(logger)
savefig("tex/plots/residual.pdf")
# conjugacyplot(A, logger)
