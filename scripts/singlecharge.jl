using Plots

using ConjugateGradient
using PoissonSchrodingerSimulator
using PoissonSchrodingerSimulator.Electrostatics: PointCharge

maxiter = 500
logger = Logger(maxiter)
L = 32
N = L + 1  # Grid size
𝛟₀ = zeros(N^2);
𝛒 = zeros(N^2);
A = DiscreteLaplacian(N);

myreshape(𝐯) = reshape(𝐯, N, N)

𝛟 = solve!(logger, A, 𝛒, 𝛟₀; maxiter=maxiter)
regionheatmap(myreshape(𝛟))
surfaceplot(myreshape(𝛟))
