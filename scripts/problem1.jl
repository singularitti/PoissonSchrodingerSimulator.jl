using SparseArrays

using Plots

using LastHomework
using LastHomework.ConjugateGradient
using LastHomework.Electrostatics

N = 128
𝛟₀ = ReshapeVector(zeros(N^2), N, N);
setbc!(𝛟₀, 0)
setsquare!(𝛟₀, 5)
𝛒₀ = ReshapeVector(zeros(N^2), N, N);
setcharges!(𝛒₀, -20)
A = sparse(DiscreteLaplacian(N))

𝛟, history = solve(A, -𝛒₀, 𝛟₀; maxiter=4)
regionheatmap(reshape(𝛟, N, N))
