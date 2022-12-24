using Plots

using LastHomework.Lanczos
using LastHomework.QuantumMechanics

include("problem1.jl")

q = 0.001
H = Hamiltonian(A, 𝛟, q)

𝛙 = loop_lanczos(H, 40)
P = probability(𝛙)
regionheatmap(myreshape(P))
savefig("tex/plots/psi_heatmap.pdf")
surfaceplot(myreshape(P))
savefig("tex/plots/psi_surface.pdf")
