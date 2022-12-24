using Plots

using LastHomework.Lanczos
using LastHomework.QuantumMechanics

include("problem1.jl")

q = 0.001
H = Hamiltonian(A, 𝛟, q)

𝛙 = loop_lanczos(H, 40)
regionheatmap(myreshape(𝛙))
savefig("tex/plots/psi_heatmap.pdf")
surfaceplot(myreshape(𝛙))
savefig("tex/plots/psi_surface.pdf")
P = probability(𝛙)
regionheatmap(myreshape(P))
savefig("tex/plots/P_heatmap.pdf")
surfaceplot(myreshape(P))
savefig("tex/plots/P_surface.pdf")
