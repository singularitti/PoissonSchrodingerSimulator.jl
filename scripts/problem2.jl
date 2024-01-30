using Plots

using PoissonSchrodingerSimulator.Lanczos
using PoissonSchrodingerSimulator.QuantumMechanics

include("problem1.jl")

q = 0.001
H = Hamiltonian(A, 𝛟, q)

ntimes = 40

𝛙 = loop_lanczos(H, ntimes)
regionheatmap(myreshape(𝛙))
savefig("tex/plots/psi_heatmap.pdf")
surfaceplot(myreshape(𝛙))
savefig("tex/plots/psi_surface.pdf")
P = probability(𝛙)
regionheatmap(myreshape(P); right_margin=(3, :mm))
savefig("tex/plots/P_heatmap.pdf")
surfaceplot(myreshape(P); right_margin=(3, :mm))
savefig("tex/plots/P_surface.pdf")

print(probability(myreshape(𝛙), 1:N, 1:(N ÷ 2)))
