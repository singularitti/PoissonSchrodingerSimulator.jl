using Plots

using LastHomework
using LastHomework.Lanczos
using LastHomework.QuantumMechanics

L = 64
N = L + 1  # Grid size
q = 0.001
A = DiscreteLaplacian(N);
H = Hamiltonian(A, ϕ, q)

lanczos(H; maxiter=30)
