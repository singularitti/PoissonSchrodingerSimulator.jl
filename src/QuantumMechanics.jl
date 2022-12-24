module QuantumMechanics

using LinearAlgebra: SymTridiagonal, norm, normalize, diagind
using SparseArrays: AbstractSparseMatrix, SparseMatrixCSC

using ..LastHomework: DiscreteLaplacian

export Hamiltonian

struct Hamiltonian{T} <: AbstractSparseMatrix{T,Int64}
    parent::SparseMatrixCSC{T,Int64}
end
Hamiltonian(A::DiscreteLaplacian, 𝛟::AbstractVector, q::Number) =
    Hamiltonian(A[diagind(A)] .+ q * 𝛟)

end
