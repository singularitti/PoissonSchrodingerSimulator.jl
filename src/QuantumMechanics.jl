module QuantumMechanics

using LinearAlgebra: SymTridiagonal, norm, normalize, diagind, ⋅
using SparseArrays: AbstractSparseMatrix, SparseMatrixCSC

using ..LastHomework: DiscreteLaplacian, Boundary, InternalSquare, validate, setvalues!

import ..Lanczos: lanczos

export Hamiltonian, probability

struct Hamiltonian{T} <: AbstractSparseMatrix{T,Int64}
    parent::SparseMatrixCSC{T,Int64}
end
Hamiltonian(parent::AbstractMatrix{T}) where {T} =
    Hamiltonian{T}(SparseMatrixCSC{T}(parent))
function Hamiltonian(A::DiscreteLaplacian, 𝛟::AbstractVector, q::Number)
    H = float(A)
    H[diagind(H)] .+= q * 𝛟
    return Hamiltonian(H)
end

function lanczos(A::Hamiltonian, 𝐪₁=normalize(rand(size(A, 1))); maxiter=30)
    N = Int(sqrt(size(A, 1)))  # A is a N² × N² matrix
    BOUNDARY = Boundary((N, N), 0)
    SQUARE = InternalSquare((N, N), 0)
    n = 1  # Initial step
    setvalues!(𝐪₁, BOUNDARY)
    setvalues!(𝐪₁, SQUARE)
    𝐪₁ = normalize(𝐪₁)
    Q = Matrix{eltype(𝐪₁)}(undef, size(A, 1), maxiter)  # N² × M
    Q[:, 1] = 𝐪₁
    𝐩₁ = A * 𝐪₁
    setvalues!(𝐩₁, BOUNDARY)
    setvalues!(𝐩₁, SQUARE)
    α₁ = 𝐪₁ ⋅ 𝐩₁  # 𝐪ₙ⊺ A 𝐪ₙ
    𝐫ₙ = 𝐩₁ - α₁ * 𝐪₁  # 𝐫₁, Gram–Schmidt process
    validate(𝐫ₙ, BOUNDARY)
    validate(𝐫ₙ, SQUARE)
    # setvalues!(𝐫ₙ, BOUNDARY)
    # setvalues!(𝐫ₙ, SQUARE)
    𝛂 = Vector{eltype(float(α₁))}(undef, maxiter)
    𝛃 = Vector{Float64}(undef, maxiter)
    𝛂[n], 𝛃[n] = α₁, 0
    for n in 2:maxiter
        𝐫ₙ₋₁ = 𝐫ₙ
        𝛃[n] = norm(𝐫ₙ₋₁)
        if iszero(𝛃[n])
            error("")
        else
            𝐪ₙ = 𝐫ₙ₋₁ / 𝛃[n]
            validate(𝐫ₙ, BOUNDARY)
            validate(𝐫ₙ, SQUARE)
            Q[:, n] = 𝐪ₙ
        end
        𝐩ₙ = A * 𝐪ₙ
        setvalues!(𝐩ₙ, BOUNDARY)
        setvalues!(𝐩ₙ, SQUARE)
        𝛂[n] = 𝐪ₙ ⋅ 𝐩ₙ  # 𝐪ₙ⊺ A 𝐪ₙ
        𝐫ₙ = 𝐩ₙ - 𝛂[n] * 𝐪ₙ - 𝛃[n] * Q[:, n - 1]
        validate(𝐪ₙ, BOUNDARY)
        validate(𝐪ₙ, SQUARE)
    end
    T = SymTridiagonal(𝛂, 𝛃)
    return T, Q
end

probability(𝛙::AbstractVector) = abs2.(normalize(𝛙))

Base.parent(S::Hamiltonian) = S.parent

Base.size(S::Hamiltonian) = size(parent(S))

Base.IndexStyle(::Type{Hamiltonian{T}}) where {T} = IndexLinear()

Base.getindex(S::Hamiltonian, i) = getindex(parent(S), i)

Base.setindex!(S::Hamiltonian, v, i) = setindex!(parent(S), v, i)

end
