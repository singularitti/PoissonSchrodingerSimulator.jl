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

function lanczos(A::Hamiltonian, 𝐯₁=normalize(rand(size(A, 1))); maxiter=30)
    N = Int(sqrt(size(A, 1)))  # A is a N² × N² matrix
    BOUNDARY = Boundary((N, N), 0)
    SQUARE = InternalSquare((N, N), 0)
    n = 1  # Initial step
    setvalues!(𝐯₁, BOUNDARY)
    setvalues!(𝐯₁, SQUARE)
    𝐯₁ = normalize(𝐯₁)
    V = Matrix{eltype(𝐯₁)}(undef, length(𝐯₁), maxiter)  # N² × maxiter
    V[:, n] = 𝐯₁
    𝐰′₁ = A * 𝐯₁
    setvalues!(𝐰′₁, BOUNDARY)
    setvalues!(𝐰′₁, SQUARE)
    α₁ = 𝐰′₁ ⋅ 𝐯₁   # 𝐯₁⊺ A 𝐯₁
    𝐰ₙ = 𝐰′₁ - α₁ * 𝐯₁  # 𝐰₁, Gram–Schmidt process
    # validate(𝐰ₙ, BOUNDARY)
    # validate(𝐰ₙ, SQUARE)
    𝛂 = Vector{eltype(float(α₁))}(undef, maxiter)
    𝛃 = Vector{Float64}(undef, maxiter)
    𝛂[n], 𝛃[n] = α₁, 0
    for n in 2:maxiter
        𝐰ₙ₋₁ = 𝐰ₙ
        𝛃[n] = norm(𝐰ₙ₋₁)
        if iszero(𝛃[n])
            error("")
        else
            𝐯ₙ = 𝐰ₙ₋₁ / 𝛃[n]
            # validate(𝐯ₙ, BOUNDARY)
            # validate(𝐯ₙ, SQUARE)
            V[:, n] = 𝐯ₙ
        end
        𝐰′ₙ = A * 𝐯ₙ
        setvalues!(𝐰′ₙ, BOUNDARY)
        setvalues!(𝐰′ₙ, SQUARE)
        𝛂[n] = 𝐰′ₙ ⋅ 𝐯ₙ  # 𝐯ₙ⊺ A 𝐯ₙ
        𝐰ₙ = 𝐰′ₙ - 𝛂[n] * 𝐯ₙ - 𝛃[n] * V[:, n - 1]
        # validate(𝐰ₙ, BOUNDARY)
        # validate(𝐰ₙ, SQUARE)
    end
    T = SymTridiagonal(𝛂, 𝛃)
    return T, V
end

probability(𝛙::AbstractVector) = abs2.(normalize(𝛙))
function probability(𝛙::AbstractMatrix, xrange=1:size(𝛙, 1), yrange=1:size(𝛙, 2))
    𝛙′ = normalize(𝛙)
    return sum(abs2.(𝛙′[xrange, yrange]))
end

Base.parent(S::Hamiltonian) = S.parent

Base.size(S::Hamiltonian) = size(parent(S))

Base.IndexStyle(::Type{Hamiltonian{T}}) where {T} = IndexLinear()

Base.getindex(S::Hamiltonian, i) = getindex(parent(S), i)

Base.setindex!(S::Hamiltonian, v, i) = setindex!(parent(S), v, i)

end
