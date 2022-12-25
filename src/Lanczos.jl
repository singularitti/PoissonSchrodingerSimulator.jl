module Lanczos

using LinearAlgebra: SymTridiagonal, norm, normalize, eigen, ⋅
using ProgressMeter: @showprogress

export lanczos, restart_lanczos, loop_lanczos

function lanczos(A::AbstractMatrix, 𝐪₁=normalize(rand(size(A, 1))), β₁=0; maxiter=30)
    n = 1  # Initial step
    𝐪₁ = normalize(𝐪₁)
    Q = Matrix{eltype(𝐪₁)}(undef, size(A, 1), maxiter)  # N × M
    Q[:, 1] = 𝐪₁
    𝐩₁ = A * 𝐪₁
    α₁ = 𝐪₁ ⋅ 𝐩₁  # 𝐪ₙ⊺ A 𝐪ₙ
    𝐫ₙ = 𝐩₁ - α₁ * 𝐪₁  # 𝐫₁, Gram–Schmidt process
    𝛂 = Vector{eltype(float(α₁))}(undef, maxiter)
    𝛃 = Vector{eltype(float(β₁))}(undef, maxiter)
    𝛂[n], 𝛃[n] = α₁, β₁
    for n in 2:maxiter
        𝐫ₙ₋₁ = 𝐫ₙ
        𝛃[n] = norm(𝐫ₙ₋₁)
        if iszero(𝛃[n])
            error("")
        else
            𝐪ₙ = 𝐫ₙ₋₁ / 𝛃[n]
            Q[:, n] = 𝐪ₙ
        end
        𝐩ₙ = A * 𝐪ₙ
        𝛂[n] = 𝐪ₙ ⋅ 𝐩ₙ  # 𝐪ₙ⊺ A 𝐪ₙ
        𝐫ₙ = 𝐩ₙ - 𝛂[n] * 𝐪ₙ - 𝛃[n] * Q[:, n - 1]
    end
    T = SymTridiagonal(𝛂, 𝛃)
    return T, Q
end

recover_eigvec(Q, 𝐰) = normalize(Q * 𝐰)

function restart_lanczos(T, Q)
    vals, vecs = eigen(T)
    if all(vals .> 0)
        index = argmin(vals)  # Index of the smallest eigenvalue
    else
        index = argmax(abs.(vals))  # Index of the smallest eigenvalue
    end
    𝐰 = vecs[:, index]  # Associated eigenvector
    return recover_eigvec(Q, 𝐰)
end

function loop_lanczos(
    A::AbstractMatrix, n, 𝐪₁=normalize(rand(size(A, 1))), β₁=0; maxiter=30
)
    @showprogress for _ in 1:n
        T, Q = lanczos(A, 𝐪₁, β₁; maxiter=maxiter)
        𝐪₁ = restart_lanczos(T, Q)
    end
    return 𝐪₁
end

end
