module Lanczos

using LinearAlgebra: SymTridiagonal, norm, normalize, eigen, ⋅
using ProgressMeter: @showprogress

export lanczos, restart_lanczos, loop_lanczos

function lanczos(A::AbstractMatrix, 𝐪₁=normalize(rand(size(A, 1))); maxiter=30)
    n = 1  # Initial step
    𝐪₁ = normalize(𝐪₁)
    V = Matrix{eltype(𝐪₁)}(undef, size(A, 1), maxiter)  # N × M
    V[:, 1] = 𝐪₁
    𝐩₁ = A * 𝐪₁
    α₁ = 𝐪₁ ⋅ 𝐩₁  # 𝐪ₙ⊺ A 𝐪ₙ
    𝐫ₙ = 𝐩₁ - α₁ * 𝐪₁  # 𝐫₁, Gram–Schmidt process
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
            V[:, n] = 𝐪ₙ
        end
        𝐩ₙ = A * 𝐪ₙ
        𝛂[n] = 𝐪ₙ ⋅ 𝐩ₙ  # 𝐪ₙ⊺ A 𝐪ₙ
        𝐫ₙ = 𝐩ₙ - 𝛂[n] * 𝐪ₙ - 𝛃[n] * V[:, n - 1]
    end
    T = SymTridiagonal(𝛂, 𝛃)
    return T, V
end

recover_eigvec(V, 𝐰) = normalize(V * 𝐰)

function restart_lanczos(T, V)
    vals, vecs = eigen(T)
    index = if all(vals .> 0)
        argmin(vals)  # Index of the smallest eigenvalue
    else
        argmax(abs.(vals))  # Index of the smallest eigenvalue
    end
    𝐰 = vecs[:, index]  # Associated eigenvector
    return recover_eigvec(V, 𝐰)
end

function loop_lanczos(A::AbstractMatrix, n, 𝐪₁=normalize(rand(size(A, 1))); maxiter=30)
    @showprogress for _ in 1:n
        T, V = lanczos(A, 𝐪₁; maxiter=maxiter)
        𝐪₁ = restart_lanczos(T, V)
    end
    return 𝐪₁
end

end
