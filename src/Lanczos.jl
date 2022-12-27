module Lanczos

using LinearAlgebra: SymTridiagonal, norm, normalize, eigen, ⋅
using ProgressMeter: @showprogress

export lanczos, restart_lanczos, loop_lanczos

function lanczos(A::AbstractMatrix, 𝐯₁=rand(size(A, 1)); maxiter=30)
    n = 1  # Initial step
    𝐯₁ = normalize(𝐯₁)
    V = Matrix{eltype(𝐯₁)}(undef, length(𝐯₁), maxiter)
    V[:, n] = 𝐯₁
    𝐰′₁ = A * 𝐯₁
    α₁ = 𝐰′₁ ⋅ 𝐯₁   # 𝐯₁⊤ A 𝐯₁
    𝐰ₙ = 𝐰′₁ - α₁ * 𝐯₁  # 𝐰₁, Gram–Schmidt process
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
            V[:, n] = 𝐯ₙ
        end
        𝐰′ₙ = A * 𝐯ₙ
        𝛂[n] = 𝐰′ₙ ⋅ 𝐯ₙ  # 𝐯ₙ⊤ A 𝐯ₙ
        𝐰ₙ = 𝐰′ₙ - 𝛂[n] * 𝐯ₙ - 𝛃[n] * V[:, n - 1]
    end
    T = SymTridiagonal(𝛂, 𝛃[2:end])
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

function loop_lanczos(A::AbstractMatrix, n, 𝐯₁=rand(size(A, 1)); maxiter=30)
    @showprogress for _ in 1:n
        T, V = lanczos(A, 𝐯₁; maxiter=maxiter)
        𝐯₁ = restart_lanczos(T, V)
    end
    return 𝐯₁
end

end
