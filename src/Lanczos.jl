module Lanczos

using LinearAlgebra: SymTridiagonal, norm, normalize, eigen, ⋅

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

function restart_lanczos(
    A::AbstractMatrix, 𝐪₁=normalize(rand(size(A, 1))), β₁=0; maxiter=30
)
    T, Q = lanczos(A, 𝐪₁, β₁; maxiter=maxiter)
    vals, vecs = eigen(T)
    index = argmin(vals)
    𝐰 = vecs[index]
    subspacedim = length(𝐰)
    return normalize(vec(𝐰' * Q[axes(𝐰, 1), begin:subspacedim]))
end

function loop_lanczos(
    A::AbstractMatrix, n, M=size(A, 2), 𝐪₁=normalize(rand(size(A, 1))), β₁=0; maxiter=30
)
    for _ in 1:n
        𝐪₁ = restart_lanczos(A, 𝐪₁, β₁; maxiter=maxiter)
    end
    return lanczos(A, 𝐪₁, β₁; maxiter=maxiter)  # Do one last Lanczos
end

end
