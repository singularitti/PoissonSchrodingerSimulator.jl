module Lanczos

using LinearAlgebra: norm, normalize, ⋅
using SparseArrays: spdiagm

function lanczos(A::AbstractMatrix, M=size(A, 2), 𝐪₁=normalize(rand(M)), β₁=0)
    n = 1  # Initial step
    𝐪₁ = normalize(𝐪₁)
    Q = Matrix{eltype(𝐪₁)}(undef, size(A, 1), M)  # N × M
    Q[:, 1] = 𝐪₁
    𝐩₁ = A * 𝐪₁
    α₁ = 𝐪₁ ⋅ 𝐩₁  # 𝐪ₙ⊺ A 𝐪ₙ
    𝐫ₙ = 𝐩₁ - α₁ * 𝐪₁  # 𝐫₁
    𝛂 = Vector{eltype(float(α₁))}(undef, M)
    𝛃 = Vector{eltype(float(β₁))}(undef, M)
    𝛂[n], 𝛃[n] = α₁, β₁
    for n in 2:M
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
    T = spdiagm(0 => 𝛂, 1 => 𝛃, -1 => 𝛃)
    return T, Q
end

end
