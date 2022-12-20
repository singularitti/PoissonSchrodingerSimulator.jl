module ConjugateGradient

using LinearAlgebra: dot, norm

struct IterationStep
    n::UInt64
    alpha::Float64
    beta::Float64
    x::Vector{Float64}
    r::Vector{Float64}
    p::Vector{Float64}
end

mutable struct ConvergenceHistory
    isconverged::Bool
    history::Vector{IterationStep}
end

function solve(A, 𝐛, 𝐱₀, ε=eps(), maxiter=2000)
    history = ConvergenceHistory(false, [])
    𝐱 = 𝐱₀
    𝐫 = 𝐛 - A * 𝐱  # Residual
    𝐩 = 𝐫  # Momentum
    α = compute_alpha(A, 𝐫, 𝐩)
    β = compute_beta(A, 𝐫, 𝐩)
    push!(history.history, IterationStep(0, α, β, 𝐱₀, 𝐫, 𝐩))
    for n in 1:maxiter
        𝐱 = 𝐱 + α * 𝐩  # Do not do in-place change!
        𝐫′ = 𝐫 - α * A * 𝐩  # Trial move
        if norm(𝐫′) < ε
            history.isconverged = true
        else
            β = compute_beta(𝐫′, 𝐫)
            𝐩 = 𝐫′ + β * 𝐩
            𝐫 = 𝐫′  # Accept the trial move
            push!(history.history, IterationStep(n, α, β, 𝐱, 𝐫, 𝐩))
        end
    end
    return 𝐱, history
end

compute_alpha(A, 𝐫, 𝐩) = dot(𝐫, 𝐫) / dot(𝐩, A, 𝐩)

compute_beta(A, 𝐫, 𝐩) = -dot(𝐩, A, 𝐫) / dot(𝐩, A, 𝐩)
compute_beta(𝐫ₙ₊₁, 𝐫ₙ) = dot(𝐫ₙ₊₁, 𝐫ₙ₊₁) / dot(𝐫ₙ, 𝐫ₙ)

isconverged(ch::ConvergenceHistory) = ch.isconverged

end
